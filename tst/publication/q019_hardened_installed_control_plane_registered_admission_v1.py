#!/usr/bin/env python3
"""Hardened registered-execution admission for one Q019 nonlinear Bell case.

The caller supplies only immutable path anchors.  This module re-reads paired
control-plane receipts, re-derives them with the installed producer, retains
the exact sealed raw bytes, validates restart publication, and independently
re-runs the Q019 raw reducer.  The result admits analysis of one execution but
does not grant nonlinear-saturation or publication authority.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import struct
import tarfile
from typing import Mapping

import numpy as np

from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design
from tst.publication import q019_registered_raw_reduction_v1 as raw_reduction
from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import (
    q023_registered_execution_linear_qualification_successor_v1 as q023,
)
from tst.publication import (
    q043_registered_execution_raw_oracle_qualification_successor_v1 as q043,
)


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_hardened_installed_control_plane_registered_admission_v1"
REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTERED_CAMPAIGN = "q019_nonlinear_bell_registered_successor_v1"
EXECUTION_RECEIPT_NAME = "q019_registered_execution_receipt.json"
TERMINAL_RECEIPT_NAME = "q019_terminal_receipt.json"
EXECUTION_RECEIPT_RECORD_TYPE = "q019_reconciled_registered_execution_receipt"
TERMINAL_RECEIPT_RECORD_TYPE = "q019_registered_execution_terminal_receipt"
PROJECT_HOME_MIRROR_NAMESPACE = Path("ledger/q019_registered_execution_receipts")
PROJECT_HOME_LEDGER_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
PRODUCER_ENTRYPOINT = "reconcile_q019_registered_execution.py"
DECK_MANIFEST_PATH = (
    "inputs/publication/q019_physics_first_nonlinear_bell_successor_v2/"
    "deck_manifest.json"
)
RUNTIME_CONTROLLER_MANIFEST_PATH = (
    "inputs/publication/q019_nonlinear_bell_runtime_controller_v1/"
    "deck_manifest.json"
)
RUNTIME_CONTROLLER_ROOT = "inputs/publication/q019_nonlinear_bell_runtime_controller_v1"
REQUIRED_SOURCE_PATHS = frozenset(
    {
        "src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.cpp",
        "src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.hpp",
        DECK_MANIFEST_PATH,
        "tst/publication/q019_physics_first_nonlinear_bell_successor_v2.py",
        "tst/publication/analyze_q019_physics_first_nonlinear_bell_successor_v2.py",
        "tst/publication/q019_finite_rigidity_early_time_physics_predecessor_v2.py",
        "tst/publication/q019_nonlinear_bell_particle_state.py",
        "tst/publication/q019_particle_state_analysis_bridge_v2.py",
        "tst/publication/q019_registered_raw_reduction_v1.py",
        "tst/publication/q019_nonlinear_bell_runtime_controller_v1.py",
        "tst/publication/q019_excluded_pilot_campaign_driver_v1.py",
        "tst/publication/q019_excluded_pilot_engineering_qualification_v1.py",
        "tst/publication/q019_excluded_pilot_launch_policy_preparation_v1.py",
        "tst/publication/q019_excluded_physical_pilot_launch_preparation_v1.py",
        "tst/publication/q019_hardened_provenance_boundary_v2.py",
        "tst/publication/q019_hardened_installed_control_plane_registered_admission_v1.py",
        "tst/publication/frontier_control_plane/reconcile_q019_registered_execution.py",
        "tst/publication/analyze_q011_section54_outputs.py",
        "tst/publication/q011_section54_restart.py",
        RUNTIME_CONTROLLER_MANIFEST_PATH,
        *{
            f"{RUNTIME_CONTROLLER_ROOT}/{overlay['artifact_id']}.athinput"
            for overlay in controller.expected_overlays()
        },
    }
)
EXECUTING_QUALIFICATION_SOURCE_PATHS = frozenset(
    {
        "tst/publication/q019_hardened_installed_control_plane_registered_admission_v1.py",
        "tst/publication/q019_physics_first_nonlinear_bell_successor_v2.py",
        "tst/publication/q019_registered_raw_reduction_v1.py",
        "tst/publication/q019_nonlinear_bell_runtime_controller_v1.py",
        "tst/publication/q019_nonlinear_bell_particle_state.py",
        "tst/publication/q019_particle_state_analysis_bridge_v2.py",
        "tst/publication/analyze_q011_section54_outputs.py",
        "tst/publication/q011_section54_restart.py",
        "tst/publication/q023_registered_execution_linear_qualification_successor_v1.py",
        "tst/publication/analyze_q023_paper_bell_linear_joverc.py",
        "tst/publication/analyze_q023_paper_bell_linear.py",
        "tst/publication/q023_registered_launch_policy_preparation_successor_v1.py",
        "tst/publication/q043_registered_launch_policy_preparation_successor_v1.py",
        "tst/publication/q043_registered_execution_raw_oracle_qualification_successor_v1.py",
        "tst/publication/q043_bell_current_volume_aware_deposited_current_oracle.py",
        "tst/publication/immutable_orion_tree.py",
        "tst/publication/frontier_control_plane/control_plane_common.py",
    }
)
AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "qualification_authorized": False,
    "raw_production_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}
_SHA256 = re.compile(r"[0-9a-f]{64}")
_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
_MARKER = re.compile(
    rb"ATHENAK_RESTART_COMPLETE_V1\nsize=([0-9]+)\nfnv1a64=([0-9a-f]{16})\n"
)
_RESTART_MESH_HEADER_FORMAT = "<ii9d19i19iddii"


class RegisteredAdmissionError(ValueError):
    """Reject substituted, mutable, incomplete, or authority-bearing evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RegisteredAdmissionError(message)


def _canonical_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise RegisteredAdmissionError("Q019 admission contains noncanonical JSON") from error


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _strict_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


def _relative_path(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected relative path")
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and path.as_posix() == value
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe relative path",
    )
    return value


def _stable_read_only(
    path: Path,
    *,
    root: Path,
    label: str,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
    executable: bool = False,
    read_only: bool = True,
) -> tuple[Path, bytes, dict[str, int]]:
    root = Path(os.path.abspath(root)).resolve(strict=True)
    lexical = Path(os.path.abspath(path))
    try:
        resolved = lexical.resolve(strict=True)
        lexical.relative_to(root)
    except (OSError, ValueError) as error:
        raise RegisteredAdmissionError(f"{label}: path escaped or is unavailable") from error
    _require(resolved == lexical, f"{label}: path contains a symlink or alias")
    descriptor = os.open(lexical, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and (not read_only or not before.st_mode & 0o222)
            and (not executable or bool(before.st_mode & 0o111)),
            f"{label}: expected one stable regular file with required mode",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = lexical.stat(follow_symlinks=False)
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
            and len(payload) == after.st_size,
            f"{label}: file changed while reading",
        )
        digest = hashlib.sha256(payload).hexdigest()
        _require(
            expected_sha256 is None or digest == expected_sha256,
            f"{label}: SHA-256 drifted",
        )
        _require(
            expected_size is None or len(payload) == expected_size,
            f"{label}: byte count drifted",
        )
        return lexical, payload, {"device": after.st_dev, "inode": after.st_ino}
    finally:
        os.close(descriptor)


def _json_object(payload: bytes, *, label: str) -> dict[str, object]:
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RegisteredAdmissionError(f"{label}: invalid UTF-8 JSON") from error
    _require(type(value) is dict, f"{label}: expected JSON object")
    return value


def _fnv1a64(payload: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in payload:
        value ^= byte
        value = (value * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return value


def _validate_marker(marker: bytes, payload: bytes, *, label: str) -> None:
    match = _MARKER.fullmatch(marker)
    _require(match is not None, f"{label}: malformed completion marker")
    _require(
        int(match.group(1)) == len(payload)
        and int(match.group(2), 16) == _fnv1a64(payload),
        f"{label}: completion marker does not bind payload",
    )


def _restart_cycle_time(payload: bytes, *, label: str) -> tuple[int, float]:
    header, marker, body = payload.partition(b"<par_end>\n")
    _require(
        marker == b"<par_end>\n" and bool(header) and bool(body),
        f"{label}: invalid restart payload",
    )
    offset = len(header) + len(marker)
    try:
        values = struct.unpack_from(_RESTART_MESH_HEADER_FORMAT, payload, offset)
    except struct.error as error:
        raise RegisteredAdmissionError(
            f"{label}: restart mesh header is truncated"
        ) from error
    nmb_total = values[0]
    root_level = values[1]
    time, timestep, cycle, original_nranks = values[-4:]
    _require(
        type(cycle) is int
        and cycle >= 0
        and type(nmb_total) is int
        and nmb_total > 0
        and type(root_level) is int
        and root_level >= 0
        and type(original_nranks) is int
        and 0 < original_nranks <= nmb_total
        and math.isfinite(time)
        and time >= 0.0
        and math.isfinite(timestep)
        and timestep > 0.0,
        f"{label}: restart mesh chronology or cardinality is invalid",
    )
    return int(cycle), float(time)


def _validate_restart_publication(
    *,
    stem: str,
    payloads: Mapping[str, bytes],
) -> tuple[int, float]:
    expected = {
        stem,
        stem + ".complete",
        stem + ".manifest",
        stem + ".manifest.complete",
    }
    _require(expected <= set(payloads), f"{stem}: restart publication is incomplete")
    restart = payloads[stem]
    marker = payloads[stem + ".complete"]
    manifest_payload = payloads[stem + ".manifest"]
    manifest_marker = payloads[stem + ".manifest.complete"]
    _validate_marker(marker, restart, label=stem + ".complete")
    _validate_marker(
        manifest_marker, manifest_payload, label=stem + ".manifest.complete"
    )
    manifest = _json_object(manifest_payload, label=stem + ".manifest")
    _require(
        set(manifest) == {"schema", "members"}
        and manifest["schema"] == "ATHENAK_RESTART_MANIFEST_V1"
        and type(manifest["members"]) is list
        and len(manifest["members"]) == 1,
        f"{stem}: restart manifest schema or shared layout drifted",
    )
    member = manifest["members"][0]
    _require(
        type(member) is dict
        and set(member) == {"path", "size", "fnv1a64"}
        and member["path"] == stem
        and type(member["size"]) is int
        and member["size"] == len(restart)
        and type(member["fnv1a64"]) is str
        and re.fullmatch(r"[0-9a-f]{16}", member["fnv1a64"]) is not None
        and int(member["fnv1a64"], 16) == _fnv1a64(restart),
        f"{stem}: restart manifest does not bind payload bytes",
    )
    return _restart_cycle_time(restart, label=stem)


def _case(case_id: str) -> dict[str, object]:
    matches = [row for row in design.expected_cases() if row["case_id"] == case_id]
    _require(len(matches) == 1, "Q019 registered case identity is unknown")
    case = dict(matches[0])
    deck_path = (
        "inputs/publication/q019_physics_first_nonlinear_bell_successor_v2/"
        f"{case_id}.athinput"
    )
    deck_payload = design.render_deck(case).encode("utf-8")
    case["path"] = deck_path
    case["sha256"] = hashlib.sha256(deck_payload).hexdigest()
    return case


def _source_archive_payloads(
    payload: bytes,
    *,
    required_paths: frozenset[str] = REQUIRED_SOURCE_PATHS,
) -> dict[str, bytes]:
    observed: dict[str, tarfile.TarInfo] = {}
    try:
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:*") as archive:
            for member in archive.getmembers():
                name = member.name.rstrip("/")
                path = PurePosixPath(name)
                _require(
                    bool(name)
                    and not path.is_absolute()
                    and "." not in path.parts
                    and ".." not in path.parts
                    and path.as_posix() == name
                    and (member.isdir() or member.isreg())
                    and name not in observed,
                    "Q019 source archive contains an unsafe or duplicate member",
                )
                observed[name] = member
            _require(
                required_paths <= set(observed),
                "Q019 source archive lacks the complete analysis/runtime closure",
            )
            retained = {}
            for relative in sorted(required_paths):
                member = observed[relative]
                _require(member.isreg(), f"Q019 source archive member is not regular: {relative}")
                extracted = archive.extractfile(member)
                _require(extracted is not None, f"Q019 source archive member is unreadable: {relative}")
                retained[relative] = extracted.read()
            return retained
    except (tarfile.TarError, OSError) as error:
        raise RegisteredAdmissionError("Q019 source archive is invalid") from error


def _payload_bindings(
    payloads: Mapping[str, bytes],
) -> dict[str, dict[str, object]]:
    return {
        relative: {
            "sha256": hashlib.sha256(payload).hexdigest(),
            "byte_count": len(payload),
        }
        for relative, payload in sorted(payloads.items())
    }


def _source_archive_bindings(payload: bytes) -> dict[str, dict[str, object]]:
    return _payload_bindings(_source_archive_payloads(payload))


def _candidate_analysis_closure(
    payloads: Mapping[str, bytes],
    *,
    case: Mapping[str, object],
    execution_deck_sha256: str | None = None,
) -> dict[str, object]:
    manifest = _json_object(
        payloads[DECK_MANIFEST_PATH],
        label="Q019 archived deck manifest",
    )
    analysis_bindings = manifest.get("analysis_bindings")
    _require(
        manifest.get("schema_version") == design.SCHEMA_VERSION
        and manifest.get("record_type")
        == "q019_physics_first_nonlinear_bell_successor_v2_deck_manifest"
        and type(analysis_bindings) is list,
        "Q019 archived deck manifest identity drifted",
    )
    expected_analysis_paths = list(design.ANALYSIS_BINDING_PATHS)
    _require(
        len(analysis_bindings) == len(expected_analysis_paths),
        "Q019 archived analysis binding inventory drifted",
    )
    normalized = []
    for index, expected_path in enumerate(expected_analysis_paths):
        item = analysis_bindings[index]
        _require(
            type(item) is dict
            and set(item) == {"path", "sha256"}
            and item["path"] == expected_path
            and type(item["sha256"]) is str
            and _SHA256.fullmatch(item["sha256"]) is not None
            and expected_path in payloads
            and item["sha256"]
            == hashlib.sha256(payloads[expected_path]).hexdigest(),
            f"Q019 archived analysis binding[{index}] drifted",
        )
        normalized.append(dict(item))
    deck_path = str(case["path"])
    decks = manifest.get("decks")
    _require(type(decks) is list, "Q019 archived deck inventory is malformed")
    matches = [
        item
        for item in decks
        if type(item) is dict and item.get("case_id") == case["case_id"]
    ]
    _require(
        len(matches) == 1
        and matches[0].get("path") == deck_path
        and matches[0].get("sha256") == case["sha256"]
        and hashlib.sha256(payloads[deck_path]).hexdigest() == case["sha256"],
        "Q019 archived registered deck binding drifted",
    )
    selected_deck_sha256 = (
        str(case["sha256"])
        if execution_deck_sha256 is None
        else execution_deck_sha256
    )
    _require(
        _SHA256.fullmatch(selected_deck_sha256) is not None,
        "Q019 execution deck digest is malformed",
    )
    if selected_deck_sha256 == case["sha256"]:
        execution_deck = {
            "kind": "base_matrix_case",
            "artifact_id": None,
            "authority": "matrix_case",
            "source_case_id": case["case_id"],
            "path": deck_path,
            "sha256": case["sha256"],
            "byte_count": len(payloads[deck_path]),
            "saturation_evidence_eligible": False,
        }
    else:
        runtime_manifest = _json_object(
            payloads[RUNTIME_CONTROLLER_MANIFEST_PATH],
            label="Q019 archived runtime-controller manifest",
        )
        _require(
            _strict_equal(runtime_manifest, controller.build_manifest()),
            "Q019 archived runtime-controller packet drifted",
        )
        overlay_matches = [
            item
            for item in runtime_manifest["artifacts"]
            if item["source_case_id"] == case["case_id"]
            and item["rendered_sha256"] == selected_deck_sha256
        ]
        _require(
            len(overlay_matches) == 1,
            "Q019 execution deck is not one exact runtime-controller overlay",
        )
        overlay = overlay_matches[0]
        overlay_path = f"{RUNTIME_CONTROLLER_ROOT}/{overlay['filename']}"
        _require(
            overlay["source_deck"] == deck_path
            and overlay["source_deck_sha256"] == case["sha256"]
            and overlay_path in payloads
            and hashlib.sha256(payloads[overlay_path]).hexdigest()
            == selected_deck_sha256,
            "Q019 runtime-controller overlay source binding drifted",
        )
        execution_deck = {
            "kind": "runtime_controller_overlay",
            "artifact_id": overlay["artifact_id"],
            "authority": overlay["authority"],
            "source_case_id": case["case_id"],
            "path": overlay_path,
            "sha256": selected_deck_sha256,
            "byte_count": len(payloads[overlay_path]),
            "controller_identity_fingerprint": overlay[
                "controller_parameters"
            ]["controller_identity_fingerprint"],
            "expected_stop_reason": overlay["expected_stop_reason"],
            "saturation_evidence_eligible": False,
        }
    executing = {}
    for relative in sorted(EXECUTING_QUALIFICATION_SOURCE_PATHS):
        _require(
            relative in payloads,
            f"Q019 archived qualification source is absent: {relative}",
        )
        _, current_payload, identity = _stable_read_only(
            REPO_ROOT / relative,
            root=REPO_ROOT,
            label=f"Q019 executing qualification source {relative}",
            read_only=False,
        )
        archived_sha256 = hashlib.sha256(payloads[relative]).hexdigest()
        _require(
            hashlib.sha256(current_payload).hexdigest() == archived_sha256,
            f"Q019 executing qualification source differs from candidate: {relative}",
        )
        executing[relative] = {
            "sha256": archived_sha256,
            "byte_count": len(current_payload),
            "filesystem_identity": identity,
        }
    return {
        "deck_manifest": {
            "path": DECK_MANIFEST_PATH,
            "sha256": hashlib.sha256(payloads[DECK_MANIFEST_PATH]).hexdigest(),
            "byte_count": len(payloads[DECK_MANIFEST_PATH]),
        },
        "analysis_bindings": normalized,
        "registered_deck": {
            "path": deck_path,
            "sha256": str(case["sha256"]),
            "byte_count": len(payloads[deck_path]),
        },
        "runtime_controller_manifest": {
            "path": RUNTIME_CONTROLLER_MANIFEST_PATH,
            "sha256": hashlib.sha256(
                payloads[RUNTIME_CONTROLLER_MANIFEST_PATH]
            ).hexdigest(),
            "byte_count": len(payloads[RUNTIME_CONTROLLER_MANIFEST_PATH]),
        },
        "execution_deck": execution_deck,
        "executing_qualification_source_bindings": executing,
    }


def _manifest_and_candidate(
    receipt: Mapping[str, object],
    *,
    case: Mapping[str, object],
    artifact_root: Path,
    authorized_orion_root: Path,
) -> dict[str, object]:
    manifest_path = Path(str(receipt.get("pre_submit_manifest_path", "")))
    manifest_path, manifest_payload, manifest_identity = _stable_read_only(
        manifest_path,
        root=authorized_orion_root,
        label="Q019 pre-submit manifest",
        expected_sha256=str(receipt.get("pre_submit_manifest_sha256", "")),
    )
    manifest = _json_object(manifest_payload, label="Q019 pre-submit manifest")
    _require(
        manifest.get("schema_version") == 1
        and manifest.get("campaign") == REGISTERED_CAMPAIGN
        and manifest.get("test_id") == case["case_id"]
        and manifest.get("submission_id") == receipt.get("submission_id")
        and manifest.get("artifact_dir") == str(artifact_root)
        and manifest.get("git_commit") == receipt.get("source_commit")
        and manifest.get("control_plane_version") == receipt.get("control_plane_version")
        and manifest.get("registered_science_authorization_id")
        == receipt.get("registered_science_authorization_id")
        and manifest.get("clean_candidate_manifest_sha256")
        == receipt.get("clean_candidate_manifest_sha256"),
        "Q019 pre-submit manifest differs from registered execution",
    )
    snapshots = manifest.get("snapshot_files")
    _require(type(snapshots) is list, "Q019 pre-submit snapshot inventory is malformed")
    by_role = {
        item["role"]: item
        for item in snapshots
        if type(item) is dict and type(item.get("role")) is str
    }
    _require(
        {"executable", "input-deck", "clean-candidate-manifest", "environment-profile"}
        <= set(by_role),
        "Q019 pre-submit manifest lacks required snapshots",
    )
    deck_snapshot = by_role["input-deck"]
    _require(
        deck_snapshot.get("sha256") == receipt.get("deck_sha256")
        and type(receipt.get("deck_sha256")) is str
        and _SHA256.fullmatch(str(receipt["deck_sha256"])) is not None,
        "Q019 exact input deck binding drifted",
    )
    candidate_manifest_path = Path(str(manifest.get("clean_candidate_manifest_path", "")))
    candidate_manifest_path, candidate_payload, candidate_identity = _stable_read_only(
        candidate_manifest_path,
        root=authorized_orion_root,
        label="Q019 clean-candidate manifest",
        expected_sha256=str(receipt.get("clean_candidate_manifest_sha256", "")),
    )
    candidate = _json_object(candidate_payload, label="Q019 clean-candidate manifest")
    source = candidate.get("source")
    build = candidate.get("build")
    candidate_root = candidate_manifest_path.parent
    _require(
        candidate.get("schema_version") == 4
        and candidate_root.parent == authorized_orion_root / "clean_candidates"
        and _UUID.fullmatch(candidate_root.name) is not None
        and type(source) is dict
        and type(build) is dict
        and source.get("worktree_status") == "clean"
        and source.get("git_commit") == receipt.get("source_commit")
        and source.get("source_bundle_sha256") == receipt.get("source_bundle_sha256")
        and source.get("archive_sha256") == receipt.get("source_archive_sha256")
        and build.get("source_archive_sha256") == receipt.get("source_archive_sha256")
        and build.get("source_bundle_sha256") == receipt.get("source_bundle_sha256")
        and build.get("executable_sha256") == receipt.get("executable_sha256"),
        "Q019 clean-candidate lineage drifted",
    )
    archive_path = candidate_root / "source.tar"
    _, archive_payload, archive_identity = _stable_read_only(
        archive_path,
        root=authorized_orion_root,
        label="Q019 clean-candidate source archive",
        expected_sha256=str(receipt.get("source_archive_sha256", "")),
    )
    executable_snapshot = manifest_path.parent / "snapshot/athena"
    _, executable_payload, executable_identity = _stable_read_only(
        executable_snapshot,
        root=authorized_orion_root,
        label="Q019 immutable executable snapshot",
        expected_sha256=str(receipt.get("executable_sha256", "")),
        executable=True,
    )
    _require(
        by_role["executable"].get("sha256") == hashlib.sha256(executable_payload).hexdigest()
        and by_role["environment-profile"].get("sha256")
        == receipt.get("environment_sha256"),
        "Q019 executable or environment snapshot binding drifted",
    )
    required_source_paths = frozenset(
        {
            *REQUIRED_SOURCE_PATHS,
            *EXECUTING_QUALIFICATION_SOURCE_PATHS,
            *design.ANALYSIS_BINDING_PATHS,
            str(case["path"]),
        }
    )
    source_payloads = _source_archive_payloads(
        archive_payload,
        required_paths=required_source_paths,
    )
    analysis_closure = _candidate_analysis_closure(
        source_payloads,
        case=case,
        execution_deck_sha256=str(receipt["deck_sha256"]),
    )
    deck_snapshot_path = Path(str(deck_snapshot.get("path", "")))
    _, deck_snapshot_payload, deck_snapshot_identity = _stable_read_only(
        deck_snapshot_path,
        root=authorized_orion_root,
        label="Q019 immutable input-deck snapshot",
        expected_sha256=str(receipt["deck_sha256"]),
    )
    _require(
        deck_snapshot_payload
        == source_payloads[analysis_closure["execution_deck"]["path"]],
        "Q019 input-deck snapshot differs from the archived execution deck",
    )
    return {
        "pre_submit_manifest": {
            "path": str(manifest_path),
            "sha256": hashlib.sha256(manifest_payload).hexdigest(),
            "byte_count": len(manifest_payload),
            "filesystem_identity": manifest_identity,
        },
        "clean_candidate_manifest": {
            "path": str(candidate_manifest_path),
            "sha256": hashlib.sha256(candidate_payload).hexdigest(),
            "byte_count": len(candidate_payload),
            "filesystem_identity": candidate_identity,
        },
        "source_archive": {
            "path": str(archive_path),
            "sha256": hashlib.sha256(archive_payload).hexdigest(),
            "byte_count": len(archive_payload),
            "filesystem_identity": archive_identity,
        },
        "source_bindings": _payload_bindings(source_payloads),
        "analysis_closure": analysis_closure,
        "executable_snapshot": {
            "path": str(executable_snapshot),
            "sha256": hashlib.sha256(executable_payload).hexdigest(),
            "byte_count": len(executable_payload),
            "filesystem_identity": executable_identity,
        },
        "deck_snapshot": {
            **dict(deck_snapshot),
            "byte_count": len(deck_snapshot_payload),
            "filesystem_identity": deck_snapshot_identity,
        },
    }


def _receipt_pair(
    *,
    artifact_root: Path,
    authorized_orion_root: Path,
    authorized_project_home_root: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    receipt_path, receipt_payload, receipt_identity = _stable_read_only(
        artifact_root / "analysis" / EXECUTION_RECEIPT_NAME,
        root=authorized_orion_root,
        label="Q019 Orion execution receipt",
    )
    receipt = _json_object(receipt_payload, label="Q019 execution receipt")
    submission_id = receipt.get("submission_id")
    _require(
        receipt.get("record_type") == EXECUTION_RECEIPT_RECORD_TYPE
        and type(submission_id) is str
        and _UUID.fullmatch(submission_id) is not None,
        "Q019 execution receipt identity drifted",
    )
    mirror_root = authorized_project_home_root / PROJECT_HOME_MIRROR_NAMESPACE / submission_id
    receipt_mirror = mirror_root / EXECUTION_RECEIPT_NAME
    terminal_mirror = mirror_root / TERMINAL_RECEIPT_NAME
    _, receipt_mirror_payload, receipt_mirror_identity = _stable_read_only(
        receipt_mirror,
        root=authorized_project_home_root,
        label="Q019 Project Home execution receipt",
    )
    _require(
        receipt_mirror_payload == receipt_payload,
        "Q019 Orion and Project Home execution receipts differ",
    )
    terminal_path, terminal_payload, terminal_identity = _stable_read_only(
        artifact_root / "analysis" / TERMINAL_RECEIPT_NAME,
        root=authorized_orion_root,
        label="Q019 Orion terminal receipt",
        expected_sha256=str(receipt.get("terminal_receipt_sha256", "")),
    )
    terminal = _json_object(terminal_payload, label="Q019 terminal receipt")
    _, terminal_mirror_payload, terminal_mirror_identity = _stable_read_only(
        terminal_mirror,
        root=authorized_project_home_root,
        label="Q019 Project Home terminal receipt",
    )
    _require(
        terminal_mirror_payload == terminal_payload
        and terminal.get("record_type") == TERMINAL_RECEIPT_RECORD_TYPE
        and terminal.get("submission_id") == submission_id
        and terminal.get("member_id") == receipt.get("member_id")
        and terminal.get("slurm_terminal_state") == "COMPLETED"
        and terminal.get("slurm_exit_code") == "0:0"
        and terminal.get("raw_inventory_sha256") == receipt.get("raw_inventory_sha256"),
        "Q019 terminal and execution receipts do not cross-link",
    )
    return receipt, {
        "execution_receipt": {
            "orion_path": str(receipt_path),
            "project_home_path": str(receipt_mirror),
            "sha256": hashlib.sha256(receipt_payload).hexdigest(),
            "byte_count": len(receipt_payload),
            "orion_filesystem_identity": receipt_identity,
            "project_home_filesystem_identity": receipt_mirror_identity,
        },
        "terminal_receipt": {
            "orion_path": str(terminal_path),
            "project_home_path": str(terminal_mirror),
            "sha256": hashlib.sha256(terminal_payload).hexdigest(),
            "byte_count": len(terminal_payload),
            "orion_filesystem_identity": terminal_identity,
            "project_home_filesystem_identity": terminal_mirror_identity,
        },
        "terminal": terminal,
        "_receipt_payload": receipt_payload,
        "_terminal_payload": terminal_payload,
    }


def _installed_rederivation(
    receipt: Mapping[str, object],
    paired: Mapping[str, object],
    *,
    artifact_root: Path,
    authorized_orion_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    """Re-derive exact receipt bytes from mirrored ledger and installed code."""
    version = receipt.get("control_plane_version")
    _require(
        type(version) is str and _SHA256.fullmatch(version) is not None,
        "Q019 installed producer version is malformed",
    )
    try:
        with q043._installed_control_plane_modules(
            version, ("ledger.py", PRODUCER_ENTRYPOINT)
        ) as (modules, pair):
            ledger = modules["ledger.py"]
            producer = modules[PRODUCER_ENTRYPOINT]
            ledger_jsonl = authorized_orion_root / "ledger/node_hours.jsonl"
            receipts_jsonl = authorized_orion_root / "ledger/mirror_receipts.jsonl"
            mirror_jsonl = PROJECT_HOME_LEDGER_ROOT / "ledger/node_hours.jsonl"
            _require(
                PROJECT_HOME_LEDGER_ROOT.resolve(strict=True)
                == authorized_project_home_root.resolve(strict=True),
                "Q019 canonical Project Home ledger alias drifted",
            )
            with ledger.validated_read_only_mirrored_state_snapshot(
                ledger_jsonl,
                receipts_jsonl,
                mirror_jsonl,
                ledger_root=authorized_orion_root,
                receipts_root=authorized_orion_root,
                mirror_root=PROJECT_HOME_LEDGER_ROOT,
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
            event_matches = [
                item
                for item in events
                if item.get("event_type") == "reconciliation"
                and item.get("event_sha256")
                == receipt.get("reconciliation_event_sha256")
            ]
            _require(
                len(event_matches) == 1,
                "Q019 receipt lacks one canonical mirrored-ledger reconciliation",
            )
            event = event_matches[0]
            acknowledgment_matches = [
                item
                for item in acknowledgments
                if item.get("mirrored_event_sha256") == event["event_sha256"]
            ]
            _require(
                len(acknowledgment_matches) == 1,
                "Q019 receipt lacks one canonical mirror acknowledgment",
            )
            mirror_ack = acknowledgment_matches[0]
            expected_producer = {
                "entrypoint": PRODUCER_ENTRYPOINT,
                "entrypoint_sha256": pair["digests"].get(PRODUCER_ENTRYPOINT),
                "launch_trampoline_sha256": pair["digests"].get(
                    "launch_trampoline.py"
                ),
                "control_plane_version": version,
            }
            _require(
                receipt.get("producer") == expected_producer
                and event.get("campaign") == REGISTERED_CAMPAIGN
                and event.get("test_id") == receipt.get("member_id")
                and event.get("submission_id") == receipt.get("submission_id")
                and event.get("reservation_id") == receipt.get("reservation_id")
                and event.get("job_id") == receipt.get("slurm_job_id")
                and event.get("artifact_dir") == str(artifact_root)
                and event.get("state") == "COMPLETED"
                and event.get("scheduler_exit_code") == "0:0"
                and event.get("reconciled") is True
                and event.get("reconciled_by_control_plane_version") == version
                and mirror_ack.get("mirror_ack_sha256")
                == receipt.get("reconciliation_mirror_ack_sha256"),
                "Q019 receipt differs from canonical reconciliation evidence",
            )
            derived = producer.derive_q019_registered_execution_evidence(
                event,
                mirror_ack,
                pair["inventory"],
                authorized_pic_root=authorized_orion_root,
                authorized_project_home_root=authorized_project_home_root,
            )
    except (AttributeError, KeyError, OSError, TypeError, ValueError) as error:
        raise RegisteredAdmissionError(
            "Q019 installed producer could not re-derive registered evidence"
        ) from error
    receipt_evidence = paired["execution_receipt"]
    terminal_evidence = paired["terminal_receipt"]
    expected = (
        Path(str(terminal_evidence["orion_path"])),
        paired["_terminal_payload"],
        Path(str(receipt_evidence["orion_path"])),
        paired["_receipt_payload"],
        Path(str(terminal_evidence["project_home_path"])),
        Path(str(receipt_evidence["project_home_path"])),
    )
    _require(
        derived == expected,
        "installed Q019 producer re-derived different evidence bytes or paths",
    )
    return {
        "control_plane_version": version,
        "entrypoint": PRODUCER_ENTRYPOINT,
        "entrypoint_sha256": expected_producer["entrypoint_sha256"],
        "launch_trampoline_sha256": expected_producer[
            "launch_trampoline_sha256"
        ],
        "reconciliation_event_sha256": event["event_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "exact_byte_rederivation_passed": True,
    }


def _dependency_binding(
    path: Path, *, root: Path, label: str
) -> tuple[dict[str, object], dict[str, object]]:
    path, payload, identity = _stable_read_only(path, root=root, label=label)
    return _json_object(payload, label=label), {
        "path": str(path),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "byte_count": len(payload),
        "filesystem_identity": identity,
    }


def _array_binding(value: object) -> dict[str, object]:
    array = np.ascontiguousarray(np.asarray(value))
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "sha256": hashlib.sha256(array.tobytes(order="C")).hexdigest(),
    }


def _reduction_binding(reduction: Mapping[str, object]) -> dict[str, object]:
    snapshots = []
    for snapshot in reduction["snapshots"]:
        fields = snapshot["fields"]
        snapshots.append(
            {
                "cycle": snapshot["cycle"],
                "time": snapshot["time"],
                "x1_faces": _array_binding(snapshot["x1_faces"]),
                "x2_faces": _array_binding(snapshot["x2_faces"]),
                "x3_faces": _array_binding(snapshot["x3_faces"]),
                "fields": {
                    name: _array_binding(fields[name])
                    for name in sorted(fields)
                },
            }
        )
    particle_states = list(reduction["particle_states"])
    return {
        "record_type": reduction["record_type"],
        "case_id": reduction["case_id"],
        "campaign_id": reduction["campaign_id"],
        "matched_checkpoint_count": reduction["matched_checkpoint_count"],
        "chronology": reduction["chronology"],
        "reference_budget": reduction["reference_budget"],
        "execution_profile": reduction["execution_profile"],
        "runtime_controller_states_sha256": canonical_sha256(
            reduction["runtime_controller_states"]
        ),
        "snapshots": snapshots,
        "particle_states_sha256": canonical_sha256(particle_states),
    }


def _raw_bundle(
    receipt: Mapping[str, object],
    *,
    case_id: str,
    artifact_root: Path,
    authorized_orion_root: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    inventory = receipt.get("raw_inventory")
    _require(type(inventory) is list and inventory, "Q019 raw inventory is absent")
    _require(
        canonical_sha256(inventory) == receipt.get("raw_inventory_sha256"),
        "Q019 raw inventory digest drifted",
    )
    raw_root = Path(str(receipt.get("raw_output_root", "")))
    _require(
        raw_root == artifact_root / "raw"
        and raw_root.resolve(strict=True) == raw_root
        and not raw_root.stat(follow_symlinks=False).st_mode & 0o222,
        "Q019 raw root is not the sealed registered raw directory",
    )
    payloads: dict[str, bytes] = {}
    bindings = []
    identities = set()
    for index, item in enumerate(inventory):
        label = f"Q019 raw inventory[{index}]"
        _require(
            type(item) is dict
            and set(item)
            == {"path", "sha256", "byte_count", "member_id", "artifact_kind"}
            and item["member_id"] == case_id
            and type(item["sha256"]) is str
            and _SHA256.fullmatch(item["sha256"]) is not None
            and type(item["byte_count"]) is int
            and item["byte_count"] > 0,
            f"{label}: schema drifted",
        )
        relative = _relative_path(item["path"], label=f"{label}/path")
        path, payload, identity = _stable_read_only(
            raw_root / relative,
            root=authorized_orion_root,
            label=label,
            expected_sha256=item["sha256"],
            expected_size=item["byte_count"],
        )
        identity_pair = (identity["device"], identity["inode"])
        _require(identity_pair not in identities, "Q019 raw filesystem object reused")
        identities.add(identity_pair)
        payloads[relative] = payload
        bindings.append(
            {
                "path": relative,
                "sha256": item["sha256"],
                "byte_count": item["byte_count"],
                "filesystem_identity": identity,
                "absolute_path": str(path),
            }
        )
    binary_pattern = re.compile(
        rf"bin/{re.escape(case_id)}[.]({'|'.join(raw_reduction.REQUIRED_BINARY_PRODUCTS)})"
        rf"[.]([0-9]{{5,}})[.]bin"
    )
    restart_pattern = re.compile(
        rf"rst/{re.escape(case_id)}[.]([0-9]{{5,}})[.]rst"
    )
    binary_by_index: dict[int, dict[str, bytes]] = {}
    restart_stems: dict[int, str] = {}
    for relative, payload in payloads.items():
        if match := binary_pattern.fullmatch(relative):
            binary_by_index.setdefault(int(match.group(2)), {})[match.group(1)] = payload
        elif match := restart_pattern.fullmatch(relative):
            restart_stems[int(match.group(1))] = relative
    _require(
        sorted(binary_by_index) == list(receipt.get("binary_output_indices", []))
        and sorted(restart_stems) == list(receipt.get("checkpoint_output_indices", [])),
        "Q019 raw chronology indices differ from reconciled receipt",
    )
    binary_metadata: dict[tuple[int, float], tuple[int, dict[str, bytes]]] = {}
    for output_index in sorted(binary_by_index):
        snapshot = raw_reduction.compose_snapshot(
            case_id, binary_by_index[output_index]
        )
        key = (int(snapshot["cycle"]), float(snapshot["time"]))
        _require(key not in binary_metadata, "Q019 binary cycle/time is duplicated")
        binary_metadata[key] = (output_index, binary_by_index[output_index])
    checkpoints = []
    checkpoint_bindings = []
    for output_index in sorted(restart_stems):
        stem = restart_stems[output_index]
        cycle, time = _validate_restart_publication(stem=stem, payloads=payloads)
        matches = [
            (key, value)
            for key, value in binary_metadata.items()
            if key[0] == cycle
            and math.isclose(key[1], time, rel_tol=0.0, abs_tol=1.0e-12)
        ]
        _require(
            len(matches) == 1,
            f"{stem}: no unique ten-product binary snapshot matches restart cycle/time",
        )
        (_, matched_time), (binary_index, products) = matches[0]
        checkpoints.append(
            {
                "cycle": cycle,
                "time": float(matched_time),
                "binary_products": products,
                "restart_path": stem,
                "restart_payload": payloads[stem],
            }
        )
        checkpoint_bindings.append(
            {
                "checkpoint_output_index": output_index,
                "binary_output_index": binary_index,
                "cycle": cycle,
                "time": float(matched_time),
                "restart_path": stem,
                "restart_sha256": hashlib.sha256(payloads[stem]).hexdigest(),
            }
        )
    reduction = raw_reduction.reduce_matched_checkpoints(case_id, checkpoints)
    terminal = (
        checkpoints[-1]["cycle"],
        checkpoints[-1]["time"],
    )
    _require(
        terminal[0] == receipt.get("terminal_cycle")
        and math.isclose(
            terminal[1], float(receipt.get("terminal_time")), rel_tol=0.0, abs_tol=1.0e-12
        ),
        "Q019 terminal matched checkpoint differs from registered receipt",
    )
    return reduction, {
        "raw_root": str(raw_root),
        "raw_inventory_sha256": receipt["raw_inventory_sha256"],
        "retained_raw_bindings": bindings,
        "checkpoint_bindings": checkpoint_bindings,
        "reduction_binding": _reduction_binding(reduction),
    }


def _execution_lineage(
    receipt: Mapping[str, object],
    *,
    paired: Mapping[str, object],
    installed_rederivation: Mapping[str, object],
    candidate: Mapping[str, object],
    raw_binding: Mapping[str, object],
) -> dict[str, object]:
    """Project the immutable execution lineage needed by campaign indexes."""
    source_archive = candidate.get("source_archive")
    executable = candidate.get("executable_snapshot")
    pre_submit = candidate.get("pre_submit_manifest")
    clean_candidate = candidate.get("clean_candidate_manifest")
    execution_receipt = paired.get("execution_receipt")
    terminal_receipt = paired.get("terminal_receipt")
    producer = receipt.get("producer")
    _require(
        all(
            type(value) is dict
            for value in (
                source_archive,
                executable,
                pre_submit,
                clean_candidate,
                execution_receipt,
                terminal_receipt,
                producer,
            )
        )
        and source_archive.get("sha256") == receipt.get("source_archive_sha256")
        and executable.get("sha256") == receipt.get("executable_sha256")
        and pre_submit.get("sha256") == receipt.get("pre_submit_manifest_sha256")
        and clean_candidate.get("sha256")
        == receipt.get("clean_candidate_manifest_sha256")
        and execution_receipt.get("sha256")
        == hashlib.sha256(paired["_receipt_payload"]).hexdigest()
        and terminal_receipt.get("sha256")
        == hashlib.sha256(paired["_terminal_payload"]).hexdigest()
        and installed_rederivation.get("control_plane_version")
        == receipt.get("control_plane_version")
        and installed_rederivation.get("entrypoint_sha256")
        == producer.get("entrypoint_sha256")
        and installed_rederivation.get("launch_trampoline_sha256")
        == producer.get("launch_trampoline_sha256")
        and installed_rederivation.get("reconciliation_event_sha256")
        == receipt.get("reconciliation_event_sha256")
        and installed_rederivation.get("reconciliation_mirror_ack_sha256")
        == receipt.get("reconciliation_mirror_ack_sha256")
        and raw_binding.get("raw_inventory_sha256")
        == receipt.get("raw_inventory_sha256")
        and type(raw_binding.get("retained_raw_bindings")) is list
        and type(raw_binding.get("reduction_binding")) is dict,
        "Q019 projected execution lineage differs from rederived evidence",
    )
    return {
        "source_commit": receipt["source_commit"],
        "source_bundle_sha256": receipt["source_bundle_sha256"],
        "source_archive_sha256": receipt["source_archive_sha256"],
        "executable_sha256": receipt["executable_sha256"],
        "environment_sha256": receipt["environment_sha256"],
        "deck_sha256": receipt["deck_sha256"],
        "pre_submit_manifest": dict(pre_submit),
        "clean_candidate_manifest": dict(clean_candidate),
        "source_archive": dict(source_archive),
        "executable_snapshot": dict(executable),
        "execution_receipt": dict(execution_receipt),
        "terminal_receipt": dict(terminal_receipt),
        "installed_producer": {
            "control_plane_version": installed_rederivation[
                "control_plane_version"
            ],
            "entrypoint": installed_rederivation["entrypoint"],
            "entrypoint_sha256": installed_rederivation["entrypoint_sha256"],
            "launch_trampoline_sha256": installed_rederivation[
                "launch_trampoline_sha256"
            ],
            "exact_byte_rederivation_passed": installed_rederivation[
                "exact_byte_rederivation_passed"
            ],
        },
        "reconciliation_event_sha256": receipt["reconciliation_event_sha256"],
        "reconciliation_mirror_ack_sha256": receipt[
            "reconciliation_mirror_ack_sha256"
        ],
        "raw_inventory_sha256": raw_binding["raw_inventory_sha256"],
        "retained_raw_bindings": list(raw_binding["retained_raw_bindings"]),
        "reduction_binding": dict(raw_binding["reduction_binding"]),
    }


def _completion_record(
    receipt: Mapping[str, object],
    *,
    execution_profile: Mapping[str, object] | None = None,
    runtime_controller_states: object = None,
) -> dict[str, object]:
    command = receipt.get("command_evidence")
    _require(type(command) is dict, "Q019 command evidence is malformed")
    wrapper = command.get("trusted_wrapper_evidence")
    _require(type(wrapper) is dict, "Q019 trusted wrapper evidence is absent")
    eligible = wrapper.get("problem_saturation_evidence_eligible")
    status = wrapper.get("problem_final_evidence_status")
    _require(
        eligible in {"true", "false"} and type(status) is str and bool(status),
        "Q019 final problem status is malformed",
    )
    user_stop = wrapper.get("termination_reason") == "Terminating on user request"
    stop_reason = wrapper.get("termination_reason")
    trigger_cycle = None
    trigger_time = None
    trigger_metric = None
    if execution_profile is not None:
        states = runtime_controller_states
        _require(
            execution_profile.get("kind") == "runtime_controller_overlay"
            and status == "completed_not_acceptance_eligible"
            and eligible == "false"
            and type(states) is list
            and bool(states)
            and all(type(item) is dict for item in states),
            "Q019 runtime-controller completion status or state is invalid",
        )
        final_state = states[-1]
        expected_reason = execution_profile.get("expected_stop_reason")
        if expected_reason is None:
            _require(
                stop_reason == "Terminating on time limit"
                and all(
                    state.get("runtime_controller_triggered") is False
                    and state.get("runtime_controller_trigger_failure") is False
                    and state.get("runtime_controller_trigger_reason") == 0
                    and state.get("runtime_controller_trigger_cycle") == 0
                    and state.get("runtime_controller_trigger_time") == 0.0
                    and state.get("runtime_controller_trigger_metric") == -1.0
                    for state in states
                ),
                "Q019 monitor-only controller completion state drifted",
            )
        else:
            _require(
                user_stop
                and final_state.get("runtime_controller_triggered") is True
                and final_state.get("runtime_controller_trigger_failure") is False
                and final_state.get("runtime_controller_trigger_reason")
                == expected_reason
                and type(expected_reason) is int
                and type(final_state.get("runtime_controller_trigger_cycle")) is int
                and final_state["runtime_controller_trigger_cycle"] > 0
                and type(final_state.get("runtime_controller_trigger_time")) is float
                and math.isfinite(final_state["runtime_controller_trigger_time"])
                and type(final_state.get("runtime_controller_trigger_metric")) is float
                and math.isfinite(final_state["runtime_controller_trigger_metric"]),
                "Q019 runtime-controller completion reason drifted",
            )
            stop_reason = str(expected_reason)
            trigger_cycle = final_state["runtime_controller_trigger_cycle"]
            trigger_time = final_state["runtime_controller_trigger_time"]
            trigger_metric = final_state["runtime_controller_trigger_metric"]
    return {
        "record_type": "q019_runtime_completion_status_v1",
        "run_completion_status": status,
        "problem_stop_requested": user_stop,
        "stop_reason_code": stop_reason,
        "runtime_controller_trigger_cycle": trigger_cycle,
        "runtime_controller_trigger_time": trigger_time,
        "runtime_controller_trigger_metric": trigger_metric,
        "process_exit_code": 0,
        "scheduler_terminal_state": receipt.get("slurm_terminal_state"),
        "trusted_execution_binding_present": True,
        "problem_saturation_evidence_eligible": eligible == "true",
    }


def derive_case_bundle(
    *,
    case_id: str,
    artifact_root: Path,
    q043_qualification_path: Path,
    q043_artifact_root: Path,
    q023_qualification_path: Path,
    authorized_orion_root: Path = q043.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = q043.AUTHORIZED_PROJECT_HOME_ROOT,
) -> tuple[dict[str, object], dict[str, object]]:
    """Rebuild one registered Q019 admission and its exact analysis inputs."""
    authorized_orion_root = Path(os.path.abspath(authorized_orion_root)).resolve(strict=True)
    authorized_project_home_root = Path(
        os.path.abspath(authorized_project_home_root)
    ).resolve(strict=True)
    artifact_root = Path(os.path.abspath(artifact_root)).resolve(strict=True)
    _require(
        artifact_root.parent
        == authorized_orion_root / "runs" / REGISTERED_CAMPAIGN
        and _UUID.fullmatch(artifact_root.name) is not None
        and not artifact_root.stat(follow_symlinks=False).st_mode & 0o222,
        "Q019 artifact root is not one sealed registered submission",
    )
    case = _case(case_id)
    q043_record, q043_binding = _dependency_binding(
        q043_qualification_path,
        root=authorized_orion_root,
        label="Q043 registered matrix qualification",
    )
    try:
        q043_dependency = q043.validate_downstream_q023_q019_prerequisite(q043_record)
    except q043.AdmissionError as error:
        raise RegisteredAdmissionError("Q019 Q043 dependency is not passing") from error
    q023_record, q023_binding = _dependency_binding(
        q023_qualification_path,
        root=authorized_orion_root,
        label="Q023 registered matrix qualification",
    )
    try:
        q023_dependency = q023.validate_downstream_q019_prerequisite(
            q023_record,
            q043_artifact_root=q043_artifact_root,
            authorized_orion_root=authorized_orion_root,
            authorized_project_home_root=authorized_project_home_root,
        )
    except q023.QualificationError as error:
        raise RegisteredAdmissionError("Q019 Q023 dependency is not passing") from error
    q043_artifact_root = Path(os.path.abspath(q043_artifact_root)).resolve(strict=True)
    try:
        q043_relative = Path(str(q043_binding["path"])).relative_to(
            q043_artifact_root
        ).as_posix()
    except ValueError as error:
        raise RegisteredAdmissionError(
            "Q019 Q043 matrix is outside its declared artifact root"
        ) from error
    embedded_q043 = q023_dependency.get("q043_dependency")
    _require(
        type(embedded_q043) is dict
        and embedded_q043.get("binding_kind") == "registered_matrix_qualification"
        and embedded_q043.get("registered_matrix_path") == q043_relative
        and embedded_q043.get("registered_matrix_sha256") == q043_binding["sha256"]
        and embedded_q043.get("registered_matrix_record_type")
        == q043_dependency["record_type"]
        and embedded_q043.get("registered_matrix_case_bindings_sha256")
        == q043_dependency["case_bindings_sha256"],
        "Q019 Q023 qualification is not bound to the selected Q043 matrix",
    )
    receipt, paired = _receipt_pair(
        artifact_root=artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    installed_rederivation = _installed_rederivation(
        receipt,
        paired,
        artifact_root=artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _require(
        receipt.get("member_id") == case_id
        and receipt.get("campaign_id") == "Q019-PHYSICS-FIRST-NONLINEAR-BELL"
        and receipt.get("slurm_terminal_state") == "COMPLETED"
        and receipt.get("slurm_exit_code") == "0:0",
        "Q019 registered execution identity or scheduler outcome drifted",
    )
    candidate = _manifest_and_candidate(
        receipt,
        case=case,
        artifact_root=artifact_root,
        authorized_orion_root=authorized_orion_root,
    )
    reduction, raw_binding = _raw_bundle(
        receipt,
        case_id=case_id,
        artifact_root=artifact_root,
        authorized_orion_root=authorized_orion_root,
    )
    execution_deck = candidate["analysis_closure"]["execution_deck"]
    execution_profile = reduction["execution_profile"]
    _require(
        execution_deck["kind"] == execution_profile["kind"]
        and execution_deck["source_case_id"]
        == execution_profile["source_case_id"]
        and execution_deck["artifact_id"] == execution_profile["artifact_id"]
        and execution_deck["authority"] == execution_profile["authority"],
        "Q019 archived execution deck and raw runtime profile differ",
    )
    if execution_deck["kind"] == "runtime_controller_overlay":
        _require(
            execution_deck["controller_identity_fingerprint"]
            == execution_profile["controller_identity_fingerprint"]
            and execution_deck["expected_stop_reason"]
            == execution_profile["expected_stop_reason"],
            "Q019 runtime-controller identity or stop contract drifted",
        )
    completion = _completion_record(
        receipt,
        execution_profile=(
            execution_profile
            if execution_profile["kind"] == "runtime_controller_overlay"
            else None
        ),
        runtime_controller_states=reduction.get("runtime_controller_states"),
    )
    public_paired = {
        key: value for key, value in paired.items() if not key.startswith("_")
    }
    execution_lineage = _execution_lineage(
        receipt,
        paired=paired,
        installed_rederivation=installed_rederivation,
        candidate=candidate,
        raw_binding=raw_binding,
    )
    admission = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": "registered_Q019_case_raw_analysis_admitted_non_authorizing",
        "case_id": case_id,
        "campaign_id": case["campaign_id"],
        "artifact_root": str(artifact_root),
        "q043_dependency": q043_binding,
        "q043_artifact_root": str(
            q043_artifact_root
        ),
        "q043_dependency_sha256": canonical_sha256(q043_dependency),
        "q023_dependency": q023_binding,
        "q023_dependency_sha256": canonical_sha256(q023_dependency),
        "registered_execution_receipt_sha256": public_paired[
            "execution_receipt"
        ]["sha256"],
        "registered_execution_identity": {
            key: receipt[key]
            for key in (
                "reservation_id",
                "submission_id",
                "reconciliation_event_sha256",
                "reconciliation_mirror_ack_sha256",
                "control_plane_version",
                "registered_science_authorization_id",
                "slurm_job_id",
                "pre_submit_manifest_sha256",
            )
        },
        "paired_control_plane_evidence": public_paired,
        "installed_reconciliation_rederivation": installed_rederivation,
        "execution_lineage": execution_lineage,
        "candidate_binding": candidate,
        "execution_deck": execution_deck,
        "execution_profile": execution_profile,
        "raw_binding": raw_binding,
        "runtime_completion": completion,
        "raw_science_admission_eligible": True,
        "problem_reported_saturation_evidence_eligible": completion[
            "problem_saturation_evidence_eligible"
        ],
        "saturation_evidence_eligible": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }
    return admission, {
        "snapshots": reduction["snapshots"],
        "particle_states": reduction["particle_states"],
        "completion_record": {
            key: completion[key]
            for key in (
                "record_type",
                "run_completion_status",
                "problem_stop_requested",
                "stop_reason_code",
                "process_exit_code",
                "scheduler_terminal_state",
                "trusted_execution_binding_present",
            )
        },
    }


def validate_admission(value: object) -> dict[str, object]:
    """Rebuild an admission from its immutable path anchors."""
    record, _ = validate_analysis_bundle(value)
    return record


def _analysis_equal(left: object, right: object) -> bool:
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        if not isinstance(left, np.ndarray) or not isinstance(right, np.ndarray):
            return False
        return (
            left.dtype == right.dtype
            and left.shape == right.shape
            and left.tobytes(order="C") == right.tobytes(order="C")
        )
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            _analysis_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _analysis_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


def validate_analysis_bundle(
    value: object,
) -> tuple[dict[str, object], dict[str, object]]:
    """Rebuild the admission and exact arrays consumed by the analyzer."""
    _require(
        type(value) is dict and value.get("record_type") == RECORD_TYPE,
        "Q019 raw analysis requires the hardened registered admission record",
    )
    try:
        rebuilt, bundle = derive_case_bundle(
            case_id=str(value.get("case_id", "")),
            artifact_root=Path(str(value.get("artifact_root", ""))),
            q043_qualification_path=Path(
                str(value.get("q043_dependency", {}).get("path", ""))
            ),
            q043_artifact_root=Path(str(value.get("q043_artifact_root", ""))),
            q023_qualification_path=Path(
                str(value.get("q023_dependency", {}).get("path", ""))
            ),
        )
    except RegisteredAdmissionError:
        raise
    except (OSError, RuntimeError, ValueError) as error:
        raise RegisteredAdmissionError(
            "Q019 registered admission immutable path anchor is invalid"
        ) from error
    _require(
        _strict_equal(value, rebuilt),
        "Q019 registered admission derived fields or authority drifted",
    )
    return rebuilt, bundle


def validate_bound_analysis_inputs(
    value: object,
    *,
    snapshots: object,
    particle_states: object,
    completion_record: object,
) -> dict[str, object]:
    """Require analyzer inputs to be the exact rederived raw-analysis bundle."""
    rebuilt, bundle = validate_analysis_bundle(value)
    _require(
        _analysis_equal(snapshots, bundle["snapshots"])
        and _analysis_equal(particle_states, bundle["particle_states"])
        and _analysis_equal(completion_record, bundle["completion_record"]),
        "Q019 analyzer inputs differ from the sealed raw-artifact reduction",
    )
    return rebuilt


__all__ = [
    "AUTHORIZATION_BOUNDARY",
    "RECORD_TYPE",
    "RegisteredAdmissionError",
    "canonical_sha256",
    "derive_case_bundle",
    "validate_admission",
    "validate_analysis_bundle",
    "validate_bound_analysis_inputs",
]
