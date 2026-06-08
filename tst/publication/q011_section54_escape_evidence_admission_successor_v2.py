#!/usr/bin/env python3
"""Bind Q011 production admissions to retained raw particle-escape evidence.

This non-authorizing successor preserves the immutable production-science v1
admission.  It verifies a complete sealed raw tree, independently reduces every
Q011 escape-event stream, and requires that reduction to equal the terminal
restart ledger.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Mapping, Sequence

if __package__:
    from . import q011_escape_event_evidence as escape
    from . import (
        q011_section54_production_science_admission_orchestration_successor_v1
        as production_v1,
    )
else:
    import q011_escape_event_evidence as escape
    import q011_section54_production_science_admission_orchestration_successor_v1 as production_v1


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q011_section54_escape_evidence_admission_successor_v2"
RECORD_TYPE = "q011_section54_escape_evidence_attempt_admission_successor_v2"
ESCAPE_PATH = re.compile(
    r"(?:[^/]+/)*[^/]+\.q011_escape_events\.cycle[0-9]{8}\.rank[0-9]{8}\.bin"
)
SHA256 = re.compile(r"[0-9a-f]{64}")
INTEGER = re.compile(r"-?[0-9]+")
REAL = re.compile(r"-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
LEDGER_FIELDS = (
    "ps_escaped_injected_cr_count_global",
    "ps_escaped_injected_cr_mass_global",
    "ps_escaped_injected_cr_momentum_x1_global",
    "ps_escaped_injected_cr_momentum_x2_global",
    "ps_escaped_injected_cr_momentum_x3_global",
    "ps_escaped_injected_cr_energy_global",
    "ps_escaped_initial_cr_count_global",
    "ps_escaped_injected_cr_term_count_global",
    "ps_escaped_injected_cr_abs_mass_global",
    "ps_escaped_injected_cr_abs_momentum_x1_global",
    "ps_escaped_injected_cr_abs_momentum_x2_global",
    "ps_escaped_injected_cr_abs_momentum_x3_global",
    "ps_escaped_injected_cr_abs_energy_global",
)


class EscapeEvidenceAdmissionError(ValueError):
    """Reject incomplete, mutable, or inconsistent retained Q011 evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EscapeEvidenceAdmissionError(message)


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _relative_path(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected path text")
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and path.as_posix() == value
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe relative path",
    )
    return value


def _sha256(value: object, *, label: str) -> str:
    _require(
        type(value) is str and SHA256.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _stable_read(path: Path, *, root: Path, expected: Mapping[str, object]) -> bytes:
    try:
        relative = path.relative_to(root).as_posix()
    except ValueError as error:
        raise EscapeEvidenceAdmissionError(
            f"{path}: retained artifact leaves raw root"
        ) from error
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222,
            f"{relative}: retained artifact is not one sealed regular file",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        stable = (
            "st_dev",
            "st_ino",
            "st_mode",
            "st_nlink",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
        _require(
            all(getattr(before, name) == getattr(after, name) for name in stable)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size
            and expected["path"] == relative
            and expected["byte_count"] == len(payload)
            and expected["sha256"] == hashlib.sha256(payload).hexdigest(),
            f"{relative}: retained artifact identity or digest drifted",
        )
        return payload
    finally:
        os.close(descriptor)


def _scan_sealed_tree(root: Path) -> set[str]:
    _require(root.is_absolute() and root.is_dir(), "raw root is unavailable")
    _require(
        not root.stat(follow_symlinks=False).st_mode & 0o222,
        "raw root is not sealed",
    )
    observed: set[str] = set()
    for directory, names, files in os.walk(root, followlinks=False):
        directory_path = Path(directory)
        metadata = directory_path.stat(follow_symlinks=False)
        _require(
            stat.S_ISDIR(metadata.st_mode) and not metadata.st_mode & 0o222,
            f"{directory_path}: retained directory is mutable",
        )
        for name in names:
            child = directory_path / name
            _require(
                not child.is_symlink() and child.is_dir(),
                f"{child}: retained tree contains a non-directory member",
            )
        for name in files:
            child = directory_path / name
            _require(not child.is_symlink(), f"{child}: retained tree contains a symlink")
            observed.add(child.relative_to(root).as_posix())
    return observed


def _inventory(
    value: object, *, root: Path
) -> tuple[list[dict[str, object]], dict[str, bytes]]:
    _require(type(value) is list and bool(value), "retained_inventory: expected array")
    normalized = []
    for index, raw in enumerate(value):
        label = f"retained_inventory[{index}]"
        _require(
            type(raw) is dict and set(raw) == {"path", "sha256", "byte_count"},
            f"{label}: schema drifted",
        )
        record = {
            "path": _relative_path(raw["path"], label=f"{label}/path"),
            "sha256": _sha256(raw["sha256"], label=f"{label}/sha256"),
            "byte_count": raw["byte_count"],
        }
        _require(
            type(record["byte_count"]) is int and record["byte_count"] > 0,
            f"{label}/byte_count: expected positive integer",
        )
        normalized.append(record)
    normalized.sort(key=lambda record: str(record["path"]))
    paths = [str(record["path"]) for record in normalized]
    _require(len(paths) == len(set(paths)), "retained_inventory: duplicate path")
    _require(
        set(paths) == _scan_sealed_tree(root),
        "retained_inventory does not exactly cover sealed raw tree",
    )
    payloads = {
        str(record["path"]): _stable_read(
            root / str(record["path"]), root=root, expected=record
        )
        for record in normalized
    }
    return normalized, payloads


def _parameter_blocks(payload: bytes, *, label: str) -> dict[str, dict[str, str]]:
    marker = payload.find(b"<par_end>")
    _require(marker > 0, f"{label}: restart header lacks <par_end>")
    try:
        text = payload[:marker].decode("ascii")
    except UnicodeDecodeError as error:
        raise EscapeEvidenceAdmissionError(
            f"{label}: restart header is not ASCII"
        ) from error
    blocks: dict[str, dict[str, str]] = {}
    active: str | None = None
    for line_number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            _require(
                line.endswith(">") and line.count("<") == line.count(">") == 1,
                f"{label}:{line_number}: malformed restart block",
            )
            active = line[1:-1]
            _require(active not in blocks, f"{label}: duplicate restart block {active}")
            blocks[active] = {}
            continue
        _require(
            active is not None and line.count("=") == 1,
            f"{label}:{line_number}: malformed restart parameter",
        )
        key, parameter = (item.strip() for item in line.split("=", 1))
        _require(
            bool(key) and key not in blocks[active],
            f"{label}: duplicate or empty restart parameter",
        )
        blocks[active][key] = parameter
    return blocks


def _terminal_restart(
    admission: Mapping[str, object], payloads: Mapping[str, bytes]
) -> dict[str, object]:
    terminal = admission["snapshot_bindings"][-1]
    bindings = terminal["artifact_bindings"]["rst"]
    paths = sorted(
        str(binding["path"])
        for binding in bindings
        if str(binding["path"]).endswith(".rst")
    )
    _require(bool(paths), "terminal snapshot lacks restart payload")
    parsed = []
    for path in paths:
        _require(path in payloads, f"terminal restart is absent from retained inventory: {path}")
        blocks = _parameter_blocks(payloads[path], label=path)
        _require(
            "time" in blocks and "problem" in blocks,
            f"{path}: restart lacks time/problem blocks",
        )
        time_block = blocks["time"]
        problem = blocks["problem"]
        _require("cycle" in time_block, f"{path}: restart lacks time/cycle")
        _require(
            all(field in problem for field in LEDGER_FIELDS)
            and all(
                field in problem
                for field in (
                    "ps_escape_ledger_schema",
                    "ps_escape_ledger_complete",
                    "ps_escape_audit_calls",
                    "ps_escape_last_audit_time",
                )
            ),
            f"{path}: restart lacks complete escape ledger",
        )
        _require(
            INTEGER.fullmatch(time_block["cycle"]) is not None
            and INTEGER.fullmatch(problem["ps_escape_ledger_schema"]) is not None
            and INTEGER.fullmatch(problem["ps_escape_audit_calls"]) is not None,
            f"{path}: restart integer escape state is malformed",
        )
        _require(
            problem["ps_escape_ledger_complete"] in {"1", "true"},
            f"{path}: restart escape ledger is incomplete",
        )
        numeric = {
            field: problem[field]
            for field in (*LEDGER_FIELDS, "ps_escape_last_audit_time")
        }
        _require(
            all(REAL.fullmatch(value) is not None for value in numeric.values()),
            f"{path}: restart escape ledger contains malformed real values",
        )
        parsed.append(
            {
                "cycle": int(time_block["cycle"]),
                "schema": int(problem["ps_escape_ledger_schema"]),
                "audit_calls": int(problem["ps_escape_audit_calls"]),
                "last_audit_time": float(problem["ps_escape_last_audit_time"]),
                "ledger": {field: float(problem[field]) for field in LEDGER_FIELDS},
            }
        )
    reference = parsed[0]
    _require(
        all(item == reference for item in parsed[1:]),
        "terminal restart ranks disagree on escape ledger",
    )
    terminal_cycle = terminal["cycle"]
    terminal_time = terminal["observed_committed_time"]
    _require(reference["schema"] == 2, "terminal restart escape ledger schema drifted")
    _require(
        reference["cycle"] == terminal_cycle
        and reference["audit_calls"] == 2 * terminal_cycle
        and math.isclose(
            reference["last_audit_time"],
            terminal_time,
            rel_tol=2.0e-13,
            abs_tol=2.0e-13,
        ),
        "terminal restart escape chronology differs from admitted VL2 chronology",
    )
    return {
        "paths": paths,
        "cycle": terminal_cycle,
        "observed_committed_time": terminal_time,
        "escape_audit_calls": reference["audit_calls"],
        "escape_last_audit_time": reference["last_audit_time"],
        "ledger": reference["ledger"],
    }


def _authorization() -> dict[str, bool]:
    return {
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "qualification_authorized": False,
        "claim_closure_authorized": False,
        "publication_authorized": False,
    }


def build_escape_evidence_admission(
    *,
    attempt_admission: Mapping[str, object],
    retained_inventory: Sequence[Mapping[str, object]],
    source_root: str | Path,
) -> dict[str, object]:
    """Build one non-authorizing, independently recomputable evidence admission."""
    admission = production_v1.validate_attempt_admission(
        attempt_admission, source_root=source_root
    )
    root = Path(str(admission["execution_binding"]["raw_output_root"]))
    inventory, payloads = _inventory(retained_inventory, root=root)
    retained_by_path = {str(record["path"]): record for record in inventory}
    for artifact in admission["raw_artifacts"]:
        path = str(artifact["path"])
        _require(
            path in retained_by_path
            and retained_by_path[path]
            == {
                "path": path,
                "sha256": artifact["sha256"],
                "byte_count": artifact["byte_count"],
            },
            f"v1 admitted raw artifact is absent or changed: {path}",
        )
    escape_paths = sorted(path for path in payloads if ESCAPE_PATH.fullmatch(path))
    _require(bool(escape_paths), "sealed raw tree contains no Q011 escape streams")
    _require(
        all(payloads[path].startswith(escape._BINARY_HEADER_MAGIC) for path in escape_paths),
        "production Q011 escape streams must use the compact binary schema",
    )
    streams = [
        escape.parse_stream_bytes(payloads[path], path=path) for path in escape_paths
    ]
    identity = admission["attempt_identity"]
    _require(
        all(stream.header["basename"] == identity["attempt_id"] for stream in streams),
        "escape stream basename differs from admitted attempt",
    )
    fingerprints = {stream.header["control_fingerprint"] for stream in streams}
    _require(len(fingerprints) == 1, "escape streams mix control fingerprints")
    reduction = escape.reduce_streams(streams)
    terminal = _terminal_restart(admission, payloads)
    escape.assert_matches_restart_ledger(reduction, terminal["ledger"])
    stream_bindings = [
        {
            **retained_by_path[path],
            "segment_start_cycle": stream.header["segment_start_cycle"],
            "segment_start_time": stream.header["segment_start_time"],
            "rank": stream.header["rank"],
            "nranks": stream.header["nranks"],
            "event_count": len(stream.events),
        }
        for path, stream in zip(escape_paths, streams)
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": admission["campaign_id"],
        "claim_id": admission["claim_id"],
        "status": "admitted_raw_escape_evidence_non_authorizing",
        "attempt_identity": identity,
        "production_v1_attempt_admission_sha256": _canonical_sha256(admission),
        "retained_inventory_sha256": _canonical_sha256(inventory),
        "retained_file_count": len(inventory),
        "escape_stream_bindings": stream_bindings,
        "escape_stream_bindings_sha256": _canonical_sha256(stream_bindings),
        "control_fingerprint": next(iter(fingerprints)),
        "independent_reduction": reduction,
        "independent_reduction_sha256": _canonical_sha256(reduction),
        "terminal_restart": terminal,
        "authorization": _authorization(),
    }


def validate_escape_evidence_admission(
    value: object,
    *,
    attempt_admission: Mapping[str, object],
    retained_inventory: Sequence[Mapping[str, object]],
    source_root: str | Path,
) -> dict[str, object]:
    _require(type(value) is dict, "escape evidence admission must be an object")
    rebuilt = build_escape_evidence_admission(
        attempt_admission=attempt_admission,
        retained_inventory=retained_inventory,
        source_root=source_root,
    )
    _require(value == rebuilt, "escape evidence admission derived fields drifted")
    return rebuilt


__all__ = [
    "EscapeEvidenceAdmissionError",
    "RECORD_TYPE",
    "SCHEMA_VERSION",
    "SUCCESSOR_ID",
    "build_escape_evidence_admission",
    "validate_escape_evidence_admission",
]
