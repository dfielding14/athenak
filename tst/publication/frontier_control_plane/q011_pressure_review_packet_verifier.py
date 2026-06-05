#!/usr/bin/env python3
"""Read-only installed consumer verification for a published Q011 review packet.

This standalone verifier checks immutable packet and aggregate receipt bindings,
the complete retained evidence tuple, sealed authoritative-reanalysis and human
reviewer attestations, and their frozen qualifying-source binding.
"""

from __future__ import annotations

import contextlib
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Iterator


WATERMARK = "ENGINEERING CALIBRATION ONLY - NOT SUN & BAI REPRODUCTION"
QUALIFICATION_EFFECT = "none_no_selection_no_execution_authorization_no_science_claim"
AGGREGATE_EVIDENCE_CLASS = "engineering_calibration_only"
AGGREGATE_QUALIFICATION_EFFECT = "none_no_sun_bai_claim_no_execution_authorization"
CONSUMPTION_RULE = (
    "receipt_plus_inode_bound_success_seal_are_required_acceptance_markers_"
    "artifacts_without_both_verified_markers_are_unpublished"
)
PACKET_RECEIPT_RECORD_TYPE = "q011_section54_pressure_pilot_review_packet_receipt"
AGGREGATE_RECEIPT_RECORD_TYPE = "q011_section54_pressure_pilot_bundle_publication_receipt"
SUCCESS_SEAL_RECORD_TYPE = "q011_receipt_inode_bound_publication_success_seal"
INVENTORY_RECORD_TYPE = "q011_section54_pressure_pilot_review_packet_inventory"
REVIEW_METRICS_RECORD_TYPE = "q011_section54_pressure_pilot_review_metrics"
INVENTORY_NAME = "packet_inventory.json"
AUTHORIZED_PRODUCTION_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING = {
    "path": (
        "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
        "q011_section54_pressure_pilot_review_packet_receipt.json"
    ),
    "sha256": "3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d",
}
AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING = {
    "path": (
        "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
        "q011_section54_pressure_pilot_bundle_receipt.json"
    ),
    "sha256": "9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9",
}
AGGREGATE_MANIFEST_NAME = "pressure_pilot_manifest.json"
AGGREGATE_MANIFEST_RECORD_TYPE = "q011_section54_pressure_pilot_bundle_manifest"
RAW_CASE_RECORD_TYPE = "q011_section54_pressure_pilot_verified_raw_case"
RAW_CASE_INVENTORY_NAME = "artifact_inventory.json"
RAW_CASE_IDS = ("ps_p0_1p00", "ps_p0_0p05", "ps_p0_0p10", "ps_p0_0p20")
RAW_CASE_PRESSURES = (1.0, 0.05, 0.1, 0.2)
RAW_CASE_ARGV_VALUES = ("1.0", "0.05", "0.10", "0.20")
RAW_CASE_TIMES = (0.0, 15.0, 30.0, 45.0, 60.0)
RAW_CASE_DESCRIPTOR_PATH = "analysis/analysis.json"
AUTHORIZED_ACTIVE_DECK_BINDING = {
    "path": "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
    "sha256": "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1",
}
AUTHORIZED_COMMON_OVERRIDES = (
    "mesh/nx1=100",
    "mesh/x1max=1200",
    "mesh/nx2=20",
    "mesh/x2max=240",
    "mesh_refinement/refinement=none",
    "mesh_refinement/num_levels=1",
    "time/tlim=60",
    "time/nlim=4096",
    "time/ndiag=50",
    "problem/ps_enable_curvature_amr=false",
    "problem/ps_feedback_diag_dcycle=50",
    "output1/variable=mhd_w_bcc",
    "output1/id=mhd_w_bcc",
    "output1/dt=15",
    "output2/dt=15",
    "output3/dt=15",
    "output4/dt=15",
    "output5/dt=15",
    "output6/dt=15",
)
AUTHORIZED_RUNTIME_PROFILE = "frontier_minimum_supported"
AUTHORIZED_PARALLEL_RANKS = 1
AUTHORIZED_ROCR_VISIBLE_DEVICE = 0
PACKET_MEMBERS = (
    "PRESSURE_REVIEW_PACKET.md",
    "figures/terminal_mhd_pic_pressure_comparison.png",
    "figures/terminal_profile_overlays.png",
    "pressure_review_metrics.json",
)
AUTHORIZED_PRODUCTION_SOURCE_BINDINGS = {
    "historical_v2_execution_preregistration": {
        "path": (
            "tst/publication/readiness/"
            "q011_section54_pressure_pilot_registered_execution_retry_successor_v2_2026-06-02.json"
        ),
        "sha256": "fcee4f99009ff10efa1e2f2d18f2a8bb4af551c64a97aa5c6f10a2dacf95f3bc",
    },
    "postrun_aggregate_source_authorization": {
        "path": (
            "tst/publication/readiness/"
            "q011_section54_pressure_pilot_postrun_aggregate_source_authorization_successor_v5_2026-06-05.json"
        ),
        "sha256": "93c2b9174d1546b883ac4910afe3f7f021ed4324f2831f80b19aee9395e5bd3e",
    },
    "registered_execution_preregistration": {
        "path": (
            "tst/publication/readiness/"
            "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json"
        ),
        "sha256": "ce0216b6852b932beae695ab53c1376133a921930c001226bf5efdb70a1c463a",
    },
    "reviewed_source_closure": [
        {
            "path": "tst/publication/frontier_q011_section54_pressure_pilot_publish_job.sh",
            "role": "aggregate_worker_wrapper",
            "sha256": "049293d07db11353b9862a2d0ba1347c5e5ea0f3985a3e30fb61d0da2793c359",
        },
        {
            "path": "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
            "role": "aggregate_publisher",
            "sha256": "56d54edef24249d03b6c1c0913cc64dd829999109e5d516fffc78f17c5328637",
        },
        {
            "path": "tst/publication/analyze_q011_section54_pressure_pilot.py",
            "role": "aggregate_analyzer",
            "sha256": "da745b3e52a86127455fd9f70d0341e1f3563b4344979cd4266df18d4695a316",
        },
        {
            "path": "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
            "role": "aggregate_case_verifier",
            "sha256": "8596d95c9b8952dcb760b10fbe00bf7ba8a2c895713cc3d36dd1efa1aac11ecc",
        },
        {
            "path": "tst/publication/analyze_q011_section54_outputs.py",
            "role": "output_primitives",
            "sha256": "20794ed45b80843b8e75accffc2e9336fcb7da58716fb818e466e193cd346356",
        },
        {
            "path": "tst/publication/frontier_f1_structured_artifacts.py",
            "role": "structured_artifact_helper",
            "sha256": "cf090115bcdfd143f67b12339115102b57e74cf3205b1ebb1521c144c3415a5a",
        },
        {
            "path": "tst/publication/pvtk_particles.py",
            "role": "particle_vtk_reader",
            "sha256": "187254c4ed10ce20ec383e710dadfa2e03e9f45317cce53e5d644c1db2be7339",
        },
        {
            "path": "tst/publication/render_q011_section54_pressure_pilot_review_packet.py",
            "role": "review_packet_renderer",
            "sha256": "5838b83733a991ffc00ee6464ce648d91fccc8e2819a44ff6726bd778d98277d",
        },
        {
            "path": "tst/publication/readiness/plotting_environment_lock_candidate_2026-05-30.json",
            "role": "plotting_environment_lock",
            "sha256": "cfb2dda28e19609193cee99913766c34eaf032f2f279d534d991af5382c96c16",
        },
        {
            "path": "tst/publication/frontier_q011_section54_pressure_pilot_review_packet_job.sh",
            "role": "review_packet_worker_wrapper",
            "sha256": "f13577a6f5f1dc557f40e2403fa210eaa0f6571d666835f36fda817987c708d1",
        },
        {
            "path": (
                "tst/publication/readiness/"
                "q011_section54_pressure_pilot_snapshot_time_compatibility_successor_2026-06-04.json"
            ),
            "role": "snapshot_time_compatibility_successor",
            "sha256": "2e608404604d85963f8371876ce0dc3243e5d0b1398820c5724c94ffc631c287",
        },
    ],
    "runtime_source_archive": {
        "archive_sha256": "d794beac1ecb169104ca60eca5284f77212464005461cba9e4d38436f668a220",
        "execution_mode": "worker_extracted_git_archive_head_verified",
        "git_commit": "ce41d4b29bc646b4f0740468e1026f7308b29026",
        "verified_source_closure_sha256": (
            "3dd9bf99b0bf67e515390860f3625d262dbef33dbbdbbb1ae29440a754fa955c"
        ),
    },
}

_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT_RE = re.compile(r"[0-9a-f]{40}")
_FRONTIER_HOST_RE = re.compile(r"frontier[0-9]{5}")
_STABLE_STAT_FIELDS = (
    "st_dev",
    "st_ino",
    "st_mode",
    "st_nlink",
    "st_size",
    "st_mtime_ns",
    "st_ctime_ns",
)
_IDENTITY_STAT_FIELDS = ("st_dev", "st_ino", "st_mode")
MAX_RETAINED_FILE_BYTES = 128 * 1024 * 1024
MAX_JSON_BYTES = 8 * 1024 * 1024
MAX_DIRECTORY_ENTRIES = 256
MAX_TREE_ENTRIES = 1024
MAX_TREE_DEPTH = 16
PRESSURE_GATE_REVIEW_NOT_BEFORE_UTC = "2026-06-05T09:09:42Z"
PRESSURE_GATE_ATTESTATION_ROOT_NAME = "pressure_gate_attestations"
PRESSURE_GATE_ATTESTATION_FILENAME = "attestation.json"
PRESSURE_GATE_CAPTURE_TO_SEAL_MAX_SECONDS = 15 * 60
PRESSURE_REANALYSIS_RECORD_TYPE = (
    "q011_section54_pressure_authoritative_reanalysis_attestation"
)
PRESSURE_REANALYSIS_QUALIFICATION_EFFECT = (
    "reanalysis_only_no_selection_no_execution_authorization_no_science_claim"
)
PRESSURE_REANALYSIS_EXECUTION_MODE = "clean_candidate_archive_reanalysis"
PRESSURE_REANALYSIS_OPERATOR_STATEMENT = (
    "Independently recomputed the exact immutable Q011 Section 5.4 pressure-pilot "
    "publication from the bound clean-candidate source archive; this attestation "
    "does not select pressure, authorize execution, or close a science claim."
)
PRESSURE_REVIEWER_RECORD_TYPE = "q011_section54_pressure_selection_reviewer_attestation"
PRESSURE_REVIEWER_QUALIFICATION_EFFECT = (
    "human_selection_review_only_no_execution_authorization_no_science_claim"
)
PRESSURE_REVIEWER_STATEMENT = (
    "Reviewed the immutable pressure-pilot packet and authoritative reanalysis, "
    "then selected exactly one registered pressure case; this attestation does "
    "not authorize execution or close a science claim."
)
PRESSURE_REANALYSIS_SOURCE_PATHS = (
    "tst/publication/q011_section54_pressure_selection.py",
    "tst/publication/q011_section54_historical_pressure_pilot_consumer.py",
    "tst/publication/frontier_control_plane/q011_pressure_review_packet_verifier.py",
    "tst/publication/analyze_q011_section54_outputs.py",
    "tst/publication/analyze_q011_section54_pressure_pilot.py",
    "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
    "tst/publication/frontier_f1_structured_artifacts.py",
    "tst/publication/pvtk_particles.py",
)
AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION = {
    "path": (
        "tst/publication/readiness/"
        "q011_section54_pressure_pilot_postrun_aggregate_source_authorization_"
        "successor_v5_2026-06-05.json"
    ),
    "sha256": "93c2b9174d1546b883ac4910afe3f7f021ed4324f2831f80b19aee9395e5bd3e",
}
_CANONICAL_UTC_RE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z")
_PRESSURE_GATE_ID_RE = re.compile(r"[a-z][a-z0-9._-]{2,63}")
_REANALYSIS_DIRECTORY_RE = re.compile(
    r"(?P<sealed>[0-9]{8}T[0-9]{6}Z)-q011-section54-pressure-reanalysis-"
    r"(?P<operator>[a-z][a-z0-9._-]{2,63})"
)
_REVIEWER_DIRECTORY_RE = re.compile(
    r"(?P<sealed>[0-9]{8}T[0-9]{6}Z)-q011-section54-pressure-selection-"
    r"(?P<reviewer>[a-z][a-z0-9._-]{2,63})"
)


class PressureReviewPacketVerificationError(ValueError):
    """Raised when a pressure review packet cannot be accepted fail-closed."""


def _fail(message: str) -> None:
    raise PressureReviewPacketVerificationError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _reject_nonfinite(value: str) -> None:
    raise ValueError(f"non-finite JSON value {value!r}")


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _decode_canonical_object(payload: bytes, label: str) -> dict[str, Any]:
    _require(len(payload) <= MAX_JSON_BYTES, f"{label} exceeds the JSON size limit")
    try:
        text = payload.decode("utf-8")
        value = json.loads(
            text,
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError, ValueError) as exc:
        _fail(f"{label} is not valid strict JSON: {exc}")
    _require(type(value) is dict, f"{label} must be a JSON object")
    try:
        canonical = _canonical_json_bytes(value)
    except (RecursionError, TypeError, ValueError) as exc:
        _fail(f"{label} cannot be canonically encoded: {exc}")
    _require(payload == canonical, f"{label} is not canonical JSON")
    return value


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _require_exact_keys(value: object, keys: set[str], label: str) -> None:
    _require(type(value) is dict, f"{label} must be a JSON object")
    actual = set(value)
    _require(actual == keys, f"{label} has unexpected keys: {sorted(actual ^ keys)}")


def _require_int(value: object, expected: int, label: str) -> None:
    _require(type(value) is int and value == expected, f"{label} must equal {expected}")


def _require_string(value: object, expected: str, label: str) -> None:
    _require(type(value) is str and value == expected, f"{label} must equal {expected!r}")


def _require_nonempty_string(value: object, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label} must be a nonempty string")
    return value


def _require_sha256(value: object, label: str) -> str:
    _require(
        type(value) is str and _SHA256_RE.fullmatch(value) is not None,
        f"{label} must be a lowercase SHA-256",
    )
    return value


def _require_nonnegative_number(value: object, label: str) -> float | int:
    _require(
        type(value) in {float, int}
        and math.isfinite(value)
        and value >= 0,
        f"{label} must be a finite nonnegative JSON number",
    )
    return value


def _stable_stat(value: os.stat_result) -> tuple[int, ...]:
    return tuple(getattr(value, field) for field in _STABLE_STAT_FIELDS)


def _identity_stat(value: os.stat_result) -> tuple[int, ...]:
    return tuple(getattr(value, field) for field in _IDENTITY_STAT_FIELDS)


def _same_json_value(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            _same_json_value(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _same_json_value(a, b) for a, b in zip(left, right)
        )
    return left == right


def _canonical_utc(value: object, label: str) -> tuple[str, datetime]:
    text = _require_nonempty_string(value, label)
    _require(
        _CANONICAL_UTC_RE.fullmatch(text) is not None,
        f"{label} must be a canonical whole-second UTC timestamp",
    )
    try:
        parsed = datetime.strptime(text, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc
        )
    except ValueError as exc:
        _fail(f"{label} must be a valid canonical UTC timestamp: {exc}")
    return text, parsed


def _current_utc(now: datetime | None) -> datetime:
    current = now or datetime.now(timezone.utc)
    _require(
        isinstance(current, datetime) and current.tzinfo is not None,
        "current time must be timezone-aware",
    )
    return current.astimezone(timezone.utc)


def _pressure_gate_id(value: object, label: str) -> str:
    text = _require_nonempty_string(value, label)
    _require(
        _PRESSURE_GATE_ID_RE.fullmatch(text) is not None,
        f"{label} must be a canonical lowercase reviewer/operator ID",
    )
    return text


def _canonical_rationale(value: object, label: str) -> str:
    text = _require_nonempty_string(value, label)
    _require(
        text == text.strip() and 16 <= len(text) <= 1024,
        f"{label} must be trimmed and contain 16-1024 characters",
    )
    _require(
        all(0x20 <= ord(character) <= 0x7E for character in text),
        f"{label} must be single-line printable ASCII",
    )
    return text


def _source_closure_digest(records: list[dict[str, str]]) -> str:
    payload = json.dumps(
        records,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    return _sha256(payload)


def _pressure_gate_binding(
    value: object,
    *,
    authorized_pic_root: Path,
    label: str,
) -> dict[str, str]:
    _require_exact_keys(value, {"path", "sha256"}, label)
    path_text = _require_nonempty_string(value["path"], f"{label}.path")
    path = _canonical_absolute_path(path_text, f"{label}.path")
    archive_root = authorized_pic_root / PRESSURE_GATE_ATTESTATION_ROOT_NAME
    _require(
        path.name == PRESSURE_GATE_ATTESTATION_FILENAME
        and path.parent.parent == archive_root,
        f"{label}.path is outside the fixed pressure-gate attestation layout",
    )
    return {
        "path": str(path),
        "sha256": _require_sha256(value["sha256"], f"{label}.sha256"),
    }


def _published_binding(value: object, authorized_pic_root: Path, label: str) -> dict[str, str]:
    _require_exact_keys(value, {"path", "sha256"}, label)
    path = _direct_publication_path(
        value["path"],
        authorized_pic_root / "publication",
        f"{label}.path",
    )
    return {
        "path": str(path),
        "sha256": _require_sha256(value["sha256"], f"{label}.sha256"),
    }


def _expected_reanalysis_result(
    *,
    packet_receipt_binding: dict[str, str],
    aggregate_receipt_binding: dict[str, str],
    pilot_bundle_manifest_sha256: str,
    aggregate_pilot_analysis_sha256: str,
) -> dict[str, str]:
    return {
        "packet_receipt_sha256": packet_receipt_binding["sha256"],
        "aggregate_receipt_sha256": aggregate_receipt_binding["sha256"],
        "manifest_sha256": pilot_bundle_manifest_sha256,
        "analysis_result_sha256": aggregate_pilot_analysis_sha256,
        "status": "pass_engineering_calibration_only",
    }


def _open_flags(*, directory: bool = False) -> int:
    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    flags |= getattr(os, "O_NONBLOCK", 0)
    if directory:
        flags |= getattr(os, "O_DIRECTORY", 0)
    return flags


def _require_readonly_file_metadata(metadata: os.stat_result, label: str) -> None:
    _require(stat.S_ISREG(metadata.st_mode), f"{label} must be a regular file")
    _require(metadata.st_nlink == 1, f"{label} must have exactly one hard link")
    _require(metadata.st_mode & 0o222 == 0, f"{label} must be read-only")


def _require_readonly_directory_metadata(metadata: os.stat_result, label: str) -> None:
    _require(stat.S_ISDIR(metadata.st_mode), f"{label} must be a directory")
    _require(metadata.st_mode & 0o222 == 0, f"{label} must be read-only")


def _require_entry_identity(
    parent_fd: int,
    name: str,
    expected: os.stat_result,
    label: str,
) -> None:
    try:
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except OSError as exc:
        _fail(f"{label} path cannot be re-read: {exc}")
    _require(
        _stable_stat(current) == _stable_stat(expected),
        f"{label} path identity or metadata changed during verification",
    )


class _RetainedFile:
    def __init__(self, parent_fd: int, name: str, label: str) -> None:
        self._parent_fd = parent_fd
        self._name = name
        self._label = label
        self._fd = -1
        self.metadata: os.stat_result

    def __enter__(self) -> "_RetainedFile":
        try:
            self._fd = os.open(self._name, _open_flags(), dir_fd=self._parent_fd)
            self.metadata = os.fstat(self._fd)
            _require_readonly_file_metadata(self.metadata, self._label)
            _require(
                self.metadata.st_size <= MAX_RETAINED_FILE_BYTES,
                f"{self._label} exceeds the retained-file size limit",
            )
            _require_entry_identity(
                self._parent_fd,
                self._name,
                self.metadata,
                self._label,
            )
        except PressureReviewPacketVerificationError:
            self.close()
            raise
        except OSError as exc:
            self.close()
            _fail(f"{self._label} cannot be opened without following links: {exc}")
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def close(self) -> None:
        if self._fd >= 0:
            os.close(self._fd)
            self._fd = -1

    def read(self) -> bytes:
        _require(self._fd >= 0, f"{self._label} retained descriptor is closed")
        before = os.fstat(self._fd)
        _require(
            _stable_stat(before) == _stable_stat(self.metadata),
            f"{self._label} changed before stable read",
        )
        try:
            os.lseek(self._fd, 0, os.SEEK_SET)
            chunks: list[bytes] = []
            total = 0
            while True:
                chunk = os.read(self._fd, 1024 * 1024)
                if not chunk:
                    break
                chunks.append(chunk)
                total += len(chunk)
                _require(
                    total <= MAX_RETAINED_FILE_BYTES,
                    f"{self._label} exceeds the retained-file size limit",
                )
        except OSError as exc:
            _fail(f"{self._label} cannot be read: {exc}")
        after = os.fstat(self._fd)
        _require(
            _stable_stat(after) == _stable_stat(self.metadata),
            f"{self._label} changed during stable read",
        )
        _require_entry_identity(self._parent_fd, self._name, self.metadata, self._label)
        payload = b"".join(chunks)
        _require(
            len(payload) == self.metadata.st_size,
            f"{self._label} stable read is incomplete",
        )
        return payload


def _canonical_absolute_path(value: str | os.PathLike[str], label: str) -> Path:
    path = Path(value)
    _require(path.is_absolute(), f"{label} must be absolute")
    lexical = Path(os.path.abspath(os.fspath(path)))
    _require(path == lexical, f"{label} must be lexically canonical")
    return path


def _direct_publication_path(value: object, publication_root: Path, label: str) -> Path:
    text = _require_nonempty_string(value, label)
    path = Path(text)
    _require(path.is_absolute(), f"{label} must be absolute")
    _require(path == Path(os.path.abspath(text)), f"{label} must be lexically canonical")
    _require(path.parent == publication_root, f"{label} must be a direct publication child")
    _require(path.name not in {"", ".", ".."}, f"{label} has an invalid basename")
    return path


def _direct_publication_argument_path(
    value: str | os.PathLike[str],
    publication_root: Path,
    label: str,
) -> Path:
    try:
        text = os.fspath(value)
    except TypeError as exc:
        _fail(f"{label} must be a filesystem path: {exc}")
    _require(type(text) is str and bool(text), f"{label} must be a nonempty path")
    path = Path(os.path.abspath(text))
    _require(path.parent == publication_root, f"{label} must be a direct publication child")
    _require(path.name not in {"", ".", ".."}, f"{label} has an invalid basename")
    return path


def _binding(value: object, publication_root: Path, label: str) -> dict[str, str]:
    _require_exact_keys(value, {"path", "sha256"}, label)
    path = _direct_publication_path(value["path"], publication_root, f"{label}.path")
    digest = _require_sha256(value["sha256"], f"{label}.sha256")
    return {"path": str(path), "sha256": digest}


def _aggregate_bundle_binding(
    value: object,
    publication_root: Path,
    label: str,
) -> dict[str, str]:
    _require_exact_keys(value, {"path", "manifest_sha256"}, label)
    path = _direct_publication_path(value["path"], publication_root, f"{label}.path")
    digest = _require_sha256(value["manifest_sha256"], f"{label}.manifest_sha256")
    return {"path": str(path), "manifest_sha256": digest}


def _publication_identity(metadata: os.stat_result) -> dict[str, int]:
    return {"device": metadata.st_dev, "inode": metadata.st_ino}


def _require_publication_identity(
    value: object,
    metadata: os.stat_result,
    label: str,
) -> None:
    _require_exact_keys(value, {"device", "inode"}, label)
    expected = _publication_identity(metadata)
    _require(
        type(value["device"]) is int
        and type(value["inode"]) is int
        and value == expected,
        f"{label} does not identify the retained publication root",
    )


def _guard_name(receipt_name: str) -> str:
    return f".{receipt_name}.publication-invalid"


def _seal_name(receipt_name: str) -> str:
    return f".{receipt_name}.publication-success"


def _require_guard_absent(publication_fd: int, receipt_name: str, label: str) -> None:
    try:
        os.stat(_guard_name(receipt_name), dir_fd=publication_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    except OSError as exc:
        _fail(f"{label} guard absence cannot be established: {exc}")
    _fail(f"{label} publication guard is present")


def _expected_seal(
    receipt_name: str,
    receipt_sha256: str,
    receipt_metadata: os.stat_result,
    publication_metadata: os.stat_result,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "record_type": SUCCESS_SEAL_RECORD_TYPE,
        "publication_root_identity": _publication_identity(publication_metadata),
        "receipt_name": receipt_name,
        "receipt_sha256": receipt_sha256,
        "receipt_identity": {
            "device": receipt_metadata.st_dev,
            "inode": receipt_metadata.st_ino,
        },
    }


def _open_and_verify_seal(
    stack: contextlib.ExitStack,
    acceptance_fd: int,
    receipt_name: str,
    receipt_sha256: str,
    receipt_metadata: os.stat_result,
    publication_metadata: os.stat_result,
    label: str,
) -> tuple[_RetainedFile, bytes]:
    retained = stack.enter_context(
        _RetainedFile(acceptance_fd, _seal_name(receipt_name), f"{label} success seal")
    )
    payload = retained.read()
    value = _decode_canonical_object(payload, f"{label} success seal")
    expected = _expected_seal(
        receipt_name,
        receipt_sha256,
        receipt_metadata,
        publication_metadata,
    )
    _require(
        _same_json_value(value, expected),
        f"{label} success seal does not bind the retained receipt",
    )
    return retained, payload


def _relative_source_path(value: object, label: str) -> str:
    text = _require_nonempty_string(value, label)
    path = PurePosixPath(text)
    _require(not path.is_absolute(), f"{label} must be relative")
    _require(path.as_posix() == text, f"{label} must be canonical POSIX syntax")
    _require(all(part not in {"", ".", ".."} for part in path.parts), f"{label} is unsafe")
    _require(len(path.parts) <= MAX_TREE_DEPTH, f"{label} exceeds the path-depth limit")
    return text


def _validate_source_receipt_binding(value: object, label: str) -> None:
    _require_exact_keys(value, {"path", "sha256"}, label)
    _relative_source_path(value["path"], f"{label}.path")
    _require_sha256(value["sha256"], f"{label}.sha256")


def _validate_source_bindings(
    value: object,
    label: str,
    *,
    require_verified_archive: bool,
) -> None:
    keys = {
        "postrun_aggregate_source_authorization",
        "registered_execution_preregistration",
        "historical_v2_execution_preregistration",
        "reviewed_source_closure",
        "runtime_source_archive",
    }
    _require_exact_keys(value, keys, label)
    for key in (
        "postrun_aggregate_source_authorization",
        "registered_execution_preregistration",
        "historical_v2_execution_preregistration",
    ):
        _validate_source_receipt_binding(value[key], f"{label}.{key}")

    closure = value["reviewed_source_closure"]
    _require(
        type(closure) is list and bool(closure),
        f"{label}.reviewed_source_closure must be a nonempty list",
    )
    seen_roles: set[str] = set()
    seen_paths: set[str] = set()
    for index, member in enumerate(closure):
        member_label = f"{label}.reviewed_source_closure[{index}]"
        _require_exact_keys(member, {"role", "path", "sha256"}, member_label)
        role = _require_nonempty_string(member["role"], f"{member_label}.role")
        path = _relative_source_path(member["path"], f"{member_label}.path")
        _require_sha256(member["sha256"], f"{member_label}.sha256")
        _require(role not in seen_roles, f"{label}.reviewed_source_closure has duplicate roles")
        _require(path not in seen_paths, f"{label}.reviewed_source_closure has duplicate paths")
        seen_roles.add(role)
        seen_paths.add(path)

    archive = value["runtime_source_archive"]
    _require_exact_keys(
        archive,
        {
            "execution_mode",
            "git_commit",
            "archive_sha256",
            "verified_source_closure_sha256",
        },
        f"{label}.runtime_source_archive",
    )
    mode = _require_nonempty_string(
        archive["execution_mode"],
        f"{label}.runtime_source_archive.execution_mode",
    )
    if mode == "direct_api_nonproduction_only":
        _require(
            not require_verified_archive,
            f"{label}.runtime_source_archive production binding must be worker verified",
        )
        for key in ("git_commit", "archive_sha256", "verified_source_closure_sha256"):
            _require(archive[key] is None, f"{label}.runtime_source_archive.{key} must be null")
    elif mode == "worker_extracted_git_archive_head_verified":
        commit = archive["git_commit"]
        _require(
            type(commit) is str and _GIT_COMMIT_RE.fullmatch(commit) is not None,
            f"{label}.runtime_source_archive.git_commit must be a lowercase Git commit",
        )
        _require_sha256(archive["archive_sha256"], f"{label}.runtime_source_archive.archive_sha256")
        _require_sha256(
            archive["verified_source_closure_sha256"],
            f"{label}.runtime_source_archive.verified_source_closure_sha256",
        )
        expected_closure_sha256 = _sha256(
            _canonical_json_bytes(
                {
                    "postrun_aggregate_source_authorization": value[
                        "postrun_aggregate_source_authorization"
                    ],
                    "reviewed_source_closure": closure,
                }
            )
        )
        _require(
            archive["verified_source_closure_sha256"] == expected_closure_sha256,
            f"{label}.runtime_source_archive verified source closure hash drifted",
        )
    else:
        _fail(f"{label}.runtime_source_archive.execution_mode is unsupported")
    if require_verified_archive:
        _require(
            _same_json_value(value, AUTHORIZED_PRODUCTION_SOURCE_BINDINGS),
            f"{label} does not equal the authorized immutable production source tuple",
        )


def _validate_raw_cases(value: object, authorized_pic_root: Path) -> None:
    _require(type(value) is list, "aggregate receipt.raw_cases must be a list")
    _require(
        [
            case.get("case_id")
            for case in value
            if type(case) is dict
        ]
        == list(RAW_CASE_IDS),
        "aggregate receipt.raw_cases must contain the exact ordered four-case set",
    )
    runs_root = authorized_pic_root / "runs"
    for index, case in enumerate(value):
        label = f"aggregate receipt.raw_cases[{index}]"
        _require_exact_keys(
            case,
            {
                "case_id",
                "artifact_dir",
                "descriptor_path",
                "descriptor_sha256",
                "artifact_inventory_sha256",
                "runtime_artifacts",
            },
            label,
        )
        _require_string(case["case_id"], RAW_CASE_IDS[index], f"{label}.case_id")
        artifact_dir_text = _require_nonempty_string(
            case["artifact_dir"],
            f"{label}.artifact_dir",
        )
        artifact_dir = Path(artifact_dir_text)
        _require(artifact_dir.is_absolute(), f"{label}.artifact_dir must be absolute")
        _require(
            artifact_dir == Path(os.path.abspath(artifact_dir_text))
            and artifact_dir != runs_root
            and runs_root in artifact_dir.parents,
            f"{label}.artifact_dir must be canonical below the authorized runs root",
        )
        _require_string(
            case["descriptor_path"],
            RAW_CASE_DESCRIPTOR_PATH,
            f"{label}.descriptor_path",
        )
        _require_sha256(case["descriptor_sha256"], f"{label}.descriptor_sha256")
        _require_sha256(
            case["artifact_inventory_sha256"],
            f"{label}.artifact_inventory_sha256",
        )
        runtime = case["runtime_artifacts"]
        _require(type(runtime) is dict, f"{label}.runtime_artifacts must be an object")
        runtime_paths = set(runtime)
        required_runtime = {
            "athena_stdout.txt",
            "athena_stdout.sha256",
            "athena_stderr.txt",
        }
        allowlists = [
            path
            for path in runtime_paths - required_runtime
            if type(path) is str and path.endswith(".environment.allowlist.txt")
        ]
        _require(
            len(runtime_paths) == 4
            and required_runtime <= runtime_paths
            and len(allowlists) == 1,
            f"{label}.runtime_artifacts has unexpected members",
        )
        for path, digest in runtime.items():
            _relative_source_path(path, f"{label}.runtime_artifacts path")
            _require_sha256(digest, f"{label}.runtime_artifacts[{path!r}]")


def _raw_case_relative_path(case: dict[str, Any], authorized_pic_root: Path, label: str) -> str:
    artifact_dir = Path(case["artifact_dir"])
    try:
        relative = artifact_dir.relative_to(authorized_pic_root / "runs")
    except ValueError:
        _fail(f"{label}.artifact_dir is outside the authorized runs root")
    return _relative_source_path(relative.as_posix(), f"{label}.artifact_dir")


def _validate_aggregate_receipt(
    value: dict[str, Any],
    publication_root: Path,
    publication_metadata: os.stat_result,
    authorized_pic_root: Path,
) -> tuple[dict[str, str], dict[str, str]]:
    _require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "evidence_class",
            "qualification_effect",
            "consumption_rule",
            "publication_root_identity",
            "aggregate_bundle",
            "aggregate_analysis",
            "source_bindings",
            "raw_cases",
        },
        "aggregate receipt",
    )
    _require_int(value["schema_version"], 1, "aggregate receipt.schema_version")
    _require_string(
        value["record_type"],
        AGGREGATE_RECEIPT_RECORD_TYPE,
        "aggregate receipt.record_type",
    )
    _require_string(
        value["evidence_class"],
        AGGREGATE_EVIDENCE_CLASS,
        "aggregate receipt.evidence_class",
    )
    _require_string(
        value["qualification_effect"],
        AGGREGATE_QUALIFICATION_EFFECT,
        "aggregate receipt.qualification_effect",
    )
    _require_string(
        value["consumption_rule"],
        CONSUMPTION_RULE,
        "aggregate receipt.consumption_rule",
    )
    _require_publication_identity(
        value["publication_root_identity"],
        publication_metadata,
        "aggregate receipt.publication_root_identity",
    )
    aggregate_bundle = _aggregate_bundle_binding(
        value["aggregate_bundle"],
        publication_root,
        "aggregate receipt.aggregate_bundle",
    )
    aggregate_analysis = _binding(
        value["aggregate_analysis"],
        publication_root,
        "aggregate receipt.aggregate_analysis",
    )
    _validate_source_bindings(
        value["source_bindings"],
        "aggregate receipt.source_bindings",
        require_verified_archive=authorized_pic_root == AUTHORIZED_PRODUCTION_PIC_ROOT,
    )
    _validate_raw_cases(value["raw_cases"], authorized_pic_root)
    return aggregate_bundle, aggregate_analysis


def _validate_packet_receipt(
    value: dict[str, Any],
    publication_root: Path,
    publication_metadata: os.stat_result,
    authorized_pic_root: Path,
) -> tuple[dict[str, str], Path, str]:
    _require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "watermark",
            "qualification_effect",
            "consumption_rule",
            "publication_root_identity",
            "aggregate_receipt",
            "packet_root",
            "inventory_sha256",
            "source_bindings",
        },
        "packet receipt",
    )
    _require_int(value["schema_version"], 1, "packet receipt.schema_version")
    _require_string(value["record_type"], PACKET_RECEIPT_RECORD_TYPE, "packet receipt.record_type")
    _require_string(value["watermark"], WATERMARK, "packet receipt.watermark")
    _require_string(
        value["qualification_effect"],
        QUALIFICATION_EFFECT,
        "packet receipt.qualification_effect",
    )
    _require_string(value["consumption_rule"], CONSUMPTION_RULE, "packet receipt.consumption_rule")
    _require_publication_identity(
        value["publication_root_identity"],
        publication_metadata,
        "packet receipt.publication_root_identity",
    )
    aggregate_binding = _binding(
        value["aggregate_receipt"],
        publication_root,
        "packet receipt.aggregate_receipt",
    )
    packet_root = _direct_publication_path(
        value["packet_root"],
        publication_root,
        "packet receipt.packet_root",
    )
    inventory_sha256 = _require_sha256(value["inventory_sha256"], "packet receipt.inventory_sha256")
    _validate_source_bindings(
        value["source_bindings"],
        "packet receipt.source_bindings",
        require_verified_archive=authorized_pic_root == AUTHORIZED_PRODUCTION_PIC_ROOT,
    )
    return aggregate_binding, packet_root, inventory_sha256


def _open_child_directory(
    parent_fd: int,
    name: str,
    observed: os.stat_result,
    label: str,
) -> int:
    fd = -1
    try:
        fd = os.open(name, _open_flags(directory=True), dir_fd=parent_fd)
        opened = os.fstat(fd)
    except OSError as exc:
        if fd >= 0:
            os.close(fd)
        _fail(f"{label} cannot be opened without following links: {exc}")
    if _stable_stat(opened) != _stable_stat(observed):
        os.close(fd)
        _fail(f"{label} changed while being opened")
    return fd


@contextlib.contextmanager
def _retained_child_directory(
    parent_fd: int,
    name: str,
    label: str,
) -> Iterator[tuple[int, os.stat_result]]:
    try:
        observed = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except OSError as exc:
        _fail(f"{label} cannot be inspected: {exc}")
    _require(stat.S_ISDIR(observed.st_mode), f"{label} must be a directory")
    fd = _open_child_directory(parent_fd, name, observed, label)
    try:
        yield fd, observed
        after = os.fstat(fd)
        _require(
            _stable_stat(after) == _stable_stat(observed),
            f"{label} changed during verification",
        )
        _require_entry_identity(parent_fd, name, observed, label)
    finally:
        os.close(fd)


@contextlib.contextmanager
def _retained_absolute_directory(
    value: str | os.PathLike[str],
    label: str,
) -> Iterator[tuple[Path, int, os.stat_result]]:
    path = _canonical_absolute_path(value, label)
    root_fd = -1
    try:
        root_fd = os.open("/", _open_flags(directory=True))
        root_metadata = os.fstat(root_fd)
        _require(stat.S_ISDIR(root_metadata.st_mode), "filesystem root must be a directory")
        with contextlib.ExitStack() as stack:
            current_fd = root_fd
            current_metadata = root_metadata
            traversed = Path("/")
            for part in path.parts[1:]:
                traversed /= part
                current_fd, current_metadata = stack.enter_context(
                    _retained_ancestry_child_directory(
                        current_fd,
                        part,
                        f"{label} component {str(traversed)!r}",
                    )
                )
            yield path, current_fd, current_metadata
            _require(
                _identity_stat(os.fstat(root_fd)) == _identity_stat(root_metadata),
                "filesystem root changed during verification",
            )
    except PressureReviewPacketVerificationError:
        raise
    except OSError as exc:
        _fail(f"{label} cannot be opened by descriptor-relative ancestry: {exc}")
    finally:
        if root_fd >= 0:
            os.close(root_fd)


@contextlib.contextmanager
def _retained_ancestry_child_directory(
    parent_fd: int,
    name: str,
    label: str,
) -> Iterator[tuple[int, os.stat_result]]:
    try:
        observed = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except OSError as exc:
        _fail(f"{label} cannot be inspected: {exc}")
    _require(stat.S_ISDIR(observed.st_mode), f"{label} must be a directory")
    fd = -1
    try:
        fd = os.open(name, _open_flags(directory=True), dir_fd=parent_fd)
        opened = os.fstat(fd)
        _require(
            _identity_stat(opened) == _identity_stat(observed),
            f"{label} changed while being opened",
        )
        yield fd, opened
        _require(
            _identity_stat(os.fstat(fd)) == _identity_stat(opened),
            f"{label} identity changed during verification",
        )
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        _require(
            _identity_stat(current) == _identity_stat(opened),
            f"{label} path identity changed during verification",
        )
    except PressureReviewPacketVerificationError:
        raise
    except OSError as exc:
        _fail(f"{label} cannot be retained without following links: {exc}")
    finally:
        if fd >= 0:
            os.close(fd)


@contextlib.contextmanager
def _retained_relative_directory(
    parent_fd: int,
    relative: str,
    label: str,
) -> Iterator[tuple[int, os.stat_result]]:
    normalized = _relative_source_path(relative, f"{label} relative path")
    with contextlib.ExitStack() as stack:
        current_fd = parent_fd
        current_metadata: os.stat_result | None = None
        traversed: list[str] = []
        for part in PurePosixPath(normalized).parts:
            traversed.append(part)
            current_fd, current_metadata = stack.enter_context(
                _retained_child_directory(
                    current_fd,
                    part,
                    f"{label} component {'/'.join(traversed)!r}",
                )
            )
        _require(current_metadata is not None, f"{label} relative path must not be empty")
        yield current_fd, current_metadata


def _enter_retained_relative_file(
    stack: contextlib.ExitStack,
    parent_fd: int,
    relative: str,
    label: str,
) -> _RetainedFile:
    normalized = _relative_source_path(relative, f"{label} relative path")
    parts = PurePosixPath(normalized).parts
    current_fd = parent_fd
    traversed: list[str] = []
    for part in parts[:-1]:
        traversed.append(part)
        current_fd, _metadata = stack.enter_context(
            _retained_child_directory(
                current_fd,
                part,
                f"{label} parent {'/'.join(traversed)!r}",
            )
        )
    return stack.enter_context(_RetainedFile(current_fd, parts[-1], label))


def _bounded_directory_names(
    directory_fd: int,
    prefix: str,
    tree_label: str,
    *,
    require_nonempty: bool,
) -> list[str]:
    names: list[str] = []
    try:
        with os.scandir(directory_fd) as entries:
            for entry in entries:
                _require(
                    len(names) < MAX_DIRECTORY_ENTRIES,
                    f"{tree_label} directory {prefix or '.'!r} exceeds the entry limit",
                )
                names.append(entry.name)
    except OSError as exc:
        _fail(f"{tree_label} directory {prefix or '.'!r} cannot be listed: {exc}")
    if require_nonempty:
        _require(bool(names), f"{tree_label} directory {prefix or '.'!r} must not be empty")
    return sorted(names)


def _scan_immutable_directory(
    directory_fd: int,
    prefix: str,
    tree_label: str,
    directories: dict[str, tuple[int, ...]],
    files: dict[str, tuple[int, ...]],
    payloads: dict[str, bytes],
    depth: int,
) -> None:
    _require(depth <= MAX_TREE_DEPTH, f"{tree_label} exceeds the tree-depth limit")
    names = _bounded_directory_names(
        directory_fd,
        prefix,
        tree_label,
        require_nonempty=True,
    )
    for name in names:
        _require(
            len(directories) + len(files) < MAX_TREE_ENTRIES,
            f"{tree_label} exceeds the tree-entry limit",
        )
        _require(
            name not in {"", ".", ".."} and "/" not in name,
            f"{tree_label} contains an unsafe entry name",
        )
        relative = f"{prefix}/{name}" if prefix else name
        try:
            observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            _fail(f"{tree_label} entry {relative!r} cannot be inspected: {exc}")
        if stat.S_ISREG(observed.st_mode):
            with _RetainedFile(directory_fd, name, f"{tree_label} member {relative!r}") as retained:
                payloads[relative] = retained.read()
                files[relative] = _stable_stat(retained.metadata)
        elif stat.S_ISDIR(observed.st_mode):
            _require_readonly_directory_metadata(observed, f"{tree_label} directory {relative!r}")
            child_fd = _open_child_directory(
                directory_fd,
                name,
                observed,
                f"{tree_label} directory {relative!r}",
            )
            try:
                directories[relative] = _stable_stat(observed)
                _scan_immutable_directory(
                    child_fd,
                    relative,
                    tree_label,
                    directories,
                    files,
                    payloads,
                    depth + 1,
                )
                after = os.fstat(child_fd)
                _require(
                    _stable_stat(after) == _stable_stat(observed),
                    f"{tree_label} directory {relative!r} changed during verification",
                )
                _require_entry_identity(
                    directory_fd,
                    name,
                    observed,
                    f"{tree_label} directory {relative!r}",
                )
            finally:
                os.close(child_fd)
        else:
            _fail(f"{tree_label} entry {relative!r} has an unsupported file type")


def _scan_immutable_tree(directory_fd: int, tree_label: str) -> dict[str, object]:
    directories: dict[str, tuple[int, ...]] = {}
    files: dict[str, tuple[int, ...]] = {}
    payloads: dict[str, bytes] = {}
    _scan_immutable_directory(
        directory_fd,
        "",
        tree_label,
        directories,
        files,
        payloads,
        0,
    )
    return {"directories": directories, "files": files, "payloads": payloads}


def _scan_packet_tree(packet_root_fd: int) -> dict[str, object]:
    return _scan_immutable_tree(packet_root_fd, "packet")


def _relative_binding(value: object, label: str) -> dict[str, str]:
    _require_exact_keys(value, {"path", "sha256"}, label)
    return {
        "path": _relative_source_path(value["path"], f"{label}.path"),
        "sha256": _require_sha256(value["sha256"], f"{label}.sha256"),
    }


def _add_expected_binding(
    expected: dict[str, str],
    value: object,
    label: str,
    *,
    expected_path: str | None = None,
) -> dict[str, str]:
    binding = _relative_binding(value, label)
    if expected_path is not None:
        _require_string(binding["path"], expected_path, f"{label}.path")
    _require(binding["path"] not in expected, f"{label}.path is duplicated")
    expected[binding["path"]] = binding["sha256"]
    return binding


def _validate_aggregate_manifest(
    payload: bytes,
    expected_manifest_sha256: str,
    source_bindings: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, str], list[set[str]]]:
    _require(
        _sha256(payload) == expected_manifest_sha256,
        "aggregate bundle manifest hash does not match the aggregate receipt",
    )
    manifest = _decode_canonical_object(payload, "aggregate bundle manifest")
    _require_exact_keys(
        manifest,
        {
            "schema_version",
            "record_type",
            "evidence_class",
            "qualification_effect",
            "active_deck_binding",
            "preregistration_binding",
            "registered_execution_preregistration_binding",
            "cases",
        },
        "aggregate bundle manifest",
    )
    _require_int(manifest["schema_version"], 1, "aggregate bundle manifest.schema_version")
    _require_string(
        manifest["record_type"],
        AGGREGATE_MANIFEST_RECORD_TYPE,
        "aggregate bundle manifest.record_type",
    )
    _require_string(
        manifest["evidence_class"],
        AGGREGATE_EVIDENCE_CLASS,
        "aggregate bundle manifest.evidence_class",
    )
    _require_string(
        manifest["qualification_effect"],
        AGGREGATE_QUALIFICATION_EFFECT,
        "aggregate bundle manifest.qualification_effect",
    )
    active_deck = _relative_binding(
        manifest["active_deck_binding"],
        "aggregate bundle manifest.active_deck_binding",
    )
    _require(
        _same_json_value(active_deck, AUTHORIZED_ACTIVE_DECK_BINDING),
        "aggregate bundle manifest active deck is not the authorized pressure-pilot deck",
    )
    preregistration = _relative_binding(
        manifest["preregistration_binding"],
        "aggregate bundle manifest.preregistration_binding",
    )
    registered = _relative_binding(
        manifest["registered_execution_preregistration_binding"],
        "aggregate bundle manifest.registered_execution_preregistration_binding",
    )
    _require(
        _same_json_value(
            preregistration,
            source_bindings["postrun_aggregate_source_authorization"],
        ),
        "aggregate bundle manifest preregistration binding differs from source bindings",
    )
    _require(
        _same_json_value(
            registered,
            source_bindings["registered_execution_preregistration"],
        ),
        "aggregate bundle manifest registered-execution binding differs from source bindings",
    )

    cases = manifest["cases"]
    _require(
        type(cases) is list and len(cases) == len(RAW_CASE_IDS),
        "aggregate bundle manifest must contain exactly four cases",
    )
    expected = {AGGREGATE_MANIFEST_NAME: expected_manifest_sha256}
    case_paths: list[set[str]] = []
    for case_index, (case, case_id, pressure, argv_value) in enumerate(
        zip(cases, RAW_CASE_IDS, RAW_CASE_PRESSURES, RAW_CASE_ARGV_VALUES)
    ):
        label = f"aggregate bundle manifest.cases[{case_index}]"
        _require_exact_keys(
            case,
            {"case_id", "ps_p0", "overrides", "snapshots", "stdout", "terminal_restart"},
            label,
        )
        _require_string(case["case_id"], case_id, f"{label}.case_id")
        _require(
            type(case["ps_p0"]) is float and case["ps_p0"] == pressure,
            f"{label}.ps_p0 differs from the registered pressure",
        )
        overrides = case["overrides"]
        expected_overrides = list(AUTHORIZED_COMMON_OVERRIDES) + [
            f"problem/ps_p0={argv_value}"
        ]
        _require(
            type(overrides) is list and overrides == expected_overrides,
            f"{label}.overrides differs from the authorized pressure-pilot launch contract",
        )
        paths: set[str] = set()
        snapshots = case["snapshots"]
        _require(
            type(snapshots) is list and len(snapshots) == len(RAW_CASE_TIMES),
            f"{label}.snapshots must contain the exact five-state schedule",
        )
        for snapshot_index, (snapshot, expected_time) in enumerate(zip(snapshots, RAW_CASE_TIMES)):
            snapshot_label = f"{label}.snapshots[{snapshot_index}]"
            _require_exact_keys(
                snapshot,
                {"time", "mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all"},
                snapshot_label,
            )
            _require(
                type(snapshot["time"]) is float and snapshot["time"] == expected_time,
                f"{snapshot_label}.time differs from the registered schedule",
            )
            suffix = f"{snapshot_index:05d}"
            expected_paths = {
                "mhd_w_bcc": f"cases/{case_id}/bin/{case_id}.mhd_w_bcc.{suffix}.bin",
                "bmag": f"cases/{case_id}/bin/{case_id}.bmag.{suffix}.bin",
                "prtcl_jx": f"cases/{case_id}/bin/{case_id}.prtcl_jx.{suffix}.bin",
                "j2": f"cases/{case_id}/bin/{case_id}.j2.{suffix}.bin",
                "prtcl_all": f"cases/{case_id}/pvtk/{case_id}.prtcl_all.{suffix}.part.vtk",
            }
            for name, expected_path in expected_paths.items():
                binding = _add_expected_binding(
                    expected,
                    snapshot[name],
                    f"{snapshot_label}.{name}",
                    expected_path=expected_path,
                )
                paths.add(binding["path"])
        stdout = _add_expected_binding(
            expected,
            case["stdout"],
            f"{label}.stdout",
            expected_path=f"cases/{case_id}/stdout.txt",
        )
        paths.add(stdout["path"])
        restart = case["terminal_restart"]
        _require_exact_keys(
            restart,
            {"time", "manifest", "manifest_complete", "members"},
            f"{label}.terminal_restart",
        )
        _require(
            type(restart["time"]) is float and restart["time"] == RAW_CASE_TIMES[-1],
            f"{label}.terminal_restart.time differs from the registered schedule",
        )
        restart_manifest_path = f"cases/{case_id}/rst/{case_id}.00004.rst.manifest"
        for name, expected_path in (
            ("manifest", restart_manifest_path),
            ("manifest_complete", restart_manifest_path + ".complete"),
        ):
            binding = _add_expected_binding(
                expected,
                restart[name],
                f"{label}.terminal_restart.{name}",
                expected_path=expected_path,
            )
            paths.add(binding["path"])
        members = restart["members"]
        _require(type(members) is list and bool(members), f"{label}.terminal_restart.members is empty")
        for member_index, member in enumerate(members):
            member_label = f"{label}.terminal_restart.members[{member_index}]"
            _require_exact_keys(member, {"artifact", "complete"}, member_label)
            artifact = _add_expected_binding(expected, member["artifact"], f"{member_label}.artifact")
            complete = _add_expected_binding(expected, member["complete"], f"{member_label}.complete")
            _require(
                artifact["path"].startswith(f"cases/{case_id}/rst/")
                and complete["path"] == artifact["path"] + ".complete",
                f"{member_label} has an unexpected restart-member path",
            )
            paths.update((artifact["path"], complete["path"]))
        case_paths.append(paths)
    return manifest, expected, case_paths


def _expected_directories(paths: set[str]) -> set[str]:
    directories: set[str] = set()
    for relative in paths:
        parent = PurePosixPath(relative).parent
        while parent.as_posix() != ".":
            directories.add(parent.as_posix())
            parent = parent.parent
    return directories


def _scan_bound_directory(
    directory_fd: int,
    prefix: str,
    tree_label: str,
    expected: dict[str, str],
    writable_directories: set[str],
    directories: dict[str, tuple[int, ...]],
    files: dict[str, tuple[tuple[int, ...], int, str]],
    depth: int,
) -> None:
    _require(depth <= MAX_TREE_DEPTH, f"{tree_label} exceeds the tree-depth limit")
    names = _bounded_directory_names(
        directory_fd,
        prefix,
        tree_label,
        require_nonempty=False,
    )
    for name in names:
        _require(
            len(directories) + len(files) < MAX_TREE_ENTRIES,
            f"{tree_label} exceeds the tree-entry limit",
        )
        _require(name not in {"", ".", ".."} and "/" not in name, f"{tree_label} has unsafe entry")
        relative = f"{prefix}/{name}" if prefix else name
        try:
            observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            _fail(f"{tree_label} entry {relative!r} cannot be inspected: {exc}")
        if stat.S_ISREG(observed.st_mode):
            _require(relative in expected, f"{tree_label} has undeclared member {relative!r}")
            with _RetainedFile(directory_fd, name, f"{tree_label} member {relative!r}") as retained:
                payload = retained.read()
                digest = _sha256(payload)
                _require(digest == expected[relative], f"{tree_label} member hash drifted: {relative}")
                files[relative] = (_stable_stat(retained.metadata), len(payload), digest)
        elif stat.S_ISDIR(observed.st_mode):
            if relative not in writable_directories:
                _require_readonly_directory_metadata(observed, f"{tree_label} directory {relative!r}")
            child_fd = _open_child_directory(
                directory_fd,
                name,
                observed,
                f"{tree_label} directory {relative!r}",
            )
            try:
                directories[relative] = _stable_stat(observed)
                _scan_bound_directory(
                    child_fd,
                    relative,
                    tree_label,
                    expected,
                    writable_directories,
                    directories,
                    files,
                    depth + 1,
                )
                _require(
                    _stable_stat(os.fstat(child_fd)) == _stable_stat(observed),
                    f"{tree_label} directory {relative!r} changed during verification",
                )
                _require_entry_identity(
                    directory_fd,
                    name,
                    observed,
                    f"{tree_label} directory {relative!r}",
                )
            finally:
                os.close(child_fd)
        else:
            _fail(f"{tree_label} entry {relative!r} has an unsupported file type")


def _scan_bound_tree(
    directory_fd: int,
    tree_label: str,
    expected: dict[str, str],
    *,
    writable_directories: set[str] | None = None,
) -> dict[str, object]:
    allowed_writable = set() if writable_directories is None else writable_directories
    _require(
        allowed_writable <= _expected_directories(set(expected)),
        f"{tree_label} writable-directory exception is outside expected closure",
    )
    directories: dict[str, tuple[int, ...]] = {}
    files: dict[str, tuple[tuple[int, ...], int, str]] = {}
    _scan_bound_directory(
        directory_fd,
        "",
        tree_label,
        expected,
        allowed_writable,
        directories,
        files,
        0,
    )
    _require(set(files) == set(expected), f"{tree_label} member closure is not exact")
    _require(
        set(directories) == _expected_directories(set(expected)),
        f"{tree_label} directory closure is not exact",
    )
    return {"directories": directories, "files": files}


def _validate_raw_inventory(payload: bytes, expected_sha256: str, label: str) -> dict[str, dict[str, object]]:
    _require(_sha256(payload) == expected_sha256, f"{label} hash differs from aggregate receipt")
    inventory = _decode_canonical_object(payload, label)
    _require_exact_keys(inventory, {"schema_version", "files"}, label)
    _require_int(inventory["schema_version"], 1, f"{label}.schema_version")
    members = inventory["files"]
    _require(type(members) is list and bool(members), f"{label}.files must be nonempty")
    result: dict[str, dict[str, object]] = {}
    observed_paths: list[str] = []
    for index, member in enumerate(members):
        member_label = f"{label}.files[{index}]"
        _require_exact_keys(member, {"path", "sha256", "size"}, member_label)
        path = _relative_source_path(member["path"], f"{member_label}.path")
        digest = _require_sha256(member["sha256"], f"{member_label}.sha256")
        size = member["size"]
        _require(type(size) is int and size >= 0, f"{member_label}.size must be nonnegative")
        _require(path not in result, f"{label}.files contains a duplicate path")
        result[path] = {"sha256": digest, "size": size}
        observed_paths.append(path)
    _require(observed_paths == sorted(observed_paths), f"{label}.files must be sorted")
    return result


def _validate_raw_descriptor(
    payload: bytes,
    case: dict[str, Any],
    case_index: int,
    manifest_case: dict[str, Any],
    aggregate_members: dict[str, tuple[tuple[int, ...], int, str]],
    case_paths: set[str],
    inventory: dict[str, dict[str, object]],
    runtime_payloads: dict[str, bytes],
) -> dict[str, Any]:
    label = f"raw case {case['case_id']!r} descriptor"
    _require(_sha256(payload) == case["descriptor_sha256"], f"{label} hash differs from aggregate receipt")
    descriptor = _decode_canonical_object(payload, label)
    _require_exact_keys(
        descriptor,
        {
            "schema_version",
            "record_type",
            "evidence_class",
            "qualification_effect",
            "launch_contract",
            "case_id",
            "ps_p0",
            "argv_value",
            "artifact_inventory_sha256",
            "runtime_artifacts",
            "runtime_profile",
            "parallel_ranks",
            "rank_gpu_bindings",
            "manifest_case",
            "bundle_members",
        },
        label,
    )
    _require_int(descriptor["schema_version"], 1, f"{label}.schema_version")
    _require_string(descriptor["record_type"], RAW_CASE_RECORD_TYPE, f"{label}.record_type")
    _require_string(descriptor["evidence_class"], AGGREGATE_EVIDENCE_CLASS, f"{label}.evidence_class")
    _require_string(
        descriptor["qualification_effect"],
        AGGREGATE_QUALIFICATION_EFFECT,
        f"{label}.qualification_effect",
    )
    _require_string(
        descriptor["launch_contract"],
        "trusted_trampoline_athena_argv_v1",
        f"{label}.launch_contract",
    )
    _require_string(descriptor["case_id"], RAW_CASE_IDS[case_index], f"{label}.case_id")
    _require(
        type(descriptor["ps_p0"]) is float and descriptor["ps_p0"] == RAW_CASE_PRESSURES[case_index],
        f"{label}.ps_p0 differs from registered pressure",
    )
    _require_string(descriptor["argv_value"], RAW_CASE_ARGV_VALUES[case_index], f"{label}.argv_value")
    _require_string(
        descriptor["artifact_inventory_sha256"],
        case["artifact_inventory_sha256"],
        f"{label}.artifact_inventory_sha256",
    )
    _require(
        _same_json_value(descriptor["runtime_artifacts"], case["runtime_artifacts"]),
        f"{label}.runtime_artifacts differs from aggregate receipt",
    )
    _require(
        _same_json_value(descriptor["manifest_case"], manifest_case),
        f"{label}.manifest_case differs from aggregate bundle manifest",
    )
    _require_string(
        descriptor["runtime_profile"],
        AUTHORIZED_RUNTIME_PROFILE,
        f"{label}.runtime_profile",
    )
    ranks = descriptor["parallel_ranks"]
    rank_bindings = descriptor["rank_gpu_bindings"]
    _require_int(ranks, AUTHORIZED_PARALLEL_RANKS, f"{label}.parallel_ranks")
    _require(
        type(rank_bindings) is list and len(rank_bindings) == AUTHORIZED_PARALLEL_RANKS,
        f"{label}.rank_gpu_bindings drifted",
    )
    for rank, binding in enumerate(rank_bindings):
        rank_label = f"{label}.rank_gpu_bindings[{rank}]"
        _require_exact_keys(binding, {"host", "rank", "rocr_visible_device"}, rank_label)
        host = _require_nonempty_string(binding["host"], f"{rank_label}.host")
        _require(
            _FRONTIER_HOST_RE.fullmatch(host) is not None,
            f"{rank_label}.host must identify a Frontier compute node",
        )
        _require_int(binding["rank"], rank, f"{rank_label}.rank")
        _require_int(
            binding["rocr_visible_device"],
            AUTHORIZED_ROCR_VISIBLE_DEVICE,
            f"{rank_label}.rocr_visible_device",
        )

    bundle_members = descriptor["bundle_members"]
    _require(type(bundle_members) is list and bool(bundle_members), f"{label}.bundle_members must be nonempty")
    observed_paths: list[str] = []
    for member_index, member in enumerate(bundle_members):
        member_label = f"{label}.bundle_members[{member_index}]"
        _require_exact_keys(member, {"path", "sha256", "size", "source_path"}, member_label)
        target = _relative_source_path(member["path"], f"{member_label}.path")
        source = _relative_source_path(member["source_path"], f"{member_label}.source_path")
        digest = _require_sha256(member["sha256"], f"{member_label}.sha256")
        size = member["size"]
        _require(type(size) is int and size >= 0, f"{member_label}.size must be nonnegative")
        _require(source in inventory, f"{member_label}.source_path is absent from raw inventory")
        _require(
            inventory[source] == {"sha256": digest, "size": size},
            f"{member_label} differs from raw inventory",
        )
        _require(target in aggregate_members, f"{member_label}.path is absent from aggregate bundle")
        _require(
            aggregate_members[target][1:] == (size, digest),
            f"{member_label} differs from retained aggregate bundle member",
        )
        observed_paths.append(target)
    _require(observed_paths == sorted(observed_paths), f"{label}.bundle_members must be sorted")
    _require(set(observed_paths) == case_paths, f"{label}.bundle_members closure differs from manifest case")

    for path, digest in case["runtime_artifacts"].items():
        _require(path in inventory, f"{label} runtime artifact {path!r} is absent from raw inventory")
        _require(
            inventory[path]["sha256"] == digest
            and _sha256(runtime_payloads[path]) == digest
            and inventory[path]["size"] == len(runtime_payloads[path]),
            f"{label} runtime artifact differs from raw inventory or retained file: {path}",
        )
    _require(
        runtime_payloads["athena_stdout.sha256"]
        == (_sha256(runtime_payloads["athena_stdout.txt"]) + "\n").encode("ascii"),
        f"{label} stdout checksum file does not bind retained stdout",
    )
    return descriptor


def _validate_review_metrics(
    payload: bytes,
    aggregate_receipt_binding: dict[str, str],
) -> dict[str, Any]:
    metrics = _decode_canonical_object(payload, "pressure review metrics")
    _require_exact_keys(
        metrics,
        {
            "schema_version",
            "record_type",
            "watermark",
            "qualification_effect",
            "aggregate_receipt",
            "cases",
        },
        "pressure review metrics",
    )
    _require_int(metrics["schema_version"], 1, "pressure review metrics.schema_version")
    _require_string(
        metrics["record_type"],
        REVIEW_METRICS_RECORD_TYPE,
        "pressure review metrics.record_type",
    )
    _require_string(metrics["watermark"], WATERMARK, "pressure review metrics.watermark")
    _require_string(
        metrics["qualification_effect"],
        QUALIFICATION_EFFECT,
        "pressure review metrics.qualification_effect",
    )
    metrics_aggregate_binding = _binding(
        metrics["aggregate_receipt"],
        Path(aggregate_receipt_binding["path"]).parent,
        "pressure review metrics.aggregate_receipt",
    )
    _require(
        _same_json_value(metrics_aggregate_binding, aggregate_receipt_binding),
        "pressure review metrics does not bind the supplied aggregate receipt",
    )

    cases = metrics["cases"]
    _require(
        type(cases) is list and len(cases) == len(RAW_CASE_IDS),
        "pressure review metrics.cases must contain the exact ordered four-case set",
    )
    for index, case in enumerate(cases):
        label = f"pressure review metrics.cases[{index}]"
        _require_exact_keys(
            case,
            {
                "case_id",
                "problem_ps_p0",
                "terminal_particle_count",
                "particle_efficiency",
                "zone_cycles_per_second",
                "particle_updates_per_second",
                "tracked_gpu_memory_high_water_bytes",
            },
            label,
        )
        _require_string(case["case_id"], RAW_CASE_IDS[index], f"{label}.case_id")
        problem_ps_p0 = _require_nonnegative_number(
            case["problem_ps_p0"],
            f"{label}.problem_ps_p0",
        )
        _require(
            problem_ps_p0 == RAW_CASE_PRESSURES[index],
            f"{label}.problem_ps_p0 differs from registered pressure",
        )
        _require(
            type(case["terminal_particle_count"]) is int
            and case["terminal_particle_count"] >= 0,
            f"{label}.terminal_particle_count must be a nonnegative integer",
        )
        for key in (
            "particle_efficiency",
            "zone_cycles_per_second",
            "particle_updates_per_second",
            "tracked_gpu_memory_high_water_bytes",
        ):
            _require_nonnegative_number(case[key], f"{label}.{key}")
    return metrics


def _validate_packet_tree(
    scan: dict[str, object],
    inventory_sha256: str,
    aggregate_receipt_binding: dict[str, str],
) -> dict[str, Any]:
    payloads = scan["payloads"]
    _require(type(payloads) is dict, "internal packet scan payload map is invalid")
    expected = set(PACKET_MEMBERS) | {INVENTORY_NAME}
    _require(set(payloads) == expected, "packet inventory/member closure is not exact")

    inventory_payload = payloads[INVENTORY_NAME]
    _require(
        _sha256(inventory_payload) == inventory_sha256,
        "packet inventory hash does not match the packet receipt",
    )
    inventory = _decode_canonical_object(inventory_payload, "packet inventory")
    _require_exact_keys(inventory, {"schema_version", "record_type", "members"}, "packet inventory")
    _require_int(inventory["schema_version"], 1, "packet inventory.schema_version")
    _require_string(inventory["record_type"], INVENTORY_RECORD_TYPE, "packet inventory.record_type")
    members = inventory["members"]
    _require(
        type(members) is list and len(members) == len(PACKET_MEMBERS),
        "packet inventory.members has the wrong length",
    )
    expected_paths = sorted(PACKET_MEMBERS)
    observed_paths: list[str] = []
    for index, member in enumerate(members):
        label = f"packet inventory.members[{index}]"
        _require_exact_keys(member, {"path", "sha256", "size"}, label)
        path = _relative_source_path(member["path"], f"{label}.path")
        _require(path in PACKET_MEMBERS, f"{label}.path is not an expected packet member")
        digest = _require_sha256(member["sha256"], f"{label}.sha256")
        size = member["size"]
        _require(type(size) is int and size >= 0, f"{label}.size must be a nonnegative integer")
        payload = payloads[path]
        _require(len(payload) == size, f"{label}.size does not match the packet member")
        _require(_sha256(payload) == digest, f"{label}.sha256 does not match the packet member")
        observed_paths.append(path)
    _require(observed_paths == expected_paths, "packet inventory.members must be sorted and exact")
    _require(set(scan["directories"]) == {"figures"}, "packet directory closure is not exact")
    _validate_review_metrics(
        payloads["pressure_review_metrics.json"],
        aggregate_receipt_binding,
    )
    return inventory


def _verify_aggregate_artifacts(
    stack: contextlib.ExitStack,
    publication_fd: int,
    runs_fd: int,
    authorized_pic_root: Path,
    aggregate_value: dict[str, Any],
    aggregate_bundle: dict[str, str],
    aggregate_analysis: dict[str, str],
) -> tuple[dict[str, Any], dict[str, object]]:
    with _retained_child_directory(
        publication_fd,
        Path(aggregate_bundle["path"]).name,
        "aggregate bundle",
    ) as (bundle_fd, bundle_metadata):
        _require_readonly_directory_metadata(bundle_metadata, "aggregate bundle")
        manifest_file = stack.enter_context(
            _RetainedFile(bundle_fd, AGGREGATE_MANIFEST_NAME, "aggregate bundle manifest")
        )
        manifest_payload = manifest_file.read()
        manifest, expected, case_paths = _validate_aggregate_manifest(
            manifest_payload,
            aggregate_bundle["manifest_sha256"],
            aggregate_value["source_bindings"],
        )
        first_bundle_scan = _scan_bound_tree(bundle_fd, "aggregate bundle", expected)
        second_bundle_scan = _scan_bound_tree(bundle_fd, "aggregate bundle", expected)
        _require(
            first_bundle_scan == second_bundle_scan,
            "aggregate bundle changed between stable verification reads",
        )
        _require(
            manifest_file.read() == manifest_payload,
            "aggregate bundle manifest changed during verification",
        )

        aggregate_members = first_bundle_scan["files"]
        _require(type(aggregate_members) is dict, "internal aggregate member map is invalid")
        for case_index, case in enumerate(aggregate_value["raw_cases"]):
            case_label = f"aggregate receipt.raw_cases[{case_index}]"
            relative = _raw_case_relative_path(case, authorized_pic_root, case_label)
            with _retained_relative_directory(
                runs_fd,
                relative,
                f"raw case {case['case_id']!r} artifact directory",
            ) as (raw_fd, raw_metadata):
                _require_readonly_directory_metadata(
                    raw_metadata,
                    f"raw case {case['case_id']!r} artifact directory",
                )
                raw_stack = contextlib.ExitStack()
                try:
                    descriptor_file = _enter_retained_relative_file(
                        raw_stack,
                        raw_fd,
                        case["descriptor_path"],
                        f"raw case {case['case_id']!r} descriptor",
                    )
                    inventory_file = _enter_retained_relative_file(
                        raw_stack,
                        raw_fd,
                        RAW_CASE_INVENTORY_NAME,
                        f"raw case {case['case_id']!r} inventory",
                    )
                    runtime_files = {
                        path: _enter_retained_relative_file(
                            raw_stack,
                            raw_fd,
                            path,
                            f"raw case {case['case_id']!r} runtime artifact {path!r}",
                        )
                        for path in case["runtime_artifacts"]
                    }
                    descriptor_payload = descriptor_file.read()
                    inventory_payload = inventory_file.read()
                    runtime_payloads = {
                        path: retained.read() for path, retained in runtime_files.items()
                    }
                    raw_inventory = _validate_raw_inventory(
                        inventory_payload,
                        case["artifact_inventory_sha256"],
                        f"raw case {case['case_id']!r} inventory",
                    )
                    expected_raw = {
                        path: str(binding["sha256"])
                        for path, binding in raw_inventory.items()
                    }
                    expected_raw[RAW_CASE_INVENTORY_NAME] = case["artifact_inventory_sha256"]
                    expected_raw[case["descriptor_path"]] = case["descriptor_sha256"]
                    first_raw_scan = _scan_bound_tree(
                        raw_fd,
                        f"raw case {case['case_id']!r}",
                        expected_raw,
                        writable_directories={"analysis"},
                    )
                    second_raw_scan = _scan_bound_tree(
                        raw_fd,
                        f"raw case {case['case_id']!r}",
                        expected_raw,
                        writable_directories={"analysis"},
                    )
                    _require(
                        first_raw_scan == second_raw_scan,
                        f"raw case {case['case_id']!r} changed between stable verification reads",
                    )
                    raw_files = first_raw_scan["files"]
                    _require(type(raw_files) is dict, "internal raw-case member map is invalid")
                    _require(
                        all(
                            raw_files[path][1:]
                            == (binding["size"], binding["sha256"])
                            for path, binding in raw_inventory.items()
                        ),
                        f"raw case {case['case_id']!r} inventory size/hash closure drifted",
                    )
                    _validate_raw_descriptor(
                        descriptor_payload,
                        case,
                        case_index,
                        manifest["cases"][case_index],
                        aggregate_members,
                        case_paths[case_index],
                        raw_inventory,
                        runtime_payloads,
                    )
                    _require(
                        descriptor_file.read() == descriptor_payload
                        and inventory_file.read() == inventory_payload
                        and all(
                            runtime_files[path].read() == payload
                            for path, payload in runtime_payloads.items()
                        ),
                        f"raw case {case['case_id']!r} provenance files changed during verification",
                    )
                finally:
                    raw_stack.close()

        analysis_file = stack.enter_context(
            _RetainedFile(
                publication_fd,
                Path(aggregate_analysis["path"]).name,
                "aggregate analysis",
            )
        )
        analysis_payload = analysis_file.read()
        _require(
            _sha256(analysis_payload) == aggregate_analysis["sha256"],
            "aggregate analysis hash does not match aggregate receipt",
        )
        _decode_canonical_object(analysis_payload, "aggregate analysis")
        _require(
            analysis_file.read() == analysis_payload,
            "aggregate analysis changed during verification",
        )
        return manifest, first_bundle_scan


def _packet_receipt_binding_from_path(
    receipt_path: str | os.PathLike[str],
    *,
    authorized_pic_root: str | os.PathLike[str],
) -> dict[str, str]:
    root = _canonical_absolute_path(authorized_pic_root, "authorized PIC root")
    publication_root = root / "publication"
    normalized = _direct_publication_argument_path(
        receipt_path,
        publication_root,
        "packet receipt path",
    )
    with _retained_absolute_directory(root, "authorized PIC root") as (
        _root,
        root_fd,
        _root_metadata,
    ), _retained_child_directory(root_fd, "publication", "publication root") as (
        publication_fd,
        _publication_metadata,
    ), _RetainedFile(publication_fd, normalized.name, "packet receipt") as receipt:
        payload = receipt.read()
    return {"path": str(normalized), "sha256": _sha256(payload)}


def _verify_pressure_review_packet_binding(
    packet_receipt_binding: object,
    aggregate_receipt_binding: object,
    *,
    authorized_pic_root: str | os.PathLike[str],
) -> dict[str, object]:
    root = _canonical_absolute_path(authorized_pic_root, "authorized PIC root")
    publication_root = root / "publication"

    packet_binding = _binding(packet_receipt_binding, publication_root, "packet receipt binding")
    aggregate_binding = _binding(
        aggregate_receipt_binding,
        publication_root,
        "aggregate receipt binding",
    )
    if root == AUTHORIZED_PRODUCTION_PIC_ROOT:
        _require(
            _same_json_value(
                packet_binding,
                AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING,
            ),
            "packet receipt binding is not the authorized immutable production receipt",
        )
        _require(
            _same_json_value(
                aggregate_binding,
                AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING,
            ),
            "aggregate receipt binding is not the authorized immutable production receipt",
        )
    packet_name = Path(packet_binding["path"]).name
    aggregate_name = Path(aggregate_binding["path"]).name

    with _retained_absolute_directory(root, "authorized PIC root") as (
        _root,
        root_fd,
        _root_metadata,
    ), contextlib.ExitStack() as stack:
        publication_fd, publication_metadata = stack.enter_context(
            _retained_child_directory(root_fd, "publication", "publication root")
        )
        acceptance_fd, _acceptance_metadata = stack.enter_context(
            _retained_child_directory(
                root_fd,
                "publication_acceptance",
                "publication acceptance root",
            )
        )
        runs_fd, _runs_metadata = stack.enter_context(
            _retained_child_directory(root_fd, "runs", "authorized runs root")
        )
        _require_guard_absent(publication_fd, packet_name, "packet receipt")
        _require_guard_absent(publication_fd, aggregate_name, "aggregate receipt")

        packet_receipt = stack.enter_context(
            _RetainedFile(publication_fd, packet_name, "packet receipt")
        )
        aggregate_receipt = stack.enter_context(
            _RetainedFile(publication_fd, aggregate_name, "aggregate receipt")
        )
        packet_payload = packet_receipt.read()
        aggregate_payload = aggregate_receipt.read()
        _require(
            _sha256(packet_payload) == packet_binding["sha256"],
            "packet receipt binding hash mismatch",
        )
        _require(
            _sha256(aggregate_payload) == aggregate_binding["sha256"],
            "aggregate receipt binding hash mismatch",
        )

        packet_seal, packet_seal_payload = _open_and_verify_seal(
            stack,
            acceptance_fd,
            packet_name,
            packet_binding["sha256"],
            packet_receipt.metadata,
            publication_metadata,
            "packet receipt",
        )
        aggregate_seal, aggregate_seal_payload = _open_and_verify_seal(
            stack,
            acceptance_fd,
            aggregate_name,
            aggregate_binding["sha256"],
            aggregate_receipt.metadata,
            publication_metadata,
            "aggregate receipt",
        )

        aggregate_value = _decode_canonical_object(aggregate_payload, "aggregate receipt")
        packet_value = _decode_canonical_object(packet_payload, "packet receipt")
        aggregate_bundle, aggregate_analysis = _validate_aggregate_receipt(
            aggregate_value,
            publication_root,
            publication_metadata,
            root,
        )
        embedded_aggregate, packet_root, inventory_sha256 = _validate_packet_receipt(
            packet_value,
            publication_root,
            publication_metadata,
            root,
        )
        _require(
            _same_json_value(embedded_aggregate, aggregate_binding),
            "packet receipt does not bind the supplied aggregate receipt",
        )
        _require(
            _same_json_value(packet_value["source_bindings"], aggregate_value["source_bindings"]),
            "packet and aggregate receipt source bindings differ",
        )
        _verify_aggregate_artifacts(
            stack,
            publication_fd,
            runs_fd,
            root,
            aggregate_value,
            aggregate_bundle,
            aggregate_analysis,
        )
        with _retained_child_directory(
            publication_fd,
            packet_root.name,
            "packet root",
        ) as (
            packet_root_fd,
            packet_root_metadata,
        ):
            _require_readonly_directory_metadata(packet_root_metadata, "packet root")
            first_scan = _scan_packet_tree(packet_root_fd)
            inventory = _validate_packet_tree(
                first_scan,
                inventory_sha256,
                aggregate_binding,
            )
            second_scan = _scan_packet_tree(packet_root_fd)
            second_inventory = _validate_packet_tree(
                second_scan,
                inventory_sha256,
                aggregate_binding,
            )
            _require(
                first_scan == second_scan,
                "packet tree changed between stable verification reads",
            )
            _require(
                _same_json_value(inventory, second_inventory),
                "packet inventory changed between stable verification reads",
            )
            _require(
                packet_receipt.read() == packet_payload,
                "packet receipt changed during verification",
            )
            _require(
                aggregate_receipt.read() == aggregate_payload,
                "aggregate receipt changed during verification",
            )
            _require(
                packet_seal.read() == packet_seal_payload,
                "packet receipt success seal changed during verification",
            )
            _require(
                aggregate_seal.read() == aggregate_seal_payload,
                "aggregate receipt success seal changed during verification",
            )
            _require_guard_absent(publication_fd, packet_name, "packet receipt")
            _require_guard_absent(publication_fd, aggregate_name, "aggregate receipt")

    return {
        "packet_receipt_binding": packet_binding,
        "aggregate_receipt_binding": aggregate_binding,
        "packet_root": str(packet_root),
        "inventory_sha256": inventory_sha256,
        "source_bindings": packet_value["source_bindings"],
        "packet_receipt": packet_value,
        "aggregate_receipt": aggregate_value,
        "aggregate_bundle": aggregate_bundle,
        "aggregate_analysis": aggregate_analysis,
        "inventory": inventory,
    }


def verify_pressure_review_packet_binding(
    packet_receipt_binding: object,
    aggregate_receipt_binding: object,
    *,
    authorized_pic_root: str | os.PathLike[str],
) -> dict[str, object]:
    """Verify installed packet and aggregate-receipt contracts without mutation."""

    try:
        verified = _verify_pressure_review_packet_binding(
            packet_receipt_binding,
            aggregate_receipt_binding,
            authorized_pic_root=authorized_pic_root,
        )
        return {
            "packet_receipt": verified["packet_receipt_binding"],
            "aggregate_receipt": verified["aggregate_receipt_binding"],
            "packet_root": verified["packet_root"],
            "inventory_sha256": verified["inventory_sha256"],
            "source_bindings": verified["source_bindings"],
        }
    except PressureReviewPacketVerificationError:
        raise
    except (OSError, RecursionError, TypeError, ValueError) as exc:
        raise PressureReviewPacketVerificationError(
            f"pressure review packet verification failed closed: {exc}"
        ) from exc


def consume_published_pressure_pilot_review_packet(
    receipt_path: str | os.PathLike[str],
    *,
    aggregate_receipt_binding: object,
    authorized_pic_root: str | os.PathLike[str],
) -> dict[str, object]:
    """Return normalized bindings and detached parsed records for one accepted packet.

    Raises PressureReviewPacketVerificationError for every rejected input or artifact.
    """

    try:
        packet_receipt_binding = _packet_receipt_binding_from_path(
            receipt_path,
            authorized_pic_root=authorized_pic_root,
        )
        verified = _verify_pressure_review_packet_binding(
            packet_receipt_binding,
            aggregate_receipt_binding,
            authorized_pic_root=authorized_pic_root,
        )
        return {
            "receipt_binding": copy.deepcopy(verified["packet_receipt_binding"]),
            "aggregate_receipt_binding": copy.deepcopy(
                verified["aggregate_receipt_binding"]
            ),
            "packet_receipt": copy.deepcopy(verified["packet_receipt"]),
            "aggregate_receipt": copy.deepcopy(verified["aggregate_receipt"]),
            "aggregate_bundle": copy.deepcopy(verified["aggregate_bundle"]),
            "aggregate_analysis": copy.deepcopy(verified["aggregate_analysis"]),
            "source_bindings": copy.deepcopy(verified["source_bindings"]),
            "inventory": copy.deepcopy(verified["inventory"]),
        }
    except PressureReviewPacketVerificationError:
        raise
    except (OSError, RecursionError, TypeError, ValueError) as exc:
        raise PressureReviewPacketVerificationError(
            f"pressure review packet consumption failed closed: {exc}"
        ) from exc


def _read_pressure_gate_attestation(
    value: object,
    *,
    authorized_pic_root: str | os.PathLike[str],
    label: str,
) -> tuple[dict[str, str], dict[str, Any]]:
    root = _canonical_absolute_path(authorized_pic_root, "authorized PIC root")
    binding = _pressure_gate_binding(value, authorized_pic_root=root, label=label)
    path = Path(binding["path"])
    with _retained_absolute_directory(root, "authorized PIC root") as (
        _,
        root_fd,
        _,
    ):
        with _retained_child_directory(
            root_fd,
            PRESSURE_GATE_ATTESTATION_ROOT_NAME,
            "pressure-gate attestation archive root",
        ) as (archive_fd, _):
            with _retained_child_directory(
                archive_fd,
                path.parent.name,
                label,
            ) as (attestation_fd, attestation_metadata):
                _require(
                    stat.S_IMODE(attestation_metadata.st_mode) == 0o500,
                    f"{label} directory mode must be 0500",
                )
                try:
                    first_names = os.listdir(attestation_fd)
                except OSError as exc:
                    _fail(f"{label} tree cannot be listed: {exc}")
                _require(
                    first_names == [PRESSURE_GATE_ATTESTATION_FILENAME],
                    f"{label} tree must contain only attestation.json",
                )
                with _RetainedFile(
                    attestation_fd,
                    PRESSURE_GATE_ATTESTATION_FILENAME,
                    f"{label} payload",
                ) as retained:
                    _require(
                        stat.S_IMODE(retained.metadata.st_mode) == 0o400,
                        f"{label} payload mode must be 0400",
                    )
                    payload = retained.read()
                    _require(
                        _sha256(payload) == binding["sha256"],
                        f"{label} binding hash mismatch",
                    )
                    attestation = _decode_canonical_object(payload, label)
                    _require(
                        retained.read() == payload,
                        f"{label} payload changed during verification",
                    )
                try:
                    second_names = os.listdir(attestation_fd)
                except OSError as exc:
                    _fail(f"{label} tree cannot be re-listed: {exc}")
                _require(
                    second_names == first_names,
                    f"{label} tree changed during verification",
                )
    return binding, attestation


def _validate_reanalysis_source_authorization(value: object) -> dict[str, object]:
    label = "authoritative reanalysis attestation.source_authorization"
    _require_exact_keys(
        value,
        {
            "execution_mode",
            "git_commit",
            "source_archive_sha256",
            "source_closure_sha256",
            "source_closure",
            "historical_production_source_authorization",
        },
        label,
    )
    _require_string(
        value["execution_mode"],
        PRESSURE_REANALYSIS_EXECUTION_MODE,
        f"{label}.execution_mode",
    )
    git_commit = _require_nonempty_string(value["git_commit"], f"{label}.git_commit")
    _require(
        _GIT_COMMIT_RE.fullmatch(git_commit) is not None,
        f"{label}.git_commit must be a lowercase Git commit",
    )
    archive_sha256 = _require_sha256(
        value["source_archive_sha256"],
        f"{label}.source_archive_sha256",
    )
    closure_sha256 = _require_sha256(
        value["source_closure_sha256"],
        f"{label}.source_closure_sha256",
    )
    closure = value["source_closure"]
    _require(
        type(closure) is list and len(closure) == len(PRESSURE_REANALYSIS_SOURCE_PATHS),
        f"{label}.source_closure must contain the exact reanalysis source set",
    )
    normalized_closure: list[dict[str, str]] = []
    for index, (record, expected_path) in enumerate(
        zip(closure, PRESSURE_REANALYSIS_SOURCE_PATHS)
    ):
        record_label = f"{label}.source_closure[{index}]"
        _require_exact_keys(record, {"path", "sha256"}, record_label)
        path = _relative_source_path(record["path"], f"{record_label}.path")
        _require(path == expected_path, f"{record_label}.path drifted")
        normalized_closure.append(
            {
                "path": path,
                "sha256": _require_sha256(
                    record["sha256"],
                    f"{record_label}.sha256",
                ),
            }
        )
    _require(
        closure_sha256 == _source_closure_digest(normalized_closure),
        f"{label}.source_closure_sha256 drifted",
    )
    historical = value["historical_production_source_authorization"]
    _require_exact_keys(historical, {"path", "sha256"}, f"{label}.historical_production_source_authorization")
    normalized_historical = {
        "path": _relative_source_path(
            historical["path"],
            f"{label}.historical_production_source_authorization.path",
        ),
        "sha256": _require_sha256(
            historical["sha256"],
            f"{label}.historical_production_source_authorization.sha256",
        ),
    }
    _require(
        _same_json_value(
            normalized_historical,
            AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION,
        ),
        f"{label}.historical_production_source_authorization drifted",
    )
    return {
        "execution_mode": PRESSURE_REANALYSIS_EXECUTION_MODE,
        "git_commit": git_commit,
        "source_archive_sha256": archive_sha256,
        "source_closure_sha256": closure_sha256,
        "source_closure": normalized_closure,
        "historical_production_source_authorization": normalized_historical,
    }


def _normalized_reanalysis_result(value: object, label: str) -> dict[str, str]:
    _require_exact_keys(
        value,
        {
            "packet_receipt_sha256",
            "aggregate_receipt_sha256",
            "manifest_sha256",
            "analysis_result_sha256",
            "status",
        },
        label,
    )
    result = {
        key: _require_sha256(value[key], f"{label}.{key}")
        for key in (
            "packet_receipt_sha256",
            "aggregate_receipt_sha256",
            "manifest_sha256",
            "analysis_result_sha256",
        )
    }
    _require_string(
        value["status"],
        "pass_engineering_calibration_only",
        f"{label}.status",
    )
    result["status"] = "pass_engineering_calibration_only"
    return result


def consume_sealed_pressure_reanalysis_attestation(
    attestation_binding: object,
    *,
    aggregate_receipt_binding: object,
    packet_receipt_binding: object,
    pilot_bundle_manifest_sha256: object,
    aggregate_pilot_analysis_sha256: object,
    authorized_pic_root: str | os.PathLike[str],
    expected_result: object | None = None,
    now: datetime | None = None,
) -> dict[str, object]:
    """Verify one sealed authoritative-reanalysis attestation without mutation."""
    try:
        root = _canonical_absolute_path(authorized_pic_root, "authorized PIC root")
        aggregate_binding = _published_binding(
            aggregate_receipt_binding,
            root,
            "aggregate receipt binding",
        )
        packet_binding = _published_binding(
            packet_receipt_binding,
            root,
            "packet receipt binding",
        )
        manifest_sha256 = _require_sha256(
            pilot_bundle_manifest_sha256,
            "pilot bundle manifest SHA-256",
        )
        analysis_sha256 = _require_sha256(
            aggregate_pilot_analysis_sha256,
            "aggregate pilot analysis SHA-256",
        )
        binding, attestation = _read_pressure_gate_attestation(
            attestation_binding,
            authorized_pic_root=root,
            label="authoritative reanalysis attestation",
        )
        _require_exact_keys(
            attestation,
            {
                "schema_version",
                "record_type",
                "qualification_effect",
                "operator_id",
                "recomputed_utc",
                "sealed_utc",
                "operator_statement",
                "evidence",
                "source_authorization",
                "result",
            },
            "authoritative reanalysis attestation",
        )
        _require_int(
            attestation["schema_version"],
            1,
            "authoritative reanalysis attestation.schema_version",
        )
        _require_string(
            attestation["record_type"],
            PRESSURE_REANALYSIS_RECORD_TYPE,
            "authoritative reanalysis attestation.record_type",
        )
        _require_string(
            attestation["qualification_effect"],
            PRESSURE_REANALYSIS_QUALIFICATION_EFFECT,
            "authoritative reanalysis attestation.qualification_effect",
        )
        operator_id = _pressure_gate_id(
            attestation["operator_id"],
            "authoritative reanalysis attestation.operator_id",
        )
        _require_string(
            attestation["operator_statement"],
            PRESSURE_REANALYSIS_OPERATOR_STATEMENT,
            "authoritative reanalysis attestation.operator_statement",
        )
        recomputed_text, recomputed = _canonical_utc(
            attestation["recomputed_utc"],
            "authoritative reanalysis attestation.recomputed_utc",
        )
        sealed_text, sealed = _canonical_utc(
            attestation["sealed_utc"],
            "authoritative reanalysis attestation.sealed_utc",
        )
        _, not_before = _canonical_utc(
            PRESSURE_GATE_REVIEW_NOT_BEFORE_UTC,
            "pressure gate review not-before",
        )
        current = _current_utc(now)
        _require(
            not_before <= recomputed <= sealed <= current,
            "authoritative reanalysis attestation timestamps are out of order",
        )
        _require(
            (sealed - recomputed).total_seconds()
            <= PRESSURE_GATE_CAPTURE_TO_SEAL_MAX_SECONDS,
            "authoritative reanalysis attestation capture-to-seal interval is stale",
        )
        expected_directory = (
            f"{sealed.strftime('%Y%m%dT%H%M%SZ')}-"
            f"q011-section54-pressure-reanalysis-{operator_id}"
        )
        _require(
            Path(binding["path"]).parent.name == expected_directory,
            "authoritative reanalysis attestation directory name differs",
        )
        evidence = attestation["evidence"]
        _require_exact_keys(
            evidence,
            {
                "published_pressure_pilot_receipt",
                "published_pressure_pilot_review_packet_receipt",
                "pilot_bundle_manifest_sha256",
                "aggregate_pilot_analysis_sha256",
            },
            "authoritative reanalysis attestation.evidence",
        )
        normalized_evidence = {
            "published_pressure_pilot_receipt": _published_binding(
                evidence["published_pressure_pilot_receipt"],
                root,
                "authoritative reanalysis attestation.evidence.published_pressure_pilot_receipt",
            ),
            "published_pressure_pilot_review_packet_receipt": _published_binding(
                evidence["published_pressure_pilot_review_packet_receipt"],
                root,
                "authoritative reanalysis attestation.evidence.published_pressure_pilot_review_packet_receipt",
            ),
            "pilot_bundle_manifest_sha256": _require_sha256(
                evidence["pilot_bundle_manifest_sha256"],
                "authoritative reanalysis attestation.evidence.pilot_bundle_manifest_sha256",
            ),
            "aggregate_pilot_analysis_sha256": _require_sha256(
                evidence["aggregate_pilot_analysis_sha256"],
                "authoritative reanalysis attestation.evidence.aggregate_pilot_analysis_sha256",
            ),
        }
        expected_evidence = {
            "published_pressure_pilot_receipt": aggregate_binding,
            "published_pressure_pilot_review_packet_receipt": packet_binding,
            "pilot_bundle_manifest_sha256": manifest_sha256,
            "aggregate_pilot_analysis_sha256": analysis_sha256,
        }
        _require(
            _same_json_value(normalized_evidence, expected_evidence),
            "authoritative reanalysis attestation evidence tuple drifted",
        )
        source_authorization = _validate_reanalysis_source_authorization(
            attestation["source_authorization"]
        )
        result = _normalized_reanalysis_result(
            attestation["result"],
            "authoritative reanalysis attestation.result",
        )
        required_result = _expected_reanalysis_result(
            packet_receipt_binding=packet_binding,
            aggregate_receipt_binding=aggregate_binding,
            pilot_bundle_manifest_sha256=manifest_sha256,
            aggregate_pilot_analysis_sha256=analysis_sha256,
        )
        if expected_result is not None:
            supplied_result = _normalized_reanalysis_result(
                expected_result,
                "expected authoritative reanalysis result",
            )
            _require(
                _same_json_value(supplied_result, required_result),
                "expected authoritative reanalysis result differs from publication evidence",
            )
        _require(
            _same_json_value(result, required_result),
            "authoritative reanalysis attestation result tuple drifted",
        )
        return {
            "binding": copy.deepcopy(binding),
            "attestation": copy.deepcopy(attestation),
            "operator_id": operator_id,
            "recomputed_utc": recomputed_text,
            "sealed_utc": sealed_text,
            "evidence": copy.deepcopy(normalized_evidence),
            "source_authorization": copy.deepcopy(source_authorization),
            "result": copy.deepcopy(result),
        }
    except PressureReviewPacketVerificationError:
        raise
    except (OSError, RecursionError, TypeError, ValueError) as exc:
        raise PressureReviewPacketVerificationError(
            f"authoritative reanalysis attestation verification failed closed: {exc}"
        ) from exc


def _reconsume_sealed_pressure_reanalysis_verification(
    value: object,
    *,
    authorized_pic_root: str | os.PathLike[str],
    aggregate_receipt_binding: object | None = None,
    packet_receipt_binding: object | None = None,
    now: datetime | None = None,
) -> dict[str, object]:
    label = "authoritative reanalysis verification result"
    _require_exact_keys(
        value,
        {
            "binding",
            "attestation",
            "operator_id",
            "recomputed_utc",
            "sealed_utc",
            "evidence",
            "source_authorization",
            "result",
        },
        label,
    )
    evidence = value["evidence"]
    _require_exact_keys(
        evidence,
        {
            "published_pressure_pilot_receipt",
            "published_pressure_pilot_review_packet_receipt",
            "pilot_bundle_manifest_sha256",
            "aggregate_pilot_analysis_sha256",
        },
        f"{label}.evidence",
    )
    root = _canonical_absolute_path(authorized_pic_root, "authorized PIC root")
    reconsumed = consume_sealed_pressure_reanalysis_attestation(
        value["binding"],
        aggregate_receipt_binding=(
            evidence["published_pressure_pilot_receipt"]
            if aggregate_receipt_binding is None
            else aggregate_receipt_binding
        ),
        packet_receipt_binding=(
            evidence["published_pressure_pilot_review_packet_receipt"]
            if packet_receipt_binding is None
            else packet_receipt_binding
        ),
        pilot_bundle_manifest_sha256=evidence["pilot_bundle_manifest_sha256"],
        aggregate_pilot_analysis_sha256=evidence[
            "aggregate_pilot_analysis_sha256"
        ],
        authorized_pic_root=root,
        now=now,
    )
    _require(
        _same_json_value(value, reconsumed),
        f"{label} differs from re-consumed sealed attestation",
    )
    return reconsumed


def consume_sealed_pressure_reviewer_attestation(
    attestation_binding: object,
    *,
    aggregate_receipt_binding: object,
    packet_receipt_binding: object,
    reanalysis_verification: object,
    selected_case: object,
    authorized_pic_root: str | os.PathLike[str],
    now: datetime | None = None,
) -> dict[str, object]:
    """Verify one sealed human pressure-selection reviewer attestation."""
    try:
        root = _canonical_absolute_path(authorized_pic_root, "authorized PIC root")
        aggregate_binding = _published_binding(
            aggregate_receipt_binding,
            root,
            "aggregate receipt binding",
        )
        packet_binding = _published_binding(
            packet_receipt_binding,
            root,
            "packet receipt binding",
        )
        reanalysis = _reconsume_sealed_pressure_reanalysis_verification(
            reanalysis_verification,
            authorized_pic_root=root,
            aggregate_receipt_binding=aggregate_binding,
            packet_receipt_binding=packet_binding,
            now=now,
        )
        reanalysis_binding = _pressure_gate_binding(
            reanalysis["binding"],
            authorized_pic_root=root,
            label="authoritative reanalysis attestation binding",
        )
        _, reanalysis_sealed = _canonical_utc(
            reanalysis["sealed_utc"],
            "authoritative reanalysis verification sealed_utc",
        )
        _require_exact_keys(selected_case, {"case_id", "problem_ps_p0"}, "selected case")
        case_id = _require_nonempty_string(selected_case["case_id"], "selected case.case_id")
        pressure = selected_case["problem_ps_p0"]
        _require(
            type(pressure) is float
            and (case_id, pressure)
            in {
                ("ps_p0_1p00", 1.0),
                ("ps_p0_0p05", 0.05),
                ("ps_p0_0p10", 0.1),
                ("ps_p0_0p20", 0.2),
            },
            "selected case is not one exact registered pressure case",
        )
        normalized_case = {"case_id": case_id, "problem_ps_p0": pressure}
        binding, attestation = _read_pressure_gate_attestation(
            attestation_binding,
            authorized_pic_root=root,
            label="pressure reviewer attestation",
        )
        _require_exact_keys(
            attestation,
            {
                "schema_version",
                "record_type",
                "qualification_effect",
                "selection_method",
                "reviewer_id",
                "reviewed_utc",
                "sealed_utc",
                "rationale",
                "reviewer_statement",
                "published_pressure_pilot_receipt",
                "published_pressure_pilot_review_packet_receipt",
                "authoritative_reanalysis_attestation",
                "selected_case",
            },
            "pressure reviewer attestation",
        )
        _require_int(
            attestation["schema_version"],
            1,
            "pressure reviewer attestation.schema_version",
        )
        _require_string(
            attestation["record_type"],
            PRESSURE_REVIEWER_RECORD_TYPE,
            "pressure reviewer attestation.record_type",
        )
        _require_string(
            attestation["qualification_effect"],
            PRESSURE_REVIEWER_QUALIFICATION_EFFECT,
            "pressure reviewer attestation.qualification_effect",
        )
        _require_string(
            attestation["selection_method"],
            "human_review_only",
            "pressure reviewer attestation.selection_method",
        )
        reviewer_id = _pressure_gate_id(
            attestation["reviewer_id"],
            "pressure reviewer attestation.reviewer_id",
        )
        rationale = _canonical_rationale(
            attestation["rationale"],
            "pressure reviewer attestation.rationale",
        )
        _require_string(
            attestation["reviewer_statement"],
            PRESSURE_REVIEWER_STATEMENT,
            "pressure reviewer attestation.reviewer_statement",
        )
        reviewed_text, reviewed = _canonical_utc(
            attestation["reviewed_utc"],
            "pressure reviewer attestation.reviewed_utc",
        )
        sealed_text, sealed = _canonical_utc(
            attestation["sealed_utc"],
            "pressure reviewer attestation.sealed_utc",
        )
        _, not_before = _canonical_utc(
            PRESSURE_GATE_REVIEW_NOT_BEFORE_UTC,
            "pressure gate review not-before",
        )
        current = _current_utc(now)
        _require(
            not_before <= reanalysis_sealed <= reviewed <= sealed <= current,
            "pressure reviewer attestation timestamps are out of order",
        )
        _require(
            (sealed - reviewed).total_seconds()
            <= PRESSURE_GATE_CAPTURE_TO_SEAL_MAX_SECONDS,
            "pressure reviewer attestation review-to-seal interval is stale",
        )
        expected_directory = (
            f"{sealed.strftime('%Y%m%dT%H%M%SZ')}-"
            f"q011-section54-pressure-selection-{reviewer_id}"
        )
        _require(
            Path(binding["path"]).parent.name == expected_directory,
            "pressure reviewer attestation directory name differs",
        )
        _require(
            _same_json_value(
                _published_binding(
                    attestation["published_pressure_pilot_receipt"],
                    root,
                    "pressure reviewer attestation.published_pressure_pilot_receipt",
                ),
                aggregate_binding,
            )
            and _same_json_value(
                _published_binding(
                    attestation["published_pressure_pilot_review_packet_receipt"],
                    root,
                    "pressure reviewer attestation.published_pressure_pilot_review_packet_receipt",
                ),
                packet_binding,
            )
            and _same_json_value(
                _pressure_gate_binding(
                    attestation["authoritative_reanalysis_attestation"],
                    authorized_pic_root=root,
                    label="pressure reviewer attestation.authoritative_reanalysis_attestation",
                ),
                reanalysis_binding,
            )
            and _same_json_value(attestation["selected_case"], normalized_case),
            "pressure reviewer attestation selection/evidence tuple drifted",
        )
        return {
            "binding": copy.deepcopy(binding),
            "attestation": copy.deepcopy(attestation),
            "reviewer_id": reviewer_id,
            "reviewed_utc": reviewed_text,
            "sealed_utc": sealed_text,
            "rationale": rationale,
            "selected_case": copy.deepcopy(normalized_case),
            "authoritative_reanalysis_attestation": copy.deepcopy(reanalysis_binding),
        }
    except PressureReviewPacketVerificationError:
        raise
    except (OSError, RecursionError, TypeError, ValueError) as exc:
        raise PressureReviewPacketVerificationError(
            f"pressure reviewer attestation verification failed closed: {exc}"
        ) from exc


def validate_pressure_reanalysis_source_snapshot(
    reanalysis_verification: object,
    *,
    git_commit: object,
    source_archive_sha256: object,
    helper_source_closure: object,
    authorized_pic_root: str | os.PathLike[str],
) -> None:
    """Bind accepted reanalysis code to the frozen qualifying-plan source archive."""
    reanalysis = _reconsume_sealed_pressure_reanalysis_verification(
        reanalysis_verification,
        authorized_pic_root=authorized_pic_root,
    )
    authorization = _validate_reanalysis_source_authorization(
        reanalysis["source_authorization"]
    )
    _require(
        type(git_commit) is str and _GIT_COMMIT_RE.fullmatch(git_commit) is not None,
        "qualifying-plan Git commit is malformed",
    )
    archive_sha256 = _require_sha256(
        source_archive_sha256,
        "qualifying-plan source archive SHA-256",
    )
    _require(
        type(helper_source_closure) is list,
        "qualifying-plan helper-source closure must be a list",
    )
    helper_by_path: dict[str, str] = {}
    for index, record in enumerate(helper_source_closure):
        label = f"qualifying-plan helper-source closure[{index}]"
        _require_exact_keys(record, {"path", "sha256"}, label)
        path = _relative_source_path(record["path"], f"{label}.path")
        digest = _require_sha256(record["sha256"], f"{label}.sha256")
        _require(path not in helper_by_path, "qualifying-plan helper-source closure has duplicate paths")
        helper_by_path[path] = digest
    expected_closure = [
        {"path": path, "sha256": helper_by_path.get(path)}
        for path in PRESSURE_REANALYSIS_SOURCE_PATHS
    ]
    _require(
        all(record["sha256"] is not None for record in expected_closure),
        "qualifying-plan helper-source closure omits a reanalysis source",
    )
    _require(
        authorization["git_commit"] == git_commit
        and authorization["source_archive_sha256"] == archive_sha256
        and _same_json_value(authorization["source_closure"], expected_closure),
        "authoritative reanalysis source authorization differs from the frozen qualifying-plan source snapshot",
    )


__all__ = [
    "AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION",
    "PRESSURE_GATE_ATTESTATION_FILENAME",
    "PRESSURE_GATE_ATTESTATION_ROOT_NAME",
    "PRESSURE_GATE_REVIEW_NOT_BEFORE_UTC",
    "PRESSURE_REANALYSIS_EXECUTION_MODE",
    "PRESSURE_REANALYSIS_OPERATOR_STATEMENT",
    "PRESSURE_REANALYSIS_QUALIFICATION_EFFECT",
    "PRESSURE_REANALYSIS_RECORD_TYPE",
    "PRESSURE_REANALYSIS_SOURCE_PATHS",
    "PRESSURE_REVIEWER_QUALIFICATION_EFFECT",
    "PRESSURE_REVIEWER_RECORD_TYPE",
    "PRESSURE_REVIEWER_STATEMENT",
    "PressureReviewPacketVerificationError",
    "consume_sealed_pressure_reanalysis_attestation",
    "consume_sealed_pressure_reviewer_attestation",
    "consume_published_pressure_pilot_review_packet",
    "validate_pressure_reanalysis_source_snapshot",
    "verify_pressure_review_packet_binding",
]
