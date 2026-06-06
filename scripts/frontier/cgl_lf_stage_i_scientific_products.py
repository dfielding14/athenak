#!/usr/bin/env python3
"""Generate replayable Stage I scientific products from authenticated raw data.

This utility intentionally does not import ``scripts/analyze_cgl_lf_paper.py``.
It authenticates an accepted whole-case bundle, derives supported products
directly from retained histories and rank-local AthenaK snapshots, and emits a
self-digested evidence record.  Verification ignores claimed products,
recomputes the complete record from its bound request, and byte-compares it
with the retained evidence.

The implementation supports the products needed for scalar scientific-
acceptance inputs, density/anisotropy PDFs, pressure-density surfaces,
pressure-anisotropy transfer, local-field eddy anisotropy, selected-shell
alignment PDFs, peak-alignment curves, and R16/R02/R17 velocity/magnetic
spectral-shape convergence.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import stat
import sys
from typing import Iterable

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[2]
BIN_CONVERT_PATH = REPOSITORY / "vis/python/bin_convert.py"
BIN_CONVERT_SPEC = importlib.util.spec_from_file_location(
    "_cgl_lf_stage_i_bin_convert", BIN_CONVERT_PATH
)
if BIN_CONVERT_SPEC is None or BIN_CONVERT_SPEC.loader is None:
    raise RuntimeError(f"cannot load exact Athena binary parser: {BIN_CONVERT_PATH}")
bin_convert = importlib.util.module_from_spec(BIN_CONVERT_SPEC)
BIN_CONVERT_SPEC.loader.exec_module(bin_convert)


CANONICAL_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
SCHEMA_VERSION = 1
RECORD_TYPE = "stage-i-deterministic-scientific-products"
ACCEPTANCE_CONTRACT_SCHEMA_VERSION = 2
ACCEPTANCE_CONTRACT_RECORD_TYPE = "stage-i-scientific-acceptance-reviewed-products"
DETERMINISTIC_REPLAY_VERIFICATION = (
    "exact-semantic-replay-from-bound-canonical-case-inputs"
)
EVIDENCE_DIGEST_METHOD = "sha256-canonical-json-without-evidence-digest"
HISTORY_LABEL = re.compile(r"\[(\d+)\]=(\S+)")
SHA256 = re.compile(r"[0-9a-f]{64}")
CASE_ID = re.compile(r"R(?:0[2-9]|1[0-7])")
RANK_DIRECTORY = re.compile(r"rank_(\d{8})")
EDDY_ANGLE_DEGREES = 15.0
EDDY_SAMPLES = 2_000_000
EDDY_BINS = 24
EDDY_SEED = 731
REQUIRED_FIELDS = (
    "dens",
    "velx",
    "vely",
    "velz",
    "eint",
    "p_perp",
    "bcc1",
    "bcc2",
    "bcc3",
)
HISTORY_METRICS = {
    "abs_dp": ("user", "abs_dp"),
    "beta": ("user", "beta"),
    "firehose_occupancy": ("user", "fire_vol"),
    "force_power": ("user", "force_pwr"),
    "hard_occupancy": ("user", "hard_vol"),
    "kinetic": ("user", "kinetic"),
    "magnetic": ("user", "magnetic"),
    "mirror_occupancy": ("user", "mirror_vol"),
    "nu_eff": ("user", "nu_eff"),
}
LF_COUNTERS = (
    "lf_hwproj",
    "lf_cpwrk",
    "lf_cawrk",
    "lf_qface",
    "lf_qprwrk",
    "lf_qpewrk",
)
SUPPORTED_REFERENCE_CURVES = {
    "history.unstable_fraction",
    "pdf.density_fluctuation",
    "pdf.beta_delta",
    "alignment_peak.cos_theta",
    "pressure_transfer.transfer",
    "pressure_transfer.transfer_normalized_by_total",
    "eddy_anisotropy.velocity_perp",
    "eddy_anisotropy.magnetic_perp",
}
SUPPORTED_REFERENCE_SURFACES = {
    "pressure_density_joint.parallel",
    "pressure_density_joint.perpendicular",
}
class ScientificProductsError(RuntimeError):
    """Raised when authentication or deterministic replay fails."""


class UnsupportedProduct(ScientificProductsError):
    """Raised when an authenticated input uses an unsupported product format."""


@dataclass(frozen=True)
class SnapshotGroup:
    """One exact rank-local snapshot inventory bound by segment inspection."""

    time: float
    representative: Path
    rank_files: tuple[dict[str, object], ...]
    segment: str


@dataclass
class BundleContext:
    """Authenticated whole-case inputs used for one deterministic generation."""

    bundle: dict[str, object]
    bundle_binding: dict[str, object]
    stage_manifest: dict[str, object]
    stage_manifest_binding: dict[str, object]
    segment_bindings: list[dict[str, object]]
    segment_history_bindings: list[dict[str, object]]
    case_id: str
    case_name: str
    model_choices: dict[str, object]
    mhd: dict[str, list[float]]
    mhd_binding: dict[str, object]
    user: dict[str, list[float]]
    user_binding: dict[str, object]
    snapshots: list[SnapshotGroup]
    reference_root: Path
    authority_scope: str


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON keys."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ScientificProductsError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> object:
    """Reject non-finite JSON constants."""

    raise ScientificProductsError(f"invalid JSON numeric constant: {value}")


def canonical_json(value: object) -> bytes:
    """Return compact deterministic JSON bytes."""

    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def stable_json(value: object) -> bytes:
    """Return human-readable deterministic JSON bytes."""

    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    """Return a lowercase SHA-256 digest."""

    return hashlib.sha256(value).hexdigest()


def require_dict(value: object, label: str) -> dict[str, object]:
    """Require one JSON object."""

    if not isinstance(value, dict):
        raise ScientificProductsError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    """Require one JSON list."""

    if not isinstance(value, list):
        raise ScientificProductsError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    """Require one nonempty string."""

    if not isinstance(value, str) or not value:
        raise ScientificProductsError(f"{label} must be a nonempty string")
    return value


def require_sha256(value: object, label: str) -> str:
    """Require one lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise ScientificProductsError(f"{label} must be a lowercase SHA-256")
    return value


def require_finite(value: object, label: str) -> float:
    """Require one finite non-boolean number."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScientificProductsError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ScientificProductsError(f"{label} must be finite")
    return result


def require_int(value: object, label: str, minimum: int = 0) -> int:
    """Require one exact integer."""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ScientificProductsError(f"{label} must be an integer >= {minimum}")
    return value


def stable_profile(value: os.stat_result) -> tuple[int, ...]:
    """Return identity and mutation metadata for descriptor-bound reads."""

    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def descriptor_sha256(descriptor: int) -> str:
    """Hash one open descriptor without changing its offset."""

    digest = hashlib.sha256()
    offset = 0
    while True:
        block = os.pread(descriptor, 1024 * 1024, offset)
        if not block:
            break
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def regular_file_binding(
    path: Path,
    label: str,
    *,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
) -> dict[str, object]:
    """Bind one unchanged regular file without following its leaf symlink."""

    absolute = path.expanduser().absolute()
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(absolute, flags)
    except OSError as error:
        raise ScientificProductsError(f"{label} cannot be opened safely: {absolute}") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ScientificProductsError(f"{label} is not a regular file: {absolute}")
        digest = descriptor_sha256(descriptor)
        after = os.fstat(descriptor)
        if stable_profile(before) != stable_profile(after):
            raise ScientificProductsError(f"{label} changed while it was hashed")
        if expected_sha256 is not None and digest != expected_sha256:
            raise ScientificProductsError(f"{label} SHA-256 differs from its authority")
        if expected_size is not None and before.st_size != expected_size:
            raise ScientificProductsError(f"{label} size differs from its authority")
        try:
            resolved = absolute.resolve(strict=True)
        except OSError as error:
            raise ScientificProductsError(f"{label} path changed while it was hashed") from error
        return {
            "path": str(resolved),
            "size_bytes": before.st_size,
            "sha256": digest,
        }
    finally:
        os.close(descriptor)


def read_regular_bytes(
    path: Path,
    label: str,
    *,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
    maximum: int = 64 * 1024 * 1024,
) -> tuple[bytes, dict[str, object]]:
    """Read one bounded unchanged regular file and return its binding."""

    binding = regular_file_binding(
        path,
        label,
        expected_sha256=expected_sha256,
        expected_size=expected_size,
    )
    size = require_int(binding["size_bytes"], f"{label} size")
    if size > maximum:
        raise ScientificProductsError(f"{label} exceeds the reviewed size limit")
    try:
        payload = Path(str(binding["path"])).read_bytes()
    except OSError as error:
        raise ScientificProductsError(f"{label} changed between binding and read") from error
    if len(payload) != size or sha256_bytes(payload) != binding["sha256"]:
        raise ScientificProductsError(f"{label} changed between binding and read")
    return payload, binding


def load_json(
    path: Path, label: str, *, expected_sha256: str | None = None
) -> tuple[dict[str, object], dict[str, object]]:
    """Load one unambiguous bound JSON object."""

    payload, binding = read_regular_bytes(path, label, expected_sha256=expected_sha256)
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ScientificProductsError(f"{label} is invalid JSON") from error
    return require_dict(value, label), binding


def evidence_digest(value: dict[str, object]) -> str:
    """Digest an evidence object without its self-referential member."""

    body = dict(value)
    body.pop("evidence_digest", None)
    return sha256_bytes(canonical_json(body))


def seal_evidence(value: dict[str, object]) -> dict[str, object]:
    """Attach the deterministic evidence self-digest."""

    if "evidence_digest" in value:
        raise ScientificProductsError("evidence is already sealed")
    result = dict(value)
    result["evidence_digest"] = {
        "method": EVIDENCE_DIGEST_METHOD,
        "sha256": evidence_digest(result),
    }
    return result


def verify_self_digest(value: dict[str, object]) -> None:
    """Verify one evidence self-digest."""

    record = require_dict(value.get("evidence_digest"), "evidence_digest")
    if record.get("method") != EVIDENCE_DIGEST_METHOD:
        raise ScientificProductsError("evidence digest method is unsupported")
    expected = require_sha256(record.get("sha256"), "evidence digest")
    if evidence_digest(value) != expected:
        raise ScientificProductsError("evidence self-digest differs")


def write_candidate(path: Path, value: dict[str, object]) -> None:
    """Write one immutable, no-clobber candidate."""

    absolute = path.expanduser().absolute()
    try:
        absolute.resolve(strict=False).relative_to(CANONICAL_ROOT.resolve(strict=False))
    except ValueError:
        pass
    else:
        raise ScientificProductsError("candidate output beneath canonical root is forbidden")
    absolute.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    try:
        descriptor = os.open(absolute, flags, 0o444)
    except FileExistsError as error:
        raise ScientificProductsError(f"candidate already exists: {absolute}") from error
    try:
        payload = stable_json(value)
        offset = 0
        while offset < len(payload):
            offset += os.write(descriptor, payload[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def parse_history(payload: bytes, label: str) -> dict[str, list[float]]:
    """Parse one finite, strictly ordered Athena history."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ScientificProductsError(f"{label} is not UTF-8") from error
    labels: list[str] | None = None
    rows: list[list[float]] = []
    for line in text.splitlines():
        if line.startswith("#"):
            found = HISTORY_LABEL.findall(line)
            if found:
                indexed = sorted((int(index), name) for index, name in found)
                indices = [index for index, _ in indexed]
                names = [name for _, name in indexed]
                if (
                    not indices
                    or indices[0] not in (0, 1)
                    or indices != list(range(indices[0], indices[0] + len(indices)))
                    or len(names) != len(set(names))
                ):
                    raise ScientificProductsError(f"{label} labels are malformed")
                if labels is not None and labels != names:
                    raise ScientificProductsError(f"{label} has conflicting headers")
                labels = names
            continue
        if not line.strip():
            continue
        try:
            row = [float(value) for value in line.split()]
        except ValueError as error:
            raise ScientificProductsError(f"{label} contains a nonnumeric row") from error
        if not all(math.isfinite(value) for value in row):
            raise ScientificProductsError(f"{label} contains a nonfinite row")
        rows.append(row)
    if labels is None or len(rows) < 2:
        raise ScientificProductsError(f"{label} lacks labels or sufficient rows")
    if any(len(row) != len(labels) for row in rows):
        raise ScientificProductsError(f"{label} row width differs from its header")
    result = {
        name: [row[index] for row in rows] for index, name in enumerate(labels)
    }
    times = result.get("time")
    if times is None or any(right <= left for left, right in zip(times, times[1:])):
        raise ScientificProductsError(f"{label} time is not strictly increasing")
    return result


def load_history(
    path: Path,
    label: str,
    *,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
) -> tuple[dict[str, list[float]], dict[str, object]]:
    """Load and bind one retained Athena history."""

    payload, binding = read_regular_bytes(
        path,
        label,
        expected_sha256=expected_sha256,
        expected_size=expected_size,
    )
    return parse_history(payload, label), binding


def merge_history_records(
    records: list[dict[str, list[float]]], label: str
) -> dict[str, list[float]]:
    """Reproduce the accepted restart-history merge on parsed values."""

    if not records or "time" not in records[0]:
        raise ScientificProductsError(f"{label} segment histories are empty")
    names = list(records[0])
    merged = {name: [] for name in names}
    last_time = -math.inf
    for record in records:
        if list(record) != names:
            raise ScientificProductsError(f"{label} segment history columns differ")
        lengths = {len(values) for values in record.values()}
        if len(lengths) != 1:
            raise ScientificProductsError(f"{label} segment history lengths differ")
        for index, time in enumerate(record["time"]):
            if time > last_time + 1.0e-12:
                for name in names:
                    merged[name].append(record[name][index])
                last_time = time
    if len(merged["time"]) < 2:
        raise ScientificProductsError(f"{label} merged history is insufficient")
    return merged


def interpolate_at(times: list[float], values: list[float], target: float) -> float:
    """Linearly interpolate one ordered finite series."""

    if target < times[0] or target > times[-1]:
        raise ScientificProductsError("analysis window lies outside history coverage")
    for index, time in enumerate(times):
        if time == target:
            return values[index]
        if time > target:
            fraction = (target - times[index - 1]) / (time - times[index - 1])
            return values[index - 1] + fraction * (values[index] - values[index - 1])
    return values[-1]


def clipped_series(
    times: list[float], values: list[float], start: float, end: float
) -> tuple[list[float], list[float]]:
    """Clip one series to exact endpoints."""

    if len(times) != len(values) or len(times) < 2 or start >= end:
        raise ScientificProductsError("invalid history series or analysis window")
    clipped_times = [start]
    clipped_values = [interpolate_at(times, values, start)]
    for time, value in zip(times, values):
        if start < time < end:
            clipped_times.append(time)
            clipped_values.append(value)
    clipped_times.append(end)
    clipped_values.append(interpolate_at(times, values, end))
    return clipped_times, clipped_values


def trapezoidal_weights(times: list[float]) -> list[float]:
    """Return normalized endpoint-clipped trapezoidal weights."""

    result = [0.0] * len(times)
    for index, (left, right) in enumerate(zip(times, times[1:])):
        width = right - left
        if width <= 0.0 or not math.isfinite(width):
            raise ScientificProductsError("time grid is not strictly increasing")
        result[index] += 0.5 * width
        result[index + 1] += 0.5 * width
    total = sum(result)
    return [value / total for value in result]


def metric_record(
    times: list[float], values: list[float], start: float, end: float
) -> dict[str, object]:
    """Return raw samples and deterministic time-weighted moments."""

    selected_times, selected_values = clipped_series(times, values, start, end)
    weights = trapezoidal_weights(selected_times)
    mean = sum(weight * value for weight, value in zip(weights, selected_values))
    variance = sum(
        weight * (value - mean) ** 2
        for weight, value in zip(weights, selected_values)
    )
    return {
        "sample_times": selected_times,
        "sample_values": selected_values,
        "time_weighted_mean": mean,
        "time_weighted_standard_deviation": math.sqrt(max(0.0, variance)),
        "uncertainty_status": "deferred_to_reviewed_scientific_acceptance_policy",
        "method": "endpoint-clipped-trapezoidal",
    }


def path_within(path: Path, root: Path, label: str) -> Path:
    """Resolve a path and require containment beneath one root."""

    try:
        resolved = path.expanduser().absolute().resolve(strict=True)
        resolved_root = root.expanduser().absolute().resolve(strict=True)
    except OSError as error:
        raise ScientificProductsError(f"{label} path or authenticated root is missing") from error
    try:
        resolved.relative_to(resolved_root)
    except ValueError as error:
        raise ScientificProductsError(f"{label} escapes its authenticated root") from error
    return resolved


def absolute_existing_path(value: object, label: str) -> Path:
    """Resolve one declared absolute path with a controlled failure."""

    path = Path(require_text(value, label)).expanduser()
    if not path.is_absolute():
        raise ScientificProductsError(f"{label} must be absolute")
    try:
        return path.resolve(strict=True)
    except OSError as error:
        raise ScientificProductsError(f"{label} does not exist: {path}") from error


def bundle_output_path(root: Path, value: object, label: str) -> Path:
    """Resolve one bundle-relative output path."""

    text = require_text(value, label)
    path = Path(text)
    if path.is_absolute():
        raise ScientificProductsError(f"{label} must be bundle-relative")
    return path_within(root / path, root, label)


def bundle_snapshot_target(root: Path, value: object, label: str) -> Path:
    """Resolve a bundle-local snapshot entry whose leaf may link outside."""

    text = require_text(value, label)
    relative = Path(text)
    if relative.is_absolute() or ".." in relative.parts:
        raise ScientificProductsError(f"{label} must be a contained bundle-relative path")
    root_resolved = root.expanduser().absolute().resolve(strict=True)
    entry = root_resolved / relative
    path_within(entry.parent, root_resolved, f"{label} parent")
    try:
        mode = entry.lstat().st_mode
    except OSError as error:
        raise ScientificProductsError(f"{label} does not exist: {entry}") from error
    if not (stat.S_ISLNK(mode) or stat.S_ISREG(mode)):
        raise ScientificProductsError(f"{label} is not a regular file or symlink: {entry}")
    try:
        resolved = entry.resolve(strict=True)
    except OSError as error:
        raise ScientificProductsError(f"{label} target does not exist: {entry}") from error
    if not stat.S_ISREG(resolved.stat().st_mode):
        raise ScientificProductsError(f"{label} target is not a regular file: {resolved}")
    return resolved


def declared_file_binding(value: object, label: str) -> dict[str, object]:
    """Validate one path/SHA/size declaration."""

    record = require_dict(value, label)
    path = Path(require_text(record.get("path"), f"{label} path")).expanduser()
    if not path.is_absolute():
        raise ScientificProductsError(f"{label} path must be absolute")
    return {
        "path": str(path.absolute()),
        "size_bytes": require_int(record.get("size_bytes"), f"{label} size", 1),
        "sha256": require_sha256(record.get("sha256"), f"{label} sha256"),
    }


def validate_snapshot_record(
    record: object, time: float, segment: str
) -> SnapshotGroup:
    """Validate one exact rank-local snapshot declaration."""

    value = require_dict(record, "snapshot record")
    rank_values = require_list(value.get("rank_files"), "snapshot rank_files")
    if value.get("storage") != "per_rank" or not rank_values:
        raise UnsupportedProduct("only exact rank-local snapshots are supported")
    rank_files = tuple(
        declared_file_binding(item, f"snapshot rank file {index}")
        for index, item in enumerate(rank_values)
    )
    paths = [Path(str(item["path"])) for item in rank_files]
    rank_ids: list[int] = []
    for path in paths:
        match = RANK_DIRECTORY.fullmatch(path.parent.name)
        if match is None:
            raise ScientificProductsError("snapshot rank file lacks rank_######## parent")
        rank_ids.append(int(match.group(1)))
    if rank_ids != list(range(len(paths))) or len(set(paths)) != len(paths):
        raise ScientificProductsError("snapshot rank inventory is not exact and contiguous")
    representative = Path(require_text(value.get("path"), "snapshot representative"))
    if representative.expanduser().absolute() != paths[0].expanduser().absolute():
        raise ScientificProductsError("snapshot representative is not rank zero")
    if value.get("sha256") != rank_files[0]["sha256"]:
        raise ScientificProductsError("snapshot representative digest differs from rank zero")
    return SnapshotGroup(
        time=time,
        representative=representative.expanduser().absolute(),
        rank_files=rank_files,
        segment=segment,
    )


def distinct_snapshot_groups(groups: list[SnapshotGroup]) -> list[SnapshotGroup]:
    """Reproduce whole-case bundle ordering and physical-time deduplication."""

    ordered = sorted(
        groups,
        key=lambda group: (
            group.time,
            tuple(str(record["path"]) for record in group.rank_files),
        ),
    )
    distinct: list[SnapshotGroup] = []
    last_time = -math.inf
    for group in ordered:
        if group.time > last_time + 1.0e-12:
            distinct.append(group)
            last_time = group.time
    return distinct


def authenticate_bundle(request: dict[str, object]) -> BundleContext:
    """Authenticate one accepted whole-case bundle and selected raw products."""

    bundle_path = Path(require_text(request.get("bundle_manifest"), "bundle_manifest"))
    expected_bundle = require_sha256(
        request.get("expected_bundle_sha256"), "expected_bundle_sha256"
    )
    bundle, bundle_binding = load_json(
        bundle_path, "whole-case bundle manifest", expected_sha256=expected_bundle
    )
    authority_mode = request.get("authority_mode")
    if authority_mode not in ("canonical", "offline"):
        raise ScientificProductsError("authority_mode must be canonical or offline")
    bundle_resolved = Path(str(bundle_binding["path"]))
    canonical_stage_root: Path | None = None
    if authority_mode == "canonical":
        canonical_stage_root = (
            CANONICAL_ROOT
            / "runs/mks24-stage-i/E03-forcing-policy"
        ).resolve(strict=True)
        path_within(bundle_resolved, canonical_stage_root / "bundles", "canonical bundle")
        authority_scope = "canonical-path-and-external-bundle-digest"
    else:
        authority_scope = "offline-replay-external-bundle-digest-non-authorizing"

    if (
        bundle.get("workflow") != "paper-mks24-stage-i-production"
        or bundle.get("status") != "accepted_for_analysis"
        or require_finite(bundle.get("required_final_time"), "required_final_time") != 10.0
        or require_finite(bundle.get("accepted_final_time"), "accepted_final_time") != 10.0
    ):
        raise ScientificProductsError("whole-case bundle identity or final time differs")
    case_id = require_text(bundle.get("production_case_id"), "production_case_id")
    if CASE_ID.fullmatch(case_id) is None:
        raise ScientificProductsError("production_case_id is outside R02-R17")
    cases = require_list(bundle.get("cases"), "bundle cases")
    if len(cases) != 1:
        raise ScientificProductsError("whole-case bundle must contain exactly one case")
    case = require_dict(cases[0], "bundle case")
    case_name = require_text(case.get("name"), "bundle case name")
    if case.get("status") != "passed":
        raise ScientificProductsError("bundle case is not passed")
    model_choices = require_dict(case.get("model_choices"), "model_choices")

    stage_path = absolute_existing_path(bundle.get("stage_i_manifest"), "stage_i_manifest")
    stage_manifest, stage_binding = load_json(stage_path, "Stage I manifest")
    stage_cases = {
        require_text(require_dict(item, "Stage I case").get("id"), "Stage I case id"):
        require_dict(item, "Stage I case")
        for item in require_list(stage_manifest.get("cases"), "Stage I cases")
    }
    selected = stage_cases.get(case_id)
    if (
        selected is None
        or selected.get("name") != case_name
        or selected.get("input") != case.get("input")
    ):
        raise ScientificProductsError("bundle case differs from the Stage I manifest")

    segment_paths = [
        absolute_existing_path(item, "production segment manifest")
        for item in require_list(
            bundle.get("production_segment_manifests"),
            "production_segment_manifests",
        )
    ]
    if not segment_paths or len(set(segment_paths)) != len(segment_paths):
        raise ScientificProductsError("accepted segment lineage is empty or duplicated")
    segment_bindings: list[dict[str, object]] = []
    segment_history_bindings: list[dict[str, object]] = []
    segment_histories: dict[str, list[dict[str, list[float]]]] = {
        "mhd_history": [],
        "user_history": [],
    }
    snapshots: list[SnapshotGroup] = []
    previous_path: Path | None = None
    previous_final = -math.inf
    executable_sha: str | None = None
    input_sha: str | None = None
    for index, segment_path in enumerate(segment_paths):
        if canonical_stage_root is not None:
            path_within(
                segment_path,
                canonical_stage_root,
                "canonical segment manifest",
            )
        segment, binding = load_json(segment_path, f"accepted segment {index}")
        segment_bindings.append(binding)
        accounting = require_dict(segment.get("accounting"), "segment accounting")
        command = require_dict(segment.get("command"), "segment command")
        inspection = require_dict(
            segment.get("scientific_inspection"), "segment scientific_inspection"
        )
        segment_name = require_text(accounting.get("segment"), "segment name")
        final_time = require_finite(inspection.get("final_time"), "segment final_time")
        segment_executable = require_sha256(
            command.get("executable_sha256"), "segment executable_sha256"
        )
        segment_input = require_sha256(command.get("input_sha256"), "segment input_sha256")
        if (
            accounting.get("result") != "accepted"
            or accounting.get("case_id") != case_id
            or accounting.get("case_name") != case_name
            or accounting.get("executable_sha256") != segment_executable
            or inspection.get("accepted") is not True
            or inspection.get("case_id") != case_id
            or Path(require_text(inspection.get("manifest"), "inspection manifest")).resolve(strict=True)
            != Path(str(binding["path"]))
            or command.get("matrix_sha256") != stage_binding["sha256"]
            or final_time <= previous_final
        ):
            raise ScientificProductsError("accepted segment lineage semantics differ")
        parent = command.get("parent_segment")
        if previous_path is None:
            if (
                parent is not None
                or not segment_name.startswith("s00")
                or command.get("restart_file") is not None
                or require_list(command.get("restart_files"), "root restart_files")
            ):
                raise ScientificProductsError("accepted lineage root is not an s00 origin")
        else:
            parent_record = require_dict(parent, "segment parent")
            if (
                Path(require_text(parent_record.get("manifest"), "parent manifest")).resolve(strict=True)
                != previous_path
                or parent_record.get("case_id") != case_id
                or parent_record.get("result") != "accepted"
                or require_finite(parent_record.get("final_time"), "parent final_time")
                != previous_final
                or parent_record.get("executable_sha256") != segment_executable
                or parent_record.get("input_sha256") != segment_input
            ):
                raise ScientificProductsError("accepted segment parent lineage differs")
        if executable_sha is None:
            executable_sha = segment_executable
            input_sha = segment_input
        elif executable_sha != segment_executable or input_sha != segment_input:
            raise ScientificProductsError("segment executable or input lineage differs")
        current_histories: dict[str, dict[str, list[float]]] = {}
        for key in ("mhd_history", "user_history"):
            declared = declared_file_binding(
                inspection.get(key), f"segment {index} {key}"
            )
            history_path = Path(str(declared["path"]))
            if canonical_stage_root is not None:
                path_within(
                    history_path,
                    canonical_stage_root,
                    f"canonical segment {key}",
                )
            history, history_binding = load_history(
                history_path,
                f"segment {index} {key}",
                expected_sha256=str(declared["sha256"]),
                expected_size=int(declared["size_bytes"]),
            )
            current_histories[key] = history
            segment_histories[key].append(history)
            segment_history_bindings.append(history_binding)
        if (
            current_histories["mhd_history"]["time"]
            != current_histories["user_history"]["time"]
            or not math.isclose(
                current_histories["mhd_history"]["time"][-1],
                final_time,
                rel_tol=0.0,
                abs_tol=1.0e-10,
            )
        ):
            raise ScientificProductsError(
                "accepted segment histories are unsynchronized or differ from final time"
            )
        times = [
            require_finite(item, "snapshot time")
            for item in require_list(inspection.get("snapshot_times"), "snapshot_times")
        ]
        records = require_list(inspection.get("snapshots"), "snapshots")
        if len(times) != len(records):
            raise ScientificProductsError("snapshot time and record counts differ")
        groups = [
            validate_snapshot_record(record, time, segment_name)
            for time, record in zip(times, records)
        ]
        if canonical_stage_root is not None:
            for group in groups:
                for rank_file in group.rank_files:
                    path_within(
                        Path(str(rank_file["path"])),
                        canonical_stage_root,
                        "canonical snapshot rank file",
                    )
        snapshots.extend(groups)
        previous_path = Path(str(binding["path"]))
        previous_final = final_time
    if previous_final != 10.0:
        raise ScientificProductsError("accepted lineage does not terminate at t=10")
    snapshots = distinct_snapshot_groups(snapshots)

    outputs = require_dict(case.get("outputs"), "bundle case outputs")
    bundle_root = bundle_resolved.parent
    mhd_path = bundle_output_path(bundle_root, outputs.get("mhd_history"), "mhd_history")
    user_path = bundle_output_path(bundle_root, outputs.get("user_history"), "user_history")
    mhd, mhd_binding = load_history(mhd_path, "whole-case MHD history")
    user, user_binding = load_history(user_path, "whole-case user history")
    if mhd != merge_history_records(segment_histories["mhd_history"], "MHD"):
        raise ScientificProductsError(
            "whole-case MHD history differs from authenticated segment merge"
        )
    if user != merge_history_records(segment_histories["user_history"], "user"):
        raise ScientificProductsError(
            "whole-case user history differs from authenticated segment merge"
        )

    snapshot_paths = [
        bundle_snapshot_target(bundle_root, item, "bundle snapshot")
        for item in require_list(outputs.get("snapshot_paths", []), "snapshot_paths")
    ]
    bundle_targets = snapshot_paths
    inspection_targets = [group.representative.resolve(strict=True) for group in snapshots]
    if bundle_targets != inspection_targets:
        raise ScientificProductsError(
            "bundle snapshot inventory differs from accepted segment inspections"
        )

    reference_root_text = request.get("reference_root")
    if reference_root_text is None:
        provenance = require_dict(bundle.get("reference_provenance"), "reference_provenance")
        reference_root = Path(
            require_text(provenance.get("reference_manifest"), "reference manifest")
        ).resolve(strict=True).parent
    else:
        reference_root = Path(require_text(reference_root_text, "reference_root")).resolve(strict=True)

    return BundleContext(
        bundle=bundle,
        bundle_binding=bundle_binding,
        stage_manifest=stage_manifest,
        stage_manifest_binding=stage_binding,
        segment_bindings=segment_bindings,
        segment_history_bindings=segment_history_bindings,
        case_id=case_id,
        case_name=case_name,
        model_choices=model_choices,
        mhd=mhd,
        mhd_binding=mhd_binding,
        user=user,
        user_binding=user_binding,
        snapshots=sorted(snapshots, key=lambda value: value.time),
        reference_root=reference_root,
        authority_scope=authority_scope,
    )


def metadata_equal(left: object, right: object) -> bool:
    """Return exact equality for binary parser metadata."""

    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        return bool(np.array_equal(np.asarray(left), np.asarray(right)))
    return left == right


def read_snapshot_group(
    group: SnapshotGroup, max_cells: int
) -> tuple[dict[str, np.ndarray], tuple[float, float, float], list[dict[str, object]]]:
    """Authenticate and reconstruct one uniform, non-AMR, 3-D rank set."""

    payloads: list[dict[str, object]] = []
    bindings: list[dict[str, object]] = []
    for rank, declared in enumerate(group.rank_files):
        path = Path(str(declared["path"]))
        before = regular_file_binding(
            path,
            f"snapshot {group.time:g} rank {rank}",
            expected_sha256=str(declared["sha256"]),
            expected_size=int(declared["size_bytes"]),
        )
        raw = bin_convert.read_binary(str(path))
        after = regular_file_binding(
            path,
            f"snapshot {group.time:g} rank {rank} after parse",
            expected_sha256=str(declared["sha256"]),
            expected_size=int(declared["size_bytes"]),
        )
        if before != after:
            raise ScientificProductsError("snapshot changed while it was parsed")
        payloads.append(raw)
        bindings.append(before)
    reference = payloads[0]
    metadata = (
        "header",
        "time",
        "cycle",
        "var_names",
        "Nx1",
        "Nx2",
        "Nx3",
        "nvars",
        "x1min",
        "x1max",
        "x2min",
        "x2max",
        "x3min",
        "x3max",
        "nx1_mb",
        "nx2_mb",
        "nx3_mb",
        "nx1_out_mb",
        "nx2_out_mb",
        "nx3_out_mb",
    )
    for rank, raw in enumerate(payloads[1:], start=1):
        for key in metadata:
            if key not in raw or key not in reference or not metadata_equal(raw[key], reference[key]):
                raise UnsupportedProduct(
                    f"rank-local snapshot metadata differs at rank={rank}, key={key}"
                )
    missing = sorted(set(REQUIRED_FIELDS) - set(reference["var_names"]))
    if missing:
        raise UnsupportedProduct(f"snapshot lacks required fields: {missing}")
    if not math.isclose(float(reference["time"]), group.time, rel_tol=0.0, abs_tol=1.0e-12):
        raise ScientificProductsError("snapshot binary time differs from inspection")
    root = tuple(int(reference[f"Nx{axis}"]) for axis in (1, 2, 3))
    block = tuple(int(reference[f"nx{axis}_mb"]) for axis in (1, 2, 3))
    output = tuple(int(reference[f"nx{axis}_out_mb"]) for axis in (1, 2, 3))
    if (
        any(value <= 1 for value in root)
        or block != output
        or any(value <= 0 for value in block)
        or any(size % width for size, width in zip(root, block))
    ):
        raise UnsupportedProduct("only full uniform non-AMR 3-D snapshots are supported")
    cell_count = math.prod(root)
    if cell_count > max_cells:
        raise UnsupportedProduct(
            f"snapshot has {cell_count} cells, above max_snapshot_cells={max_cells}"
        )
    logical = np.concatenate(
        [np.asarray(raw["mb_logical"], dtype=np.int64) for raw in payloads], axis=0
    )
    if logical.ndim != 2 or logical.shape[1] != 4 or np.any(logical[:, 3] != 0):
        raise UnsupportedProduct("snapshot is AMR or has malformed logical locations")
    locations = [tuple(int(value) for value in row) for row in logical[:, :3]]
    counts = tuple(size // width for size, width in zip(root, block))
    expected = {
        (i, j, k)
        for k in range(counts[2])
        for j in range(counts[1])
        for i in range(counts[0])
    }
    if set(locations) != expected or len(locations) != len(expected):
        raise UnsupportedProduct("snapshot meshblock coverage is incomplete or duplicated")
    fields = {
        name: np.empty((root[2], root[1], root[0]), dtype=np.float64)
        for name in REQUIRED_FIELDS
    }
    row = 0
    for raw in payloads:
        raw_logical = np.asarray(raw["mb_logical"], dtype=np.int64)
        for local, location in enumerate(raw_logical[:, :3]):
            i, j, k = (int(value) for value in location)
            slices = (
                slice(k * block[2], (k + 1) * block[2]),
                slice(j * block[1], (j + 1) * block[1]),
                slice(i * block[0], (i + 1) * block[0]),
            )
            for name in REQUIRED_FIELDS:
                values = np.asarray(raw["mb_data"][name][local], dtype=np.float64)
                if values.shape != (block[2], block[1], block[0]):
                    raise UnsupportedProduct(f"snapshot field shape differs: {name}")
                fields[name][slices] = values
            row += 1
    if row != len(expected) or any(not np.isfinite(value).all() for value in fields.values()):
        raise UnsupportedProduct("snapshot reconstruction is incomplete or nonfinite")
    if np.any(fields["dens"] <= 0.0):
        raise UnsupportedProduct("snapshot density is nonpositive")
    lengths = tuple(
        float(reference[f"x{axis}max"] - reference[f"x{axis}min"])
        for axis in (1, 2, 3)
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in lengths):
        raise UnsupportedProduct("snapshot domain lengths are invalid")
    return fields, lengths, bindings


def histogram(
    values: np.ndarray, bins: int, value_range: tuple[float, float]
) -> dict[str, object]:
    """Return one density-normalized finite histogram."""

    density, edges = np.histogram(values.ravel(), bins=bins, range=value_range, density=True)
    if not np.isfinite(density).all():
        raise UnsupportedProduct("histogram density is nonfinite")
    return {"edges": edges.tolist(), "density": density.tolist()}


def joint_histogram(
    x: np.ndarray,
    y: np.ndarray,
    bins: int,
    value_range: tuple[tuple[float, float], tuple[float, float]],
) -> dict[str, object]:
    """Return one density-normalized two-dimensional histogram."""

    density, x_edges, y_edges = np.histogram2d(
        x.ravel(), y.ravel(), bins=bins, range=value_range, density=True
    )
    if not np.isfinite(density).all():
        raise UnsupportedProduct("joint histogram density is nonfinite")
    return {
        "x_edges": x_edges.tolist(),
        "y_edges": y_edges.tolist(),
        "density": density.tolist(),
    }


def snapshot_scalar_fields(fields: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Return supported scalar fields for PDF products."""

    density = fields["dens"]
    p_parallel = fields["eint"]
    p_perp = fields["p_perp"]
    b2 = fields["bcc1"] ** 2 + fields["bcc2"] ** 2 + fields["bcc3"] ** 2
    return {
        "density_fluctuation": density / np.mean(density) - 1.0,
        "beta_delta": 2.0 * (p_perp - p_parallel) / np.maximum(b2, np.finfo(float).tiny),
    }


def pressure_density_fields(
    fields: dict[str, np.ndarray],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return Figure 2(a)-style pressure-density coordinates."""

    density = fields["dens"]
    p_parallel = fields["eint"]
    p_perp = fields["p_perp"]
    mean_pressure = float(np.mean((2.0 * p_perp + p_parallel) / 3.0))
    density_coordinate = mean_pressure * (density / np.mean(density) - 1.0)
    return {
        "parallel": (density_coordinate, p_parallel - np.mean(p_parallel)),
        "perpendicular": (density_coordinate, p_perp - np.mean(p_perp)),
    }


def shell_indices(shape: tuple[int, int, int], lengths: tuple[float, float, float], dk: float) -> np.ndarray:
    """Return perpendicular Fourier-shell indices without allocating 3-D k grids."""

    _, ny, nx = shape
    lx, ly, _ = lengths
    kx = 2.0 * math.pi * np.fft.fftfreq(nx, d=lx / nx)
    ky = 2.0 * math.pi * np.fft.fftfreq(ny, d=ly / ny)
    perpendicular = np.sqrt(ky[:, None] ** 2 + kx[None, :] ** 2)
    shell_2d = np.floor(perpendicular / dk + 1.0e-12).astype(np.int64)
    return np.broadcast_to(shell_2d[None, :, :], shape)


def shell_spectrum(
    fields: Iterable[np.ndarray], lengths: tuple[float, float, float], dk: float
) -> dict[str, object]:
    """Return a perpendicular shell-summed spectrum."""

    values = list(fields)
    shell = shell_indices(values[0].shape, lengths, dk)
    power = np.zeros(values[0].shape, dtype=np.float64)
    normalizer = float(values[0].size)
    for value in values:
        transformed = np.fft.fftn(value - np.mean(value)) / normalizer
        power += np.abs(transformed) ** 2
    binned = np.bincount(shell.ravel(), weights=power.ravel())
    return {
        "dk": dk,
        "k": (np.arange(len(binned), dtype=float) * dk).tolist(),
        "power_per_dk": (binned / dk).tolist(),
        "perpendicular": True,
        "normalization_definition": (
            "sum_shell |FFT(field - mean(field))/N|^2 / dk"
        ),
    }


def periodic_gradient(
    value: np.ndarray, lengths: tuple[float, float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return second-order periodic gradients in x/y/z order."""

    nz, ny, nx = value.shape
    spacings = (lengths[0] / nx, lengths[1] / ny, lengths[2] / nz)
    return tuple(
        (np.roll(value, -1, axis=axis) - np.roll(value, 1, axis=axis))
        / (2.0 * spacing)
        for axis, spacing in zip((2, 1, 0), spacings)
    )


def pressure_transfer(
    density: np.ndarray,
    velocity: list[np.ndarray],
    magnetic: list[np.ndarray],
    delta_p: np.ndarray,
    lengths: tuple[float, float, float],
    dk: float,
) -> dict[str, object]:
    """Return the MKS24 CGL pressure-stress transfer shell partition.

    This is integral <sqrt(rho) u>_k dot [(B/sqrt(rho)) dot grad
    ((Delta p/B^2) B)], filtered in k_perp shells and normalized by the
    MKS24 Kolmogorov estimate T_total ~= E_K (2 pi u_rms / L_perp).
    """

    if len(velocity) != 3 or len(magnetic) != 3:
        raise UnsupportedProduct("pressure transfer requires three-vector fields")
    arrays = [density, *velocity, *magnetic, delta_p]
    if any(value.shape != density.shape for value in arrays[1:]):
        raise UnsupportedProduct("pressure-transfer fields have inconsistent shapes")
    if any(not np.isfinite(value).all() for value in arrays):
        raise UnsupportedProduct("pressure-transfer fields are nonfinite")
    if np.any(density <= 0.0):
        raise UnsupportedProduct("pressure-transfer normalization requires positive density")
    if any(not math.isfinite(value) or value <= 0.0 for value in lengths):
        raise UnsupportedProduct("pressure-transfer domain lengths are invalid")
    if not math.isfinite(dk) or dk <= 0.0:
        raise UnsupportedProduct("pressure-transfer shell spacing is invalid")

    b2 = sum(component * component for component in magnetic)
    safe_b2 = np.maximum(b2, np.finfo(float).tiny)
    root_density = np.sqrt(density)
    weighted_velocity = [root_density * component for component in velocity]
    stress_vector = [delta_p * component / safe_b2 for component in magnetic]
    directional = []
    for component in stress_vector:
        gradient = periodic_gradient(component, lengths)
        directional.append(
            sum(magnetic[index] * gradient[index] for index in range(3))
            / root_density
        )
    volume = math.prod(lengths)
    direct = volume * float(np.mean(sum(
        weighted_velocity[index] * directional[index] for index in range(3)
    )))
    shells = shell_indices(density.shape, lengths, dk)
    shell_2d = np.asarray(shells[0])
    shell_count = int(shell_2d.max()) + 1
    shell_cross_power = np.zeros(shell_count, dtype=np.float64)
    for weighted_component, directional_component in zip(
        weighted_velocity, directional
    ):
        weighted_fourier = np.fft.fftn(weighted_component)
        directional_fourier = np.fft.fftn(directional_component)
        np.conj(directional_fourier, out=directional_fourier)
        weighted_fourier *= directional_fourier
        shell_cross_power += np.bincount(
            shell_2d.ravel(),
            weights=np.sum(weighted_fourier.real, axis=0).ravel(),
            minlength=shell_count,
        )
    # Parseval's theorem is exactly the analyzer's shell-filtered real-space
    # product, without hundreds of inverse transforms on production grids.
    transfer = (
        volume * shell_cross_power / float(density.size ** 2)
    ).tolist()
    kinetic_energy = volume * float(np.mean(
        0.5 * density * sum(component * component for component in velocity)
    ))
    velocity_rms = float(np.sqrt(np.mean(
        sum(component * component for component in velocity)
    )))
    lperp = math.sqrt(lengths[0] * lengths[1])
    total_transfer_rate = kinetic_energy * (2.0 * math.pi * velocity_rms / lperp)
    normalization_available = bool(
        math.isfinite(total_transfer_rate)
        and total_transfer_rate > np.finfo(float).tiny
    )
    shell_sum = float(sum(transfer))
    if not all(math.isfinite(value) for value in (
        *transfer,
        kinetic_energy,
        velocity_rms,
        total_transfer_rate,
        direct,
        shell_sum,
        shell_sum - direct,
    )):
        raise UnsupportedProduct("pressure-transfer reduction is nonfinite")
    return {
        "definition": (
            "integral <sqrt(rho) u>_k dot [(B/sqrt(rho)) dot grad "
            "((Delta p/B^2) B)] over each perpendicular Fourier shell"
        ),
        "delta_p_definition": "Delta p = p_perp - p_parallel",
        "discretization": (
            "second-order centered periodic real-space gradients; full three-dimensional "
            "FFT Parseval cross-spectrum binned by k_perp shell"
        ),
        "dk": dk,
        "k_perp": (np.arange(len(transfer), dtype=np.float64) * dk).tolist(),
        "transfer": transfer,
        "normalization_available": normalization_available,
        "normalization_definition": (
            "T_total ~= E_K (2 pi u_rms / L_perp), with "
            "E_K = integral[0.5 rho |u|^2] dV, "
            "u_rms = sqrt(<|u|^2>), and L_perp = sqrt(Lx Ly)"
        ),
        "kinetic_energy": kinetic_energy,
        "velocity_rms": velocity_rms,
        "perpendicular_outer_scale": lperp,
        "total_transfer_rate": total_transfer_rate,
        "transfer_normalized_by_total": (
            [value / total_transfer_rate for value in transfer]
            if normalization_available else None
        ),
        "direct_real_space": direct,
        "shell_sum": shell_sum,
        "closure_error": shell_sum - direct,
    }


def eddy_anisotropy_curve(
    bin_centers: np.ndarray, perpendicular: np.ndarray, parallel: np.ndarray
) -> dict[str, object]:
    """Invert perpendicular and parallel structure functions at common power."""

    if (
        bin_centers.ndim != 1
        or perpendicular.shape != bin_centers.shape
        or parallel.shape != bin_centers.shape
        or not np.isfinite(bin_centers).all()
        or np.any(bin_centers <= 0.0)
        or np.any(np.diff(bin_centers) <= 0.0)
    ):
        raise UnsupportedProduct("eddy-anisotropy structure-function grid is invalid")
    parallel_valid = np.isfinite(parallel) & (parallel > 0.0)
    perpendicular_valid = np.isfinite(perpendicular) & (perpendicular > 0.0)
    parallel_length = bin_centers[parallel_valid]
    parallel_power = parallel[parallel_valid]
    monotonic: list[int] = []
    maximum = -math.inf
    for index, value in enumerate(parallel_power):
        if value > maximum:
            monotonic.append(index)
            maximum = float(value)
    if len(monotonic) < 2:
        return {
            "available": False,
            "reason": "parallel structure function has fewer than two increasing bins",
        }
    parallel_length = parallel_length[monotonic]
    parallel_power = parallel_power[monotonic]
    selected = perpendicular_valid & (perpendicular >= parallel_power[0]) & (
        perpendicular <= parallel_power[-1]
    )
    if np.count_nonzero(selected) < 2:
        return {
            "available": False,
            "reason": "structure functions do not overlap on at least two bins",
        }
    ell_perp = bin_centers[selected]
    ell_parallel = np.exp(np.interp(
        np.log(perpendicular[selected]),
        np.log(parallel_power),
        np.log(parallel_length),
    ))
    return {
        "available": True,
        "ell_perp_over_lperp": ell_perp.tolist(),
        "ell_parallel_over_lperp": ell_parallel.tolist(),
    }


def local_field_eddy_anisotropy(
    velocity: list[np.ndarray],
    magnetic: list[np.ndarray],
    lengths: tuple[float, float, float],
    samples: int,
    bins: int,
    seed: int,
) -> dict[str, object]:
    """Return deterministic local-field-conditioned three-point eddy scales."""

    definition = (
        "solve S2(phi; ell_perp) = S2(phi; ell_parallel), where "
        "S2 = <|phi(x+ell) - 2 phi(x) + phi(x-ell)|^2>"
    )
    if len(velocity) != 3 or len(magnetic) != 3:
        raise UnsupportedProduct("eddy anisotropy requires three-vector fields")
    shape = velocity[0].shape
    vectors = [*velocity, *magnetic]
    if any(value.shape != shape for value in vectors[1:]):
        raise UnsupportedProduct("eddy-anisotropy fields have inconsistent shapes")
    if len(shape) != 3 or any(value <= 1 for value in shape):
        raise UnsupportedProduct("eddy anisotropy requires a three-dimensional grid")
    if any(not np.isfinite(value).all() for value in vectors):
        raise UnsupportedProduct("eddy-anisotropy fields are nonfinite")
    if samples <= 0 or bins < 2 or seed < 0:
        raise UnsupportedProduct("eddy-anisotropy sampling controls are invalid")
    if any(not math.isfinite(value) or value <= 0.0 for value in lengths):
        raise UnsupportedProduct("eddy-anisotropy domain lengths are invalid")

    nz, ny, nx = shape
    lx, ly, lz = lengths
    lperp = math.sqrt(lx * ly)
    spacing_xyz = np.asarray((lx / nx, ly / ny, lz / nz), dtype=np.float64)
    minimum = max(float(np.min(spacing_xyz)), np.finfo(float).tiny)
    maximum = 0.5 * min(lx, ly)
    if maximum <= minimum:
        return {
            "computed": False,
            "available": False,
            "definition": definition,
            "reason": "snapshot has insufficient perpendicular scale separation",
        }
    edges = np.geomspace(minimum, maximum, bins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:]) / lperp
    per_bin = max(1, int(math.ceil(samples / bins)))
    generated = per_bin * bins
    generator = np.random.Generator(np.random.PCG64(seed))
    source_bins = np.repeat(np.arange(bins), per_bin)
    radius = np.exp(generator.uniform(
        np.log(edges[source_bins]), np.log(edges[source_bins + 1])
    ))
    directions = generator.normal(size=(generated, 3))
    direction_norm = np.linalg.norm(directions, axis=1)
    if np.any(direction_norm <= np.finfo(float).tiny):
        raise UnsupportedProduct("eddy-anisotropy direction sampling is degenerate")
    directions /= direction_norm[:, None]
    offsets_xyz = np.rint(
        directions * radius[:, None] / spacing_xyz[None, :]
    ).astype(int)
    separation_xyz = offsets_xyz * spacing_xyz[None, :]
    separation = np.linalg.norm(separation_xyz, axis=1)
    retained = (separation >= edges[0]) & (separation <= edges[-1])
    retained &= np.any(offsets_xyz != 0, axis=1)
    offsets_xyz = offsets_xyz[retained]
    separation_xyz = separation_xyz[retained]
    separation = separation[retained]
    sample_bins = np.searchsorted(edges, separation, side="right") - 1
    valid_bins = (sample_bins >= 0) & (sample_bins < bins)
    offsets_xyz = offsets_xyz[valid_bins]
    separation_xyz = separation_xyz[valid_bins]
    separation = separation[valid_bins]
    sample_bins = sample_bins[valid_bins]
    count = len(sample_bins)
    center_z = generator.integers(0, nz, size=count)
    center_y = generator.integers(0, ny, size=count)
    center_x = generator.integers(0, nx, size=count)
    offset_x = offsets_xyz[:, 0]
    offset_y = offsets_xyz[:, 1]
    offset_z = offsets_xyz[:, 2]
    plus = (
        (center_z + offset_z) % nz,
        (center_y + offset_y) % ny,
        (center_x + offset_x) % nx,
    )
    center = (center_z, center_y, center_x)
    minus = (
        (center_z - offset_z) % nz,
        (center_y - offset_y) % ny,
        (center_x - offset_x) % nx,
    )
    local_field = [
        (component[plus] + component[center] + component[minus]) / 3.0
        for component in magnetic
    ]
    field_norm = np.sqrt(sum(component * component for component in local_field))
    valid_field = field_norm > np.finfo(float).tiny
    bhat = [
        component / np.maximum(field_norm, np.finfo(float).tiny)
        for component in local_field
    ]
    separation_hat = separation_xyz / separation[:, None]
    cosine = np.abs(sum(
        separation_hat[:, index] * bhat[index] for index in range(3)
    ))
    radians = math.radians(EDDY_ANGLE_DEGREES)
    parallel_selected = valid_field & (cosine >= math.cos(radians))
    perpendicular_selected = valid_field & (cosine <= math.sin(radians))

    def second_order_perp(vector: list[np.ndarray]) -> np.ndarray:
        sampled = [
            [component[location] for component in vector]
            for location in (plus, center, minus)
        ]
        perpendicular_values = []
        for values in sampled:
            field_parallel = sum(values[index] * bhat[index] for index in range(3))
            perpendicular_values.append([
                values[index] - field_parallel * bhat[index] for index in range(3)
            ])
        return sum(
            (
                perpendicular_values[0][index]
                - 2.0 * perpendicular_values[1][index]
                + perpendicular_values[2][index]
            ) ** 2
            for index in range(3)
        )

    def conditioned_mean(
        values: np.ndarray, selected: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        counts = np.bincount(sample_bins[selected], minlength=bins)
        sums = np.bincount(
            sample_bins[selected], weights=values[selected], minlength=bins
        )
        mean = np.zeros(bins, dtype=np.float64)
        populated = counts > 0
        mean[populated] = sums[populated] / counts[populated]
        return mean, counts

    output: dict[str, object] = {
        "computed": True,
        "available": True,
        "definition": definition,
        "conditioning": (
            "three-point local mean magnetic field; separation vectors within "
            f"{EDDY_ANGLE_DEGREES:g} degrees of parallel or perpendicular"
        ),
        "sampling": (
            "deterministic NumPy PCG64 random lattice separations, logarithmically "
            "balanced over normalized separation bins"
        ),
        "bit_generator": "numpy.random.PCG64",
        "separation_coordinate": (
            "|ell|/L_perp binned within each angular cone; the selected "
            "parallel or perpendicular projection differs by at most "
            "1 - cos(15 degrees)"
        ),
        "normalization_definition": "lengths divided by L_perp = sqrt(Lx Ly)",
        "lperp": lperp,
        "angle_degrees": EDDY_ANGLE_DEGREES,
        "samples_requested": samples,
        "samples_generated": generated,
        "samples_retained": int(count),
        "seed": seed,
        "bins": bins,
        "bin_centers_over_lperp": centers.tolist(),
    }
    for name, vector in (
        ("velocity_perp", velocity),
        ("magnetic_perp", magnetic),
    ):
        values = second_order_perp(vector)
        perpendicular, perpendicular_counts = conditioned_mean(
            values, perpendicular_selected
        )
        parallel, parallel_counts = conditioned_mean(values, parallel_selected)
        if not np.isfinite(perpendicular).all() or not np.isfinite(parallel).all():
            raise UnsupportedProduct("eddy-anisotropy structure functions are nonfinite")
        product = eddy_anisotropy_curve(centers, perpendicular, parallel)
        product.update({
            "perpendicular_structure_function": perpendicular.tolist(),
            "parallel_structure_function": parallel.tolist(),
            "perpendicular_sample_counts": perpendicular_counts.tolist(),
            "parallel_sample_counts": parallel_counts.tolist(),
        })
        output[name] = product
    output["available"] = bool(
        output["velocity_perp"]["available"]
        or output["magnetic_perp"]["available"]
    )
    if not output["available"]:
        output["reason"] = (
            "no eddy-anisotropy product has overlapping structure functions"
        )
    return output


def mean_pressure_transfer(records: list[dict[str, object]]) -> dict[str, object]:
    """Average compatible pressure-transfer records over retained snapshots."""

    if not records:
        raise UnsupportedProduct("pressure-transfer ensemble is empty")
    first = records[0]
    for record in records[1:]:
        for key in (
            "definition",
            "delta_p_definition",
            "discretization",
            "dk",
            "k_perp",
            "normalization_definition",
            "perpendicular_outer_scale",
        ):
            if record.get(key) != first.get(key):
                raise UnsupportedProduct("pressure-transfer snapshot coordinates differ")
    normalization_available = all(
        bool(record.get("normalization_available")) for record in records
    )
    return {
        "definition": first["definition"],
        "delta_p_definition": first["delta_p_definition"],
        "discretization": first["discretization"],
        "temporal_averaging": (
            "arithmetic mean of per-snapshot shell transfer and per-snapshot "
            "MKS24-normalized transfer"
        ),
        "dk": first["dk"],
        "k_perp": first["k_perp"],
        "transfer": weighted_array(
            [record["transfer"] for record in records],
            [1.0 / len(records)] * len(records),
        ),
        "normalization_available": normalization_available,
        "normalization_definition": first["normalization_definition"],
        "kinetic_energy_mean": float(np.mean([
            require_finite(record.get("kinetic_energy"), "pressure-transfer kinetic energy")
            for record in records
        ])),
        "velocity_rms_mean": float(np.mean([
            require_finite(record.get("velocity_rms"), "pressure-transfer velocity RMS")
            for record in records
        ])),
        "perpendicular_outer_scale": first["perpendicular_outer_scale"],
        "total_transfer_rate_mean": float(np.mean([
            require_finite(record.get("total_transfer_rate"), "pressure-transfer total rate")
            for record in records
        ])),
        "transfer_normalized_by_total": (
            weighted_array(
                [record["transfer_normalized_by_total"] for record in records],
                [1.0 / len(records)] * len(records),
            )
            if normalization_available else None
        ),
        "direct_real_space_mean": float(np.mean([
            require_finite(record.get("direct_real_space"), "pressure-transfer direct value")
            for record in records
        ])),
        "shell_sum_mean": float(np.mean([
            require_finite(record.get("shell_sum"), "pressure-transfer shell sum")
            for record in records
        ])),
        "closure_error_mean": float(np.mean([
            require_finite(record.get("closure_error"), "pressure-transfer closure error")
            for record in records
        ])),
    }


def mean_eddy_anisotropy(records: list[dict[str, object]]) -> dict[str, object]:
    """Combine sampled structure functions and invert the ensemble result."""

    computed = [record for record in records if bool(record.get("computed"))]
    if not computed:
        return {
            "computed": False,
            "available": False,
            "reason": "eddy anisotropy was not computed for selected snapshots",
        }
    metadata = (
        "definition",
        "conditioning",
        "sampling",
        "separation_coordinate",
        "normalization_definition",
        "lperp",
        "angle_degrees",
        "samples_requested",
        "samples_generated",
        "bit_generator",
        "seed",
        "bins",
        "bin_centers_over_lperp",
    )
    if any(
        any(record.get(key) != computed[0].get(key) for key in metadata)
        for record in computed[1:]
    ):
        raise UnsupportedProduct("eddy-anisotropy snapshot sampling contracts differ")
    centers = np.asarray(computed[0]["bin_centers_over_lperp"], dtype=np.float64)
    output = {name: computed[0][name] for name in metadata if name != "bin_centers_over_lperp"}
    output.update({
        "computed": True,
        "snapshot_count": len(computed),
        "ensemble_aggregation": (
            "sample-count-weighted parallel and perpendicular structure functions "
            "followed by equal-power scale inversion"
        ),
        "samples_retained": int(sum(
            require_int(record.get("samples_retained"), "eddy samples retained")
            for record in computed
        )),
        "bin_centers_over_lperp": centers.tolist(),
    })
    for name in ("velocity_perp", "magnetic_perp"):
        product: dict[str, np.ndarray] = {}
        for direction in ("perpendicular", "parallel"):
            count_name = f"{direction}_sample_counts"
            value_name = f"{direction}_structure_function"
            counts = np.sum([
                np.asarray(require_dict(record.get(name), f"eddy {name}").get(count_name),
                           dtype=np.float64)
                for record in computed
            ], axis=0)
            values = [
                np.asarray(require_dict(record.get(name), f"eddy {name}").get(value_name),
                           dtype=np.float64)
                for record in computed
            ]
            if (
                counts.shape != centers.shape
                or any(value.shape != centers.shape for value in values)
                or not np.isfinite(counts).all()
                or np.any(counts < 0.0)
                or any(not np.isfinite(value).all() for value in values)
            ):
                raise UnsupportedProduct("eddy-anisotropy snapshot samples are invalid")
            weighted = np.sum([
                value * np.asarray(
                    require_dict(record.get(name), f"eddy {name}").get(count_name),
                    dtype=np.float64,
                )
                for value, record in zip(values, computed)
            ], axis=0)
            mean = np.zeros(counts.shape, dtype=np.float64)
            populated = counts > 0.0
            mean[populated] = weighted[populated] / counts[populated]
            product[value_name] = mean
            product[count_name] = counts
        curve = eddy_anisotropy_curve(
            centers,
            product["perpendicular_structure_function"],
            product["parallel_structure_function"],
        )
        curve.update({key: value.tolist() for key, value in product.items()})
        output[name] = curve
    output["available"] = bool(
        output["velocity_perp"]["available"]
        or output["magnetic_perp"]["available"]
    )
    if not output["available"]:
        output["reason"] = "no ensemble eddy-anisotropy curve could be inverted"
    return output


def alignment_histograms(
    fields: dict[str, np.ndarray],
    lengths: tuple[float, float, float],
    shells_selected: list[int],
    bins: int,
) -> dict[str, object]:
    """Return selected-shell stretching-eigenvector alignment PDFs."""

    velocity = [fields["velx"], fields["vely"], fields["velz"]]
    magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
    b2 = sum(value * value for value in magnetic)
    bhat = [value / np.sqrt(np.maximum(b2, np.finfo(float).tiny)) for value in magnetic]
    dk = 2.0 * math.pi / lengths[2]
    shell = shell_indices(velocity[0].shape, lengths, dk)
    result: dict[str, object] = {}
    for selected in shells_selected:
        if not np.any(shell == selected):
            continue
        filtered = [
            np.fft.ifftn(np.fft.fftn(value) * (shell == selected)).real
            for value in velocity
        ]
        gradients = [periodic_gradient(value, lengths) for value in filtered]
        strain = np.empty(velocity[0].shape + (3, 3), dtype=np.float64)
        for row in range(3):
            for column in range(3):
                strain[..., row, column] = 0.5 * (
                    gradients[row][column] + gradients[column][row]
                )
        _, vectors = np.linalg.eigh(strain)
        stretching = vectors[..., :, 2]
        cosine = np.abs(
            sum(stretching[..., index] * bhat[index] for index in range(3))
        )
        result[str(selected)] = histogram(cosine, bins, (0.0, 1.0))
    return result


def padded_range(low: float, high: float) -> tuple[float, float]:
    """Return a nondegenerate deterministic histogram range."""

    if high <= low:
        delta = max(abs(low), 1.0) * 1.0e-12
        return low - delta, high + delta
    return low, high


def weighted_array(records: list[object], weights: list[float]) -> list[float]:
    """Return one deterministic weighted vector."""

    arrays = [np.asarray(value, dtype=np.float64) for value in records]
    if any(value.shape != arrays[0].shape for value in arrays[1:]):
        raise UnsupportedProduct("snapshot product arrays have inconsistent shapes")
    result = sum(weight * value for weight, value in zip(weights, arrays))
    return np.asarray(result, dtype=np.float64).tolist()


def admitted_curve_products(context: BundleContext) -> frozenset[str]:
    """Return curve products admitted for the authenticated Stage I case."""

    status = require_dict(context.stage_manifest.get("panel_status"), "panel_status")
    aliases = {
        require_text(key, "analysis case alias"): require_text(value, "analysis case")
        for key, value in require_dict(
            status.get("analysis_case_aliases"), "analysis_case_aliases"
        ).items()
    }
    products: set[str] = set()
    for product_id, value in require_dict(
        status.get("reference_product_bindings"), "reference_product_bindings"
    ).items():
        binding = require_dict(value, f"reference binding {product_id}")
        case = require_text(binding.get("case"), f"{product_id} case")
        if aliases.get(case, case) != context.case_name:
            continue
        if require_text(binding.get("kind"), f"{product_id} kind") == "curve":
            products.add(require_text(binding.get("product"), f"{product_id} product"))
    return frozenset(products)


def analyze_snapshots(
    context: BundleContext,
    start: float,
    end: float,
    bins: int,
    shells: list[int],
    max_cells: int,
    mode: str,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    """Derive supported snapshot products or return an explicit inconclusive record."""

    selected = [group for group in context.snapshots if start <= group.time <= end]
    declared = [
        {
            "time": group.time,
            "segment": group.segment,
            "representative_path": str(group.representative),
            "rank_files": list(group.rank_files),
        }
        for group in selected
    ]
    if mode == "skip":
        return {
            "result": "inconclusive",
            "reason": "snapshot processing was explicitly skipped",
            "selected_snapshot_count": len(selected),
            "declared_snapshot_inventory": declared,
        }, [], declared
    if mode != "auto":
        raise ScientificProductsError("snapshot_mode must be auto or skip")
    if not selected:
        return {
            "result": "inconclusive",
            "reason": "no authenticated snapshots lie inside the exact window",
            "selected_snapshot_count": 0,
            "declared_snapshot_inventory": [],
        }, [], []

    admitted_products = admitted_curve_products(context)
    pressure_transfer_required = any(
        product == "pressure_transfer"
        or product.startswith("pressure_transfer.")
        for product in admitted_products
    )
    eddy_anisotropy_required = any(
        product == "eddy_anisotropy"
        or product.startswith("eddy_anisotropy.")
        for product in admitted_products
    )
    bindings: list[dict[str, object]] = []
    extrema: dict[str, tuple[float, float]] = {}
    joint_extrema: dict[str, tuple[tuple[float, float], tuple[float, float]]] = {}
    first_pass: list[
        tuple[SnapshotGroup, tuple[float, float, float], tuple[int, int, int]]
    ] = []
    common_lengths: tuple[float, float, float] | None = None
    common_shape: tuple[int, int, int] | None = None
    try:
        for group in selected:
            fields, lengths, current_bindings = read_snapshot_group(group, max_cells)
            bindings.extend(current_bindings)
            shape = tuple(int(value) for value in fields["dens"].shape)
            if common_lengths is None:
                common_lengths = lengths
                common_shape = shape
            elif lengths != common_lengths or shape != common_shape:
                raise UnsupportedProduct(
                    "authenticated snapshots do not share one uniform analysis grid"
                )
            first_pass.append((group, lengths, shape))
            for name, value in snapshot_scalar_fields(fields).items():
                low = float(np.min(value))
                high = float(np.max(value))
                previous = extrema.get(name)
                extrema[name] = (
                    low if previous is None else min(previous[0], low),
                    high if previous is None else max(previous[1], high),
                )
            for name, coordinates in pressure_density_fields(fields).items():
                ranges = tuple(
                    (float(np.min(value)), float(np.max(value)))
                    for value in coordinates
                )
                previous = joint_extrema.get(name)
                joint_extrema[name] = (
                    ranges[0] if previous is None else (
                        min(previous[0][0], ranges[0][0]),
                        max(previous[0][1], ranges[0][1]),
                    ),
                    ranges[1] if previous is None else (
                        min(previous[1][0], ranges[1][0]),
                        max(previous[1][1], ranges[1][1]),
                    ),
                )
            del fields
        records: list[dict[str, object]] = []
        second_bindings: list[dict[str, object]] = []
        for group, expected_lengths, expected_shape in first_pass:
            fields, lengths, current_bindings = read_snapshot_group(group, max_cells)
            second_bindings.extend(current_bindings)
            if (
                lengths != expected_lengths
                or tuple(int(value) for value in fields["dens"].shape) != expected_shape
            ):
                raise ScientificProductsError("snapshot grid changed between analysis passes")
            scalar = snapshot_scalar_fields(fields)
            pressure = pressure_density_fields(fields)
            dk = 2.0 * math.pi / lengths[2]
            velocity = [fields["velx"], fields["vely"], fields["velz"]]
            magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
            record: dict[str, object] = {
                "time": group.time,
                "pdf": {
                    name: histogram(value, bins, padded_range(*extrema[name]))
                    for name, value in scalar.items()
                },
                "pressure_density_joint": {
                    name: joint_histogram(
                        x,
                        y,
                        bins,
                        (
                            padded_range(*joint_extrema[name][0]),
                            padded_range(*joint_extrema[name][1]),
                        ),
                    )
                    for name, (x, y) in pressure.items()
                },
                "spectra": {
                    "velocity": shell_spectrum(
                        velocity, lengths, dk
                    ),
                    "magnetic_fluctuation": shell_spectrum(
                        magnetic, lengths, dk
                    ),
                },
                "alignment": alignment_histograms(fields, lengths, shells, bins),
            }
            if pressure_transfer_required:
                record["pressure_transfer"] = pressure_transfer(
                    fields["dens"],
                    velocity,
                    magnetic,
                    fields["p_perp"] - fields["eint"],
                    lengths,
                    dk,
                )
            if eddy_anisotropy_required:
                record["eddy_anisotropy"] = local_field_eddy_anisotropy(
                    velocity,
                    magnetic,
                    lengths,
                    EDDY_SAMPLES,
                    EDDY_BINS,
                    EDDY_SEED,
                )
            records.append(record)
            del fields
        if bindings != second_bindings:
            raise ScientificProductsError("snapshot bindings changed between analysis passes")
    except (UnsupportedProduct, MemoryError, ValueError, np.linalg.LinAlgError) as error:
        return {
            "result": "inconclusive",
            "reason": str(error),
            "selected_snapshot_count": len(selected),
            "declared_snapshot_inventory": declared,
            "verified_rank_file_count_before_inconclusive": len(bindings),
        }, bindings, declared

    times = [float(record["time"]) for record in records]
    weights = [1.0 / len(times)] * len(times)
    pdf_names = records[0]["pdf"].keys()
    joint_names = records[0]["pressure_density_joint"].keys()
    spectra_names = records[0]["spectra"].keys()
    alignment_names = set(records[0]["alignment"].keys())
    for record in records[1:]:
        alignment_names.intersection_update(record["alignment"].keys())
    try:
        pressure_transfer_ensemble = (
            mean_pressure_transfer([
                require_dict(record.get("pressure_transfer"), "pressure transfer")
                for record in records
            ])
            if pressure_transfer_required else None
        )
        eddy_anisotropy_ensemble = (
            mean_eddy_anisotropy([
                require_dict(record.get("eddy_anisotropy"), "eddy anisotropy")
                for record in records
            ])
            if eddy_anisotropy_required else None
        )
    except (UnsupportedProduct, MemoryError, ValueError, np.linalg.LinAlgError) as error:
        return {
            "result": "inconclusive",
            "reason": str(error),
            "selected_snapshot_count": len(selected),
            "declared_snapshot_inventory": declared,
            "verified_rank_file_count_before_inconclusive": len(bindings),
        }, bindings, declared
    ensemble = {
        "result": "pass",
        "reason": "supported snapshot products derived from authenticated raw rank files",
        "snapshot_count": len(records),
        "snapshot_times": times,
        "temporal_averaging": "arithmetic-mean-over-retained-snapshots",
        "pdf": {
            name: {
                "edges": records[0]["pdf"][name]["edges"],
                "density": weighted_array(
                    [record["pdf"][name]["density"] for record in records], weights
                ),
            }
            for name in pdf_names
        },
        "pressure_density_joint": {
            name: {
                "x_edges": records[0]["pressure_density_joint"][name]["x_edges"],
                "y_edges": records[0]["pressure_density_joint"][name]["y_edges"],
                "density": weighted_array(
                    [record["pressure_density_joint"][name]["density"] for record in records],
                    weights,
                ),
            }
            for name in joint_names
        },
        "spectra": {
            name: {
                "dk": records[0]["spectra"][name]["dk"],
                "k": records[0]["spectra"][name]["k"],
                "power_per_dk": weighted_array(
                    [record["spectra"][name]["power_per_dk"] for record in records],
                    weights,
                ),
                "perpendicular": True,
            }
            for name in spectra_names
        },
        "alignment": {
            name: {
                "edges": records[0]["alignment"][name]["edges"],
                "density": weighted_array(
                    [record["alignment"][name]["density"] for record in records], weights
                ),
            }
            for name in sorted(alignment_names, key=int)
        },
        **(
            {"pressure_transfer": pressure_transfer_ensemble}
            if pressure_transfer_required else {}
        ),
        **(
            {"eddy_anisotropy": eddy_anisotropy_ensemble}
            if eddy_anisotropy_required else {}
        ),
        "admitted_curve_products": sorted(admitted_products),
        "declared_snapshot_inventory": declared,
    }
    return ensemble, bindings, declared


def curve_from_product(
    context: BundleContext, ensemble: dict[str, object], product: str
) -> tuple[np.ndarray, np.ndarray]:
    """Return one supported derived curve."""

    if product == "history.unstable_fraction":
        user = context.user
        required = ("time", "volume", "mirror_vol", "fire_vol")
        if any(name not in user for name in required):
            raise UnsupportedProduct("user history lacks unstable-fraction columns")
        values = [
            (mirror + fire) / volume
            for mirror, fire, volume in zip(
                user["mirror_vol"], user["fire_vol"], user["volume"]
            )
        ]
        return np.asarray(user["time"]), np.asarray(values)
    if ensemble.get("result") != "pass":
        raise UnsupportedProduct("snapshot products are inconclusive")
    family, _, name = product.partition(".")
    if family == "pdf":
        record = require_dict(
            require_dict(ensemble.get("pdf"), "ensemble pdf").get(name),
            f"pdf {name}",
        )
        edges = np.asarray(record["edges"], dtype=np.float64)
        return 0.5 * (edges[1:] + edges[:-1]), np.asarray(record["density"], dtype=np.float64)
    if family == "alignment":
        record = require_dict(
            require_dict(ensemble.get("alignment"), "ensemble alignment").get(name),
            f"alignment {name}",
        )
        edges = np.asarray(record["edges"], dtype=np.float64)
        return 0.5 * (edges[1:] + edges[:-1]), np.asarray(record["density"], dtype=np.float64)
    if product == "alignment_peak.cos_theta":
        return alignment_peak_curve(ensemble)
    if family == "pressure_transfer":
        record = require_dict(
            ensemble.get("pressure_transfer"), "ensemble pressure transfer"
        )
        if name == "transfer":
            values = record.get("transfer")
        elif name == "transfer_normalized_by_total":
            if not bool(record.get("normalization_available")):
                raise UnsupportedProduct(
                    "pressure-transfer normalization is unavailable"
                )
            values = record.get("transfer_normalized_by_total")
        else:
            raise UnsupportedProduct(f"unsupported reference curve product: {product}")
        if not isinstance(values, list):
            raise UnsupportedProduct(f"derived reference curve is missing: {product}")
        return (
            np.asarray(record.get("k_perp"), dtype=np.float64),
            np.asarray(values, dtype=np.float64),
        )
    if family == "eddy_anisotropy":
        ensemble_eddy = require_dict(
            ensemble.get("eddy_anisotropy"), "ensemble eddy anisotropy"
        )
        record = require_dict(ensemble_eddy.get(name), f"eddy anisotropy {name}")
        if not bool(record.get("available")):
            raise UnsupportedProduct(require_text(
                record.get("reason", ensemble_eddy.get("reason", "eddy curve is unavailable")),
                f"eddy anisotropy {name} reason",
            ))
        return (
            np.asarray(record.get("ell_perp_over_lperp"), dtype=np.float64),
            np.asarray(record.get("ell_parallel_over_lperp"), dtype=np.float64),
        )
    raise UnsupportedProduct(f"unsupported reference curve product: {product}")


def alignment_peak_curve(ensemble: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """Return peak stretching-eigenvector alignment against k_perp."""

    if ensemble.get("result") != "pass":
        raise UnsupportedProduct("snapshot products are inconclusive")
    alignment = require_dict(ensemble.get("alignment"), "ensemble alignment")
    spectra = require_dict(ensemble.get("spectra"), "ensemble spectra")
    velocity = require_dict(spectra.get("velocity"), "velocity spectrum")
    dk = require_finite(velocity.get("dk"), "velocity dk")
    peaks: list[tuple[float, float]] = []
    for shell, value in sorted(alignment.items(), key=lambda item: int(item[0])):
        record = require_dict(value, f"alignment shell {shell}")
        edges = np.asarray(record["edges"], dtype=np.float64)
        density = np.asarray(record["density"], dtype=np.float64)
        centers = 0.5 * (edges[1:] + edges[:-1])
        peaks.append((int(shell) * dk, float(centers[np.argmax(density)])))
    if len(peaks) < 2:
        raise UnsupportedProduct("peak alignment requires at least two shells")
    return (
        np.asarray([value[0] for value in peaks]),
        np.asarray([value[1] for value in peaks]),
    )


def surface_from_product(
    ensemble: dict[str, object], product: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return one supported derived surface."""

    if ensemble.get("result") != "pass":
        raise UnsupportedProduct("snapshot products are inconclusive")
    prefix, _, name = product.partition(".")
    if prefix != "pressure_density_joint":
        raise UnsupportedProduct(f"unsupported reference surface product: {product}")
    record = require_dict(
        require_dict(ensemble.get("pressure_density_joint"), "pressure-density products").get(name),
        f"pressure-density {name}",
    )
    x_edges = np.asarray(record["x_edges"], dtype=np.float64)
    y_edges = np.asarray(record["y_edges"], dtype=np.float64)
    return (
        0.5 * (x_edges[1:] + x_edges[:-1]),
        0.5 * (y_edges[1:] + y_edges[:-1]),
        np.asarray(record["density"], dtype=np.float64),
    )


def interpolate_curve(
    target: np.ndarray, source_x: np.ndarray, source_y: np.ndarray, method: str
) -> np.ndarray:
    """Interpolate one derived curve onto reference coordinates."""

    finite = np.isfinite(source_x) & np.isfinite(source_y)
    if method == "loglog":
        finite &= (source_x > 0.0) & (source_y > 0.0)
        if np.any(target <= 0.0):
            raise UnsupportedProduct("loglog reference coordinates are nonpositive")
    source_x = source_x[finite]
    source_y = source_y[finite]
    if source_x.size < 2 or np.any(np.diff(source_x) <= 0.0):
        raise UnsupportedProduct("derived reference curve is not finite and ordered")
    if np.any(target < source_x[0]) or np.any(target > source_x[-1]):
        raise UnsupportedProduct("reference coordinates lie outside derived product range")
    if method == "linear":
        return np.interp(target, source_x, source_y)
    if method == "loglog":
        return np.exp(np.interp(np.log(target), np.log(source_x), np.log(source_y)))
    raise UnsupportedProduct(f"unsupported interpolation method: {method}")


def interpolate_surface(
    x: np.ndarray,
    y: np.ndarray,
    source_x: np.ndarray,
    source_y: np.ndarray,
    source_z: np.ndarray,
) -> np.ndarray:
    """Bilinearly interpolate one rectangular derived surface."""

    if (
        source_z.shape != (source_x.size, source_y.size)
        or np.any(np.diff(source_x) <= 0.0)
        or np.any(np.diff(source_y) <= 0.0)
        or np.any(x < source_x[0])
        or np.any(x > source_x[-1])
        or np.any(y < source_y[0])
        or np.any(y > source_y[-1])
    ):
        raise UnsupportedProduct("reference surface lies outside a valid derived grid")
    values = []
    for x_value, y_value in zip(x, y):
        row = np.asarray([
            np.interp(x_value, source_x, source_z[:, index])
            for index in range(source_y.size)
        ])
        values.append(np.interp(y_value, source_y, row))
    return np.asarray(values)


def csv_rows(payload: bytes, label: str) -> list[dict[str, float]]:
    """Parse one finite numeric reference CSV."""

    try:
        reader = csv.DictReader(io.StringIO(payload.decode("utf-8")))
        rows = [
            {str(key): float(value) for key, value in row.items()}
            for row in reader
        ]
    except (UnicodeDecodeError, ValueError, TypeError) as error:
        raise ScientificProductsError(f"{label} is not a finite numeric CSV") from error
    if not rows or any(not all(math.isfinite(value) for value in row.values()) for row in rows):
        raise ScientificProductsError(f"{label} is empty or nonfinite")
    return rows


def reference_products_for_case(
    context: BundleContext, ensemble: dict[str, object]
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    """Compute exact residual vectors for every admitted supported case product."""

    status = require_dict(context.stage_manifest.get("panel_status"), "panel_status")
    aliases = {
        str(key): str(value)
        for key, value in require_dict(
            status.get("analysis_case_aliases"), "analysis_case_aliases"
        ).items()
    }
    manifest_bindings = require_dict(
        status.get("reference_manifests"), "reference_manifests"
    )
    product_bindings = require_dict(
        status.get("reference_product_bindings"), "reference_product_bindings"
    )
    panel_by_product: dict[str, str] = {}
    for panel_value in require_list(status.get("panels"), "panel_status panels"):
        panel = require_dict(panel_value, "panel")
        panel_id = require_text(panel.get("id"), "panel id")
        for product in require_list(panel.get("reference_products"), "panel reference_products"):
            product_id = require_text(product, "panel reference product")
            if product_id in panel_by_product:
                raise ScientificProductsError("reference product appears in multiple panels")
            panel_by_product[product_id] = panel_id

    comparisons: dict[str, object] = {}
    bindings: list[dict[str, object]] = []
    deferred: list[dict[str, object]] = []
    for product_id, binding_value in sorted(product_bindings.items()):
        binding = require_dict(binding_value, f"reference binding {product_id}")
        reference_case = require_text(binding.get("case"), f"{product_id} case")
        if aliases.get(reference_case, reference_case) != context.case_name:
            continue
        kind = require_text(binding.get("kind"), f"{product_id} kind")
        product = require_text(binding.get("product"), f"{product_id} product")
        manifest_id = require_text(binding.get("reference_manifest"), f"{product_id} manifest")
        manifest_declared = require_dict(
            manifest_bindings.get(manifest_id), f"reference manifest {manifest_id}"
        )
        manifest_path = path_within(
            context.reference_root / require_text(manifest_declared.get("path"), "reference manifest path"),
            context.reference_root,
            "reference manifest",
        )
        manifest, manifest_binding = load_json(
            manifest_path,
            f"reference manifest {manifest_id}",
            expected_sha256=require_sha256(
                manifest_declared.get("sha256"), f"reference manifest {manifest_id} sha256"
            ),
        )
        bindings.append(manifest_binding)
        collection_name = "curves" if kind == "curve" else "surfaces" if kind == "surface" else None
        if collection_name is None:
            raise ScientificProductsError(f"reference product {product_id} has invalid kind")
        entries = [
            require_dict(item, f"reference {kind}")
            for item in require_list(manifest.get(collection_name), collection_name)
            if require_dict(item, f"reference {kind}").get("id") == product_id
        ]
        if len(entries) != 1:
            raise ScientificProductsError(f"reference product {product_id} is absent or duplicated")
        entry = entries[0]
        expected = {
            "case": reference_case,
            "product": product,
            "data_file": binding.get("data_file"),
            "data_sha256": binding.get("data_sha256"),
        }
        if any(entry.get(key) != value for key, value in expected.items()):
            raise ScientificProductsError(f"reference product {product_id} differs from Stage I binding")
        data_path = path_within(
            manifest_path.parent / require_text(entry.get("data_file"), f"{product_id} data_file"),
            manifest_path.parent,
            f"reference CSV {product_id}",
        )
        payload, data_binding = read_regular_bytes(
            data_path,
            f"reference CSV {product_id}",
            expected_sha256=require_sha256(entry.get("data_sha256"), f"{product_id} data_sha256"),
        )
        bindings.append(data_binding)
        rows = csv_rows(payload, f"reference CSV {product_id}")
        try:
            if kind == "curve":
                required = {"x", "y", "y_uncertainty"}
                if not required.issubset(rows[0]):
                    raise ScientificProductsError(f"reference curve {product_id} lacks required columns")
                x = np.asarray([row["x"] for row in rows])
                reference = np.asarray([row["y"] for row in rows])
                uncertainty = np.asarray([row["y_uncertainty"] for row in rows])
                if np.any(uncertainty <= 0.0) or np.any(np.diff(x) <= 0.0):
                    raise ScientificProductsError(f"reference curve {product_id} is invalid")
                source_x, source_y = curve_from_product(context, ensemble, product)
                simulated = interpolate_curve(
                    x, source_x, source_y, str(entry.get("interpolation", "linear"))
                )
                coordinate_record: dict[str, object] = {"x": x.tolist()}
            else:
                required = {"x", "y", "z", "z_uncertainty"}
                if not required.issubset(rows[0]):
                    raise ScientificProductsError(f"reference surface {product_id} lacks required columns")
                x = np.asarray([row["x"] for row in rows])
                y = np.asarray([row["y"] for row in rows])
                reference = np.asarray([row["z"] for row in rows])
                uncertainty = np.asarray([row["z_uncertainty"] for row in rows])
                if np.any(uncertainty <= 0.0):
                    raise ScientificProductsError(f"reference surface {product_id} is invalid")
                source_x, source_y, source_z = surface_from_product(ensemble, product)
                simulated = interpolate_surface(x, y, source_x, source_y, source_z)
                coordinate_record = {"x": x.tolist(), "y": y.tolist()}
            residual = simulated - reference
            normalized = residual / uncertainty
            interpolation = str(
                entry.get("interpolation", "bilinear" if kind == "surface" else "linear")
            )
            value_record = (
                {
                    "reference_y": reference.tolist(),
                    "reference_y_uncertainty": uncertainty.tolist(),
                    "simulated_y": simulated.tolist(),
                }
                if kind == "curve"
                else {
                    "reference_z": reference.tolist(),
                    "reference_z_uncertainty": uncertainty.tolist(),
                    "simulated_z": simulated.tolist(),
                }
            )
            comparisons[str(product_id)] = {
                "available": True,
                "panel_id": panel_by_product.get(str(product_id)),
                "kind": kind,
                "case": reference_case,
                "analysis_case": context.case_name,
                "product": product,
                "stage_i_binding_validated": True,
                "reference_data_file": require_text(
                    entry.get("data_file"), f"{product_id} data_file"
                ),
                "reference_manifest_sha256": manifest_binding["sha256"],
                "data_file": data_binding["path"],
                "data_sha256": data_binding["sha256"],
                "interpolation": interpolation,
                "sample_count": int(reference.size),
                "reference_manifest": manifest_binding,
                "reference_data": data_binding,
                **coordinate_record,
                **value_record,
                "reference_values": reference.tolist(),
                "reference_uncertainty": uncertainty.tolist(),
                "simulated_values": simulated.tolist(),
                "residual": residual.tolist(),
                "normalized_residual": normalized.tolist(),
                "normalized_residual_rms": float(np.sqrt(np.mean(normalized ** 2))),
                "maximum_absolute_normalized_residual": float(np.max(np.abs(normalized))),
                "rms_residual": float(np.sqrt(np.mean(residual ** 2))),
                "maximum_absolute_residual": float(np.max(np.abs(residual))),
                "rms_normalized_by_reported_uncertainty": float(
                    np.sqrt(np.mean(normalized ** 2))
                ),
            }
        except UnsupportedProduct as error:
            record = {
                "product_id": str(product_id),
                "panel_id": panel_by_product.get(str(product_id)),
                "product": product,
                "reason": str(error),
            }
            deferred.append(record)
            comparisons[str(product_id)] = {"available": False, **record}
    return comparisons, bindings, deferred


def scientific_metrics(
    context: BundleContext, start: float, end: float
) -> tuple[dict[str, object], dict[str, object]]:
    """Derive raw scientific-acceptance metrics and LF increments."""

    metrics: dict[str, object] = {}
    for metric, (source, column) in HISTORY_METRICS.items():
        history = context.user if source == "user" else context.mhd
        if column not in history:
            raise ScientificProductsError(f"{source} history lacks required column {column}")
        metrics[metric] = {
            "source_history": source,
            "source_column": column,
            **metric_record(history["time"], history[column], start, end),
        }
    if all(name in context.user for name in ("mirror_vol", "fire_vol", "volume")):
        unstable = [
            (mirror + fire) / volume
            for mirror, fire, volume in zip(
                context.user["mirror_vol"],
                context.user["fire_vol"],
                context.user["volume"],
            )
        ]
        metrics["unstable_occupancy"] = {
            "source_history": "user",
            "source_columns": ["mirror_vol", "fire_vol", "volume"],
            "normalization": "(mirror_vol + fire_vol) / volume",
            **metric_record(context.user["time"], unstable, start, end),
        }
    if all(name in context.user for name in ("force_prp2", "force_prl2")):
        parallel_fraction = [
            parallel / max(parallel + perpendicular, 1.0e-300)
            for perpendicular, parallel in zip(
                context.user["force_prp2"], context.user["force_prl2"]
            )
        ]
        metrics["parallel_forcing_fraction"] = {
            "source_history": "user",
            "source_columns": ["force_prp2", "force_prl2"],
            **metric_record(context.user["time"], parallel_fraction, start, end),
        }
    increments: dict[str, object] = {}
    for column in LF_COUNTERS:
        if column not in context.mhd:
            increments[column] = {"available": False, "reason": "column is absent"}
            continue
        increments[column] = {
            "available": True,
            "start_value": interpolate_at(context.mhd["time"], context.mhd[column], start),
            "end_value": interpolate_at(context.mhd["time"], context.mhd[column], end),
        }
        increments[column]["increment"] = (
            increments[column]["end_value"] - increments[column]["start_value"]
        )
    return metrics, increments


def convergence_products(
    ensemble: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Return R16/R02/R17-ready products in reviewed physical-k units."""

    if ensemble.get("result") != "pass":
        return {}, [{
            "product": "R16/R02/R17 convergence",
            "reason": "snapshot products are inconclusive",
        }]
    deferred: list[dict[str, object]] = []
    result: dict[str, object] = {}
    try:
        x, y = alignment_peak_curve(ensemble)
        result["peak_alignment"] = {
            "x_coordinate": "physical_k_perp",
            "coordinate_definition": "physical_k_perp = (k_perp_over_pi) * exact pi",
            "x": x.tolist(),
            "x_over_pi": (x / math.pi).tolist(),
            "y": y.tolist(),
        }
    except UnsupportedProduct as error:
        deferred.append({"product": "peak_alignment", "reason": str(error)})
    spectra = require_dict(ensemble.get("spectra"), "ensemble spectra")
    for source, target in (
        ("velocity", "velocity_spectrum_shape"),
        ("magnetic_fluctuation", "magnetic_fluctuation_spectrum_shape"),
    ):
        record = require_dict(spectra.get(source), f"{source} spectrum")
        k = np.asarray(record["k"], dtype=np.float64)
        power = np.asarray(record["power_per_dk"], dtype=np.float64)
        coordinate = k / math.pi
        selected = (coordinate >= 4.0) & (coordinate <= 24.0)
        if np.count_nonzero(selected) < 3:
            deferred.append({
                "product": target,
                "reason": "spectrum lacks at least three samples on 4 <= k_perp/pi <= 24",
            })
            continue
        x = k[selected]
        y = power[selected]
        integral = float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))
        if not math.isfinite(integral) or integral <= 0.0:
            deferred.append({"product": target, "reason": "spectrum normalization is nonpositive"})
            continue
        result[target] = {
            "x_coordinate": "physical_k_perp",
            "coordinate_definition": "physical_k_perp = (k_perp_over_pi) * exact pi",
            "normalization": "unit integral over 4*pi <= physical_k_perp <= 24*pi",
            "x": x.tolist(),
            "x_over_pi": coordinate[selected].tolist(),
            "y": (y / integral).tolist(),
        }
    return result, deferred


def deduplicate_bindings(values: list[dict[str, object]]) -> list[dict[str, object]]:
    """Return stable path-sorted unique bindings."""

    by_path: dict[str, dict[str, object]] = {}
    for value in values:
        path = require_text(value.get("path"), "binding path")
        previous = by_path.get(path)
        if previous is not None and previous != value:
            raise ScientificProductsError(f"input binding changed within generation: {path}")
        by_path[path] = value
    return [by_path[path] for path in sorted(by_path)]


def scientific_acceptance_contract(
    context: BundleContext,
    generator_binding: dict[str, object],
    start: float,
    end: float,
) -> dict[str, object]:
    """Return the exact nested replay contract consumed by scientific acceptance."""

    return {
        "schema_version": ACCEPTANCE_CONTRACT_SCHEMA_VERSION,
        "record_type": ACCEPTANCE_CONTRACT_RECORD_TYPE,
        "case_id": context.case_id,
        "case_name": context.case_name,
        "analysis_window": {"time_start": start, "time_end": end},
        "accepted_bundle_manifest": context.bundle_binding,
        "generator": generator_binding,
        "deterministic_replay_verification": DETERMINISTIC_REPLAY_VERIFICATION,
        "stage_i_manifest_sha256": context.stage_manifest_binding["sha256"],
        "mhd_history": context.mhd_binding,
        "user_history": context.user_binding,
    }


def scientific_acceptance_case_record(
    start: float,
    end: float,
    increments: dict[str, object],
    ensemble: dict[str, object],
    convergence: dict[str, object],
) -> dict[str, object]:
    """Return one analyzer-compatible case record from deterministic products."""

    return {
        "analysis_window": {"time_start": start, "time_end": end},
        "lf_counter_increments": increments,
        "snapshot_ensemble": ensemble,
        "scientific_acceptance_convergence": convergence,
    }


def build_evidence(request: dict[str, object]) -> dict[str, object]:
    """Recompute a complete deterministic scientific-products evidence record."""

    start = require_finite(request.get("time_start"), "time_start")
    end = require_finite(request.get("time_end"), "time_end")
    if start >= end:
        raise ScientificProductsError("time_start must be less than time_end")
    bins = require_int(request.get("pdf_bins"), "pdf_bins", 8)
    max_cells = require_int(request.get("max_snapshot_cells"), "max_snapshot_cells", 1)
    shells = [
        require_int(item, "alignment shell", 1)
        for item in require_list(request.get("alignment_shells"), "alignment_shells")
    ]
    if not shells or len(shells) != len(set(shells)):
        raise ScientificProductsError("alignment_shells must be unique and nonempty")
    context = authenticate_bundle(request)
    for history, label in ((context.mhd, "MHD"), (context.user, "user")):
        interpolate_at(history["time"], history["time"], start)
        interpolate_at(history["time"], history["time"], end)
    metrics, increments = scientific_metrics(context, start, end)
    ensemble, snapshot_bindings, declared_snapshots = analyze_snapshots(
        context,
        start,
        end,
        bins,
        shells,
        max_cells,
        require_text(request.get("snapshot_mode"), "snapshot_mode"),
    )
    references, reference_bindings, reference_deferred = reference_products_for_case(
        context, ensemble
    )
    convergence, convergence_deferred = convergence_products(ensemble)
    deferred = [
        *reference_deferred,
        *convergence_deferred,
    ]
    inputs = deduplicate_bindings([
        context.bundle_binding,
        context.stage_manifest_binding,
        *context.segment_bindings,
        *context.segment_history_bindings,
        context.mhd_binding,
        context.user_binding,
        *snapshot_bindings,
        *reference_bindings,
    ])
    generator_binding = regular_file_binding(Path(__file__), "scientific products generator")
    parser_binding = regular_file_binding(BIN_CONVERT_PATH, "Athena binary parser")
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "result": "complete" if not deferred and ensemble.get("result") == "pass" else "inconclusive",
        "request": request,
        "generator": generator_binding,
        "runtime": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "byteorder": sys.byteorder,
            "binary_parser": parser_binding,
        },
        "authentication": {
            "authority_scope": context.authority_scope,
            "case_id": context.case_id,
            "case_name": context.case_name,
            "accepted_final_time": 10.0,
            "bundle_manifest": context.bundle_binding,
            "stage_i_manifest": context.stage_manifest_binding,
            "accepted_segment_count": len(context.segment_bindings),
            "authenticated_segment_history_file_count": len(
                context.segment_history_bindings
            ),
            "selected_snapshot_count": len(declared_snapshots),
            "verified_snapshot_rank_file_count": len(snapshot_bindings),
        },
        "analysis_window": {"time_start": start, "time_end": end},
        "scientific_acceptance_contract": scientific_acceptance_contract(
            context, generator_binding, start, end
        ),
        "cases": {
            context.case_name: scientific_acceptance_case_record(
                start, end, increments, ensemble, convergence
            )
        },
        "scientific_acceptance_metrics": metrics,
        "lf_counter_increments": increments,
        "snapshot_ensemble": ensemble,
        "reference_curve_comparisons": {
            "available": bool(references),
            "comparisons": references,
        },
        "scientific_acceptance_convergence": convergence,
        "deferred_products": deferred,
        "implementation_scope": {
            "implemented_product_families": [
                {
                    "product": "pressure_transfer",
                    "definition": (
                        "MKS24 perpendicular-shell partition of the CGL "
                        "pressure-anisotropy stress transfer"
                    ),
                    "selection": "computed for admitted Stage I case roles",
                },
                {
                    "product": "eddy_anisotropy",
                    "definition": (
                        "three-point local-field-conditioned perpendicular-vector "
                        "structure functions with equal-power scale inversion"
                    ),
                    "selection": "computed for admitted Stage I case roles",
                    "sampling": {
                        "angle_degrees": EDDY_ANGLE_DEGREES,
                        "samples_per_snapshot": EDDY_SAMPLES,
                        "bins": EDDY_BINS,
                        "seed": EDDY_SEED,
                    },
                },
            ],
            "deferred_product_families": [],
            "status_effect": (
                "An admitted product is a case blocker only when its deterministic "
                "derivation or exact reference comparison appears in deferred_products."
            ),
        },
        "input_bindings": inputs,
        "non_authorizing_statement": (
            "This deterministic product evidence is non-authorizing until bound "
            "to independently reviewed scientific-acceptance policy and canonical "
            "campaign authority."
        ),
    }
    return seal_evidence(evidence)


def normalized_request(args: argparse.Namespace) -> dict[str, object]:
    """Return the exact deterministic request represented by CLI arguments."""

    shells = [int(value) for value in args.alignment_shells.split(",") if value.strip()]
    return {
        "authority_mode": args.authority_mode,
        "bundle_manifest": str(args.bundle_manifest.expanduser().absolute().resolve(strict=True)),
        "expected_bundle_sha256": require_sha256(
            args.expected_bundle_sha256, "expected_bundle_sha256"
        ),
        "time_start": float(args.time_start),
        "time_end": float(args.time_end),
        "reference_root": (
            str(args.reference_root.expanduser().absolute().resolve(strict=True))
            if args.reference_root is not None else None
        ),
        "snapshot_mode": args.snapshot_mode,
        "pdf_bins": int(args.pdf_bins),
        "alignment_shells": shells,
        "max_snapshot_cells": int(args.max_snapshot_cells),
    }


def load_evidence(path: Path, expected_sha256: str) -> tuple[dict[str, object], bytes]:
    """Load one exact evidence file and verify its external and self digests."""

    payload, _ = read_regular_bytes(
        path, "scientific products evidence", expected_sha256=expected_sha256
    )
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ScientificProductsError("scientific products evidence is invalid JSON") from error
    evidence = require_dict(value, "scientific products evidence")
    if evidence.get("schema_version") != SCHEMA_VERSION or evidence.get("record_type") != RECORD_TYPE:
        raise ScientificProductsError("scientific products evidence identity differs")
    verify_self_digest(evidence)
    return evidence, payload


def replay_evidence(path: Path, expected_sha256: str) -> tuple[dict[str, object], bool]:
    """Recompute one evidence record and compare exact serialized bytes."""

    evidence, payload = load_evidence(path, expected_sha256)
    request = require_dict(evidence.get("request"), "evidence request")
    recomputed = build_evidence(request)
    return recomputed, stable_json(recomputed) == payload


def command_parser() -> argparse.ArgumentParser:
    """Return the CLI parser."""

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--bundle-manifest", type=Path, required=True)
    generate.add_argument("--expected-bundle-sha256", required=True)
    generate.add_argument("--authority-mode", choices=("canonical", "offline"), required=True)
    generate.add_argument("--time-start", type=float, required=True)
    generate.add_argument("--time-end", type=float, required=True)
    generate.add_argument("--reference-root", type=Path)
    generate.add_argument("--snapshot-mode", choices=("auto", "skip"), default="auto")
    generate.add_argument("--pdf-bins", type=int, default=64)
    generate.add_argument(
        "--alignment-shells",
        default="1,2,3,4,6,8,12,16,24,32,64,128",
    )
    generate.add_argument("--max-snapshot-cells", type=int, default=150_000_000)
    generate.add_argument("--output", type=Path, required=True)

    for name in ("verify", "replay"):
        command = subparsers.add_parser(name)
        command.add_argument("--evidence", type=Path, required=True)
        command.add_argument("--expected-evidence-sha256", required=True)
        if name == "replay":
            command.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line interface."""

    args = command_parser().parse_args(argv)
    try:
        if args.command == "generate":
            evidence = build_evidence(normalized_request(args))
            write_candidate(args.output, evidence)
            result = {
                "result": evidence["result"],
                "output": str(args.output.expanduser().absolute()),
                "output_sha256": sha256_bytes(stable_json(evidence)),
                "evidence_digest": evidence["evidence_digest"]["sha256"],
                "deferred_product_count": len(evidence["deferred_products"]),
            }
        elif args.command == "verify":
            recomputed, identical = replay_evidence(
                args.evidence,
                require_sha256(args.expected_evidence_sha256, "expected_evidence_sha256"),
            )
            if not identical:
                raise ScientificProductsError(
                    "recomputed evidence is not byte-identical to retained evidence"
                )
            result = {
                "verified": True,
                "result": recomputed["result"],
                "evidence_sha256": args.expected_evidence_sha256,
            }
        else:
            recomputed, identical = replay_evidence(
                args.evidence,
                require_sha256(args.expected_evidence_sha256, "expected_evidence_sha256"),
            )
            if not identical:
                raise ScientificProductsError(
                    "retained evidence differs from deterministic replay"
                )
            write_candidate(args.output, recomputed)
            result = {
                "replayed": True,
                "byte_identical_to_source": True,
                "output": str(args.output.expanduser().absolute()),
                "output_sha256": sha256_bytes(stable_json(recomputed)),
            }
    except (ScientificProductsError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    sys.stdout.write(json.dumps(result, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
