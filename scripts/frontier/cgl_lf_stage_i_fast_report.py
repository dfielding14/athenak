#!/usr/bin/env python3
"""Analyze the direct-fast CGL-LF Stage I campaign without controller gates.

The adapter treats production run directories as read-only inputs.  It
deterministically reconstructs one restart-linked lineage per case, merges
histories, indexes exact rank-local snapshots, reuses the paper analyzer's pure
scientific kernels, and writes case, campaign, table, figure, and manuscript
artifacts beneath a separate output directory.

Unexpected scientific behavior is reported as a warning.  It never prevents
analysis of the remaining available data.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import gc
import hashlib
import importlib.util
import json
import math
import multiprocessing
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import types
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
DEFAULT_FROZEN_SOURCE = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422"
)
MATRIX_RELATIVE = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
FAST_RUNS_RELATIVE = Path("runs/mks24-stage-i-fast/E03-forcing-policy")
RACE_RUNS_RELATIVE = Path("runs/mks24-stage-i-fast-races/E03-forcing-policy")
RELAXED_RUNS_RELATIVE = Path("runs/mks24-stage-i-fast-relaxed/E03-forcing-policy")
FAST_RUNS_RELATIVES = (
    FAST_RUNS_RELATIVE,
    RACE_RUNS_RELATIVE,
    RELAXED_RUNS_RELATIVE,
)
LEGACY_RUNS_RELATIVE = Path("runs/mks24-stage-i/E03-forcing-policy")
R02_BUNDLE_RELATIVE = LEGACY_RUNS_RELATIVE / "bundles/R02"
DEFAULT_OUTPUT_RELATIVE = Path(
    "analysis/mks24-stage-i-fast/E03-forcing-policy/R02-R17-final"
)
REFERENCE_ROOT_RELATIVE = Path(
    "build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2"
)
WRITING_GUIDE = DEFAULT_ROOT / "writing_guide.md"
TARGET_TIME = 10.0
TIME_TOLERANCE = 1.0e-12
RESTART_TIME_TOLERANCE = 5.0e-6
DEFAULT_SNAPSHOT_WORKERS = 4
DEFAULT_SNAPSHOT_MEMORY_BUDGET_GIB = 384.0
SNAPSHOT_PEAK_BYTES_PER_INPUT_BYTE = 24.0
GIBIBYTE = 1024 ** 3
ACTIVE_SLURM_STATES = {
    "CONFIGURING",
    "COMPLETING",
    "PENDING",
    "RUNNING",
    "SUSPENDED",
}
FAILED_SLURM_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "TIMEOUT",
}
STRICT_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)
FATAL_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
)
NONFATAL_HARD_BOUND_VARIANT = "finite_limiter_hard_bound_diagnostic_nonfatal"
NONFATAL_HARD_BOUND_CASES = ("R14", "R15")
MONOTONIC_LF_COUNTERS = (
    "lf_nstage",
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_mirror",
    "lf_firehs",
    "lf_hardbd",
    "lf_hwproj",
    "lf_qface",
    "lf_qprcap",
    "lf_qpr10",
    "lf_qpecap",
    "lf_qpe10",
)
WINDOWS = {
    "full": (0.0, 10.0),
    "steady": (4.0, 10.0),
    "early": (4.0, 7.0),
    "late": (7.0, 10.0),
}
ACTIVE_PASSIVE_PAIRS = (
    ("R02", "R06"),
    ("R03", "R07"),
    ("R04", "R08"),
    ("R05", "R09"),
)
COMPARISON_GROUPS = {
    "active_passive": [list(pair) for pair in ACTIVE_PASSIVE_PAIRS],
    "forcing_geometry": [
        ["R02", "R04"],
        ["R03", "R05"],
        ["R06", "R08"],
        ["R07", "R09"],
    ],
    "beta": [
        ["R02", "R03"],
        ["R04", "R05"],
        ["R06", "R07"],
        ["R08", "R09"],
        ["R10", "R04"],
    ],
    "correlation_time": [["R05", "R11"]],
    "heat_flux": [["R12", "R02", "R06", "R13"]],
    "limiter": [["R14", "R15", "R03", "R07"]],
    "resolution": [["R16", "R02", "R17"]],
}
STEADY_SCALARS = (
    "kinetic_mean",
    "magnetic_mean",
    "therm_cgl_mean",
    "abs_dp_mean",
    "delta_p_mean",
    "beta_mean",
    "nu_eff_mean",
    "force_pwr_mean",
    "mirror_vol_fraction_mean",
    "fire_vol_fraction_mean",
    "hard_vol_fraction_mean",
    "c_b2_mean",
)
HISTORY_LABEL = re.compile(r"\[(\d+)\]=(\S+)")
FAST_SEGMENT = re.compile(r"fast_s(\d{3})_")


class ReportError(RuntimeError):
    """An unrecoverable adapter configuration or output error."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(value: object) -> str:
    payload = json.dumps(
        json_safe(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return sha256_bytes(payload)


def json_safe(value: object) -> object:
    """Return a deterministic JSON representation with nonfinite values nulled."""

    if isinstance(value, dict):
        return {
            str(key): json_safe(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as stream:
            staged = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(staged, path)
        staged = None
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def write_json(path: Path, value: object) -> None:
    atomic_write_bytes(
        path,
        (
            json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False)
            + "\n"
        ).encode("utf-8"),
    )


def write_text(path: Path, value: str) -> None:
    atomic_write_bytes(path, value.encode("utf-8"))


def load_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ReportError(f"expected JSON object: {path}")
    return value


def artifact_binding(path: Path, *, digest: bool = True) -> dict[str, object]:
    stat = path.stat()
    record: dict[str, object] = {
        "path": str(path.absolute()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if digest:
        record["sha256"] = sha256_file(path)
    return record


def read_stable_bytes(path: Path, attempts: int = 3) -> tuple[bytes, dict[str, object]]:
    """Read one possibly live file, retaining whether the observed bytes were stable."""

    last_payload = b""
    last_before = None
    last_after = None
    for _ in range(attempts):
        before = path.stat()
        payload = path.read_bytes()
        after = path.stat()
        last_payload, last_before, last_after = payload, before, after
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        ) == (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            return payload, {
                "path": str(path.absolute()),
                "size_bytes": len(payload),
                "sha256": sha256_bytes(payload),
                "stable_read": True,
                "mtime_ns": after.st_mtime_ns,
            }
    assert last_before is not None and last_after is not None
    return last_payload, {
        "path": str(path.absolute()),
        "size_bytes": len(last_payload),
        "sha256": sha256_bytes(last_payload),
        "stable_read": False,
        "mtime_ns_before": last_before.st_mtime_ns,
        "mtime_ns_after": last_after.st_mtime_ns,
    }


def expand_cases(values: Iterable[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        for item in value.split(","):
            item = item.strip()
            if not item:
                continue
            match = re.fullmatch(r"R(\d{2})-R(\d{2})", item)
            if match:
                start, end = (int(number) for number in match.groups())
                if start > end:
                    raise ReportError(f"descending case range is invalid: {item}")
                result.extend(f"R{number:02d}" for number in range(start, end + 1))
            else:
                result.append(item)
    invalid = [case for case in result if not re.fullmatch(r"R(?:0[2-9]|1[0-7])", case)]
    if invalid:
        raise ReportError(f"unsupported cases: {', '.join(invalid)}")
    return list(dict.fromkeys(result))


def safe_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(parent.resolve(strict=False))
        return True
    except ValueError:
        return False


def require_output_is_separate(root: Path, frozen_source: Path, output: Path) -> None:
    source_roots = (
        *(root / relative for relative in FAST_RUNS_RELATIVES),
        root / LEGACY_RUNS_RELATIVE,
        frozen_source,
        REPO_ROOT,
    )
    for source in source_roots:
        if safe_relative_to(output, source):
            raise ReportError(f"output directory must not be inside source tree: {source}")


def matrix_cases(matrix: Path) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    data = load_json(matrix)
    values = data.get("cases")
    if not isinstance(values, list):
        raise ReportError(f"matrix has no cases list: {matrix}")
    cases = {
        str(value["id"]): value
        for value in values
        if isinstance(value, dict) and isinstance(value.get("id"), str)
    }
    return cases, data


def parse_history_payload(
    path: Path, payload: bytes, stable_read: bool
) -> tuple[dict[str, object], list[str]]:
    warnings: list[str] = []
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ReportError(f"history is not UTF-8 text: {path}") from error
    labels: list[str] | None = None
    header: list[str] = []
    rows: list[list[float]] = []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith("#"):
            header.append(line)
            found = HISTORY_LABEL.findall(line)
            if found:
                labels = [
                    name for _, name in sorted(found, key=lambda item: int(item[0]))
                ]
            continue
        if not line.strip():
            continue
        try:
            row = [float(value) for value in line.split()]
        except ValueError:
            if index == len(lines) - 1 and not payload.endswith(b"\n"):
                warnings.append(f"ignored incomplete trailing history row: {path}")
                continue
            raise ReportError(f"history contains a nonnumeric row: {path}")
        if labels is not None and len(row) != len(labels):
            if index == len(lines) - 1 and not payload.endswith(b"\n"):
                warnings.append(f"ignored incomplete trailing history row: {path}")
                continue
            raise ReportError(f"history row width differs from header: {path}")
        if not all(math.isfinite(value) for value in row):
            warnings.append(f"history contains nonfinite row values: {path}")
        rows.append(row)
    if labels is None or not rows:
        raise ReportError(f"history lacks a labeled header or rows: {path}")
    if any(len(row) != len(labels) for row in rows):
        raise ReportError(f"history row width differs from header: {path}")
    if not stable_read:
        warnings.append(f"history changed while being read: {path}")
    return {
        "path": str(path.absolute()),
        "labels": labels,
        "header": header,
        "rows": rows,
    }, warnings


def read_history_source(path: Path) -> tuple[dict[str, object], list[str]]:
    payload, binding = read_stable_bytes(path)
    record, warnings = parse_history_payload(path, payload, bool(binding["stable_read"]))
    record["binding"] = binding
    return record, warnings


def history_paths(output: Path) -> tuple[Path | None, Path | None]:
    mhd = sorted(output.glob("*.mhd.hst"))
    user = sorted(output.glob("*.user.hst"))
    return (
        mhd[0] if len(mhd) == 1 else None,
        user[0] if len(user) == 1 else None,
    )


def peek_history_final(output: Path) -> float:
    mhd, _ = history_paths(output)
    if mhd is None:
        return -math.inf
    try:
        record, _ = read_history_source(mhd)
        return float(record["rows"][-1][0])  # type: ignore[index]
    except (OSError, ReportError, ValueError, IndexError):
        return -math.inf


def rows_close(first: list[float], second: list[float]) -> bool:
    return len(first) == len(second) and all(
        math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-12)
        for left, right in zip(first, second)
    )


def merge_histories(
    sources: list[tuple[str, Path]], destination: Path, label: str
) -> dict[str, object]:
    warnings: list[str] = []
    errors: list[str] = []
    parsed: list[tuple[str, dict[str, object]]] = []
    for segment, path in sources:
        try:
            record, source_warnings = read_history_source(path)
            parsed.append((segment, record))
            warnings.extend(source_warnings)
        except (OSError, ReportError) as error:
            errors.append(str(error))
    if not parsed:
        destination.unlink(missing_ok=True)
        return {
            "available": False,
            "label": label,
            "warnings": warnings,
            "errors": errors or [f"no {label} history sources were available"],
            "sources": [],
        }
    labels = list(parsed[0][1]["labels"])  # type: ignore[arg-type]
    header = list(parsed[0][1]["header"])  # type: ignore[arg-type]
    retained: list[list[float]] = []
    retained_segment: list[str] = []
    duplicate_rows = 0
    conflicting_duplicates = 0
    skipped_overlap_rows = 0
    for segment, record in parsed:
        if record["labels"] != labels:
            errors.append(f"{label} history columns differ in segment {segment}")
            continue
        for row in record["rows"]:  # type: ignore[union-attr]
            row = list(row)
            if not retained or row[0] > retained[-1][0] + TIME_TOLERANCE:
                retained.append(row)
                retained_segment.append(segment)
            elif math.isclose(
                row[0], retained[-1][0], rel_tol=0.0, abs_tol=TIME_TOLERANCE
            ):
                duplicate_rows += 1
                if not rows_close(row, retained[-1]):
                    conflicting_duplicates += 1
                    warnings.append(
                        f"{label} duplicate time {row[0]:.16g} differed; "
                        f"retained later lineage row from {segment}"
                    )
                    retained[-1] = row
                    retained_segment[-1] = segment
            else:
                skipped_overlap_rows += 1
                warnings.append(
                    f"{label} skipped out-of-order overlap row t={row[0]:.16g} "
                    f"from {segment}"
                )
    if not retained:
        destination.unlink(missing_ok=True)
        return {
            "available": False,
            "label": label,
            "warnings": warnings,
            "errors": errors or [f"no {label} rows survived merge"],
            "sources": [record["binding"] for _, record in parsed],
        }
    output_lines = [
        *header,
        *[" ".join(f"{value:24.16e}" for value in row) for row in retained],
    ]
    write_text(destination, "\n".join(output_lines) + "\n")
    intervals = [
        right[0] - left[0] for left, right in zip(retained, retained[1:])
        if right[0] > left[0]
    ]
    median_interval = (
        sorted(intervals)[len(intervals) // 2] if intervals else None
    )
    maximum_interval = max(intervals) if intervals else None
    return {
        "available": True,
        "label": label,
        "path": str(destination.absolute()),
        "binding": artifact_binding(destination),
        "columns": labels,
        "rows": len(retained),
        "time_first": retained[0][0],
        "time_final": retained[-1][0],
        "duplicate_rows_removed": duplicate_rows,
        "conflicting_duplicate_rows": conflicting_duplicates,
        "overlap_rows_skipped": skipped_overlap_rows,
        "median_positive_interval": median_interval,
        "maximum_positive_interval": maximum_interval,
        "source_segments": retained_segment,
        "sources": [record["binding"] for _, record in parsed],
        "warnings": warnings,
        "errors": errors,
    }


def binary_snapshot_time(path: Path) -> float:
    with path.open("rb") as stream:
        code_header = stream.readline().split()
        if not code_header or code_header[0] != b"Athena":
            raise ReportError(f"invalid Athena binary snapshot: {path}")
        count_line = stream.readline()
        try:
            count = int(count_line.split(b"=")[-1])
        except ValueError as error:
            raise ReportError(f"invalid Athena binary preheader: {path}") from error
        values: dict[str, str] = {}
        for _ in range(count - 1):
            line = stream.readline().decode("utf-8")
            if "=" not in line:
                continue
            key, value = (token.strip() for token in line.split("=", 1))
            values[key] = value
    if "time" not in values:
        raise ReportError(f"snapshot has no physical time: {path}")
    value = float(values["time"])
    if not math.isfinite(value):
        raise ReportError(f"snapshot time is nonfinite: {path}")
    return value


def snapshot_group(
    representative: Path, expected_ranks: int | None, segment: str, order: int
) -> tuple[dict[str, object], list[str]]:
    warnings: list[str] = []
    rank_root = representative.parent.parent
    rank_dirs = sorted(
        path for path in rank_root.glob("rank_*")
        if path.is_dir() and re.fullmatch(r"rank_\d{8}", path.name)
    )
    if expected_ranks is None:
        expected_ranks = len(rank_dirs)
    expected_names = [f"rank_{rank:08d}" for rank in range(expected_ranks)]
    actual_names = [path.name for path in rank_dirs]
    if actual_names != expected_names:
        warnings.append(
            f"inexact rank directory inventory for snapshot {representative}"
        )
    members = [rank_root / name / representative.name for name in expected_names]
    missing = [str(path) for path in members if not path.is_file()]
    empty = [str(path) for path in members if path.is_file() and path.stat().st_size == 0]
    try:
        time = binary_snapshot_time(representative)
    except (OSError, ReportError, ValueError) as error:
        time = None
        warnings.append(str(error))
    return {
        "segment": segment,
        "lineage_order": order,
        "time": time,
        "representative": str(representative.absolute()),
        "expected_ranks": expected_ranks,
        "rank_files": [
            {
                "path": str(path.absolute()),
                "size_bytes": path.stat().st_size if path.is_file() else None,
            }
            for path in members
        ],
        "complete": not missing and not empty and time is not None,
        "missing_rank_files": missing,
        "empty_rank_files": empty,
    }, warnings


def index_snapshots(
    sources: list[tuple[str, Path, int]], destination: Path
) -> dict[str, object]:
    warnings: list[str] = []
    groups: list[dict[str, object]] = []
    for order, (segment, output, expected_ranks) in enumerate(sources):
        rank_zero = output / "bin/rank_00000000"
        if not rank_zero.is_dir():
            continue
        for representative in sorted(rank_zero.glob("*.bin")):
            group, group_warnings = snapshot_group(
                representative, expected_ranks, segment, order
            )
            groups.append(group)
            warnings.extend(group_warnings)
    groups.sort(
        key=lambda item: (
            float(item["time"]) if item["time"] is not None else math.inf,
            int(item["lineage_order"]),
            str(item["representative"]),
        )
    )
    distinct: list[dict[str, object]] = []
    duplicates: list[dict[str, object]] = []
    for group in groups:
        if group["time"] is None:
            duplicates.append({**group, "deduplication_reason": "unreadable_time"})
            continue
        if distinct and math.isclose(
            float(group["time"]), float(distinct[-1]["time"]),
            rel_tol=0.0, abs_tol=TIME_TOLERANCE,
        ):
            replaced = distinct[-1]
            distinct[-1] = group
            duplicates.append({
                **replaced,
                "deduplication_reason": "later_lineage_snapshot_at_same_time_retained",
            })
        else:
            distinct.append(group)
    result = {
        "schema_version": 1,
        "generated_utc": utc_now(),
        "snapshot_count": len(distinct),
        "complete_snapshot_count": sum(bool(item["complete"]) for item in distinct),
        "time_first": distinct[0]["time"] if distinct else None,
        "time_last": distinct[-1]["time"] if distinct else None,
        "snapshots": distinct,
        "duplicates": duplicates,
        "warnings": warnings,
    }
    write_json(destination, result)
    return result


def fast_manifest_path(segment: Path) -> Path:
    return segment / "manifest/fast_run.json"


def normalized_slurm_state(value: str) -> str:
    """Normalize one Slurm state while preserving its operational meaning."""

    stripped = value.strip()
    return stripped.split("+", 1)[0].split()[0].upper() if stripped else "UNKNOWN"


def slurm_job_evidence(job_id: str) -> dict[str, object]:
    """Capture current scheduler evidence without making it scientific provenance."""

    try:
        queued = subprocess.run(
            ["/usr/bin/squeue", "-h", "-j", job_id, "-o", "%T"],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        queued = subprocess.CompletedProcess([], 1, stdout="", stderr="")
    queued_state = queued.stdout.strip()
    if queued_state:
        return {
            "job_id": job_id,
            "source": "squeue",
            "state": normalized_slurm_state(queued_state.splitlines()[0]),
            "exit_code": None,
        }
    try:
        accounted = subprocess.run(
            [
                "/usr/bin/sacct",
                "-X",
                "-n",
                "-P",
                "-j",
                job_id,
                "-o",
                "State,ExitCode",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        accounted = subprocess.CompletedProcess([], 1, stdout="", stderr="")
    row = accounted.stdout.strip().splitlines()
    if row:
        fields = row[0].split("|")
        return {
            "job_id": job_id,
            "source": "sacct",
            "state": normalized_slurm_state(fields[0]),
            "exit_code": fields[1] if len(fields) > 1 and fields[1] else None,
        }
    return {
        "job_id": job_id,
        "source": "unavailable",
        "state": "UNKNOWN",
        "exit_code": None,
    }


def fast_candidate_state(item: dict[str, object]) -> str:
    manifest = item["manifest"]
    assert isinstance(manifest, dict)
    final = float(item["observed_final_time"])
    exit_code = segment_exit_code(Path(item["segment"]))
    if exit_code is not None and exit_code != 0:
        return "failed"
    if exit_code == 0 and math.isfinite(final) and final >= (
        float(manifest.get("target_time", TARGET_TIME)) - TIME_TOLERANCE
    ):
        return "complete"
    if exit_code == 0:
        return "exited_success_partial"
    scheduler = item.get("scheduler_evidence")
    scheduler_state = (
        str(scheduler.get("state"))
        if isinstance(scheduler, dict) and scheduler.get("state")
        else "UNKNOWN"
    )
    if scheduler_state in FAILED_SLURM_STATES:
        return "failed"
    if scheduler_state in ACTIVE_SLURM_STATES:
        return "in_progress"
    if scheduler_state == "COMPLETED":
        return "completed_without_exit_artifact"
    mhd, user = history_paths(Path(item["output"]))
    if manifest.get("job_id") is None:
        if mhd is not None or user is not None:
            return "unsubmitted_partial"
        return "prepared"
    if mhd is not None or user is not None:
        return "submitted_unmarked"
    return "submitted_or_pending"


def exact_command_line_overrides(manifest: dict[str, object]) -> tuple[str, ...]:
    """Return exact effective override order, preserving repeated assignments."""

    overrides = manifest.get("command_line_overrides")
    if overrides is None:
        return ()
    if not isinstance(overrides, list) or not all(
        isinstance(value, str) for value in overrides
    ):
        raise ReportError("command_line_overrides must be a list of strings")
    return tuple(overrides)


def fast_candidate_configuration(item: dict[str, object]) -> tuple[str, tuple[str, ...]]:
    manifest = item["manifest"]
    assert isinstance(manifest, dict)
    return (
        str(manifest.get("variant") or "standard"),
        exact_command_line_overrides(manifest),
    )


def fast_candidate_summary(item: dict[str, object]) -> dict[str, object]:
    manifest = item["manifest"]
    assert isinstance(manifest, dict)
    final = float(item["observed_final_time"])
    return {
        "segment": str(item["segment"]),
        "sequence": item["sequence"],
        "source_family": item["source_family"],
        "source_root": str(item["source_root"]),
        "variant": manifest.get("variant"),
        "command_line_overrides": manifest.get("command_line_overrides", []),
        "claim_scope": manifest.get("claim_scope"),
        "job_id": manifest.get("job_id"),
        "scheduler_evidence": item.get("scheduler_evidence"),
        "state": fast_candidate_state(item),
        "run_exit_code": segment_exit_code(Path(item["segment"])),
        "observed_final_time": final if math.isfinite(final) else None,
        "restart_link_valid": item.get("restart_link_valid"),
    }


def fast_candidates(
    root: Path, case_id: str, case_name: str,
    rejections: list[dict[str, object]] | None = None,
) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    seen: set[Path] = set()
    rejected = rejections if rejections is not None else []
    for relative in FAST_RUNS_RELATIVES:
        source_root = root / relative
        if not source_root.is_dir():
            continue
        manifests = source_root.glob(f"**/{case_id}/fast_s*/manifest/fast_run.json")
        for manifest_path in sorted(manifests):
            manifest_path = manifest_path.absolute()
            if manifest_path in seen:
                continue
            seen.add(manifest_path)
            segment = manifest_path.parent.parent
            match = FAST_SEGMENT.match(segment.name)
            if not segment.is_dir() or match is None:
                rejected.append({
                    "reason": "invalid_fast_segment_path",
                    "manifest": str(manifest_path),
                })
                continue
            try:
                manifest = load_json(manifest_path)
            except (OSError, ValueError, json.JSONDecodeError, ReportError) as error:
                rejected.append({
                    "reason": "unreadable_fast_manifest",
                    "manifest": str(manifest_path),
                    "error": str(error),
                })
                continue
            sequence_value = manifest.get("sequence")
            if not isinstance(sequence_value, int) or isinstance(sequence_value, bool):
                rejected.append({
                    "reason": "manifest_sequence_is_not_an_integer",
                    "manifest": str(manifest_path),
                    "observed_sequence": sequence_value,
                })
                continue
            sequence = sequence_value
            path_sequence = int(match.group(1))
            if sequence != path_sequence:
                rejected.append({
                    "reason": "manifest_sequence_mismatch",
                    "manifest": str(manifest_path),
                    "observed_sequence": sequence,
                    "expected_sequence": path_sequence,
                })
                continue
            if str(manifest.get("case_id")) != case_id:
                rejected.append({
                    "reason": "manifest_case_id_mismatch",
                    "manifest": str(manifest_path),
                    "observed_case_id": manifest.get("case_id"),
                    "expected_case_id": case_id,
                })
                continue
            if str(manifest.get("case_name")) != case_name:
                rejected.append({
                    "reason": "manifest_case_name_mismatch",
                    "manifest": str(manifest_path),
                    "observed_case_name": manifest.get("case_name"),
                    "expected_case_name": case_name,
                })
                continue
            try:
                exact_command_line_overrides(manifest)
            except ReportError as error:
                rejected.append({
                    "reason": "manifest_command_line_overrides_invalid",
                    "manifest": str(manifest_path),
                    "error": str(error),
                })
                continue
            missing_identities = [
                name for name in (
                    "input_sha256", "matrix_sha256", "executable_sha256"
                )
                if not isinstance(manifest.get(name), str)
                or re.fullmatch(r"[0-9a-f]{64}", str(manifest.get(name))) is None
            ]
            if missing_identities:
                rejected.append({
                    "reason": "manifest_execution_identity_missing_or_invalid",
                    "manifest": str(manifest_path),
                    "fields": missing_identities,
                })
                continue
            declared_run = Path(str(manifest.get("run_dir", segment))).absolute()
            declared_output = Path(
                str(manifest.get("output_dir", segment / "output"))
            ).absolute()
            if declared_run != segment.absolute() or declared_output != (
                segment / "output"
            ).absolute():
                rejected.append({
                    "reason": "manifest_run_or_output_path_mismatch",
                    "manifest": str(manifest_path),
                    "declared_run_dir": str(declared_run),
                    "expected_run_dir": str(segment.absolute()),
                    "declared_output_dir": str(declared_output),
                    "expected_output_dir": str((segment / "output").absolute()),
                })
                continue
            output = declared_output
            final = peek_history_final(output)
            job_id = manifest.get("job_id")
            candidates.append({
                "segment": segment,
                "manifest_path": manifest_path,
                "manifest": manifest,
                "sequence": sequence,
                "output": output,
                "observed_final_time": final,
                "source_root": source_root.absolute(),
                "source_family": (
                    "original" if relative == FAST_RUNS_RELATIVE
                    else "race" if relative == RACE_RUNS_RELATIVE
                    else "relaxed"
                ),
                "scheduler_evidence": (
                    slurm_job_evidence(str(job_id))
                    if isinstance(job_id, str) and job_id else None
                ),
                "restart_link_valid": sequence == 0 and not bool(manifest.get("restart")),
            })
    return candidates


def fast_restart_time(path: Path) -> float:
    with path.open("rb") as stream:
        payload = stream.read(128 * 1024)
    end = payload.find(b"<par_end>")
    if end < 0:
        raise ReportError(f"restart lacks parameter terminator: {path}")
    text = payload[:end].decode("utf-8")
    block = ""
    values: list[float] = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = (item.strip() for item in line.split("=", 1))
            if key in {"time", "restart_time"}:
                values.append(float(value))
    if len(values) != 1 or not math.isfinite(values[0]):
        raise ReportError(f"restart has ambiguous physical time: {path}")
    return values[0]


def validate_fast_restart_link(
    parent: dict[str, object], child: dict[str, object], restart: Path
) -> str | None:
    parent_manifest = parent["manifest"]
    child_manifest = child["manifest"]
    assert isinstance(parent_manifest, dict)
    assert isinstance(child_manifest, dict)
    if int(child["sequence"]) != int(parent["sequence"]) + 1:
        return "child sequence is not parent sequence plus one"
    if not restart.is_file():
        return f"restart rank-zero file is missing: {restart}"
    expected_sha = child_manifest.get("restart_sha256")
    if not isinstance(expected_sha, str) or not expected_sha:
        return "child manifest has no restart_sha256"
    if sha256_file(restart) != expected_sha:
        return f"restart rank-zero checksum differs: {restart}"
    try:
        restart_time = fast_restart_time(restart)
    except (OSError, UnicodeDecodeError, ValueError, ReportError) as error:
        return str(error)
    child_start = float(child_manifest.get("start_time", math.nan))
    parent_final = float(parent["observed_final_time"])
    if not math.isclose(
        child_start, restart_time, rel_tol=0.0, abs_tol=RESTART_TIME_TOLERANCE
    ):
        return (
            f"child start_time {child_start:.16g} differs from restart "
            f"time {restart_time:.16g}"
        )
    if not math.isfinite(parent_final) or not math.isclose(
        parent_final, restart_time, rel_tol=0.0, abs_tol=RESTART_TIME_TOLERANCE
    ):
        return (
            f"parent history final time {parent_final:.16g} differs from restart "
            f"time {restart_time:.16g}"
        )
    return None


def validate_fast_seed_restart(child: dict[str, object], restart: Path) -> str | None:
    child_manifest = child["manifest"]
    assert isinstance(child_manifest, dict)
    if not restart.is_file():
        return f"historical seed restart rank-zero file is missing: {restart}"
    expected_sha = child_manifest.get("restart_sha256")
    if not isinstance(expected_sha, str) or not expected_sha:
        return "historical seed manifest has no restart_sha256"
    if sha256_file(restart) != expected_sha:
        return f"historical seed restart rank-zero checksum differs: {restart}"
    try:
        restart_time = fast_restart_time(restart)
    except (OSError, UnicodeDecodeError, ValueError, ReportError) as error:
        return str(error)
    child_start = float(child_manifest.get("start_time", math.nan))
    if not math.isclose(
        child_start, restart_time, rel_tol=0.0, abs_tol=RESTART_TIME_TOLERANCE
    ):
        return (
            f"historical seed start_time {child_start:.16g} differs from restart "
            f"time {restart_time:.16g}"
        )
    return None


def select_fast_lineage(
    root: Path, case_id: str, case_name: str
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    rejections: list[dict[str, object]] = []
    candidates = fast_candidates(root, case_id, case_name, rejections)
    if not candidates:
        return [], rejections, []
    warnings: list[str] = []
    parent_by_path: dict[Path, dict[str, object]] = {}
    children: dict[Path, list[Path]] = {
        Path(item["segment"]): [] for item in candidates
    }
    for child in candidates:
        restart_value = child["manifest"].get("restart")  # type: ignore[union-attr]
        if not isinstance(restart_value, str) or not restart_value:
            continue
        restart = Path(restart_value).absolute()
        if int(child["sequence"]) == 0:
            invalid = validate_fast_seed_restart(child, restart)
            if invalid is not None:
                warnings.append(
                    f"invalid historical seed restart for {child['segment']}: {invalid}"
                )
            else:
                child["restart_link_valid"] = True
            continue
        parents = [
            parent for parent in candidates
            if safe_relative_to(restart, Path(parent["output"]) / "rst")
            and fast_candidate_configuration(parent) == fast_candidate_configuration(child)
        ]
        if len(parents) == 1:
            parent = parents[0]
            invalid = validate_fast_restart_link(parent, child, restart)
            if invalid is not None:
                warnings.append(
                    f"invalid fast restart link for {child['segment']}: {invalid}"
                )
                continue
            child["restart_link_valid"] = True
            parent_by_path[Path(child["segment"])] = parent
            children[Path(parent["segment"])].append(Path(child["segment"]))
        else:
            warnings.append(
                f"fast restart for {child['segment']} matched {len(parents)} parent(s)"
            )
    terminals = [
        item for item in candidates if not children[Path(item["segment"])]
    ]
    def lineage_for_terminal(terminal: dict[str, object]) -> list[dict[str, object]]:
        lineage: list[dict[str, object]] = []
        seen: set[Path] = set()
        current = terminal
        while True:
            path = Path(current["segment"])
            if path in seen:
                warnings.append(f"cycle detected in fast lineage at {path}")
                break
            seen.add(path)
            lineage.append(current)
            if path not in parent_by_path:
                break
            current = parent_by_path[path]
        lineage.reverse()
        return lineage

    lineages = [lineage_for_terminal(terminal) for terminal in terminals]

    def lineage_score(lineage: list[dict[str, object]]) -> tuple[object, ...]:
        terminal = lineage[-1]
        manifest = terminal["manifest"]
        assert isinstance(manifest, dict)
        final = float(terminal["observed_final_time"])
        target = float(manifest.get("target_time", TARGET_TIME))
        state = fast_candidate_state(terminal)
        variant, overrides = fast_candidate_configuration(terminal)
        allowed_variant = (
            variant == "standard"
            or (
                case_id in NONFATAL_HARD_BOUND_CASES
                and variant == NONFATAL_HARD_BOUND_VARIANT
            )
        )
        variant_priority = (
            2
            if case_id in NONFATAL_HARD_BOUND_CASES
            and variant == NONFATAL_HARD_BOUND_VARIANT
            else int(allowed_variant)
        )
        mhd, user = history_paths(Path(terminal["output"]))
        synchronized_final = False
        if mhd is not None and user is not None:
            mhd_final = peek_history_final(Path(terminal["output"]))
            try:
                user_record, _ = read_history_source(user)
                user_final = float(user_record["rows"][-1][0])  # type: ignore[index]
                synchronized_final = math.isclose(
                    mhd_final, user_final, rel_tol=0.0, abs_tol=TIME_TOLERANCE
                )
            except (OSError, ReportError, ValueError, IndexError):
                synchronized_final = False
        return (
            int(all(bool(item.get("restart_link_valid")) for item in lineage)),
            variant_priority,
            int(math.isfinite(final) and final >= target - TIME_TOLERANCE),
            {
                "in_progress": 5,
                "submitted_or_pending": 4,
                "submitted_unmarked": 3,
                "exited_success_partial": 2,
                "completed_without_exit_artifact": 2,
                "unsubmitted_partial": 1,
                "prepared": 1,
                "failed": 0,
            }.get(state, 0),
            int(synchronized_final),
            final,
            int(terminal["sequence"]),
            variant,
            overrides,
            str(terminal["segment"]),
        )

    eligible = [
        lineage for lineage in lineages
        if all(bool(item.get("restart_link_valid")) for item in lineage)
    ]
    selected = max(eligible, key=lineage_score) if eligible else []
    selected_paths = {Path(item["segment"]) for item in selected}
    unselected = [
        {"reason": "rejected_manifest_candidate", **item} for item in rejections
    ]
    for lineage in lineages:
        if selected and lineage is selected:
            continue
        unselected.append({
            "reason": (
                "lower_ranked_restart_linked_lineage"
                if all(bool(item.get("restart_link_valid")) for item in lineage)
                else "invalid_restart_linked_lineage"
            ),
            "selection_score": list(lineage_score(lineage)),
            "terminal": fast_candidate_summary(lineage[-1]),
            "segments": [fast_candidate_summary(item) for item in lineage],
        })
    orphaned = [
        item for item in candidates if Path(item["segment"]) not in selected_paths
        and all(Path(item["segment"]) not in {
            Path(member["segment"]) for member in lineage
        } for lineage in lineages)
    ]
    for item in orphaned:
        unselected.append({
            "reason": "orphaned_fast_segment_candidate",
            "selection_score": None,
            "terminal": fast_candidate_summary(item),
            "segments": [fast_candidate_summary(item)],
        })
    if unselected:
        warnings.append(
            f"{case_id} has {len(unselected)} unselected fast lineage(s)"
        )
    if not selected:
        warnings.append(f"{case_id} has no provenance-valid fast lineage")
    return selected, unselected, warnings


def historical_chain_from_restart(restart: Path) -> tuple[list[Path], list[str]]:
    warnings: list[str] = []
    try:
        segment = restart.parent.parent.parent.parent
    except IndexError:
        return [], [f"cannot locate historical segment for seed restart: {restart}"]
    manifest = segment / "manifest/prepared_run.json"
    if not manifest.is_file():
        return [], [f"seed restart has no historical prepared manifest: {restart}"]
    chain: list[Path] = []
    seen: set[Path] = set()
    current = manifest
    while current.is_file():
        current = current.absolute()
        if current in seen:
            warnings.append(f"cycle detected in historical parent manifests: {current}")
            break
        seen.add(current)
        chain.append(current)
        try:
            data = load_json(current)
        except (OSError, ReportError, json.JSONDecodeError) as error:
            warnings.append(str(error))
            break
        command = data.get("command")
        parent = command.get("parent_segment") if isinstance(command, dict) else None
        parent_manifest = parent.get("manifest") if isinstance(parent, dict) else None
        if not isinstance(parent_manifest, str) or not parent_manifest:
            break
        current = Path(parent_manifest)
    chain.reverse()
    return chain, warnings


def segment_exit_code(segment: Path) -> int | None:
    path = segment / "manifest/run_exit_code"
    if not path.is_file():
        return None
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except ValueError:
        return None


def fast_segment_record(item: dict[str, object], order: int) -> dict[str, object]:
    segment = Path(item["segment"])
    manifest = item["manifest"]
    assert isinstance(manifest, dict)
    output = Path(item["output"])
    mhd, user = history_paths(output)
    exit_code = segment_exit_code(segment)
    analysis_path = segment / "manifest/fast_analysis.json"
    final = float(item["observed_final_time"])
    state = fast_candidate_state(item)
    return {
        "kind": "fast",
        "order": order,
        "segment": segment.name,
        "segment_dir": str(segment.absolute()),
        "source_family": item["source_family"],
        "source_root": str(item["source_root"]),
        "manifest": artifact_binding(Path(item["manifest_path"])),
        "start_time": manifest.get("start_time"),
        "target_time": manifest.get("target_time"),
        "observed_final_time": final if math.isfinite(final) else None,
        "state": state,
        "job_id": manifest.get("job_id"),
        "scheduler_evidence": item.get("scheduler_evidence"),
        "run_exit_code": exit_code,
        "fast_analysis": (
            artifact_binding(analysis_path) if analysis_path.is_file() else None
        ),
        "restart": manifest.get("restart"),
        "restart_sha256": manifest.get("restart_sha256"),
        "case_id": manifest.get("case_id"),
        "case_name": manifest.get("case_name"),
        "input": manifest.get("input"),
        "input_sha256": manifest.get("input_sha256"),
        "matrix_sha256": manifest.get("matrix_sha256"),
        "executable": manifest.get("executable"),
        "executable_sha256": manifest.get("executable_sha256"),
        "fast_script": manifest.get("fast_script"),
        "fast_script_sha256": manifest.get("fast_script_sha256"),
        "nodes": int(manifest.get("nodes", 0)),
        "ranks": int(manifest.get("ranks", 0)),
        "ranks_per_node": int(manifest.get("ranks_per_node", 0)),
        "variant": manifest.get("variant"),
        "command_line_overrides": manifest.get("command_line_overrides", []),
        "claim_scope": manifest.get("claim_scope"),
        "prepared_utc": manifest.get("prepared_utc"),
        "submitted_utc": manifest.get("submitted_utc"),
        "run_dir": manifest.get("run_dir"),
        "run_basename": manifest.get("run_basename"),
        "run_sbatch": (
            artifact_binding(segment / "manifest/run.sbatch")
            if (segment / "manifest/run.sbatch").is_file() else None
        ),
        "run_environment": (
            artifact_binding(segment / "manifest/run_environment.txt")
            if (segment / "manifest/run_environment.txt").is_file() else None
        ),
        "run_exit_code_artifact": (
            artifact_binding(segment / "manifest/run_exit_code")
            if (segment / "manifest/run_exit_code").is_file() else None
        ),
        "output": str(output.absolute()),
        "mhd_history": str(mhd.absolute()) if mhd is not None else None,
        "user_history": str(user.absolute()) if user is not None else None,
    }


def historical_segment_record(manifest_path: Path, order: int) -> dict[str, object]:
    manifest = load_json(manifest_path)
    segment = manifest_path.parent.parent
    paths = manifest.get("paths")
    output_value = paths.get("output_dir") if isinstance(paths, dict) else None
    output = Path(str(output_value)) if output_value else segment / "output"
    mhd, user = history_paths(output)
    allocation = manifest.get("allocation")
    ranks = 0
    if isinstance(allocation, dict):
        if isinstance(allocation.get("ranks"), int):
            ranks = int(allocation["ranks"])
        else:
            ranks = int(allocation.get("nodes", 0)) * int(
                allocation.get("ranks_per_node", 1)
            )
    inspection = manifest.get("scientific_inspection")
    final = (
        float(inspection["final_time"])
        if isinstance(inspection, dict) and isinstance(inspection.get("final_time"), (int, float))
        else peek_history_final(output)
    )
    state = str(manifest.get("state", "unknown"))
    if isinstance(inspection, dict):
        state = (
            "accepted" if inspection.get("accepted")
            else "clean_partial" if inspection.get("clean_for_continuation")
            else state
        )
    command = manifest.get("command")
    return {
        "kind": "historical_seed_prefix",
        "order": order,
        "segment": segment.name,
        "segment_dir": str(segment.absolute()),
        "manifest": artifact_binding(manifest_path),
        "observed_final_time": final if math.isfinite(final) else None,
        "state": state,
        "input": command.get("input_file") if isinstance(command, dict) else None,
        "input_sha256": command.get("input_sha256") if isinstance(command, dict) else None,
        "matrix_sha256": command.get("matrix_sha256") if isinstance(command, dict) else None,
        "executable": command.get("executable") if isinstance(command, dict) else None,
        "executable_sha256": (
            command.get("executable_sha256") if isinstance(command, dict) else None
        ),
        "ranks": ranks,
        "output": str(output.absolute()),
        "mhd_history": str(mhd.absolute()) if mhd is not None else None,
        "user_history": str(user.absolute()) if user is not None else None,
    }


def load_workflow_module() -> object:
    path = REPO_ROOT / "scripts/cgl_lf_workflow.py"
    module_name = "_cgl_lf_fast_report_workflow"
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ReportError(f"cannot load workflow module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        raise
    return module


def model_choices_for_input(
    path: Path, overrides: list[str] | None = None
) -> dict[str, str]:
    workflow = load_workflow_module()
    source = path.read_text(encoding="utf-8")
    effective_overrides = overrides or []
    choices = workflow.model_choices(source, effective_overrides)
    if not isinstance(choices, dict):
        raise ReportError(f"workflow returned invalid model choices: {path}")
    result = {str(key): str(value) for key, value in choices.items()}
    strict = workflow.input_block_value(
        source, "mhd", "cgl_lf_strict_admissibility"
    ) or "false"
    for override in effective_overrides:
        prefix = "mhd/cgl_lf_strict_admissibility="
        if override.startswith(prefix):
            strict = override[len(prefix):]
    result["cgl_lf_strict_admissibility"] = strict
    return result


def case_status(lineage: list[dict[str, object]], final_time: float | None) -> str:
    if not lineage:
        return "not_started"
    terminal = lineage[-1]
    if terminal.get("state") == "failed":
        return "failed_partial"
    if (
        final_time is not None
        and final_time >= TARGET_TIME - TIME_TOLERANCE
        and terminal.get("state") in ("accepted", "complete")
    ):
        return "complete"
    if terminal.get("state") in (
        "in_progress",
        "prepared",
        "submitted_or_pending",
        "submitted_unmarked",
    ):
        return "in_progress"
    return "partial"


def assemble_r02(
    root: Path, output: Path, case: dict[str, object]
) -> dict[str, object]:
    case_id = "R02"
    case_dir = output / "cases" / case_id
    bundle = root / R02_BUNDLE_RELATIVE
    warnings: list[str] = []
    errors: list[str] = []
    if not (bundle / "manifest.json").is_file():
        return {
            "case_id": case_id,
            "case_name": case["name"],
            "status": "not_started",
            "warnings": [],
            "errors": [f"R02 bundle is missing: {bundle}"],
        }
    manifest = load_json(bundle / "manifest.json")
    cases = manifest.get("cases")
    if not isinstance(cases, list) or len(cases) != 1 or not isinstance(cases[0], dict):
        raise ReportError(f"R02 bundle has invalid cases: {bundle}")
    bundle_case = cases[0]
    outputs = bundle_case.get("outputs")
    if not isinstance(outputs, dict):
        raise ReportError(f"R02 bundle case lacks outputs: {bundle}")
    mhd_source = bundle / str(outputs["mhd_history"])
    user_source = bundle / str(outputs["user_history"])
    history_dir = case_dir / "history"
    mhd = merge_histories(
        [("R02 accepted bundle", mhd_source)],
        history_dir / f"{case['name']}.mhd.hst",
        "MHD",
    )
    user = merge_histories(
        [("R02 accepted bundle", user_source)],
        history_dir / f"{case['name']}.user.hst",
        "user",
    )
    warnings.extend(mhd["warnings"])
    warnings.extend(user["warnings"])
    errors.extend(mhd["errors"])
    errors.extend(user["errors"])
    snapshots: list[dict[str, object]] = []
    snapshot_warnings: list[str] = []
    for index, relative in enumerate(outputs.get("snapshot_paths", [])):
        representative = bundle / str(relative)
        group, group_warnings = snapshot_group(
            representative, None, "R02 accepted bundle", index
        )
        snapshots.append(group)
        snapshot_warnings.extend(group_warnings)
    snapshots.sort(key=lambda item: float(item["time"]) if item["time"] is not None else math.inf)
    snapshot_index = {
        "schema_version": 1,
        "generated_utc": utc_now(),
        "snapshot_count": len(snapshots),
        "complete_snapshot_count": sum(bool(item["complete"]) for item in snapshots),
        "time_first": snapshots[0]["time"] if snapshots else None,
        "time_last": snapshots[-1]["time"] if snapshots else None,
        "snapshots": snapshots,
        "duplicates": [],
        "warnings": snapshot_warnings,
    }
    write_json(case_dir / "snapshots.json", snapshot_index)
    final_time = float(mhd["time_final"]) if mhd.get("available") else None
    input_path = bundle / str(bundle_case["execution_input"])
    model = (
        bundle_case.get("model_choices")
        if isinstance(bundle_case.get("model_choices"), dict)
        else model_choices_for_input(input_path)
    )
    status = case_status([{"state": "accepted"}], final_time)
    if errors:
        status = "assembly_error"
    record = {
        "schema_version": 1,
        "assembled_utc": utc_now(),
        "case_id": case_id,
        "case_name": case["name"],
        "matrix_case": case,
        "status": status,
        "target_time": TARGET_TIME,
        "final_time": final_time,
        "model_choices": model,
        "input": artifact_binding(input_path),
        "lineage": [{
            "kind": "accepted_r02_bundle",
            "order": 0,
            "segment": "R02 accepted bundle",
            "state": str(manifest.get("status", "unknown")),
            "manifest": artifact_binding(bundle / "manifest.json"),
            "output": str(bundle.absolute()),
            "ranks": snapshots[0]["expected_ranks"] if snapshots else None,
            "mhd_history": str(mhd_source.absolute()),
            "user_history": str(user_source.absolute()),
        }],
        "selected_fast_lineage": None,
        "unselected_lineages": [],
        "unselected_segments": [],
        "histories": {"mhd": mhd, "user": user},
        "snapshots": {
            "path": str((case_dir / "snapshots.json").absolute()),
            "snapshot_count": snapshot_index["snapshot_count"],
            "complete_snapshot_count": snapshot_index["complete_snapshot_count"],
        },
        "warnings": [*warnings, *snapshot_warnings],
        "errors": errors,
    }
    write_json(case_dir / "lineage.json", record)
    return record


def assemble_fast_case(
    root: Path, frozen_source: Path, output: Path, case_id: str,
    case: dict[str, object],
) -> dict[str, object]:
    case_dir = output / "cases" / case_id
    selected, unselected, warnings = select_fast_lineage(
        root, case_id, str(case["name"])
    )
    errors: list[str] = []
    historical_manifests: list[Path] = []
    if selected:
        first_manifest = selected[0]["manifest"]
        assert isinstance(first_manifest, dict)
        restart = first_manifest.get("restart")
        if isinstance(restart, str) and restart:
            historical_manifests, historical_warnings = historical_chain_from_restart(
                Path(restart)
            )
            warnings.extend(historical_warnings)
    lineage: list[dict[str, object]] = []
    for manifest in historical_manifests:
        try:
            lineage.append(historical_segment_record(manifest, len(lineage)))
        except (OSError, ReportError, json.JSONDecodeError, ValueError) as error:
            errors.append(str(error))
    for item in selected:
        lineage.append(fast_segment_record(item, len(lineage)))
    variants = sorted({
        str(item["variant"]) for item in lineage if item.get("variant")
    })
    claim_scopes = sorted({
        str(item["claim_scope"]) for item in lineage if item.get("claim_scope")
    })
    command_line_overrides = (
        list(exact_command_line_overrides(selected[-1]["manifest"]))
        if selected and isinstance(selected[-1].get("manifest"), dict)
        else []
    )
    expected_input = frozen_source / str(case["input"])
    try:
        model = model_choices_for_input(expected_input, command_line_overrides)
        input_binding = artifact_binding(expected_input)
    except (OSError, ReportError) as error:
        model = {}
        input_binding = {"path": str(expected_input.absolute()), "available": False}
        errors.append(str(error))
    identities = {
        "input_sha256": sorted({
            str(item["input_sha256"]) for item in lineage if item.get("input_sha256")
        }),
        "matrix_sha256": sorted({
            str(item["matrix_sha256"]) for item in lineage if item.get("matrix_sha256")
        }),
        "executable_sha256": sorted({
            str(item["executable_sha256"])
            for item in lineage if item.get("executable_sha256")
        }),
    }
    for name, values in identities.items():
        if len(values) != 1:
            errors.append(
                f"{case_id} lineage must have exactly one {name} value: {values}"
            )
    expected_input_sha = input_binding.get("sha256")
    selected_input_shas = identities["input_sha256"]
    if expected_input_sha and selected_input_shas and selected_input_shas != [expected_input_sha]:
        errors.append(
            f"{case_id} selected input hash differs from frozen matrix input: "
            f"selected={selected_input_shas}, expected={expected_input_sha}"
        )
    mhd_sources = [
        (str(item["segment"]), Path(str(item["mhd_history"])))
        for item in lineage if item.get("mhd_history")
    ]
    user_sources = [
        (str(item["segment"]), Path(str(item["user_history"])))
        for item in lineage if item.get("user_history")
    ]
    history_dir = case_dir / "history"
    mhd = merge_histories(
        mhd_sources, history_dir / f"{case['name']}.mhd.hst", "MHD"
    )
    user = merge_histories(
        user_sources, history_dir / f"{case['name']}.user.hst", "user"
    )
    warnings.extend(mhd["warnings"])
    warnings.extend(user["warnings"])
    errors.extend(mhd["errors"])
    errors.extend(user["errors"])
    snapshot_sources = [
        (str(item["segment"]), Path(str(item["output"])), int(item.get("ranks") or 0))
        for item in lineage if item.get("output") and int(item.get("ranks") or 0) > 0
    ]
    snapshots = index_snapshots(snapshot_sources, case_dir / "snapshots.json")
    warnings.extend(snapshots["warnings"])
    final_time = float(mhd["time_final"]) if mhd.get("available") else None
    status = case_status(lineage, final_time)
    if errors:
        status = "assembly_error"
    record = {
        "schema_version": 1,
        "assembled_utc": utc_now(),
        "case_id": case_id,
        "case_name": case["name"],
        "matrix_case": case,
        "status": status,
        "target_time": TARGET_TIME,
        "final_time": final_time,
        "model_choices": model,
        "input": input_binding,
        "lineage_identities": identities,
        "lineage_variants": variants,
        "lineage_claim_scopes": claim_scopes,
        "lineage_command_line_overrides": command_line_overrides,
        "lineage": lineage,
        "selected_fast_lineage": (
            {
                "reason": "highest_ranked_restart_linked_lineage",
                "terminal": fast_candidate_summary(selected[-1]),
                "segments": [fast_candidate_summary(item) for item in selected],
            }
            if selected else None
        ),
        "unselected_lineages": unselected,
        "unselected_segments": unselected,
        "histories": {"mhd": mhd, "user": user},
        "snapshots": {
            "path": str((case_dir / "snapshots.json").absolute()),
            "snapshot_count": snapshots["snapshot_count"],
            "complete_snapshot_count": snapshots["complete_snapshot_count"],
        },
        "warnings": warnings,
        "errors": errors,
    }
    write_json(case_dir / "lineage.json", record)
    return record


def command_assemble(args: argparse.Namespace) -> int:
    cases, _ = matrix_cases(args.matrix)
    selected = expand_cases(args.cases)
    records: dict[str, object] = {}
    for case_id in selected:
        case = cases.get(case_id)
        if case is None:
            raise ReportError(f"matrix lacks requested case: {case_id}")
        try:
            if case_id == "R02":
                record = assemble_r02(args.root, args.output, case)
            else:
                record = assemble_fast_case(
                    args.root, args.frozen_source, args.output, case_id, case
                )
        except Exception as error:  # retain other cases when one source is malformed
            record = {
                "schema_version": 1,
                "assembled_utc": utc_now(),
                "case_id": case_id,
                "case_name": case.get("name"),
                "status": "assembly_error",
                "warnings": [],
                "errors": [f"{type(error).__name__}: {error}"],
            }
            write_json(args.output / "cases" / case_id / "lineage.json", record)
        records[case_id] = record
        print(
            f"{case_id}: {record.get('status')} "
            f"t={record.get('final_time')} errors={len(record.get('errors', []))}"
        )
    inventory = {
        "schema_version": 1,
        "assembled_utc": utc_now(),
        "root": str(args.root.absolute()),
        "output": str(args.output.absolute()),
        "matrix": artifact_binding(args.matrix),
        "adapter": artifact_binding(Path(__file__)),
        "searched_fast_run_roots": [
            str((args.root / relative).absolute()) for relative in FAST_RUNS_RELATIVES
        ],
        "cases": records,
    }
    write_json(args.output / "inventory.json", inventory)
    write_json(args.output / "manifest.json", {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_direct_fast_report",
        "created_utc": utc_now(),
        "root": str(args.root.absolute()),
        "output": str(args.output.absolute()),
        "matrix": inventory["matrix"],
        "adapter": inventory["adapter"],
        "commands": [
            "assemble", "analyze-case", "aggregate", "render", "manuscript", "verify"
        ],
    })
    return 0


def load_pure_analyzer() -> object:
    """Load paper-analysis kernels while replacing its validation import with a stub."""

    path = REPO_ROOT / "scripts/analyze_cgl_lf_paper.py"
    module_name = "_cgl_lf_fast_report_analyzer"
    if module_name in sys.modules:
        return sys.modules[module_name]
    dependency_name = "cgl_lf_stage_i_validate_segment"
    previous = sys.modules.get(dependency_name)
    stub = types.ModuleType(dependency_name)
    stub.QUALIFIED_RESTART_BINARY_ABIS = ()
    sys.modules[dependency_name] = stub
    try:
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            raise ReportError(f"cannot load analyzer: {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        if previous is None:
            sys.modules.pop(dependency_name, None)
        else:
            sys.modules[dependency_name] = previous


def snapshot_worker_plan(
    requested_workers: int,
    selected: list[Path],
    exact_rank_sets: dict[str, list[Path]],
    memory_budget_gib: float,
) -> dict[str, object]:
    """Bound independent snapshot workers by affinity and conservative memory use."""

    if requested_workers < 1:
        raise ReportError("--snapshot-workers must be positive")
    if not math.isfinite(memory_budget_gib) or memory_budget_gib <= 0.0:
        raise ReportError("--snapshot-memory-budget-gib must be positive and finite")
    if not selected:
        return {
            "requested_workers": requested_workers,
            "workers": 1,
            "snapshot_count": 0,
            "available_cpus": 1,
            "memory_budget_bytes": int(memory_budget_gib * GIBIBYTE),
            "maximum_snapshot_input_bytes": 0,
            "estimated_peak_bytes_per_process": 0,
        }
    try:
        available_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        available_cpus = os.cpu_count() or 1
    maximum_snapshot_input_bytes = max(
        sum(path.stat().st_size for path in exact_rank_sets[str(snapshot)])
        for snapshot in selected
    )
    estimated_peak_bytes = max(
        1,
        int(math.ceil(
            maximum_snapshot_input_bytes * SNAPSHOT_PEAK_BYTES_PER_INPUT_BYTE
        )),
    )
    memory_budget_bytes = int(memory_budget_gib * GIBIBYTE)
    # The coordinator performs the common-range pass before workers start. Reserve
    # one snapshot-sized process so retained allocator pages cannot exhaust a node.
    memory_limited_workers = max(
        1, memory_budget_bytes // estimated_peak_bytes - 1
    )
    workers = min(
        requested_workers,
        len(selected),
        max(1, available_cpus),
        memory_limited_workers,
    )
    return {
        "requested_workers": requested_workers,
        "workers": max(1, workers),
        "snapshot_count": len(selected),
        "available_cpus": available_cpus,
        "memory_budget_bytes": memory_budget_bytes,
        "maximum_snapshot_input_bytes": maximum_snapshot_input_bytes,
        "estimated_peak_bytes_per_process": estimated_peak_bytes,
    }


def selected_snapshot_inventory(
    analyzer: object,
    paths: list[Path],
    time_start: float | None,
    time_end: float | None,
    expected_ranks_by_path: dict[str, int],
) -> tuple[list[Path], dict[str, list[Path]]]:
    """Select snapshots and freeze their exact rank-local sibling sets."""

    candidate_rank_sets = {
        str(path): analyzer.snapshot_sibling_paths(
            path, expected_ranks_by_path.get(str(path))
        )
        for path in paths
    }
    selected = [
        path
        for path in paths
        if analyzer.time_mask(
            analyzer.np.asarray([analyzer.snapshot_time(path)]),
            time_start,
            time_end,
        )[0]
    ]
    if len({str(path) for path in selected}) != len(selected):
        raise ValueError("selected snapshot paths must be unique")
    return selected, {
        str(path): candidate_rank_sets[str(path)] for path in selected
    }


def read_selected_snapshot(
    analyzer: object, path: Path, exact_rank_set: list[Path]
) -> tuple[dict[str, object], tuple[float, float, float], float]:
    """Read one shared or exact rank-local snapshot using analyzer semantics."""

    if path.parent.name == "rank_00000000":
        return analyzer.read_snapshot(path, exact_rank_set)
    return analyzer.read_snapshot(path)


def snapshot_range_record(
    analyzer: object, path: Path, exact_rank_set: list[Path]
) -> dict[str, object]:
    """Return exact local PDF extrema for one authenticated snapshot."""

    fields, lengths, _ = read_selected_snapshot(analyzer, path, exact_rank_set)
    extrema: dict[str, list[float]] = {}
    joint_extrema: dict[str, list[list[float]]] = {}
    for name, values in analyzer.pdf_fields(fields, lengths).items():
        finite = values[analyzer.np.isfinite(values)]
        if finite.size == 0:
            continue
        extrema[name] = [
            float(analyzer.np.min(finite)),
            float(analyzer.np.max(finite)),
        ]
    for name, (x_values, y_values) in analyzer.pressure_density_fields(fields).items():
        joint_extrema[name] = [
            [float(analyzer.np.min(values)), float(analyzer.np.max(values))]
            for values in (x_values, y_values)
        ]
    del fields
    gc.collect()
    return {"extrema": extrema, "joint_extrema": joint_extrema}


def merge_snapshot_range_records(
    completed: list[tuple[str, dict[str, object]]],
) -> tuple[
    dict[str, tuple[float, float]],
    dict[str, tuple[tuple[float, float], tuple[float, float]]],
]:
    """Reduce per-snapshot extrema in submitted snapshot order."""

    extrema: dict[str, list[float]] = {}
    joint_extrema: dict[str, list[list[float]]] = {}
    for _path, record in completed:
        local_extrema = record["extrema"]
        local_joint_extrema = record["joint_extrema"]
        assert isinstance(local_extrema, dict)
        assert isinstance(local_joint_extrema, dict)
        for name, bounds in local_extrema.items():
            low, high = float(bounds[0]), float(bounds[1])
            if name in extrema:
                extrema[name][0] = min(extrema[name][0], low)
                extrema[name][1] = max(extrema[name][1], high)
            else:
                extrema[name] = [low, high]
        for name, coordinates in local_joint_extrema.items():
            if name not in joint_extrema:
                joint_extrema[name] = [
                    [float(bounds[0]), float(bounds[1])] for bounds in coordinates
                ]
                continue
            for index, bounds in enumerate(coordinates):
                joint_extrema[name][index][0] = min(
                    joint_extrema[name][index][0], float(bounds[0])
                )
                joint_extrema[name][index][1] = max(
                    joint_extrema[name][index][1], float(bounds[1])
                )
    ranges: dict[str, tuple[float, float]] = {}
    for name, (low, high) in extrema.items():
        if high <= low:
            delta = max(abs(low), 1.0) * 1.0e-12
            low -= delta
            high += delta
        ranges[name] = (low, high)
    joint_ranges: dict[
        str, tuple[tuple[float, float], tuple[float, float]]
    ] = {}
    for name, coordinates in joint_extrema.items():
        padded: list[tuple[float, float]] = []
        for low, high in coordinates:
            if high <= low:
                delta = max(abs(low), 1.0) * 1.0e-12
                low -= delta
                high += delta
            padded.append((low, high))
        joint_ranges[name] = (padded[0], padded[1])
    return ranges, joint_ranges


def common_snapshot_ranges(
    analyzer: object,
    selected: list[Path],
    exact_rank_sets: dict[str, list[Path]],
) -> tuple[
    dict[str, tuple[float, float]],
    dict[str, tuple[tuple[float, float], tuple[float, float]]],
]:
    """Discover shared histogram ranges in original deterministic snapshot order."""

    return merge_snapshot_range_records([
        (
            str(path),
            snapshot_range_record(analyzer, path, exact_rank_sets[str(path)]),
        )
        for path in selected
    ])


def configure_snapshot_worker() -> None:
    """Prevent nested BLAS/OpenMP pools from oversubscribing snapshot processes."""

    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def analyze_snapshot_worker(task: dict[str, object]) -> tuple[str, dict[str, object]]:
    """Analyze one authenticated snapshot in an isolated bounded-memory process."""

    analyzer = load_pure_analyzer()
    path = Path(str(task["path"]))
    exact_rank_set = [Path(str(value)) for value in task["exact_rank_set"]]
    fields, lengths, time = read_selected_snapshot(analyzer, path, exact_rank_set)
    record = analyzer.analyze_fields(
        fields,
        lengths,
        time,
        int(task["bins"]),
        [int(value) for value in task["alignment_shells"]],
        task["ranges"],
        task["model_choices"],
        int(task["eddy_samples"]),
        int(task["eddy_bins"]),
        int(task["eddy_seed"]),
        task["joint_ranges"],
    )
    del fields
    gc.collect()
    return str(path), record


def snapshot_range_worker(task: dict[str, object]) -> tuple[str, dict[str, object]]:
    """Scan one authenticated snapshot for common-range extrema."""

    analyzer = load_pure_analyzer()
    path = Path(str(task["path"]))
    exact_rank_set = [Path(str(value)) for value in task["exact_rank_set"]]
    return str(path), snapshot_range_record(analyzer, path, exact_rank_set)


def run_snapshot_workers(
    tasks: list[dict[str, object]], workers: int
) -> list[tuple[str, dict[str, object]]]:
    """Run one task per child lifetime and preserve submitted task order."""

    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=configure_snapshot_worker,
        max_tasks_per_child=1,
    ) as executor:
        return list(executor.map(analyze_snapshot_worker, tasks, chunksize=1))


def run_snapshot_range_workers(
    tasks: list[dict[str, object]], workers: int
) -> list[tuple[str, dict[str, object]]]:
    """Run one common-range task per child lifetime in submitted order."""

    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=configure_snapshot_worker,
        max_tasks_per_child=1,
    ) as executor:
        return list(executor.map(snapshot_range_worker, tasks, chunksize=1))


def analyze_snapshot_paths_bounded(
    analyzer: object,
    paths: list[Path],
    bins: int,
    alignment_shells: list[int],
    time_start: float | None,
    time_end: float | None,
    model_choices: dict[str, object] | None,
    eddy_samples: int,
    eddy_bins: int,
    eddy_seed: int,
    expected_ranks_by_path: dict[str, int],
    requested_workers: int,
    memory_budget_gib: float,
    worker_runner=run_snapshot_workers,
    range_runner=run_snapshot_range_workers,
) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    """Analyze independent snapshots concurrently without changing reductions."""

    selected, exact_rank_sets = selected_snapshot_inventory(
        analyzer, paths, time_start, time_end, expected_ranks_by_path
    )
    plan = snapshot_worker_plan(
        requested_workers, selected, exact_rank_sets, memory_budget_gib
    )
    workers = int(plan["workers"])
    print(
        "snapshot analysis plan: "
        f"snapshots={len(selected)} workers={workers}/{requested_workers} "
        f"max_input_gib={float(plan['maximum_snapshot_input_bytes']) / GIBIBYTE:.3f} "
        f"estimated_peak_gib_per_process="
        f"{float(plan['estimated_peak_bytes_per_process']) / GIBIBYTE:.3f} "
        f"memory_budget_gib={memory_budget_gib:g}"
    )
    if workers == 1:
        return analyzer.analyze_snapshot_paths(
            paths,
            bins,
            alignment_shells,
            time_start,
            time_end,
            model_choices,
            eddy_samples,
            eddy_bins,
            eddy_seed,
            expected_ranks_by_path,
        )
    snapshot_provenance = {
        str(path): analyzer.snapshot_digest_provenance(
            path,
            expected_ranks_by_path.get(str(path)),
            exact_rank_sets[str(path)],
        )
        for path in selected
    }
    range_tasks = [
        {
            "path": str(path),
            "exact_rank_set": [str(value) for value in exact_rank_sets[str(path)]],
        }
        for path in selected
    ]
    range_completed = range_runner(range_tasks, workers)
    if [path for path, _record in range_completed] != [
        str(path) for path in selected
    ]:
        raise ValueError("parallel snapshot range scan returned snapshots out of order")
    ranges, joint_ranges = merge_snapshot_range_records(range_completed)
    tasks = [
        {
            "path": str(path),
            "exact_rank_set": [str(value) for value in exact_rank_sets[str(path)]],
            "bins": bins,
            "alignment_shells": alignment_shells,
            "ranges": ranges,
            "model_choices": model_choices,
            "eddy_samples": eddy_samples,
            "eddy_bins": eddy_bins,
            "eddy_seed": eddy_seed,
            "joint_ranges": joint_ranges,
        }
        for path in selected
    ]
    completed = worker_runner(tasks, workers)
    if [path for path, _record in completed] != [str(path) for path in selected]:
        raise ValueError("parallel snapshot analysis returned snapshots out of order")
    records: dict[str, dict[str, object]] = {}
    for path, record in completed:
        record["snapshot_provenance"] = snapshot_provenance[path]
        records[path] = record
    for path in selected:
        if analyzer.snapshot_digest_provenance(
            path,
            expected_ranks_by_path.get(str(path)),
            exact_rank_sets[str(path)],
        ) != snapshot_provenance[str(path)]:
            raise ValueError(f"snapshot changed while being analyzed: {path}")
    ensemble = analyzer.average_snapshot_records(records)
    ensemble["time_start"] = time_start
    ensemble["time_end"] = time_end
    if "firehose_threshold_occupancy" in ensemble:
        ensemble["firehose_threshold_occupancy"]["analysis_window"].update({
            "requested_time_start": time_start,
            "requested_time_end": time_end,
        })
    return records, ensemble


def parsed_history_columns(path: Path) -> dict[str, list[float]]:
    record, _ = read_history_source(path)
    labels = list(record["labels"])  # type: ignore[arg-type]
    return {
        label: [float(row[index]) for row in record["rows"]]  # type: ignore[union-attr]
        for index, label in enumerate(labels)
    }


def relative_mass_drift(values: list[float]) -> float:
    if not values:
        return math.inf
    reference = values[0]
    return max(abs(value - reference) for value in values) / max(abs(reference), 1.0)


def max_parallel_forcing_fraction(user: dict[str, list[float]]) -> float | None:
    if "force_prp2" not in user or "force_prl2" not in user:
        return None
    return max(
        parallel / max(parallel + perpendicular, 1.0e-300)
        for perpendicular, parallel in zip(user["force_prp2"], user["force_prl2"])
    )


def compute_health(
    lineage: dict[str, object], mhd: dict[str, list[float]],
    user: dict[str, list[float]], model: dict[str, object],
) -> dict[str, object]:
    structural_warnings = list(lineage.get("warnings", []))
    structural_errors = list(lineage.get("errors", []))
    science_warnings: list[str] = []
    numerical_warnings: list[str] = []
    mhd_times = mhd.get("time", [])
    user_times = user.get("time", [])
    synchronized = len(mhd_times) == len(user_times) and all(
        math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-12)
        for left, right in zip(mhd_times, user_times)
    )
    if not synchronized:
        structural_errors.append("merged MHD and user history times are not synchronized")
    strict_maxima = {
        name: max((abs(value) for value in mhd.get(name, [])), default=math.inf)
        for name in STRICT_FAILURE_COLUMNS
    }
    fatal_maxima = {
        name: strict_maxima[name] for name in FATAL_FAILURE_COLUMNS
    }
    if any(value != 0.0 for value in fatal_maxima.values()):
        numerical_warnings.append(f"fatal LF failure counters are nonzero: {fatal_maxima}")
    variants = lineage.get("lineage_variants", [])
    nonfatal_hard_bound = (
        isinstance(variants, list) and NONFATAL_HARD_BOUND_VARIANT in variants
    )
    hard_bound_maximum = strict_maxima["lf_hardbd"]
    if hard_bound_maximum != 0.0:
        if nonfatal_hard_bound:
            case_id = str(lineage.get("case_id", "finite-limiter"))
            science_warnings.append(
                f"nonfatal {case_id} diagnostic retained hard-bound events: "
                f"lf_hardbd={hard_bound_maximum:.6g}"
            )
        else:
            numerical_warnings.append(
                f"strict LF hard-bound counter is nonzero: {hard_bound_maximum:.6g}"
            )
    mhd_mass_drift = relative_mass_drift(mhd.get("mass", []))
    user_mass_drift = relative_mass_drift(user.get("mass", []))
    if mhd_mass_drift > 1.0e-8 or user_mass_drift > 1.0e-8:
        numerical_warnings.append(
            f"mass drift exceeds 1e-8: mhd={mhd_mass_drift:.6g}, "
            f"user={user_mass_drift:.6g}"
        )
    mass_mismatch = math.inf
    if synchronized and mhd.get("mass") and user.get("mass"):
        scale = max(abs(mhd["mass"][0]), 1.0)
        mass_mismatch = max(
            abs(left - right) for left, right in zip(mhd["mass"], user["mass"])
        ) / scale
        if mass_mismatch > 1.0e-8:
            numerical_warnings.append(
                f"MHD/user mass mismatch exceeds 1e-8: {mass_mismatch:.6g}"
            )
    monotonic_violations: dict[str, int] = {}
    for name in MONOTONIC_LF_COUNTERS:
        values = mhd.get(name)
        if values is None:
            continue
        count = sum(
            right < left - max(abs(left), abs(right), 1.0) * 1.0e-13
            for left, right in zip(values, values[1:])
        )
        if count:
            monotonic_violations[name] = count
    if monotonic_violations:
        numerical_warnings.append(
            f"LF cumulative counters decreased: {monotonic_violations}"
        )
    hard_volume_max = max(user.get("hard_vol", [0.0]))
    if hard_volume_max != 0.0:
        science_warnings.append(f"hard-bound volume is nonzero: {hard_volume_max:.6g}")
    parallel_fraction = max_parallel_forcing_fraction(user)
    forcing_mode = str(model.get("forcing_mode", "unspecified"))
    if parallel_fraction is not None:
        if forcing_mode == "alfvenic_z_perpendicular" and parallel_fraction > 1.0e-10:
            science_warnings.append(
                f"Alfvenic forcing has parallel fraction {parallel_fraction:.6g}"
            )
        if forcing_mode == "isotropic_random" and parallel_fraction <= 0.05:
            science_warnings.append(
                f"random forcing parallel fraction is only {parallel_fraction:.6g}"
            )
    hardwall = str(model.get("limiter_hardwall", "false")).lower() in ("true", "1")
    hwproj = mhd.get("lf_hwproj", [])
    hw_increment = hwproj[-1] - hwproj[0] if len(hwproj) >= 2 else None
    complete = lineage.get("status") == "complete"
    if hardwall and complete and hw_increment is not None and hw_increment <= 0.0:
        science_warnings.append("hard-wall case has no accumulated hard-wall projections")
    if not hardwall and hw_increment is not None and hw_increment != 0.0:
        science_warnings.append(
            f"finite-limiter case has nonzero hard-wall projection increment {hw_increment}"
        )
    final_time = mhd_times[-1] if mhd_times else None
    return {
        "schema_version": 1,
        "case_status": lineage.get("status"),
        "final_time": final_time,
        "target_time": TARGET_TIME,
        "complete": final_time is not None and final_time >= TARGET_TIME - TIME_TOLERANCE,
        "history_rows": {"mhd": len(mhd_times), "user": len(user_times)},
        "histories_synchronized": synchronized,
        "strict_lf_failure_maxima": strict_maxima,
        "fatal_lf_failure_maxima": fatal_maxima,
        "hard_bound_diagnostic_maximum": hard_bound_maximum,
        "nonfatal_hard_bound_variant": nonfatal_hard_bound,
        "mass_relative_drift": {"mhd": mhd_mass_drift, "user": user_mass_drift},
        "mhd_user_mass_relative_mismatch": mass_mismatch,
        "lf_counter_monotonicity_violations": monotonic_violations,
        "hard_bound_volume_maximum": hard_volume_max,
        "parallel_forcing_fraction_maximum": parallel_fraction,
        "hardwall_projection_increment": hw_increment,
        "structural_errors": structural_errors,
        "structural_warnings": structural_warnings,
        "numerical_warnings": numerical_warnings,
        "science_warnings": science_warnings,
        "result": (
            "structural_error" if structural_errors
            else "warnings" if numerical_warnings or science_warnings or structural_warnings
            else "clean"
        ),
    }


def window_summaries(
    analyzer: object, user_path: Path, mhd_path: Path, model: dict[str, object]
) -> dict[str, object]:
    result: dict[str, object] = {}
    for name, (start, end) in WINDOWS.items():
        result[name] = {
            "time_start": start,
            "time_end": end,
            "history": analyzer.summarize_history(user_path, start, end),
            "lf_history": analyzer.summarize_lf_history(mhd_path, start, end),
            "forcing_energy_budget": analyzer.summarize_forcing_energy_budget(
                user_path, mhd_path, model, start, end
            ),
        }
    return result


def case_report_markdown(diagnostics: dict[str, object]) -> str:
    health = diagnostics["health"]
    assert isinstance(health, dict)
    lines = [
        f"# {diagnostics['case_id']}: {diagnostics['case_name']}",
        "",
        f"- Assembly status: `{diagnostics['assembly_status']}`",
        f"- Analysis status: `{diagnostics['analysis_status']}`",
        f"- Final time: `{health.get('final_time')}`",
        f"- Numerical/health result: `{health.get('result')}`",
        f"- Retained late snapshots: `{diagnostics.get('selected_snapshot_count', 0)}`",
        "",
        "## Warnings",
        "",
    ]
    warnings = [
        *health.get("structural_warnings", []),
        *health.get("numerical_warnings", []),
        *health.get("science_warnings", []),
        *diagnostics.get("analysis_warnings", []),
    ]
    lines.extend(f"- {value}" for value in warnings)
    if not warnings:
        lines.append("- None.")
    errors = [*health.get("structural_errors", []), *diagnostics.get("analysis_errors", [])]
    lines.extend(["", "## Structural Errors", ""])
    lines.extend(f"- {value}" for value in errors)
    if not errors:
        lines.append("- None.")
    lines.extend([
        "",
        "## Interpretation Boundary",
        "",
        "Unexpected scientific behavior is retained for review and does not "
        "invalidate or suppress the available analysis products.",
    ])
    return "\n".join(lines) + "\n"


def command_analyze_case(args: argparse.Namespace) -> int:
    case_id = expand_cases([args.case])[0]
    case_dir = args.output / "cases" / case_id
    lineage_path = case_dir / "lineage.json"
    if not lineage_path.is_file():
        raise ReportError(f"assemble the case before analysis: {case_id}")
    lineage = load_json(lineage_path)
    histories = lineage.get("histories")
    if not isinstance(histories, dict):
        raise ReportError(f"assembled case has no history inventory: {case_id}")
    mhd_record = histories.get("mhd")
    user_record = histories.get("user")
    if not isinstance(mhd_record, dict) or not isinstance(user_record, dict):
        raise ReportError(f"assembled case history inventory is malformed: {case_id}")
    if not mhd_record.get("available") or not user_record.get("available"):
        diagnostics = {
            "schema_version": 1,
            "analyzed_utc": utc_now(),
            "case_id": case_id,
            "case_name": lineage.get("case_name"),
            "assembly_status": lineage.get("status"),
            "analysis_status": "unavailable",
            "analysis_warnings": [],
            "analysis_errors": ["merged MHD and user histories are unavailable"],
            "health": {
                "result": "structural_error",
                "structural_errors": ["merged MHD and user histories are unavailable"],
                "structural_warnings": lineage.get("warnings", []),
                "numerical_warnings": [],
                "science_warnings": [],
            },
        }
        write_json(case_dir / "diagnostics.json", diagnostics)
        write_json(case_dir / "health.json", diagnostics["health"])
        write_text(case_dir / "report.md", case_report_markdown(diagnostics))
        print(f"{case_id}: analysis unavailable")
        return 0
    mhd_path = Path(str(mhd_record["path"]))
    user_path = Path(str(user_record["path"]))
    mhd = parsed_history_columns(mhd_path)
    user = parsed_history_columns(user_path)
    model = lineage.get("model_choices")
    if not isinstance(model, dict):
        model = {}
    health = compute_health(lineage, mhd, user, model)
    warnings: list[str] = []
    errors: list[str] = []
    analyzer = load_pure_analyzer()
    try:
        windows = window_summaries(analyzer, user_path, mhd_path, model)
    except Exception as error:
        windows = {}
        errors.append(f"history analysis failed: {type(error).__name__}: {error}")
    snapshots = load_json(case_dir / "snapshots.json")
    selected = [
        item for item in snapshots.get("snapshots", [])
        if isinstance(item, dict)
        and item.get("complete")
        and isinstance(item.get("time"), (int, float))
        and args.snapshot_time_start - TIME_TOLERANCE
        <= float(item["time"])
        <= args.snapshot_time_end + TIME_TOLERANCE
    ]
    snapshot_records: dict[str, object] = {}
    snapshot_ensemble: dict[str, object] = {"snapshot_count": 0}
    snapshot_status = "skipped" if args.skip_snapshots else "unavailable"
    if not args.skip_snapshots and selected:
        paths = [Path(str(item["representative"])) for item in selected]
        expected = {
            str(path): int(item["expected_ranks"])
            for path, item in zip(paths, selected)
        }
        eddy_samples = (
            args.eddy_samples
            if args.eddy_samples is not None
            else 2_000_000 if case_id in ("R02", "R04") else 0
        )
        try:
            snapshot_records, snapshot_ensemble = analyze_snapshot_paths_bounded(
                analyzer,
                paths,
                args.pdf_bins,
                [int(value) for value in args.alignment_shells.split(",") if value],
                args.snapshot_time_start,
                args.snapshot_time_end,
                model,
                eddy_samples,
                args.eddy_bins,
                args.eddy_seed,
                expected,
                args.snapshot_workers,
                args.snapshot_memory_budget_gib,
            )
            snapshot_status = "complete"
        except Exception as error:
            snapshot_status = "failed"
            warnings.append(
                f"snapshot analysis failed but history analysis was retained: "
                f"{type(error).__name__}: {error}"
            )
    elif not args.skip_snapshots:
        snapshot_status = "not_yet_available"
        warnings.append(
            f"no complete snapshots are available in "
            f"t={args.snapshot_time_start:g}..{args.snapshot_time_end:g}"
        )
    full_history = None
    steady_lf = None
    steady_budget = None
    if isinstance(windows, dict) and isinstance(windows.get("full"), dict):
        full_history = windows["full"].get("history")
    if isinstance(windows, dict) and isinstance(windows.get("steady"), dict):
        steady_lf = windows["steady"].get("lf_history")
        steady_budget = windows["steady"].get("forcing_energy_budget")
    passive = str(model.get("passive_delta", "false")).lower() in ("true", "1")
    if passive and isinstance(steady_lf, dict):
        work = steady_lf.get("applied_pressure_work")
        if isinstance(work, dict) and (
            float(work.get("total", 0.0)) != 0.0
            or float(work.get("anisotropic", 0.0)) != 0.0
        ):
            health["science_warnings"].append(
                f"passive case has nonzero applied pressure work: {work}"
            )
    if not passive and isinstance(steady_budget, dict) and steady_budget.get("available"):
        residual = float(steady_budget.get("relative_residual", math.inf))
        if residual > 1.0e-8:
            health["numerical_warnings"].append(
                f"forcing-energy closure relative residual is {residual:.6g}"
            )
    health["result"] = (
        "structural_error" if health["structural_errors"]
        else "warnings" if (
            health["structural_warnings"]
            or health["numerical_warnings"]
            or health["science_warnings"]
        )
        else "clean"
    )
    diagnostics = {
        "schema_version": 1,
        "analyzed_utc": utc_now(),
        "case_id": case_id,
        "case_name": lineage.get("case_name"),
        "assembly_status": lineage.get("status"),
        "analysis_status": (
            "complete" if not errors and lineage.get("status") == "complete"
            else "partial" if not errors
            else "history_analysis_error"
        ),
        "model_choices": model,
        "health": health,
        "windows": windows,
        "selected_snapshot_count": len(selected),
        "snapshot_analysis_status": snapshot_status,
        "snapshots": snapshot_records,
        "snapshot_ensemble": snapshot_ensemble,
        "analysis_warnings": warnings,
        "analysis_errors": errors,
        "provenance": {
            "lineage": artifact_binding(lineage_path),
            "snapshot_index": artifact_binding(case_dir / "snapshots.json"),
            "merged_mhd_history": artifact_binding(mhd_path),
            "merged_user_history": artifact_binding(user_path),
            "analyzer": artifact_binding(REPO_ROOT / "scripts/analyze_cgl_lf_paper.py"),
            "adapter": artifact_binding(Path(__file__)),
        },
        "compat": {
            "analysis_window": {
                "time_start": args.snapshot_time_start,
                "time_end": args.snapshot_time_end,
            },
            "histories": [full_history] if isinstance(full_history, dict) else [],
            "lf_histories": [steady_lf] if isinstance(steady_lf, dict) else [],
            "forcing_energy_budgets": (
                [steady_budget] if isinstance(steady_budget, dict) else []
            ),
            "snapshot_ensemble": snapshot_ensemble,
        },
    }
    write_json(case_dir / "diagnostics.json", diagnostics)
    write_json(case_dir / "health.json", health)
    write_text(case_dir / "report.md", case_report_markdown(diagnostics))
    print(
        f"{case_id}: {diagnostics['analysis_status']}, health={health['result']}, "
        f"snapshots={snapshot_ensemble.get('snapshot_count', 0)}"
    )
    return 0


def nested_value(value: object, path: str) -> object | None:
    current = value
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def steady_values(diagnostics: dict[str, object]) -> dict[str, float | None]:
    window = nested_value(diagnostics, "windows.steady.history.analysis_window")
    if not isinstance(window, dict):
        return {name: None for name in STEADY_SCALARS}
    return {
        name: float(window[name]) if isinstance(window.get(name), (int, float)) else None
        for name in STEADY_SCALARS
    }


def pair_comparison(
    cases: dict[str, dict[str, object]], left: str, right: str
) -> dict[str, object]:
    if left not in cases or right not in cases:
        return {"available": False, "left": left, "right": right, "reason": "case missing"}
    first = steady_values(cases[left])
    second = steady_values(cases[right])
    scalars: dict[str, object] = {}
    for name in STEADY_SCALARS:
        left_value, right_value = first[name], second[name]
        if left_value is None or right_value is None:
            continue
        scale = max(abs(left_value), abs(right_value), 1.0e-300)
        scalars[name] = {
            "left": left_value,
            "right": right_value,
            "difference_right_minus_left": right_value - left_value,
            "absolute_relative_difference": abs(right_value - left_value) / scale,
            "ratio_right_over_left": (
                right_value / left_value if left_value != 0.0 else None
            ),
        }
    return {
        "available": bool(scalars),
        "left": left,
        "right": right,
        "steady_scalar_comparisons": scalars,
    }


def stationarity_record(diagnostics: dict[str, object]) -> dict[str, object]:
    early = nested_value(diagnostics, "windows.early.history.analysis_window")
    late = nested_value(diagnostics, "windows.late.history.analysis_window")
    if not isinstance(early, dict) or not isinstance(late, dict):
        return {"available": False}
    metrics: dict[str, object] = {}
    for name in STEADY_SCALARS:
        if not isinstance(early.get(name), (int, float)) or not isinstance(
            late.get(name), (int, float)
        ):
            continue
        first, second = float(early[name]), float(late[name])
        metrics[name] = {
            "early": first,
            "late": second,
            "difference": second - first,
            "absolute_relative_change": abs(second - first)
            / max(abs(first), abs(second), 1.0e-300),
        }
    return {"available": bool(metrics), "metrics": metrics}


def resolution_comparison(
    analyzer: object, cases: dict[str, dict[str, object]]
) -> dict[str, object]:
    required = ("R16", "R02", "R17")
    if any(case not in cases for case in required):
        return {"available": False, "reason": "R16/R02/R17 diagnostics are incomplete"}
    try:
        np = analyzer.np
        products: dict[str, object] = {}
        ensembles = {
            case: cases[case]["compat"]["snapshot_ensemble"] for case in required
        }
        for product in ("velocity", "magnetic_fluctuation"):
            curves: dict[str, tuple[object, object]] = {}
            for case in required:
                record = ensembles[case]["spectra"][product]
                k = np.asarray(record["k"], dtype=float) / math.pi
                power = np.asarray(record["power_per_dk"], dtype=float)
                selected = (k >= 4.0) & (k <= 24.0) & (power > 0.0)
                curves[case] = (k[selected], power[selected])
            common = curves["R02"][0]
            selected = (common >= max(curves[case][0][0] for case in required)) & (
                common <= min(curves[case][0][-1] for case in required)
            )
            common = common[selected]
            normalized: dict[str, object] = {}
            for case in required:
                k, power = curves[case]
                values = np.exp(np.interp(np.log(common), np.log(k), np.log(power)))
                values /= np.trapz(values, common)
                normalized[case] = values
            def distance(left: str, right: str) -> float:
                return float(np.sqrt(np.mean(
                    (np.log(normalized[left]) - np.log(normalized[right])) ** 2
                )))
            low = distance("R16", "R02")
            high = distance("R02", "R17")
            products[product] = {
                "k_perp_over_pi": common.tolist(),
                "log_rms_R16_R02": low,
                "log_rms_R02_R17": high,
                "high_to_low_distance_ratio": high / low if low > 0.0 else None,
            }
        return {"available": True, "products": products}
    except Exception as error:
        return {
            "available": False,
            "reason": f"{type(error).__name__}: {error}",
        }


def reference_manifest_paths(matrix_data: dict[str, object]) -> list[Path]:
    panel = matrix_data.get("panel_status")
    values = panel.get("reference_manifests") if isinstance(panel, dict) else None
    if not isinstance(values, dict):
        return []
    root = REPO_ROOT / REFERENCE_ROOT_RELATIVE
    result = []
    for record in values.values():
        if isinstance(record, dict) and isinstance(record.get("path"), str):
            result.append(root / str(record["path"]))
    return result


def table_write(
    directory: Path, name: str, columns: list[str], rows: list[dict[str, object]]
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    csv_path = directory / f"{name}.csv"
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=directory,
        prefix=f".{name}.", delete=False
    ) as stream:
        staged = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(staged, csv_path)
    def tex(value: object) -> str:
        return str(value).replace("\\", r"\textbackslash{}").replace("_", r"\_").replace(
            "%", r"\%"
        ).replace("&", r"\&")
    lines = [
        r"\begin{tabular}{" + "l" * len(columns) + "}",
        r"\hline",
        " & ".join(tex(value) for value in columns) + r" \\",
        r"\hline",
    ]
    lines.extend(
        " & ".join(tex(row.get(column, "")) for column in columns) + r" \\"
        for row in rows
    )
    lines.extend([r"\hline", r"\end{tabular}", ""])
    write_text(directory / f"{name}.tex", "\n".join(lines))


def flatten_reference_rows(reference: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for collection in ("comparisons", "surface_comparisons"):
        values = reference.get(collection)
        if not isinstance(values, dict):
            continue
        for product_id, record in values.items():
            if not isinstance(record, dict):
                continue
            rows.append({
                "product_id": product_id,
                "case": record.get("case"),
                "product": record.get("product"),
                "available": record.get("available"),
                "rms_residual": record.get("rms_residual"),
                "maximum_absolute_residual": record.get("maximum_absolute_residual"),
                "rms_normalized_by_reported_uncertainty": record.get(
                    "rms_normalized_by_reported_uncertainty"
                ),
            })
    for record in reference.get("omitted_products", []):
        if isinstance(record, dict):
            rows.append({
                "product_id": record.get("id"),
                "case": record.get("case"),
                "product": record.get("product"),
                "available": False,
                "reason": record.get("reason"),
            })
    return rows


def write_campaign_tables(
    output: Path, cases: dict[str, dict[str, object]],
    comparisons: dict[str, object], reference: dict[str, object],
) -> None:
    directory = output / "tables"
    matrix_rows = []
    health_rows = []
    steady_rows = []
    stationarity_rows = []
    for case_id, diagnostics in sorted(cases.items()):
        model = diagnostics.get("model_choices")
        if not isinstance(model, dict):
            model = {}
        health = diagnostics.get("health")
        if not isinstance(health, dict):
            health = {}
        matrix_rows.append({
            "case_id": case_id,
            "case_name": diagnostics.get("case_name"),
            "status": diagnostics.get("assembly_status"),
            "resolution": "x".join(
                str(model.get(name, "?")) for name in ("mesh_nx1", "mesh_nx2", "mesh_nx3")
            ),
            "beta0": model.get("beta0"),
            "passive_delta": model.get("passive_delta"),
            "forcing_mode": model.get("forcing_mode"),
            "lf_k_parallel": model.get("lf_k_parallel"),
            "limiter_hardwall": model.get("limiter_hardwall"),
            "limiter_nu_coll": model.get("limiter_nu_coll"),
        })
        health_rows.append({
            "case_id": case_id,
            "status": diagnostics.get("assembly_status"),
            "final_time": health.get("final_time"),
            "health_result": health.get("result"),
            "strict_failure_maximum": max(
                health.get("strict_lf_failure_maxima", {}).values(), default=None
            ) if isinstance(health.get("strict_lf_failure_maxima"), dict) else None,
            "fatal_failure_maximum": max(
                health.get("fatal_lf_failure_maxima", {}).values(), default=None
            ) if isinstance(health.get("fatal_lf_failure_maxima"), dict) else None,
            "hard_bound_diagnostic_maximum": health.get(
                "hard_bound_diagnostic_maximum"
            ),
            "nonfatal_hard_bound_variant": health.get(
                "nonfatal_hard_bound_variant"
            ),
            "mhd_mass_drift": nested_value(health, "mass_relative_drift.mhd"),
            "mass_mismatch": health.get("mhd_user_mass_relative_mismatch"),
            "structural_errors": len(health.get("structural_errors", [])),
            "numerical_warnings": len(health.get("numerical_warnings", [])),
            "science_warnings": len(health.get("science_warnings", [])),
        })
        steady_rows.append({"case_id": case_id, **steady_values(diagnostics)})
        stationarity = stationarity_record(diagnostics)
        if stationarity.get("available"):
            for metric, record in stationarity["metrics"].items():
                stationarity_rows.append({"case_id": case_id, "metric": metric, **record})
    table_write(directory, "case_matrix", list(matrix_rows[0]) if matrix_rows else [], matrix_rows)
    table_write(
        directory, "numerical_health", list(health_rows[0]) if health_rows else [], health_rows
    )
    table_write(
        directory, "steady_state_summary", list(steady_rows[0]) if steady_rows else [], steady_rows
    )
    table_write(
        directory, "stationarity",
        list(stationarity_rows[0]) if stationarity_rows else ["case_id", "metric"],
        stationarity_rows,
    )
    active_rows = []
    for key, record in comparisons.get("active_passive", {}).items():
        if isinstance(record, dict):
            for metric, values in record.get("steady_scalar_comparisons", {}).items():
                active_rows.append({"pair": key, "metric": metric, **values})
    table_write(
        directory, "active_passive",
        list(active_rows[0]) if active_rows else ["pair", "metric"],
        active_rows,
    )
    for name, case_ids in (
        ("heat_flux_scan", ["R12", "R02", "R06", "R13"]),
        ("limiter_scan", ["R14", "R15", "R03", "R07"]),
        ("resolution_scan", ["R16", "R02", "R17"]),
    ):
        rows = [
            {"case_id": case_id, **steady_values(cases[case_id])}
            for case_id in case_ids if case_id in cases
        ]
        table_write(directory, name, list(rows[0]) if rows else ["case_id"], rows)
    reference_rows = flatten_reference_rows(reference)
    table_write(
        directory, "reference_comparisons",
        list(reference_rows[0]) if reference_rows else ["product_id", "case", "available"],
        reference_rows,
    )


def campaign_report_markdown(
    cases: dict[str, dict[str, object]], comparisons: dict[str, object],
    reference: dict[str, object],
) -> str:
    complete = [case for case, value in cases.items() if value.get("assembly_status") == "complete"]
    partial = [case for case, value in cases.items() if value.get("assembly_status") != "complete"]
    lines = [
        "# CGL-LF Stage I Direct-Fast Campaign Report",
        "",
        f"- Cases analyzed: `{len(cases)}`",
        f"- Complete through t=10: `{', '.join(complete) or 'none'}`",
        f"- Partial or unavailable: `{', '.join(partial) or 'none'}`",
        "",
        "## Case Status",
        "",
        "| Case | Assembly | Analysis | Health | Final time |",
        "| --- | --- | --- | --- | ---: |",
    ]
    for case_id, value in sorted(cases.items()):
        health = value.get("health", {})
        lines.append(
            f"| {case_id} | {value.get('assembly_status')} | "
            f"{value.get('analysis_status')} | "
            f"{health.get('result') if isinstance(health, dict) else None} | "
            f"{health.get('final_time') if isinstance(health, dict) else None} |"
        )
    lines.extend([
        "",
        "## Comparison Boundary",
        "",
        "The comparison products are descriptive evidence. Unexpected or unresolved "
        "behavior remains visible as warnings and is not silently converted into "
        "a failed or passing scientific claim.",
        "",
        f"- Reference products available: `{len(flatten_reference_rows(reference))}`",
        f"- Comparison families: `{', '.join(sorted(comparisons))}`",
        "",
    ])
    return "\n".join(lines)


def command_aggregate(args: argparse.Namespace) -> int:
    selected = expand_cases(args.cases)
    case_diagnostics: dict[str, dict[str, object]] = {}
    compat_cases: dict[str, object] = {}
    for case_id in selected:
        path = args.output / "cases" / case_id / "diagnostics.json"
        if not path.is_file():
            continue
        diagnostics = load_json(path)
        case_diagnostics[case_id] = diagnostics
        compat = diagnostics.get("compat")
        ensemble = compat.get("snapshot_ensemble") if isinstance(compat, dict) else None
        if (
            isinstance(compat, dict)
            and isinstance(ensemble, dict)
            and int(ensemble.get("snapshot_count", 0)) > 0
        ):
            compat_cases[str(diagnostics["case_name"])] = compat
    comparisons: dict[str, object] = {}
    for group, groups in COMPARISON_GROUPS.items():
        records: dict[str, object] = {}
        for values in groups:
            if len(values) == 2:
                records[f"{values[0]}_{values[1]}"] = pair_comparison(
                    case_diagnostics, values[0], values[1]
                )
            else:
                records["_".join(values)] = {
                    "available_cases": [
                        case for case in values if case in case_diagnostics
                    ],
                    "steady_values": {
                        case: steady_values(case_diagnostics[case])
                        for case in values if case in case_diagnostics
                    },
                }
        comparisons[group] = records
    comparisons["stationarity"] = {
        case_id: stationarity_record(diagnostics)
        for case_id, diagnostics in case_diagnostics.items()
    }
    analyzer = load_pure_analyzer()
    comparisons["resolution_detail"] = resolution_comparison(analyzer, case_diagnostics)
    _, matrix_data = matrix_cases(args.matrix)
    reference: dict[str, object] = {
        "available": False,
        "reason": "no snapshot-analyzed cases are available",
        "comparisons": {},
        "surface_comparisons": {},
        "omitted_products": [],
    }
    if compat_cases:
        try:
            configuration = analyzer.stage_i_panels_configuration(args.matrix)
            reference = analyzer.combined_reference_curve_comparisons(
                {"cases": compat_cases},
                reference_manifest_paths(matrix_data),
                allow_missing_cases=True,
                analysis_case_aliases=configuration["analysis_case_aliases"],
                stage_i_reference_bindings=configuration["reference_product_bindings"],
            )
        except Exception as error:
            reference = {
                "available": False,
                "reason": f"{type(error).__name__}: {error}",
                "comparisons": {},
                "surface_comparisons": {},
                "omitted_products": [],
            }
    campaign = args.output / "campaign"
    campaign.mkdir(parents=True, exist_ok=True)
    diagnostics = {
        "schema_version": 1,
        "aggregated_utc": utc_now(),
        "cases": compat_cases,
        "case_status": {
            case_id: {
                "case_name": value.get("case_name"),
                "assembly_status": value.get("assembly_status"),
                "analysis_status": value.get("analysis_status"),
                "health": value.get("health"),
            }
            for case_id, value in case_diagnostics.items()
        },
        "reference_curve_comparisons": reference,
    }
    write_json(campaign / "diagnostics.json", diagnostics)
    write_json(campaign / "comparisons.json", comparisons)
    write_json(campaign / "reference_comparisons.json", reference)
    health = {
        "complete_cases": [
            case for case, value in case_diagnostics.items()
            if value.get("assembly_status") == "complete"
        ],
        "partial_or_unavailable_cases": [
            case for case, value in case_diagnostics.items()
            if value.get("assembly_status") != "complete"
        ],
        "structural_error_cases": [
            case for case, value in case_diagnostics.items()
            if nested_value(value, "health.result") == "structural_error"
        ],
        "warning_cases": [
            case for case, value in case_diagnostics.items()
            if nested_value(value, "health.result") == "warnings"
        ],
    }
    write_json(campaign / "health.json", health)
    write_text(campaign / "report.md", campaign_report_markdown(
        case_diagnostics, comparisons, reference
    ))
    write_campaign_tables(args.output, case_diagnostics, comparisons, reference)
    print(
        f"aggregated {len(case_diagnostics)} cases; "
        f"complete={len(health['complete_cases'])}, "
        f"partial={len(health['partial_or_unavailable_cases'])}"
    )
    return 0


def load_plotter() -> object:
    path = REPO_ROOT / "scripts/plot_cgl_lf_paper.py"
    spec = importlib.util.spec_from_file_location("_cgl_lf_fast_report_plotter", path)
    if spec is None or spec.loader is None:
        raise ReportError(f"cannot load plotter: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def render_history_dashboard(
    output: Path, case_ids: list[str], generated: list[str]
) -> list[str]:
    warnings: list[str] = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:
        return [f"history dashboard unavailable: {type(error).__name__}: {error}"]
    fields = (
        ("kinetic", "Kinetic energy / volume"),
        ("magnetic", "Magnetic energy / volume"),
        ("beta", "Mean beta"),
        ("abs_dp", "Mean |Delta p|"),
        ("nu_eff", "Mean effective collision rate"),
        ("hard_vol", "Hard-bound volume fraction"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(11.0, 9.0), sharex=True)
    for case_id in case_ids:
        lineage_path = output / "cases" / case_id / "lineage.json"
        if not lineage_path.is_file():
            continue
        lineage = load_json(lineage_path)
        user = nested_value(lineage, "histories.user.path")
        if not isinstance(user, str) or not Path(user).is_file():
            continue
        try:
            data = parsed_history_columns(Path(user))
        except Exception as error:
            warnings.append(f"{case_id} dashboard history failed: {error}")
            continue
        volume = data.get("volume", [1.0] * len(data["time"]))
        for axis, (name, label) in zip(axes.flat, fields):
            if name not in data:
                continue
            values = [
                value / vol if vol != 0.0 else math.nan
                for value, vol in zip(data[name], volume)
            ]
            axis.plot(data["time"], values, lw=0.8, label=case_id)
            axis.set_ylabel(label)
            axis.grid(True, alpha=0.25)
    for axis in axes[-1]:
        axis.set_xlabel(r"$t/(L_\perp/v_A)$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="outside lower center", ncol=8, fontsize=7)
    fig.suptitle("CGL-LF Stage I global histories")
    path = output / "figures/history_dashboard.pdf"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    generated.append(str(path.relative_to(output)))
    return warnings


def command_render(args: argparse.Namespace) -> int:
    diagnostics_path = args.output / "campaign/diagnostics.json"
    if not diagnostics_path.is_file():
        raise ReportError("aggregate the campaign before rendering")
    data = load_json(diagnostics_path)
    figure_dir = args.output / "figures/paper"
    figure_dir.mkdir(parents=True, exist_ok=True)
    plotter = load_plotter()
    generated: list[str] = []
    warnings: list[str] = []
    ensembles = plotter.case_ensembles(data)
    operations = (
        ("pdfs", lambda: plotter.plot_pdfs(ensembles, figure_dir, generated)),
        (
            "pressure density",
            lambda: plotter.plot_pressure_density_joint(ensembles, figure_dir, generated),
        ),
        (
            "compressive spectra",
            lambda: plotter.plot_compressive_pressure_spectra(
                ensembles, figure_dir, generated
            ),
        ),
        (
            "projected spectra",
            lambda: plotter.plot_spatial_spectra(ensembles, figure_dir, generated),
        ),
        (
            "transfer alignment",
            lambda: plotter.plot_transfer_and_alignment(ensembles, figure_dir, generated),
        ),
        (
            "heat flux",
            lambda: plotter.plot_heat_flux_proxy(ensembles, figure_dir, generated),
        ),
        (
            "pressure work",
            lambda: plotter.plot_pressure_work(ensembles, figure_dir, generated),
        ),
        (
            "eddy anisotropy",
            lambda: plotter.plot_eddy_anisotropy(ensembles, figure_dir, generated),
        ),
        (
            "reference curves",
            lambda: plotter.plot_reference_comparisons(data, figure_dir, generated),
        ),
        (
            "reference surfaces",
            lambda: plotter.plot_reference_surface_comparisons(
                data, figure_dir, generated
            ),
        ),
    )
    for label, operation in operations:
        try:
            operation()
        except Exception as error:
            warnings.append(f"{label} rendering failed: {type(error).__name__}: {error}")
    case_ids = sorted(load_json(args.output / "inventory.json").get("cases", {}))
    generated = [
        str(Path("figures/paper") / name) if "/" not in name else name
        for name in generated
    ]
    warnings.extend(render_history_dashboard(args.output, case_ids, generated))
    index = {
        "generated_utc": utc_now(),
        "figures": sorted(set(generated)),
        "warnings": warnings,
        "plotter": artifact_binding(REPO_ROOT / "scripts/plot_cgl_lf_paper.py"),
    }
    write_json(args.output / "figures/figure_manifest.json", index)
    print(f"rendered {len(index['figures'])} figures with {len(warnings)} warning(s)")
    return 0


def latex_macro_name(case_id: str, metric: str) -> str:
    pieces = re.findall(r"[A-Za-z0-9]+", metric)
    suffix = "".join(piece[:1].upper() + piece[1:] for piece in pieces)
    return f"CGL{case_id}{suffix}"


def command_manuscript(args: argparse.Namespace) -> int:
    campaign = args.output / "campaign"
    diagnostics = load_json(campaign / "diagnostics.json")
    comparisons = load_json(campaign / "comparisons.json")
    health = load_json(campaign / "health.json")
    manuscript = args.output / "manuscript"
    manuscript.mkdir(parents=True, exist_ok=True)
    inventory = load_json(args.output / "inventory.json")
    case_results: dict[str, object] = {}
    macros: list[str] = [
        "% Generated by cgl_lf_stage_i_fast_report.py; do not edit by hand."
    ]
    warnings: list[str] = []
    for case_id in sorted(inventory.get("cases", {})):
        path = args.output / "cases" / case_id / "diagnostics.json"
        if not path.is_file():
            continue
        case = load_json(path)
        values = steady_values(case)
        case_results[case_id] = {
            "case_name": case.get("case_name"),
            "assembly_status": case.get("assembly_status"),
            "analysis_status": case.get("analysis_status"),
            "health_result": nested_value(case, "health.result"),
            "steady_state": values,
            "stationarity": stationarity_record(case),
        }
        for metric, value in values.items():
            if value is not None:
                macros.append(
                    rf"\newcommand{{\{latex_macro_name(case_id, metric)}}}"
                    rf"{{{value:.8g}}}"
                )
        case_health = case.get("health")
        if isinstance(case_health, dict):
            warnings.extend(
                f"{case_id}: {warning}"
                for key in ("structural_warnings", "numerical_warnings", "science_warnings")
                for warning in case_health.get(key, [])
            )
    results = {
        "schema_version": 1,
        "generated_utc": utc_now(),
        "case_results": case_results,
        "comparisons": comparisons,
        "campaign_health": health,
        "reference_comparisons": diagnostics.get("reference_curve_comparisons"),
        "interpretation_boundary": (
            "Quantitative products are evidence inputs. Claims require author review "
            "and must distinguish direct results from interpretation."
        ),
    }
    write_json(manuscript / "results.json", results)
    write_text(manuscript / "results_macros.tex", "\n".join(macros) + "\n")
    all_complete = not health.get("partial_or_unavailable_cases")
    no_structural = not health.get("structural_error_cases")
    claims = [
        {
            "claim_id": "campaign_complete",
            "claim": "All Stage I cases R02-R17 reach t=10.",
            "status": "supported" if all_complete else "pending",
            "evidence": ["campaign/health.json", "tables/numerical_health.csv"],
            "scope": "direct-fast Stage I campaign",
        },
        {
            "claim_id": "numerical_integrity",
            "claim": "The analyzed campaign has no structural data-integrity errors.",
            "status": "supported" if no_structural else "needs_review",
            "evidence": ["campaign/health.json", "tables/numerical_health.csv"],
            "scope": "retained histories and indexed snapshots",
        },
        {
            "claim_id": "active_passive_feedback",
            "claim": "Active/passive comparisons measure the effect of pressure feedback.",
            "status": "evidence_available",
            "evidence": ["campaign/comparisons.json", "tables/active_passive.csv"],
            "scope": "matched R02-R09 pairs; interpretation requires review",
        },
        {
            "claim_id": "resolution",
            "claim": "R16/R02/R17 provide a resolution-convergence test.",
            "status": (
                "evidence_available"
                if nested_value(comparisons, "resolution_detail.available")
                else "pending"
            ),
            "evidence": ["campaign/comparisons.json", "tables/resolution_scan.csv"],
            "scope": "common retained scales",
        },
    ]
    write_json(manuscript / "claim_evidence.json", {"claims": claims})
    claim_lines = [
        "# Claim-Evidence Ledger",
        "",
        "| Claim | Status | Scope | Evidence |",
        "| --- | --- | --- | --- |",
    ]
    claim_lines.extend(
        f"| {claim['claim']} | {claim['status']} | {claim['scope']} | "
        f"{', '.join(claim['evidence'])} |"
        for claim in claims
    )
    write_text(manuscript / "claim_evidence.md", "\n".join(claim_lines) + "\n")
    figure_manifest = args.output / "figures/figure_manifest.json"
    write_json(
        manuscript / "figure_manifest.json",
        load_json(figure_manifest) if figure_manifest.is_file()
        else {"figures": [], "warnings": ["render has not been run"]},
    )
    provenance_inputs = [
        args.output / "manifest.json",
        args.output / "inventory.json",
        campaign / "diagnostics.json",
        campaign / "comparisons.json",
        campaign / "health.json",
        Path(__file__),
        REPO_ROOT / "scripts/analyze_cgl_lf_paper.py",
        REPO_ROOT / "scripts/plot_cgl_lf_paper.py",
        args.matrix,
    ]
    if WRITING_GUIDE.is_file():
        provenance_inputs.append(WRITING_GUIDE)
    write_json(manuscript / "provenance_manifest.json", {
        "schema_version": 1,
        "generated_utc": utc_now(),
        "inputs": [artifact_binding(path) for path in provenance_inputs if path.is_file()],
        "case_inputs": {
            case_id: record.get("input")
            for case_id, record in inventory.get("cases", {}).items()
            if isinstance(record, dict)
        },
    })
    open_lines = [
        "# Open Questions And Required Review",
        "",
        "- Select the target journal and finalize its formatting requirements.",
        "- Review every evidence-available claim before converting it to manuscript prose.",
        "- Keep CGL-LF conclusions distinct from kinetic validation.",
        "- Resolve any remaining MKS24 normalization-blocked comparisons.",
        "",
        "## Retained Warnings",
        "",
    ]
    open_lines.extend(f"- {warning}" for warning in warnings)
    if not warnings:
        open_lines.append("- None.")
    write_text(manuscript / "open_questions.md", "\n".join(open_lines) + "\n")
    print(f"wrote manuscript-ready outputs for {len(case_results)} cases")
    return 0


def verify_snapshot_index(index: dict[str, object]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    times: list[float] = []
    for item in index.get("snapshots", []):
        if not isinstance(item, dict):
            errors.append("snapshot index contains a non-object entry")
            continue
        time = item.get("time")
        if isinstance(time, (int, float)):
            times.append(float(time))
        members = item.get("rank_files")
        if not isinstance(members, list) or len(members) != item.get("expected_ranks"):
            errors.append(f"snapshot rank inventory differs: {item.get('representative')}")
            continue
        missing = [
            str(member.get("path"))
            for member in members if isinstance(member, dict)
            and not Path(str(member.get("path"))).is_file()
        ]
        if missing:
            errors.append(
                f"snapshot now lacks {len(missing)} rank file(s): {item.get('representative')}"
            )
    if any(right <= left for left, right in zip(times, times[1:])):
        errors.append("snapshot index physical times are not strictly increasing")
    if int(index.get("complete_snapshot_count", 0)) < int(index.get("snapshot_count", 0)):
        warnings.append("snapshot index contains incomplete groups")
    return errors, warnings


def command_verify(args: argparse.Namespace) -> int:
    errors: list[str] = []
    warnings: list[str] = []
    inventory_path = args.output / "inventory.json"
    if not inventory_path.is_file():
        raise ReportError("assemble before verification")
    inventory = load_json(inventory_path)
    selected = expand_cases(args.cases)
    case_results: dict[str, object] = {}
    for case_id in selected:
        case_errors: list[str] = []
        case_warnings: list[str] = []
        lineage_path = args.output / "cases" / case_id / "lineage.json"
        snapshot_path = args.output / "cases" / case_id / "snapshots.json"
        if not lineage_path.is_file():
            case_errors.append("lineage.json is missing")
            case_results[case_id] = {"errors": case_errors, "warnings": case_warnings}
            errors.extend(f"{case_id}: {value}" for value in case_errors)
            continue
        lineage = load_json(lineage_path)
        retained_lineage_errors = lineage.get("errors")
        if isinstance(retained_lineage_errors, list):
            case_errors.extend(
                f"lineage error: {value}" for value in retained_lineage_errors
            )
        elif retained_lineage_errors is not None:
            case_errors.append("lineage errors field is malformed")
        retained_segments = lineage.get("lineage")
        if not isinstance(retained_segments, list) or not retained_segments:
            case_errors.append("selected lineage is missing or malformed")
        else:
            for index, segment in enumerate(retained_segments):
                if not isinstance(segment, dict):
                    case_errors.append(
                        f"selected lineage segment {index} is malformed"
                    )
                    continue
                exit_code = segment.get("run_exit_code")
                if isinstance(exit_code, int) and not isinstance(exit_code, bool):
                    if exit_code != 0:
                        case_errors.append(
                            f"selected lineage segment {index} has nonzero "
                            f"run exit code: {exit_code}"
                        )
                elif exit_code is not None:
                    case_errors.append(
                        f"selected lineage segment {index} run exit code is malformed"
                    )
                if segment.get("kind") in {"fast", "historical_seed_prefix"}:
                    for identity in (
                        "input_sha256", "matrix_sha256", "executable_sha256"
                    ):
                        value = segment.get(identity)
                        if not isinstance(value, str) or re.fullmatch(
                            r"[0-9a-f]{64}", value
                        ) is None:
                            case_errors.append(
                                f"selected lineage segment {index} lacks valid "
                                f"{identity}"
                            )
        if args.require_complete and lineage.get("status") != "complete":
            case_errors.append(f"case is not complete: {lineage.get('status')}")
        elif lineage.get("status") != "complete":
            case_warnings.append(f"case is not complete: {lineage.get('status')}")
        for kind in ("mhd", "user"):
            record = nested_value(lineage, f"histories.{kind}")
            if not isinstance(record, dict) or not record.get("available"):
                case_errors.append(f"merged {kind} history is unavailable")
                continue
            path = Path(str(record["path"]))
            if not path.is_file():
                case_errors.append(f"merged {kind} history is missing: {path}")
                continue
            binding = record.get("binding")
            if isinstance(binding, dict) and binding.get("sha256") != sha256_file(path):
                case_errors.append(f"merged {kind} history digest changed: {path}")
            try:
                data = parsed_history_columns(path)
                if any(
                    right <= left for left, right in zip(data["time"], data["time"][1:])
                ):
                    case_errors.append(f"merged {kind} history time is not increasing")
            except Exception as error:
                case_errors.append(f"merged {kind} history parse failed: {error}")
        if snapshot_path.is_file():
            snapshot_errors, snapshot_warnings = verify_snapshot_index(load_json(snapshot_path))
            case_errors.extend(snapshot_errors)
            case_warnings.extend(snapshot_warnings)
        else:
            case_errors.append("snapshots.json is missing")
        diagnostics_path = args.output / "cases" / case_id / "diagnostics.json"
        if diagnostics_path.is_file():
            diagnostics = load_json(diagnostics_path)
            for field in ("errors", "analysis_errors"):
                diagnostic_errors = diagnostics.get(field)
                if isinstance(diagnostic_errors, list):
                    case_errors.extend(
                        f"diagnostic error: {value}" for value in diagnostic_errors
                    )
                elif diagnostic_errors is not None:
                    case_errors.append(f"diagnostics {field} field is malformed")
            structural_errors = nested_value(diagnostics, "health.structural_errors")
            if isinstance(structural_errors, list):
                case_errors.extend(
                    f"diagnostic structural error: {value}"
                    for value in structural_errors
                )
            elif structural_errors is not None:
                case_errors.append("diagnostic structural errors field is malformed")
            case_warnings.extend(diagnostics.get("analysis_warnings", []))
            case_warnings.extend(nested_value(diagnostics, "health.numerical_warnings") or [])
            case_warnings.extend(nested_value(diagnostics, "health.science_warnings") or [])
        else:
            case_warnings.append("case diagnostics have not been generated")
        case_results[case_id] = {"errors": case_errors, "warnings": case_warnings}
        errors.extend(f"{case_id}: {value}" for value in case_errors)
        warnings.extend(f"{case_id}: {value}" for value in case_warnings)
    for path in (
        args.output / "campaign/diagnostics.json",
        args.output / "campaign/comparisons.json",
        args.output / "campaign/health.json",
    ):
        if not path.is_file():
            warnings.append(f"campaign output has not been generated: {path.name}")
    result = {
        "schema_version": 1,
        "verified_utc": utc_now(),
        "result": "fail" if errors else "warnings" if warnings else "pass",
        "require_complete": args.require_complete,
        "cases": case_results,
        "errors": errors,
        "warnings": warnings,
        "inventory": artifact_binding(inventory_path),
        "adapter": artifact_binding(Path(__file__)),
    }
    write_json(args.output / "verify.json", result)
    print(f"verify: {result['result']}; errors={len(errors)}, warnings={len(warnings)}")
    return 1 if errors else 0


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    command.add_argument(
        "--output", type=Path, default=DEFAULT_ROOT / DEFAULT_OUTPUT_RELATIVE
    )
    command.add_argument("--frozen-source", type=Path, default=DEFAULT_FROZEN_SOURCE)
    command.add_argument("--matrix", type=Path)
    subcommands = command.add_subparsers(dest="command", required=True)

    assemble = subcommands.add_parser(
        "assemble", help="assemble deterministic read-only case lineages"
    )
    assemble.add_argument("cases", nargs="*", default=["R02-R17"])

    analyze = subcommands.add_parser(
        "analyze-case", help="analyze one assembled complete or partial case"
    )
    analyze.add_argument("case")
    analyze.add_argument("--skip-snapshots", action="store_true")
    analyze.add_argument("--snapshot-time-start", type=float, default=8.0)
    analyze.add_argument("--snapshot-time-end", type=float, default=10.0)
    analyze.add_argument(
        "--snapshot-workers",
        type=int,
        default=DEFAULT_SNAPSHOT_WORKERS,
        help=(
            "maximum independent snapshot-analysis processes; bounded by CPU "
            "affinity and --snapshot-memory-budget-gib"
        ),
    )
    analyze.add_argument(
        "--snapshot-memory-budget-gib",
        type=float,
        default=DEFAULT_SNAPSHOT_MEMORY_BUDGET_GIB,
        help=(
            "conservative aggregate memory budget for coordinator and snapshot "
            "workers"
        ),
    )
    analyze.add_argument("--pdf-bins", type=int, default=64)
    analyze.add_argument(
        "--alignment-shells", default="2,4,6,8,12,16,24,32,64,128"
    )
    analyze.add_argument("--eddy-samples", type=int)
    analyze.add_argument("--eddy-bins", type=int, default=24)
    analyze.add_argument("--eddy-seed", type=int, default=731)

    aggregate = subcommands.add_parser(
        "aggregate", help="aggregate available case analyses and write tables"
    )
    aggregate.add_argument("cases", nargs="*", default=["R02-R17"])

    subcommands.add_parser("render", help="render available campaign figures")
    subcommands.add_parser(
        "manuscript", help="write manuscript-ready data and provenance outputs"
    )

    verify = subcommands.add_parser(
        "verify", help="verify generated outputs and retained source inventory"
    )
    verify.add_argument("cases", nargs="*", default=["R02-R17"])
    verify.add_argument("--require-complete", action="store_true")
    return command


def main(argv: list[str] | None = None) -> int:
    command = parser()
    args = command.parse_args(argv)
    args.root = args.root.expanduser().absolute()
    args.output = args.output.expanduser().absolute()
    args.frozen_source = args.frozen_source.expanduser().absolute()
    args.matrix = (
        args.matrix.expanduser().absolute()
        if args.matrix is not None
        else args.frozen_source / MATRIX_RELATIVE
    )
    try:
        require_output_is_separate(args.root, args.frozen_source, args.output)
        if not args.matrix.is_file():
            raise ReportError(f"matrix does not exist: {args.matrix}")
        if args.command == "assemble":
            return command_assemble(args)
        if args.command == "analyze-case":
            return command_analyze_case(args)
        if args.command == "aggregate":
            return command_aggregate(args)
        if args.command == "render":
            return command_render(args)
        if args.command == "manuscript":
            return command_manuscript(args)
        if args.command == "verify":
            return command_verify(args)
        raise ReportError(f"unsupported command: {args.command}")
    except (ReportError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
