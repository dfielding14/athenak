#!/usr/bin/env python3
"""Render deterministic publication summaries for direct-fast Stage I evidence.

This renderer is intentionally downstream-only.  It reads the direct-fast
report output, discovers optional scientific-acceptance evidence, and writes a
fixed set of summary figures and tables.  Missing or partial campaign evidence
is shown explicitly rather than treated as a failed result.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
import hashlib
import io
import json
import math
import os
import platform
from pathlib import Path
import re
import sys
import tempfile
import textwrap
from typing import Any, Iterable


CASE_IDS = tuple(f"R{number:02d}" for number in range(2, 18))
ACTIVE_PASSIVE_PAIRS = (
    ("R02", "R06"),
    ("R03", "R07"),
    ("R04", "R08"),
    ("R05", "R09"),
)
ROBUSTNESS_CONTRASTS = (
    ("forcing A, beta=10", "R02", "R04"),
    ("forcing A, beta=100", "R03", "R05"),
    ("forcing P, beta=10", "R06", "R08"),
    ("forcing P, beta=100", "R07", "R09"),
    ("beta, A Alfvenic", "R02", "R03"),
    ("beta, A random", "R04", "R05"),
    ("beta, P Alfvenic", "R06", "R07"),
    ("beta, P random", "R08", "R09"),
    ("forcing correlation", "R05", "R11"),
)
HEAT_FLUX_CASES = ("R12", "R02", "R06", "R13")
LIMITER_CASES = ("R14", "R15", "R03", "R07")
RESOLUTION_CASES = ("R16", "R02", "R17")
TARGET_TIME = 10.0
TIME_TOLERANCE = 1.0e-10
HISTORY_LABEL = re.compile(r"\[(\d+)\]=(\S+)")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
CASE_PATH_PATTERN = re.compile(r"(?:^|/)(R(?:0[2-9]|1[0-7]))(?:/|$)")
CASE_DIRECTORY_PATTERN = re.compile(
    r"(?:^|/)cases/(R(?:0[2-9]|1[0-7]))(?:/|$)"
)
EVIDENCE_DIGEST_METHOD = "sha256-canonical-json-without-evidence_digest-v1"
CT_EVIDENCE_DIGEST_METHOD = "sha256-canonical-json-without-evidence-digest"
SCIENCE_RECORD_TYPE = "cgl-lf-stage-i-direct-fast-reviewed-science-comparisons"
SCIENCE_PROVENANCE_RECORD_TYPE = (
    "cgl-lf-stage-i-direct-fast-reviewed-science-provenance"
)
SCIENCE_AUTHORITY = "non-authorizing-direct-fast-scientific-assessment"
CT_AUDIT_RECORD_TYPE = "stage-i-direct-fast-ct-audit"
EVIDENCE_RESULTS = {"pass", "fail", "inconclusive"}
SCIENCE_CONTRAST_RESULTS = {*EVIDENCE_RESULTS, "available"}
R15_STRICT_FAILURE_RECORD_TYPE = "cgl-lf-stage-i-retained-strict-failure-evidence"
RENDERER_PATH = Path(__file__).resolve()
ACCEPTANCE_RECORD_TYPES = {
    "stage-i-scientific-criteria-validation",
    "stage-i-scientific-case-evidence",
    "stage-i-scientific-campaign-evidence",
    "stage-i-scientific-ct-divb-evidence",
    "cgl-lf-stage-i-direct-fast-case-acceptance",
    "cgl-lf-stage-i-direct-fast-campaign-evidence",
}

STATUS_COLORS = {
    "pass": "#3a923a",
    "clean": "#3a923a",
    "warning": "#e6a83a",
    "warnings": "#e6a83a",
    "restricted": "#8c6bb1",
    "configuration": "#b8c4cc",
    "incomplete": "#4c78a8",
    "inconclusive": "#9e9e9e",
    "blocked_out_of_scope": "#9e9e9e",
    "fail": "#c43c39",
    "structural_error": "#c43c39",
    "unknown": "#e6e6e6",
}
CASE_COLORS = {
    case_id: color
    for case_id, color in zip(
        CASE_IDS,
        (
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            "#bcbd22", "#17becf", "#4e79a7", "#f28e2b",
            "#59a14f", "#e15759", "#b07aa1", "#76b7b2",
        ),
    )
}


class PublicationError(RuntimeError):
    """An unrecoverable publication-renderer input or output error."""


@dataclass
class CaseRecord:
    """All available downstream evidence for one Stage I case."""

    case_id: str
    diagnostics: dict[str, Any] | None = None
    lineage: dict[str, Any] | None = None
    acceptance: dict[str, Any] | None = None
    direct_acceptance: dict[str, Any] | None = None
    model: dict[str, Any] = field(default_factory=dict)
    user_history: dict[str, list[float]] | None = None
    mhd_history: dict[str, list[float]] | None = None
    lineage_path: Path | None = None
    history_paths: dict[str, Path] = field(default_factory=dict)


@dataclass
class PublicationData:
    """Normalized fast-report and acceptance evidence."""

    analysis: Path
    cases: dict[str, CaseRecord]
    aggregate: dict[str, Any] | None
    comparisons: dict[str, Any] | None
    campaign_acceptance: dict[str, Any] | None
    science_record: dict[str, Any] | None
    ct_audit_record: dict[str, Any] | None
    acceptance_records: list[dict[str, Any]]
    audit_records: list[dict[str, Any]]
    source_paths: set[Path]
    ingestion_warnings: list[str]


def integrated_audit_records(data: PublicationData) -> list[dict[str, Any]]:
    """Return authenticated audit records represented in publication products."""

    records = list(data.audit_records)
    if isinstance(data.ct_audit_record, dict):
        records.append(data.ct_audit_record)
    return records


def nested(value: object, path: str) -> object | None:
    """Return one dotted-path value from nested dictionaries."""

    current = value
    for part in path.split("."):
        if isinstance(current, dict) and part in current:
            current = current[part]
        elif isinstance(current, (list, tuple)) and part.isdigit():
            index = int(part)
            if index >= len(current):
                return None
            current = current[index]
        else:
            return None
    return current


def as_float(value: object) -> float | None:
    """Return one finite floating-point value, or None."""

    if not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def text_value(value: object) -> str:
    """Format one deterministic table cell."""

    if value is None:
        return "--"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "--"
        return f"{value:.8g}"
    if isinstance(value, (list, tuple)):
        return "; ".join(text_value(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return str(value)


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: object) -> bytes:
    """Return the canonical JSON representation used by scientific evidence."""

    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def evidence_digest(value: dict[str, Any]) -> str:
    """Return the scientific-evidence digest excluding its self-digest member."""

    body = dict(value)
    body.pop("evidence_digest", None)
    return hashlib.sha256(canonical_json(body)).hexdigest()


def binding_path(binding: object) -> Path | None:
    """Return the absolute path declared by one file binding."""

    if not isinstance(binding, dict) or not isinstance(binding.get("path"), str):
        return None
    return Path(binding["path"]).expanduser().absolute()


def verify_file_binding(binding: object, label: str) -> list[str]:
    """Return validation errors for one required path/SHA-256/size binding."""

    if not isinstance(binding, dict):
        return [f"{label} is not a file binding"]
    path = binding_path(binding)
    expected_sha = binding.get("sha256")
    expected_size = binding.get("size_bytes")
    if path is None:
        return [f"{label} lacks a path"]
    if not isinstance(expected_sha, str) or SHA256_PATTERN.fullmatch(expected_sha) is None:
        return [f"{label} lacks a valid SHA-256"]
    if not path.is_file():
        return [f"{label} path is unavailable: {path}"]
    errors: list[str] = []
    try:
        observed_size = path.stat().st_size
        observed_sha = sha256_file(path)
    except OSError as error:
        return [f"{label} could not be bound: {type(error).__name__}: {error}"]
    if isinstance(expected_size, int) and expected_size != observed_size:
        errors.append(
            f"{label} size differs: declared={expected_size} observed={observed_size}"
        )
    elif expected_size is not None and not isinstance(expected_size, int):
        errors.append(f"{label} size_bytes is not an integer")
    if expected_sha != observed_sha:
        errors.append(f"{label} SHA-256 differs")
    return errors


def recursive_file_bindings(value: object) -> list[dict[str, Any]]:
    """Return every nested object that declares a path and SHA-256."""

    bindings: list[dict[str, Any]] = []
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and "sha256" in value:
            bindings.append(value)
        for child in value.values():
            bindings.extend(recursive_file_bindings(child))
    elif isinstance(value, list):
        for child in value:
            bindings.extend(recursive_file_bindings(child))
    return bindings


def verify_record_digest(
    record: dict[str, Any], label: str, expected_method: str
) -> list[str]:
    """Return validation errors for one canonical self-digest."""

    digest = record.get("evidence_digest")
    if not isinstance(digest, dict):
        return [f"{label} lacks a scientific-evidence self-digest"]
    if digest.get("method") != expected_method:
        return [f"{label} uses an unsupported evidence-digest method"]
    expected = digest.get("sha256")
    if not isinstance(expected, str) or SHA256_PATTERN.fullmatch(expected) is None:
        return [f"{label} evidence digest is not a valid SHA-256"]
    try:
        observed = evidence_digest(record)
    except (TypeError, ValueError) as error:
        return [f"{label} cannot be canonically digested: {error}"]
    return [] if observed == expected else [f"{label} evidence self-digest differs"]


def verify_self_digest(record: dict[str, Any], label: str) -> list[str]:
    """Return validation errors for one scientific-evidence self-digest."""

    return verify_record_digest(record, label, EVIDENCE_DIGEST_METHOD)


def binding_freshness_errors(
    binding: object, current_path: Path | None, label: str
) -> list[str]:
    """Return errors when evidence does not bind the selected current case input."""

    if current_path is None:
        return [f"{label} cannot be matched to the selected current case input"]
    declared_path = binding_path(binding)
    if declared_path is None:
        return [f"{label} lacks a bound path"]
    try:
        declared_resolved = declared_path.resolve(strict=True)
        current_resolved = current_path.resolve(strict=True)
    except OSError as error:
        return [f"{label} freshness path is unavailable: {error}"]
    if declared_resolved != current_resolved:
        return [
            f"{label} is stale: evidence path {declared_resolved} differs from "
            f"selected path {current_resolved}"
        ]
    return verify_file_binding(binding, label)


def direct_fast_output_provenance_errors(
    path: Path, record: dict[str, Any], label: str
) -> list[str]:
    """Return errors when a direct-fast record lacks its external output binding."""

    case_id = record.get("case_id")
    if record.get("record_type") == "cgl-lf-stage-i-direct-fast-case-acceptance":
        provenance_path = path.parents[2] / "provenance.json"
        output_path = f"outputs.case_acceptance.{case_id}"
    else:
        provenance_path = path.parent / "provenance.json"
        output_path = "outputs.campaign_evidence"
    if not provenance_path.is_file():
        return [f"{label} lacks direct-fast external provenance: {provenance_path}"]
    try:
        provenance = load_json(provenance_path)
    except (OSError, json.JSONDecodeError, PublicationError) as error:
        return [f"{label} external provenance is unreadable: {error}"]
    binding = nested(provenance, output_path)
    return binding_freshness_errors(binding, path, f"{label} external output binding")


def evidence_record_errors(
    path: Path, record: dict[str, Any], cases: dict[str, CaseRecord]
) -> list[str]:
    """Return provenance, digest, and freshness errors for acceptance evidence."""

    label = f"acceptance evidence {path}"
    record_type = record.get("record_type")
    errors: list[str] = []
    if record_type != "cgl-lf-stage-i-direct-fast-case-acceptance":
        errors.extend(verify_self_digest(record, label))
    if record_type in {
        "cgl-lf-stage-i-direct-fast-case-acceptance",
        "cgl-lf-stage-i-direct-fast-campaign-evidence",
    }:
        errors.extend(direct_fast_output_provenance_errors(path, record, label))

    provenance = record.get("provenance")
    if not isinstance(provenance, dict):
        errors.append(f"{label} lacks provenance")
    else:
        if isinstance(record_type, str) and record_type.startswith("stage-i-scientific-"):
            for name in ("criteria", "criteria_review", "acceptance_utility"):
                if not isinstance(provenance.get(name), dict):
                    errors.append(f"{label} lacks required provenance binding {name}")
            inventory_name = (
                "case_evidence"
                if record_type == "stage-i-scientific-campaign-evidence"
                else "inputs"
            )
            if not isinstance(provenance.get(inventory_name), list):
                errors.append(
                    f"{label} lacks required provenance {inventory_name} inventory"
                )
        seen: set[tuple[str, str]] = set()
        for index, binding in enumerate(recursive_file_bindings(provenance)):
            key = (str(binding.get("path")), str(binding.get("sha256")))
            if key in seen:
                continue
            seen.add(key)
            errors.extend(verify_file_binding(binding, f"{label} provenance binding {index}"))

    case_id = record.get("case_id")
    case = cases.get(case_id) if isinstance(case_id, str) else None
    if record_type in {
        "stage-i-scientific-case-evidence",
        "cgl-lf-stage-i-direct-fast-case-acceptance",
    }:
        if case is None:
            errors.append(f"{label} does not identify a selected Stage I case")
            return errors
        if record_type == "stage-i-scientific-case-evidence":
            evaluation = record.get("evaluation_inputs")
            if not isinstance(evaluation, dict):
                errors.append(f"{label} lacks evaluation_inputs")
            else:
                for kind in ("mhd", "user"):
                    binding = evaluation.get(f"{kind}_history")
                    errors.extend(
                        binding_freshness_errors(
                            binding, case.history_paths.get(kind),
                            f"{label} {kind} history",
                        )
                    )
        else:
            if not isinstance(provenance, dict):
                return errors
            errors.extend(
                binding_freshness_errors(
                    provenance.get("lineage"), case.lineage_path, f"{label} lineage"
                )
            )
            histories = provenance.get("histories")
            if not isinstance(histories, dict):
                errors.append(f"{label} lacks history provenance")
            else:
                for kind in ("mhd", "user"):
                    errors.extend(
                        binding_freshness_errors(
                            histories.get(kind), case.history_paths.get(kind),
                            f"{label} {kind} history",
                        )
                    )
    return errors


def referenced_case_ids(path: Path, record: dict[str, Any]) -> set[str]:
    """Return case IDs unambiguously referenced by an audit path or inventory."""

    declared = record.get("case_id")
    if isinstance(declared, str) and declared in CASE_IDS:
        return {declared}
    strings = [str(path)]
    provenance = record.get("provenance")
    if isinstance(provenance, dict):
        patterns = provenance.get("input_patterns")
        if isinstance(patterns, list):
            strings.extend(str(value) for value in patterns if isinstance(value, str))
    snapshots = record.get("snapshots")
    if isinstance(snapshots, list):
        for snapshot in snapshots:
            if not isinstance(snapshot, dict):
                continue
            for key in ("root", "basename"):
                if isinstance(snapshot.get(key), str):
                    strings.append(str(snapshot[key]))
            rank_files = snapshot.get("rank_files")
            if isinstance(rank_files, list):
                strings.extend(
                    str(value["path"])
                    for value in rank_files
                    if isinstance(value, dict) and isinstance(value.get("path"), str)
                )
    directory_matches = {
        match.group(1)
        for value in strings
        for match in CASE_DIRECTORY_PATTERN.finditer(value)
    }
    if directory_matches:
        return directory_matches
    return {
        match.group(1)
        for value in strings
        for match in CASE_PATH_PATTERN.finditer(value)
    }


def audit_record_errors(path: Path, record: dict[str, Any]) -> list[str]:
    """Return provenance and retained-input freshness errors for one audit."""

    label = f"hyperbolicity audit {path}"
    provenance = record.get("provenance")
    if not isinstance(provenance, dict):
        return [f"{label} lacks provenance"]
    errors: list[str] = []
    for prefix in ("script", "bin_convert"):
        binding = {
            "path": provenance.get(f"{prefix}_path"),
            "sha256": provenance.get(f"{prefix}_sha256"),
        }
        errors.extend(verify_file_binding(binding, f"{label} {prefix}"))
    snapshots = record.get("snapshots")
    if not isinstance(snapshots, list) or not snapshots:
        errors.append(f"{label} lacks snapshot evidence")
        return errors
    for snapshot_index, snapshot in enumerate(snapshots):
        if not isinstance(snapshot, dict):
            errors.append(f"{label} snapshot {snapshot_index} is not an object")
            continue
        rank_files = snapshot.get("rank_files")
        if not isinstance(rank_files, list) or not rank_files:
            errors.append(f"{label} snapshot {snapshot_index} lacks rank-file inventory")
            continue
        try:
            inventory_digest = hashlib.sha256(canonical_json(rank_files)).hexdigest()
        except (TypeError, ValueError) as error:
            errors.append(f"{label} snapshot {snapshot_index} inventory is invalid: {error}")
            continue
        if inventory_digest != snapshot.get("input_inventory_sha256"):
            errors.append(f"{label} snapshot {snapshot_index} inventory digest differs")
        for rank_index, rank_file in enumerate(rank_files):
            rank_label = f"{label} snapshot {snapshot_index} rank file {rank_index}"
            if not isinstance(rank_file, dict):
                errors.append(f"{rank_label} is not an object")
                continue
            rank_path = binding_path(rank_file)
            if rank_path is None or not rank_path.is_file():
                errors.append(f"{rank_label} path is unavailable")
                continue
            try:
                stat = rank_path.stat()
            except OSError as error:
                errors.append(f"{rank_label} could not be inspected: {error}")
                continue
            if rank_file.get("size_bytes") != stat.st_size:
                errors.append(f"{rank_label} size differs")
            if rank_file.get("mtime_ns") != stat.st_mtime_ns:
                errors.append(f"{rank_label} mtime differs")
            if "sha256" in rank_file:
                errors.extend(verify_file_binding(rank_file, rank_label))
    if not referenced_case_ids(path, record):
        errors.append(f"{label} cannot be associated with a Stage I case")
    return errors


def validated_result(value: object) -> str:
    """Return one admitted evidence result without promoting unknown values."""

    return str(value) if value in EVIDENCE_RESULTS else "inconclusive"


def validated_science_contrast_result(value: object) -> str:
    """Return one admitted reviewed-science contrast disposition."""

    return str(value) if value in SCIENCE_CONTRAST_RESULTS else "inconclusive"


def science_record_errors(
    path: Path, record: dict[str, Any], analysis: Path
) -> list[str]:
    """Authenticate one reviewed direct-fast science aggregate."""

    label = f"reviewed science aggregate {path}"
    errors: list[str] = []
    if record.get("schema_version") != 1:
        errors.append(f"{label} has unsupported schema_version")
    if record.get("record_type") != SCIENCE_RECORD_TYPE:
        errors.append(f"{label} has unexpected record_type")
    if record.get("authority") != SCIENCE_AUTHORITY:
        errors.append(f"{label} has unexpected authority")
    if record.get("release_authorizing") is not False:
        errors.append(f"{label} must be explicitly non-release-authorizing")
    if record.get("result") not in EVIDENCE_RESULTS:
        errors.append(f"{label} has an invalid aggregate result")
    errors.extend(verify_self_digest(record, label))

    selected = record.get("selected_cases")
    selected_cases = selected if isinstance(selected, list) else []
    if (
        not isinstance(selected, list)
        or len(selected_cases) != len(set(selected_cases))
        or any(case_id not in CASE_IDS for case_id in selected_cases)
    ):
        errors.append(f"{label} has an invalid selected_cases inventory")
    if record.get("result") == "pass" and set(selected_cases) != set(CASE_IDS):
        errors.append(f"{label} partial selected_cases cannot pass")
    dispositions = record.get("case_dispositions")
    if not isinstance(dispositions, dict):
        errors.append(f"{label} lacks case_dispositions")
    elif set(dispositions) != set(selected_cases):
        errors.append(f"{label} case_dispositions differ from selected_cases")
    else:
        for case_id, disposition in dispositions.items():
            if not isinstance(disposition, dict):
                errors.append(f"{label} case disposition {case_id} is not an object")
                continue
            claim_eligible = disposition.get("claim_eligible")
            acceptance_result = disposition.get("acceptance_result")
            if not isinstance(claim_eligible, bool):
                errors.append(
                    f"{label} case disposition {case_id} lacks boolean claim eligibility"
                )
            if (
                acceptance_result is not None
                and acceptance_result not in EVIDENCE_RESULTS
            ):
                errors.append(
                    f"{label} case disposition {case_id} has an invalid acceptance result"
                )
            if claim_eligible is True and acceptance_result != "pass":
                errors.append(
                    f"{label} case disposition {case_id} is eligible without "
                    "acceptance pass"
                )

    provenance = record.get("provenance")
    if not isinstance(provenance, dict):
        errors.append(f"{label} lacks provenance")
    else:
        for name in (
            "inventory",
            "acceptance_provenance",
            "acceptance_campaign_evidence",
            "criteria",
            "criteria_review",
            "reviewed_acceptance_utility",
            "fast_report_utility",
            "paper_analyzer",
            "aggregator",
        ):
            if not isinstance(provenance.get(name), dict):
                errors.append(f"{label} lacks required provenance binding {name}")
        seen: set[tuple[str, str]] = set()
        for index, file_binding in enumerate(recursive_file_bindings(provenance)):
            key = (str(file_binding.get("path")), str(file_binding.get("sha256")))
            if key in seen:
                continue
            seen.add(key)
            errors.extend(
                verify_file_binding(file_binding, f"{label} provenance binding {index}")
            )
        errors.extend(
            binding_freshness_errors(
                provenance.get("inventory"),
                analysis / "inventory.json",
                f"{label} inventory",
            )
        )

    external_path = path.parent / "provenance.json"
    if not external_path.is_file():
        errors.append(f"{label} lacks external provenance: {external_path}")
    else:
        try:
            external = load_json(external_path)
        except (OSError, json.JSONDecodeError, PublicationError) as error:
            errors.append(f"{label} external provenance is unreadable: {error}")
        else:
            if (
                external.get("schema_version") != 1
                or external.get("record_type") != SCIENCE_PROVENANCE_RECORD_TYPE
            ):
                errors.append(f"{label} external provenance has unexpected schema")
            errors.extend(
                binding_freshness_errors(
                    nested(external, "outputs.science"),
                    path,
                    f"{label} external output binding",
                )
            )
            if isinstance(provenance, dict) and external.get("inputs") != provenance:
                errors.append(f"{label} external provenance inputs differ")

    families = record.get("families")
    if not isinstance(families, dict):
        errors.append(f"{label} lacks comparison families")
    else:
        for family, contrasts in families.items():
            if not isinstance(contrasts, dict):
                errors.append(f"{label} family {family} is not an object")
                continue
            for name, contrast in contrasts.items():
                if not isinstance(contrast, dict):
                    errors.append(f"{label} contrast {family}.{name} is not an object")
                    continue
                result = contrast.get("result")
                if result not in SCIENCE_CONTRAST_RESULTS:
                    errors.append(
                        f"{label} contrast {family}.{name} has an invalid result"
                    )
                references = [
                    contrast.get(key)
                    for key in ("active", "passive", "left", "right")
                    if contrast.get(key) is not None
                ]
                if any(case_id not in CASE_IDS for case_id in references):
                    errors.append(
                        f"{label} contrast {family}.{name} references an invalid case"
                    )
                if (
                    result == "pass"
                    and contrast.get("claim_eligible") is not True
                ):
                    errors.append(
                        f"{label} contrast {family}.{name} passes without "
                        "claim eligibility"
                    )
                if (
                    contrast.get("claim_eligible") is True
                    and isinstance(dispositions, dict)
                    and any(
                        not isinstance(dispositions.get(case_id), dict)
                        or dispositions[case_id].get("claim_eligible") is not True
                        for case_id in references
                    )
                ):
                    errors.append(
                        f"{label} contrast {family}.{name} is eligible with "
                        "ineligible contributing cases"
                    )
                metrics = contrast.get("metrics")
                if not isinstance(metrics, list):
                    errors.append(f"{label} contrast {family}.{name} lacks metrics")
                else:
                    for metric in metrics:
                        if (
                            isinstance(metric, dict)
                            and metric.get("holm_significant") is True
                            and metric.get("available") is not True
                        ):
                            errors.append(
                                f"{label} contrast {family}.{name} has a Holm result "
                                "without an available metric"
                            )
    gates = record.get("gates")
    if not isinstance(gates, list) or not gates:
        errors.append(f"{label} lacks required science gates")
    else:
        gate_results: list[str] = []
        for index, gate in enumerate(gates):
            if not isinstance(gate, dict):
                errors.append(f"{label} gate {index} is not an object")
                continue
            result = gate.get("result")
            if result not in {*EVIDENCE_RESULTS, "blocked_out_of_scope"}:
                errors.append(f"{label} gate {index} has an invalid result")
                continue
            gate_results.append(str(result))
            observations = gate.get("observations")
            if result == "pass" and isinstance(observations, list) and any(
                isinstance(observation, dict)
                and (
                    observation.get("available") is False
                    or observation.get("passed") is False
                    or observation.get("result") in {
                        "fail", "inconclusive", "blocked_out_of_scope"
                    }
                )
                for observation in observations
            ):
                errors.append(f"{label} gate {index} passes with adverse observations")
        expected = (
            "fail" if "fail" in gate_results
            else "inconclusive" if "inconclusive" in gate_results
            else "inconclusive" if (
                not gate_results
                or all(result == "blocked_out_of_scope" for result in gate_results)
            )
            else "pass"
        )
        if record.get("result") != expected:
            errors.append(f"{label} aggregate result differs from required gates")
    return errors


def ct_audit_record_errors(
    path: Path, record: dict[str, Any], analysis: Path
) -> list[str]:
    """Authenticate one non-authorizing direct CT audit aggregate."""

    label = f"direct CT audit {path}"
    errors: list[str] = []
    if record.get("schema_version") != 2:
        errors.append(f"{label} has unsupported schema_version")
    if record.get("record_type") != CT_AUDIT_RECORD_TYPE:
        errors.append(f"{label} has unexpected record_type")
    if record.get("result") not in EVIDENCE_RESULTS:
        errors.append(f"{label} has an invalid aggregate result")
    errors.extend(verify_record_digest(record, label, CT_EVIDENCE_DIGEST_METHOD))

    claim_boundary = record.get("claim_boundary")
    if not isinstance(claim_boundary, dict):
        errors.append(f"{label} lacks claim_boundary")
    elif (
        claim_boundary.get("campaign_authority_eligible") is not False
        or claim_boundary.get("release_authorizing") is not False
    ):
        errors.append(f"{label} must be explicitly non-authorizing")
    errors.extend(
        binding_freshness_errors(
            record.get("inventory"),
            analysis / "inventory.json",
            f"{label} inventory",
        )
    )
    source_bindings = record.get("source_bindings")
    if not isinstance(source_bindings, dict):
        errors.append(f"{label} lacks source_bindings")
    else:
        for index, file_binding in enumerate(recursive_file_bindings(source_bindings)):
            errors.extend(
                verify_file_binding(file_binding, f"{label} source binding {index}")
            )

    selection = record.get("selection")
    selected = selection.get("cases") if isinstance(selection, dict) else None
    selected_cases = selected if isinstance(selected, list) else []
    if (
        not isinstance(selected, list)
        or not selected_cases
        or len(selected_cases) != len(set(selected_cases))
        or any(case_id not in CASE_IDS for case_id in selected_cases)
    ):
        errors.append(f"{label} has an invalid selected-case inventory")
    cases = record.get("cases")
    if not isinstance(cases, dict):
        errors.append(f"{label} lacks cases")
        cases = {}
    elif set(cases) != set(selected_cases):
        errors.append(f"{label} cases differ from selected-case inventory")
    results: list[str] = []
    for case_id in selected_cases:
        case = cases.get(case_id)
        if not isinstance(case, dict):
            errors.append(f"{label} lacks selected case {case_id}")
            results.append("inconclusive")
            continue
        raw_result = case.get("ct_result")
        if raw_result not in EVIDENCE_RESULTS:
            errors.append(f"{label} case {case_id} has an invalid CT result")
        result = validated_result(raw_result)
        results.append(result)
        native = case.get("native_restart_ct")
        if not isinstance(native, dict):
            errors.append(f"{label} case {case_id} lacks native_restart_ct")
            continue
        if (
            native.get("campaign_authority_eligible") is not False
            or native.get("release_authorizing") is not False
        ):
            errors.append(f"{label} case {case_id} must remain non-authorizing")
        if result in {"pass", "fail"} and not (
            case.get("provenance_authenticated") is True
            and case.get("ct_evidence_available") is True
            and case.get("ct_claim_supported") is True
            and native.get("ct_evidence_available") is True
            and native.get("ct_claim_supported") is True
            and native.get("result") == result
        ):
            errors.append(
                f"{label} case {case_id} has a conclusive result without "
                "consistent CT evidence"
            )
        if result == "pass" and native.get("coverage_complete") is not True:
            errors.append(f"{label} case {case_id} passes without complete CT coverage")
        if result in {"pass", "fail"}:
            maximum = as_float(native.get("maximum_normalized_ct_divb"))
            threshold = as_float(native.get("normalized_ct_divb_lt"))
            numerically_consistent = (
                maximum is not None
                and threshold is not None
                and threshold > 0.0
                and native.get("numerical_result") == result
                and (
                    (maximum < threshold)
                    if result == "pass"
                    else (maximum >= threshold)
                )
            )
            if not numerically_consistent:
                errors.append(
                    f"{label} case {case_id} has a CT result inconsistent with "
                    "its numerical evidence"
                )
    expected = (
        "fail" if "fail" in results
        else "pass" if results and all(result == "pass" for result in results)
        else "inconclusive"
    )
    if record.get("result") != expected:
        errors.append(
            f"{label} aggregate result differs from selected-case CT results"
        )
    return errors


def atomic_write(path: Path, payload: bytes) -> None:
    """Atomically replace one output file."""

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
    """Write canonical human-readable JSON."""

    atomic_write(
        path,
        (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False)
            + "\n"
        ).encode("utf-8"),
    )


def write_text(path: Path, value: str) -> None:
    """Write one UTF-8 text artifact atomically."""

    atomic_write(path, value.encode("utf-8"))


def load_json(path: Path) -> dict[str, Any]:
    """Load one JSON object."""

    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise PublicationError(f"expected a JSON object: {path}")
    return value


def load_optional_json(
    path: Path, source_paths: set[Path], warnings: list[str]
) -> dict[str, Any] | None:
    """Load one optional JSON object while retaining a readable warning."""

    if not path.is_file():
        return None
    try:
        value = load_json(path)
    except (OSError, json.JSONDecodeError, PublicationError) as error:
        warnings.append(f"could not read {path}: {type(error).__name__}: {error}")
        return None
    source_paths.add(path.absolute())
    return value


def parse_history(path: Path) -> dict[str, list[float]]:
    """Parse a labeled Athena history file, ignoring a partial trailing row."""

    payload = path.read_bytes()
    labels: list[str] | None = None
    rows: list[list[float]] = []
    lines = payload.decode("utf-8").splitlines()
    for index, line in enumerate(lines):
        if line.startswith("#"):
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
                continue
            raise PublicationError(f"nonnumeric history row: {path}")
        if labels is not None and len(row) != len(labels):
            if index == len(lines) - 1 and not payload.endswith(b"\n"):
                continue
            raise PublicationError(f"history row width differs from header: {path}")
        rows.append(row)
    if labels is None or not rows or any(len(row) != len(labels) for row in rows):
        raise PublicationError(f"history lacks a complete labeled data table: {path}")
    return {
        label: [float(row[index]) for row in rows]
        for index, label in enumerate(labels)
    }


def history_path(
    case_dir: Path, lineage: dict[str, Any] | None, kind: str
) -> Path | None:
    """Select one merged MHD or user history from direct-fast report output."""

    declared = nested(lineage, f"histories.{kind}.path")
    if isinstance(declared, str) and Path(declared).is_file():
        return Path(declared)
    matches = sorted((case_dir / "history").glob(f"*.{kind}.hst"))
    return matches[0] if matches else None


def discover_acceptance_paths(analysis: Path, extras: Iterable[Path]) -> list[Path]:
    """Return deterministic candidate acceptance and audit JSON paths."""

    candidates: set[Path] = set()
    roots = [
        analysis / "scientific-acceptance",
        analysis / "acceptance",
        analysis / "audits",
        analysis / "scientific-audits",
        *extras,
    ]
    for root in roots:
        if root.is_file():
            candidates.add(root.absolute())
        elif root.is_dir():
            candidates.update(path.absolute() for path in root.rglob("*.json"))
    return sorted(candidates)


def contains_hyperbolicity_evidence(record: dict[str, Any]) -> bool:
    """Return whether one JSON record carries retained-state hyperbolicity data."""

    snapshots = record.get("snapshots")
    if isinstance(snapshots, list) and any(
        isinstance(snapshot, dict) and isinstance(snapshot.get("aggregate"), dict)
        for snapshot in snapshots
    ):
        return True
    paths = (
        "hyperbolicity",
        "retained_state_hyperbolicity",
        "negative_discriminant_fraction",
        "minimum_discriminant",
        "summary.negative_discriminant_fraction",
        "summary.minimum_discriminant",
    )
    return any(nested(record, path) is not None for path in paths)


def select_acceptance_record(
    records: list[tuple[Path, dict[str, Any]]],
    record_type: str,
    *,
    trusted_only: bool = True,
) -> dict[str, Any] | None:
    """Select the most information-rich deterministic trusted record of one type."""

    eligible = [
        (path, record)
        for path, record in records
        if record.get("record_type") == record_type
        and (
            record.get("_publication_evidence_validated") is True
            or not trusted_only
        )
    ]
    if not eligible:
        return None
    _, selected = max(
        eligible,
        key=lambda item: (
            len(item[1].get("gates", []))
            if isinstance(item[1].get("gates"), list) else 0,
            len(item[1].get("metrics", {}))
            if isinstance(item[1].get("metrics"), dict) else 0,
            str(item[0]),
        ),
    )
    return selected


def select_integrated_record(
    records: list[tuple[Path, dict[str, Any]]], *, kind: str
) -> dict[str, Any] | None:
    """Select one deterministic authenticated science or CT aggregate."""

    if not records:
        return None
    _, selected = max(
        records,
        key=lambda item: (
            len(item[1].get("selected_cases", []))
            if kind == "science" and isinstance(item[1].get("selected_cases"), list)
            else len(nested(item[1], "selection.cases") or [])
            if kind == "ct" and isinstance(nested(item[1], "selection.cases"), list)
            else 0,
            len(item[1].get("gates", []))
            if isinstance(item[1].get("gates"), list) else 0,
            str(item[0]),
        ),
    )
    return selected


def discover_data(analysis: Path, acceptance_paths: Iterable[Path]) -> PublicationData:
    """Load all available downstream evidence without requiring campaign completion."""

    source_paths: set[Path] = {RENDERER_PATH}
    warnings: list[str] = []
    aggregate = load_optional_json(
        analysis / "campaign/diagnostics.json", source_paths, warnings
    )
    comparisons = load_optional_json(
        analysis / "campaign/comparisons.json", source_paths, warnings
    )
    load_optional_json(analysis / "campaign/health.json", source_paths, warnings)
    load_optional_json(analysis / "inventory.json", source_paths, warnings)
    load_optional_json(analysis / "manifest.json", source_paths, warnings)

    cases: dict[str, CaseRecord] = {}
    for case_id in CASE_IDS:
        case_dir = analysis / "cases" / case_id
        diagnostics = load_optional_json(
            case_dir / "diagnostics.json", source_paths, warnings
        )
        lineage = load_optional_json(case_dir / "lineage.json", source_paths, warnings)
        model = nested(diagnostics, "model_choices")
        if not isinstance(model, dict):
            model = nested(lineage, "model_choices")
        if not isinstance(model, dict):
            model = {}
        record = CaseRecord(
            case_id=case_id,
            diagnostics=diagnostics,
            lineage=lineage,
            model=dict(model),
            lineage_path=case_dir / "lineage.json" if lineage is not None else None,
        )
        for kind in ("user", "mhd"):
            path = history_path(case_dir, lineage, kind)
            if path is None:
                continue
            try:
                parsed = parse_history(path)
                source_paths.add(path.absolute())
                if kind == "user":
                    record.user_history = parsed
                else:
                    record.mhd_history = parsed
                record.history_paths[kind] = path
            except (OSError, UnicodeDecodeError, PublicationError) as error:
                warnings.append(
                    f"{case_id} {kind} history unavailable: "
                    f"{type(error).__name__}: {error}"
                )
        cases[case_id] = record

    loaded_acceptance: list[tuple[Path, dict[str, Any]]] = []
    loaded_audits: list[tuple[Path, dict[str, Any]]] = []
    loaded_science: list[tuple[Path, dict[str, Any]]] = []
    loaded_ct_audits: list[tuple[Path, dict[str, Any]]] = []
    for path in discover_acceptance_paths(analysis, acceptance_paths):
        try:
            record = load_json(path)
        except (OSError, json.JSONDecodeError, PublicationError) as error:
            warnings.append(
                f"acceptance candidate unreadable {path}: "
                f"{type(error).__name__}: {error}"
            )
            continue
        source_paths.add(path.absolute())
        if record.get("record_type") == SCIENCE_RECORD_TYPE:
            validation_errors = science_record_errors(path, record, analysis)
            if validation_errors:
                warnings.extend(
                    f"rejected {path}: {error}" for error in validation_errors
                )
                continue
            loaded = dict(record)
            loaded["_publication_evidence_validated"] = True
            loaded["_publication_source_path"] = str(path)
            loaded_science.append((path, loaded))
            external = path.parent / "provenance.json"
            if external.is_file():
                source_paths.add(external.absolute())
            continue
        if record.get("record_type") == CT_AUDIT_RECORD_TYPE:
            validation_errors = ct_audit_record_errors(path, record, analysis)
            if validation_errors:
                warnings.extend(
                    f"rejected {path}: {error}" for error in validation_errors
                )
                continue
            loaded = dict(record)
            loaded["_publication_evidence_validated"] = True
            loaded["_publication_source_path"] = str(path)
            loaded_ct_audits.append((path, loaded))
            continue
        if record.get("record_type") in ACCEPTANCE_RECORD_TYPES:
            validation_errors = evidence_record_errors(path, record, cases)
            loaded = dict(record)
            loaded["_publication_evidence_validated"] = not validation_errors
            loaded["_publication_source_path"] = str(path)
            if validation_errors:
                warnings.extend(
                    f"rejected {path}: {error}" for error in validation_errors
                )
            loaded_acceptance.append((path, loaded))
            continue
        if contains_hyperbolicity_evidence(record):
            validation_errors = audit_record_errors(path, record)
            if validation_errors:
                warnings.extend(
                    f"rejected {path}: {error}" for error in validation_errors
                )
                continue
            loaded = dict(record)
            loaded["_publication_evidence_validated"] = True
            loaded["_publication_source_path"] = str(path)
            loaded["_publication_case_ids"] = sorted(referenced_case_ids(path, record))
            loaded_audits.append((path, loaded))
    for case_id, case in cases.items():
        matching = [
            (path, record)
            for path, record in loaded_acceptance
            if record.get("record_type") == "stage-i-scientific-case-evidence"
            and record.get("case_id") == case_id
            and record.get("_publication_evidence_validated") is True
        ]
        if len(matching) > 1:
            warnings.append(
                f"{case_id} has {len(matching)} acceptance records; "
                "selected the most information-rich deterministic record"
            )
        case.acceptance = select_acceptance_record(
            matching, "stage-i-scientific-case-evidence"
        )
        direct_matching = [
            (path, record)
            for path, record in loaded_acceptance
            if record.get("record_type")
            == "cgl-lf-stage-i-direct-fast-case-acceptance"
            and record.get("case_id") == case_id
            and record.get("_publication_evidence_validated") is True
        ]
        case.direct_acceptance = select_acceptance_record(
            direct_matching, "cgl-lf-stage-i-direct-fast-case-acceptance"
        )
    for path, record in loaded_acceptance:
        record_type = record.get("record_type")
        if record_type not in {
            "cgl-lf-stage-i-direct-fast-campaign-evidence",
            "stage-i-scientific-campaign-evidence",
        } or record.get("_publication_evidence_validated") is not True:
            continue
        case_results = record.get("case_results")
        if not isinstance(case_results, dict):
            record["_publication_evidence_validated"] = False
            warnings.append(f"rejected {path}: campaign evidence lacks case_results")
            continue
        selected_name = (
            "direct_acceptance"
            if record_type == "cgl-lf-stage-i-direct-fast-campaign-evidence"
            else "acceptance"
        )
        stale = [
            case_id
            for case_id, result in case_results.items()
            if case_id not in cases
            or not isinstance(getattr(cases[case_id], selected_name), dict)
            or getattr(cases[case_id], selected_name).get("result") != result
        ]
        if stale:
            record["_publication_evidence_validated"] = False
            warnings.append(
                f"rejected {path}: campaign evidence is stale for "
                f"{', '.join(sorted(stale))}"
            )
    campaign_acceptance = (
        select_acceptance_record(
            loaded_acceptance, "cgl-lf-stage-i-direct-fast-campaign-evidence"
        )
        or select_acceptance_record(
            loaded_acceptance, "stage-i-scientific-campaign-evidence"
        )
    )
    if campaign_acceptance is None:
        campaign_acceptance = (
            select_acceptance_record(
                loaded_acceptance,
                "cgl-lf-stage-i-direct-fast-campaign-evidence",
                trusted_only=False,
            )
            or select_acceptance_record(
                loaded_acceptance,
                "stage-i-scientific-campaign-evidence",
                trusted_only=False,
            )
        )
    if len(loaded_science) > 1:
        warnings.append(
            f"found {len(loaded_science)} authenticated reviewed science aggregates; "
            "selected the most information-rich deterministic record"
        )
    if len(loaded_ct_audits) > 1:
        warnings.append(
            f"found {len(loaded_ct_audits)} authenticated direct CT audits; "
            "selected the most information-rich deterministic record"
        )
    return PublicationData(
        analysis=analysis,
        cases=cases,
        aggregate=aggregate,
        comparisons=comparisons,
        campaign_acceptance=campaign_acceptance,
        science_record=select_integrated_record(loaded_science, kind="science"),
        ct_audit_record=select_integrated_record(loaded_ct_audits, kind="ct"),
        acceptance_records=[record for _, record in loaded_acceptance],
        audit_records=[record for _, record in loaded_audits],
        source_paths=source_paths,
        ingestion_warnings=warnings,
    )


def normalized_history_series(
    case: CaseRecord, name: str
) -> tuple[list[float], list[float]] | None:
    """Return a normalized user-history time series."""

    history = case.user_history
    if not isinstance(history, dict) or "time" not in history:
        return None
    times = history["time"]
    volume = history.get("volume", [1.0] * len(times))
    if name == "unstable":
        if "mirror_vol" not in history or "fire_vol" not in history:
            return None
        values = [
            (mirror + firehose) / max(abs(vol), 1.0e-300)
            for mirror, firehose, vol in zip(
                history["mirror_vol"], history["fire_vol"], volume
            )
        ]
    elif name == "hard_bound":
        if "hard_vol" not in history:
            return None
        values = [
            value / max(abs(vol), 1.0e-300)
            for value, vol in zip(history["hard_vol"], volume)
        ]
    elif name == "parallel_forcing_fraction":
        if "force_prp2" not in history or "force_prl2" not in history:
            return None
        values = [
            parallel / max(parallel + perpendicular, 1.0e-300)
            for perpendicular, parallel in zip(
                history["force_prp2"], history["force_prl2"]
            )
        ]
    elif name in history:
        values = [
            value / max(abs(vol), 1.0e-300)
            for value, vol in zip(history[name], volume)
        ]
    else:
        return None
    return list(times), values


def interpolate_at(times: list[float], values: list[float], target: float) -> float:
    """Linearly interpolate one ordered time series."""

    if target < times[0] - TIME_TOLERANCE or target > times[-1] + TIME_TOLERANCE:
        raise PublicationError("requested time is outside the available history")
    for index, time in enumerate(times):
        if math.isclose(time, target, rel_tol=0.0, abs_tol=TIME_TOLERANCE):
            return values[index]
        if time > target and index > 0:
            left_t, right_t = times[index - 1], time
            fraction = (target - left_t) / (right_t - left_t)
            return values[index - 1] + fraction * (values[index] - values[index - 1])
    return values[-1]


def time_weighted_mean(
    times: list[float], values: list[float], start: float = 4.0, end: float = 10.0
) -> float | None:
    """Return a trapezoidal exact-window mean when the complete window exists."""

    if len(times) != len(values) or len(times) < 2:
        return None
    if times[0] > start + TIME_TOLERANCE or times[-1] < end - TIME_TOLERANCE:
        return None
    selected_times = [start]
    selected_values = [interpolate_at(times, values, start)]
    for time, value in zip(times, values):
        if start + TIME_TOLERANCE < time < end - TIME_TOLERANCE:
            selected_times.append(time)
            selected_values.append(value)
    selected_times.append(end)
    selected_values.append(interpolate_at(times, values, end))
    integral = sum(
        0.5 * (right_v + left_v) * (right_t - left_t)
        for left_t, right_t, left_v, right_v in zip(
            selected_times,
            selected_times[1:],
            selected_values,
            selected_values[1:],
        )
    )
    return integral / (end - start)


def diagnostic_steady_metric(case: CaseRecord, metric: str) -> float | None:
    """Return one fast-report developed-window metric."""

    window = nested(case.diagnostics, "windows.steady.history.analysis_window")
    if not isinstance(window, dict):
        return None
    aliases = {
        "kinetic": "kinetic_mean",
        "magnetic": "magnetic_mean",
        "abs_dp": "abs_dp_mean",
        "beta": "beta_mean",
        "nu_eff": "nu_eff_mean",
        "hard_bound": "hard_vol_fraction_mean",
        "mirror": "mirror_vol_fraction_mean",
        "firehose": "fire_vol_fraction_mean",
    }
    if metric == "unstable":
        mirror = as_float(window.get("mirror_vol_fraction_mean"))
        firehose = as_float(window.get("fire_vol_fraction_mean"))
        return mirror + firehose if mirror is not None and firehose is not None else None
    return as_float(window.get(aliases.get(metric, metric)))


def acceptance_metric(case: CaseRecord, metric: str) -> float | None:
    """Return one exact-window scientific-acceptance metric when present."""

    if not isinstance(case.acceptance, dict):
        return None
    aliases = {
        "hard_bound": "hard_bound_occupancy",
        "unstable": "unstable_occupancy",
    }
    names = (metric, aliases.get(metric, metric))
    for name in names:
        for root in ("metrics", "analyzer_metrics"):
            record = nested(case.acceptance, f"{root}.{name}")
            if not isinstance(record, dict):
                continue
            for window in ("full", "late"):
                value = nested(record, f"{window}.mean")
                parsed = as_float(value)
                if parsed is not None:
                    return parsed
            parsed = as_float(record.get("mean"))
            if parsed is not None:
                return parsed
    return None


def metric_value(case: CaseRecord, metric: str) -> float | None:
    """Return the best available developed-window scalar for one case."""

    value = diagnostic_steady_metric(case, metric)
    if value is not None:
        return value
    value = acceptance_metric(case, metric)
    if value is not None:
        return value
    series = normalized_history_series(case, metric)
    return time_weighted_mean(*series) if series is not None else None


def scientific_response_eligible(data: PublicationData, case: CaseRecord) -> bool:
    """Return whether one case may populate scientific-response products."""

    if completion_status(case) != "pass":
        return False
    direct_health = nested(case.direct_acceptance, "health.result")
    if isinstance(direct_health, str) and direct_health != "pass":
        return False
    fast_health = fast_health_status(case)
    if fast_health in {"structural_error", "fail"}:
        return False
    return direct_health == "pass" or fast_health in {"clean", "warnings"}


def response_metric_value(
    data: PublicationData, case: CaseRecord, metric: str
) -> float | None:
    """Return a developed-window scalar only for a numerically healthy case."""

    return metric_value(case, metric) if scientific_response_eligible(data, case) else None


def lf_ledger_metric(case: CaseRecord, kind: str) -> float | None:
    """Return one developed-window applied-work ledger magnitude."""

    if kind == "heat_flux_work":
        paths = (
            "windows.steady.lf_history.applied_heat_flux_work.total",
            "windows.steady.lf_history.applied_heat_flux_work.parallel",
        )
    elif kind == "pressure_work":
        paths = (
            "windows.steady.lf_history.applied_pressure_work.total",
            "windows.steady.lf_history.applied_pressure_work.anisotropic",
        )
    else:
        raise PublicationError(f"unsupported LF ledger metric: {kind}")
    values = [
        abs(value)
        for value in (as_float(nested(case.diagnostics, path)) for path in paths)
        if value is not None
    ]
    return max(values) if values else None


def response_lf_ledger_metric(
    data: PublicationData, case: CaseRecord, kind: str
) -> float | None:
    """Return an LF ledger response only for a numerically healthy case."""

    return (
        lf_ledger_metric(case, kind)
        if scientific_response_eligible(data, case)
        else None
    )


def final_time(case: CaseRecord) -> float | None:
    """Return the best available physical final time."""

    for value in (
        nested(case.diagnostics, "health.final_time"),
        nested(case.lineage, "final_time"),
        (
            case.user_history.get("time", [None])[-1]
            if isinstance(case.user_history, dict) else None
        ),
    ):
        parsed = as_float(value)
        if parsed is not None:
            return parsed
    return None


def assembly_status(case: CaseRecord) -> str:
    """Return the best available direct-fast assembly status."""

    for value in (
        nested(case.diagnostics, "assembly_status"),
        nested(case.lineage, "status"),
    ):
        if isinstance(value, str):
            return value
    return "unavailable"


def completion_status(case: CaseRecord) -> str:
    """Return a normalized completion status."""

    time = final_time(case)
    if time is not None and time >= TARGET_TIME - TIME_TOLERANCE:
        return "pass"
    status = assembly_status(case)
    if status in {"failed", "error", "structural_error"}:
        return "fail"
    return "incomplete" if time is not None else "unknown"


def fast_health_status(case: CaseRecord) -> str:
    """Return a normalized fast-report health status."""

    result = nested(case.diagnostics, "health.result")
    if result in ("clean", "warnings", "structural_error"):
        return str(result)
    return "unknown"


def case_warning_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return distinct case-level structural, numerical, and scientific warnings."""

    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    paths = (
        ("structural", "health.structural_warnings"),
        ("numerical", "health.numerical_warnings"),
        ("scientific", "health.science_warnings"),
        ("analysis", "analysis_warnings"),
        ("lineage", "warnings"),
    )
    for case_id in CASE_IDS:
        case = data.cases[case_id]
        for category, path in paths:
            for root in (case.diagnostics, case.lineage):
                values = nested(root, path)
                if not isinstance(values, list):
                    continue
                for value in values:
                    if not isinstance(value, str):
                        continue
                    key = (case_id, category, value)
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append({
                        "case_id": case_id,
                        "category": category,
                        "warning": value,
                    })
    return sorted(
        rows,
        key=lambda row: (
            str(row["case_id"]), str(row["category"]), str(row["warning"])
        ),
    )


def fatal_counter_status(case: CaseRecord) -> str:
    """Return whether fatal LF counters are observed zero."""

    maxima = nested(case.diagnostics, "health.fatal_lf_failure_maxima")
    if isinstance(maxima, dict):
        values = [as_float(value) for value in maxima.values()]
        finite = [value for value in values if value is not None]
        if finite:
            return "pass" if max(finite) == 0.0 else "fail"
    history = case.mhd_history
    if isinstance(history, dict):
        columns = ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos")
        if all(column in history for column in columns):
            maximum = max(
                abs(value)
                for column in columns
                for value in history[column]
            )
            return "pass" if maximum == 0.0 else "fail"
    return "unknown"


def acceptance_status(data: PublicationData, case: CaseRecord) -> str:
    """Return the selected non-authorizing scientific-acceptance result."""

    direct_result = (
        case.direct_acceptance.get("result")
        if isinstance(case.direct_acceptance, dict)
        else None
    )
    if isinstance(direct_result, str):
        return direct_result
    if (
        isinstance(case.acceptance, dict)
        and isinstance(case.acceptance.get("result"), str)
    ):
        return str(case.acceptance["result"])
    result = (
        nested(data.campaign_acceptance, f"case_results.{case.case_id}")
        if isinstance(data.campaign_acceptance, dict)
        and data.campaign_acceptance.get("_publication_evidence_validated") is True
        else None
    )
    return str(result) if isinstance(result, str) else "unknown"


def scope_status(case_id: str) -> str:
    """Return the fixed campaign claim scope."""

    return "restricted" if case_id in {"R10", "R14", "R15"} else "configuration"


def science_case_status(data: PublicationData, case_id: str) -> str:
    """Return reviewed-science eligibility for one selected case."""

    disposition = nested(data.science_record, f"case_dispositions.{case_id}")
    if not isinstance(disposition, dict):
        return "unknown"
    if (
        case_id == "R15"
        and selected_nonfatal_hard_bound_variant(data.cases["R15"])
    ):
        return "restricted"
    if (
        case_id == "R15"
        and r15_strict_failure_disposition(data) != "unavailable/inconclusive"
    ):
        return "fail"
    if disposition.get("claim_eligible") is True:
        return "pass"
    return (
        "fail"
        if disposition.get("acceptance_result") == "fail"
        else "inconclusive"
    )


def ct_case_record(data: PublicationData, case_id: str) -> dict[str, Any] | None:
    """Return one authenticated direct CT case record."""

    record = nested(data.ct_audit_record, f"cases.{case_id}")
    return record if isinstance(record, dict) else None


def ct_case_status(data: PublicationData, case_id: str) -> str:
    """Return one authenticated numerical CT result without granting authority."""

    record = ct_case_record(data, case_id)
    return validated_result(record.get("ct_result")) if record is not None else "unknown"


def aggregate_science_status(data: PublicationData) -> str:
    """Return the authenticated reviewed-science aggregate result."""

    return (
        validated_result(data.science_record.get("result"))
        if isinstance(data.science_record, dict) else "unknown"
    )


def aggregate_ct_status(data: PublicationData) -> str:
    """Return the authenticated direct CT aggregate numerical result."""

    return (
        validated_result(data.ct_audit_record.get("result"))
        if isinstance(data.ct_audit_record, dict) else "unknown"
    )


def ct_campaign_coverage_complete(data: PublicationData) -> bool:
    """Return whether every Stage I case has complete authenticated CT coverage."""

    selected = nested(data.ct_audit_record, "selection.cases")
    return (
        isinstance(selected, list)
        and set(selected) == set(CASE_IDS)
        and all(
            nested(
                data.ct_audit_record,
                f"cases.{case_id}.native_restart_ct.coverage_complete",
            )
            is True
            for case_id in CASE_IDS
        )
    )


def r15_science_scope(data: PublicationData, *labels: object) -> str:
    """Return the publication scope for results involving selected diagnostic R15."""

    if (
        any("R15" in str(label) for label in labels)
        and selected_nonfatal_hard_bound_variant(data.cases["R15"])
    ):
        return "restricted nonfatal-hard-bound diagnostic; not strict R15 success"
    return "standard reviewed-science scope"


def r15_strict_failure_disposition(data: PublicationData) -> str:
    """Return authenticated strict-R15 failure evidence without inventing details."""

    case = data.cases["R15"]
    roots: list[dict[str, Any]] = []
    for record in (case.direct_acceptance, case.acceptance):
        if (
            isinstance(record, dict)
            and record.get("_publication_evidence_validated") is True
        ):
            roots.append(record)
    if (
        isinstance(data.science_record, dict)
        and data.science_record.get("_publication_evidence_validated") is True
    ):
        binding = nested(data.science_record, "provenance.case_acceptance.R15")
        path = binding_path(binding)
        if path is not None and not verify_file_binding(
            binding, "R15 science-bound case acceptance"
        ):
            try:
                roots.append(load_json(path))
                data.source_paths.add(path.absolute())
            except (OSError, json.JSONDecodeError, PublicationError):
                pass
    candidates: list[dict[str, Any]] = []
    for root in roots:
        if root.get("record_type") == R15_STRICT_FAILURE_RECORD_TYPE:
            candidates.append(root)
        for path in (
            "retained_strict_failure_evidence",
            "scope.retained_strict_failure_evidence",
            "claim_scope.retained_strict_failure_evidence",
        ):
            values = nested(root, path)
            if isinstance(values, list):
                candidates.extend(value for value in values if isinstance(value, dict))
    authenticated: list[dict[str, Any]] = []
    for record in candidates:
        counters = record.get("failure_counters")
        provenance = record.get("provenance")
        if (
            record.get("schema_version") != 1
            or record.get("record_type") != R15_STRICT_FAILURE_RECORD_TYPE
            or record.get("case_id") != "R15"
            or record.get("result") != "fail"
            or record.get("strict_admissibility_evidence") is not True
            or as_float(record.get("failure_time")) is None
            or not isinstance(counters, dict)
            or (as_float(counters.get("lf_hardbd")) or 0.0) <= 0.0
            or not isinstance(provenance, dict)
            or any(
                not isinstance(provenance.get(name), dict)
                for name in ("manifest", "run_exit_code", "slurm_log")
            )
            or any(
                verify_file_binding(
                    provenance.get(name), f"R15 strict-failure {name} provenance"
                )
                for name in ("manifest", "run_exit_code", "slurm_log")
            )
        ):
            continue
        authenticated.append(record)
    if not authenticated:
        return "unavailable/inconclusive"
    selected = min(
        authenticated,
        key=lambda record: (
            float(record["failure_time"]),
            str(record.get("job_id", "")),
        ),
    )
    time = float(selected["failure_time"])
    counters = selected["failure_counters"]
    hard_bound = float(counters["lf_hardbd"])
    job_id = selected.get("job_id")
    details = ["fail", f"t={time:.8g}", f"hard_bound={hard_bound:.8g}"]
    if isinstance(job_id, (str, int)):
        details.append(f"job={job_id}")
    return "; ".join(details)


def publication_evidence_state(data: PublicationData) -> str:
    """Return an explicit renderer-wide completeness label."""

    if any(completion_status(data.cases[case_id]) != "pass" for case_id in CASE_IDS):
        return "partial/transient"
    if (
        data.ingestion_warnings
        or aggregate_science_status(data) in {"unknown", "inconclusive"}
        or aggregate_ct_status(data) in {"unknown", "inconclusive"}
        or not ct_campaign_coverage_complete(data)
        or any(
        acceptance_status(data, data.cases[case_id]) in {"unknown", "inconclusive"}
        for case_id in CASE_IDS
        )
    ):
        return "complete integration / partial evidence"
    return "complete"


def display_status(status: str) -> str:
    """Return a compact status label for figure cells."""

    labels = {
        "pass": "pass",
        "clean": "clean",
        "warnings": "warn",
        "warning": "warn",
        "restricted": "restricted",
        "configuration": "configured",
        "incomplete": "partial",
        "inconclusive": "inconclusive",
        "blocked_out_of_scope": "blocked",
        "fail": "fail",
        "structural_error": "error",
        "unknown": "--",
    }
    return labels.get(status, status)


def status_color(status: str) -> str:
    """Return the fixed color for one status."""

    return STATUS_COLORS.get(status, STATUS_COLORS["unknown"])


def tex_escape(value: object) -> str:
    """Escape one table value for a LaTeX tabular."""

    text = text_value(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "_": r"\_",
        "#": r"\#",
        "$": r"\$",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(character, character) for character in text)


def write_table(
    directory: Path,
    name: str,
    columns: list[str],
    rows: list[dict[str, object]],
) -> list[Path]:
    """Write deterministic CSV and LaTeX versions of one table."""

    directory.mkdir(parents=True, exist_ok=True)
    csv_buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        csv_buffer, fieldnames=columns, extrasaction="ignore", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({column: text_value(row.get(column)) for column in columns})
    csv_path = directory / f"{name}.csv"
    write_text(csv_path, csv_buffer.getvalue())

    alignment = "l" * len(columns)
    lines = [
        r"\begin{tabular}{" + alignment + "}",
        r"\hline",
        " & ".join(tex_escape(column) for column in columns) + r" \\",
        r"\hline",
    ]
    lines.extend(
        " & ".join(tex_escape(row.get(column)) for column in columns) + r" \\"
        for row in rows
    )
    lines.extend([r"\hline", r"\end{tabular}", ""])
    tex_path = directory / f"{name}.tex"
    write_text(tex_path, "\n".join(lines))
    return [csv_path, tex_path]


def configure_plotting() -> tuple[Any, Any, Any]:
    """Load Matplotlib with deterministic manuscript-oriented defaults."""

    os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "figure.dpi": 120,
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "legend.fontsize": 7,
        "lines.linewidth": 1.25,
        "pdf.compression": 9,
        "pdf.fonttype": 42,
        "savefig.transparent": False,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    })
    import matplotlib.colors as colors
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    return plt, colors, patches


def save_figure(fig: Any, path: Path) -> None:
    """Atomically save one deterministic PDF figure."""

    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=".pdf", dir=path.parent, prefix=f".{path.stem}.", delete=False
        ) as stream:
            staged = Path(stream.name)
        fig.savefig(
            staged,
            format="pdf",
            bbox_inches="tight",
            metadata={
                "Title": path.stem,
                "Author": "AthenaK CGL-LF Stage I",
                "Creator": "cgl_lf_stage_i_fast_publication.py",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        os.replace(staged, path)
        staged = None
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def status_matrix_figure(
    plt: Any,
    colors: Any,
    patches: Any,
    row_labels: list[str],
    column_labels: list[str],
    statuses: list[list[str]],
    cell_text: list[list[str]],
    title: str,
    path: Path,
    footnote: str | None = None,
) -> None:
    """Render one categorical status matrix."""

    palette = [
        STATUS_COLORS["pass"],
        STATUS_COLORS["warning"],
        STATUS_COLORS["restricted"],
        STATUS_COLORS["configuration"],
        STATUS_COLORS["incomplete"],
        STATUS_COLORS["inconclusive"],
        STATUS_COLORS["fail"],
        STATUS_COLORS["unknown"],
    ]
    categories = {
        "pass": 0, "clean": 0,
        "warning": 1, "warnings": 1,
        "restricted": 2,
        "configuration": 3,
        "incomplete": 4,
        "inconclusive": 5, "blocked_out_of_scope": 5,
        "fail": 6, "structural_error": 6,
        "unknown": 7,
    }
    matrix = [
        [categories.get(status, categories["unknown"]) for status in row]
        for row in statuses
    ]
    height = max(3.0, 0.34 * len(row_labels) + 1.8)
    width = max(8.0, 1.45 * len(column_labels))
    fig, axis = plt.subplots(figsize=(width, height))
    cmap = colors.ListedColormap(palette)
    norm = colors.BoundaryNorm([value - 0.5 for value in range(9)], cmap.N)
    axis.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    axis.set_xticks(range(len(column_labels)), column_labels)
    axis.set_yticks(range(len(row_labels)), row_labels)
    axis.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
    for row, values in enumerate(cell_text):
        for column, value in enumerate(values):
            axis.text(column, row, value, ha="center", va="center", fontsize=6.8)
    axis.set_xticks(
        [value - 0.5 for value in range(1, len(column_labels))], minor=True
    )
    axis.set_yticks([value - 0.5 for value in range(1, len(row_labels))], minor=True)
    axis.grid(which="minor", color="white", linewidth=1.0)
    axis.tick_params(which="minor", bottom=False, left=False)
    axis.set_title(title, pad=30)
    legend = [
        patches.Patch(facecolor=STATUS_COLORS[name], label=label)
        for name, label in (
            ("pass", "pass / clean"),
            ("warning", "warning"),
            ("restricted", "restricted scope"),
            ("configuration", "configuration / policy"),
            ("incomplete", "partial"),
            ("inconclusive", "inconclusive"),
            ("fail", "fail / error"),
            ("unknown", "not available"),
        )
    ]
    fig.legend(
        handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.065),
        ncol=4, frameon=False,
    )
    if footnote:
        wrapped_footnote = textwrap.fill(
            footnote, width=max(82, 20 * len(column_labels))
        )
        fig.text(
            0.02, 0.012, wrapped_footnote, ha="left", va="bottom", fontsize=6.8
        )
        footnote_lines = wrapped_footnote.count("\n") + 1
    else:
        footnote_lines = 0
    fig.subplots_adjust(
        bottom=(0.22 + 0.025 * footnote_lines) if footnote else 0.18,
        top=0.86,
    )
    save_figure(fig, path)
    plt.close(fig)


def science_contrast_rows(data: PublicationData) -> list[dict[str, object]]:
    """Flatten authenticated reviewed-science contrasts, including Holm results."""

    rows: list[dict[str, object]] = []
    families = (
        data.science_record.get("families")
        if isinstance(data.science_record, dict) else None
    )
    if not isinstance(families, dict):
        return rows
    for family, contrasts in sorted(families.items()):
        if not isinstance(contrasts, dict):
            continue
        for name, contrast in sorted(contrasts.items()):
            if not isinstance(contrast, dict):
                continue
            left = contrast.get("active", contrast.get("left"))
            right = contrast.get("passive", contrast.get("right"))
            metrics = contrast.get("metrics")
            if not isinstance(metrics, list) or not metrics:
                metrics = [{}]
            for metric in metrics:
                if not isinstance(metric, dict):
                    continue
                rows.append({
                    "family": family,
                    "contrast": name,
                    "left": left,
                    "right": right,
                    "result": validated_science_contrast_result(
                        contrast.get("result")
                    ),
                    "claim_eligible": contrast.get("claim_eligible") is True,
                    "metric": metric.get("metric"),
                    "available": metric.get("available") is True,
                    "left_mean": metric.get("left_mean"),
                    "right_mean": metric.get("right_mean"),
                    "difference": metric.get(
                        "difference", metric.get("difference_left_minus_right")
                    ),
                    "combined_standard_error": metric.get("combined_standard_error"),
                    "z_score": metric.get("z_score"),
                    "two_sided_p": metric.get("two_sided_p"),
                    "standardized_effect": metric.get("standardized_effect"),
                    "holm_threshold": metric.get("holm_threshold"),
                    "holm_significant": (
                        metric.get("holm_significant")
                        if metric.get("available") is True else None
                    ),
                    "reason": metric.get("reason", contrast.get("reason")),
                    "claim_scope": r15_science_scope(data, name, left, right),
                    "authority": SCIENCE_AUTHORITY,
                    "release_authorizing": False,
                })
    return rows


def science_gate_rows(data: PublicationData) -> list[dict[str, object]]:
    """Flatten authenticated reviewed-science gates."""

    gates = (
        data.science_record.get("gates")
        if isinstance(data.science_record, dict) else None
    )
    if not isinstance(gates, list):
        return []
    return [
        {
            "gate": gate.get("name"),
            "result": validated_result(gate.get("result")),
            "reason": gate.get("reason"),
            "claim_scope": r15_science_scope(data, gate.get("name")),
            "observations": gate.get("observations"),
            "limits": gate.get("limits"),
            "authority": SCIENCE_AUTHORITY,
            "release_authorizing": False,
        }
        for gate in gates if isinstance(gate, dict)
    ]


def science_resolution_rows(data: PublicationData) -> list[dict[str, object]]:
    """Flatten authenticated reviewed common-range resolution criteria."""

    resolution = (
        data.science_record.get("resolution")
        if isinstance(data.science_record, dict) else None
    )
    observations = (
        resolution.get("observations") if isinstance(resolution, dict) else None
    )
    if not isinstance(observations, list):
        return []
    return [
        {
            "result": validated_result(resolution.get("result")),
            "kind": value.get("kind"),
            "name": value.get("product", value.get("metric")),
            "available": value.get("available") is True,
            "passed": value.get("passed") if value.get("available") is True else None,
            "observations": value,
            "limits": resolution.get("limits"),
            "reason": value.get("reason", resolution.get("reason")),
            "authority": SCIENCE_AUTHORITY,
            "release_authorizing": False,
        }
        for value in observations if isinstance(value, dict)
    ]


def science_mks24_rows(data: PublicationData) -> list[dict[str, object]]:
    """Flatten authenticated admitted MKS24 residual and drift criteria."""

    mks24 = (
        data.science_record.get("mks24")
        if isinstance(data.science_record, dict) else None
    )
    panels = mks24.get("panels") if isinstance(mks24, dict) else None
    if not isinstance(panels, dict):
        return []
    rows: list[dict[str, object]] = []
    for panel_id, panel in sorted(panels.items()):
        products = panel.get("products") if isinstance(panel, dict) else None
        if not isinstance(products, list):
            continue
        for product in products:
            if not isinstance(product, dict):
                continue
            observations = product.get("observations")
            rows.append({
                "panel": panel_id,
                "panel_result": validated_result(panel.get("result")),
                "product_id": product.get("product_id"),
                "case_id": product.get("case_id"),
                "source": product.get("source"),
                "result": validated_result(product.get("result")),
                "normalized_residual_rms": nested(
                    observations, "normalized_residual_rms"
                ),
                "maximum_absolute_normalized_residual": nested(
                    observations, "maximum_absolute_normalized_residual"
                ),
                "early_late_vector_drift_rms": nested(
                    observations, "early_late_vector_drift_rms"
                ),
                "limits": product.get("limits"),
                "reason": product.get("reason"),
                "authority": SCIENCE_AUTHORITY,
                "release_authorizing": False,
            })
    return rows


def ct_health_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return authenticated direct CT numerical health without authority promotion."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        record = ct_case_record(data, case_id)
        native = record.get("native_restart_ct") if isinstance(record, dict) else None
        rows.append({
            "case_id": case_id,
            "ct_result": ct_case_status(data, case_id),
            "ct_evidence_available": (
                record.get("ct_evidence_available") is True
                if isinstance(record, dict) else False
            ),
            "ct_claim_supported": (
                record.get("ct_claim_supported") is True
                if isinstance(record, dict) else False
            ),
            "coverage_complete": (
                native.get("coverage_complete") is True
                if isinstance(native, dict) else False
            ),
            "maximum_normalized_ct_divb": (
                native.get("maximum_normalized_ct_divb")
                if isinstance(native, dict) else None
            ),
            "normalized_ct_divb_lt": (
                native.get("normalized_ct_divb_lt")
                if isinstance(native, dict) else None
            ),
            "campaign_authority_eligible": False,
            "release_authorizing": False,
            "reason": record.get("reason") if isinstance(record, dict) else (
                "case not selected by the authenticated direct CT audit"
                if isinstance(data.ct_audit_record, dict)
                else "no authenticated direct CT audit selected"
            ),
        })
    return rows


def health_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return the complete health/completion table."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        case = data.cases[case_id]
        health = nested(case.diagnostics, "health")
        if not isinstance(health, dict):
            health = {}
        rows.append({
            "case_id": case_id,
            "assembly_status": assembly_status(case),
            "final_time": final_time(case),
            "target_time": TARGET_TIME,
            "completion": completion_status(case),
            "fast_health": fast_health_status(case),
            "fatal_lf_counters": fatal_counter_status(case),
            "acceptance": acceptance_status(data, case),
            "reviewed_science": science_case_status(data, case_id),
            "direct_ct_numerical": ct_case_status(data, case_id),
            "direct_ct_release_authorizing": False,
            "claim_scope": (
                "restricted" if scope_status(case_id) == "restricted" else "standard"
            ),
            "strict_admissibility": case.model.get("cgl_lf_strict_admissibility"),
            "hard_bound_maximum": health.get("hard_bound_diagnostic_maximum"),
            "hard_bound_volume_maximum": health.get("hard_bound_volume_maximum"),
            "structural_error_count": len(health.get("structural_errors", []))
            if isinstance(health.get("structural_errors"), list) else None,
            "structural_warning_count": len(health.get("structural_warnings", []))
            if isinstance(health.get("structural_warnings"), list) else None,
            "numerical_warning_count": len(health.get("numerical_warnings", []))
            if isinstance(health.get("numerical_warnings"), list) else None,
            "science_warning_count": len(health.get("science_warnings", []))
            if isinstance(health.get("science_warnings"), list) else None,
        })
    return rows


def render_health(
    data: PublicationData, plt: Any, colors: Any, patches: Any, path: Path
) -> None:
    """Render the R02-R17 completion and health strip."""

    statuses: list[list[str]] = []
    labels: list[list[str]] = []
    for case_id in CASE_IDS:
        case = data.cases[case_id]
        completion = completion_status(case)
        health = fast_health_status(case)
        fatal = fatal_counter_status(case)
        acceptance = acceptance_status(data, case)
        science = science_case_status(data, case_id)
        ct = ct_case_status(data, case_id)
        scope = scope_status(case_id)
        statuses.append([completion, health, fatal, acceptance, science, ct, scope])
        time = final_time(case)
        labels.append([
            f"{time:.2f}/10" if time is not None else "--",
            display_status(health),
            "zero" if fatal == "pass" else display_status(fatal),
            display_status(acceptance),
            display_status(science),
            display_status(ct),
            "restricted" if scope == "restricted" else "standard",
        ])
    status_matrix_figure(
        plt, colors, patches, list(CASE_IDS),
        [
            "Completion", "Fast health", "Fatal LF", "Acceptance",
            "Reviewed science", "Direct CT", "Claim scope",
        ],
        statuses, labels,
        (
            "Stage I completion, numerical health, and claim eligibility "
            f"[{publication_evidence_state(data)} evidence]"
        ),
        path,
        (
            "Acceptance, reviewed science, and direct CT are non-authorizing "
            "direct-fast evidence. R10, R14, and R15 remain restricted even when "
            "completion and numerical-health cells pass. "
            "Claim scope and strict-policy settings are configurations, not passes; "
            "Fatal LF excludes selected nonfatal hard-bound diagnostic variants. "
            "Any R15 strict-run failure details are displayed only when carried by "
            "authenticated evidence."
        ),
    )


def render_active_passive(data: PublicationData, plt: Any, path: Path) -> None:
    """Render matched active/passive history comparisons."""

    fig, axes = plt.subplots(4, 2, figsize=(9.0, 10.2), sharex=True)
    metrics = (("abs_dp", r"$\langle|\Delta p|\rangle$"),
               ("unstable", "mirror + firehose volume fraction"))
    for row, (active, passive) in enumerate(ACTIVE_PASSIVE_PAIRS):
        for column, (metric, ylabel) in enumerate(metrics):
            axis = axes[row, column]
            available = 0
            for case_id, style, label in (
                (active, "-", f"{active} active"),
                (passive, "--", f"{passive} passive"),
            ):
                case = data.cases[case_id]
                if not scientific_response_eligible(data, case):
                    continue
                series = normalized_history_series(case, metric)
                if series is None:
                    continue
                times, values = series
                partial = completion_status(case) != "pass"
                display_label = (
                    f"{label} (partial to t={times[-1]:.2f})"
                    if partial else label
                )
                axis.plot(
                    times, values, style, color=CASE_COLORS[case_id],
                    label=display_label, linewidth=1.2,
                )
                if partial:
                    axis.scatter(
                        [times[-1]], [values[-1]], s=22, marker="o",
                        facecolors="none", edgecolors=CASE_COLORS[case_id],
                        linewidths=0.9, zorder=4,
                    )
                available += 1
            axis.axvspan(4.0, 10.0, color="#eeeeee", alpha=0.5, zorder=-10)
            axis.grid(True, alpha=0.25)
            axis.set_xlim(0.0, 10.0)
            axis.set_ylabel(ylabel)
            axis.set_title(f"{active} / {passive}: {ylabel}")
            if available:
                axis.legend(frameon=False, loc="best")
            else:
                axis.text(
                    0.5, 0.5, "history not yet available",
                    transform=axis.transAxes, ha="center", va="center", color="#666666",
                )
    for axis in axes[-1]:
        axis.set_xlabel(r"$t/(L_\perp/v_A)$")
    fig.suptitle(
        "Matched-design active/passive histories (health-eligible cases only)",
        y=0.995,
    )
    fig.text(
        0.5, 0.008,
        "Failed, incomplete, and numerically inconclusive cases are excluded. "
        "Interpret a pair only when both curves cover the shaded t=4--10 interval.",
        ha="center", va="bottom", fontsize=6.8,
    )
    fig.tight_layout(rect=(0.0, 0.025, 1.0, 0.985))
    save_figure(fig, path)
    plt.close(fig)


def active_passive_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return active/passive developed-window summary rows."""

    rows: list[dict[str, object]] = []
    for active, passive in ACTIVE_PASSIVE_PAIRS:
        for metric in ("kinetic", "magnetic", "abs_dp", "unstable", "nu_eff"):
            left = response_metric_value(data, data.cases[active], metric)
            right = response_metric_value(data, data.cases[passive], metric)
            rows.append({
                "pair": f"{active}/{passive}",
                "metric": metric,
                "active": left,
                "passive": right,
                "passive_minus_active": (
                    right - left if left is not None and right is not None else None
                ),
                "passive_over_active": (
                    right / left
                    if left is not None and right is not None and left != 0.0 else None
                ),
                "active_acceptance": acceptance_status(data, data.cases[active]),
                "passive_acceptance": acceptance_status(data, data.cases[passive]),
            })
    return rows


def robustness_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return the fixed robustness-contrast table."""

    rows: list[dict[str, object]] = []
    for label, reference, variant in ROBUSTNESS_CONTRASTS:
        for metric in ("kinetic", "magnetic", "abs_dp", "unstable", "nu_eff"):
            left = response_metric_value(data, data.cases[reference], metric)
            right = response_metric_value(data, data.cases[variant], metric)
            scale = (
                max(abs(left), abs(right), 1.0e-300)
                if left is not None and right is not None else None
            )
            rows.append({
                "contrast": label,
                "reference": reference,
                "variant": variant,
                "metric": metric,
                "reference_value": left,
                "variant_value": right,
                "variant_minus_reference": (
                    right - left if left is not None and right is not None else None
                ),
                "signed_relative_difference": (
                    (right - left) / scale
                    if left is not None and right is not None and scale is not None
                    else None
                ),
            })
    return rows


def render_robustness(data: PublicationData, plt: Any, colors: Any, path: Path) -> None:
    """Render a fixed signed-relative-effect robustness heatmap."""

    metrics = ("kinetic", "magnetic", "abs_dp", "unstable", "nu_eff")
    rows = robustness_rows(data)
    values: list[list[float]] = []
    for label, _, _ in ROBUSTNESS_CONTRASTS:
        by_metric = {
            str(row["metric"]): row["signed_relative_difference"]
            for row in rows if row["contrast"] == label
        }
        values.append([
            float(by_metric[metric])
            if isinstance(by_metric.get(metric), (int, float)) else math.nan
            for metric in metrics
        ])
    has_values = any(math.isfinite(value) for row in values for value in row)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(STATUS_COLORS["unknown"])
    fig, axis = plt.subplots(figsize=(8.4, 5.8))
    image = axis.imshow(values, cmap=cmap, vmin=-1.0, vmax=1.0, aspect="auto")
    axis.set_xticks(range(len(metrics)), metrics)
    axis.set_yticks(
        range(len(ROBUSTNESS_CONTRASTS)),
        [label for label, _, _ in ROBUSTNESS_CONTRASTS],
    )
    axis.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    for row, values_row in enumerate(values):
        for column, value in enumerate(values_row):
            label = "--" if not math.isfinite(value) else f"{value:+.2f}"
            axis.text(column, row, label, ha="center", va="center", fontsize=7)
    if has_values:
        colorbar = fig.colorbar(image, ax=axis, pad=0.02)
        colorbar.set_label(
            r"$(\mathrm{variant}-\mathrm{reference})/"
            r"\max(|\mathrm{variant}|,|\mathrm{reference}|)$"
        )
    else:
        axis.text(
            0.5, -0.07,
            "No paired complete developed-window contrasts are available yet.",
            transform=axis.transAxes, ha="center", va="top", color="#555555",
        )
    axis.set_title(
        (
            "Developed-window robustness summary "
            "(no paired contrasts available yet)"
            if not has_values else
            "Developed-window robustness summary "
            "(R10 excluded from strict-CGL trends)"
        ),
        pad=28,
    )
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def limiter_heat_flux_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return limiter and heat-flux scan summaries."""

    rows: list[dict[str, object]] = []
    for family, case_ids in (
        ("heat_flux", HEAT_FLUX_CASES),
        ("limiter", LIMITER_CASES),
    ):
        for case_id in case_ids:
            case = data.cases[case_id]
            rows.append({
                "family": family,
                "case_id": case_id,
                "completion": completion_status(case),
                "acceptance": acceptance_status(data, case),
                "claim_scope": (
                    "restricted nonfatal-hard-bound diagnostic"
                    if case_id in {"R14", "R15"}
                    and selected_nonfatal_hard_bound_variant(case)
                    else "restricted" if scope_status(case_id) == "restricted"
                    else "standard"
                ),
                "lf_k_parallel": case.model.get("lf_k_parallel"),
                "limiter_hardwall": case.model.get("limiter_hardwall"),
                "limiter_nu_coll": case.model.get("limiter_nu_coll"),
                "strict_admissibility": case.model.get("cgl_lf_strict_admissibility"),
                "abs_dp": response_metric_value(data, case, "abs_dp"),
                "unstable_fraction": response_metric_value(data, case, "unstable"),
                "nu_eff": response_metric_value(data, case, "nu_eff"),
                "hard_bound_fraction": response_metric_value(data, case, "hard_bound"),
                "applied_heat_flux_work_abs": response_lf_ledger_metric(
                    data, case, "heat_flux_work"
                ),
                "applied_pressure_work_abs": response_lf_ledger_metric(
                    data, case, "pressure_work"
                ),
            })
    return rows


def compact_number(value: object) -> str:
    """Return a compact numeric model label."""

    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return "--"
    return f"{parsed:.3g}" if math.isfinite(parsed) else "--"


def scan_tick_labels(
    data: PublicationData, case_ids: tuple[str, ...], family: str
) -> list[str]:
    """Return self-contained categorical labels for one parameter scan."""

    labels: list[str] = []
    for case_id in case_ids:
        case = data.cases[case_id]
        passive = str(case.model.get("passive_delta", "")).lower() in {"true", "1"}
        if family == "heat_flux":
            detail = rf"$k_\parallel$={compact_number(case.model.get('lf_k_parallel'))}"
        elif str(case.model.get("limiter_hardwall", "")).lower() in {"true", "1"}:
            detail = "hard wall"
        else:
            detail = (
                rf"$\nu_{{lim}}$={compact_number(case.model.get('limiter_nu_coll'))}"
            )
        qualifiers = [detail]
        if passive:
            qualifiers.append("passive")
        if str(case.model.get("cgl_lf_strict_admissibility", "")).lower() in {
            "false", "0"
        }:
            qualifiers.append("strict off")
        if selected_nonfatal_hard_bound_variant(case):
            qualifiers.append("diagnostic scope")
        labels.append("\n".join([case_id, *qualifiers]))
    return labels


def render_scan_row(
    axes: Iterable[Any],
    data: PublicationData,
    case_ids: tuple[str, ...],
    metrics: tuple[tuple[str, str], ...],
    tick_labels: list[str],
) -> None:
    """Render one categorical scan row."""

    for axis, (metric, label) in zip(axes, metrics):
        values: list[float | None] = []
        for case_id in case_ids:
            case = data.cases[case_id]
            if metric in {"heat_flux_work", "pressure_work"}:
                value = response_lf_ledger_metric(data, case, metric)
            else:
                value = response_metric_value(data, case, metric)
            values.append(value)
        finite = [value for value in values if value is not None]
        for index, (case_id, value) in enumerate(zip(case_ids, values)):
            if value is not None:
                axis.scatter(
                    index, value, s=34, color=CASE_COLORS[case_id],
                    edgecolor="black", linewidth=0.35, zorder=3,
                )
        if finite and min(finite) > 0.0 and max(finite) / min(finite) > 100.0:
            axis.set_yscale("log")
        if not finite:
            axis.text(
                0.5, 0.5, "not yet available",
                transform=axis.transAxes, ha="center", va="center", color="#666666",
            )
        axis.set_xticks(range(len(case_ids)), tick_labels)
        axis.set_ylabel(label)
        axis.grid(True, axis="y", alpha=0.25)


def render_limiter_heat_flux(data: PublicationData, plt: Any, path: Path) -> None:
    """Render limiter and heat-flux scan summaries."""

    fig, axes = plt.subplots(4, 2, figsize=(7.6, 10.2))
    heat_axes = tuple(axes[0]) + tuple(axes[1])
    limiter_axes = tuple(axes[2]) + tuple(axes[3])
    render_scan_row(
        heat_axes, data, HEAT_FLUX_CASES,
        (
            ("abs_dp", r"$\langle|\Delta p|\rangle$"),
            ("unstable", "unstable fraction"),
            ("heat_flux_work", r"$|W_q|$"),
            ("pressure_work", r"$|W_{\Delta p}|$"),
        ),
        scan_tick_labels(data, HEAT_FLUX_CASES, "heat_flux"),
    )
    render_scan_row(
        limiter_axes, data, LIMITER_CASES,
        (
            ("abs_dp", r"$\langle|\Delta p|\rangle$"),
            ("unstable", "unstable fraction"),
            ("nu_eff", r"$\langle\nu_{\rm eff}\rangle$"),
            ("hard_bound", "hard-bound fraction"),
        ),
        scan_tick_labels(data, LIMITER_CASES, "limiter"),
    )
    axes[0, 0].text(
        -0.12, 1.28, "Heat-flux scan", transform=axes[0, 0].transAxes,
        fontsize=10, fontweight="bold",
    )
    axes[2, 0].text(
        -0.12, 1.28, "Limiter / transport-coupled scan",
        transform=axes[2, 0].transAxes, fontsize=10, fontweight="bold",
    )
    heat_available = sum(
        response_metric_value(data, data.cases[case_id], "abs_dp") is not None
        for case_id in HEAT_FLUX_CASES
    )
    limiter_available = sum(
        response_metric_value(data, data.cases[case_id], "abs_dp") is not None
        for case_id in LIMITER_CASES
    )
    if heat_available < 2:
        axes[0, 1].text(
            1.0, 1.28, "fewer than two developed-window cases; no scan trend",
            transform=axes[0, 1].transAxes, ha="right", color="#666666",
        )
    if limiter_available < 2:
        axes[2, 1].text(
            1.0, 1.28, "fewer than two developed-window cases; no scan trend",
            transform=axes[2, 1].transAxes, ha="right", color="#666666",
        )
    fig.suptitle(
        "Developed-window limiter and Landau-fluid summaries", y=0.995
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.975), h_pad=2.5)
    save_figure(fig, path)
    plt.close(fig)


def resolution_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return resolution-case scalars and available convergence distances."""

    rows: list[dict[str, object]] = []
    for case_id in RESOLUTION_CASES:
        case = data.cases[case_id]
        for metric in ("kinetic", "magnetic", "abs_dp", "unstable", "nu_eff"):
            rows.append({
                "record_type": "case_metric",
                "case_or_product": case_id,
                "metric": metric,
                "value": response_metric_value(data, case, metric),
                "result": acceptance_status(data, case),
            })
    detail = (
        nested(data.comparisons, "resolution_detail.products")
        if all(
            scientific_response_eligible(data, data.cases[case_id])
            for case_id in RESOLUTION_CASES
        )
        else None
    )
    if isinstance(detail, dict):
        for product, record in sorted(detail.items()):
            if not isinstance(record, dict):
                continue
            rows.append({
                "record_type": "fast_report_distance",
                "case_or_product": product,
                "metric": "log_rms_R16_R02",
                "value": record.get("log_rms_R16_R02"),
                "result": None,
            })
            rows.append({
                "record_type": "fast_report_distance",
                "case_or_product": product,
                "metric": "log_rms_R02_R17",
                "value": record.get("log_rms_R02_R17"),
                "result": None,
            })
    gate = find_gate(data.campaign_acceptance, "R16_R02_R17_resolution_convergence")
    if isinstance(gate, dict) and isinstance(gate.get("observations"), list):
        for index, record in enumerate(gate["observations"]):
            if isinstance(record, dict):
                rows.append({
                    "record_type": "acceptance_convergence",
                    "case_or_product": record.get("product", record.get("metric", index)),
                    "metric": "campaign_gate",
                    "value": record.get(
                        "R02_R17_distance",
                        record.get("relative_R02_R17_difference"),
                    ),
                    "result": record.get("passed"),
                })
    for record in science_resolution_rows(data):
        rows.append({
            "record_type": "reviewed_science_convergence",
            "case_or_product": record.get("name"),
            "metric": record.get("kind"),
            "value": record.get("observations"),
            "result": record.get("passed")
            if record.get("available") is True else record.get("result"),
        })
    return rows


def find_gate(evidence: dict[str, Any] | None, name: str) -> dict[str, Any] | None:
    """Return one named acceptance gate."""

    if (
        isinstance(evidence, dict)
        and evidence.get("_publication_evidence_validated") is not True
    ):
        return None
    gates = evidence.get("gates") if isinstance(evidence, dict) else None
    if not isinstance(gates, list):
        return None
    for gate in gates:
        if isinstance(gate, dict) and gate.get("name") == name:
            return gate
    return None


def render_resolution(data: PublicationData, plt: Any, path: Path) -> None:
    """Render scalar resolution ratios and available convergence distances."""

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4))
    metrics = ("kinetic", "magnetic", "abs_dp", "unstable")
    width = 0.22
    reference = data.cases["R02"]
    any_ratio = False
    comparison_cases = ("R16", "R17")
    for case_index, case_id in enumerate(comparison_cases):
        ratios: list[float] = []
        positions: list[float] = []
        for metric_index, metric in enumerate(metrics):
            baseline = response_metric_value(data, reference, metric)
            value = response_metric_value(data, data.cases[case_id], metric)
            if baseline is None or value is None or baseline == 0.0:
                continue
            ratios.append(value / baseline)
            positions.append(metric_index + (case_index - 0.5) * width)
        if ratios:
            axes[0].bar(
                positions, ratios, width=width, label=case_id,
                color=CASE_COLORS[case_id], edgecolor="black", linewidth=0.35,
            )
            any_ratio = True
    axes[0].axhline(1.0, color="black", linewidth=0.8, linestyle=":")
    axes[0].set_xticks(range(len(metrics)), metrics)
    axes[0].set_ylabel("ratio to R02")
    axes[0].set_title("Developed-window scalar sensitivity")
    axes[0].grid(True, axis="y", alpha=0.25)
    if any_ratio:
        axes[0].legend(frameon=False)
    else:
        axes[0].text(
            0.5, 0.54,
            "no non-reference developed-window metrics yet",
            transform=axes[0].transAxes, ha="center", va="center", color="#666666",
        )
        axes[0].text(
            0.5, 0.43,
            "R02/R02 self-ratios are intentionally suppressed",
            transform=axes[0].transAxes, ha="center", va="center",
        )

    distance_records: list[tuple[str, float, float]] = []
    detail = (
        nested(data.comparisons, "resolution_detail.products")
        if all(
            scientific_response_eligible(data, data.cases[case_id])
            for case_id in RESOLUTION_CASES
        )
        else None
    )
    if isinstance(detail, dict):
        for product, record in sorted(detail.items()):
            if not isinstance(record, dict):
                continue
            low = as_float(record.get("log_rms_R16_R02"))
            high = as_float(record.get("log_rms_R02_R17"))
            if low is not None and high is not None:
                distance_records.append((str(product), low, high))
    if distance_records:
        positions = list(range(len(distance_records)))
        axes[1].bar(
            [value - 0.18 for value in positions],
            [record[1] for record in distance_records],
            width=0.36, label="R16/R02", color=CASE_COLORS["R16"],
        )
        axes[1].bar(
            [value + 0.18 for value in positions],
            [record[2] for record in distance_records],
            width=0.36, label="R02/R17", color=CASE_COLORS["R17"],
        )
        axes[1].set_xticks(
            positions, [record[0] for record in distance_records], rotation=25
        )
        axes[1].set_ylabel("log-RMS curve distance")
        axes[1].legend(frameon=False)
    else:
        gate = (
            find_gate(data.science_record, "R16_R02_R17_resolution_convergence")
            or find_gate(
                data.campaign_acceptance, "R16_R02_R17_resolution_convergence"
            )
        )
        gate_result = gate.get("result") if isinstance(gate, dict) else "not available"
        axes[1].text(
            0.5, 0.56, "curve-distance evidence not yet available",
            transform=axes[1].transAxes, ha="center", va="center", color="#666666",
        )
        axes[1].text(
            0.5, 0.43, f"reviewed convergence gate: {gate_result}",
            transform=axes[1].transAxes, ha="center", va="center",
        )
    axes[1].set_title("Common-scale convergence evidence")
    axes[1].grid(True, axis="y", alpha=0.25)
    fig.suptitle("R16 / R02 / R17 resolution summary", y=0.995)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def case_evidence_roots(data: PublicationData, case_id: str) -> list[dict[str, Any]]:
    """Return all case-specific diagnostic, acceptance, and audit roots."""

    case = data.cases[case_id]
    roots: list[dict[str, Any]] = []
    for record in data.audit_records:
        case_ids = record.get("_publication_case_ids")
        if record.get("case_id") == case_id or (
            isinstance(case_ids, list) and case_id in case_ids
        ):
            roots.append(record)
        case_record = nested(record, f"cases.{case_id}")
        if isinstance(case_record, dict):
            roots.append(case_record)
    roots.extend(
        root for root in (case.acceptance, case.diagnostics) if isinstance(root, dict)
    )
    for record in data.acceptance_records:
        if record.get("_publication_evidence_validated") is not True:
            continue
        if record.get("case_id") == case_id:
            roots.append(record)
        case_record = nested(record, f"cases.{case_id}")
        if isinstance(case_record, dict):
            roots.append(case_record)
    return roots


def first_evidence_value(
    roots: Iterable[dict[str, Any]], paths: Iterable[str]
) -> object | None:
    """Return the first available value among ordered evidence paths."""

    for root in roots:
        for path in paths:
            value = nested(root, path)
            if value is not None:
                return value
    return None


def hyperbolicity_diagnostics(
    data: PublicationData, case_id: str
) -> dict[str, object]:
    """Return normalized retained-state hyperbolicity diagnostics when available."""

    roots = case_evidence_roots(data, case_id)
    snapshot_records: list[dict[str, object]] = []
    for root in roots:
        snapshots = root.get("snapshots")
        if not isinstance(snapshots, list):
            continue
        evaluated = 0.0
        negative = 0.0
        nonfinite = 0.0
        minimum: float | None = None
        for snapshot in snapshots:
            aggregate = snapshot.get("aggregate") if isinstance(snapshot, dict) else None
            if not isinstance(aggregate, dict):
                continue
            evaluated_value = as_float(aggregate.get("evaluated"))
            negative_value = as_float(aggregate.get("negative"))
            nonfinite_value = as_float(aggregate.get("nonfinite_discriminant"))
            minimum_value = as_float(aggregate.get("minimum"))
            if evaluated_value is not None:
                evaluated += evaluated_value
            if negative_value is not None:
                negative += negative_value
            if nonfinite_value is not None:
                nonfinite += nonfinite_value
            if minimum_value is not None and (
                minimum is None or minimum_value < minimum
            ):
                minimum = minimum_value
        if evaluated > 0.0:
            snapshot_records.append({
                "result": (
                    "fail" if negative > 0.0 or nonfinite > 0.0 or (
                        minimum is not None and minimum < 0.0
                    ) else "pass"
                ),
                "negative_discriminant_fraction": negative / evaluated,
                "negative_discriminant_count": negative,
                "cell_direction_evaluations": evaluated,
                "minimum_discriminant": minimum,
                "nonfinite_discriminant_count": nonfinite,
            })
    if snapshot_records:
        adverse = [
            record for record in snapshot_records if record["result"] == "fail"
        ]
        return max(
            adverse or snapshot_records,
            key=lambda record: (
                as_float(record["negative_discriminant_fraction"]) or 0.0,
                as_float(record["cell_direction_evaluations"]) or 0.0,
            ),
        )
    status = first_evidence_value(
        roots,
        (
            "hyperbolicity.result",
            "retained_state_hyperbolicity.result",
        ),
    )
    normalized_status = str(status) if status in {"pass", "fail", "inconclusive"} else None
    negative_fraction = as_float(first_evidence_value(
        roots,
        (
            "hyperbolicity.negative_discriminant_fraction",
            "hyperbolicity.negative_fraction",
            "retained_state_hyperbolicity.negative_discriminant_fraction",
            "retained_state_hyperbolicity.negative_fraction",
            "summary.negative_discriminant_fraction",
            "summary.negative_fraction",
            "negative_discriminant_fraction",
            "negative_fraction",
        ),
    ))
    negative_count = as_float(first_evidence_value(
        roots,
        (
            "hyperbolicity.negative_discriminant_count",
            "hyperbolicity.negative_count",
            "retained_state_hyperbolicity.negative_discriminant_count",
            "retained_state_hyperbolicity.negative_count",
            "summary.negative_discriminant_count",
            "summary.negative_count",
            "negative_discriminant_count",
            "negative_count",
        ),
    ))
    evaluation_count = as_float(first_evidence_value(
        roots,
        (
            "hyperbolicity.cell_direction_evaluations",
            "hyperbolicity.evaluation_count",
            "retained_state_hyperbolicity.cell_direction_evaluations",
            "retained_state_hyperbolicity.evaluation_count",
            "summary.cell_direction_evaluations",
            "summary.evaluation_count",
            "cell_direction_evaluations",
            "evaluation_count",
        ),
    ))
    minimum = as_float(first_evidence_value(
        roots,
        (
            "hyperbolicity.minimum_discriminant",
            "hyperbolicity.minimum",
            "retained_state_hyperbolicity.minimum_discriminant",
            "retained_state_hyperbolicity.minimum",
            "summary.minimum_discriminant",
            "summary.minimum",
            "minimum_discriminant",
            "minimum",
        ),
    ))
    if normalized_status is None:
        if (
            (negative_fraction is not None and negative_fraction > 0.0)
            or (negative_count is not None and negative_count > 0.0)
            or (minimum is not None and minimum < 0.0)
        ):
            normalized_status = "fail"
        elif negative_count == 0.0 and evaluation_count is not None and evaluation_count > 0.0:
            normalized_status = "pass"
        else:
            normalized_status = "unknown"
    return {
        "result": normalized_status,
        "negative_discriminant_fraction": negative_fraction,
        "negative_discriminant_count": negative_count,
        "cell_direction_evaluations": evaluation_count,
        "minimum_discriminant": minimum,
    }


def hyperbolicity_status(data: PublicationData, case_id: str) -> str:
    """Return the normalized retained-state hyperbolicity disposition."""

    return str(hyperbolicity_diagnostics(data, case_id)["result"])


def scope_observed_diagnostic(
    data: PublicationData, case_id: str
) -> tuple[str, str]:
    """Return one compact observed scope diagnostic and its status color."""

    case = data.cases[case_id]
    hard_bound = as_float(nested(case.diagnostics, "health.hard_bound_diagnostic_maximum"))
    if hard_bound is None:
        hard_bound = as_float(nested(case.direct_acceptance, "health.hard_bound_maximum"))
    if case_id in {"R14", "R15"} and hard_bound is not None:
        status = "warning" if hard_bound > 0.0 else "pass"
        return f"hard-bound={hard_bound:.3g}", status
    diagnostics = hyperbolicity_diagnostics(data, case_id)
    fraction = as_float(diagnostics["negative_discriminant_fraction"])
    minimum = as_float(diagnostics["minimum_discriminant"])
    if fraction is not None:
        status = "fail" if fraction > 0.0 else "pass"
        return f"negative D={100.0 * fraction:.3g}%", status
    if minimum is not None and minimum < 0.0:
        return f"min D={minimum:.3g}", "fail"
    if hard_bound is not None:
        status = "warning" if hard_bound > 0.0 else "pass"
        return f"hard-bound={hard_bound:.3g}", status
    return "--", "unknown"


def selected_nonfatal_hard_bound_variant(case: CaseRecord) -> bool:
    """Return whether the selected lineage declares a nonfatal hard-bound variant."""

    variants = nested(case.lineage, "lineage_variants")
    return isinstance(variants, list) and any(
        isinstance(value, str)
        and "nonfatal" in value.lower()
        and "hard_bound" in value.lower()
        for value in variants
    )


def scope_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return explicit R10/R14/R15 claim-scope disclosures."""

    rows: list[dict[str, object]] = []
    for case_id in ("R10", "R14", "R15"):
        case = data.cases[case_id]
        health = nested(case.diagnostics, "health")
        if not isinstance(health, dict):
            health = {}
        hyper = hyperbolicity_diagnostics(data, case_id)
        observed, _ = scope_observed_diagnostic(data, case_id)
        if case_id == "R10":
            scope = (
                "Exploratory stress test only; exclude from strict-hyperbolic CGL "
                "and beta-trend claims."
            )
        elif case_id == "R14":
            scope = (
                "Finite-limiter hard-bound diagnostic; strict admissibility is "
                "nonfatal and the limiter scan also changes local LF transport."
            )
        else:
            scope = (
                "Any selected nonfatal-hard-bound R15 lineage is a restricted "
                "diagnostic variant only, never a strict-admissibility success. "
                "Strict R15 failure details require authenticated evidence."
            )
        variants = nested(case.lineage, "lineage_variants")
        variant = (
            "; ".join(value for value in variants if isinstance(value, str))
            if isinstance(variants, list) and any(
                isinstance(value, str) for value in variants
            )
            else None
        )
        if case_id == "R14" and not selected_nonfatal_hard_bound_variant(case):
            scope = (
                "R14 is restricted, but its intended nonfatal-hard-bound variant "
                "cannot be verified from the selected lineage variant evidence."
            )
        if case_id == "R15" and not selected_nonfatal_hard_bound_variant(case):
            scope += (
                " The independent nonfatal-hard-bound diagnostic variant is not yet "
                "verified as the selected lineage."
            )
        hard_bound_maximum = as_float(
            nested(case.diagnostics, "health.hard_bound_diagnostic_maximum")
        )
        if hard_bound_maximum is None:
            hard_bound_maximum = as_float(
                nested(case.direct_acceptance, "health.hard_bound_maximum")
            )
        rows.append({
            "case_id": case_id,
            "variant": variant,
            "completion": completion_status(case),
            "fast_health": fast_health_status(case),
            "acceptance": acceptance_status(data, case),
            "strict_admissibility": case.model.get("cgl_lf_strict_admissibility"),
            "strict_run_disposition": (
                r15_strict_failure_disposition(data) if case_id == "R15" else None
            ),
            "retained_state_hyperbolicity": hyper["result"],
            "negative_discriminant_fraction": hyper["negative_discriminant_fraction"],
            "negative_discriminant_count": hyper["negative_discriminant_count"],
            "cell_direction_evaluations": hyper["cell_direction_evaluations"],
            "minimum_discriminant": hyper["minimum_discriminant"],
            "hard_bound_maximum": hard_bound_maximum,
            "hard_bound_volume_maximum": health.get("hard_bound_volume_maximum"),
            "observed_diagnostic": observed,
            "claim_scope": scope,
        })
    return rows


def render_scope(
    data: PublicationData, plt: Any, colors: Any, patches: Any, path: Path
) -> None:
    """Render explicit R10/R14/R15 scope restrictions."""

    statuses: list[list[str]] = []
    labels: list[list[str]] = []
    for case_id in ("R10", "R14", "R15"):
        case = data.cases[case_id]
        strict = str(case.model.get("cgl_lf_strict_admissibility", "")).lower()
        r15_strict = r15_strict_failure_disposition(data)
        strict_status = (
            "fail" if case_id == "R15" and r15_strict != "unavailable/inconclusive"
            else "inconclusive" if case_id == "R15"
            else (
            "configuration" if strict in {"true", "1", "false", "0"} else "unknown"
            )
        )
        hyper = hyperbolicity_status(data, case_id)
        observed, observed_status = scope_observed_diagnostic(data, case_id)
        statuses.append([
            completion_status(case),
            acceptance_status(data, case),
            strict_status,
            hyper,
            observed_status,
            "restricted",
        ])
        labels.append([
            display_status(completion_status(case)),
            display_status(acceptance_status(data, case)),
            display_status(strict_status) if case_id == "R15" else (
                "enabled" if strict in {"true", "1"} else (
                "disabled" if strict in {"false", "0"} else "--"
                )
            ),
            display_status(hyper),
            observed,
            "restricted",
        ])
    status_matrix_figure(
        plt, colors, patches, ["R10", "R14", "R15"],
        [
            "Completion", "Acceptance", "Strict policy",
            "Retained-state\nhyperbolicity", "Observed\ndiagnostic", "Claim scope",
        ],
        statuses, labels, "Explicit scope restrictions for R10, R14, and R15", path,
        (
            "R10: exploratory stress test; exclude from strict-hyperbolic and beta-trend "
            "claims. R14: nonfatal hard-bound diagnostic; disclose coupled limiter/LF "
            "transport semantics. Any selected R15 nonfatal variant is diagnostic "
            "only; strict-R15 failure details require authenticated evidence. "
            "Strict-policy cells report configuration only except the explicit R15 fail."
        ),
    )


def render_science_ct_summary(
    data: PublicationData, plt: Any, colors: Any, patches: Any, path: Path
) -> None:
    """Render reviewed science gates and direct CT numerical health."""

    rows: list[tuple[str, str]] = [
        ("reviewed science aggregate", aggregate_science_status(data))
    ]
    rows.extend(
        (
            str(row["gate"]) + (
                " [diagnostic scope]"
                if str(row.get("claim_scope", "")).startswith("restricted")
                else ""
            ),
            str(row["result"]),
        )
        for row in science_gate_rows(data)
    )
    resolution = nested(data.science_record, "resolution.result")
    mks24 = nested(data.science_record, "mks24.result")
    rows.extend([
        (
            "resolution aggregate",
            validated_result(resolution) if resolution is not None else "unknown",
        ),
        (
            "MKS24 residual/drift aggregate",
            validated_result(mks24) if mks24 is not None else "unknown",
        ),
        ("direct CT numerical aggregate", aggregate_ct_status(data)),
    ])
    status_matrix_figure(
        plt, colors, patches,
        [name for name, _ in rows],
        ["Result", "Release authority"],
        [[result, "configuration"] for _, result in rows],
        [[display_status(result), "non-authorizing"] for _, result in rows],
        "Authenticated reviewed-science gates and direct CT health",
        path,
        (
            "All displayed results are authenticated direct-fast evidence and remain "
            "explicitly non-authorizing. Inconclusive or unavailable partial evidence "
            "is never promoted to pass. A direct CT pass is a sampled numerical result, "
            "not campaign or release authority."
        ),
    )


def acceptance_gate_rows(data: PublicationData) -> list[dict[str, object]]:
    """Flatten all discovered scientific-acceptance gates."""

    rows: list[dict[str, object]] = []
    for record in data.acceptance_records:
        record_type = record.get("record_type")
        case_id = record.get("case_id", "campaign")
        if record.get("_publication_evidence_validated") is not True:
            rows.append({
                "record_type": record_type,
                "case_id": case_id,
                "gate": "evidence_authentication",
                "result": "fail",
                "reason": "record failed publication-renderer provenance or freshness validation",
            })
            continue
        gates = record.get("gates")
        if not isinstance(gates, list):
            rows.append({
                "record_type": record_type,
                "case_id": case_id,
                "gate": "--",
                "result": record.get("result", record.get("valid")),
                "reason": record.get("non_authorizing_statement"),
            })
            continue
        for gate in gates:
            if isinstance(gate, dict):
                rows.append({
                    "record_type": record_type,
                    "case_id": case_id,
                    "gate": gate.get("name"),
                    "result": gate.get("result"),
                    "reason": gate.get("reason"),
                })
    return rows


def source_binding(path: Path) -> dict[str, object]:
    """Return a deterministic input or product binding."""

    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def source_identity_summary(data: PublicationData) -> dict[str, list[str]]:
    """Return selected-lineage executable, matrix, input, and launcher identities."""

    identities: dict[str, set[str]] = {
        "executable_sha256": set(),
        "matrix_sha256": set(),
        "input_sha256": set(),
        "fast_script_sha256": set(),
    }
    for case in data.cases.values():
        for name in tuple(identities):
            values = nested(case.lineage, f"lineage_identities.{name}")
            if isinstance(values, list):
                identities[name].update(
                    value for value in values if isinstance(value, str)
                )
            lineage = nested(case.lineage, "lineage")
            if isinstance(lineage, list):
                identities[name].update(
                    str(segment[name])
                    for segment in lineage
                    if isinstance(segment, dict) and isinstance(segment.get(name), str)
                )
    return {
        name: sorted(values)
        for name, values in identities.items()
        if values
    }


def normalized_invocation(
    analysis: Path, output: Path, acceptance_paths: Iterable[Path]
) -> list[str]:
    """Return a deterministic normalized renderer invocation."""

    invocation = [
        sys.executable,
        str(RENDERER_PATH),
        str(analysis),
        "--output",
        str(output),
    ]
    for path in sorted(path.absolute() for path in acceptance_paths):
        invocation.extend(["--acceptance", str(path)])
    return invocation


def report_markdown(data: PublicationData, products: list[Path], output: Path) -> str:
    """Return a compact renderer report."""

    complete = [
        case_id for case_id in CASE_IDS
        if completion_status(data.cases[case_id]) == "pass"
    ]
    analyzed = [
        case_id for case_id in CASE_IDS
        if data.cases[case_id].diagnostics is not None
    ]
    accepted = [
        case_id for case_id in CASE_IDS
        if acceptance_status(data, data.cases[case_id]) == "pass"
    ]
    case_warnings = case_warning_rows(data)
    paired_contrasts = sum(
        row["signed_relative_difference"] is not None
        for row in robustness_rows(data)
    )
    nonreference_resolution_cases = [
        case_id for case_id in ("R16", "R17")
        if any(
            response_metric_value(data, data.cases[case_id], metric) is not None
            for metric in ("kinetic", "magnetic", "abs_dp", "unstable")
        )
    ]
    science_cases = nested(data.science_record, "selected_cases")
    ct_cases = nested(data.ct_audit_record, "selection.cases")
    lines = [
        "# Stage I Publication Products",
        "",
        "This downstream renderer is descriptive and non-authorizing. Missing evidence "
        "remains visible in every product.",
        "",
        f"- Evidence state: **{publication_evidence_state(data)}**",
        f"- Complete cases: `{', '.join(complete) or 'none'}`",
        f"- Fast-report diagnostics available: `{', '.join(analyzed) or 'none'}`",
        f"- Scientific-acceptance passes: `{', '.join(accepted) or 'none'}`",
        f"- Acceptance records discovered: `{len(data.acceptance_records)}`",
        f"- Audit records discovered: `{len(integrated_audit_records(data))}`",
        f"- Reviewed-science aggregate: **{aggregate_science_status(data)}** "
        f"(`{SCIENCE_AUTHORITY}`, release_authorizing=false)",
        f"- Reviewed-science selected cases: "
        f"`{', '.join(science_cases) if isinstance(science_cases, list) else 'none'}`",
        f"- Direct CT numerical aggregate: **{aggregate_ct_status(data)}** "
        "(campaign_authority_eligible=false, release_authorizing=false)",
        f"- Direct CT selected cases: "
        f"`{', '.join(ct_cases) if isinstance(ct_cases, list) else 'none'}`; "
        f"full Stage I coverage={text_value(ct_campaign_coverage_complete(data))}",
        "",
        "## Products",
        "",
    ]
    lines.extend(
        f"- `{path.relative_to(output)}`" for path in sorted(products)
    )
    lines.extend(["", "## Renderer Ingestion Warnings", ""])
    lines.extend(
        f"- {warning}" for warning in data.ingestion_warnings
    )
    if not data.ingestion_warnings:
        lines.append("- none")
    lines.extend(["", "## Case Lineage, Health, And Scientific Warnings", ""])
    lines.extend(
        f"- `{row['case_id']}` `{row['category']}`: {row['warning']}"
        for row in case_warnings
    )
    if not case_warnings:
        lines.append("- none")
    lines.extend([
        "",
        "## Partial-Product Limitations",
        "",
        f"- Populated paired developed-window robustness cells: `{paired_contrasts}`.",
        "- Active/passive history curves marked as partial are transient-only and "
        "must not be interpreted as developed-window comparisons.",
        "- Resolution self-ratios for R02/R02 are suppressed. Non-reference cases "
        f"with developed-window scalars: `{', '.join(nonreference_resolution_cases) or 'none'}`.",
        "- Parameter-scan trends require at least two populated developed-window cases.",
    ])
    lines.extend([
        "",
        "## Scope",
        "",
        "- R10 remains exploratory and is excluded from strict-hyperbolic and "
        "beta-trend claims.",
        "- R14 remains a disclosed nonfatal hard-bound diagnostic and not a "
        "uniform strict-admissibility case.",
        f"- R15 strict-run disposition is `{r15_strict_failure_disposition(data)}`. "
        "Any selected "
        "nonfatal-hard-bound R15 variant is a restricted diagnostic and never "
        "strict-admissibility success.",
        "- Standard claim scope and strict-policy enabled/disabled cells describe "
        "configuration only; they are not scientific passes.",
        "",
    ])
    return "\n".join(lines)


def render_products(data: PublicationData, output: Path) -> list[Path]:
    """Write all fixed publication figures, tables, and report products."""

    plt, colors, patches = configure_plotting()
    figures = output / "figures"
    tables = output / "tables"
    products: list[Path] = []

    figure_paths = {
        "health": figures / "fig01_health_completion.pdf",
        "active_passive": figures / "fig02_active_passive_history.pdf",
        "robustness": figures / "fig03_robustness_summary.pdf",
        "limiter_heat_flux": figures / "fig04_limiter_heat_flux_summary.pdf",
        "resolution": figures / "fig05_resolution_summary.pdf",
        "scope": figures / "fig06_r10_r14_r15_scope.pdf",
        "science_ct": figures / "fig07_reviewed_science_ct_summary.pdf",
    }
    render_health(data, plt, colors, patches, figure_paths["health"])
    render_active_passive(data, plt, figure_paths["active_passive"])
    render_robustness(data, plt, colors, figure_paths["robustness"])
    render_limiter_heat_flux(data, plt, figure_paths["limiter_heat_flux"])
    render_resolution(data, plt, figure_paths["resolution"])
    render_scope(data, plt, colors, patches, figure_paths["scope"])
    render_science_ct_summary(
        data, plt, colors, patches, figure_paths["science_ct"]
    )
    products.extend(figure_paths.values())

    table_specs = (
        (
            "health_completion",
            [
                "case_id", "assembly_status", "final_time", "target_time",
                "completion", "fast_health", "fatal_lf_counters", "acceptance",
                "reviewed_science", "direct_ct_numerical",
                "direct_ct_release_authorizing", "claim_scope",
                "strict_admissibility", "hard_bound_maximum",
                "hard_bound_volume_maximum", "structural_error_count",
                "structural_warning_count", "numerical_warning_count",
                "science_warning_count",
            ],
            health_rows(data),
        ),
        (
            "active_passive_summary",
            [
                "pair", "metric", "active", "passive", "passive_minus_active",
                "passive_over_active", "active_acceptance", "passive_acceptance",
            ],
            active_passive_rows(data),
        ),
        (
            "robustness_summary",
            [
                "contrast", "reference", "variant", "metric", "reference_value",
                "variant_value", "variant_minus_reference",
                "signed_relative_difference",
            ],
            robustness_rows(data),
        ),
        (
            "limiter_heat_flux_summary",
            [
                "family", "case_id", "completion", "acceptance", "lf_k_parallel",
                "claim_scope", "limiter_hardwall", "limiter_nu_coll",
                "strict_admissibility",
                "abs_dp", "unstable_fraction", "nu_eff", "hard_bound_fraction",
                "applied_heat_flux_work_abs", "applied_pressure_work_abs",
            ],
            limiter_heat_flux_rows(data),
        ),
        (
            "resolution_summary",
            ["record_type", "case_or_product", "metric", "value", "result"],
            resolution_rows(data),
        ),
        (
            "r10_r14_r15_scope",
            [
                "case_id", "variant", "completion", "fast_health", "acceptance",
                "strict_admissibility", "strict_run_disposition",
                "retained_state_hyperbolicity",
                "negative_discriminant_fraction", "negative_discriminant_count",
                "cell_direction_evaluations", "minimum_discriminant",
                "hard_bound_maximum", "hard_bound_volume_maximum",
                "observed_diagnostic", "claim_scope",
            ],
            scope_rows(data),
        ),
        (
            "case_health_scientific_warnings",
            ["case_id", "category", "warning"],
            case_warning_rows(data),
        ),
        (
            "acceptance_gates",
            ["record_type", "case_id", "gate", "result", "reason"],
            acceptance_gate_rows(data),
        ),
        (
            "reviewed_science_contrasts",
            [
                "family", "contrast", "left", "right", "result",
                "claim_eligible", "metric", "available", "left_mean",
                "right_mean", "difference", "combined_standard_error", "z_score",
                "two_sided_p", "standardized_effect", "holm_threshold",
                "holm_significant", "reason", "authority", "release_authorizing",
                "claim_scope",
            ],
            science_contrast_rows(data),
        ),
        (
            "reviewed_science_gates",
            [
                "gate", "result", "reason", "claim_scope", "observations", "limits",
                "authority", "release_authorizing",
            ],
            science_gate_rows(data),
        ),
        (
            "reviewed_science_resolution",
            [
                "result", "kind", "name", "available", "passed", "observations",
                "limits", "reason", "authority", "release_authorizing",
            ],
            science_resolution_rows(data),
        ),
        (
            "reviewed_science_mks24",
            [
                "panel", "panel_result", "product_id", "case_id", "source",
                "result", "normalized_residual_rms",
                "maximum_absolute_normalized_residual",
                "early_late_vector_drift_rms", "limits", "reason", "authority",
                "release_authorizing",
            ],
            science_mks24_rows(data),
        ),
        (
            "direct_ct_health",
            [
                "case_id", "ct_result", "ct_evidence_available",
                "ct_claim_supported", "coverage_complete",
                "maximum_normalized_ct_divb", "normalized_ct_divb_lt",
                "campaign_authority_eligible", "release_authorizing", "reason",
            ],
            ct_health_rows(data),
        ),
    )
    for name, columns, rows in table_specs:
        products.extend(write_table(tables, name, columns, rows))

    captions = f"""# Suggested Figure Captions

**Evidence state: {publication_evidence_state(data)}.** Scientific-response products
exclude failed, incomplete, and numerically inconclusive cases.

1. **Health and completion.** Direct-fast completion, numerical-health diagnostics,
   non-authorizing scientific-acceptance disposition, and claim scope for R02--R17.
   Green denotes an evidence-backed pass or clean result. Blue denotes incomplete
   integration; light blue-gray denotes configuration or policy, not eligibility;
   gray denotes missing or inconclusive evidence.
2. **Active/passive histories.** Matched-design active and passive CGL-pressure
   histories for Alfvenic and random forcing at beta 10 and 100. Failed, incomplete,
   and numerically inconclusive cases are excluded. The shaded interval is the
   developed-turbulence analysis window.
3. **Robustness.** Signed developed-window response to forcing geometry, beta, and
   forcing correlation time where both cases have complete-window evidence. Missing
   cells are not zero response. R10 is excluded from strict-CGL trend inference.
4. **Limiter and heat flux.** Developed-window summaries for the Landau-fluid
   strength scan and the limiter/transport-coupled scan. Configuration is printed
   below each case. Fewer than two populated cases do not establish a scan trend.
5. **Resolution.** R16/R02/R17 scalar sensitivity and common-scale convergence
   evidence when available. R02/R02 self-ratios are suppressed.
6. **Restricted scope.** Claim restrictions and available observed diagnostics for
   R10, R14, and R15. Any selected nonfatal-hard-bound R15 variant is diagnostic
   only; strict-R15 failure details appear only from authenticated evidence.
7. **Reviewed science and direct CT.** Authenticated Holm-corrected comparison
   gates, common-range convergence, admitted MKS24 residual/drift criteria, and
   sampled direct CT numerical health. All are explicitly non-authorizing.
"""
    captions_path = output / "captions.md"
    write_text(captions_path, captions)
    products.append(captions_path)

    report_path = output / "report.md"
    write_text(report_path, report_markdown(data, products, output))
    products.append(report_path)
    return products


def parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument(
        "analysis_output", type=Path,
        help="direct-fast report output containing cases/ and optional campaign/",
    )
    command.add_argument(
        "--output", type=Path,
        help=(
            "publication product directory "
            "(default: ANALYSIS_OUTPUT/publication-products)"
        ),
    )
    command.add_argument(
        "--acceptance", type=Path, action="append", default=[],
        help="additional scientific-acceptance JSON file or directory; repeatable",
    )
    return command


def main(argv: list[str] | None = None) -> int:
    """Render deterministic partial- or complete-campaign publication products."""

    args = parser().parse_args(argv)
    analysis = args.analysis_output.absolute()
    if not analysis.is_dir():
        raise SystemExit(f"analysis output does not exist: {analysis}")
    output = (
        args.output.absolute()
        if args.output is not None
        else analysis / "publication-products"
    )
    data = discover_data(analysis, args.acceptance)
    products = render_products(data, output)
    import matplotlib

    manifest = {
        "schema_version": 2,
        "record_type": "cgl_lf_stage_i_fast_publication_products",
        "evidence_state": publication_evidence_state(data),
        "analysis_output": str(analysis),
        "renderer": source_binding(RENDERER_PATH),
        "normalized_invocation": normalized_invocation(
            analysis, output, args.acceptance
        ),
        "runtime": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "matplotlib_version": matplotlib.__version__,
        },
        "determinism_contract": (
            "No wall-clock fields; sorted inputs and rows; fixed figure definitions; "
            "PDF creation/modification dates suppressed."
        ),
        "source_identity_summary": source_identity_summary(data),
        "case_ids": list(CASE_IDS),
        "acceptance_record_types": sorted({
            str(record.get("record_type")) for record in data.acceptance_records
        }),
        "audit_record_types": sorted({
            str(record.get("record_type"))
            for record in integrated_audit_records(data)
        }),
        "reviewed_science": {
            "record_type": data.science_record.get("record_type"),
            "result": aggregate_science_status(data),
            "selected_cases": data.science_record.get("selected_cases"),
            "authority": data.science_record.get("authority"),
            "release_authorizing": False,
        } if isinstance(data.science_record, dict) else None,
        "direct_ct_audit": {
            "record_type": data.ct_audit_record.get("record_type"),
            "numerical_result": aggregate_ct_status(data),
            "selected_cases": nested(data.ct_audit_record, "selection.cases"),
            "full_stage_i_coverage": ct_campaign_coverage_complete(data),
            "campaign_authority_eligible": False,
            "release_authorizing": False,
        } if isinstance(data.ct_audit_record, dict) else None,
        "sources": [
            source_binding(path) for path in sorted(data.source_paths)
            if path.is_file()
        ],
        "products": [
            source_binding(path) for path in sorted(products)
            if path.is_file()
        ],
        "renderer_ingestion_warnings": data.ingestion_warnings,
        "case_warning_records": case_warning_rows(data),
    }
    write_json(output / "manifest.json", manifest)
    print(
        f"rendered {len(products)} publication products in {output}; "
        f"acceptance_records={len(data.acceptance_records)}, "
        f"audit_records={len(integrated_audit_records(data))}, "
        f"ingestion_warnings={len(data.ingestion_warnings)}, "
        f"case_warnings={len(case_warning_rows(data))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
