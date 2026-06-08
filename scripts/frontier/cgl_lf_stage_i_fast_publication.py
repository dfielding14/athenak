#!/usr/bin/env python3
"""Render deterministic publication summaries for direct-fast Stage I evidence.

This renderer is intentionally downstream-only.  It reads the direct-fast
report output, discovers optional scientific-acceptance evidence, and writes a
fixed set of summary figures and tables.  Missing or partial campaign evidence
is shown explicitly rather than treated as a failed result.  Products are
rendered in a sibling staging directory and promoted before ``manifest.json``;
the canonical manifest is therefore the last-published authority.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass, field
import fcntl
import hashlib
import io
import json
import math
import os
import platform
from pathlib import Path
import re
import secrets
import stat
import sys
import tempfile
import textwrap
from typing import Any, Iterable, Iterator


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
ACTIVE_ENERGY_CASES = (
    "R02", "R03", "R04", "R05", "R10", "R11",
    "R12", "R13", "R14", "R15", "R16", "R17",
)
PRIMARY_SCALAR_METRICS = (
    "kinetic",
    "magnetic",
    "abs_dp",
    "beta",
    "mirror_occupancy",
    "firehose_occupancy",
    "nu_eff",
)
FATAL_LF_COUNTERS = ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos")
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
HYPERBOLICITY_JOB_RECORD_TYPE = "cgl_lf_stage_i_direct_fast_hyperbolicity_job"
HYPERBOLICITY_FORMULA_IDS = frozenset(
    ("qualified-legacy", "literature-correct")
)
HYPERBOLICITY_FORMULA_FAMILIES = {
    "qualified-legacy": "legacy_implementation",
    "literature-correct": "literature_correct",
}
PUBLICATION_AUTHORITY_TOKEN_NAME = ".cgl-lf-publication-authority-token"
PUBLICATION_MANIFEST_QUARANTINE_PREFIX = (
    ".cgl-lf-publication-manifest-quarantine-"
)
PUBLICATION_OWNERSHIP_RECORD_TYPE = (
    "cgl_lf_stage_i_fast_publication_output_ownership"
)
ACTIVE_PASSIVE_INTERVENTION_SCOPE = {
    "estimand": "total_effect_of_enabling_active_cgl",
    "enabled_components": [
        "pressure_feedback",
        "thermodynamic_evolution",
        "characteristic_speeds_and_fluxes",
        "realized_forcing_after_trajectory_divergence",
    ],
    "excluded_interpretation": "anisotropic_stress_alone",
    "declaration": (
        "The active/passive comparison estimates the total effect of enabling active-CGL "
        "pressure feedback, thermodynamic evolution, characteristic speeds and fluxes, "
        "and realized forcing after trajectory divergence; it is not an "
        "anisotropic-stress-only comparison."
    ),
    "claim_scope": "descriptive_within_realization",
}
CURRENT_SCIENCE_SCOPE_LIMITATION = {
    "binding_semantics": (
        "Exact current-byte bindings identify the implementation that carries both "
        "historically independently approved scope and current reviewed production-"
        "science corrections; they do not expand any approval beyond its explicit scope."
    ),
    "current_reviewed_production_science_corrections": [
        "Aggressive family-gate revisions, including descriptive active/passive, LF-strength, finite-limiter, and forcing decisions.",
        "Exactly-once reduction of raw intensive history totals to volume means or fractions.",
        "The robust inertial-range alignment scalar and its descriptive within-realization acceptance statistics.",
    ],
    "disposition": (
        "accepted_reviewed_production_science_corrections_without_new_independent_"
        "plasma_or_statistical_approval"
    ),
    "historical_independent_approval_scope": [
        "The physical-time-stationarity-r03-r17-v3 criteria change only, as declared by each historical_independent_reviews record.",
        "Pressure-transfer methods explicitly enumerated in scientific_products_method_review.",
        "Local-field eddy-anisotropy methods explicitly enumerated in scientific_products_method_review.",
        "Replay security only for the exact checks and bindings enumerated in replay_tool_promotion_review and the scientific-replay-security method approval.",
    ],
    "full_scope_independent_review_complete": False,
    "independent_reviewer_identity_claimed_for_current_corrections": False,
}
DEPRECATED_INFERENTIAL_SCIENCE_FIELDS = frozenset({
    "holm_significant",
    "holm_threshold",
    "p_value",
    "pvalue",
    "two_sided_p",
    "z_score",
})
RETAINED_DISCRIMINANT_STATE_SCOPE = "retained_cell_centered_snapshots"
RETAINED_DISCRIMINANT_DIRECTION_SCOPE = "three_coordinate_normal_directions"
# This exact authenticated audit-script identity evaluates the frozen defective
# implementation formula.  It predates explicit formula_id result provenance.
KNOWN_LEGACY_HYPERBOLICITY_AUDIT_SCRIPTS = frozenset({
    "1b1682fece477445867b9d84ef3a1678c7346e099b927817341aa74bc8660e96",
})
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
    "hyperbolic": "#3a923a",
    "nonnegative_discriminant": "#3a923a",
    "warning": "#e6a83a",
    "warnings": "#e6a83a",
    "restricted": "#8c6bb1",
    "configuration": "#b8c4cc",
    "incomplete": "#4c78a8",
    "inconclusive": "#9e9e9e",
    "blocked_out_of_scope": "#9e9e9e",
    "fail": "#c43c39",
    "negative": "#c43c39",
    "nonfinite": "#c43c39",
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


def csv_text_value(value: object) -> str:
    """Format one deterministic, round-trip-safe CSV table cell."""

    if isinstance(value, float):
        return repr(value) if math.isfinite(value) else "--"
    if isinstance(value, (list, tuple)):
        return "; ".join(csv_text_value(item) for item in value)
    return text_value(value)


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


def hyperbolicity_manifest_errors(
    path: Path, record: dict[str, Any]
) -> list[str]:
    """Authenticate one retained-state selection and any completed audit result."""

    label = f"hyperbolicity manifest {path}"
    errors: list[str] = []
    if record.get("schema_version") != 1:
        errors.append(f"{label} has unsupported schema_version")
    case_id = record.get("case_id")
    if case_id not in CASE_IDS:
        errors.append(f"{label} does not identify a Stage I case")
    for name in ("case_lineage", "snapshot_index", "audit_script", "launcher"):
        errors.extend(verify_file_binding(record.get(name), f"{label} {name}"))

    selected_value = record.get("selected_snapshots")
    selected = (
        selected_value
        if isinstance(selected_value, list)
        else [record["selected_snapshot"]]
        if isinstance(record.get("selected_snapshot"), dict)
        else []
    )
    coverage = record.get("snapshot_coverage")
    if coverage is None:
        return errors
    if not isinstance(coverage, dict) or not selected or not all(
        isinstance(snapshot, dict) for snapshot in selected
    ):
        errors.append(f"{label} all-snapshot coverage declaration is malformed")
        return errors
    policy = record.get("snapshot_policy")
    if policy not in {"all", "latest"} or coverage.get("snapshot_policy") != policy:
        errors.append(f"{label} snapshot policy differs from coverage")
    selected_times = [
        as_float(snapshot.get("time")) for snapshot in selected
        if isinstance(snapshot, dict)
    ]
    selected_positions = [
        snapshot.get("index_position") for snapshot in selected
        if isinstance(snapshot, dict)
    ]
    selected_digest = hashlib.sha256(canonical_json(selected)).hexdigest()
    complete_count = coverage.get("snapshot_index_complete_count")
    if (
        coverage.get("selected_snapshot_count") != len(selected)
        or coverage.get("selected_snapshot_positions") != selected_positions
        or selected_times.count(None) > 0
        or coverage.get("selected_snapshot_times") != selected_times
        or coverage.get("selected_snapshots_sha256") != selected_digest
        or not isinstance(complete_count, int)
        or complete_count < len(selected)
        or not isinstance(
            coverage.get("all_complete_retained_snapshots_selected"), bool
        )
    ):
        errors.append(f"{label} selected retained-snapshot coverage differs")
    if policy == "all" and (
        coverage.get("all_complete_retained_snapshots_selected") is not True
        or complete_count != len(selected)
    ):
        errors.append(f"{label} does not select every complete retained snapshot")

    result_record = record.get("result")
    if not isinstance(result_record, dict):
        return errors
    if (
        result_record.get("required_coverage")
        != "exactly_once_per_selected_snapshot"
        or result_record.get("expected_snapshot_count") != len(selected)
        or result_record.get("expected_selected_snapshots_sha256") != selected_digest
    ):
        errors.append(f"{label} result coverage requirement differs")
        return errors
    result_path_value = result_record.get("path")
    result_sha_value = result_record.get("sha256_path")
    if not isinstance(result_path_value, str) or not isinstance(result_sha_value, str):
        errors.append(f"{label} result paths are malformed")
        return errors
    result_path = Path(result_path_value).expanduser().absolute()
    result_sha_path = Path(result_sha_value).expanduser().absolute()
    if (
        result_path != path.parent / "result.json"
        or result_sha_path != path.parent / "result.sha256"
    ):
        errors.append(f"{label} result paths differ from manifest attempt directory")
        return errors
    if not result_path.is_file() and not result_sha_path.is_file():
        return errors
    if not result_path.is_file() or not result_sha_path.is_file():
        errors.append(f"{label} has an incomplete result/sidecar pair")
        return errors
    try:
        declared_sha = result_sha_path.read_text(encoding="utf-8").split()[0]
        result = load_json(result_path)
    except (OSError, IndexError, json.JSONDecodeError, PublicationError) as error:
        errors.append(f"{label} result is unreadable: {error}")
        return errors
    if (
        SHA256_PATTERN.fullmatch(declared_sha) is None
        or sha256_file(result_path) != declared_sha
    ):
        errors.append(f"{label} result SHA-256 sidecar differs")
        return errors
    if result.get("record_type") == HYPERBOLICITY_JOB_RECORD_TYPE:
        errors.append(f"{label} result recursively identifies as a job manifest")
        return errors
    result_errors = audit_record_errors(result_path, result)
    if result_errors:
        errors.extend(result_errors)
        return errors
    result_snapshots = result.get("snapshots")
    if not isinstance(result_snapshots, list) or len(result_snapshots) != len(selected):
        errors.append(f"{label} result snapshot count differs")
        return errors
    result_provenance = result.get("provenance")
    audit_binding = record.get("audit_script")
    manifest_formula_id = record.get("formula_id", "qualified-legacy")
    result_formula_id = (
        result_provenance.get("formula_id", "qualified-legacy")
        if isinstance(result_provenance, dict) else None
    )
    selected_patterns = [
        snapshot.get("audit_input_pattern")
        for snapshot in selected if isinstance(snapshot, dict)
    ]
    if (
        not isinstance(result_provenance, dict)
        or not isinstance(audit_binding, dict)
        or result_provenance.get("script_path") != audit_binding.get("path")
        or result_provenance.get("script_sha256") != audit_binding.get("sha256")
        or result_provenance.get("input_patterns") != selected_patterns
        or result_provenance.get("hash_inputs") is not True
        or manifest_formula_id not in HYPERBOLICITY_FORMULA_IDS
        or result_formula_id != manifest_formula_id
    ):
        errors.append(f"{label} result provenance differs from selected coverage")
        return errors

    def profile(snapshot: dict[str, Any]) -> tuple[tuple[object, ...], ...] | None:
        rank_files = snapshot.get("rank_files")
        if not isinstance(rank_files, list) or not all(
            isinstance(rank_file, dict) for rank_file in rank_files
        ):
            return None
        if not all(
            isinstance(rank_file.get("path"), str)
            and isinstance(rank_file.get("size_bytes"), int)
            and isinstance(rank_file.get("mtime_ns"), int)
            for rank_file in rank_files
        ):
            return None
        return tuple(sorted(
            (
                rank_file.get("path"),
                rank_file.get("size_bytes"),
                rank_file.get("mtime_ns"),
            )
            for rank_file in rank_files
        ))

    selected_by_profile: dict[tuple[tuple[object, ...], ...], float] = {}
    for snapshot in selected:
        assert isinstance(snapshot, dict)
        snapshot_profile = profile(snapshot)
        time = as_float(snapshot.get("time"))
        if (
            snapshot_profile is None
            or time is None
            or snapshot_profile in selected_by_profile
        ):
            errors.append(f"{label} selected snapshot inventory is ambiguous")
            return errors
        selected_by_profile[snapshot_profile] = time
    for snapshot in result_snapshots:
        if not isinstance(snapshot, dict):
            errors.append(f"{label} result snapshot is malformed")
            return errors
        rank_files = snapshot.get("rank_files")
        if (
            snapshot.get("active_cgl_signal_speed") is not True
            or snapshot.get("ranks_contiguous_from_zero") is not True
            or not isinstance(rank_files, list)
            or any(
                not isinstance(rank_file, dict)
                or SHA256_PATTERN.fullmatch(str(rank_file.get("sha256"))) is None
                for rank_file in rank_files
            )
        ):
            errors.append(f"{label} result lacks active hash-bound rank coverage")
            return errors
        snapshot_profile = profile(snapshot)
        selected_time = selected_by_profile.pop(snapshot_profile, None)
        observed_time = as_float(snapshot.get("time"))
        if (
            selected_time is None
            or observed_time is None
            or not math.isclose(
                selected_time, observed_time, rel_tol=0.0, abs_tol=1.0e-12
            )
        ):
            errors.append(f"{label} result differs from selected retained coverage")
            return errors
    if selected_by_profile:
        errors.append(f"{label} result omits selected retained snapshots")
        return errors
    record["_publication_hyperbolicity_result"] = result
    record["_publication_hyperbolicity_formula_id"] = manifest_formula_id
    record["_publication_hyperbolicity_result_path"] = str(result_path)
    record["_publication_hyperbolicity_result_sha_path"] = str(result_sha_path)
    return errors


def audit_record_errors(path: Path, record: dict[str, Any]) -> list[str]:
    """Return provenance and retained-input freshness errors for one audit."""

    if record.get("record_type") == HYPERBOLICITY_JOB_RECORD_TYPE:
        return hyperbolicity_manifest_errors(path, record)
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


def deprecated_inferential_science_fields(value: object) -> set[str]:
    """Return deprecated population-inference fields found in reviewed science."""

    found: set[str] = set()
    if isinstance(value, dict):
        found.update(set(value) & DEPRECATED_INFERENTIAL_SCIENCE_FIELDS)
        for nested_value in value.values():
            found.update(deprecated_inferential_science_fields(nested_value))
    elif isinstance(value, list):
        for nested_value in value:
            found.update(deprecated_inferential_science_fields(nested_value))
    return found


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
    if record.get("current_science_scope_limitation") != CURRENT_SCIENCE_SCOPE_LIMITATION:
        errors.append(f"{label} lacks the exact current science scope limitation")
    if record.get("active_passive_intervention_scope") != ACTIVE_PASSIVE_INTERVENTION_SCOPE:
        errors.append(f"{label} lacks the exact active/passive intervention scope")
    deprecated_fields = deprecated_inferential_science_fields(record)
    if deprecated_fields:
        errors.append(
            f"{label} contains deprecated population-inference fields: "
            + ", ".join(sorted(deprecated_fields))
        )
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
                if (
                    family == "active_passive"
                    and contrast.get("intervention_scope")
                    != ACTIVE_PASSIVE_INTERVENTION_SCOPE
                ):
                    errors.append(
                        f"{label} contrast {family}.{name} lacks the exact "
                        "active/passive intervention scope"
                    )
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


def json_payload(value: object) -> bytes:
    """Return canonical human-readable JSON bytes."""

    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def write_json(path: Path, value: object) -> None:
    """Write canonical human-readable JSON."""

    atomic_write(path, json_payload(value))


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


def builtin_acceptance_roots(analysis: Path) -> tuple[Path, ...]:
    """Return built-in acceptance roots that publication output may not overlap."""

    return (
        analysis / "scientific-acceptance",
        analysis / "acceptance",
        analysis / "audits",
        analysis / "scientific-audits",
    )


def discover_acceptance_paths(analysis: Path, extras: Iterable[Path]) -> list[Path]:
    """Return deterministic candidate acceptance and audit JSON paths."""

    candidates: set[Path] = set()
    roots = [*builtin_acceptance_roots(analysis), *extras]
    for root in roots:
        if root.is_file():
            candidates.add(root.absolute())
        elif root.is_dir():
            candidates.update(path.absolute() for path in root.rglob("*.json"))
    return sorted(candidates)


def contains_hyperbolicity_evidence(record: dict[str, Any]) -> bool:
    """Return whether one JSON record carries retained-state hyperbolicity data."""

    if record.get("record_type") == HYPERBOLICITY_JOB_RECORD_TYPE:
        return True
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
            for private_name in (
                "_publication_hyperbolicity_result_path",
                "_publication_hyperbolicity_result_sha_path",
            ):
                private_path = loaded.get(private_name)
                if isinstance(private_path, str):
                    source_paths.add(Path(private_path).absolute())
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
    elif name == "c_b2":
        if "b2" not in history or "b4" not in history:
            return None
        values = [
            fourth * vol / max(second * second, 1.0e-300) - 1.0
            for second, fourth, vol in zip(history["b2"], history["b4"], volume)
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


def authenticated_strict_failure_records(
    data: PublicationData, case_id: str
) -> list[dict[str, Any]]:
    """Return authenticated retained strict-failure records for one scoped case."""

    if case_id not in {"R14", "R15"}:
        return []
    case = data.cases[case_id]
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
        binding = nested(data.science_record, f"provenance.case_acceptance.{case_id}")
        path = binding_path(binding)
        if path is not None and not verify_file_binding(
            binding, f"{case_id} science-bound case acceptance"
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
            or record.get("case_id") != case_id
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
                    provenance.get(name),
                    f"{case_id} strict-failure {name} provenance",
                )
                for name in ("manifest", "run_exit_code", "slurm_log")
            )
        ):
            continue
        authenticated.append(record)
    return sorted(
        authenticated,
        key=lambda record: (
            float(record["failure_time"]),
            str(record.get("job_id", "")),
        ),
    )


def strict_failure_disposition(data: PublicationData, case_id: str) -> str:
    """Return one compact authenticated strict-failure disposition."""

    authenticated = authenticated_strict_failure_records(data, case_id)
    if not authenticated:
        return "unavailable/inconclusive"
    selected = authenticated[0]
    time = float(selected["failure_time"])
    counters = selected["failure_counters"]
    hard_bound = float(counters["lf_hardbd"])
    job_id = selected.get("job_id")
    details = ["fail", f"t={time:.8g}", f"hard_bound={hard_bound:.8g}"]
    if isinstance(job_id, (str, int)):
        details.append(f"job={job_id}")
    return "; ".join(details)


def r15_strict_failure_disposition(data: PublicationData) -> str:
    """Return authenticated strict-R15 failure evidence without inventing details."""

    return strict_failure_disposition(data, "R15")


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
        "hyperbolic": "hyperbolic",
        "nonnegative_discriminant": "D >= 0",
        "warnings": "warn",
        "warning": "warn",
        "restricted": "restricted",
        "configuration": "configured",
        "incomplete": "partial",
        "inconclusive": "inconclusive",
        "blocked_out_of_scope": "blocked",
        "fail": "fail",
        "negative": "negative",
        "nonfinite": "nonfinite",
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
        writer.writerow({
            column: csv_text_value(row.get(column)) for column in columns
        })
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
        "pass": 0, "clean": 0, "hyperbolic": 0, "nonnegative_discriminant": 0,
        "warning": 1, "warnings": 1,
        "restricted": 2,
        "configuration": 3,
        "incomplete": 4,
        "inconclusive": 5, "blocked_out_of_scope": 5,
        "fail": 6, "structural_error": 6, "negative": 6, "nonfinite": 6,
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
    """Flatten authenticated descriptive reviewed-science contrasts."""

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
            intervention_scope = (
                contrast.get("intervention_scope")
                if family == "active_passive"
                else None
            )
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
                    "standardized_effect": metric.get("standardized_effect"),
                    "intervention_estimand": (
                        intervention_scope.get("estimand")
                        if isinstance(intervention_scope, dict) else None
                    ),
                    "intervention_enabled_components": (
                        intervention_scope.get("enabled_components")
                        if isinstance(intervention_scope, dict) else None
                    ),
                    "excluded_interpretation": (
                        intervention_scope.get("excluded_interpretation")
                        if isinstance(intervention_scope, dict) else None
                    ),
                    "intervention_declaration": (
                        intervention_scope.get("declaration")
                        if isinstance(intervention_scope, dict) else None
                    ),
                    "reason": metric.get("reason", contrast.get("reason")),
                    "claim_scope": r15_science_scope(data, name, left, right),
                    "inference_scope": (
                        intervention_scope.get("claim_scope")
                        if isinstance(intervention_scope, dict) else
                        "descriptive_within_realization"
                    ),
                    "current_science_scope_disposition": (
                        CURRENT_SCIENCE_SCOPE_LIMITATION["disposition"]
                    ),
                    "full_scope_independent_review_complete": False,
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


def authenticated_case_diagnostics(
    data: PublicationData, case_id: str
) -> dict[str, Any] | None:
    """Return current diagnostics only when bound by authenticated reviewed science."""

    if (
        not isinstance(data.science_record, dict)
        or data.science_record.get("_publication_evidence_validated") is not True
    ):
        return None
    binding = nested(data.science_record, f"provenance.case_diagnostics.{case_id}")
    path = data.analysis / "cases" / case_id / "diagnostics.json"
    if binding_freshness_errors(
        binding, path, f"{case_id} reviewed-science diagnostics"
    ):
        return None
    try:
        diagnostics = load_json(path)
    except (OSError, json.JSONDecodeError, PublicationError):
        return None
    data.source_paths.add(path.absolute())
    return diagnostics


def authenticated_case_lineage(
    data: PublicationData, case_id: str
) -> dict[str, Any] | None:
    """Return current lineage only when bound by authenticated evidence."""

    case = data.cases[case_id]
    path = case.lineage_path
    if path is None:
        return None
    bindings: list[object] = []
    direct = authenticated_direct_acceptance(case)
    if isinstance(direct, dict):
        bindings.append(nested(direct, "provenance.lineage"))
    if (
        isinstance(data.science_record, dict)
        and data.science_record.get("_publication_evidence_validated") is True
    ):
        bindings.append(
            nested(data.science_record, f"provenance.case_lineages.{case_id}")
        )
    if all(
        binding_freshness_errors(binding, path, f"{case_id} authenticated lineage")
        for binding in bindings
    ):
        return None
    try:
        lineage = load_json(path)
    except (OSError, json.JSONDecodeError, PublicationError):
        return None
    data.source_paths.add(path.absolute())
    return lineage


def authenticated_direct_acceptance(case: CaseRecord) -> dict[str, Any] | None:
    """Return one selected direct-fast case record only after provenance validation."""

    record = case.direct_acceptance
    return (
        record
        if isinstance(record, dict)
        and record.get("_publication_evidence_validated") is True
        else None
    )


def case_gate(case: CaseRecord, name: str) -> dict[str, Any] | None:
    """Return one authenticated reviewed per-case gate."""

    evidence = case.acceptance
    if (
        not isinstance(evidence, dict)
        or evidence.get("_publication_evidence_validated") is not True
    ):
        return None
    gates = evidence.get("gates")
    if not isinstance(gates, list):
        return None
    for gate in gates:
        if isinstance(gate, dict) and gate.get("name") == name:
            return gate
    return None


def energy_closure_summary(case: CaseRecord) -> dict[str, object]:
    """Return authenticated active-energy closure evidence without inventing values."""

    if case.case_id not in ACTIVE_ENERGY_CASES:
        return {
            "result": "not_applicable",
            "maximum_increment_normalized_residual": None,
            "maximum_state_normalized_mismatch": None,
        }
    gate = case_gate(case, "active_energy_closure")
    if gate is None:
        return {
            "result": "inconclusive",
            "maximum_increment_normalized_residual": None,
            "maximum_state_normalized_mismatch": None,
        }
    windows = nested(gate, "observations.windows")
    values = (
        [windows[name] for name in ("whole_lineage", "developed")]
        if isinstance(windows, dict)
        and all(isinstance(windows.get(name), dict) for name in (
            "whole_lineage", "developed"
        ))
        else []
    )
    increment = [
        value
        for value in (
            as_float(record.get("increment_normalized_residual"))
            for record in values if isinstance(record, dict)
        )
        if value is not None
    ]
    state = [
        value
        for value in (
            as_float(record.get("state_normalized_mismatch"))
            for record in values if isinstance(record, dict)
        )
        if value is not None
    ]
    return {
        "result": validated_result(gate.get("result")),
        "maximum_increment_normalized_residual": max(increment) if increment else None,
        "maximum_state_normalized_mismatch": max(state) if state else None,
    }


def numerical_health_provenance_rows(
    data: PublicationData,
) -> list[dict[str, object]]:
    """Return compact authenticated all-case numerical-health and provenance rows."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        case = data.cases[case_id]
        direct = authenticated_direct_acceptance(case)
        diagnostics = authenticated_case_diagnostics(data, case_id)
        health = diagnostics.get("health") if isinstance(diagnostics, dict) else None
        direct_health = direct.get("health") if isinstance(direct, dict) else None
        fatal = (
            direct_health.get("fatal_counter_maxima")
            if isinstance(direct_health, dict) else None
        )
        fatal_values = [
            value
            for value in (
                as_float(raw) for raw in fatal.values()
            )
            if value is not None
        ] if isinstance(fatal, dict) and set(fatal) == set(FATAL_LF_COUNTERS) else []
        mass = health.get("mass_relative_drift") if isinstance(health, dict) else None
        mass_values = [
            value
            for value in (
                as_float(raw) for raw in mass.values()
            )
            if value is not None
        ] if isinstance(mass, dict) and set(mass) == {"mhd", "user"} else []
        energy = energy_closure_summary(case)
        strict_records = authenticated_strict_failure_records(data, case_id)
        strict_record = strict_records[0] if strict_records else None
        direct_scope = direct.get("scope") if isinstance(direct, dict) else None
        ct_record = ct_case_record(data, case_id)
        ct_native = (
            ct_record.get("native_restart_ct")
            if isinstance(ct_record, dict) else None
        )
        direct_complete = (
            direct_health.get("complete_to_target")
            if isinstance(direct_health, dict) else None
        )
        direct_final = (
            as_float(direct_health.get("observed_final_time"))
            if isinstance(direct_health, dict) else None
        )
        fatal_status = (
            "pass" if fatal_values and max(fatal_values) == 0.0
            else "fail" if fatal_values
            else "inconclusive"
        )
        reporter_provenance = (
            "authenticated" if diagnostics is not None else "inconclusive"
        )
        acceptance_provenance = (
            "authenticated" if direct is not None else "inconclusive"
        )
        science_provenance = (
            "authenticated"
            if isinstance(data.science_record, dict)
            and data.science_record.get("_publication_evidence_validated") is True
            and case_id in (data.science_record.get("selected_cases") or [])
            else "inconclusive"
        )
        ct_provenance = (
            "authenticated" if isinstance(ct_record, dict) else "inconclusive"
        )
        rows.append({
            "case_id": case_id,
            "completion": (
                "pass" if direct_complete is True
                else "incomplete" if direct_complete is False
                else "inconclusive"
            ),
            "final_time": direct_final,
            "fatal_lf_counters": fatal_status,
            "fatal_counter_maximum": max(fatal_values) if fatal_values else None,
            "mass_relative_drift_maximum": max(mass_values) if mass_values else None,
            "mhd_user_mass_relative_mismatch": (
                health.get("mhd_user_mass_relative_mismatch")
                if isinstance(health, dict) else None
            ),
            "active_energy_closure": energy["result"],
            "energy_increment_residual_maximum": energy[
                "maximum_increment_normalized_residual"
            ],
            "energy_state_mismatch_maximum": energy[
                "maximum_state_normalized_mismatch"
            ],
            "direct_ct_numerical": ct_case_status(data, case_id),
            "maximum_normalized_ct_divb": (
                ct_native.get("maximum_normalized_ct_divb")
                if isinstance(ct_native, dict) else None
            ),
            "claim_scope": (
                direct_scope.get("classification")
                if isinstance(direct_scope, dict)
                else "restricted" if scope_status(case_id) == "restricted"
                else "standard"
            ),
            "acceptance": acceptance_status(data, case),
            "reviewed_science": science_case_status(data, case_id),
            "strict_failure": (
                strict_failure_disposition(data, case_id)
                if case_id in {"R14", "R15"} else "not_applicable"
            ),
            "strict_failure_counters": (
                strict_record.get("failure_counters")
                if isinstance(strict_record, dict) else None
            ),
            "reporter_diagnostics_provenance": reporter_provenance,
            "direct_acceptance_provenance": acceptance_provenance,
            "reviewed_science_provenance": science_provenance,
            "direct_ct_provenance": ct_provenance,
            "provenance_summary": (
                f"reporter={reporter_provenance}; "
                f"acceptance={acceptance_provenance}; "
                f"science={science_provenance}; CT={ct_provenance}"
            ),
        })
    return rows


def primary_full_window_scalar_rows(
    data: PublicationData,
) -> list[dict[str, object]]:
    """Return compact authenticated primary full-window scalar estimates."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        case = data.cases[case_id]
        direct = authenticated_direct_acceptance(case)
        statistics = (
            direct.get("history_statistics") if isinstance(direct, dict) else None
        )
        direct_scope = direct.get("scope") if isinstance(direct, dict) else None
        for metric in PRIMARY_SCALAR_METRICS:
            record = statistics.get(metric) if isinstance(statistics, dict) else None
            full = nested(record, "windows.full.statistics")
            stationarity = record.get("stationarity") if isinstance(record, dict) else None
            confidence = (
                full.get("confidence_interval_95") if isinstance(full, dict) else None
            )
            available = (
                isinstance(full, dict)
                and as_float(full.get("mean")) is not None
                and as_float(full.get("standard_error")) is not None
                and as_float(full.get("effective_sample_count")) is not None
                and isinstance(confidence, list)
                and len(confidence) == 2
                and all(as_float(value) is not None for value in confidence)
            )
            rows.append({
                "case_id": case_id,
                "metric": metric,
                "history": record.get("history") if isinstance(record, dict) else None,
                "column": record.get("column") if isinstance(record, dict) else None,
                "availability": "available" if available else "inconclusive",
                "mean": full.get("mean") if available else None,
                "standard_error": (
                    full.get("standard_error") if available else None
                ),
                "ci95_lower": (
                    confidence[0] if available else None
                ),
                "ci95_upper": (
                    confidence[1] if available else None
                ),
                "effective_sample_count": (
                    full.get("effective_sample_count") if available else None
                ),
                "stationarity": (
                    validated_result(stationarity.get("result"))
                    if isinstance(stationarity, dict) else "inconclusive"
                ),
                "sampling_adequacy": (
                    record.get("sampling_adequacy")
                    if isinstance(record, dict) else "inconclusive"
                ),
                "acceptance": acceptance_status(data, case),
                "claim_scope": (
                    direct_scope.get("classification")
                    if isinstance(direct_scope, dict)
                    else "restricted" if scope_status(case_id) == "restricted"
                    else "standard"
                ),
            })
    return rows


def exact_nonnegative_integer(value: object) -> int | None:
    """Return one exact nonnegative integer without converting it to float."""

    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 0 else None
    if isinstance(value, float) and math.isfinite(value) and value.is_integer():
        return int(value) if value >= 0.0 else None
    return None


def snapshot_hyperbolicity_summary(
    snapshots: object,
) -> dict[str, object] | None:
    """Aggregate retained coordinate-direction discriminant dispositions."""

    if not isinstance(snapshots, list) or not snapshots:
        return None
    evaluated = 0
    negative = 0
    nonfinite = 0
    minimum: float | None = None
    aggregate_count = 0
    for snapshot in snapshots:
        aggregate = snapshot.get("aggregate") if isinstance(snapshot, dict) else None
        if not isinstance(aggregate, dict):
            continue
        evaluated_value = exact_nonnegative_integer(aggregate.get("evaluated"))
        negative_value = exact_nonnegative_integer(aggregate.get("negative"))
        nonfinite_value = exact_nonnegative_integer(
            aggregate.get("nonfinite_discriminant")
        )
        minimum_value = as_float(aggregate.get("minimum"))
        if (
            evaluated_value is None
            or negative_value is None
            or nonfinite_value is None
            or negative_value > evaluated_value
            or nonfinite_value > evaluated_value
        ):
            continue
        aggregate_count += 1
        evaluated += evaluated_value
        negative += negative_value
        nonfinite += nonfinite_value
        if minimum_value is not None and (
            minimum is None or minimum_value < minimum
        ):
            minimum = minimum_value
    if aggregate_count == 0 or evaluated <= 0:
        return None
    disposition = (
        "nonfinite"
        if nonfinite > 0
        else "negative"
        if negative > 0 or (minimum is not None and minimum < 0.0)
        else "nonnegative_discriminant"
    )
    return {
        "result": disposition,
        "negative_discriminant_fraction": negative / evaluated,
        "negative_discriminant_count": negative,
        "cell_direction_evaluations": evaluated,
        "minimum_discriminant": minimum,
        "nonfinite_discriminant_count": nonfinite,
    }


def authenticated_formula_summary(
    result: object, manifest: object = None
) -> dict[str, object]:
    """Return formula identity and executable compatibility from one bound result."""

    provenance = (
        result.get("provenance") if isinstance(result, dict) else None
    )
    provenance = provenance if isinstance(provenance, dict) else {}

    def result_or_provenance(name: str) -> object | None:
        value = result.get(name) if isinstance(result, dict) else None
        return value if value is not None else provenance.get(name)

    explicit_formula_id = result_or_provenance("formula_id")
    formula_id = (
        explicit_formula_id
        if isinstance(explicit_formula_id, str) and explicit_formula_id
        else None
    )
    formula_id_provenance = (
        "authenticated_result" if formula_id is not None else "inconclusive"
    )
    manifest_formula_id = (
        manifest.get("_publication_hyperbolicity_formula_id")
        if isinstance(manifest, dict) else None
    )
    if formula_id is None and manifest_formula_id in HYPERBOLICITY_FORMULA_IDS:
        formula_id = manifest_formula_id
        formula_id_provenance = "authenticated_legacy_default_contract"
    script_sha = provenance.get("script_sha256")
    if (
        formula_id is None
        and isinstance(script_sha, str)
        and script_sha in KNOWN_LEGACY_HYPERBOLICITY_AUDIT_SCRIPTS
    ):
        formula_id = "qualified-legacy"
        formula_id_provenance = "authenticated_known_legacy_audit_script"
    formula_family = HYPERBOLICITY_FORMULA_FAMILIES.get(str(formula_id))

    explicit_formula_provenance = result_or_provenance("formula_provenance")
    if explicit_formula_provenance not in (None, "", [], {}):
        formula_provenance: object | None = explicit_formula_provenance
    else:
        formula_provenance_members = {
            name: provenance[name]
            for name in (
                "formula_reference",
                "script_path",
                "script_version",
                "script_sha256",
            )
            if provenance.get(name) is not None
        }
        formula_provenance = formula_provenance_members or None
    formula_provenance_status = (
        "authenticated" if formula_provenance is not None else "inconclusive"
    )

    executable_formula_id_value = result_or_provenance("executable_formula_id")
    executable_formula_id = (
        executable_formula_id_value
        if isinstance(executable_formula_id_value, str)
        and executable_formula_id_value
        else None
    )
    executable_formula_id_provenance = (
        "authenticated_result"
        if executable_formula_id is not None else "inconclusive"
    )

    compatibility_value = result_or_provenance("formula_executable_compatibility")
    compatibility_reason = result_or_provenance(
        "formula_executable_compatibility_reason"
    )
    compatibility_provenance = "inconclusive"
    executable_formula_family = HYPERBOLICITY_FORMULA_FAMILIES.get(
        str(executable_formula_id)
    )
    identity_compatibility = (
        "compatible" if formula_family == executable_formula_family else "incompatible"
    ) if formula_family is not None and executable_formula_family is not None else None
    if isinstance(compatibility_value, dict):
        compatibility_reason = (
            compatibility_value.get("reason")
            if compatibility_value.get("reason") is not None
            else compatibility_reason
        )
        compatibility_value = (
            compatibility_value.get("status")
            if compatibility_value.get("status") is not None
            else compatibility_value.get("result")
        )
    if isinstance(compatibility_value, bool):
        compatibility = "compatible" if compatibility_value else "incompatible"
        compatibility_provenance = "authenticated_result"
    elif compatibility_value in {"compatible", "incompatible", "inconclusive"}:
        compatibility = str(compatibility_value)
        compatibility_provenance = "authenticated_result"
    elif identity_compatibility is not None:
        compatibility = identity_compatibility
        compatibility_provenance = "derived_from_authenticated_formula_ids"
        compatibility_reason = (
            "authenticated audit and executable formula identities agree"
            if compatibility == "compatible"
            else "authenticated audit and executable formula identities differ"
        )
    else:
        compatibility = "inconclusive"
    if compatibility == "compatible" and identity_compatibility == "incompatible":
        compatibility = "incompatible"
        compatibility_provenance = "authenticated_result_conflicts_with_formula_ids"
        compatibility_reason = (
            "authenticated compatibility declaration conflicts with differing audit "
            "and executable formula identities"
        )
    if not isinstance(compatibility_reason, str) or not compatibility_reason:
        compatibility_reason = (
            "authenticated formula/executable compatibility is unavailable"
            if compatibility == "inconclusive"
            else "authenticated result declares formula/executable compatibility"
        )

    numerical = snapshot_hyperbolicity_summary(
        result.get("snapshots") if isinstance(result, dict) else None
    )
    numerical_result = (
        str(numerical["result"]) if isinstance(numerical, dict) else "inconclusive"
    )
    return {
        "formula_id": formula_id,
        "formula_disposition_family": formula_family,
        "formula_id_provenance": formula_id_provenance,
        "formula_provenance": formula_provenance,
        "formula_provenance_status": formula_provenance_status,
        "executable_formula_id": executable_formula_id,
        "executable_formula_id_provenance": executable_formula_id_provenance,
        "formula_executable_compatibility": compatibility,
        "formula_executable_compatibility_provenance": compatibility_provenance,
        "formula_executable_compatibility_reason": compatibility_reason,
        "legacy_implementation_disposition": (
            numerical_result
            if formula_family == "legacy_implementation" else "inconclusive"
        ),
        "literature_correct_disposition": (
            numerical_result
            if formula_family == "literature_correct" else "inconclusive"
        ),
    }


def strict_hyperbolic_claim_summary(
    experiment_scope: str,
    coverage_result: str,
    numerical_result: str,
    formula_id: object = None,
    formula_provenance_status: str = "inconclusive",
    formula_executable_compatibility: str = "inconclusive",
) -> dict[str, object]:
    """Gate strict claims without overextending retained directional audits."""

    formula_family = HYPERBOLICITY_FORMULA_FAMILIES.get(str(formula_id))
    if experiment_scope == "passive_delta":
        return {
            "status": "not_applicable",
            "eligible": False,
            "reason": "passive-delta case has no active CGL signal-speed claim",
        }
    if coverage_result != "pass":
        return {
            "status": "inconclusive",
            "eligible": None,
            "reason": (
                "authenticated all-retained-snapshot coverage and result are required"
            ),
        }
    if formula_family == "legacy_implementation":
        return {
            "status": "inconclusive_legacy_implementation",
            "eligible": None,
            "reason": (
                "legacy-implementation disposition neither supports nor refutes the "
                "literature-correct strict-hyperbolic claim"
            ),
        }
    if formula_family != "literature_correct":
        return {
            "status": "inconclusive_formula_identity",
            "eligible": None,
            "reason": (
                "authenticated literature-correct audit formula identity is required"
            ),
        }
    if formula_provenance_status != "authenticated":
        return {
            "status": "inconclusive_formula_provenance",
            "eligible": None,
            "reason": "authenticated literature-correct formula provenance is required",
        }
    if formula_executable_compatibility == "incompatible":
        return {
            "status": "excluded_formula_executable_incompatible",
            "eligible": False,
            "reason": (
                "literature-correct audit formula is incompatible with the audited "
                "executable or campaign"
            ),
        }
    if formula_executable_compatibility != "compatible":
        return {
            "status": "inconclusive_formula_executable_compatibility",
            "eligible": None,
            "reason": (
                "authenticated literature-correct formula/executable compatibility "
                "is required"
            ),
        }
    if numerical_result == "negative":
        return {
            "status": "excluded_negative",
            "eligible": False,
            "reason": (
                "authenticated all-snapshot audit contains negative discriminants"
            ),
        }
    if numerical_result == "nonfinite":
        return {
            "status": "excluded_nonfinite",
            "eligible": False,
            "reason": (
                "authenticated all-snapshot audit contains nonfinite discriminants"
            ),
        }
    if numerical_result not in {"nonnegative_discriminant", "hyperbolic"}:
        return {
            "status": "inconclusive",
            "eligible": None,
            "reason": (
                "authenticated coordinate-direction discriminant disposition is "
                "unavailable"
            ),
        }
    if experiment_scope == "restricted":
        return {
            "status": "excluded_experiment_scope",
            "eligible": False,
            "reason": (
                "retained cell-centered coordinate-direction discriminants are "
                "nonnegative, but experiment scope is restricted"
            ),
        }
    return {
        "status": "inconclusive_retained_coordinate_discriminant_scope",
        "eligible": None,
        "reason": (
            "compatible authenticated literature-correct discriminants are "
            "nonnegative only for retained cell-centered snapshots in the three "
            "coordinate-normal directions; this does not test reconstructed faces, "
            "intermediate states, oblique directions, full or strict hyperbolicity, "
            "or absence of sqrt(|D|) fallback"
        ),
    }


def authenticated_hyperbolicity_manifest(
    data: PublicationData, case_id: str
) -> dict[str, Any] | None:
    """Select the most complete authenticated current hyperbolicity manifest."""

    case = data.cases[case_id]

    def current(record: dict[str, Any]) -> bool:
        snapshot_path_value = nested(case.lineage, "snapshots.path")
        snapshot_path = (
            Path(snapshot_path_value)
            if isinstance(snapshot_path_value, str) else None
        )
        return (
            not binding_freshness_errors(
                record.get("case_lineage"), case.lineage_path,
                f"{case_id} hyperbolicity manifest lineage",
            )
            and not binding_freshness_errors(
                record.get("snapshot_index"), snapshot_path,
                f"{case_id} hyperbolicity manifest snapshot index",
            )
            and isinstance(record.get("_publication_source_path"), str)
            and not audit_record_errors(
                Path(str(record["_publication_source_path"])), record
            )
        )

    records = [
        record for record in data.audit_records
        if record.get("_publication_evidence_validated") is True
        and record.get("record_type") == HYPERBOLICITY_JOB_RECORD_TYPE
        and record.get("case_id") == case_id
        and current(record)
    ]
    if not records:
        return None
    return max(
        records,
        key=lambda record: (
            int(record.get("snapshot_policy") == "all"),
            int(isinstance(record.get("_publication_hyperbolicity_result"), dict)),
            int(nested(record, "snapshot_coverage.selected_snapshot_count") or 0),
            int(record.get("attempt") or 0),
            str(record.get("_publication_source_path", "")),
        ),
    )


def hyperbolicity_coverage_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return authenticated retained coordinate-direction discriminant coverage."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        if case_id not in ACTIVE_ENERGY_CASES:
            claim = strict_hyperbolic_claim_summary(
                "passive_delta", "not_applicable", "not_applicable"
            )
            rows.append({
                "case_id": case_id,
                "coverage_result": "not_applicable",
                "numerical_result": "not_applicable",
                "coverage_reason": (
                    "passive-delta case has no active CGL signal-speed audit"
                ),
                "experiment_scope": "passive_delta",
                "formula_id": None,
                "formula_disposition_family": None,
                "formula_id_provenance": "not_applicable",
                "formula_provenance": None,
                "formula_provenance_status": "not_applicable",
                "executable_formula_id": None,
                "executable_formula_id_provenance": "not_applicable",
                "formula_executable_compatibility": "not_applicable",
                "formula_executable_compatibility_provenance": "not_applicable",
                "formula_executable_compatibility_reason": (
                    "passive-delta case has no active CGL signal-speed audit"
                ),
                "legacy_implementation_disposition": "not_applicable",
                "literature_correct_disposition": "not_applicable",
                "strict_hyperbolic_claim_status": claim["status"],
                "strict_hyperbolic_claim_eligible": claim["eligible"],
                "strict_hyperbolic_claim_reason": claim["reason"],
                "audit_state_scope": "not_applicable",
                "audit_direction_scope": "not_applicable",
                "selection_provenance": "not_applicable",
                "result_provenance": "not_applicable",
            })
            continue
        manifest = authenticated_hyperbolicity_manifest(data, case_id)
        coverage = (
            manifest.get("snapshot_coverage") if isinstance(manifest, dict) else None
        )
        result = (
            manifest.get("_publication_hyperbolicity_result")
            if isinstance(manifest, dict) else None
        )
        formula = authenticated_formula_summary(result, manifest)
        snapshots = result.get("snapshots") if isinstance(result, dict) else None
        numerical = snapshot_hyperbolicity_summary(snapshots)
        snapshot_times = [
            value
            for value in (
                as_float(snapshot.get("time"))
                for snapshot in snapshots
                if isinstance(snapshot, dict)
            )
            if value is not None
        ] if isinstance(snapshots, list) else []
        policy = manifest.get("snapshot_policy") if isinstance(manifest, dict) else None
        all_selected = (
            coverage.get("all_complete_retained_snapshots_selected")
            if isinstance(coverage, dict) else None
        )
        complete_count = (
            coverage.get("snapshot_index_complete_count")
            if isinstance(coverage, dict) else None
        )
        selected_count = (
            coverage.get("selected_snapshot_count")
            if isinstance(coverage, dict) else None
        )
        result_count = len(snapshots) if isinstance(snapshots, list) else None
        coverage_pass = (
            policy == "all"
            and all_selected is True
            and isinstance(complete_count, int)
            and complete_count > 0
            and complete_count == selected_count == result_count
        )
        experiment_scope = (
            "restricted" if scope_status(case_id) == "restricted" else "standard"
        )
        numerical_result = (
            str(numerical["result"])
            if isinstance(numerical, dict) else "inconclusive"
        )
        claim = strict_hyperbolic_claim_summary(
            experiment_scope,
            "pass" if coverage_pass else "inconclusive",
            numerical_result,
            formula["formula_id"],
            str(formula["formula_provenance_status"]),
            str(formula["formula_executable_compatibility"]),
        )
        coverage_reason = (
            "authenticated result covers every complete retained cell-centered "
            "snapshot in the three coordinate-normal directions"
            if coverage_pass
            else "no authenticated current directional-discriminant selection"
            if not isinstance(manifest, dict)
            else "authenticated manifest lacks all-snapshot coverage metadata"
            if not isinstance(coverage, dict)
            else "authenticated selection is not snapshot_policy=all"
            if policy != "all"
            else "authenticated all-snapshot selection lacks a completed bound result"
        )
        rows.append({
            "case_id": case_id,
            "coverage_result": "pass" if coverage_pass else "inconclusive",
            "numerical_result": numerical_result,
            "snapshot_policy": policy,
            "complete_retained_snapshot_count": complete_count,
            "selected_snapshot_count": selected_count,
            "audited_snapshot_count": result_count,
            "all_complete_retained_snapshots_selected": all_selected,
            "time_first": min(snapshot_times) if snapshot_times else None,
            "time_last": max(snapshot_times) if snapshot_times else None,
            "negative_discriminant_count": (
                numerical["negative_discriminant_count"]
                if isinstance(numerical, dict) else None
            ),
            "nonfinite_discriminant_count": (
                numerical["nonfinite_discriminant_count"]
                if isinstance(numerical, dict) else None
            ),
            "cell_direction_evaluations": (
                numerical["cell_direction_evaluations"]
                if isinstance(numerical, dict) else None
            ),
            "minimum_discriminant": (
                numerical["minimum_discriminant"]
                if isinstance(numerical, dict) else None
            ),
            "coverage_reason": coverage_reason,
            "experiment_scope": experiment_scope,
            **formula,
            "strict_hyperbolic_claim_status": claim["status"],
            "strict_hyperbolic_claim_eligible": claim["eligible"],
            "strict_hyperbolic_claim_reason": claim["reason"],
            "audit_state_scope": RETAINED_DISCRIMINANT_STATE_SCOPE,
            "audit_direction_scope": RETAINED_DISCRIMINANT_DIRECTION_SCOPE,
            "selection_provenance": (
                "authenticated" if isinstance(manifest, dict) else "inconclusive"
            ),
            "result_provenance": (
                "authenticated" if isinstance(result, dict) else "inconclusive"
            ),
        })
    return rows


def signed_lf_cap_work_ledger_rows(
    data: PublicationData,
) -> list[dict[str, object]]:
    """Return signed applied ledgers and distinct snapshot reconstructions."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        diagnostics = authenticated_case_diagnostics(data, case_id)
        lf = nested(diagnostics, "windows.steady.lf_history")
        applied_heat = (
            lf.get("applied_heat_flux_work") if isinstance(lf, dict) else None
        )
        applied_pressure = (
            lf.get("applied_pressure_work") if isinstance(lf, dict) else None
        )
        caps = lf.get("heat_flux_cap_fractions") if isinstance(lf, dict) else None
        ensemble = authenticated_snapshot_ensemble(data, case_id)
        pressure = nested(ensemble, "pressure_work_decomposition")
        heat_proxy = nested(ensemble, "heat_flux_transport_proxy")
        pressure_integral = nested(pressure, "time_integral_estimate")
        heat_integral = nested(heat_proxy, "time_integral_estimate")
        applied_heat_available = (
            isinstance(applied_heat, dict)
            and applied_heat.get("signed") is True
        )
        applied_pressure_available = (
            isinstance(applied_pressure, dict)
            and applied_pressure.get("signed") is True
        )
        pressure_available = (
            isinstance(pressure, dict) and pressure.get("available") is True
        )
        heat_proxy_available = (
            isinstance(heat_proxy, dict) and heat_proxy.get("available") is True
        )
        rows.append({
            "case_id": case_id,
            "availability": (
                "available"
                if applied_heat_available or applied_pressure_available
                or pressure_available or heat_proxy_available
                else "inconclusive"
            ),
            "diagnostics_provenance": (
                "authenticated" if isinstance(diagnostics, dict) else "inconclusive"
            ),
            "applied_heat_flux_availability": (
                "available" if applied_heat_available else "inconclusive"
            ),
            "applied_pressure_work_availability": (
                "available" if applied_pressure_available else "inconclusive"
            ),
            "reconstructed_pressure_availability": (
                "available" if pressure_available else "inconclusive"
            ),
            "reconstructed_heat_flux_availability": (
                "available" if heat_proxy_available else "inconclusive"
            ),
            "applied_ledgers_signed": (
                True if applied_heat_available and applied_pressure_available else None
            ),
            "applied_heat_flux_parallel": (
                as_float(applied_heat.get("parallel"))
                if applied_heat_available else None
            ),
            "applied_heat_flux_perpendicular": (
                as_float(applied_heat.get("perpendicular"))
                if applied_heat_available else None
            ),
            "applied_heat_flux_total": (
                as_float(applied_heat.get("total")) if applied_heat_available else None
            ),
            "applied_pressure_work_total": (
                as_float(applied_pressure.get("total"))
                if applied_pressure_available else None
            ),
            "applied_pressure_work_anisotropic": (
                as_float(applied_pressure.get("anisotropic"))
                if applied_pressure_available else None
            ),
            "cap_parallel_over_1": (
                as_float(caps.get("parallel_over_1"))
                if applied_heat_available and isinstance(caps, dict) else None
            ),
            "cap_parallel_over_10": (
                as_float(caps.get("parallel_over_10"))
                if applied_heat_available and isinstance(caps, dict) else None
            ),
            "cap_perpendicular_over_1": (
                as_float(caps.get("perpendicular_over_1"))
                if applied_heat_available and isinstance(caps, dict) else None
            ),
            "cap_perpendicular_over_10": (
                as_float(caps.get("perpendicular_over_10"))
                if applied_heat_available and isinstance(caps, dict) else None
            ),
            "reconstructed_pressure_snapshot_count": (
                pressure.get("snapshot_count") if pressure_available else None
            ),
            "reconstructed_pressure_applied_to_flow": (
                pressure.get("applied_to_flow") if pressure_available else None
            ),
            "reconstructed_isotropic_perpendicular_pressure_power_mean": (
                as_float(pressure.get("isotropic_perpendicular_pressure_power_mean"))
                if pressure_available else None
            ),
            "reconstructed_anisotropic_stress_power_mean": (
                as_float(pressure.get("anisotropic_stress_power_mean"))
                if pressure_available else None
            ),
            "reconstructed_total_cgl_pressure_power_mean": (
                as_float(pressure.get("total_cgl_pressure_power_mean"))
                if pressure_available else None
            ),
            "reconstructed_anisotropic_stress_power_integral": (
                as_float(pressure_integral.get("anisotropic_stress_power_integral"))
                if isinstance(pressure_integral, dict)
                and pressure_integral.get("available") is True else None
            ),
            "reconstructed_heat_flux_snapshot_count": (
                heat_proxy.get("snapshot_count") if heat_proxy_available else None
            ),
            "reconstructed_regularized_heat_flux_power_mean": (
                as_float(heat_proxy.get("regularized_total_power_mean"))
                if heat_proxy_available else None
            ),
            "reconstructed_unlimited_heat_flux_power_mean": (
                as_float(heat_proxy.get("unlimited_total_power_mean"))
                if heat_proxy_available else None
            ),
            "reconstructed_parallel_cap_active_volume_fraction_mean": (
                as_float(heat_proxy.get("parallel_cap_active_volume_fraction_mean"))
                if heat_proxy_available else None
            ),
            "reconstructed_perpendicular_cap_active_volume_fraction_mean": (
                as_float(heat_proxy.get("perpendicular_cap_active_volume_fraction_mean"))
                if heat_proxy_available else None
            ),
            "reconstructed_regularized_heat_flux_power_integral": (
                as_float(heat_integral.get("regularized_total_power_integral"))
                if isinstance(heat_integral, dict)
                and heat_integral.get("available") is True else None
            ),
            "reconstructed_unlimited_heat_flux_power_integral": (
                as_float(heat_integral.get("unlimited_total_power_integral"))
                if isinstance(heat_integral, dict)
                and heat_integral.get("available") is True else None
            ),
            "semantics": (
                "applied columns are signed stage ledgers; reconstructed columns are "
                "sparse retained-snapshot estimates and are not applied accounting"
            ),
        })
    return rows


def mks24_panel_disposition_rows(
    data: PublicationData,
) -> list[dict[str, object]]:
    """Return admitted and explicitly blocked/external MKS24 panel dispositions."""

    rows: list[dict[str, object]] = []
    panels = (
        nested(data.science_record, "mks24.panels")
        if isinstance(data.science_record, dict)
        and data.science_record.get("_publication_evidence_validated") is True
        else None
    )
    if isinstance(panels, dict):
        for panel_id, panel in sorted(panels.items()):
            products = panel.get("products") if isinstance(panel, dict) else None
            product_records = [
                product for product in products if isinstance(product, dict)
            ] if isinstance(products, list) else []
            results = [
                validated_result(product.get("result")) for product in product_records
            ]
            rows.append({
                "panel": panel_id,
                "disposition": "admitted",
                "result": (
                    validated_result(panel.get("result"))
                    if isinstance(panel, dict) else "inconclusive"
                ),
                "product_count": len(product_records),
                "pass_count": results.count("pass"),
                "fail_count": results.count("fail"),
                "inconclusive_count": results.count("inconclusive"),
                "sources": sorted({
                    str(product.get("source"))
                    for product in product_records if product.get("source") is not None
                }),
                "reason": panel.get("reason") if isinstance(panel, dict) else None,
                "evidence": "authenticated reviewed science",
            })
    for gate in acceptance_gate_rows(data):
        name = str(gate.get("gate") or "")
        if gate.get("result") != "blocked_out_of_scope" or "panel" not in name.lower():
            continue
        rows.append({
            "panel": name.split(":", 1)[1] if ":" in name else name,
            "disposition": "blocked_or_external",
            "result": "blocked_out_of_scope",
            "product_count": None,
            "pass_count": None,
            "fail_count": None,
            "inconclusive_count": None,
            "sources": None,
            "reason": gate.get("reason"),
            "evidence": gate.get("record_type"),
        })
    if not rows:
        rows.append({
            "panel": None,
            "disposition": "inconclusive",
            "result": "inconclusive",
            "reason": "no authenticated admitted or blocked/external panel evidence",
            "evidence": "inconclusive",
        })
    return rows


def lineage_disposition_rows(data: PublicationData) -> list[dict[str, object]]:
    """Return authenticated selected and unselected fast-lineage dispositions."""

    rows: list[dict[str, object]] = []
    for case_id in CASE_IDS:
        lineage = authenticated_case_lineage(data, case_id)
        if not isinstance(lineage, dict):
            rows.append({
                "case_id": case_id,
                "lineage_index": None,
                "selected": None,
                "disposition": "inconclusive",
                "reason": "current lineage lacks authenticated binding",
                "provenance": "inconclusive",
            })
            continue
        candidates: list[tuple[bool, str, dict[str, Any], int]] = []
        selected = lineage.get("selected_fast_lineage")
        if isinstance(selected, dict):
            terminal = selected.get("terminal")
            if isinstance(terminal, dict):
                candidates.append((
                    True, str(selected.get("reason") or "selected_fast_lineage"),
                    terminal, len(selected.get("segments", []))
                    if isinstance(selected.get("segments"), list) else 0,
                ))
        elif isinstance(lineage.get("lineage"), list) and lineage["lineage"]:
            terminal = lineage["lineage"][-1]
            if isinstance(terminal, dict):
                candidates.append((
                    True, "selected_reporter_lineage", terminal, len(lineage["lineage"])
                ))
        unselected = lineage.get("unselected_lineages")
        if isinstance(unselected, list):
            for record in unselected:
                if not isinstance(record, dict):
                    continue
                terminal = record.get("terminal")
                candidates.append((
                    False, str(record.get("reason") or "unselected"),
                    terminal if isinstance(terminal, dict) else {},
                    len(record.get("segments", []))
                    if isinstance(record.get("segments"), list) else 0,
                ))
        if not candidates:
            rows.append({
                "case_id": case_id,
                "lineage_index": None,
                "selected": None,
                "disposition": "inconclusive",
                "reason": "authenticated lineage contains no lineage candidates",
                "provenance": "authenticated",
            })
            continue
        for index, (
            is_selected, reason, terminal, segment_count
        ) in enumerate(candidates):
            state = terminal.get("state")
            disposition = (
                "selected" if is_selected
                else "failed" if state == "failed"
                else "superseded" if reason == "lower_ranked_restart_linked_lineage"
                else "unselected"
            )
            rows.append({
                "case_id": case_id,
                "lineage_index": index,
                "selected": is_selected,
                "disposition": disposition,
                "reason": reason,
                "terminal_state": state,
                "source_family": terminal.get("source_family"),
                "variant": terminal.get("variant"),
                "job_id": terminal.get("job_id"),
                "observed_final_time": terminal.get("observed_final_time"),
                "run_exit_code": terminal.get("run_exit_code"),
                "restart_link_valid": terminal.get("restart_link_valid"),
                "segment_count": segment_count,
                "terminal_segment": terminal.get("segment"),
                "provenance": "authenticated",
            })
    return rows


def mechanism_metric_value(
    data: PublicationData, case_id: str, metric: str
) -> float | None:
    """Return one authenticated descriptive mechanism quantity."""

    if metric == "c_b2_full_window_mean":
        case = data.cases[case_id]
        direct = authenticated_direct_acceptance(case)
        if (
            not authenticated_science_response_case(data, case_id)
            or not isinstance(direct, dict)
            or binding_freshness_errors(
                nested(direct, "provenance.histories.user"),
                case.history_paths.get("user"),
                f"{case_id} coherent-direction user history",
            )
        ):
            return None
        series = normalized_history_series(case, "c_b2")
        return time_weighted_mean(*series) if series is not None else None
    diagnostics = authenticated_case_diagnostics(data, case_id)
    applied_paths = {
        "applied_pressure_work_total": (
            "windows.steady.lf_history.applied_pressure_work.total"
        ),
        "applied_pressure_work_anisotropic": (
            "windows.steady.lf_history.applied_pressure_work.anisotropic"
        ),
    }
    if metric in applied_paths:
        return as_float(nested(diagnostics, applied_paths[metric]))
    ensemble = authenticated_snapshot_ensemble(data, case_id)
    reconstructed_paths = {
        "reconstructed_anisotropic_stress_power_mean": (
            "pressure_work_decomposition.anisotropic_stress_power_mean"
        ),
        "parallel_strain_rms_mean": (
            "pressure_work_decomposition.parallel_strain_rms_mean"
        ),
    }
    if metric in reconstructed_paths:
        return as_float(nested(ensemble, reconstructed_paths[metric]))
    if metric == "reviewed_abs_dp_standardized_effect":
        effects = [
            row for active, passive in ACTIVE_PASSIVE_PAIRS
            if active == case_id
            for row in reviewed_pair_effect_rows(data, active, passive)
            if row.get("metric") == "abs_dp"
        ]
        return (
            as_float(effects[0].get("standardized_effect")) if len(effects) == 1
            else None
        )
    raise PublicationError(f"unsupported mechanism metric: {metric}")


def coherent_direction_mechanism_rows(
    data: PublicationData,
) -> list[dict[str, object]]:
    """Summarize descriptive active-minus-passive directions without a pass gate."""

    metrics = (
        "c_b2_full_window_mean",
        "applied_pressure_work_total",
        "applied_pressure_work_anisotropic",
        "reconstructed_anisotropic_stress_power_mean",
        "parallel_strain_rms_mean",
        "reviewed_abs_dp_standardized_effect",
    )
    rows: list[dict[str, object]] = []
    for metric in metrics:
        differences: dict[str, float | None] = {}
        for active, passive in ACTIVE_PASSIVE_PAIRS:
            if metric == "reviewed_abs_dp_standardized_effect":
                active_value = mechanism_metric_value(data, active, metric)
                passive_value = 0.0 if active_value is not None else None
            else:
                active_value = mechanism_metric_value(data, active, metric)
                passive_value = mechanism_metric_value(data, passive, metric)
            differences[f"{active}_{passive}"] = (
                active_value - passive_value
                if active_value is not None and passive_value is not None else None
            )
        available = [value for value in differences.values() if value is not None]
        positive = sum(value > 0.0 for value in available)
        negative = sum(value < 0.0 for value in available)
        equal = sum(value == 0.0 for value in available)
        direction = (
            "active_gt_passive"
            if len(available) == len(ACTIVE_PASSIVE_PAIRS) and positive == len(available)
            else "active_lt_passive"
            if len(available) == len(ACTIVE_PASSIVE_PAIRS) and negative == len(available)
            else "equal"
            if len(available) == len(ACTIVE_PASSIVE_PAIRS) and equal == len(available)
            else "mixed"
            if len(available) == len(ACTIVE_PASSIVE_PAIRS)
            else "inconclusive"
        )
        rows.append({
            "metric": metric,
            "available_pair_count": len(available),
            "positive_active_minus_passive_count": positive,
            "negative_active_minus_passive_count": negative,
            "equal_count": equal,
            "inconclusive_pair_count": len(ACTIVE_PASSIVE_PAIRS) - len(available),
            "descriptive_direction": direction,
            "pair_active_minus_passive": differences,
            "inference_scope": "descriptive_only_no_preregistered_pass_gate",
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


def authenticated_science_response_case(
    data: PublicationData, case_id: str
) -> bool:
    """Return whether one case may populate authenticated response figures."""

    case = data.cases[case_id]
    direct = authenticated_direct_acceptance(case)
    return (
        direct is not None
        and direct.get("result") == "pass"
        and scientific_response_eligible(data, case)
        and science_case_status(data, case_id) == "pass"
    )


def finite_curve(
    record: object, x_name: str, y_name: str
) -> tuple[list[float], list[float]] | None:
    """Return one finite, strictly ordered curve from an evidence record."""

    if not isinstance(record, dict):
        return None
    x_values = record.get(x_name)
    y_values = record.get(y_name)
    if (
        not isinstance(x_values, list)
        or not isinstance(y_values, list)
        or len(x_values) != len(y_values)
        or len(x_values) < 2
    ):
        return None
    x = [as_float(value) for value in x_values]
    y = [as_float(value) for value in y_values]
    if any(value is None for value in (*x, *y)):
        return None
    finite_x = [float(value) for value in x if value is not None]
    finite_y = [float(value) for value in y if value is not None]
    if any(right <= left for left, right in zip(finite_x, finite_x[1:])):
        return None
    return finite_x, finite_y


def histogram_curve(record: object) -> tuple[list[float], list[float]] | None:
    """Return finite histogram centers and nonnegative density."""

    if not isinstance(record, dict):
        return None
    edges = record.get("edges")
    density = record.get("density")
    if (
        not isinstance(edges, list)
        or not isinstance(density, list)
        or len(edges) != len(density) + 1
        or len(density) < 2
    ):
        return None
    parsed_edges = [as_float(value) for value in edges]
    parsed_density = [as_float(value) for value in density]
    if any(value is None for value in (*parsed_edges, *parsed_density)):
        return None
    edge_values = [float(value) for value in parsed_edges if value is not None]
    density_values = [float(value) for value in parsed_density if value is not None]
    if (
        any(right <= left for left, right in zip(edge_values, edge_values[1:]))
        or any(value < 0.0 for value in density_values)
    ):
        return None
    return (
        [0.5 * (left + right) for left, right in zip(edge_values, edge_values[1:])],
        density_values,
    )


def authenticated_snapshot_ensemble(
    data: PublicationData, case_id: str
) -> dict[str, Any] | None:
    """Return one authenticated complete snapshot ensemble."""

    diagnostics = authenticated_case_diagnostics(data, case_id)
    ensemble = diagnostics.get("snapshot_ensemble") if isinstance(diagnostics, dict) else None
    if (
        not isinstance(ensemble, dict)
        or diagnostics.get("snapshot_analysis_status") != "complete"
        or not isinstance(ensemble.get("snapshot_count"), int)
        or ensemble.get("snapshot_count", 0) <= 0
    ):
        return None
    return ensemble


def reviewed_pair_effect_rows(
    data: PublicationData, active: str, passive: str
) -> list[dict[str, object]]:
    """Return authenticated standardized active/passive effects for one matched pair."""

    return [
        row
        for row in science_contrast_rows(data)
        if row.get("family") == "active_passive"
        and row.get("left") == active
        and row.get("right") == passive
        and row.get("claim_eligible") is True
        and row.get("available") is True
        and as_float(row.get("standardized_effect")) is not None
    ]


def render_causal_mechanism(data: PublicationData, plt: Any, path: Path) -> None:
    """Render matched C_B2 histories, strain PDFs, and reviewed causal effects."""

    fig, axes = plt.subplots(4, 3, figsize=(11.2, 10.6))
    for row_index, (active, passive) in enumerate(ACTIVE_PASSIVE_PAIRS):
        history_axis, strain_axis, effect_axis = axes[row_index]
        history_count = 0
        strain_count = 0
        for case_id, style, label in (
            (active, "-", f"{active} active"),
            (passive, "--", f"{passive} passive"),
        ):
            if not authenticated_science_response_case(data, case_id):
                continue
            history = normalized_history_series(data.cases[case_id], "c_b2")
            if history is not None:
                history_axis.plot(
                    history[0], history[1], style, color=CASE_COLORS[case_id], label=label
                )
                history_count += 1
            ensemble = authenticated_snapshot_ensemble(data, case_id)
            distribution = (
                histogram_curve(nested(ensemble, "pdf.bb_grad_velocity"))
                if isinstance(ensemble, dict) else None
            )
            if distribution is not None:
                strain_axis.plot(
                    distribution[0], distribution[1], style,
                    color=CASE_COLORS[case_id], label=label,
                )
                strain_count += 1
        history_axis.axvspan(4.0, 10.0, color="#eeeeee", alpha=0.5, zorder=-10)
        history_axis.set_xlim(0.0, 10.0)
        history_axis.set_ylabel(r"$C_{B^2}$")
        history_axis.grid(True, alpha=0.25)
        if history_count:
            history_axis.legend(frameon=False)
        else:
            history_axis.text(
                0.5, 0.5, "inconclusive: authenticated pair history unavailable",
                transform=history_axis.transAxes, ha="center", va="center",
                color="#666666",
            )

        strain_axis.set_ylabel(r"PDF of $\hat{b}\hat{b}:\nabla u$")
        strain_axis.grid(True, alpha=0.25)
        if strain_count:
            strain_axis.legend(frameon=False)
        else:
            strain_axis.text(
                0.5, 0.5, "inconclusive: authenticated strain PDFs unavailable",
                transform=strain_axis.transAxes, ha="center", va="center",
                color="#666666",
            )

        effects = reviewed_pair_effect_rows(data, active, passive)
        if effects:
            effect_axis.barh(
                list(range(len(effects))),
                [float(row["standardized_effect"]) for row in effects],
                color=CASE_COLORS[active],
                edgecolor="black",
                linewidth=0.35,
            )
            effect_axis.set_yticks(
                list(range(len(effects))),
                [str(row.get("metric")) for row in effects],
            )
            effect_axis.axvline(0.0, color="black", linewidth=0.8, linestyle=":")
        else:
            effect_axis.text(
                0.5, 0.5, "inconclusive: reviewed effect evidence unavailable",
                transform=effect_axis.transAxes, ha="center", va="center",
                color="#666666",
            )
        effect_axis.set_xlabel("standardized active-minus-passive effect")
        effect_axis.grid(True, axis="x", alpha=0.25)
        history_axis.set_title(f"{active} / {passive}: $C_{{B^2}}(t)$")
        strain_axis.set_title(f"{active} / {passive}: parallel strain")
        effect_axis.set_title(f"{active} / {passive}: reviewed effects")
    for axis in axes[-1, :2]:
        axis.set_xlabel(r"$t/(L_\perp/v_A)$" if axis is axes[-1, 0] else "strain")
    fig.suptitle(
        "Matched active/passive causal-mechanism evidence "
        "(authenticated complete-case products only)",
        y=0.997,
    )
    fig.text(
        0.5, 0.006,
        "Effects are descriptive within-realization active-minus-passive contrasts; "
        "missing evidence is reported as inconclusive and is never plotted as zero.",
        ha="center", va="bottom", fontsize=6.8,
    )
    fig.tight_layout(rect=(0.0, 0.02, 1.0, 0.985))
    save_figure(fig, path)
    plt.close(fig)


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


def resolution_common_range(data: PublicationData) -> tuple[float, float] | None:
    """Return the authenticated preregistered common k_perp/pi interval."""

    if (
        not isinstance(data.science_record, dict)
        or data.science_record.get("_publication_evidence_validated") is not True
    ):
        return None
    values = nested(data.science_record, "resolution.limits.common_k_perp_over_pi")
    if not isinstance(values, list) or len(values) != 2:
        return None
    low, high = as_float(values[0]), as_float(values[1])
    if low is None or high is None or low <= 0.0 or high <= low:
        return None
    return low, high


def normalized_resolution_spectrum(
    data: PublicationData, case_id: str, product: str
) -> tuple[list[float], list[float]] | None:
    """Return one authenticated spectrum normalized over the common range."""

    common_range = resolution_common_range(data)
    ensemble = authenticated_snapshot_ensemble(data, case_id)
    curve = (
        finite_curve(nested(ensemble, f"spectra.{product}"), "k", "power_per_dk")
        if isinstance(ensemble, dict) else None
    )
    if common_range is None or curve is None:
        return None
    low, high = common_range
    selected = [
        (x / math.pi, y)
        for x, y in zip(*curve)
        if low - TIME_TOLERANCE <= x / math.pi <= high + TIME_TOLERANCE and y > 0.0
    ]
    if len(selected) < 2:
        return None
    x_values = [value[0] for value in selected]
    y_values = [value[1] for value in selected]
    if x_values[0] > low + TIME_TOLERANCE or x_values[-1] < high - TIME_TOLERANCE:
        return None
    integral = sum(
        0.5 * (left_y + right_y) * (right_x - left_x)
        for left_x, right_x, left_y, right_y in zip(
            x_values, x_values[1:], y_values, y_values[1:]
        )
    )
    if not math.isfinite(integral) or integral <= 0.0:
        return None
    return x_values, [value / integral for value in y_values]


def resolution_alignment_curve(
    data: PublicationData, case_id: str
) -> tuple[list[float], list[float]] | None:
    """Return authenticated peak alignment on the preregistered common range."""

    common_range = resolution_common_range(data)
    ensemble = authenticated_snapshot_ensemble(data, case_id)
    alignment = ensemble.get("alignment") if isinstance(ensemble, dict) else None
    dk = as_float(nested(ensemble, "spectra.velocity.dk"))
    if common_range is None or not isinstance(alignment, dict) or dk is None or dk <= 0.0:
        return None
    peaks: list[tuple[float, float]] = []
    try:
        ordered = sorted(alignment.items(), key=lambda item: int(item[0]))
    except (TypeError, ValueError):
        return None
    for shell, record in ordered:
        curve = histogram_curve(record)
        if curve is None:
            return None
        x = int(shell) * dk / math.pi
        if common_range[0] - TIME_TOLERANCE <= x <= common_range[1] + TIME_TOLERANCE:
            maximum = max(range(len(curve[1])), key=lambda index: curve[1][index])
            peaks.append((x, curve[0][maximum]))
    if len(peaks) < 2:
        return None
    return [value[0] for value in peaks], [value[1] for value in peaks]


def render_resolution_curves(data: PublicationData, plt: Any, path: Path) -> None:
    """Render actual authenticated R16/R02/R17 spectra and alignment criteria."""

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.4))
    common_range = resolution_common_range(data)
    resolution = (
        data.science_record.get("resolution")
        if isinstance(data.science_record, dict) else None
    )
    result = (
        validated_result(resolution.get("result"))
        if isinstance(resolution, dict) else "inconclusive"
    )
    products = (
        ("velocity", "Velocity spectrum", axes[0, 0]),
        ("magnetic_fluctuation", "Magnetic-fluctuation spectrum", axes[0, 1]),
    )
    for product, title, axis in products:
        curves = {
            case_id: normalized_resolution_spectrum(data, case_id, product)
            for case_id in RESOLUTION_CASES
            if authenticated_science_response_case(data, case_id)
        }
        if len(curves) == len(RESOLUTION_CASES) and all(
            curve is not None for curve in curves.values()
        ):
            for case_id in RESOLUTION_CASES:
                curve = curves[case_id]
                assert curve is not None
                axis.loglog(
                    curve[0], curve[1], color=CASE_COLORS[case_id], label=case_id
                )
            axis.legend(frameon=False)
        else:
            axis.text(
                0.5, 0.5,
                "inconclusive: all three authenticated common-range curves required",
                transform=axis.transAxes, ha="center", va="center", color="#666666",
            )
        axis.set_xlabel(r"$k_\perp/\pi$")
        axis.set_ylabel("common-range normalized power")
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.25)

    alignment_axis = axes[1, 0]
    alignment = {
        case_id: resolution_alignment_curve(data, case_id)
        for case_id in RESOLUTION_CASES
        if authenticated_science_response_case(data, case_id)
    }
    if len(alignment) == len(RESOLUTION_CASES) and all(
        curve is not None for curve in alignment.values()
    ):
        for case_id in RESOLUTION_CASES:
            curve = alignment[case_id]
            assert curve is not None
            alignment_axis.plot(
                curve[0], curve[1], marker="o", markersize=3.0,
                color=CASE_COLORS[case_id], label=case_id,
            )
        alignment_axis.legend(frameon=False)
    else:
        alignment_axis.text(
            0.5, 0.5,
            "inconclusive: all three authenticated alignment curves required",
            transform=alignment_axis.transAxes, ha="center", va="center",
            color="#666666",
        )
    alignment_axis.set_xlabel(r"$k_\perp/\pi$")
    alignment_axis.set_ylabel(r"peak $|\cos\theta|$")
    alignment_axis.set_title("Peak velocity/magnetic alignment")
    alignment_axis.grid(True, alpha=0.25)

    criteria_axis = axes[1, 1]
    criteria_axis.axis("off")
    observations = (
        resolution.get("observations") if isinstance(resolution, dict) else None
    )
    lines = [
        f"Reviewed convergence result: {result}",
        (
            f"Preregistered common range: {common_range[0]:g} <= k_perp/pi "
            f"<= {common_range[1]:g}"
            if common_range is not None else
            "Preregistered common range: inconclusive/unavailable"
        ),
        "",
    ]
    if isinstance(observations, list) and observations:
        for observation in observations:
            if not isinstance(observation, dict):
                continue
            name = observation.get("product", observation.get("metric", "--"))
            if observation.get("available") is not True:
                lines.append(f"{name}: inconclusive")
                continue
            distance = as_float(observation.get("R02_R17_distance"))
            limit = as_float(observation.get("R02_R17_limit"))
            relative = as_float(observation.get("relative_R02_R17_difference"))
            measured = distance if distance is not None else relative
            passed = observation.get("passed")
            decision = (
                "pass" if passed is True
                else "fail" if passed is False
                else "inconclusive"
            )
            lines.append(
                f"{name}: {decision}; "
                f"R02/R17={text_value(measured)}"
                + (f"; limit={limit:.4g}" if limit is not None else "")
                + (
                    f"; improved={text_value(observation.get('improved'))}"
                    if "improved" in observation else ""
                )
            )
    else:
        lines.append("Reviewed preregistered criteria: inconclusive/unavailable")
    criteria_axis.text(
        0.02, 0.98, "\n".join(lines), transform=criteria_axis.transAxes,
        ha="left", va="top", fontsize=7.2,
    )
    fig.suptitle(
        "R16 / R02 / R17 authenticated common-range resolution evidence", y=0.995
    )
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
            result = record.get("_publication_hyperbolicity_result")
            if isinstance(result, dict):
                roots.append(result)
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
    """Return normalized retained-state directional-discriminant diagnostics."""

    roots = case_evidence_roots(data, case_id)
    snapshot_records: list[dict[str, object]] = []
    for root in roots:
        summary = snapshot_hyperbolicity_summary(root.get("snapshots"))
        if summary is not None:
            snapshot_records.append(summary)
    if snapshot_records:
        severity = {
            "nonnegative_discriminant": 0, "hyperbolic": 0,
            "negative": 1, "nonfinite": 2,
        }
        return max(
            snapshot_records,
            key=lambda record: (
                severity.get(str(record["result"]), -1),
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
    normalized_status = (
        "nonnegative_discriminant" if status in {"pass", "hyperbolic"}
        else str(status)
        if status in {
            "nonnegative_discriminant", "negative", "nonfinite", "inconclusive",
        }
        else None
    )
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
    negative_count = exact_nonnegative_integer(first_evidence_value(
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
    evaluation_count = exact_nonnegative_integer(first_evidence_value(
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
    nonfinite_count = exact_nonnegative_integer(first_evidence_value(
        roots,
        (
            "hyperbolicity.nonfinite_discriminant_count",
            "retained_state_hyperbolicity.nonfinite_discriminant_count",
            "summary.nonfinite_discriminant_count",
            "nonfinite_discriminant_count",
        ),
    ))
    if nonfinite_count is not None and nonfinite_count > 0:
        normalized_status = "nonfinite"
    elif (
            (negative_fraction is not None and negative_fraction > 0.0)
            or (negative_count is not None and negative_count > 0.0)
            or (minimum is not None and minimum < 0.0)
    ):
        normalized_status = "negative"
    elif (
        normalized_status is None
        and negative_count == 0
        and nonfinite_count in {None, 0}
        and evaluation_count is not None
        and evaluation_count > 0
    ):
        normalized_status = "nonnegative_discriminant"
    elif normalized_status is None:
        normalized_status = "inconclusive"
    return {
        "result": normalized_status,
        "negative_discriminant_fraction": negative_fraction,
        "negative_discriminant_count": negative_count,
        "cell_direction_evaluations": evaluation_count,
        "minimum_discriminant": minimum,
        "nonfinite_discriminant_count": nonfinite_count,
    }


def hyperbolicity_status(data: PublicationData, case_id: str) -> str:
    """Return the normalized retained-state directional-discriminant disposition."""

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
            "retained_state_coordinate_discriminant": hyper["result"],
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
            "Retained coordinate\nD coverage", "Observed\ndiagnostic", "Claim scope",
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


def render_hyperbolicity_coverage(
    data: PublicationData, plt: Any, colors: Any, patches: Any, path: Path
) -> None:
    """Render compact retained coordinate-direction discriminant coverage."""

    rows = [
        row for row in hyperbolicity_coverage_rows(data)
        if row["case_id"] in ACTIVE_ENERGY_CASES
    ]
    statuses: list[list[str]] = []
    labels: list[list[str]] = []
    for row in rows:
        coverage = str(row["coverage_result"])
        formula_id = row.get("formula_id")
        formula_family = row.get("formula_disposition_family")
        formula_status = (
            "warning" if formula_family == "legacy_implementation"
            else "configuration" if formula_family == "literature_correct"
            else "inconclusive"
        )
        legacy = str(row["legacy_implementation_disposition"])
        literature = str(row["literature_correct_disposition"])
        experiment_scope = str(row["experiment_scope"])
        experiment_status = (
            "restricted" if experiment_scope == "restricted" else "configuration"
        )
        claim = str(row["strict_hyperbolic_claim_status"])
        claim_label = {
            "eligible": "eligible",
            "inconclusive_legacy_implementation": "legacy-only",
            "inconclusive_formula_identity": "formula unknown",
            "inconclusive_formula_provenance": "formula provenance?",
            "inconclusive_formula_executable_compatibility": "compatibility?",
            "excluded_formula_executable_incompatible": "incompatible",
            "excluded_negative": "negative",
            "excluded_nonfinite": "nonfinite",
            "excluded_experiment_scope": "restricted",
            "inconclusive_retained_coordinate_discriminant_scope": (
                "retained coordinate directions only"
            ),
            "not_applicable": "not applicable",
            "inconclusive": "inconclusive",
        }.get(claim, claim)
        claim_status = (
            "pass" if claim == "eligible"
            else "fail" if claim in {
                "excluded_negative",
                "excluded_nonfinite",
                "excluded_formula_executable_incompatible",
            }
            else "restricted" if claim == "excluded_experiment_scope"
            else "configuration" if claim == "not_applicable"
            else "inconclusive"
        )
        selected = row.get("selected_snapshot_count")
        complete = row.get("complete_retained_snapshot_count")
        audited = row.get("audited_snapshot_count")
        statuses.append([
            coverage, formula_status, legacy, literature, experiment_status,
            claim_status,
        ])
        labels.append([
            (
                f"{text_value(audited)}/{text_value(complete)} audited"
                if audited is not None or complete is not None else "inconclusive"
            ),
            text_value(formula_id),
            display_status(legacy),
            display_status(literature),
            experiment_scope,
            claim_label,
        ])
        if selected is not None and audited is None:
            labels[-1][0] = f"{text_value(selected)}/{text_value(complete)} selected"
    status_matrix_figure(
        plt, colors, patches,
        [str(row["case_id"]) for row in rows],
        [
            "All retained snapshots", "Audit formula",
            "Legacy implementation", "Literature correct",
            "Experiment scope", "Strict-hyperbolic claim",
        ],
        statuses,
        labels,
        "Authenticated retained-snapshot directional-discriminant coverage",
        path,
        (
            "Coverage passes only when an authenticated snapshot_policy=all manifest "
            "and completed result cover every complete retained snapshot exactly once. "
            "Dispositions are assigned only to the authenticated audit formula. "
            "Nonnegative values cover retained cell-centered states in three "
            "coordinate-normal directions only; they do not establish face, "
            "intermediate, oblique, full, or strict hyperbolicity and do not exclude "
            "sqrt(|D|) fallback. Experiment scope remains separate; passive-delta "
            "cases are not applicable."
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


def canonical_product_records(
    products: Iterable[Path], staging: Path, output: Path
) -> list[tuple[Path, Path, dict[str, object]]]:
    """Bind staged products to their final canonical publication paths."""

    staging = staging.absolute()
    output = output.absolute()
    if staging.is_symlink() or not staging.is_dir():
        raise PublicationError(
            f"publication staging root is missing, invalid, or a symlink: {staging}"
        )
    resolved_staging = staging.resolve(strict=True)
    records: list[tuple[Path, Path, dict[str, object]]] = []
    seen: set[Path] = set()
    for staged in sorted(path.absolute() for path in products):
        try:
            relative = staged.relative_to(staging)
        except ValueError as error:
            raise PublicationError(
                f"staged publication product escapes staging root: {staged}"
            ) from error
        if not relative.parts or ".." in relative.parts:
            raise PublicationError(
                f"staged publication product escapes staging root: {staged}"
            )
        if (
            relative.parts[0] in {
                "manifest.json", PUBLICATION_AUTHORITY_TOKEN_NAME,
            }
            or relative.parts[0].startswith(
                PUBLICATION_MANIFEST_QUARANTINE_PREFIX
            )
        ):
            raise PublicationError(
                "reserved publication authority paths may not be promoted as products "
                "or product directories"
            )
        if staged.is_symlink() or not staged.is_file():
            raise PublicationError(
                f"staged publication product is missing, invalid, or a symlink: {staged}"
            )
        try:
            resolved_staged = staged.resolve(strict=True)
            resolved_staged.relative_to(resolved_staging)
        except (OSError, ValueError) as error:
            raise PublicationError(
                f"staged publication product escapes staging root: {staged}"
            ) from error
        canonical = output / relative
        try:
            canonical.relative_to(output)
        except ValueError as error:
            raise PublicationError(
                f"canonical publication product escapes output root: {canonical}"
            ) from error
        if canonical in seen:
            raise PublicationError(f"duplicate canonical publication product: {canonical}")
        seen.add(canonical)
        binding = source_binding(staged)
        binding["path"] = str(canonical)
        records.append((staged, canonical, binding))
    if not records:
        raise PublicationError("publication renderer produced no products")
    return records


def canonical_product_bindings(
    products: Iterable[Path], staging: Path, output: Path
) -> list[dict[str, object]]:
    """Return deterministic bindings for staged bytes at canonical paths."""

    return [
        binding
        for _, _, binding in canonical_product_records(products, staging, output)
    ]


def fsync_regular_file(path: Path) -> None:
    """Durably flush one no-follow regular file."""

    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise PublicationError(
            f"cannot open staged publication regular file for fsync: {path}"
        ) from error
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise PublicationError(
                f"staged publication path is not a regular file: {path}"
            )
        os.fsync(descriptor)
    except OSError as error:
        raise PublicationError(
            f"cannot fsync staged publication regular file: {path}"
        ) from error
    finally:
        os.close(descriptor)


def directory_open_flags() -> int:
    """Return no-follow flags for a directory capability."""

    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    return flags


def open_output_directory(output: Path) -> int:
    """Open one canonical output root without following its final component."""

    try:
        return os.open(output, directory_open_flags())
    except OSError as error:
        raise PublicationError(
            f"cannot open publication output directory without following symlinks: "
            f"{output}"
        ) from error


def open_output_directory_at(parent_descriptor: int, output_name: str) -> int:
    """Open one output root relative to its held, locked parent capability."""

    try:
        return os.open(output_name, directory_open_flags(), dir_fd=parent_descriptor)
    except OSError as error:
        raise PublicationError(
            "cannot open publication output directory relative to its held parent: "
            f"{output_name}"
        ) from error


@contextmanager
def relative_parent_directory(
    output_descriptor: int, relative: Path, *, create: bool
) -> Iterator[tuple[int, str]]:
    """Yield the no-follow parent capability and basename for one relative path."""

    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise PublicationError(
            f"canonical publication product escapes output root: {relative}"
        )
    current = os.dup(output_descriptor)
    try:
        for part in relative.parent.parts:
            if part in {"", "."}:
                continue
            if create:
                try:
                    os.mkdir(part, mode=0o755, dir_fd=current)
                    os.fsync(current)
                except FileExistsError:
                    pass
                except OSError as error:
                    raise PublicationError(
                        f"cannot create canonical publication directory: "
                        f"{relative.parent}"
                    ) from error
            try:
                child = os.open(part, directory_open_flags(), dir_fd=current)
            except OSError as error:
                raise PublicationError(
                    f"canonical publication parent is missing, invalid, or a symlink: "
                    f"{relative.parent}"
                ) from error
            os.close(current)
            current = child
        yield current, relative.name
    finally:
        os.close(current)


def promote_file(staged: Path, relative: Path, output_descriptor: int) -> None:
    """Atomically replace one product relative to the held output capability."""

    with relative_parent_directory(
        output_descriptor, relative, create=True
    ) as (parent_descriptor, name):
        try:
            os.replace(staged, name, dst_dir_fd=parent_descriptor)
            os.fsync(parent_descriptor)
        except OSError as error:
            raise PublicationError(
                f"cannot promote canonical publication product: {relative}"
            ) from error


def write_durable_regular_file_at(
    parent_descriptor: int, name: str, payload: bytes
) -> None:
    """Durably replace one regular file through a held parent capability."""

    if not name or name in {".", ".."} or Path(name).name != name:
        raise PublicationError(f"invalid publication sibling filename: {name}")
    staged_name = f".{name}.staging-{secrets.token_hex(32)}"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(staged_name, flags, 0o600, dir_fd=parent_descriptor)
    except OSError as error:
        raise PublicationError(
            f"cannot create staged publication sibling file: {staged_name}"
        ) from error
    try:
        offset = 0
        while offset < len(payload):
            offset += os.write(descriptor, payload[offset:])
        os.fsync(descriptor)
    except OSError as error:
        raise PublicationError(
            f"cannot durably write staged publication sibling file: {staged_name}"
        ) from error
    finally:
        os.close(descriptor)
    try:
        os.replace(
            staged_name,
            name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        os.fsync(parent_descriptor)
    except OSError as error:
        # Failures intentionally preserve the unique staged path. Removing it after
        # a pathname race could delete bytes not created by this publisher.
        raise PublicationError(
            f"cannot promote durable publication sibling file: {name}"
        ) from error


def publication_lock_path(output: Path) -> Path:
    """Return the persistent sibling lock serializing canonical publication."""

    return output.parent / f".{output.name}.publication.lock"


def publication_ownership_path(output: Path) -> Path:
    """Return the persistent sibling ownership record for canonical output."""

    return output.parent / f".{output.name}.publication-owner.json"


def normalized_filesystem_path(path: Path) -> Path:
    """Return one absolute path with existing symlinked ancestors resolved."""

    return path.expanduser().resolve(strict=False)


def path_contains(root: Path, path: Path) -> bool:
    """Return whether one normalized path equals or contains another."""

    normalized_root = normalized_filesystem_path(root)
    normalized_path = normalized_filesystem_path(path)
    return (
        normalized_path == normalized_root
        or normalized_path.is_relative_to(normalized_root)
    )


def validate_publication_output_overlap(
    analysis: Path,
    output: Path,
    evidence_paths: Iterable[Path],
    acceptance_roots: Iterable[Path] = (),
) -> None:
    """Reject output roots that could delete analysis or discovered evidence."""

    if path_contains(output, analysis):
        raise PublicationError(
            f"publication output may not equal or contain the analysis root: {output}"
        )
    overlapping_roots = sorted({
        normalized_filesystem_path(path)
        for path in acceptance_roots
        if path_contains(path, output) or path_contains(output, path)
    })
    if overlapping_roots:
        raise PublicationError(
            "publication output overlaps an acceptance root: "
            + ", ".join(str(path) for path in overlapping_roots)
        )
    overlapping = sorted({
        normalized_filesystem_path(path)
        for path in evidence_paths
        if path_contains(output, path)
    })
    if overlapping:
        raise PublicationError(
            "publication output contains discovered source or acceptance evidence: "
            + ", ".join(str(path) for path in overlapping)
        )


def read_regular_file_at(output_descriptor: int, relative: Path) -> bytes:
    """Read one regular file relative to a held output capability."""

    with relative_parent_directory(
        output_descriptor, relative, create=False
    ) as (parent_descriptor, name):
        flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            descriptor = os.open(name, flags, dir_fd=parent_descriptor)
        except OSError as error:
            raise PublicationError(
                f"cannot open canonical publication regular file: {relative}"
            ) from error
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise PublicationError(
                    f"canonical publication path is not a regular file: {relative}"
                )
            chunks = []
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    return b"".join(chunks)
                chunks.append(chunk)
        finally:
            os.close(descriptor)


def fsync_exact_regular_file_at(
    output_descriptor: int, relative: Path, expected_payload: bytes
) -> None:
    """Require and durably flush exact bytes through a held output capability."""

    with relative_parent_directory(
        output_descriptor, relative, create=False
    ) as (parent_descriptor, name):
        flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            descriptor = os.open(name, flags, dir_fd=parent_descriptor)
        except OSError as error:
            raise PublicationError(
                f"cannot open exact publication regular file for fsync: {relative}"
            ) from error
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise PublicationError(
                    f"exact publication path is not a regular file: {relative}"
                )
            chunks: list[bytes] = []
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                chunks.append(chunk)
            if not secrets.compare_digest(b"".join(chunks), expected_payload):
                raise PublicationError(
                    f"publication regular file differs from exact expected bytes: "
                    f"{relative}"
                )
            os.fsync(descriptor)
        except OSError as error:
            raise PublicationError(
                f"cannot fsync exact publication regular file: {relative}"
            ) from error
        finally:
            os.close(descriptor)
        try:
            observed = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
        except OSError as error:
            raise PublicationError(
                f"cannot revalidate exact publication regular file: {relative}"
            ) from error
        if not stat.S_ISREG(observed.st_mode):
            raise PublicationError(
                f"exact publication path changed during fsync: {relative}"
            )
    if not secrets.compare_digest(
        read_regular_file_at(output_descriptor, relative), expected_payload
    ):
        raise PublicationError(
            f"publication regular file changed after exact fsync: {relative}"
        )


def publication_authority_token(output_descriptor: int) -> bytes | None:
    """Return and validate the persistent random token inside output."""

    try:
        token = read_regular_file_at(
            output_descriptor, Path(PUBLICATION_AUTHORITY_TOKEN_NAME)
        )
    except PublicationError as error:
        try:
            os.stat(
                PUBLICATION_AUTHORITY_TOKEN_NAME,
                dir_fd=output_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            return None
        raise error
    try:
        text = token.decode("ascii")
    except UnicodeDecodeError as error:
        raise PublicationError("publication authority token is malformed") from error
    if (
        len(token) != 65
        or not text.endswith("\n")
        or SHA256_PATTERN.fullmatch(text[:-1]) is None
    ):
        raise PublicationError("publication authority token is malformed")
    return token


def create_publication_authority_token(output_descriptor: int) -> bytes:
    """Create one cryptographically random persistent token inside output."""

    token = (secrets.token_hex(32) + "\n").encode("ascii")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(
            PUBLICATION_AUTHORITY_TOKEN_NAME,
            flags,
            0o600,
            dir_fd=output_descriptor,
        )
    except OSError as error:
        raise PublicationError("cannot create publication authority token") from error
    try:
        offset = 0
        while offset < len(token):
            offset += os.write(descriptor, token[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.fsync(output_descriptor)
    return token


def publication_ownership_record(
    output: Path, token: bytes
) -> dict[str, object]:
    """Return token-bound persistent ownership proof for one canonical output."""

    return {
        "schema_version": 3,
        "record_type": PUBLICATION_OWNERSHIP_RECORD_TYPE,
        "output": str(normalized_filesystem_path(output)),
        "token_file": PUBLICATION_AUTHORITY_TOKEN_NAME,
        "token_sha256": hashlib.sha256(token).hexdigest(),
        "authority": "cgl_lf_stage_i_fast_publication.py",
    }


def load_publication_ownership(output: Path) -> dict[str, Any] | None:
    """Load a regular sibling ownership record, if present."""

    path = publication_ownership_path(output)
    if path.is_symlink() or not path.is_file():
        return None
    try:
        return load_json(path)
    except (OSError, json.JSONDecodeError, PublicationError):
        return None


def load_publication_ownership_at(
    output: Path, parent_descriptor: int
) -> dict[str, Any] | None:
    """Load a regular sibling ownership record through a held parent capability."""

    name = publication_ownership_path(output).name
    try:
        observed = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return None
    except OSError as error:
        raise PublicationError("cannot inspect publication ownership record") from error
    if not stat.S_ISREG(observed.st_mode):
        return None
    try:
        value = json.loads(read_regular_file_at(parent_descriptor, Path(name)))
    except (UnicodeDecodeError, json.JSONDecodeError, PublicationError):
        return None
    return value if isinstance(value, dict) else None


def valid_publication_ownership_at(
    output: Path,
    output_descriptor: int,
    token: bytes | None = None,
    parent_descriptor: int | None = None,
) -> bool:
    """Return whether the sibling record binds the token in held output."""

    try:
        observed_token = (
            token
            if token is not None
            else publication_authority_token(output_descriptor)
        )
    except PublicationError:
        return False
    return (
        observed_token is not None
        and (
            load_publication_ownership_at(output, parent_descriptor)
            if parent_descriptor is not None
            else load_publication_ownership(output)
        )
        == publication_ownership_record(output, observed_token)
    )


def valid_publication_ownership(output: Path) -> bool:
    """Return whether token-backed persistent ownership binds canonical output."""

    if output.is_symlink() or not output.is_dir():
        return False
    try:
        descriptor = open_output_directory(output)
    except PublicationError:
        return False
    try:
        return valid_publication_ownership_at(output, descriptor)
    finally:
        os.close(descriptor)


def fsync_directory(path: Path) -> None:
    """Durably publish prior directory-entry changes."""

    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise PublicationError(f"cannot open directory for fsync: {path}") from error
    try:
        os.fsync(descriptor)
    except OSError as error:
        raise PublicationError(f"cannot fsync directory: {path}") from error
    finally:
        os.close(descriptor)


def canonical_path_resolves_to_held_output(
    output: Path,
    output_descriptor: int,
    token: bytes | None = None,
    *,
    require_ownership: bool,
    parent_descriptor: int | None = None,
) -> bool:
    """Return whether canonical path still names the held token-bearing instance."""

    try:
        observed_descriptor = open_output_directory(output)
    except PublicationError:
        return False
    try:
        # This is only a live held-descriptor race check. Persistent authority is
        # the random token, never a stored device/inode tuple.
        if not os.path.samestat(
            os.fstat(output_descriptor), os.fstat(observed_descriptor)
        ):
            return False
        if token is not None and publication_authority_token(
            observed_descriptor
        ) != token:
            return False
        return (
            not require_ownership
            or token is not None
            and valid_publication_ownership_at(
                output, observed_descriptor, token, parent_descriptor
            )
        )
    except (OSError, PublicationError):
        return False
    finally:
        os.close(observed_descriptor)


def canonical_parent_resolves_to_held_parent(
    output: Path, parent_descriptor: int
) -> bool:
    """Return whether the parent pathname still names the held locked directory."""

    try:
        observed_descriptor = os.open(
            normalized_filesystem_path(output.parent), directory_open_flags()
        )
    except OSError:
        return False
    try:
        return os.path.samestat(
            os.fstat(parent_descriptor), os.fstat(observed_descriptor)
        )
    except OSError:
        return False
    finally:
        os.close(observed_descriptor)


def require_stable_publication_parent(
    output: Path, parent_descriptor: int, phase: str
) -> None:
    """Require the trusted locked parent pathname to remain stable."""

    if not canonical_parent_resolves_to_held_parent(output, parent_descriptor):
        raise PublicationError(
            "publication output parent no longer resolves to the held trusted "
            f"stable parent boundary {phase}: {output.parent}"
        )


def write_durable_publication_ownership(
    output: Path, parent_descriptor: int, output_descriptor: int, token: bytes
) -> None:
    """Durably bind sibling ownership to the token in the held output."""

    if not canonical_path_resolves_to_held_output(
        output,
        output_descriptor,
        token,
        require_ownership=False,
        parent_descriptor=parent_descriptor,
    ):
        raise PublicationError(
            "canonical publication path no longer resolves to the authorized "
            f"output directory instance: {output}"
        )
    path = publication_ownership_path(output)
    require_stable_publication_parent(
        output, parent_descriptor, "before publication ownership write"
    )
    write_durable_regular_file_at(
        parent_descriptor,
        path.name,
        json_payload(publication_ownership_record(output, token)),
    )
    require_stable_publication_parent(
        output, parent_descriptor, "after publication ownership write"
    )
    if not canonical_path_resolves_to_held_output(
        output,
        output_descriptor,
        token,
        require_ownership=True,
        parent_descriptor=parent_descriptor,
    ):
        raise PublicationError(
            f"publication ownership record failed post-write validation: {path}"
        )


def manifest_output_path(record: dict[str, Any]) -> Path | None:
    """Return the output path bound by one normalized publication invocation."""

    invocation = record.get("normalized_invocation")
    if not isinstance(invocation, list):
        return None
    indices = [
        index
        for index, value in enumerate(invocation)
        if value == "--output" and index + 1 < len(invocation)
    ]
    if len(indices) != 1 or not isinstance(invocation[indices[0] + 1], str):
        return None
    return normalized_filesystem_path(Path(invocation[indices[0] + 1]))


def valid_existing_publication_manifest(output: Path) -> bool:
    """Return whether the existing canonical manifest proves output ownership."""

    manifest_path = output / "manifest.json"
    if manifest_path.is_symlink() or not manifest_path.is_file():
        return False
    try:
        manifest = load_json(manifest_path)
    except (OSError, json.JSONDecodeError, PublicationError):
        return False
    if (
        manifest.get("record_type") != "cgl_lf_stage_i_fast_publication_products"
        or manifest_output_path(manifest) != normalized_filesystem_path(output)
    ):
        return False
    products = manifest.get("products")
    if not isinstance(products, list) or not products:
        return False
    seen: set[Path] = set()
    for index, binding in enumerate(products):
        canonical = binding_path(binding)
        if canonical is None or ".." in canonical.parts:
            return False
        canonical = normalized_filesystem_path(canonical)
        if (
            canonical == normalized_filesystem_path(output / "manifest.json")
            or not canonical.is_relative_to(normalized_filesystem_path(output))
            or canonical in seen
            or canonical.is_symlink()
            or verify_file_binding(binding, f"existing publication product {index}")
        ):
            return False
        seen.add(canonical)
    return True


def valid_existing_publication_manifest_at(
    output: Path, output_descriptor: int
) -> bool:
    """Validate a legacy manifest while proving canonical path still names held output."""

    return (
        canonical_path_resolves_to_held_output(
            output, output_descriptor, require_ownership=False
        )
        and valid_existing_publication_manifest(output)
        and canonical_path_resolves_to_held_output(
            output, output_descriptor, require_ownership=False
        )
    )


def require_or_create_publication_ownership(
    output: Path, parent_descriptor: int, output_descriptor: int
) -> bytes:
    """Require token authority, safely creating or migrating it before withdrawal."""

    ownership_path = publication_ownership_path(output)
    try:
        ownership_stat = os.stat(
            ownership_path.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
    except FileNotFoundError:
        ownership_stat = None
    except OSError as error:
        raise PublicationError(
            f"cannot inspect publication ownership record: {ownership_path}"
        ) from error
    if ownership_stat is not None and stat.S_ISLNK(ownership_stat.st_mode):
        raise PublicationError(
            f"publication ownership record may not be a symlink: {ownership_path}"
        )
    token = publication_authority_token(output_descriptor)
    if (
        token is not None
        and valid_publication_ownership_at(
            output, output_descriptor, token, parent_descriptor
        )
        and canonical_path_resolves_to_held_output(
            output,
            output_descriptor,
            token,
            require_ownership=True,
            parent_descriptor=parent_descriptor,
        )
    ):
        return token

    legacy_manifest = valid_existing_publication_manifest_at(
        output, output_descriptor
    )
    names = set(os.listdir(output_descriptor))
    token_only = names == {PUBLICATION_AUTHORITY_TOKEN_NAME}
    if ownership_stat is not None and not legacy_manifest:
        raise PublicationError(
            "publication ownership record does not bind the current "
            f"output directory instance: {ownership_path}"
        )
    if token is None:
        if names and not legacy_manifest:
            raise PublicationError(
                "refusing publication commit to nonempty unowned output: "
                f"{output}"
            )
        token = create_publication_authority_token(output_descriptor)
    elif not token_only and not legacy_manifest:
        raise PublicationError(
            "refusing publication commit to nonempty unowned output containing "
            f"an unbound authority token: {output}"
        )
    write_durable_publication_ownership(
        output, parent_descriptor, output_descriptor, token
    )
    return token


@contextmanager
def exclusive_publication_lock(output: Path) -> Iterator[int]:
    """Hold parent-bound and legacy locks for one publication commit."""

    lock_path = publication_lock_path(output)
    try:
        parent_descriptor = os.open(
            normalized_filesystem_path(output.parent), directory_open_flags()
        )
    except OSError as error:
        raise PublicationError(
            f"cannot open publication parent directory for locking: {output.parent}"
        ) from error
    parent_locked = False
    descriptor: int | None = None
    locked = False
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        # Locking the held parent directory makes replacement of the legacy sibling
        # lock pathname irrelevant to publishers using this implementation.
        fcntl.flock(parent_descriptor, fcntl.LOCK_EX)
        parent_locked = True
        require_stable_publication_parent(
            output, parent_descriptor, "after acquiring the parent lock"
        )
        try:
            descriptor = os.open(
                lock_path.name, flags, 0o600, dir_fd=parent_descriptor
            )
        except OSError as error:
            raise PublicationError(
                "cannot open publication sibling lock without following symlinks: "
                f"{lock_path}"
            ) from error
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise PublicationError(
                f"publication sibling lock is not a regular file: {lock_path}"
            )
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        locked = True
        require_stable_publication_parent(
            output, parent_descriptor, "before entering the publication commit"
        )
        yield parent_descriptor
    finally:
        if locked and descriptor is not None:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        if descriptor is not None:
            os.close(descriptor)
        if parent_locked:
            fcntl.flock(parent_descriptor, fcntl.LOCK_UN)
        os.close(parent_descriptor)


def validate_canonical_output_tree(output: Path) -> None:
    """Reject canonical output roots and descendants that use symlinks."""

    if output.is_symlink():
        raise PublicationError(
            f"publication output root may not be a symlink: {output}"
        )
    if not output.exists():
        return
    if not output.is_dir():
        raise PublicationError(f"publication output root is not a directory: {output}")
    for root, directories, files in os.walk(output, topdown=True, followlinks=False):
        root_path = Path(root)
        for name in sorted([*directories, *files]):
            path = root_path / name
            if path.is_symlink():
                raise PublicationError(
                    f"canonical publication descendant may not be a symlink: {path}"
                )


def canonical_tree_entries(
    output_descriptor: int,
) -> tuple[set[Path], set[Path]]:
    """Return regular files and directories beneath held output without symlinks."""

    files: set[Path] = set()
    directories: set[Path] = {Path(".")}

    def visit(descriptor: int, prefix: Path) -> None:
        try:
            names = sorted(os.listdir(descriptor))
        except OSError as error:
            raise PublicationError("cannot enumerate canonical publication tree") from error
        for name in names:
            relative = Path(name) if prefix == Path(".") else prefix / name
            try:
                observed = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            except OSError as error:
                raise PublicationError(
                    f"cannot inspect canonical publication descendant: {relative}"
                ) from error
            if stat.S_ISLNK(observed.st_mode):
                raise PublicationError(
                    f"canonical publication descendant may not be a symlink: {relative}"
                )
            if stat.S_ISREG(observed.st_mode):
                files.add(relative)
                continue
            if not stat.S_ISDIR(observed.st_mode):
                raise PublicationError(
                    f"canonical publication descendant has unsupported type: {relative}"
                )
            directories.add(relative)
            try:
                child = os.open(name, directory_open_flags(), dir_fd=descriptor)
            except OSError as error:
                raise PublicationError(
                    f"cannot open canonical publication directory: {relative}"
                ) from error
            try:
                visit(child, relative)
            finally:
                os.close(child)

    visit(output_descriptor, Path("."))
    return files, directories


def expected_relative_directories(products: Iterable[Path]) -> set[Path]:
    """Return exact relative directories required by canonical products."""

    expected = {Path(".")}
    for relative in products:
        if relative.is_absolute() or ".." in relative.parts:
            raise PublicationError(
                f"canonical publication product escapes output root: {relative}"
            )
        parent = relative.parent
        while parent != Path("."):
            expected.add(parent)
            parent = parent.parent
    return expected


def require_recognized_canonical_tree(
    output_descriptor: int,
    products: Iterable[Path],
    directories: set[Path],
) -> None:
    """Fail closed if canonical output contains stale or unrecognized paths."""

    product_paths = set(products)
    allowed_files = {
        *product_paths,
        Path("manifest.json"),
        Path(PUBLICATION_AUTHORITY_TOKEN_NAME),
    }
    files, actual_directories = canonical_tree_entries(output_descriptor)
    unexpected = sorted(
        [path.as_posix() for path in files - allowed_files]
        + [path.as_posix() + "/" for path in actual_directories - directories]
    )
    if unexpected:
        raise PublicationError(
            "canonical publication tree contains stale or unrecognized paths; "
            "refusing to delete them: " + ", ".join(unexpected)
        )


def withdraw_authenticated_canonical_manifest(
    output_descriptor: int,
    expected_payload: bytes,
    *,
    parent_descriptor: int | None = None,
    output_name: str = "publication",
) -> bool:
    """Quarantine and authenticate authority without deleting raced pathnames."""

    try:
        observed = os.stat(
            "manifest.json", dir_fd=output_descriptor, follow_symlinks=False
        )
    except FileNotFoundError:
        return False
    except OSError as error:
        raise PublicationError("cannot inspect canonical publication manifest") from error
    if not stat.S_ISREG(observed.st_mode):
        return False

    owned_parent_descriptor = parent_descriptor is None
    if parent_descriptor is None:
        try:
            parent_descriptor = os.open("..", directory_open_flags(), dir_fd=output_descriptor)
        except OSError as error:
            raise PublicationError(
                "cannot open publication parent for manifest quarantine"
            ) from error

    quarantine_name: str | None = None
    for _ in range(16):
        candidate = (
            f".{output_name}{PUBLICATION_MANIFEST_QUARANTINE_PREFIX}"
            f"{secrets.token_hex(32)}"
        )
        try:
            os.mkdir(candidate, mode=0o700, dir_fd=parent_descriptor)
            os.fsync(parent_descriptor)
            quarantine_name = candidate
            break
        except FileExistsError:
            continue
        except OSError as error:
            raise PublicationError(
                "cannot create publication-manifest quarantine directory"
            ) from error
    if quarantine_name is None:
        raise PublicationError(
            "cannot allocate a unique publication-manifest quarantine directory"
        )
    try:
        quarantine_descriptor = os.open(
            quarantine_name, directory_open_flags(), dir_fd=parent_descriptor
        )
    except OSError as error:
        raise PublicationError(
            "cannot open publication-manifest quarantine directory"
        ) from error

    try:
        try:
            os.rename(
                "manifest.json",
                "manifest.json",
                src_dir_fd=output_descriptor,
                dst_dir_fd=quarantine_descriptor,
            )
            os.fsync(quarantine_descriptor)
            os.fsync(output_descriptor)
            os.fsync(parent_descriptor)
        except FileNotFoundError:
            return False
        except OSError as error:
            raise PublicationError(
                "cannot quarantine canonical publication manifest"
            ) from error

        try:
            quarantined_payload = read_regular_file_at(
                quarantine_descriptor, Path("manifest.json")
            )
        except PublicationError:
            return False
        if not secrets.compare_digest(quarantined_payload, expected_payload):
            try:
                # Restore an unrecognized manifest without overwriting a canonical
                # replacement and without deleting either pathname afterward.
                os.link(
                    "manifest.json",
                    "manifest.json",
                    src_dir_fd=quarantine_descriptor,
                    dst_dir_fd=output_descriptor,
                    follow_symlinks=False,
                )
                os.fsync(output_descriptor)
            except FileExistsError:
                pass
            except OSError as error:
                raise PublicationError(
                    "cannot restore unrecognized quarantined publication manifest"
                ) from error
            return False
        # The authenticated prior authority remains as non-canonical sibling evidence.
        # POSIX has no compare-and-unlink primitive; deleting its mutable pathname after
        # authentication could delete an unknown same-user replacement.
        return True
    finally:
        os.close(quarantine_descriptor)
        if owned_parent_descriptor:
            os.close(parent_descriptor)


def source_binding_at(
    output: Path, output_descriptor: int, relative: Path
) -> dict[str, object]:
    """Return a canonical binding read through the held output capability."""

    payload = read_regular_file_at(output_descriptor, relative)
    return {
        "path": str(output / relative),
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def recognized_canonical_manifest_payload_at(
    output: Path, output_descriptor: int
) -> bytes | None:
    """Return exact recognized canonical manifest bytes or fail closed."""

    try:
        payload = read_regular_file_at(output_descriptor, Path("manifest.json"))
    except PublicationError as error:
        try:
            os.stat(
                "manifest.json", dir_fd=output_descriptor, follow_symlinks=False
            )
        except FileNotFoundError:
            return None
        raise PublicationError(
            "canonical publication manifest is unrecognized; refusing to replace it"
        ) from error
    try:
        manifest = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PublicationError(
            "canonical publication manifest is unrecognized; refusing to replace it"
        ) from error
    if (
        not isinstance(manifest, dict)
        or manifest.get("record_type")
        != "cgl_lf_stage_i_fast_publication_products"
        or manifest_output_path(manifest) != normalized_filesystem_path(output)
    ):
        raise PublicationError(
            "canonical publication manifest is unrecognized; refusing to replace it"
        )
    products = manifest.get("products")
    if not isinstance(products, list) or not products:
        raise PublicationError(
            "canonical publication manifest is unrecognized; refusing to replace it"
        )
    seen: set[Path] = set()
    for binding in products:
        canonical = binding_path(binding)
        if canonical is None or canonical.is_symlink() or ".." in canonical.parts:
            raise PublicationError(
                "canonical publication manifest is unrecognized; refusing to replace it"
            )
        try:
            relative = canonical.relative_to(output)
        except ValueError as error:
            raise PublicationError(
                "canonical publication manifest is unrecognized; refusing to replace it"
            ) from error
        if (
            not relative.parts
            or relative in seen
            or relative == Path("manifest.json")
            or relative == Path(PUBLICATION_AUTHORITY_TOKEN_NAME)
            or relative.parts[0].startswith(PUBLICATION_MANIFEST_QUARANTINE_PREFIX)
        ):
            raise PublicationError(
                "canonical publication manifest is unrecognized; refusing to replace it"
            )
        seen.add(relative)
        try:
            observed = source_binding_at(output, output_descriptor, relative)
        except PublicationError as error:
            raise PublicationError(
                "canonical publication manifest is unrecognized; refusing to replace it"
            ) from error
        if observed != binding:
            raise PublicationError(
                "canonical publication manifest is unrecognized; refusing to replace it"
            )
    return payload


def verify_exact_canonical_products(
    output: Path,
    output_descriptor: int,
    products: Iterable[Path],
    directories: set[Path],
    bindings: Iterable[dict[str, object]],
    *,
    manifest_published: bool,
) -> None:
    """Verify the exact token-bearing canonical tree through held output."""

    expected_products = set(products)
    expected_files = {
        *expected_products,
        Path(PUBLICATION_AUTHORITY_TOKEN_NAME),
    }
    if manifest_published:
        expected_files.add(Path("manifest.json"))
    actual_files, actual_directories = canonical_tree_entries(output_descriptor)
    if actual_files != expected_files or actual_directories != directories:
        raise PublicationError(
            "canonical publication tree differs from the staged publication"
        )
    for binding in bindings:
        canonical = Path(str(binding["path"]))
        try:
            relative = canonical.relative_to(output)
        except ValueError as error:
            raise PublicationError(
                f"promoted publication product escapes output root: {canonical}"
            ) from error
        if source_binding_at(output, output_descriptor, relative) != binding:
            raise PublicationError(
                f"promoted publication product differs from staged binding: {canonical}"
            )


def promote_staged_publication(
    staging: Path,
    output: Path,
    products: Iterable[Path],
    manifest: dict[str, object],
    *,
    analysis: Path,
    evidence_paths: Iterable[Path],
    acceptance_roots: Iterable[Path] = (),
) -> Path:
    """Promote products atomically and publish their canonical manifest last."""

    staging = staging.absolute()
    output = output.absolute()
    analysis = analysis.absolute()
    evidence_paths = tuple(path.absolute() for path in evidence_paths)
    acceptance_roots = (
        *builtin_acceptance_roots(analysis),
        *(path.absolute() for path in acceptance_roots),
    )
    validate_publication_output_overlap(
        analysis, output, evidence_paths, acceptance_roots
    )
    if staging.parent.resolve(strict=True) != output.parent.resolve(strict=True):
        raise PublicationError("publication staging directory must be a sibling of output")
    records = canonical_product_records(products, staging, output)
    expected_bindings = [binding for _, _, binding in records]
    if manifest.get("products") != expected_bindings:
        raise PublicationError(
            "publication manifest products differ from staged canonical bindings"
        )

    staged_manifest = staging / "manifest.json"
    staged_manifest_payload = json_payload(manifest)
    atomic_write(staged_manifest, staged_manifest_payload)
    canonical_products = [canonical for _, canonical, _ in records]
    relative_products = [
        canonical.relative_to(output) for canonical in canonical_products
    ]
    canonical_directories = expected_relative_directories(relative_products)
    canonical_manifest = output / "manifest.json"

    with exclusive_publication_lock(output) as parent_descriptor:
        require_stable_publication_parent(
            output, parent_descriptor, "before canonical output validation"
        )
        validate_publication_output_overlap(
            analysis, output, evidence_paths, acceptance_roots
        )
        validate_canonical_output_tree(output)
        if not output.exists():
            if publication_ownership_path(output).exists():
                raise PublicationError(
                    "publication ownership record exists while canonical output is "
                    f"absent: {publication_ownership_path(output)}"
                )
            try:
                os.mkdir(output.name, mode=0o755, dir_fd=parent_descriptor)
                os.fsync(parent_descriptor)
            except OSError as error:
                raise PublicationError(
                    f"cannot create publication output under held parent: {output}"
                ) from error
        require_stable_publication_parent(
            output, parent_descriptor, "before opening canonical output"
        )
        output_descriptor = open_output_directory_at(parent_descriptor, output.name)
        manifest_promotion_attempted = False
        try:
            require_recognized_canonical_tree(
                output_descriptor, relative_products, canonical_directories
            )
            prior_manifest_payload = recognized_canonical_manifest_payload_at(
                output, output_descriptor
            )
            token = require_or_create_publication_ownership(
                output, parent_descriptor, output_descriptor
            )
            require_recognized_canonical_tree(
                output_descriptor, relative_products, canonical_directories
            )
            if not canonical_path_resolves_to_held_output(
                output,
                output_descriptor,
                token,
                require_ownership=True,
                parent_descriptor=parent_descriptor,
            ):
                raise PublicationError(
                    "publication ownership no longer binds the canonical output "
                    f"directory instance: {output}"
                )
            if recognized_canonical_manifest_payload_at(
                output, output_descriptor
            ) != prior_manifest_payload:
                raise PublicationError(
                    "canonical publication manifest changed while establishing "
                    "publication ownership; unrecognized manifest bytes were preserved"
                )

            # Flush every staged product, including PDFs, before withdrawing prior
            # authority or mutating any canonical product.
            for staged, _, _ in records:
                fsync_regular_file(staged)
            fsync_regular_file(staged_manifest)

            # Withdraw authority before replacing recognized products. No stale or
            # unrecognized path is ever removed.
            require_stable_publication_parent(
                output, parent_descriptor, "before prior-authority withdrawal"
            )
            if (
                prior_manifest_payload is not None
                and not withdraw_authenticated_canonical_manifest(
                    output_descriptor,
                    prior_manifest_payload,
                    parent_descriptor=parent_descriptor,
                    output_name=output.name,
                )
            ):
                raise PublicationError(
                    "canonical publication manifest changed before authority "
                    "withdrawal; unrecognized manifest bytes were preserved"
                )
            for (staged, _, _), relative in zip(records, relative_products):
                require_stable_publication_parent(
                    output, parent_descriptor, "before product promotion"
                )
                promote_file(staged, relative, output_descriptor)
            verify_exact_canonical_products(
                output,
                output_descriptor,
                relative_products,
                canonical_directories,
                expected_bindings,
                manifest_published=False,
            )
            if not canonical_path_resolves_to_held_output(
                output,
                output_descriptor,
                token,
                require_ownership=True,
                parent_descriptor=parent_descriptor,
            ):
                raise PublicationError(
                    "canonical publication path changed during product promotion"
                )
            require_stable_publication_parent(
                output, parent_descriptor, "before manifest promotion"
            )
            if not secrets.compare_digest(
                staged_manifest.read_bytes(), staged_manifest_payload
            ):
                raise PublicationError(
                    "staged publication manifest differs from exact expected bytes"
                )
            manifest_promotion_attempted = True
            promote_file(staged_manifest, Path("manifest.json"), output_descriptor)
            verify_exact_canonical_products(
                output,
                output_descriptor,
                relative_products,
                canonical_directories,
                expected_bindings,
                manifest_published=True,
            )
            fsync_exact_regular_file_at(
                output_descriptor, Path("manifest.json"), staged_manifest_payload
            )
            if not canonical_path_resolves_to_held_output(
                output,
                output_descriptor,
                token,
                require_ownership=True,
                parent_descriptor=parent_descriptor,
            ):
                raise PublicationError(
                    "canonical publication path changed before commit return"
                )
            require_stable_publication_parent(
                output, parent_descriptor, "before publication commit return"
            )
        except BaseException as error:
            if manifest_promotion_attempted:
                try:
                    withdrawn = withdraw_authenticated_canonical_manifest(
                        output_descriptor,
                        staged_manifest_payload,
                        parent_descriptor=parent_descriptor,
                        output_name=output.name,
                    )
                    if not withdrawn and hasattr(error, "add_note"):
                        error.add_note(
                            "canonical manifest was absent or no longer matched this "
                            "commit; no unrecognized path was deleted"
                        )
                except BaseException as withdrawal_error:
                    if hasattr(error, "add_note"):
                        error.add_note(
                            "failed to withdraw authenticated canonical manifest: "
                            f"{withdrawal_error}"
                        )
            raise
        finally:
            os.close(output_descriptor)
    return canonical_manifest


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
        "- The compact numerical-health/provenance table excludes unsupported retained-"
        "state floor margins. Missing authenticated mass, energy, CT, or strict-failure "
        "evidence remains inconclusive.",
        "- Primary full-window scalar rows reproduce authenticated direct-acceptance "
        "source-history statistics without additional publication-layer normalization.",
        "- All-snapshot directional-discriminant coverage passes only for an authenticated "
        "snapshot_policy=all selection with an exactly once completed result. "
        "Dispositions remain nonnegative-discriminant, negative, or nonfinite and are "
        "assigned separately to legacy-implementation or literature-correct formula "
        "evidence. Latest-only, missing, stale, or unidentified-formula results remain "
        "inconclusive. Static experiment scope, formula/executable compatibility, and "
        "strict-hyperbolic claim eligibility are reported separately. Retained "
        "cell-centered coordinate-normal coverage cannot establish face, intermediate, "
        "oblique, full, or strict hyperbolicity or absence of sqrt(|D|) fallback.",
        "- Signed LF applied-stage ledgers remain distinct from sparse retained-"
        "snapshot pressure-work and heat-flux reconstructions; signs are preserved.",
        "- MKS24 panel and lineage disposition tables report authenticated admitted, "
        "blocked/external, selected, superseded, failed, and unselected evidence "
        "without promoting absent evidence.",
        "- The coherent-direction mechanism summary is descriptive only and does not "
        "define or imply a new preregistered pass gate.",
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
        "causal_mechanism": figures / "fig08_causal_mechanism.pdf",
        "resolution_curves": figures / "fig09_resolution_curves.pdf",
        "hyperbolicity_coverage": figures / "fig10_hyperbolicity_coverage.pdf",
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
    render_causal_mechanism(data, plt, figure_paths["causal_mechanism"])
    render_resolution_curves(data, plt, figure_paths["resolution_curves"])
    render_hyperbolicity_coverage(
        data, plt, colors, patches, figure_paths["hyperbolicity_coverage"]
    )
    products.extend(figure_paths.values())

    table_specs = (
        (
            "numerical_health_provenance",
            [
                "case_id", "completion", "final_time", "fatal_lf_counters",
                "fatal_counter_maximum", "mass_relative_drift_maximum",
                "mhd_user_mass_relative_mismatch", "active_energy_closure",
                "energy_increment_residual_maximum",
                "direct_ct_numerical", "maximum_normalized_ct_divb",
                "claim_scope", "strict_failure", "strict_failure_counters",
                "provenance_summary",
            ],
            numerical_health_provenance_rows(data),
        ),
        (
            "primary_full_window_scalars",
            [
                "case_id", "metric", "history", "column", "availability", "mean",
                "standard_error", "ci95_lower", "ci95_upper", "effective_sample_count",
                "stationarity", "sampling_adequacy", "acceptance", "claim_scope",
            ],
            primary_full_window_scalar_rows(data),
        ),
        (
            "hyperbolicity_all_snapshot_coverage",
            [
                "case_id", "coverage_result", "numerical_result", "snapshot_policy",
                "complete_retained_snapshot_count", "selected_snapshot_count",
                "audited_snapshot_count",
                "all_complete_retained_snapshots_selected", "time_first", "time_last",
                "negative_discriminant_count", "nonfinite_discriminant_count",
                "cell_direction_evaluations", "minimum_discriminant",
                "coverage_reason", "experiment_scope",
                "formula_id", "formula_disposition_family",
                "formula_id_provenance", "formula_provenance",
                "formula_provenance_status", "executable_formula_id",
                "executable_formula_id_provenance",
                "formula_executable_compatibility",
                "formula_executable_compatibility_provenance",
                "formula_executable_compatibility_reason",
                "legacy_implementation_disposition",
                "literature_correct_disposition",
                "strict_hyperbolic_claim_status",
                "strict_hyperbolic_claim_eligible",
                "strict_hyperbolic_claim_reason", "audit_state_scope",
                "audit_direction_scope", "selection_provenance",
                "result_provenance",
            ],
            hyperbolicity_coverage_rows(data),
        ),
        (
            "signed_lf_cap_work_ledger",
            [
                "case_id", "availability", "diagnostics_provenance",
                "applied_heat_flux_availability",
                "applied_pressure_work_availability",
                "reconstructed_pressure_availability",
                "reconstructed_heat_flux_availability",
                "applied_ledgers_signed", "applied_heat_flux_parallel",
                "applied_heat_flux_perpendicular", "applied_heat_flux_total",
                "applied_pressure_work_total", "applied_pressure_work_anisotropic",
                "cap_parallel_over_1", "cap_parallel_over_10",
                "cap_perpendicular_over_1", "cap_perpendicular_over_10",
                "reconstructed_pressure_snapshot_count",
                "reconstructed_pressure_applied_to_flow",
                "reconstructed_isotropic_perpendicular_pressure_power_mean",
                "reconstructed_anisotropic_stress_power_mean",
                "reconstructed_total_cgl_pressure_power_mean",
                "reconstructed_anisotropic_stress_power_integral",
                "reconstructed_heat_flux_snapshot_count",
                "reconstructed_regularized_heat_flux_power_mean",
                "reconstructed_unlimited_heat_flux_power_mean",
                "reconstructed_parallel_cap_active_volume_fraction_mean",
                "reconstructed_perpendicular_cap_active_volume_fraction_mean",
                "reconstructed_regularized_heat_flux_power_integral",
                "reconstructed_unlimited_heat_flux_power_integral", "semantics",
            ],
            signed_lf_cap_work_ledger_rows(data),
        ),
        (
            "mks24_panel_dispositions",
            [
                "panel", "disposition", "result", "product_count", "pass_count",
                "fail_count", "inconclusive_count", "sources", "reason", "evidence",
            ],
            mks24_panel_disposition_rows(data),
        ),
        (
            "lineage_dispositions",
            [
                "case_id", "lineage_index", "selected", "disposition", "reason",
                "terminal_state", "source_family", "variant", "job_id",
                "observed_final_time", "run_exit_code", "restart_link_valid",
                "segment_count", "terminal_segment", "provenance",
            ],
            lineage_disposition_rows(data),
        ),
        (
            "coherent_direction_mechanism",
            [
                "metric", "available_pair_count",
                "positive_active_minus_passive_count",
                "negative_active_minus_passive_count", "equal_count",
                "inconclusive_pair_count", "descriptive_direction",
                "pair_active_minus_passive", "inference_scope",
            ],
            coherent_direction_mechanism_rows(data),
        ),
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
                "retained_state_coordinate_discriminant",
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
                "right_mean", "difference", "combined_standard_error",
                "standardized_effect", "intervention_estimand",
                "intervention_enabled_components", "excluded_interpretation",
                "intervention_declaration", "inference_scope",
                "current_science_scope_disposition",
                "full_scope_independent_review_complete", "reason", "authority",
                "release_authorizing", "claim_scope",
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
8. **Causal mechanism.** Matched active/passive magnetic-intermittency
   \\(C_{{B^2}}\\) histories, parallel-strain PDFs, and reviewed standardized effects.
   Curves require
   complete claim-eligible cases; strain PDFs additionally require diagnostics bytes
   bound by authenticated reviewed science. Missing evidence remains inconclusive.
9. **Resolution curves.** Actual R16/R02/R17 velocity and magnetic-fluctuation
   spectra and peak-alignment curves over the preregistered common range, accompanied
   by the reviewed convergence decisions and limits. All three authenticated curves
   are required in each panel.
10. **All-snapshot coordinate-direction discriminant coverage.** Authenticated
    active-CGL retained cell-centered coverage in the three coordinate-normal
    directions, numerical disposition by authenticated formula identity, static
    experiment scope, and formula/executable compatibility. Coverage passes only when
    every complete retained snapshot is represented exactly once in a completed bound
    result. Negative or nonfinite compatible literature-correct evidence can exclude a
    strict claim. Nonnegative retained-state coordinate-direction evidence cannot
    establish face, intermediate, oblique, full, or strict hyperbolicity or absence of
    sqrt(|D|) fallback.
"""
    captions_path = output / "captions.md"
    write_text(captions_path, captions)
    products.append(captions_path)

    report_path = output / "report.md"
    write_text(report_path, report_markdown(data, products, output))
    products.append(report_path)
    return products


def build_publication_manifest(
    data: PublicationData,
    analysis: Path,
    output: Path,
    acceptance_paths: Iterable[Path],
    staged_products: Iterable[Path],
    staging: Path,
) -> dict[str, object]:
    """Build the final-path authority for one fully staged publication."""

    import matplotlib

    return {
        "schema_version": 2,
        "record_type": "cgl_lf_stage_i_fast_publication_products",
        "evidence_state": publication_evidence_state(data),
        "analysis_output": str(analysis),
        "renderer": source_binding(RENDERER_PATH),
        "normalized_invocation": normalized_invocation(
            analysis, output, acceptance_paths
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
        "publication_commit_contract": {
            "authority": "manifest.json is published last and exact bytes are withdrawn",
            "manifest_withdrawal": (
                "rename-to-random-quarantine then authenticate; unknown bytes are "
                "restored or preserved without overwrite"
            ),
            "stable_parent_boundary": (
                "the output parent pathname must continue to name the held locked "
                "directory through commit return"
            ),
        },
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
        "products": canonical_product_bindings(staged_products, staging, output),
        "renderer_ingestion_warnings": data.ingestion_warnings,
        "case_warning_records": case_warning_rows(data),
    }


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
    acceptance_roots = (
        *builtin_acceptance_roots(analysis),
        *(path.absolute() for path in args.acceptance),
    )
    validate_publication_output_overlap(
        analysis, output, (), acceptance_roots
    )
    acceptance_candidates = discover_acceptance_paths(analysis, args.acceptance)
    data = discover_data(analysis, args.acceptance)
    evidence_paths = set(data.source_paths) | set(acceptance_candidates)
    validate_publication_output_overlap(
        analysis, output, evidence_paths, acceptance_roots
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output.parent, prefix=f".{output.name}.staging-"
    ) as staging_value:
        staging = Path(staging_value).absolute()
        products = render_products(data, staging)
        evidence_paths.update(data.source_paths)
        validate_publication_output_overlap(
            analysis, output, evidence_paths, acceptance_roots
        )
        manifest = build_publication_manifest(
            data, analysis, output, args.acceptance, products, staging
        )
        promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=analysis,
            evidence_paths=evidence_paths,
            acceptance_roots=acceptance_roots,
        )
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
