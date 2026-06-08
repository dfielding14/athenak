#!/usr/bin/env python3
"""Aggregate provenance-bound direct-fast Stage I scientific comparisons.

This downstream-only utility combines the direct-fast report inventory, case
diagnostics and snapshot products with direct-fast acceptance evidence.  It
uses the reviewed scientific-acceptance kernels for scalar extraction, active
/ passive descriptive contrasts, convergence distances, and MKS24 residual checks.

Partial campaigns remain useful but cannot silently become claim-grade:
observations are retained when available, while a gate can pass only when all
contributing direct-fast case-acceptance records themselves pass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
REVIEWED_UTILITY = REPO_ROOT / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
FAST_ACCEPTANCE_UTILITY = REPO_ROOT / "scripts/frontier/cgl_lf_stage_i_fast_acceptance.py"
FAST_REPORT_UTILITY = REPO_ROOT / "scripts/frontier/cgl_lf_stage_i_fast_report.py"
PAPER_ANALYZER = REPO_ROOT / "scripts/analyze_cgl_lf_paper.py"
DEFAULT_CAMPAIGN_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
DEFAULT_INVENTORY = DEFAULT_CAMPAIGN_ROOT / (
    "analysis/mks24-stage-i-fast/E03-forcing-policy/R02-R17-live/inventory.json"
)
DEFAULT_CRITERIA = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
DEFAULT_CRITERIA_REVIEW = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)
CASE_ID_PATTERN = re.compile(r"R(?:0[2-9]|1[0-7])")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
ACTIVE_PASSIVE_PAIRS = (
    ("R02", "R06"),
    ("R03", "R07"),
    ("R04", "R08"),
    ("R05", "R09"),
)
CONTRAST_FAMILIES = {
    "forcing": (
        ("R02", "R04"),
        ("R03", "R05"),
        ("R06", "R08"),
        ("R07", "R09"),
    ),
    "beta": (
        ("R02", "R03"),
        ("R04", "R05"),
        ("R06", "R07"),
        ("R08", "R09"),
    ),
    "tcorr": (("R05", "R11"),),
    "lf_strength": (
        ("R12", "R02"),
        ("R13", "R02"),
    ),
    "limiter": (
        ("R15", "R14"),
        ("R14", "R03"),
        ("R15", "R03"),
        ("R03", "R07"),
    ),
}
VALID_RESULTS = {"pass", "fail", "inconclusive"}
WINDOW_TIME_TOLERANCE = 1.0e-12


class ScienceError(RuntimeError):
    """Raised when an input is malformed, stale, or provenance-inconsistent."""


def reject_json_constant(value: str) -> None:
    """Reject non-standard NaN and infinity constants."""

    raise ValueError(f"non-finite JSON constant is forbidden: {value}")


def sha256_file(path: Path) -> str:
    """Return a lowercase SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: object) -> bytes:
    """Return canonical JSON bytes used for evidence and inventory digests."""

    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def json_safe(value: object) -> object:
    """Return a deterministic JSON-safe representation."""

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


def stable_json_bytes(value: object) -> bytes:
    """Return readable deterministic JSON bytes."""

    return (
        json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")


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
    """Write deterministic JSON."""

    atomic_write(path, stable_json_bytes(value))


def write_text(path: Path, value: str) -> None:
    """Write deterministic UTF-8 text."""

    atomic_write(path, value.encode("utf-8"))


def load_json(path: Path, label: str) -> dict[str, object]:
    """Load one finite JSON object."""

    try:
        payload = path.read_bytes()
        value = json.loads(
            payload.decode("utf-8"), parse_constant=reject_json_constant
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ScienceError(f"cannot load {label} {path}: {error}") from error
    if not isinstance(value, dict):
        raise ScienceError(f"{label} is not a JSON object: {path}")
    return value


def binding(path: Path) -> dict[str, object]:
    """Bind one existing regular file."""

    resolved = path.expanduser().absolute().resolve(strict=True)
    if not resolved.is_file():
        raise ScienceError(f"expected a regular file: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def binding_identity(value: object, label: str) -> tuple[str, int, str]:
    """Return one verified binding's normalized identity."""

    if not isinstance(value, dict):
        raise ScienceError(f"{label} is not a file binding")
    path_value = value.get("path")
    digest = value.get("sha256")
    size = value.get("size_bytes")
    if not isinstance(path_value, str) or not path_value:
        raise ScienceError(f"{label} lacks a path")
    if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
        raise ScienceError(f"{label} lacks a valid SHA-256")
    if not isinstance(size, int) or isinstance(size, bool) or size < 0:
        raise ScienceError(f"{label} lacks an integer size_bytes")
    path = Path(path_value).expanduser().absolute().resolve(strict=True)
    observed = binding(path)
    if observed["sha256"] != digest or observed["size_bytes"] != size:
        raise ScienceError(f"{label} differs from its declared binding")
    return str(path), size, digest


def require_same_binding(
    value: object, expected: object, label: str
) -> dict[str, object]:
    """Verify and require two bindings to identify the same current bytes."""

    observed_identity = binding_identity(value, label)
    expected_identity = binding_identity(expected, f"expected {label}")
    if observed_identity != expected_identity:
        raise ScienceError(f"{label} does not bind the selected current artifact")
    return {
        "path": observed_identity[0],
        "size_bytes": observed_identity[1],
        "sha256": observed_identity[2],
    }


def require_same_content_binding(
    value: object, expected: object, label: str
) -> dict[str, object]:
    """Verify two bindings and require identical bytes, allowing distinct copies."""

    observed_identity = binding_identity(value, label)
    expected_identity = binding_identity(expected, f"expected {label}")
    if observed_identity[1:] != expected_identity[1:]:
        raise ScienceError(f"{label} differs from the selected authoritative bytes")
    return {
        "path": observed_identity[0],
        "size_bytes": observed_identity[1],
        "sha256": observed_identity[2],
    }


def recursive_file_bindings(value: object) -> list[dict[str, object]]:
    """Return every nested path/SHA/size binding."""

    result: list[dict[str, object]] = []
    if isinstance(value, dict):
        if (
            isinstance(value.get("path"), str)
            and "sha256" in value
            and "size_bytes" in value
        ):
            result.append(value)
        for child in value.values():
            result.extend(recursive_file_bindings(child))
    elif isinstance(value, list):
        for child in value:
            result.extend(recursive_file_bindings(child))
    return result


def verify_recursive_bindings(value: object, label: str) -> None:
    """Verify every unique nested binding."""

    seen: set[tuple[str, str]] = set()
    for index, item in enumerate(recursive_file_bindings(value)):
        key = (str(item.get("path")), str(item.get("sha256")))
        if key in seen:
            continue
        seen.add(key)
        binding_identity(item, f"{label} binding {index}")


def load_module(name: str, path: Path) -> object:
    """Load one local Python utility by exact path."""

    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ScienceError(f"cannot import utility: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_reviewed_module() -> object:
    """Load the reviewed scientific-acceptance utility."""

    return load_module("_cgl_lf_fast_science_reviewed", REVIEWED_UTILITY)


def load_report_module() -> object:
    """Load the direct-fast report utility."""

    return load_module("_cgl_lf_fast_science_report", FAST_REPORT_UTILITY)


def path_contains(parent: Path, child: Path) -> bool:
    """Return whether child is parent or is located below parent."""

    return child == parent or parent in child.parents


def reject_nested_output_symlink_escapes(output: Path) -> None:
    """Reject existing nested symlinks that resolve outside the output root."""

    if not output.exists():
        return
    root = output.resolve(strict=True)
    for directory, names, files in os.walk(root, followlinks=False):
        parent = Path(directory)
        for name in [*names, *files]:
            candidate = parent / name
            if candidate.is_symlink() and not path_contains(
                root, candidate.resolve(strict=False)
            ):
                raise ScienceError(
                    f"science output contains a nested symlink escape: {candidate}"
                )


def require_separate_output(output: Path, report_root: Path, acceptance_root: Path) -> None:
    """Require an output tree separate from both input trees and the repository."""

    resolved_output = output.resolve(strict=False)
    for source in (report_root, acceptance_root, REPO_ROOT):
        resolved = source.resolve(strict=True)
        if path_contains(resolved, resolved_output) or path_contains(resolved_output, resolved):
            raise ScienceError(f"science output must be separate from input tree: {source}")
    reject_nested_output_symlink_escapes(resolved_output)


def parse_case_selection(values: list[str] | None, required: list[str]) -> list[str]:
    """Parse optional repeated comma-separated case selections."""

    if not values:
        return list(required)
    selected: set[str] = set()
    for value in values:
        for token in value.replace(",", " ").split():
            if CASE_ID_PATTERN.fullmatch(token) is None:
                raise ScienceError(f"invalid Stage I case ID: {token}")
            selected.add(token)
    return sorted(selected)


def history_binding_from_lineage(
    lineage: dict[str, object], kind: str, label: str
) -> dict[str, object]:
    """Return and verify one selected merged-history binding."""

    histories = lineage.get("histories")
    record = histories.get(kind) if isinstance(histories, dict) else None
    if not isinstance(record, dict) or record.get("available") is not True:
        raise ScienceError(f"{label} selected {kind} history is unavailable")
    declared = record.get("binding")
    path_value = record.get("path")
    if not isinstance(path_value, str):
        raise ScienceError(f"{label} selected {kind} history lacks a path")
    current = binding(Path(path_value))
    require_same_binding(declared, current, f"{label} selected {kind} history")
    return current


def validate_inventory(
    inventory_path: Path,
    policy: dict[str, object],
    selected: list[str],
) -> tuple[dict[str, object], Path, dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    """Authenticate the report inventory and selected case lineages."""

    inventory = load_json(inventory_path, "direct-fast report inventory")
    if inventory.get("schema_version") != 1:
        raise ScienceError("direct-fast report inventory schema differs")
    output_value = inventory.get("output")
    if not isinstance(output_value, str):
        raise ScienceError("direct-fast report inventory lacks an output root")
    report_root = Path(output_value).expanduser().absolute().resolve(strict=True)
    if inventory_path.resolve(strict=True) != (report_root / "inventory.json").resolve(
        strict=True
    ):
        raise ScienceError("direct-fast report inventory is outside its declared output root")
    stage_manifest = policy["verified_sources"]["stage_i_manifest"]
    matrix = require_same_content_binding(
        inventory.get("matrix"), stage_manifest, "direct-fast report matrix"
    )
    adapter = require_same_binding(
        inventory.get("adapter"), binding(FAST_REPORT_UTILITY), "direct-fast report adapter"
    )
    cases = inventory.get("cases")
    if not isinstance(cases, dict):
        raise ScienceError("direct-fast report inventory cases are malformed")
    lineages: dict[str, dict[str, object]] = {}
    lineage_bindings: dict[str, dict[str, object]] = {}
    for case_id in selected:
        inline = cases.get(case_id)
        lineage_path = report_root / "cases" / case_id / "lineage.json"
        if inline is None and not lineage_path.exists():
            continue
        if not isinstance(inline, dict) or not lineage_path.is_file():
            raise ScienceError(f"{case_id} inventory and lineage file coverage differ")
        lineage = load_json(lineage_path, f"{case_id} lineage")
        if canonical_json(inline) != canonical_json(lineage):
            raise ScienceError(f"{case_id} lineage differs from the authenticated inventory")
        if lineage.get("case_id") != case_id:
            raise ScienceError(f"{case_id} lineage identity differs")
        lineages[case_id] = lineage
        lineage_bindings[case_id] = binding(lineage_path)
    inventory["_science_verified_matrix"] = matrix
    inventory["_science_verified_adapter"] = adapter
    return inventory, report_root, lineages, lineage_bindings


def verify_policy_provenance(
    value: dict[str, object],
    policy: dict[str, object],
    label: str,
    utility_key: str,
) -> None:
    """Require exact current criteria, review, and reviewed-utility provenance."""

    require_same_binding(value.get("criteria"), policy["criteria_binding"], f"{label} criteria")
    require_same_binding(
        value.get("criteria_review"), policy["review_binding"], f"{label} criteria review"
    )
    require_same_binding(
        value.get(utility_key),
        binding(REVIEWED_UTILITY),
        f"{label} reviewed acceptance utility",
    )


def exact_science_scope_records(
    reviewed: object, policy: dict[str, object]
) -> tuple[dict[str, object], dict[str, object]]:
    """Return the exact current review limitation and total-intervention scope."""

    limitation = policy.get("current_science_scope_limitation")
    if (
        not isinstance(limitation, dict)
        or limitation != reviewed.CURRENT_SCIENCE_SCOPE_LIMITATION
    ):
        raise ScienceError("current science scope limitation differs")
    try:
        intervention = policy["criteria"]["family_gates"]["active_passive"][
            "intervention_scope"
        ]
    except (KeyError, TypeError) as error:
        raise ScienceError("active/passive intervention scope is unavailable") from error
    if (
        not isinstance(intervention, dict)
        or intervention != reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
    ):
        raise ScienceError("active/passive intervention scope differs")
    return limitation, intervention


def require_exact_science_scope_records(
    value: dict[str, object],
    limitation: dict[str, object],
    intervention: dict[str, object],
    label: str,
) -> None:
    """Require one downstream record to retain both exact scope declarations."""

    if value.get("current_science_scope_limitation") != limitation:
        raise ScienceError(f"{label} current science scope limitation differs")
    if value.get("active_passive_intervention_scope") != intervention:
        raise ScienceError(f"{label} active/passive intervention scope differs")


def validate_acceptance_root(
    acceptance_root: Path,
    inventory_path: Path,
    policy: dict[str, object],
) -> tuple[
    dict[str, dict[str, object]],
    dict[str, object],
    dict[str, object],
    dict[str, dict[str, object]],
]:
    """Authenticate direct-fast acceptance provenance and its outputs."""

    provenance_path = acceptance_root / "provenance.json"
    provenance = load_json(provenance_path, "direct-fast acceptance provenance")
    if provenance.get("record_type") != "cgl-lf-stage-i-direct-fast-acceptance-provenance":
        raise ScienceError("direct-fast acceptance provenance record type differs")
    inputs = provenance.get("inputs")
    outputs = provenance.get("outputs")
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise ScienceError("direct-fast acceptance provenance inputs or outputs are malformed")
    require_same_binding(
        inputs.get("inventory"), binding(inventory_path), "acceptance input inventory"
    )
    require_same_binding(
        inputs.get("driver"), binding(FAST_ACCEPTANCE_UTILITY), "fast-acceptance driver"
    )
    verify_policy_provenance(inputs, policy, "fast-acceptance input", "reviewed_acceptance_utility")
    verify_recursive_bindings(inputs, "fast-acceptance inputs")

    campaign_binding = outputs.get("campaign_evidence")
    campaign_path_text = (
        campaign_binding.get("path") if isinstance(campaign_binding, dict) else None
    )
    if not isinstance(campaign_path_text, str):
        raise ScienceError("acceptance provenance lacks campaign evidence output")
    campaign_path = Path(campaign_path_text)
    if campaign_path.resolve(strict=True) != (
        acceptance_root / "campaign_evidence.json"
    ).resolve(strict=True):
        raise ScienceError("campaign evidence output escapes its acceptance directory")
    require_same_binding(campaign_binding, binding(campaign_path), "campaign evidence output")
    campaign = load_json(campaign_path, "direct-fast campaign evidence")
    reviewed = load_reviewed_module()
    limitation, intervention = exact_science_scope_records(reviewed, policy)
    try:
        reviewed.verify_evidence_digest(campaign, "direct-fast campaign evidence")
    except Exception as error:
        raise ScienceError(f"direct-fast campaign evidence is forged: {error}") from error
    if campaign.get("record_type") != "cgl-lf-stage-i-direct-fast-campaign-evidence":
        raise ScienceError("direct-fast campaign evidence record type differs")
    campaign_provenance = campaign.get("provenance")
    if not isinstance(campaign_provenance, dict):
        raise ScienceError("direct-fast campaign evidence lacks provenance")
    verify_policy_provenance(
        campaign_provenance,
        policy,
        "direct-fast campaign evidence",
        "reviewed_acceptance_utility",
    )
    require_exact_science_scope_records(
        campaign, limitation, intervention, "direct-fast campaign evidence"
    )

    output_cases = outputs.get("case_acceptance")
    case_results = campaign.get("case_results")
    if not isinstance(output_cases, dict) or not isinstance(case_results, dict):
        raise ScienceError(
            "acceptance case output inventory or campaign case results are malformed"
        )
    if set(output_cases) != set(case_results):
        raise ScienceError("acceptance case output inventory differs from campaign case results")
    cases: dict[str, dict[str, object]] = {}
    case_bindings: dict[str, dict[str, object]] = {}
    for case_id, declared in sorted(output_cases.items()):
        if CASE_ID_PATTERN.fullmatch(str(case_id)) is None:
            raise ScienceError(f"acceptance output has invalid case ID: {case_id}")
        path_text = declared.get("path") if isinstance(declared, dict) else None
        if not isinstance(path_text, str):
            raise ScienceError(f"{case_id} acceptance output lacks a path")
        path = Path(path_text)
        if path.resolve(strict=True) != (
            acceptance_root / "cases" / str(case_id) / "case_acceptance.json"
        ).resolve(strict=True):
            raise ScienceError(f"{case_id} acceptance output escapes its acceptance directory")
        current = binding(path)
        require_same_binding(declared, current, f"{case_id} acceptance output")
        record = load_json(path, f"{case_id} direct-fast acceptance")
        if (
            record.get("record_type") != "cgl-lf-stage-i-direct-fast-case-acceptance"
            or record.get("case_id") != case_id
            or record.get("result") not in VALID_RESULTS
            or record.get("result") != case_results[case_id]
        ):
            raise ScienceError(f"{case_id} direct-fast acceptance identity differs")
        require_exact_science_scope_records(
            record, limitation, intervention, f"{case_id} direct-fast acceptance"
        )
        cases[str(case_id)] = record
        case_bindings[str(case_id)] = current
    return cases, campaign, provenance, case_bindings


def validate_reviewed_case_evidence(
    case_id: str,
    summary: dict[str, object],
    policy: dict[str, object],
    lineage: dict[str, object],
) -> dict[str, object] | None:
    """Load optional sealed reviewed case evidence retained by fast acceptance."""

    reviewed_record = summary.get("reviewed_case_evidence")
    declared = reviewed_record.get("binding") if isinstance(reviewed_record, dict) else None
    if declared is None:
        return None
    path_text = declared.get("path") if isinstance(declared, dict) else None
    if not isinstance(path_text, str):
        raise ScienceError(f"{case_id} reviewed case evidence binding is malformed")
    path = Path(path_text)
    require_same_binding(declared, binding(path), f"{case_id} reviewed case evidence")
    value = load_json(path, f"{case_id} reviewed case evidence")
    reviewed = load_reviewed_module()
    try:
        reviewed.verify_evidence_digest(value, f"{case_id} reviewed case evidence")
    except Exception as error:
        raise ScienceError(f"{case_id} reviewed case evidence is forged: {error}") from error
    if (
        value.get("record_type") != "stage-i-scientific-case-evidence"
        or value.get("case_id") != case_id
    ):
        raise ScienceError(f"{case_id} reviewed case evidence identity differs")
    provenance = value.get("provenance")
    if not isinstance(provenance, dict):
        raise ScienceError(f"{case_id} reviewed case evidence lacks provenance")
    verify_policy_provenance(
        provenance, policy, f"{case_id} reviewed case evidence", "acceptance_utility"
    )
    evaluation = value.get("evaluation_inputs")
    if not isinstance(evaluation, dict):
        raise ScienceError(f"{case_id} reviewed case evidence lacks evaluation inputs")
    for kind in ("mhd", "user"):
        current = history_binding_from_lineage(lineage, kind, case_id)
        require_same_binding(
            evaluation.get(f"{kind}_history"),
            current,
            f"{case_id} reviewed {kind} history",
        )
    verify_recursive_bindings(evaluation, f"{case_id} reviewed evaluation inputs")
    verify_recursive_bindings(provenance, f"{case_id} reviewed case provenance")
    return value


def validate_case_acceptance(
    case_id: str,
    summary: dict[str, object] | None,
    lineage: dict[str, object] | None,
    lineage_binding: dict[str, object] | None,
    policy: dict[str, object],
) -> tuple[bool, str, dict[str, object] | None]:
    """Authenticate one case summary and return its gate eligibility."""

    if summary is None:
        return False, "fast-acceptance case record is unavailable", None
    if lineage is None or lineage_binding is None:
        raise ScienceError(f"{case_id} acceptance exists without a selected report lineage")
    provenance = summary.get("provenance")
    if not isinstance(provenance, dict):
        raise ScienceError(f"{case_id} acceptance lacks provenance")
    require_same_binding(
        provenance.get("lineage"), lineage_binding, f"{case_id} acceptance lineage"
    )
    histories = provenance.get("histories")
    if not isinstance(histories, dict):
        raise ScienceError(f"{case_id} acceptance history provenance is malformed")
    for kind in ("mhd", "user"):
        current = history_binding_from_lineage(lineage, kind, case_id)
        require_same_binding(
            histories.get(kind), current, f"{case_id} acceptance {kind} history"
        )
    reviewed_evidence = validate_reviewed_case_evidence(case_id, summary, policy, lineage)
    reviewed_summary = summary.get("reviewed_case_evidence")
    reviewed_summary_result = (
        reviewed_summary.get("result") if isinstance(reviewed_summary, dict) else None
    )
    if (
        reviewed_evidence is not None
        and reviewed_evidence.get("result") != reviewed_summary_result
    ):
        raise ScienceError(f"{case_id} reviewed case result differs from its summary")
    health = summary.get("health")
    scope = summary.get("scope")
    comparison = summary.get("comparison_evidence")
    eligible = (
        summary.get("result") == "pass"
        and isinstance(health, dict)
        and health.get("result") == "pass"
        and isinstance(scope, dict)
        and scope.get("campaign_interpretation_eligible") is True
        and isinstance(comparison, dict)
        and comparison.get("result") == "pass"
        and reviewed_summary_result == "pass"
        and isinstance(reviewed_evidence, dict)
        and reviewed_evidence.get("result") == "pass"
    )
    if summary.get("result") == "pass" and not eligible:
        raise ScienceError(f"{case_id} acceptance pass has inconsistent supporting gates")
    reason = (
        "direct-fast case acceptance passed"
        if eligible
        else f"direct-fast case acceptance is {summary.get('result', 'unavailable')}"
    )
    return eligible, reason, reviewed_evidence


def snapshot_provenance_errors(diagnostics: dict[str, object], case_id: str) -> list[str]:
    """Return errors for raw snapshot bindings used by analyzed products."""

    errors: list[str] = []
    snapshots = diagnostics.get("snapshots")
    if not isinstance(snapshots, dict):
        return [f"{case_id} snapshot records are malformed"]
    ensemble = diagnostics.get("snapshot_ensemble")
    expected_count = (
        ensemble.get("snapshot_count") if isinstance(ensemble, dict) else None
    )
    if (
        not isinstance(expected_count, int)
        or isinstance(expected_count, bool)
        or expected_count <= 0
        or len(snapshots) != expected_count
    ):
        errors.append(f"{case_id} snapshot record count differs from its ensemble")
    for source, record in sorted(snapshots.items()):
        if not isinstance(record, dict):
            errors.append(f"{case_id} snapshot record is malformed: {source}")
            continue
        provenance = record.get("snapshot_provenance")
        if not isinstance(provenance, dict):
            errors.append(f"{case_id} snapshot lacks provenance: {source}")
            continue
        files = provenance.get("files")
        if not isinstance(files, list) or not files:
            errors.append(f"{case_id} snapshot provenance lacks files: {source}")
            continue
        if provenance.get("expected_rank_count") != len(files):
            errors.append(f"{case_id} snapshot rank count differs: {source}")
        try:
            observed_digest = hashlib.sha256(canonical_json(files)).hexdigest()
        except (TypeError, ValueError) as error:
            errors.append(f"{case_id} snapshot inventory is not canonical: {error}")
            continue
        if provenance.get("aggregate_sha256") != observed_digest:
            errors.append(f"{case_id} snapshot aggregate digest differs: {source}")
        try:
            verify_recursive_bindings(files, f"{case_id} snapshot {source}")
        except ScienceError as error:
            errors.append(str(error))
    return errors


def validate_diagnostics(
    case_id: str,
    report_root: Path,
    lineage: dict[str, object] | None,
    lineage_binding: dict[str, object] | None,
    inventory_adapter: dict[str, object],
) -> tuple[dict[str, object] | None, bool, str, dict[str, object] | None]:
    """Authenticate optional direct-fast case diagnostics and snapshot products."""

    path = report_root / "cases" / case_id / "diagnostics.json"
    if not path.is_file():
        return None, False, "case diagnostics are unavailable", None
    if lineage is None or lineage_binding is None:
        raise ScienceError(f"{case_id} diagnostics exist without an inventoried lineage")
    diagnostics = load_json(path, f"{case_id} diagnostics")
    if (
        diagnostics.get("schema_version") != 1
        or diagnostics.get("case_id") != case_id
        or diagnostics.get("case_name") != lineage.get("case_name")
    ):
        raise ScienceError(f"{case_id} diagnostics identity differs")
    provenance = diagnostics.get("provenance")
    if not isinstance(provenance, dict):
        raise ScienceError(f"{case_id} diagnostics lack provenance")
    require_same_binding(
        provenance.get("lineage"), lineage_binding, f"{case_id} diagnostics lineage"
    )
    require_same_binding(
        provenance.get("snapshot_index"),
        binding(report_root / "cases" / case_id / "snapshots.json"),
        f"{case_id} diagnostics snapshot index",
    )
    require_same_binding(
        provenance.get("adapter"), inventory_adapter, f"{case_id} diagnostics adapter"
    )
    require_same_binding(
        provenance.get("analyzer"), binding(PAPER_ANALYZER), f"{case_id} diagnostics analyzer"
    )
    for kind in ("mhd", "user"):
        current = history_binding_from_lineage(lineage, kind, case_id)
        require_same_binding(
            provenance.get(f"merged_{kind}_history"),
            current,
            f"{case_id} diagnostics merged {kind} history",
        )
    verify_recursive_bindings(provenance, f"{case_id} diagnostics provenance")
    ensemble = diagnostics.get("snapshot_ensemble")
    snapshot_ready = (
        diagnostics.get("snapshot_analysis_status") == "complete"
        and isinstance(ensemble, dict)
        and isinstance(ensemble.get("snapshot_count"), int)
        and int(ensemble["snapshot_count"]) > 0
    )
    if snapshot_ready:
        errors = snapshot_provenance_errors(diagnostics, case_id)
        if errors:
            raise ScienceError("; ".join(errors))
        compat = diagnostics.get("compat")
        compat_ensemble = (
            compat.get("snapshot_ensemble") if isinstance(compat, dict) else None
        )
        if canonical_json(compat_ensemble) != canonical_json(ensemble):
            raise ScienceError(f"{case_id} compat and primary snapshot ensembles differ")
        try:
            analyzer = load_report_module().load_pure_analyzer()
            recomputed = analyzer.average_snapshot_records(diagnostics["snapshots"])
            recomputed["time_start"] = ensemble.get("time_start")
            recomputed["time_end"] = ensemble.get("time_end")
            if "firehose_threshold_occupancy" in recomputed:
                recomputed["firehose_threshold_occupancy"]["analysis_window"].update({
                    "requested_time_start": ensemble.get("time_start"),
                    "requested_time_end": ensemble.get("time_end"),
                })
            if canonical_json(recomputed) != canonical_json(ensemble):
                raise ScienceError(
                    f"{case_id} snapshot ensemble differs from its retained records"
                )
        except ScienceError:
            raise
        except Exception as error:
            raise ScienceError(
                f"{case_id} snapshot ensemble replay failed: {error}"
            ) from error
        reason = "snapshot diagnostics and raw snapshot provenance authenticated"
    else:
        reason = "authenticated diagnostics lack complete snapshot analysis"
    return diagnostics, snapshot_ready, reason, binding(path)


def case_metrics_from_summary(summary: dict[str, object] | None) -> dict[str, object]:
    """Convert persisted fast-acceptance statistics to reviewed scalar schema."""

    result: dict[str, object] = {}
    values = summary.get("history_statistics") if isinstance(summary, dict) else None
    if not isinstance(values, dict):
        return result
    for metric, record in sorted(values.items()):
        windows = record.get("windows") if isinstance(record, dict) else None
        if not isinstance(windows, dict):
            continue
        selected: dict[str, object] = {}
        for window, window_record in sorted(windows.items()):
            statistics = (
                window_record.get("statistics")
                if isinstance(window_record, dict)
                else None
            )
            if isinstance(statistics, dict):
                selected[str(window)] = statistics
        if selected:
            result[str(metric)] = selected
    return result


def direct_convergence_products(
    diagnostics: dict[str, object] | None,
    snapshot_ready: bool,
    case_name: str | None,
) -> dict[str, object]:
    """Extract convergence curves from authenticated direct-fast products."""

    if not snapshot_ready or diagnostics is None or not isinstance(case_name, str):
        return {}
    reviewed = load_reviewed_module()
    report = load_report_module()
    ensemble = diagnostics.get("snapshot_ensemble")
    wrapped = {"cases": {case_name: {"snapshot_ensemble": ensemble}}}
    try:
        products = reviewed.convergence_products(wrapped, case_name)
        if isinstance(ensemble, dict) and isinstance(ensemble.get("alignment"), dict):
            analyzer = report.load_pure_analyzer()
            x, y = analyzer.alignment_peak_curve(ensemble)
            products["peak_alignment"] = reviewed.finite_curve(
                {"x": x.tolist(), "y": y.tolist()}, "direct peak alignment"
            )
        return products
    except Exception:
        return {}


def build_case_science(
    case_id: str,
    summary: dict[str, object] | None,
    reviewed_evidence: dict[str, object] | None,
    diagnostics: dict[str, object] | None,
    snapshot_ready: bool,
    case_name: str | None,
) -> dict[str, object]:
    """Build the normalized case object consumed by reviewed comparison kernels."""

    evidence = {
        "case_id": case_id,
        "metrics": case_metrics_from_summary(summary),
        "analyzer_metrics": {},
        "convergence_products": direct_convergence_products(
            diagnostics, snapshot_ready, case_name
        ),
        "panel_products": [],
        "gates": [],
    }
    if isinstance(reviewed_evidence, dict):
        for key in (
            "metrics",
            "analyzer_metrics",
            "convergence_products",
            "panel_products",
            "gates",
        ):
            value = reviewed_evidence.get(key)
            if isinstance(value, dict) and key not in ("panel_products", "gates"):
                evidence[key] = value
            elif isinstance(value, list) and key in ("panel_products", "gates"):
                evidence[key] = value
    return evidence


def scalar_metric_names(policy: dict[str, object]) -> list[str]:
    """Return the complete declared scalar comparison inventory."""

    case_metrics = policy["criteria"].get("case_metrics")
    names = set(str(value) for value in case_metrics) if isinstance(case_metrics, dict) else set()
    names.update(("unstable_occupancy", "peak_alignment"))
    return sorted(names)


def generic_contrast(
    reviewed: object,
    left_id: str,
    right_id: str,
    cases: dict[str, dict[str, object]],
    eligible: dict[str, bool],
    metrics: list[str],
) -> dict[str, object]:
    """Build a descriptive matched contrast with reviewed scalar extraction."""

    left_case = cases.get(left_id, {})
    right_case = cases.get(right_id, {})
    records: list[dict[str, object]] = []
    for metric in metrics:
        try:
            left = reviewed.scalar_from_case(left_case, metric)
            right = reviewed.scalar_from_case(right_case, metric)
        except Exception:
            left = right = None
        if left is None or right is None:
            records.append({
                "metric": metric,
                "available": False,
                "reason": "one or both exact-window scalar estimates are unavailable",
            })
            continue
        difference = float(left["mean"]) - float(right["mean"])
        combined_se = math.hypot(
            float(left["standard_error"]), float(right["standard_error"])
        )
        pooled = math.sqrt(
            0.5
            * (
                float(left["standard_deviation"]) ** 2
                + float(right["standard_deviation"]) ** 2
            )
        )
        records.append({
            "metric": metric,
            "available": True,
            "left_mean": float(left["mean"]),
            "right_mean": float(right["mean"]),
            "difference_left_minus_right": difference,
            "combined_standard_error": combined_se,
            "pooled_within_realization_standard_deviation": pooled,
            "standardized_effect": reviewed.finite_ratio(difference, pooled),
            "standardized_effect_scope": "descriptive_within_realization",
            "claim_scope": "descriptive_within_realization",
            "population_inference": reviewed.POPULATION_INFERENCE_LIMITATION,
        })
    both_eligible = eligible.get(left_id, False) and eligible.get(right_id, False)
    return {
        "left": left_id,
        "right": right_id,
        "result": "available" if any(record["available"] for record in records) else "inconclusive",
        "claim_eligible": both_eligible,
        "population_inference": reviewed.POPULATION_INFERENCE_LIMITATION,
        "reason": (
            "descriptive exact-window contrast; no preregistered family pass threshold"
            if both_eligible
            else "one or both contributing cases have not passed fast acceptance"
        ),
        "metrics": records,
    }


def active_passive_contrasts(
    reviewed: object,
    policy: dict[str, object],
    cases: dict[str, dict[str, object]],
    eligible: dict[str, bool],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Evaluate active/passive pairs with reviewed descriptive-direction criteria."""

    records: dict[str, object] = {}
    gates: list[dict[str, object]] = []
    intervention_scope = policy["criteria"]["family_gates"]["active_passive"].get(
        "intervention_scope"
    )
    if (
        not isinstance(intervention_scope, dict)
        or intervention_scope != reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
    ):
        raise ScienceError("active/passive intervention scope differs")
    for active, passive in ACTIVE_PASSIVE_PAIRS:
        try:
            contrast = reviewed.pair_contrast(
                cases.get(active, {}), cases.get(passive, {}), policy
            )
        except Exception as error:
            contrast = {
                "result": "inconclusive",
                "reason": f"reviewed pair contrast unavailable: {error}",
                "metrics": [],
                "intervention_scope": intervention_scope,
                "population_inference": reviewed.POPULATION_INFERENCE_LIMITATION,
            }
        if contrast.get("intervention_scope") != intervention_scope:
            raise ScienceError("reviewed pair contrast intervention scope differs")
        claim_eligible = eligible.get(active, False) and eligible.get(passive, False)
        provisional = str(contrast.get("result", "inconclusive"))
        if not claim_eligible:
            contrast["result"] = "inconclusive"
            contrast["reason"] = "one or both contributing cases have not passed fast acceptance"
        contrast["active"] = active
        contrast["passive"] = passive
        contrast["claim_eligible"] = claim_eligible
        contrast["reviewed_kernel_result_before_case_gate"] = provisional
        key = f"{active}_{passive}"
        records[key] = contrast
        gates.append(reviewed.gate(
            f"active_passive_pair:{active}:{passive}",
            str(contrast["result"]),
            reason=str(contrast["reason"]),
            observations=contrast,
        ))
    return records, gates


def limiter_ordering_gate(
    reviewed: object,
    cases: dict[str, dict[str, object]],
    eligible: dict[str, bool],
    summaries: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Evaluate the preregistered R15 greater-than R14 late-nu_eff ordering."""

    scope = summaries.get("R14", {}).get("scope")
    admitted = (
        isinstance(scope, dict)
        and scope.get("classification") == "scoped_nonfatal_hard_bound_variant"
    )
    if not admitted or not eligible.get("R14", False) or not eligible.get("R15", False):
        return reviewed.gate(
            "finite_limiter_ordering:R15_gt_R14",
            "inconclusive",
            reason="R14/R15 passing admitted fast-acceptance evidence is unavailable",
        )
    lower = reviewed.scalar_from_case(cases["R14"], "nu_eff", "late")
    upper = reviewed.scalar_from_case(cases["R15"], "nu_eff", "late")
    if lower is None or upper is None:
        return reviewed.gate(
            "finite_limiter_ordering:R15_gt_R14",
            "inconclusive",
            reason="late-window R14/R15 nu_eff estimates are unavailable",
        )
    difference = float(upper["mean"]) - float(lower["mean"])
    descriptive_block_lower_95 = difference - 1.96 * math.hypot(
        float(upper["standard_error"]), float(lower["standard_error"])
    )
    return reviewed.gate(
        "finite_limiter_ordering:R15_gt_R14",
        "pass" if descriptive_block_lower_95 > 0.0 else "fail",
        reason=(
            "descriptive within-trajectory block bound supports R15 greater than R14"
            if descriptive_block_lower_95 > 0.0
            else "finite-limiter descriptive ordering is unsupported"
        ),
        observations={
            "R15_minus_R14": difference,
            "descriptive_block_lower_95": descriptive_block_lower_95,
            "claim_scope": "descriptive_within_trajectory",
        },
    )


def lf_strength_gate(
    reviewed: object,
    policy: dict[str, object],
    cases: dict[str, dict[str, object]],
    eligible: dict[str, bool],
) -> dict[str, object]:
    """Report active LF-strength activity and descriptive response."""

    assessment = reviewed.lf_strength_assessment(policy, cases, eligible)
    return reviewed.gate(
        "lf_strength_descriptive_response_activity",
        str(assessment["result"]),
        reason=str(assessment["reason"]),
        observations=assessment["observations"],
        limits=assessment["limits"],
    )


def resolution_gate(
    reviewed: object,
    policy: dict[str, object],
    cases: dict[str, dict[str, object]],
    eligible: dict[str, bool],
) -> dict[str, object]:
    """Evaluate the declared R16/R02/R17 convergence criteria."""

    criteria = policy["criteria"]["resolution_convergence"]
    required = ("R16", "R02", "R17")
    interval_over_pi = tuple(float(value) for value in criteria["common_k_perp_over_pi"])
    interval = tuple(value * math.pi for value in interval_over_pi)
    records: list[dict[str, object]] = []
    product_specs = {
        "peak_alignment": ("alignment", float(criteria["alignment_max_abs_lte"])),
        "velocity_spectrum_shape": ("spectrum", float(criteria["spectral_log_rms_lte"])),
        "magnetic_fluctuation_spectrum_shape": (
            "spectrum",
            float(criteria["spectral_log_rms_lte"]),
        ),
    }
    for product, (kind, limit) in product_specs.items():
        curves = [
            cases.get(case_id, {}).get("convergence_products", {}).get(product)
            if isinstance(cases.get(case_id, {}).get("convergence_products"), dict)
            else None
            for case_id in required
        ]
        if not all(isinstance(curve, dict) for curve in curves):
            records.append({
                "kind": "curve",
                "product": product,
                "available": False,
                "reason": "one or more authenticated convergence curves are unavailable",
            })
            continue
        try:
            samples = (
                [float(value) * math.pi for value in criteria["alignment_shells"]]
                if kind == "alignment"
                else None
            )
            low_mid = reviewed.curve_distance_resolution(
                curves[0],
                curves[1],
                interval,
                kind,
                samples,
                float(criteria["resolved_disagreement_sigma_gt"]),
            )
            mid_high = reviewed.curve_distance(
                curves[1], curves[2], interval, kind, samples
            )
        except Exception as error:
            records.append({
                "kind": "curve",
                "product": product,
                "available": False,
                "reason": str(error),
            })
            continue
        improved = (
            mid_high
            <= float(criteria["improvement_ratio_lte"]) * float(low_mid["distance"])
            if low_mid["resolved"]
            else True
        )
        records.append({
            "kind": "curve",
            "product": product,
            "available": True,
            "common_k_perp_over_pi": list(interval_over_pi),
            "R16_R02_distance": low_mid["distance"],
            "R16_R02_distance_standard_error": low_mid["distance_standard_error"],
            "resolved_R16_R02_difference": low_mid["resolved"],
            "R02_R17_distance": mid_high,
            "R02_R17_limit": limit,
            "improvement_ratio_lte": criteria["improvement_ratio_lte"],
            "improvement_applicability": low_mid["resolution_applicability"],
            "improved": improved,
            "passed": mid_high <= limit and improved,
        })
    for metric in criteria["scalar_metrics"]:
        values = [reviewed.scalar_from_case(cases.get(case_id, {}), metric) for case_id in required]
        if any(value is None for value in values):
            records.append({
                "kind": "scalar",
                "metric": metric,
                "available": False,
                "reason": "one or more exact-window scalar estimates are unavailable",
            })
            continue
        low, mid, high = values
        mid_high = abs(float(high["mean"]) - float(mid["mean"]))
        relative = reviewed.finite_ratio(mid_high, abs(float(mid["mean"])))
        within_uncertainty = mid_high <= 2.0 * math.hypot(
            float(high["standard_error"]), float(mid["standard_error"])
        )
        low_mid = abs(float(mid["mean"]) - float(low["mean"]))
        resolved_low_mid = low_mid > 2.0 * math.hypot(
            float(mid["standard_error"]), float(low["standard_error"])
        )
        improved = (
            not resolved_low_mid
            or mid_high <= float(criteria["improvement_ratio_lte"]) * low_mid
        )
        records.append({
            "kind": "scalar",
            "metric": metric,
            "available": True,
            "relative_R02_R17_difference": relative,
            "within_two_combined_standard_errors": within_uncertainty,
            "resolved_R16_R02_difference": resolved_low_mid,
            "improved": improved,
            "passed": (
                (
                    relative <= float(criteria["scalar_relative_difference_lte"])
                    or within_uncertainty
                )
                and improved
            ),
        })
    all_eligible = all(eligible.get(case_id, False) for case_id in required)
    if not all_eligible or any(record.get("available") is not True for record in records):
        result = "inconclusive"
        reason = "passing cases or required convergence products are unavailable"
    else:
        result = "pass" if all(record.get("passed") is True for record in records) else "fail"
        reason = (
            "resolution convergence gates passed"
            if result == "pass"
            else "one or more resolution convergence criteria failed"
        )
    return reviewed.gate(
        "R16_R02_R17_resolution_convergence",
        result,
        reason=reason,
        observations=records,
        limits=criteria,
    )


def recompute_direct_reference_products(
    policy: dict[str, object],
    diagnostics: dict[str, dict[str, object]],
    snapshot_ready: dict[str, bool],
) -> tuple[dict[str, dict[str, dict[str, object]]], str | None]:
    """Recompute full, early, and late MKS24 products from authenticated records."""

    report = load_report_module()
    analyzer = report.load_pure_analyzer()
    windows = policy["criteria"].get("analysis_windows")
    if not isinstance(windows, dict):
        return {}, "reviewed policy lacks analysis windows"
    compat_by_window: dict[str, dict[str, object]] = {
        str(window): {} for window in windows
    }
    for case_id, value in sorted(diagnostics.items()):
        if not snapshot_ready.get(case_id, False):
            continue
        compat = value.get("compat")
        name = value.get("case_name")
        records = value.get("snapshots")
        if (
            not isinstance(compat, dict)
            or not isinstance(name, str)
            or not isinstance(records, dict)
        ):
            continue
        for window, bounds in sorted(windows.items()):
            if (
                not isinstance(bounds, list)
                or len(bounds) != 2
                or not all(isinstance(bound, (int, float)) for bound in bounds)
            ):
                return {}, f"reviewed analysis window {window} is malformed"
            start, end = (float(bounds[0]), float(bounds[1]))
            selected = {
                str(path): record
                for path, record in sorted(records.items())
                if isinstance(record, dict)
                and isinstance(record.get("time"), (int, float))
                and start - WINDOW_TIME_TOLERANCE
                <= float(record["time"])
                <= end + WINDOW_TIME_TOLERANCE
            }
            if not selected:
                continue
            ensemble = analyzer.average_snapshot_records(selected)
            ensemble["time_start"] = start
            ensemble["time_end"] = end
            occupancy = ensemble.get("firehose_threshold_occupancy")
            if isinstance(occupancy, dict):
                analysis_window = occupancy.get("analysis_window")
                if isinstance(analysis_window, dict):
                    analysis_window.update({
                        "requested_time_start": start,
                        "requested_time_end": end,
                    })
            window_compat = dict(compat)
            window_compat["analysis_window"] = {
                "time_start": start,
                "time_end": end,
            }
            window_compat["snapshot_ensemble"] = ensemble
            compat_by_window[str(window)][name] = window_compat
    if not any(compat_by_window.values()):
        return {}, "no authenticated complete snapshot analyses are available"
    try:
        manifest_path = Path(str(policy["verified_sources"]["stage_i_manifest"]["path"]))
        configuration = analyzer.stage_i_panels_configuration(manifest_path)
        reviewed = load_reviewed_module()
        products: dict[str, dict[str, dict[str, object]]] = {}
        for window, compat_cases in sorted(compat_by_window.items()):
            if not compat_cases:
                continue
            reference = analyzer.combined_reference_curve_comparisons(
                {"cases": compat_cases},
                report.reference_manifest_paths(policy["manifest"]),
                allow_missing_cases=True,
                analysis_case_aliases=configuration["analysis_case_aliases"],
                stage_i_reference_bindings=configuration[
                    "reference_product_bindings"
                ],
            )
            products[window] = reviewed.reference_comparison_records(
                {"reference_curve_comparisons": reference}
            )
        return products, None
    except Exception as error:
        return {}, f"direct MKS24 recomputation unavailable: {type(error).__name__}: {error}"


def mks24_assessment(
    reviewed: object,
    policy: dict[str, object],
    cases: dict[str, dict[str, object]],
    diagnostics: dict[str, dict[str, object]],
    snapshot_ready: dict[str, bool],
    eligible: dict[str, bool],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Assess admitted MKS24 products from reviewed or direct authenticated evidence."""

    direct_windows, direct_error = recompute_direct_reference_products(
        policy, diagnostics, snapshot_ready
    )
    direct = direct_windows.get("full", {})
    direct_early = direct_windows.get("early", {})
    direct_late = direct_windows.get("late", {})
    panel_criteria = {
        str(item["id"]): item for item in policy["criteria"]["comparison_panels"]
    }
    manifest_status = policy["manifest"]["panel_status"]
    manifest_panels = {
        str(item["id"]): item
        for item in manifest_status["panels"]
        if item.get("disposition") == "comparison"
    }
    bindings = manifest_status["reference_product_bindings"]
    aliases = manifest_status["analysis_case_aliases"]
    name_to_case = {
        str(item["name"]): str(item["id"]) for item in policy["manifest"]["cases"]
    }
    reviewed_products: dict[str, tuple[str, dict[str, object]]] = {}
    for case_id, case in cases.items():
        values = case.get("panel_products")
        if not isinstance(values, list):
            continue
        for value in values:
            if isinstance(value, dict) and isinstance(value.get("product_id"), str):
                reviewed_products[str(value["product_id"])] = (case_id, value)
    panels: dict[str, object] = {}
    gates: list[dict[str, object]] = []
    for panel_id, criteria in sorted(panel_criteria.items()):
        expected = [str(value) for value in manifest_panels[panel_id]["reference_products"]]
        products: list[dict[str, object]] = []
        for product_id in expected:
            binding_record = bindings[product_id]
            expected_name = str(aliases.get(binding_record["case"], binding_record["case"]))
            case_id = name_to_case.get(expected_name)
            if product_id in reviewed_products:
                source_case, record = reviewed_products[product_id]
                result = str(record.get("result", "inconclusive"))
                if not eligible.get(source_case, False):
                    result = "inconclusive"
                products.append({
                    "panel_id": panel_id,
                    "product_id": product_id,
                    "case_id": source_case,
                    "source": "sealed_reviewed_case_evidence",
                    "result": result,
                    "reason": (
                        str(record.get("reason"))
                        if eligible.get(source_case, False)
                        else "contributing case has not passed fast acceptance"
                    ),
                    "observations": record.get("observations"),
                    "limits": record.get("limits"),
                })
                continue
            record = direct.get(product_id)
            if not isinstance(record, dict):
                products.append({
                    "panel_id": panel_id,
                    "product_id": product_id,
                    "case_id": case_id,
                    "source": "unavailable",
                    "result": "inconclusive",
                    "reason": (
                        direct_error
                        or "authenticated direct comparison product is unavailable"
                    ),
                })
                continue
            try:
                rms, maximum, _ = reviewed.normalized_reference_metrics(
                    policy, product_id, record, binding_record, expected_name
                )
            except Exception as error:
                products.append({
                    "panel_id": panel_id,
                    "product_id": product_id,
                    "case_id": case_id,
                    "source": "authenticated_direct_fast_recomputation",
                    "result": (
                        "fail"
                        if case_id and eligible.get(case_id, False)
                        else "inconclusive"
                    ),
                    "reason": f"reference product validation failed: {error}",
                })
                continue
            limits = {
                "normalized_residual_rms_lte": criteria["normalized_residual_rms_lte"],
                "maximum_absolute_normalized_residual_lte": criteria[
                    "maximum_absolute_normalized_residual_lte"
                ],
                "early_late_vector_drift_rms_lte": criteria[
                    "early_late_vector_drift_rms_lte"
                ],
            }
            residual_pass = (
                rms <= float(limits["normalized_residual_rms_lte"])
                and maximum <= float(limits["maximum_absolute_normalized_residual_lte"])
            )
            drift: float | None = None
            early_record = direct_early.get(product_id)
            late_record = direct_late.get(product_id)
            if (
                isinstance(early_record, dict)
                and isinstance(late_record, dict)
                and early_record.get("available") is True
                and late_record.get("available") is True
            ):
                try:
                    _, _, early_values = reviewed.normalized_reference_metrics(
                        policy, product_id, early_record, binding_record, expected_name
                    )
                    _, _, late_values = reviewed.normalized_reference_metrics(
                        policy, product_id, late_record, binding_record, expected_name
                    )
                    uncertainty = record.get("reference_y_uncertainty")
                    if binding_record.get("kind") == "surface":
                        uncertainty = record.get("reference_z_uncertainty")
                    if (
                        not isinstance(uncertainty, list)
                        or not early_values
                        or not (
                            len(early_values)
                            == len(late_values)
                            == len(uncertainty)
                        )
                    ):
                        raise ValueError("early/late reference vector lengths differ")
                    drift = math.sqrt(sum(
                        (
                            (float(early) - float(late)) / float(error)
                        ) ** 2
                        for early, late, error in zip(
                            early_values, late_values, uncertainty
                        )
                    ) / len(early_values))
                    if not math.isfinite(drift) or any(
                        not math.isfinite(float(error)) or float(error) <= 0.0
                        for error in uncertainty
                    ):
                        raise ValueError(
                            "early/late drift has nonfinite values or uncertainty"
                        )
                except Exception as error:
                    products.append({
                        "panel_id": panel_id,
                        "product_id": product_id,
                        "case_id": case_id,
                        "source": "authenticated_direct_fast_recomputation",
                        "result": (
                            "fail"
                            if case_id and eligible.get(case_id, False)
                            else "inconclusive"
                        ),
                        "reason": f"early/late reference product validation failed: {error}",
                    })
                    continue
            if not case_id or not eligible.get(case_id, False):
                result = "inconclusive"
                reason = "contributing case has not passed fast acceptance"
            elif not residual_pass:
                result = "fail"
                reason = "available normalized residual criterion failed"
            elif drift is None:
                result = "inconclusive"
                reason = "early/late panel vectors are unavailable"
            elif drift > float(limits["early_late_vector_drift_rms_lte"]):
                result = "fail"
                reason = "early/late panel stationarity criterion failed"
            else:
                result = "pass"
                reason = "reference residual and panel stationarity gates passed"
            products.append({
                "panel_id": panel_id,
                "product_id": product_id,
                "case_id": case_id,
                "source": "authenticated_direct_fast_recomputation",
                "result": result,
                "reason": reason,
                "observations": {
                    "normalized_residual_rms": rms,
                    "maximum_absolute_normalized_residual": maximum,
                    "early_late_vector_drift_rms": drift,
                },
                "limits": limits,
            })
        product_results = [str(item["result"]) for item in products]
        result = (
            "fail"
            if "fail" in product_results
            else "pass"
            if products and all(value == "pass" for value in product_results)
            else "inconclusive"
        )
        reason = (
            "all admitted panel products passed"
            if result == "pass"
            else "one or more admitted panel products failed"
            if result == "fail"
            else "one or more admitted panel products are unavailable or inconclusive"
        )
        panels[panel_id] = {
            "result": result,
            "reason": reason,
            "products": products,
        }
        gates.append(reviewed.gate(
            f"mks24_panel:{panel_id}",
            result,
            reason=reason,
            observations=products,
            limits=criteria,
        ))
    return {"result": reviewed.aggregate_gate_result(gates), "panels": panels}, gates


def aggregate_science(
    inventory_path: Path,
    acceptance_root: Path,
    output: Path,
    criteria: Path,
    criteria_review: Path,
    selected_values: list[str] | None,
) -> tuple[dict[str, object], list[Path]]:
    """Build and write the complete deterministic science aggregation."""

    reviewed = load_reviewed_module()
    try:
        policy = reviewed.load_validated_policy(criteria, criteria_review)
    except Exception as error:
        raise ScienceError(f"cannot load reviewed scientific policy: {error}") from error
    limitation, intervention = exact_science_scope_records(reviewed, policy)
    required = [str(value) for value in policy["criteria"]["required_cases"]]
    selected = parse_case_selection(selected_values, required)
    inventory, report_root, lineages, lineage_bindings = validate_inventory(
        inventory_path, policy, selected
    )
    summaries, _campaign, acceptance_provenance, summary_bindings = validate_acceptance_root(
        acceptance_root, inventory_path, policy
    )
    require_separate_output(output, report_root, acceptance_root)

    eligible: dict[str, bool] = {}
    diagnostics: dict[str, dict[str, object]] = {}
    diagnostic_bindings: dict[str, dict[str, object]] = {}
    snapshot_ready: dict[str, bool] = {}
    case_science: dict[str, dict[str, object]] = {}
    case_rows: list[dict[str, object]] = []
    for case_id in selected:
        summary = summaries.get(case_id)
        lineage = lineages.get(case_id)
        lineage_binding = lineage_bindings.get(case_id)
        is_eligible, reason, reviewed_case = validate_case_acceptance(
            case_id, summary, lineage, lineage_binding, policy
        )
        diagnostic, snapshots_ok, diagnostic_reason, diagnostic_binding = validate_diagnostics(
            case_id,
            report_root,
            lineage,
            lineage_binding,
            inventory["_science_verified_adapter"],
        )
        eligible[case_id] = is_eligible
        snapshot_ready[case_id] = snapshots_ok
        if diagnostic is not None:
            diagnostics[case_id] = diagnostic
        if diagnostic_binding is not None:
            diagnostic_bindings[case_id] = diagnostic_binding
        case_name = (
            str(lineage.get("case_name"))
            if isinstance(lineage, dict) and isinstance(lineage.get("case_name"), str)
            else None
        )
        case_science[case_id] = build_case_science(
            case_id, summary, reviewed_case, diagnostic, snapshots_ok, case_name
        )
        case_rows.append({
            "case_id": case_id,
            "case_name": case_name,
            "inventory_status": lineage.get("status") if isinstance(lineage, dict) else None,
            "acceptance_result": summary.get("result") if isinstance(summary, dict) else None,
            "claim_eligible": is_eligible,
            "eligibility_reason": reason,
            "diagnostics_available": diagnostic is not None,
            "snapshot_products_authenticated": snapshots_ok,
            "diagnostics_reason": diagnostic_reason,
        })

    active_passive, active_gates = active_passive_contrasts(
        reviewed, policy, case_science, eligible
    )
    metric_names = scalar_metric_names(policy)
    descriptive_families: dict[str, object] = {}
    for family, pairs in CONTRAST_FAMILIES.items():
        descriptive_families[family] = {
            f"{left}_{right}": generic_contrast(
                reviewed, left, right, case_science, eligible, metric_names
            )
            for left, right in pairs
        }
    limiter_gate = limiter_ordering_gate(reviewed, case_science, eligible, summaries)
    lf_gate = lf_strength_gate(reviewed, policy, case_science, eligible)
    convergence = resolution_gate(reviewed, policy, case_science, eligible)
    mks24, mks24_gates = mks24_assessment(
        reviewed, policy, case_science, diagnostics, snapshot_ready, eligible
    )
    gates = [*active_gates, limiter_gate, lf_gate, convergence, *mks24_gates]
    result = reviewed.aggregate_gate_result(gates)
    case_dispositions = {
        row["case_id"]: {
            key: value for key, value in row.items() if key != "case_id"
        }
        for row in case_rows
    }
    science = reviewed.seal_evidence({
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-comparisons",
        "authority": "non-authorizing-direct-fast-scientific-assessment",
        "release_authorizing": False,
        "result": result,
        "selected_cases": selected,
        "case_dispositions": case_dispositions,
        "active_passive_intervention_scope": intervention,
        "current_science_scope_limitation": limitation,
        "families": {
            "active_passive": active_passive,
            **descriptive_families,
        },
        "resolution": convergence,
        "mks24": mks24,
        "gates": gates,
        "provenance": {
            "inventory": binding(inventory_path),
            "acceptance_provenance": binding(acceptance_root / "provenance.json"),
            "acceptance_campaign_evidence": binding(
                Path(str(acceptance_provenance["outputs"]["campaign_evidence"]["path"]))
            ),
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "reviewed_acceptance_utility": binding(REVIEWED_UTILITY),
            "fast_report_utility": binding(FAST_REPORT_UTILITY),
            "paper_analyzer": binding(PAPER_ANALYZER),
            "aggregator": binding(Path(__file__)),
            "case_acceptance": {
                case_id: summary_bindings[case_id]
                for case_id in selected
                if case_id in summary_bindings
            },
            "case_lineages": {
                case_id: lineage_bindings[case_id]
                for case_id in selected
                if case_id in lineage_bindings
            },
            "case_diagnostics": diagnostic_bindings,
        },
    })
    output.mkdir(parents=True, exist_ok=True)
    science_path = output / "science.json"
    write_json(science_path, science)
    table_paths = write_tables(output, case_rows, science)
    provenance = {
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-provenance",
        "inputs": science["provenance"],
        "outputs": {
            "science": binding(science_path),
            "tables": [binding(path) for path in table_paths],
        },
    }
    write_json(output / "provenance.json", provenance)
    return science, table_paths


def table_value(value: object) -> str:
    """Format one deterministic table cell."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.10g}" if math.isfinite(value) else ""
    if isinstance(value, (list, tuple)):
        return "; ".join(table_value(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(json_safe(value), sort_keys=True, separators=(",", ":"))
    return str(value)


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    """Write deterministic CSV."""

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: table_value(row.get(field)) for field in fields})
    write_text(path, stream.getvalue())


def write_markdown_table(
    path: Path,
    title: str,
    fields: list[str],
    rows: Iterable[dict[str, object]],
) -> None:
    """Write deterministic Markdown."""

    lines = [
        f"# {title}",
        "",
        "| " + " | ".join(fields) + " |",
        "| " + " | ".join("---" for _ in fields) + " |",
    ]
    for row in rows:
        cells = [
            table_value(row.get(field)).replace("|", "\\|").replace("\n", " ")
            for field in fields
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    write_text(path, "\n".join(lines))


def contrast_rows(science: dict[str, object]) -> list[dict[str, object]]:
    """Flatten all contrast families."""

    rows: list[dict[str, object]] = []
    families = science.get("families")
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
            intervention_scope = contrast.get("intervention_scope")
            if not isinstance(intervention_scope, dict):
                intervention_scope = {}
            intervention_fields = {
                "intervention_estimand": intervention_scope.get("estimand"),
                "intervention_enabled_components": intervention_scope.get(
                    "enabled_components"
                ),
                "excluded_interpretation": intervention_scope.get(
                    "excluded_interpretation"
                ),
                "intervention_declaration": intervention_scope.get("declaration"),
            }
            metrics = contrast.get("metrics")
            if not isinstance(metrics, list) or not metrics:
                rows.append({
                    "family": family,
                    "contrast": name,
                    "left": left,
                    "right": right,
                    "result": contrast.get("result"),
                    "claim_eligible": contrast.get("claim_eligible"),
                    "reason": contrast.get("reason"),
                    **intervention_fields,
                })
                continue
            for metric in metrics:
                if not isinstance(metric, dict):
                    continue
                rows.append({
                    "family": family,
                    "contrast": name,
                    "left": left,
                    "right": right,
                    "result": contrast.get("result"),
                    "claim_eligible": contrast.get("claim_eligible"),
                    "metric": metric.get("metric"),
                    "available": metric.get("available"),
                    "left_mean": metric.get("left_mean"),
                    "right_mean": metric.get("right_mean"),
                    "difference": metric.get(
                        "difference", metric.get("difference_left_minus_right")
                    ),
                    "combined_standard_error": metric.get("combined_standard_error"),
                    "standardized_effect": metric.get("standardized_effect"),
                    "pooled_within_realization_standard_deviation": metric.get(
                        "pooled_within_realization_standard_deviation"
                    ),
                    "expected_direction": metric.get("expected_direction"),
                    "direction_coherent": metric.get("direction_coherent"),
                    "large_direction_coherent_effect": metric.get(
                        "large_direction_coherent_effect"
                    ),
                    "claim_scope": metric.get("claim_scope"),
                    "reason": metric.get("reason", contrast.get("reason")),
                    **intervention_fields,
                })
    return rows


def gate_rows(science: dict[str, object]) -> list[dict[str, object]]:
    """Flatten science gates."""

    result: list[dict[str, object]] = []
    for gate in science.get("gates", []):
        if isinstance(gate, dict):
            result.append({
                "gate": gate.get("name"),
                "result": gate.get("result"),
                "reason": gate.get("reason"),
                "observations": gate.get("observations"),
                "limits": gate.get("limits"),
            })
    return result


def resolution_rows(science: dict[str, object]) -> list[dict[str, object]]:
    """Flatten resolution observations."""

    rows: list[dict[str, object]] = []
    resolution = science.get("resolution")
    values = resolution.get("observations") if isinstance(resolution, dict) else None
    if not isinstance(values, list):
        return rows
    for value in values:
        if isinstance(value, dict):
            rows.append({
                "result": resolution.get("result"),
                "kind": value.get("kind"),
                "name": value.get("product", value.get("metric")),
                "available": value.get("available"),
                "passed": value.get("passed"),
                "observations": value,
            })
    return rows


def mks24_rows(science: dict[str, object]) -> list[dict[str, object]]:
    """Flatten admitted MKS24 products."""

    rows: list[dict[str, object]] = []
    mks24 = science.get("mks24")
    panels = mks24.get("panels") if isinstance(mks24, dict) else None
    if not isinstance(panels, dict):
        return rows
    for panel_id, panel in sorted(panels.items()):
        products = panel.get("products") if isinstance(panel, dict) else None
        if not isinstance(products, list):
            continue
        for product in products:
            if isinstance(product, dict):
                observations = product.get("observations")
                limits = product.get("limits")
                rows.append({
                    "panel": panel_id,
                    "panel_result": panel.get("result"),
                    "product_id": product.get("product_id"),
                    "case_id": product.get("case_id"),
                    "source": product.get("source"),
                    "result": product.get("result"),
                    "normalized_residual_rms": (
                        observations.get("normalized_residual_rms")
                        if isinstance(observations, dict)
                        else None
                    ),
                    "maximum_absolute_normalized_residual": (
                        observations.get("maximum_absolute_normalized_residual")
                        if isinstance(observations, dict)
                        else None
                    ),
                    "early_late_vector_drift_rms": (
                        observations.get("early_late_vector_drift_rms")
                        if isinstance(observations, dict)
                        else None
                    ),
                    "limits": limits,
                    "reason": product.get("reason"),
                })
    return rows


def write_tables(
    output: Path,
    cases: list[dict[str, object]],
    science: dict[str, object],
) -> list[Path]:
    """Write deterministic CSV and Markdown comparison tables."""

    tables = output / "tables"
    specifications = [
        (
            "cases",
            "Case Science Eligibility",
            [
                "case_id", "case_name", "inventory_status", "acceptance_result",
                "claim_eligible", "eligibility_reason", "diagnostics_available",
                "snapshot_products_authenticated", "diagnostics_reason",
            ],
            cases,
        ),
        (
            "contrasts",
            "Final Scientific Contrasts",
            [
                "family", "contrast", "left", "right", "result", "claim_eligible",
                "metric", "available", "left_mean", "right_mean", "difference",
                "combined_standard_error", "pooled_within_realization_standard_deviation",
                "standardized_effect", "expected_direction", "direction_coherent",
                "large_direction_coherent_effect", "claim_scope",
                "intervention_estimand", "intervention_enabled_components",
                "excluded_interpretation", "intervention_declaration",
                "reason",
            ],
            contrast_rows(science),
        ),
        (
            "gates",
            "Final Scientific Gates",
            ["gate", "result", "reason", "observations", "limits"],
            gate_rows(science),
        ),
        (
            "resolution",
            "Resolution Criteria",
            ["result", "kind", "name", "available", "passed", "observations"],
            resolution_rows(science),
        ),
        (
            "mks24",
            "Admitted MKS24 Criteria",
            [
                "panel", "panel_result", "product_id", "case_id", "source", "result",
                "normalized_residual_rms", "maximum_absolute_normalized_residual",
                "early_late_vector_drift_rms", "limits", "reason",
            ],
            mks24_rows(science),
        ),
    ]
    paths: list[Path] = []
    for stem, title, fields, rows in specifications:
        csv_path = tables / f"{stem}.csv"
        md_path = tables / f"{stem}.md"
        write_csv(csv_path, fields, rows)
        write_markdown_table(md_path, title, fields, rows)
        paths.extend((csv_path, md_path))
    return paths


def output_default(inventory_path: Path) -> Path:
    """Return the deterministic sibling science-output default."""

    parent = inventory_path.expanduser().absolute().parent
    return parent.with_name(f"{parent.name}-science")


def acceptance_default(inventory_path: Path) -> Path:
    """Return the deterministic sibling fast-acceptance default."""

    parent = inventory_path.expanduser().absolute().parent
    return parent.with_name(f"{parent.name}-acceptance")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument(
        "--acceptance",
        type=Path,
        help="fast-acceptance output directory; defaults to INVENTORY_PARENT-acceptance",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="separate science output directory; defaults to INVENTORY_PARENT-science",
    )
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--criteria-review", type=Path, default=DEFAULT_CRITERIA_REVIEW)
    parser.add_argument(
        "--cases",
        action="append",
        help="optional repeated comma-separated R02-R17 selection",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the reviewed direct-fast science aggregation."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        inventory = args.inventory.expanduser().absolute().resolve(strict=True)
        acceptance = (
            args.acceptance.expanduser().absolute().resolve(strict=True)
            if args.acceptance is not None
            else acceptance_default(inventory).resolve(strict=True)
        )
        output = (
            args.output.expanduser().absolute()
            if args.output is not None
            else output_default(inventory)
        ).resolve(strict=False)
        science, _ = aggregate_science(
            inventory,
            acceptance,
            output,
            args.criteria.expanduser().absolute().resolve(strict=True),
            args.criteria_review.expanduser().absolute().resolve(strict=True),
            args.cases,
        )
        eligible = sum(
            value.get("claim_eligible") is True
            for value in science["case_dispositions"].values()
        )
        print(
            f"science: result={science['result']} "
            f"claim_eligible_cases={eligible}/{len(science['selected_cases'])} "
            f"output={output}"
        )
        return 0
    except (ScienceError, OSError, ValueError, KeyError, TypeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
