#!/usr/bin/env python3
"""Build lean scientific-acceptance evidence from direct-fast Stage I reports.

This driver is read-only with respect to simulations and fast-report outputs.
It imports the reviewed Stage I scientific-acceptance utility for policy
validation, exact-window statistics, complete-case evaluation, pair contrasts,
and evidence sealing.  Partial cases are retained as inconclusive evidence
instead of aborting the campaign assessment.

R10 is always exploratory-only because the qualified executable regularizes a
non-hyperbolic CGL fast-speed discriminant.  R14 is admitted only as the
explicit finite-limiter, nonfatal-hard-bound variant when that variant is
recorded in the selected fast-report lineage.
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
DEFAULT_CAMPAIGN_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
DEFAULT_INVENTORY = DEFAULT_CAMPAIGN_ROOT / (
    "analysis/mks24-stage-i-fast/E03-forcing-policy/R02-R17-final/inventory.json"
)
ACCEPTANCE_UTILITY = (
    REPO_ROOT / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
)
WORKFLOW_UTILITY = REPO_ROOT / "scripts/cgl_lf_workflow.py"
DEFAULT_CRITERIA = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
DEFAULT_CRITERIA_REVIEW = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)
CASE_ID_PATTERN = re.compile(r"R(?:0[2-9]|1[0-7])")
TIME_TOLERANCE = 1.0e-12
FATAL_COUNTERS = ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos")
HARD_BOUND_COUNTER = "lf_hardbd"
R14_NONFATAL_VARIANT = "finite_limiter_hard_bound_diagnostic_nonfatal"
R14_OVERRIDE = "mhd/cgl_lf_strict_admissibility=false"


class FastAcceptanceError(RuntimeError):
    """Raised for a malformed driver invocation or fast-report inventory."""


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a regular file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    """Write UTF-8 text atomically."""

    atomic_write(path, value.encode("utf-8"))


def load_json(path: Path) -> dict[str, object]:
    """Load one JSON object."""

    try:
        with path.open(encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise FastAcceptanceError(f"cannot load JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise FastAcceptanceError(f"expected a JSON object: {path}")
    return value


def binding(path: Path) -> dict[str, object]:
    """Bind one existing regular file without nondeterministic metadata."""

    resolved = path.expanduser().absolute().resolve(strict=True)
    if not resolved.is_file():
        raise FastAcceptanceError(f"expected a regular file: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def load_acceptance_module() -> object:
    """Import the reviewed scientific-acceptance utility by exact path."""

    module_name = "_cgl_lf_stage_i_reviewed_acceptance_for_direct_fast"
    existing = sys.modules.get(module_name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(module_name, ACCEPTANCE_UTILITY)
    if spec is None or spec.loader is None:
        raise FastAcceptanceError(
            f"cannot import reviewed acceptance utility: {ACCEPTANCE_UTILITY}"
        )
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def load_workflow_module() -> object:
    """Import the input-policy parser used to derive effective model choices."""

    module_name = "_cgl_lf_stage_i_workflow_for_direct_fast_acceptance"
    existing = sys.modules.get(module_name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(module_name, WORKFLOW_UTILITY)
    if spec is None or spec.loader is None:
        raise FastAcceptanceError(f"cannot import workflow utility: {WORKFLOW_UTILITY}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def parse_case_selection(values: list[str] | None) -> list[str] | None:
    """Parse repeated comma-separated case selections."""

    if not values:
        return None
    selected: set[str] = set()
    for value in values:
        for token in value.replace(",", " ").split():
            if CASE_ID_PATTERN.fullmatch(token) is None:
                raise FastAcceptanceError(f"invalid Stage I case ID: {token}")
            selected.add(token)
    return sorted(selected)


def inventory_output_root(inventory_path: Path, inventory: dict[str, object]) -> Path:
    """Return the fast-report output root declared by an inventory."""

    output = inventory.get("output")
    if isinstance(output, str):
        return Path(output).expanduser().resolve(strict=True)
    return inventory_path.expanduser().resolve(strict=True).parent


def path_contains(parent: Path, child: Path) -> bool:
    """Return whether child is parent or is located below parent."""

    return child == parent or parent in child.parents


def reject_nested_output_symlink_escapes(output: Path) -> None:
    """Reject pre-existing nested symlinks that resolve outside the output root."""

    if not output.exists():
        return
    root = output.resolve(strict=True)
    for directory, names, files in os.walk(root, followlinks=False):
        parent = Path(directory)
        for name in [*names, *files]:
            candidate = parent / name
            if not candidate.is_symlink():
                continue
            resolved = candidate.resolve(strict=False)
            if not path_contains(root, resolved):
                raise FastAcceptanceError(
                    "acceptance output contains a nested symlink escape: "
                    f"{candidate} -> {resolved}"
                )


def discover_case_ids(
    output_root: Path,
    inventory: dict[str, object],
    requested: list[str] | None,
    inventory_only: bool,
) -> list[str]:
    """Discover requested or available fast-report case records."""

    if requested is not None:
        return requested
    cases = inventory.get("cases")
    found = {
        str(case_id)
        for case_id in cases
        if isinstance(cases, dict) and CASE_ID_PATTERN.fullmatch(str(case_id))
    } if isinstance(cases, dict) else set()
    if not inventory_only:
        for path in (output_root / "cases").glob("R??/lineage.json"):
            if CASE_ID_PATTERN.fullmatch(path.parent.name):
                found.add(path.parent.name)
    return sorted(found)


def case_record(
    case_id: str,
    output_root: Path,
    inventory: dict[str, object],
) -> tuple[dict[str, object] | None, dict[str, object] | None, str]:
    """Load one case lineage, preferring the explicit case artifact."""

    lineage_path = output_root / "cases" / case_id / "lineage.json"
    if lineage_path.is_file():
        return load_json(lineage_path), binding(lineage_path), "case_directory"
    cases = inventory.get("cases")
    inline = cases.get(case_id) if isinstance(cases, dict) else None
    if isinstance(inline, dict):
        return inline, None, "inventory_inline"
    return None, None, "missing"


def history_path(lineage: dict[str, object], kind: str) -> Path | None:
    """Return one assembled MHD or user history path."""

    histories = lineage.get("histories")
    record = histories.get(kind) if isinstance(histories, dict) else None
    if not isinstance(record, dict) or record.get("available") is not True:
        return None
    value = record.get("path")
    return Path(value) if isinstance(value, str) else None


def selected_variant_metadata(lineage: dict[str, object]) -> dict[str, list[str]]:
    """Return exact selected-lineage variant metadata."""

    variants = lineage.get("lineage_variants")
    overrides = lineage.get("lineage_command_line_overrides")
    scopes = lineage.get("lineage_claim_scopes")
    return {
        "variants": [str(value) for value in variants]
        if isinstance(variants, list)
        else [],
        "command_line_overrides": [str(value) for value in overrides]
        if isinstance(overrides, list)
        else [],
        "declared_claim_scopes": [str(value) for value in scopes]
        if isinstance(scopes, list)
        else [],
    }


def strict_model_value(lineage: dict[str, object]) -> str | None:
    """Return the normalized effective strict-admissibility model choice."""

    model = lineage.get("model_choices")
    value = model.get("cgl_lf_strict_admissibility") if isinstance(model, dict) else None
    return str(value).strip().lower() if value is not None else None


def comparable_binding(value: object) -> dict[str, object] | None:
    """Return the stable identity fields from one file binding."""

    if not isinstance(value, dict):
        return None
    return {key: value.get(key) for key in ("path", "size_bytes", "sha256")}


def current_declared_binding(
    value: object, label: str
) -> tuple[dict[str, object] | None, list[str]]:
    """Verify a declared regular-file binding against its current bytes."""

    declared = comparable_binding(value)
    if declared is None or not isinstance(declared.get("path"), str):
        return None, [f"{label} binding is missing or malformed"]
    try:
        current = binding(Path(str(declared["path"])))
    except (OSError, FastAcceptanceError) as error:
        return None, [f"{label} binding cannot be verified: {type(error).__name__}: {error}"]
    if comparable_binding(current) != declared:
        return current, [f"{label} declared binding differs from current bytes"]
    return current, []


def authoritative_matrix_case(
    policy: dict[str, object], case_id: str
) -> tuple[dict[str, object] | None, list[str]]:
    """Return the unique reviewed-manifest case record."""

    manifest = policy.get("manifest")
    cases = manifest.get("cases") if isinstance(manifest, dict) else None
    matches = [
        value
        for value in cases
        if isinstance(cases, list)
        and isinstance(value, dict)
        and value.get("id") == case_id
    ] if isinstance(cases, list) else []
    if len(matches) != 1:
        return None, [f"reviewed Stage I manifest lacks unique case {case_id}"]
    return dict(matches[0]), []


def authoritative_input_binding(
    policy: dict[str, object], matrix_case: dict[str, object]
) -> tuple[dict[str, object] | None, dict[str, object] | None, list[str]]:
    """Bind the reviewed matrix and its authoritative case input."""

    verified = policy.get("verified_sources")
    declared_matrix = (
        verified.get("stage_i_manifest") if isinstance(verified, dict) else None
    )
    matrix_binding, errors = current_declared_binding(
        declared_matrix, "reviewed Stage I manifest"
    )
    if matrix_binding is None:
        return None, None, errors
    relative_value = matrix_case.get("input")
    relative = Path(str(relative_value)) if isinstance(relative_value, str) else None
    if relative is None or relative.is_absolute() or ".." in relative.parts:
        return matrix_binding, None, [*errors, "matrix case input path is invalid"]
    matrix_path = Path(str(matrix_binding["path"]))
    matches = {
        candidate.resolve(strict=True)
        for parent in matrix_path.parents
        if (candidate := parent / relative).is_file()
    }
    if len(matches) != 1:
        return matrix_binding, None, [
            *errors,
            f"matrix-authoritative input is not uniquely available: {relative}",
        ]
    try:
        input_binding = binding(matches.pop())
    except (OSError, FastAcceptanceError) as error:
        return matrix_binding, None, [
            *errors,
            f"matrix-authoritative input cannot be bound: {type(error).__name__}: {error}",
        ]
    return matrix_binding, input_binding, errors


def recomputed_model_choices(
    lineage: dict[str, object],
) -> tuple[dict[str, str] | None, dict[str, object] | None, list[str]]:
    """Recompute effective model choices from the currently bound input bytes."""

    current_input, errors = current_declared_binding(lineage.get("input"), "case input")
    if current_input is None:
        return None, None, errors
    try:
        workflow = load_workflow_module()
        source = Path(str(current_input["path"])).read_text(encoding="utf-8")
        overrides = selected_variant_metadata(lineage)["command_line_overrides"]
        choices = {
            str(key): str(value)
            for key, value in workflow.model_choices(source, overrides).items()
        }
        strict = workflow.input_block_value(
            source, "mhd", "cgl_lf_strict_admissibility"
        ) or "false"
        prefix = "mhd/cgl_lf_strict_admissibility="
        for override in overrides:
            if override.startswith(prefix):
                strict = override[len(prefix):]
        choices["cgl_lf_strict_admissibility"] = str(strict)
    except (OSError, UnicodeDecodeError, FastAcceptanceError, ValueError) as error:
        errors.append(
            f"effective model choices cannot be recomputed: {type(error).__name__}: {error}"
        )
        return None, current_input, errors
    return choices, current_input, errors


def configuration_identity_errors(
    policy: dict[str, object], case_id: str, lineage: dict[str, object]
) -> list[str]:
    """Authenticate matrix identity, bound input bytes, and effective model choices."""

    errors: list[str] = []
    matrix_case, matrix_errors = authoritative_matrix_case(policy, case_id)
    errors.extend(matrix_errors)
    if matrix_case is None:
        return errors
    if lineage.get("matrix_case") != matrix_case:
        errors.append("lineage matrix_case differs from reviewed Stage I manifest")
    matrix_binding, matrix_input, binding_errors = authoritative_input_binding(
        policy, matrix_case
    )
    errors.extend(binding_errors)
    choices, current_input, choice_errors = recomputed_model_choices(lineage)
    errors.extend(choice_errors)
    if matrix_input is not None and current_input is not None and (
        current_input.get("sha256") != matrix_input.get("sha256")
        or current_input.get("size_bytes") != matrix_input.get("size_bytes")
    ):
        errors.append("case bound input differs from matrix-authoritative input")

    identities = lineage.get("lineage_identities")
    if case_id != "R02" or identities is not None:
        if not isinstance(identities, dict):
            errors.append("lineage identities are missing or malformed")
        else:
            expected_input_sha = matrix_input.get("sha256") if matrix_input else None
            expected_matrix_sha = matrix_binding.get("sha256") if matrix_binding else None
            if identities.get("input_sha256") != [expected_input_sha]:
                errors.append("lineage input identity differs from matrix-authoritative input")
            if identities.get("matrix_sha256") != [expected_matrix_sha]:
                errors.append("lineage matrix identity differs from reviewed Stage I manifest")

    declared_model = lineage.get("model_choices")
    if choices is None:
        return errors
    if not isinstance(declared_model, dict):
        errors.append("lineage model choices are missing or malformed")
        return errors
    declared = {str(key): str(value) for key, value in declared_model.items()}
    if case_id == "R02" and "cgl_lf_strict_admissibility" not in declared:
        expected = dict(choices)
        expected.pop("cgl_lf_strict_admissibility")
    else:
        expected = choices
    if declared != expected:
        errors.append("lineage model choices differ from bound-input effective choices")
    return errors


def lineage_identity_errors(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    lineage: dict[str, object],
    lineage_binding: dict[str, object] | None,
) -> list[str]:
    """Authenticate direct-fast case identity and effective configuration."""

    errors: list[str] = []
    expected_name = acceptance.case_name(policy, case_id)
    if lineage.get("case_id") != case_id:
        errors.append(
            f"lineage case_id differs: {lineage.get('case_id')!r} != {case_id!r}"
        )
    if lineage.get("case_name") != expected_name:
        errors.append(
            "lineage case_name differs: "
            f"{lineage.get('case_name')!r} != {expected_name!r}"
        )
    if lineage_binding is None:
        errors.append("case lineage is not bound to an explicit lineage.json artifact")
    retained = lineage.get("errors")
    if not isinstance(retained, list):
        errors.append("lineage errors field is malformed")
    elif retained:
        errors.extend(f"retained lineage error: {value}" for value in retained)
    errors.extend(configuration_identity_errors(policy, case_id, lineage))

    selected = selected_variant_metadata(lineage)
    strict = strict_model_value(lineage)
    if case_id == "R14":
        if selected["variants"] != [R14_NONFATAL_VARIANT]:
            errors.append("R14 selected variant list is not the exact admitted variant")
        if selected["command_line_overrides"] != [R14_OVERRIDE]:
            errors.append("R14 override list is not exactly the admitted false override")
        if strict != "false":
            errors.append("R14 effective strict-admissibility model choice is not false")
    else:
        if selected["variants"]:
            errors.append(f"standard case has unexpected variants: {selected['variants']}")
        if selected["command_line_overrides"]:
            errors.append(
                "standard case has unexpected command-line overrides: "
                f"{selected['command_line_overrides']}"
            )
        # Only the accepted legacy R02 bundle may omit this explicit model-choice field.
        if strict != "true" and not (case_id == "R02" and strict is None):
            errors.append("standard case effective strict-admissibility choice is not true")
    return errors


def scientific_scope(case_id: str, lineage: dict[str, object]) -> dict[str, object]:
    """Return the explicit direct-fast claim scope for one case."""

    selected = selected_variant_metadata(lineage)
    if case_id == "R10":
        return {
            "classification": "exploratory_only",
            "campaign_interpretation_eligible": False,
            "uniform_strict_diagnostics_eligible": False,
            "hard_bound_is_fatal": True,
            "reason": (
                "R10 is retained only as an exploratory compressive stress test because "
                "the qualified executable regularizes negative CGL fast-speed "
                "discriminants with an absolute value."
            ),
            "allowed_claims": [
                "exploratory behavior of the qualified regularized implementation"
            ],
            "excluded_claims": [
                "strict-hyperbolic CGL evidence",
                "claim-grade beta trend",
                "uniform campaign acceptance",
            ],
            **selected,
        }
    if case_id == "R14":
        admitted = (
            selected["variants"] == [R14_NONFATAL_VARIANT]
            and selected["command_line_overrides"] == [R14_OVERRIDE]
            and strict_model_value(lineage) == "false"
        )
        return {
            "classification": (
                "scoped_nonfatal_hard_bound_variant"
                if admitted
                else "r14_variant_not_authenticated"
            ),
            "campaign_interpretation_eligible": admitted,
            "uniform_strict_diagnostics_eligible": False,
            "hard_bound_is_fatal": not admitted,
            "reason": (
                "R14 is admitted for the finite-rate limiter comparison with every "
                "hard-bound event retained as a nonfatal science diagnostic."
                if admitted
                else (
                    "R14 lacks the exact recorded nonfatal-hard-bound variant "
                    "and override."
                )
            ),
            "allowed_claims": (
                [
                    "finite-rate limiter response with hard-bound occupancy disclosed",
                    "numerical health excluding hard-bound occupancy as a fatal counter",
                ]
                if admitted
                else []
            ),
            "excluded_claims": [
                "uniform strict-admissibility compliance",
                "hard-bound-free evolution",
                "pure pressure-limiter comparison independent of LF transport",
            ],
            **selected,
        }
    return {
        "classification": "standard_claim_scope",
        "campaign_interpretation_eligible": True,
        "uniform_strict_diagnostics_eligible": True,
        "hard_bound_is_fatal": True,
        "reason": "standard direct-fast Stage I claim scope",
        "allowed_claims": ["declared Stage I comparison after all required audits pass"],
        "excluded_claims": [],
        **selected,
    }


def load_histories(
    acceptance: object,
    lineage: dict[str, object],
) -> tuple[
    dict[str, dict[str, list[float]]],
    dict[str, dict[str, object]],
    list[str],
]:
    """Load and bind available assembled histories through reviewed parsing."""

    histories: dict[str, dict[str, list[float]]] = {}
    bindings: dict[str, dict[str, object]] = {}
    errors: list[str] = []
    for kind in ("mhd", "user"):
        path = history_path(lineage, kind)
        if path is None:
            errors.append(f"evidence gap: assembled {kind} history is unavailable")
            continue
        try:
            data, reviewed_binding = acceptance.load_history(
                path, f"direct-fast {kind} history"
            )
            histories_value = lineage.get("histories")
            declared = (
                histories_value.get(kind)
                if isinstance(histories_value, dict)
                else None
            )
            declared_binding = (
                declared.get("binding") if isinstance(declared, dict) else None
            )
            if not isinstance(declared_binding, dict):
                errors.append(f"{kind} history lacks a declared binding")
                continue
            comparable = {
                key: reviewed_binding.get(key)
                for key in ("path", "size_bytes", "sha256")
            }
            declared_comparable = {
                key: declared_binding.get(key)
                for key in ("path", "size_bytes", "sha256")
            }
            if comparable != declared_comparable:
                errors.append(f"{kind} history declared binding differs from current bytes")
                continue
            histories[kind] = data
            bindings[kind] = dict(reviewed_binding)
        except (acceptance.AcceptanceError, OSError) as error:
            errors.append(f"{kind} history: {type(error).__name__}: {error}")
    return histories, bindings, errors


def maximum_counter(history: dict[str, list[float]], column: str) -> float | None:
    """Return the maximum finite value of one cumulative counter."""

    values = history.get(column)
    return max(values) if values else None


def numerical_health(
    case_id: str,
    lineage: dict[str, object],
    scope: dict[str, object],
    histories: dict[str, dict[str, list[float]]],
    history_errors: list[str],
    identity_errors: list[str],
) -> dict[str, object]:
    """Assess lean direct-fast completion and numerical health."""

    target = float(lineage.get("target_time", 10.0))
    status = str(lineage.get("status", "unknown"))
    final_times = {
        kind: data["time"][-1] for kind, data in sorted(histories.items())
    }
    synchronized = (
        len(final_times) == 2
        and math.isclose(
            final_times["mhd"],
            final_times["user"],
            rel_tol=0.0,
            abs_tol=TIME_TOLERANCE,
        )
    )
    observed_final = min(final_times.values()) if final_times else None
    complete = (
        status == "complete"
        and synchronized
        and observed_final is not None
        and observed_final >= target - TIME_TOLERANCE
    )
    mhd = histories.get("mhd", {})
    counters = {
        column: maximum_counter(mhd, column)
        for column in (*FATAL_COUNTERS, HARD_BOUND_COUNTER)
    }
    fatal_nonzero = {
        column: value
        for column, value in counters.items()
        if column in FATAL_COUNTERS and value is not None and value != 0.0
    }
    hard_bound = counters[HARD_BOUND_COUNTER]
    hard_bound_fatal = (
        hard_bound is not None
        and hard_bound != 0.0
        and scope.get("hard_bound_is_fatal") is True
    )
    evidence_gaps = [
        value for value in history_errors if value.startswith("evidence gap:")
    ]
    structural_errors = [
        *identity_errors,
        *(value for value in history_errors if not value.startswith("evidence gap:")),
    ]
    if status == "complete":
        structural_errors.extend(evidence_gaps)
        evidence_gaps = []
    if len(histories) == 2 and not synchronized:
        structural_errors.append(
            "assembled MHD and user histories have different final times"
        )
    missing_counters = sorted(
        column for column, value in counters.items() if value is None
    )
    if structural_errors or fatal_nonzero or hard_bound_fatal:
        result = "fail"
    elif complete and not missing_counters:
        result = "pass"
    else:
        result = "inconclusive"
    hard_bound_disposition = (
        "retained_nonfatal_science_diagnostic"
        if hard_bound is not None
        and hard_bound != 0.0
        and scope.get("hard_bound_is_fatal") is False
        else "fatal_counter_nonzero"
        if hard_bound_fatal
        else "zero_or_unavailable"
    )
    return {
        "result": result,
        "assembly_status": status,
        "complete_to_target": complete,
        "target_time": target,
        "observed_final_time": observed_final,
        "history_final_times": final_times,
        "history_times_synchronized": synchronized,
        "structural_errors": structural_errors,
        "evidence_gaps": evidence_gaps,
        "fatal_counter_maxima": {
            column: counters[column] for column in FATAL_COUNTERS
        },
        "missing_counter_columns": missing_counters,
        "fatal_nonzero_counters": fatal_nonzero,
        "hard_bound_maximum": hard_bound,
        "hard_bound_disposition": hard_bound_disposition,
    }


def forcing_tcorr(lineage: dict[str, object]) -> float:
    """Return forcing correlation time recomputed from the bound input."""

    model, _input_binding, errors = recomputed_model_choices(lineage)
    value = model.get("forcing_tcorr") if isinstance(model, dict) and not errors else None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return 0.0
    return result if math.isfinite(result) and result >= 0.0 else 0.0


def window_is_covered(times: list[float], start: float, end: float) -> bool:
    """Return whether one history covers exact analysis endpoints."""

    return times[0] <= start and times[-1] >= end


def available_history_statistics(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    lineage: dict[str, object],
    histories: dict[str, dict[str, list[float]]],
) -> tuple[dict[str, object], list[str]]:
    """Compute every exact reviewed window currently covered by a case."""

    criteria = policy["criteria"]
    windows = criteria["analysis_windows"]
    statistics_policy = criteria["statistics_policy"]
    gap_policy = statistics_policy["gap_policy"]
    metrics: dict[str, object] = {}
    errors: list[str] = []
    for metric, spec in sorted(criteria["case_metrics"].items()):
        source = str(spec["history"])
        column = str(spec["column"])
        history = histories.get(source)
        if history is None or column not in history:
            errors.append(f"{metric}: source {source}/{column} is unavailable")
            continue
        metric_windows: dict[str, object] = {}
        for window_name, limits in sorted(windows.items()):
            start, end = (float(value) for value in limits)
            if not window_is_covered(history["time"], start, end):
                metric_windows[window_name] = {
                    "result": "inconclusive",
                    "reason": "history does not yet cover the exact reviewed window",
                    "window": {"start": start, "end": end},
                }
                continue
            try:
                metric_windows[window_name] = {
                    "result": "available",
                    "statistics": acceptance.window_statistics(
                        history["time"],
                        history[column],
                        start,
                        end,
                        replicates=int(statistics_policy["bootstrap_replicates"]),
                        seed_text=(
                            f"direct-fast:{case_id}:{metric}:{window_name}:"
                            f"{policy['criteria_binding']['sha256']}"
                        ),
                        minimum_block_duration=forcing_tcorr(lineage),
                        expected_cadence=float(gap_policy["expected_history_cadence"]),
                        maximum_gap_expected_cadence_multiplier=float(
                            gap_policy["maximum_gap_expected_cadence_multiplier"]
                        ),
                        maximum_gap_minimum_block_duration_fraction=float(
                            gap_policy["maximum_gap_forcing_tcorr_fraction"]
                        ),
                    ),
                }
            except acceptance.AcceptanceError as error:
                metric_windows[window_name] = {
                    "result": "inconclusive",
                    "reason": f"{type(error).__name__}: {error}",
                    "window": {"start": start, "end": end},
                }
        metrics[str(metric)] = {
            "history": source,
            "column": column,
            "windows": metric_windows,
        }
    return metrics, errors


def reviewed_complete_case_evidence(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    lineage: dict[str, object],
    health: dict[str, object],
) -> tuple[dict[str, object] | None, str | None]:
    """Evaluate one target-complete case only with an authenticated canonical bundle."""

    if health.get("result") != "pass":
        return None, "case does not pass direct-fast numerical health"
    mhd_path = history_path(lineage, "mhd")
    user_path = history_path(lineage, "user")
    if mhd_path is None or user_path is None:
        return None, "assembled MHD or user history is unavailable"
    canonical = acceptance.canonical_bundle_path(case_id)
    retained_segments = lineage.get("lineage")
    retained_manifests: list[Path] = []
    if isinstance(retained_segments, list):
        for segment in retained_segments:
            manifest = segment.get("manifest") if isinstance(segment, dict) else None
            if (
                isinstance(segment, dict)
                and segment.get("kind") == "accepted_r02_bundle"
                and isinstance(manifest, dict)
                and isinstance(manifest.get("path"), str)
            ):
                retained_manifests.append(
                    Path(str(manifest["path"])).resolve(strict=True)
                )
    if canonical.resolve(strict=False) not in retained_manifests:
        return None, (
            "no authenticated canonical accepted bundle is retained; "
            "direct-fast exact-window evidence remains non-authorizing"
        )
    try:
        evidence = acceptance.evaluate_case(
            policy,
            case_id,
            mhd_path,
            user_path,
            None,
            None,
            canonical,
        )
    except (acceptance.AcceptanceError, OSError) as error:
        return None, f"{type(error).__name__}: {error}"
    return evidence, None


def comparison_case_evidence(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    lineage: dict[str, object],
    scope: dict[str, object],
    health: dict[str, object],
    histories: dict[str, dict[str, list[float]]],
    reviewed: dict[str, object] | None = None,
) -> tuple[dict[str, object] | None, str | None]:
    """Build scoped comparison evidence with the authenticated forcing time."""

    if health.get("result") != "pass":
        return None, "case does not pass direct-fast numerical health"
    if scope.get("campaign_interpretation_eligible") is not True:
        return None, "case scope is not eligible for campaign interpretation"
    tcorr = forcing_tcorr(lineage)
    if tcorr <= 0.0:
        return None, "effective forcing correlation time is unavailable"
    if reviewed is not None and isinstance(reviewed.get("metrics"), dict):
        metrics = reviewed["metrics"]
        adequate = all(
            isinstance(record, dict)
            and record.get("sampling_adequacy") == "pass"
            and isinstance(record.get("stationarity"), dict)
            and record["stationarity"].get("result") == "pass"
            for record in metrics.values()
        )
        return {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-comparison-evidence",
            "authority": "non-authorizing-direct-fast-scientific-assessment",
            "case_id": case_id,
            "result": "pass" if adequate else "inconclusive",
            "minimum_block_duration": tcorr,
            "metrics": metrics,
            "analyzer_metrics": reviewed.get("analyzer_metrics", {}),
        }, None
    metrics: dict[str, object] = {}
    criteria = policy["criteria"]
    try:
        for metric, spec in sorted(criteria["case_metrics"].items()):
            source = str(spec["history"])
            column = str(spec["column"])
            history = histories.get(source)
            if history is None or column not in history:
                return None, f"{metric}: required history column is unavailable"
            metrics[str(metric)] = acceptance.metric_statistics(
                history,
                history[column],
                str(metric),
                policy,
                kind=str(spec["stationarity_kind"]),
                minimum_block_duration=tcorr,
            )
    except acceptance.AcceptanceError as error:
        return None, f"{type(error).__name__}: {error}"
    adequate = all(
        isinstance(record, dict)
        and record.get("sampling_adequacy") == "pass"
        and isinstance(record.get("stationarity"), dict)
        and record["stationarity"].get("result") == "pass"
        for record in metrics.values()
    )
    return {
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-direct-fast-comparison-evidence",
        "authority": "non-authorizing-direct-fast-scientific-assessment",
        "case_id": case_id,
        "result": "pass" if adequate else "inconclusive",
        "minimum_block_duration": tcorr,
        "metrics": metrics,
        "analyzer_metrics": {},
    }, None


def history_statistics_from_reviewed(
    policy: dict[str, object],
    evidence: dict[str, object],
) -> dict[str, object]:
    """Normalize reviewed complete-case metrics to the partial-case table schema."""

    reviewed_metrics = evidence.get("metrics")
    if not isinstance(reviewed_metrics, dict):
        return {}
    result: dict[str, object] = {}
    for metric, spec in sorted(policy["criteria"]["case_metrics"].items()):
        record = reviewed_metrics.get(metric)
        if not isinstance(record, dict):
            continue
        windows: dict[str, object] = {}
        for window in sorted(policy["criteria"]["analysis_windows"]):
            statistics = record.get(window)
            if isinstance(statistics, dict):
                windows[window] = {
                    "result": "available",
                    "statistics": statistics,
                }
        result[str(metric)] = {
            "history": spec["history"],
            "column": spec["column"],
            "windows": windows,
            "sampling_adequacy": record.get("sampling_adequacy"),
            "stationarity": record.get("stationarity"),
        }
    return result


def direct_fast_case_result(
    case_id: str,
    health: dict[str, object],
    scope: dict[str, object],
    reviewed: dict[str, object] | None,
) -> tuple[str, str]:
    """Return the scoped direct-fast case disposition."""

    if health["result"] == "fail":
        return "fail", "direct-fast structural or numerical health failed"
    if case_id == "R10":
        return "inconclusive", "R10 is deliberately exploratory-only"
    if case_id == "R14" and scope["campaign_interpretation_eligible"] is not True:
        return "fail", "R14 lacks its exact admitted nonfatal-hard-bound variant"
    if health["result"] != "pass":
        return "inconclusive", "case is incomplete or required histories are unavailable"
    if reviewed is None:
        return "inconclusive", "reviewed complete-case evaluation is unavailable"
    if reviewed.get("result") == "fail":
        return "fail", "reviewed per-case scientific evaluation failed"
    if reviewed.get("result") == "pass":
        return "pass", "direct-fast health and reviewed per-case evaluation passed"
    return "inconclusive", "reviewed per-case evaluation retains inconclusive gates"


def summarize_case(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    lineage: dict[str, object] | None,
    lineage_binding: dict[str, object] | None,
    source: str,
    case_output: Path,
) -> tuple[dict[str, object], dict[str, object] | None]:
    """Build and write one direct-fast case assessment."""

    reviewed_path = case_output / "reviewed_case_evidence.json"
    reviewed_path.unlink(missing_ok=True)
    if lineage is None:
        summary = {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-case-acceptance",
            "authority": "non-authorizing-direct-fast-scientific-assessment",
            "case_id": case_id,
            "case_name": None,
            "result": "inconclusive",
            "reason": "case is absent from the fast-report inventory and case directory",
            "source": source,
            "scope": {
                "classification": (
                    "exploratory_only" if case_id == "R10" else "unknown_missing_case"
                ),
                "campaign_interpretation_eligible": False,
            },
            "health": {
                "result": "inconclusive",
                "complete_to_target": False,
                "structural_errors": ["case fast-report record is unavailable"],
            },
            "history_statistics": {},
            "reviewed_case_evidence": None,
            "provenance": {"lineage": None, "histories": {}},
        }
        write_json(case_output / "case_acceptance.json", summary)
        return summary, None

    identity_errors = lineage_identity_errors(
        acceptance, policy, case_id, lineage, lineage_binding
    )
    scope = scientific_scope(case_id, lineage)
    histories, history_bindings, history_errors = load_histories(acceptance, lineage)
    health = numerical_health(
        case_id, lineage, scope, histories, history_errors, identity_errors
    )
    reviewed, reviewed_error = reviewed_complete_case_evidence(
        acceptance, policy, case_id, lineage, health
    )
    comparison, comparison_error = comparison_case_evidence(
        acceptance, policy, case_id, lineage, scope, health, histories, reviewed
    )
    if comparison is not None:
        history_statistics = history_statistics_from_reviewed(policy, comparison)
        statistics_errors: list[str] = []
    else:
        history_statistics, statistics_errors = available_history_statistics(
            acceptance, policy, case_id, lineage, histories
        )
    reviewed_binding: dict[str, object] | None = None
    if reviewed is not None:
        write_json(reviewed_path, reviewed)
        reviewed_binding = binding(reviewed_path)
    result, reason = direct_fast_case_result(case_id, health, scope, reviewed)
    summary = {
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-direct-fast-case-acceptance",
        "authority": "non-authorizing-direct-fast-scientific-assessment",
        "case_id": case_id,
        "case_name": lineage.get("case_name"),
        "result": result,
        "reason": reason,
        "source": source,
        "scope": scope,
        "health": health,
        "history_statistics": history_statistics,
        "history_statistics_errors": statistics_errors,
        "reviewed_case_evidence": (
            {
                "result": reviewed.get("result"),
                "binding": reviewed_binding,
            }
            if reviewed is not None
            else {
                "result": "inconclusive",
                "reason": reviewed_error,
                "binding": None,
            }
        ),
        "comparison_evidence": (
            {
                "result": comparison.get("result"),
                "minimum_block_duration": comparison.get("minimum_block_duration"),
            }
            if comparison is not None
            else {
                "result": "inconclusive",
                "reason": comparison_error,
            }
        ),
        "provenance": {
            "lineage": lineage_binding,
            "histories": history_bindings,
        },
    }
    write_json(case_output / "case_acceptance.json", summary)
    return summary, comparison


def campaign_gate(
    acceptance: object,
    name: str,
    result: str,
    reason: str,
    observations: object = None,
) -> dict[str, object]:
    """Create one campaign gate through the reviewed gate constructor."""

    return acceptance.gate(
        name, result, reason=reason, observations=observations
    )


def reviewed_pair_gates(
    acceptance: object,
    policy: dict[str, object],
    comparison_cases: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    """Build reviewed-kernel active/passive gates from scoped comparison evidence."""

    result: list[dict[str, object]] = []
    pairs = policy["criteria"]["family_gates"]["active_passive"]["pairs"]
    for active, passive in pairs:
        active_evidence = comparison_cases.get(active)
        passive_evidence = comparison_cases.get(passive)
        if (
            isinstance(active_evidence, dict)
            and active_evidence.get("result") == "pass"
            and isinstance(passive_evidence, dict)
            and passive_evidence.get("result") == "pass"
        ):
            contrast = acceptance.pair_contrast(
                active_evidence, passive_evidence, policy
            )
        else:
            contrast = {
                "result": "inconclusive",
                "reason": (
                    "one or both case comparison-evidence records are unavailable "
                    "or do not pass"
                ),
                "metrics": None,
            }
        result.append(campaign_gate(
            acceptance,
            f"active_passive_pair:{active}:{passive}",
            str(contrast["result"]),
            str(contrast["reason"]),
            contrast.get("metrics"),
        ))
    return result


def finite_limiter_gate(
    acceptance: object,
    comparison_cases: dict[str, dict[str, object]],
    case_summaries: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Build the scoped R15 greater-than R14 finite-limiter ordering gate."""

    scope = case_summaries.get("R14", {}).get("scope")
    admitted = (
        isinstance(scope, dict)
        and scope.get("classification") == "scoped_nonfatal_hard_bound_variant"
    )
    if not admitted:
        return campaign_gate(
            acceptance,
            "finite_limiter_ordering:R15_gt_R14",
            "inconclusive",
            "R14 exact nonfatal-hard-bound variant is unavailable",
        )
    r14_evidence = comparison_cases.get("R14")
    r15_evidence = comparison_cases.get("R15")
    if (
        not isinstance(r14_evidence, dict)
        or r14_evidence.get("result") != "pass"
        or not isinstance(r15_evidence, dict)
        or r15_evidence.get("result") != "pass"
    ):
        return campaign_gate(
            acceptance,
            "finite_limiter_ordering:R15_gt_R14",
            "inconclusive",
            "R14 or R15 comparison evidence is unavailable or does not pass",
        )
    lower = acceptance.scalar_from_case(r14_evidence, "nu_eff", "late")
    upper = acceptance.scalar_from_case(r15_evidence, "nu_eff", "late")
    if lower is None or upper is None:
        return campaign_gate(
            acceptance,
            "finite_limiter_ordering:R15_gt_R14",
            "inconclusive",
            "late-time R14 or R15 effective-collisionality statistics are unavailable",
        )
    difference = upper["mean"] - lower["mean"]
    lower_95 = difference - 1.96 * math.hypot(
        upper["standard_error"], lower["standard_error"]
    )
    return campaign_gate(
        acceptance,
        "finite_limiter_ordering:R15_gt_R14",
        "pass" if lower_95 > 0.0 else "fail",
        (
            "R15 effective collisionality exceeds R14 with 95% confidence"
            if lower_95 > 0.0
            else "finite-limiter ordering is not resolved"
        ),
        {
            "R15_minus_R14": difference,
            "lower_95": lower_95,
            "R14_scope": scope,
        },
    )


def build_campaign_evidence(
    acceptance: object,
    policy: dict[str, object],
    case_summaries: dict[str, dict[str, object]],
    comparison_cases: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Build lean scoped campaign evidence with reviewed comparison kernels."""

    required = [str(value) for value in policy["criteria"]["required_cases"]]
    missing = sorted(set(required) - set(case_summaries))
    incomplete = sorted(
        case_id
        for case_id in required
        if case_id not in case_summaries
        or case_summaries[case_id].get("health", {}).get("complete_to_target") is not True
    )
    numerical_failures = sorted(
        case_id
        for case_id, value in case_summaries.items()
        if value.get("health", {}).get("result") == "fail"
    )
    claim_grade_cases = sorted(
        case_id
        for case_id, value in case_summaries.items()
        if value.get("scope", {}).get("campaign_interpretation_eligible") is True
    )
    exploratory = sorted(
        case_id
        for case_id, value in case_summaries.items()
        if value.get("scope", {}).get("classification") == "exploratory_only"
    )
    scoped_variants = sorted(
        case_id
        for case_id, value in case_summaries.items()
        if value.get("scope", {}).get("classification")
        == "scoped_nonfatal_hard_bound_variant"
    )
    case_results = {
        case_id: str(value.get("result"))
        for case_id, value in sorted(case_summaries.items())
    }
    gates = [
        campaign_gate(
            acceptance,
            "required_case_inventory",
            "pass" if not missing else "inconclusive",
            (
                "all required cases are inventoried"
                if not missing
                else "required cases are missing"
            ),
            {"missing_cases": missing},
        ),
        campaign_gate(
            acceptance,
            "completion_to_target",
            "pass" if not missing and not incomplete else "inconclusive",
            (
                "all required cases reach the target time"
                if not missing and not incomplete
                else "one or more required cases remain partial or unavailable"
            ),
            {"incomplete_cases": incomplete},
        ),
        campaign_gate(
            acceptance,
            "direct_fast_numerical_health",
            "fail" if numerical_failures else (
                "pass" if not missing and not incomplete else "inconclusive"
            ),
            (
                "one or more cases failed direct-fast numerical health"
                if numerical_failures
                else "all completed cases pass direct-fast numerical health"
            ),
            {"failed_cases": numerical_failures},
        ),
        campaign_gate(
            acceptance,
            "scoped_exact_window_comparison_evidence",
            (
                "pass"
                if len(comparison_cases) == len(required) - 1
                and all(
                    value.get("result") == "pass"
                    for value in comparison_cases.values()
                )
                else "fail"
                if any(value.get("result") == "fail" for value in comparison_cases.values())
                else "inconclusive"
            ),
            "tcorr-aware exact-window comparison evidence is retained for every "
            "campaign-eligible target-complete case; R10 is excluded",
            {
                "available_cases": sorted(comparison_cases),
                "unavailable_cases": sorted(
                    set(required) - {"R10"} - set(comparison_cases)
                ),
                "results": {
                    case_id: value.get("result")
                    for case_id, value in sorted(comparison_cases.items())
                },
            },
        ),
        campaign_gate(
            acceptance,
            "explicit_claim_scope",
            "pass",
            "R10 exploratory and R14 nonfatal-hard-bound scopes are explicit",
            {
                "claim_grade_cases": claim_grade_cases,
                "exploratory_cases": exploratory,
                "scoped_variants": scoped_variants,
            },
        ),
        *reviewed_pair_gates(acceptance, policy, comparison_cases),
        finite_limiter_gate(acceptance, comparison_cases, case_summaries),
    ]
    result = acceptance.aggregate_gate_result(gates)
    return acceptance.seal_evidence({
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-direct-fast-campaign-evidence",
        "authority": "non-authorizing-direct-fast-scientific-assessment",
        "release_authorizing": False,
        "result": result,
        "case_results": case_results,
        "claim_scope": {
            "claim_grade_cases": claim_grade_cases,
            "exploratory_cases": exploratory,
            "scoped_variants": scoped_variants,
            "R10": case_summaries.get("R10", {}).get("scope"),
            "R14": case_summaries.get("R14", {}).get("scope"),
        },
        "gates": gates,
        "provenance": {
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "reviewed_acceptance_utility": binding(ACCEPTANCE_UTILITY),
        },
    })


def table_value(value: object) -> str:
    """Format one scalar for publication tables."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.10g}"
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


def markdown_escape(value: object) -> str:
    """Format and escape one Markdown table cell."""

    return table_value(value).replace("|", "\\|").replace("\n", " ")


def write_markdown_table(
    path: Path,
    title: str,
    fields: list[str],
    rows: Iterable[dict[str, object]],
) -> None:
    """Write one deterministic publication-ready Markdown table."""

    row_values = list(rows)
    lines = [
        f"# {title}",
        "",
        "| " + " | ".join(fields) + " |",
        "| " + " | ".join("---" for _ in fields) + " |",
    ]
    for row in row_values:
        lines.append(
            "| " + " | ".join(markdown_escape(row.get(field)) for field in fields) + " |"
        )
    lines.append("")
    write_text(path, "\n".join(lines))


def case_table_rows(
    required: list[str],
    cases: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    """Return one publication case-disposition row per required case."""

    rows: list[dict[str, object]] = []
    for case_id in required:
        value = cases.get(case_id, {})
        health = value.get("health") if isinstance(value.get("health"), dict) else {}
        scope = value.get("scope") if isinstance(value.get("scope"), dict) else {}
        reviewed = (
            value.get("reviewed_case_evidence")
            if isinstance(value.get("reviewed_case_evidence"), dict)
            else {}
        )
        counters = (
            health.get("fatal_counter_maxima")
            if isinstance(health.get("fatal_counter_maxima"), dict)
            else {}
        )
        rows.append({
            "case_id": case_id,
            "case_name": value.get("case_name"),
            "scope": scope.get("classification"),
            "assembly_status": health.get("assembly_status"),
            "complete_to_target": health.get("complete_to_target"),
            "final_time": health.get("observed_final_time"),
            "target_time": health.get("target_time"),
            "direct_fast_health": health.get("result"),
            "reviewed_result": reviewed.get("result"),
            "case_result": value.get("result"),
            "dfloor_max": counters.get("lf_dfloor"),
            "pfloor_max": counters.get("lf_pfloor"),
            "nonfinite_max": counters.get("lf_nonfin"),
            "nonpositive_max": counters.get("lf_nonpos"),
            "hard_bound_max": health.get("hard_bound_maximum"),
            "hard_bound_disposition": health.get("hard_bound_disposition"),
            "reason": value.get("reason"),
        })
    return rows


def statistics_table_rows(
    cases: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    """Flatten available reviewed-window statistics."""

    rows: list[dict[str, object]] = []
    for case_id, case in sorted(cases.items()):
        metrics = case.get("history_statistics")
        if not isinstance(metrics, dict):
            continue
        for metric, metric_record in sorted(metrics.items()):
            if not isinstance(metric_record, dict):
                continue
            windows = metric_record.get("windows")
            if not isinstance(windows, dict):
                continue
            for window, window_record in sorted(windows.items()):
                if not isinstance(window_record, dict):
                    continue
                stats = window_record.get("statistics")
                if not isinstance(stats, dict):
                    rows.append({
                        "case_id": case_id,
                        "metric": metric,
                        "window": window,
                        "availability": window_record.get("result"),
                        "reason": window_record.get("reason"),
                    })
                    continue
                confidence = stats.get("confidence_interval_95")
                lower = (
                    confidence[0]
                    if isinstance(confidence, list) and confidence
                    else None
                )
                upper = (
                    confidence[1]
                    if isinstance(confidence, list) and len(confidence) > 1
                    else None
                )
                rows.append({
                    "case_id": case_id,
                    "metric": metric,
                    "window": window,
                    "availability": "available",
                    "mean": stats.get("mean"),
                    "standard_deviation": stats.get("standard_deviation"),
                    "standard_error": stats.get("standard_error"),
                    "ci95_lower": lower,
                    "ci95_upper": upper,
                    "effective_sample_count": stats.get("effective_sample_count"),
                    "independent_time_block_count": stats.get(
                        "independent_time_block_count"
                    ),
                    "gap_adequacy": stats.get("gap_adequacy"),
                    "sample_count": stats.get("sample_count"),
                    "reason": "",
                })
    return rows


def gate_table_rows(campaign: dict[str, object]) -> list[dict[str, object]]:
    """Flatten campaign gates for publication."""

    rows: list[dict[str, object]] = []
    gates = campaign.get("gates")
    if not isinstance(gates, list):
        return rows
    for gate in gates:
        if not isinstance(gate, dict):
            continue
        rows.append({
            "gate": gate.get("name"),
            "result": gate.get("result"),
            "reason": gate.get("reason"),
            "observations": gate.get("observations"),
        })
    return rows


def write_tables(
    output: Path,
    required: list[str],
    cases: dict[str, dict[str, object]],
    campaign: dict[str, object],
) -> list[Path]:
    """Write publication-ready CSV and Markdown tables."""

    tables = output / "tables"
    case_fields = [
        "case_id",
        "case_name",
        "scope",
        "assembly_status",
        "complete_to_target",
        "final_time",
        "target_time",
        "direct_fast_health",
        "reviewed_result",
        "case_result",
        "dfloor_max",
        "pfloor_max",
        "nonfinite_max",
        "nonpositive_max",
        "hard_bound_max",
        "hard_bound_disposition",
        "reason",
    ]
    statistics_fields = [
        "case_id",
        "metric",
        "window",
        "availability",
        "mean",
        "standard_deviation",
        "standard_error",
        "ci95_lower",
        "ci95_upper",
        "effective_sample_count",
        "independent_time_block_count",
        "gap_adequacy",
        "sample_count",
        "reason",
    ]
    gate_fields = ["gate", "result", "reason", "observations"]
    case_rows = case_table_rows(required, cases)
    statistic_rows = statistics_table_rows(cases)
    gate_rows = gate_table_rows(campaign)
    paths = [
        tables / "case_acceptance.csv",
        tables / "case_acceptance.md",
        tables / "history_statistics.csv",
        tables / "history_statistics.md",
        tables / "campaign_gates.csv",
        tables / "campaign_gates.md",
    ]
    write_csv(paths[0], case_fields, case_rows)
    write_markdown_table(paths[1], "Direct-Fast Case Acceptance", case_fields, case_rows)
    write_csv(paths[2], statistics_fields, statistic_rows)
    write_markdown_table(
        paths[3], "Reviewed-Window History Statistics", statistics_fields, statistic_rows
    )
    write_csv(paths[4], gate_fields, gate_rows)
    write_markdown_table(paths[5], "Direct-Fast Campaign Gates", gate_fields, gate_rows)
    return paths


def output_default(inventory_path: Path) -> Path:
    """Return a separate sibling output directory for one inventory."""

    parent = inventory_path.expanduser().absolute().parent
    return parent.with_name(f"{parent.name}-acceptance")


def remove_stale_adapter_outputs(output: Path) -> None:
    """Remove only known generated records before a deterministic rerun."""

    if not output.exists():
        return
    for path in output.glob("cases/R??/reviewed_case_evidence.json"):
        path.unlink(missing_ok=True)
    for path in output.glob("cases/R??/case_acceptance.json"):
        path.unlink(missing_ok=True)
    for name in ("campaign_evidence.json", "summary.json", "provenance.json"):
        (output / name).unlink(missing_ok=True)
    for path in (output / "tables").glob("*"):
        if path.is_file() and path.suffix in {".csv", ".md"}:
            path.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument(
        "--output",
        type=Path,
        help="separate acceptance output; defaults to INVENTORY_PARENT-acceptance",
    )
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--criteria-review", type=Path, default=DEFAULT_CRITERIA_REVIEW)
    parser.add_argument(
        "--cases",
        action="append",
        help="optional repeated comma-separated R02-R17 selection",
    )
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="do not discover additional case lineage files under the inventory output",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the direct-fast scientific-acceptance adapter."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        inventory_path = args.inventory.expanduser().absolute().resolve(strict=True)
        inventory = load_json(inventory_path)
        output_root = inventory_output_root(inventory_path, inventory)
        output = (
            args.output.expanduser().absolute()
            if args.output is not None
            else output_default(inventory_path)
        ).resolve(strict=False)
        if path_contains(output_root, output) or path_contains(output, output_root):
            raise FastAcceptanceError(
                "acceptance output must be separate from the fast-report output"
            )
        reject_nested_output_symlink_escapes(output)
        remove_stale_adapter_outputs(output)
        acceptance = load_acceptance_module()
        policy = acceptance.load_validated_policy(args.criteria, args.criteria_review)
        requested = parse_case_selection(args.cases)
        selected = discover_case_ids(
            output_root, inventory, requested, args.inventory_only
        )
        if not selected:
            raise FastAcceptanceError("no direct-fast cases were selected or discovered")

        case_summaries: dict[str, dict[str, object]] = {}
        comparison_cases: dict[str, dict[str, object]] = {}
        for case_id in selected:
            lineage, lineage_binding, source = case_record(
                case_id, output_root, inventory
            )
            summary, comparison = summarize_case(
                acceptance,
                policy,
                case_id,
                lineage,
                lineage_binding,
                source,
                output / "cases" / case_id,
            )
            case_summaries[case_id] = summary
            if comparison is not None:
                comparison_cases[case_id] = comparison
            print(
                f"{case_id}: result={summary['result']} "
                f"health={summary['health']['result']} "
                f"t={summary['health'].get('observed_final_time')}"
            )

        campaign = build_campaign_evidence(
            acceptance, policy, case_summaries, comparison_cases
        )
        campaign_path = output / "campaign_evidence.json"
        write_json(campaign_path, campaign)
        required = [str(value) for value in policy["criteria"]["required_cases"]]
        table_paths = write_tables(output, required, case_summaries, campaign)

        summary = {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-acceptance-summary",
            "authority": "non-authorizing-direct-fast-scientific-assessment",
            "result": campaign["result"],
            "selected_cases": selected,
            "required_cases": required,
            "case_results": {
                case_id: value["result"]
                for case_id, value in sorted(case_summaries.items())
            },
            "campaign_evidence": binding(campaign_path),
            "tables": [binding(path) for path in table_paths],
        }
        summary_path = output / "summary.json"
        write_json(summary_path, summary)
        provenance = {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-acceptance-provenance",
            "inputs": {
                "inventory": binding(inventory_path),
                "criteria": policy["criteria_binding"],
                "criteria_review": policy["review_binding"],
                "reviewed_acceptance_utility": binding(ACCEPTANCE_UTILITY),
                "driver": binding(Path(__file__)),
                "case_lineages": {
                    case_id: value.get("provenance", {}).get("lineage")
                    for case_id, value in sorted(case_summaries.items())
                },
                "case_histories": {
                    case_id: value.get("provenance", {}).get("histories")
                    for case_id, value in sorted(case_summaries.items())
                },
            },
            "outputs": {
                "summary": binding(summary_path),
                "campaign_evidence": binding(campaign_path),
                "case_acceptance": {
                    case_id: binding(output / "cases" / case_id / "case_acceptance.json")
                    for case_id in selected
                },
                "reviewed_case_evidence": {
                    case_id: binding(
                        output / "cases" / case_id / "reviewed_case_evidence.json"
                    )
                    for case_id in selected
                    if (output / "cases" / case_id / "reviewed_case_evidence.json").is_file()
                },
                "tables": [binding(path) for path in table_paths],
            },
        }
        write_json(output / "provenance.json", provenance)
        print(
            f"campaign: result={campaign['result']} cases={len(selected)} "
            f"comparison_ready={len(comparison_cases)} output={output}"
        )
    except (
        FastAcceptanceError,
        OSError,
        ValueError,
    ) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
