#!/usr/bin/env python3
"""Validate corrected/composite Stage I products and write a manuscript marker.

This adapter is intentionally validation-only.  The pinned final supervisor
must render the corrected report and exact four-root publication, then invoke
corrected downstream completion before calling this adapter.  ``run`` only
revalidates those retained products and exclusively creates one immutable
``manuscript_ready`` marker.  ``validate`` rechecks the marker and all bindings.
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
import stat
import sys
import tempfile
from typing import Iterable


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
CORRECTED_REPORT_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_corrected_report.py")
CORRECTED_DOWNSTREAM_TOOL = SCRIPT_PATH.with_name(
    "cgl_lf_stage_i_fast_corrected_downstream.py"
)
CORRECTED_SCIENCE_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_corrected_science.py")
PUBLICATION_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_publication.py")

ACTIVE_CASES = (
    "R02",
    "R03",
    "R04",
    "R05",
    "R10",
    "R11",
    "R12",
    "R13",
    "R14",
    "R15",
    "R16",
    "R17",
)
PASSIVE_CASES = ("R06", "R07", "R08", "R09")
ALL_CASES = tuple(f"R{number:02d}" for number in range(2, 18))
ACTIVE_PASSIVE_PAIRS = (
    ("R02", "R06"),
    ("R03", "R07"),
    ("R04", "R08"),
    ("R05", "R09"),
)
ACTIVE_EVIDENCE_CLASS = "corrected_active_production"
PASSIVE_EVIDENCE_CLASS = "authenticated_legacy_passive_control"
VALID_SCIENCE_RESULTS = {"pass", "fail", "inconclusive"}
VALID_CT_RESULTS = {"pass", "fail", "inconclusive"}
SCIENCE_AUTHORITY = "non-authorizing-direct-fast-scientific-assessment"
SIGNED_WORK_SEMANTICS = (
    "applied columns are signed stage ledgers; reconstructed columns are "
    "sparse retained-snapshot estimates and are not applied accounting"
)
SIGNED_WORK_REPORT_STATEMENT = (
    "Signed LF applied-stage ledgers remain distinct from sparse retained-"
    "snapshot pressure-work and heat-flux reconstructions; signs are preserved."
)
REVIEWED_CONTRAST_COLUMNS = (
    "family",
    "contrast",
    "left",
    "right",
    "result",
    "claim_eligible",
    "metric",
    "available",
    "left_mean",
    "right_mean",
    "left_minus_right",
    "right_minus_left",
    "active_mean",
    "passive_mean",
    "active_minus_passive",
    "passive_minus_active",
    "combined_standard_error",
    "pooled_within_realization_standard_deviation",
    "standardized_effect",
    "standardized_effect_scope",
    "signed_standardized_left_minus_right_effect",
    "signed_standardized_active_minus_passive_effect",
    "expected_direction",
    "direction_coherent",
    "large_direction_coherent_effect",
    "intervention_estimand",
    "intervention_enabled_components",
    "excluded_interpretation",
    "intervention_declaration",
    "inference_scope",
    "current_science_scope_disposition",
    "full_scope_independent_review_complete",
    "reason",
    "authority",
    "release_authorizing",
    "claim_scope",
)
ACTIVE_PASSIVE_SUMMARY_COLUMNS = (
    "pair",
    "metric",
    "active",
    "passive",
    "active_minus_passive",
    "passive_minus_active",
    "combined_standard_error",
    "pooled_within_realization_standard_deviation",
    "standardized_effect",
    "standardized_effect_scope",
    "signed_standardized_active_minus_passive_effect",
    "expected_direction",
    "direction_coherent",
    "large_direction_coherent_effect",
    "result",
    "claim_eligible",
    "claim_scope",
    "inference_scope",
    "authority",
    "release_authorizing",
    "reason",
)
ROBUSTNESS_CONTRASTS = (
    ("forcing A, beta=10", "forcing", "R02", "R04"),
    ("forcing A, beta=100", "forcing", "R03", "R05"),
    ("forcing P, beta=10", "forcing", "R06", "R08"),
    ("forcing P, beta=100", "forcing", "R07", "R09"),
    ("beta, A Alfvenic", "beta", "R02", "R03"),
    ("beta, A random", "beta", "R04", "R05"),
    ("beta, P Alfvenic", "beta", "R06", "R07"),
    ("beta, P random", "beta", "R08", "R09"),
    ("forcing correlation", "tcorr", "R05", "R11"),
)
ROBUSTNESS_METRICS = (
    ("kinetic", "kinetic"),
    ("magnetic", "magnetic"),
    ("abs_dp", "abs_dp"),
    ("unstable", "unstable_occupancy"),
    ("nu_eff", "nu_eff"),
)
ROBUSTNESS_SUMMARY_COLUMNS = (
    "contrast",
    "family",
    "reference",
    "variant",
    "metric",
    "reference_value",
    "variant_value",
    "variant_minus_reference",
    "signed_relative_difference",
    "reviewed_left_minus_right",
    "combined_standard_error",
    "pooled_within_realization_standard_deviation",
    "standardized_effect",
    "standardized_effect_scope",
    "signed_standardized_reference_minus_variant_effect",
    "result",
    "claim_eligible",
    "claim_scope",
    "inference_scope",
    "authority",
    "release_authorizing",
    "reason",
)
STANDARD_REVIEWED_CLAIM_SCOPE = "standard reviewed-science scope"
R15_RESTRICTED_CLAIM_SCOPE = (
    "restricted nonfatal-hard-bound diagnostic; not strict R15 success"
)
SIGNED_WORK_COLUMNS = (
    "case_id",
    "availability",
    "diagnostics_provenance",
    "applied_heat_flux_availability",
    "applied_pressure_work_availability",
    "reconstructed_pressure_availability",
    "reconstructed_heat_flux_availability",
    "applied_ledgers_signed",
    "applied_heat_flux_parallel",
    "applied_heat_flux_perpendicular",
    "applied_heat_flux_total",
    "applied_pressure_work_total",
    "applied_pressure_work_anisotropic",
    "cap_parallel_over_1",
    "cap_parallel_over_10",
    "cap_perpendicular_over_1",
    "cap_perpendicular_over_10",
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
    "reconstructed_unlimited_heat_flux_power_integral",
    "semantics",
)

SCHEMA = "athenak-cgl-corrected-composite-manuscript-ready"
SCHEMA_VERSION = 1
RECORD_TYPE = "cgl-lf-stage-i-corrected-composite-manuscript-ready"
RECORD_NAME = "manuscript-ready.json"
DOWNSTREAM_POINTER_SCHEMA = "athenak-cgl-corrected-downstream-pointer"
DOWNSTREAM_COMPLETION_SCHEMA = "athenak-cgl-corrected-downstream-completion"
SCIENCE_SCHEMA = "athenak-cgl-corrected-composite-reviewed-science"
SCIENCE_RECORD_TYPE = "cgl-lf-stage-i-corrected-composite-reviewed-science"
DIRECT_SCIENCE_RECORD_TYPE = "cgl-lf-stage-i-direct-fast-reviewed-science-comparisons"
PUBLICATION_RECORD_TYPE = "cgl_lf_stage_i_fast_publication_products"
CT_RECORD_TYPE = "stage-i-direct-fast-ct-audit"

PUBLICATION_FIGURES = (
    "figures/fig01_health_completion.pdf",
    "figures/fig02_active_passive_history.pdf",
    "figures/fig03_robustness_summary.pdf",
    "figures/fig04_limiter_heat_flux_summary.pdf",
    "figures/fig05_resolution_summary.pdf",
    "figures/fig06_r10_r14_r15_scope.pdf",
    "figures/fig07_reviewed_science_ct_summary.pdf",
    "figures/fig08_causal_mechanism.pdf",
    "figures/fig09_resolution_curves.pdf",
    "figures/fig10_hyperbolicity_coverage.pdf",
)
PUBLICATION_TABLES = (
    "numerical_health_provenance",
    "primary_full_window_scalars",
    "hyperbolicity_all_snapshot_coverage",
    "signed_lf_cap_work_ledger",
    "mks24_panel_dispositions",
    "lineage_dispositions",
    "coherent_direction_mechanism",
    "health_completion",
    "active_passive_summary",
    "robustness_summary",
    "limiter_heat_flux_summary",
    "resolution_summary",
    "r10_r14_r15_scope",
    "case_health_scientific_warnings",
    "acceptance_gates",
    "reviewed_science_contrasts",
    "reviewed_science_gates",
    "reviewed_science_resolution",
    "reviewed_science_mks24",
    "direct_ct_health",
)
REQUIRED_PUBLICATION_PRODUCTS = tuple(
    (
        *PUBLICATION_FIGURES,
        "captions.md",
        "report.md",
        *(
            f"tables/{name}.{suffix}"
            for name in PUBLICATION_TABLES
            for suffix in ("csv", "tex")
        ),
    )
)
MANUSCRIPT_PRODUCTS = (
    "manuscript/results.json",
    "manuscript/results_macros.tex",
    "manuscript/claim_evidence.json",
    "manuscript/claim_evidence.md",
    "manuscript/figure_manifest.json",
    "manuscript/provenance_manifest.json",
    "manuscript/open_questions.md",
    "verify.json",
)

DEFAULT_CAMPAIGN = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/campaigns/"
    "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
)
DEFAULT_IDENTITY = DEFAULT_CAMPAIGN / "campaign-identity.json"
DEFAULT_INVENTORY = DEFAULT_CAMPAIGN / (
    "analysis/mks24-stage-i-composite/E03-forcing-policy/R02-R17-final/inventory.json"
)
DEFAULT_ANALYSIS = DEFAULT_INVENTORY.parent
DEFAULT_ACCEPTANCE = DEFAULT_ANALYSIS.with_name(f"{DEFAULT_ANALYSIS.name}-acceptance")
DEFAULT_SCIENCE = DEFAULT_ANALYSIS.with_name(f"{DEFAULT_ANALYSIS.name}-science")
DEFAULT_WORKFLOW = DEFAULT_CAMPAIGN / (
    "analysis/mks24-stage-i-corrected-downstream/E03-forcing-policy/R02-R17-final"
)
DEFAULT_MARKER = DEFAULT_WORKFLOW / RECORD_NAME
DEFAULT_CRITERIA = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
DEFAULT_CRITERIA_REVIEW = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)


class ManuscriptReadyError(RuntimeError):
    """The corrected/composite manuscript-ready data gate failed."""


def load_module(name: str, path: Path):
    """Import one reviewed local validation dependency by exact path."""

    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ManuscriptReadyError(f"cannot import validation dependency: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def require_dict(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ManuscriptReadyError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list):
        raise ManuscriptReadyError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ManuscriptReadyError(f"{label} must be a nonempty string")
    return value


def require_bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ManuscriptReadyError(f"{label} must be a boolean")
    return value


def require_finite(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ManuscriptReadyError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ManuscriptReadyError(f"{label} must be a finite number")
    return result


def require_close(observed: object, expected: float, label: str) -> float:
    current = require_finite(observed, label)
    if not math.isclose(current, expected, rel_tol=1.0e-12, abs_tol=1.0e-15):
        raise ManuscriptReadyError(
            f"{label} differs: observed={current!r}, expected={expected!r}"
        )
    return current


def nested(value: object, path: str) -> object:
    current = value
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def require_exact_scope_records(
    value: dict[str, object],
    limitation: dict[str, object],
    intervention: dict[str, object],
    label: str,
) -> None:
    """Require one release input to carry both exact reviewed-science scopes."""

    if value.get("current_science_scope_limitation") != limitation:
        raise ManuscriptReadyError(f"{label} current science scope limitation differs")
    if value.get("active_passive_intervention_scope") != intervention:
        raise ManuscriptReadyError(f"{label} active/passive intervention scope differs")


def publication_claim_scope(left: object, right: object) -> str:
    """Return the exact reviewed-science scope rendered for this campaign."""

    return (
        R15_RESTRICTED_CLAIM_SCOPE
        if "R15" in {str(left), str(right)}
        else STANDARD_REVIEWED_CLAIM_SCOPE
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_binding(path: Path) -> dict[str, object]:
    """Return a stable binding for one current regular file."""

    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise ManuscriptReadyError(f"artifact does not resolve: {path}") from error
    before = resolved.stat()
    if not stat.S_ISREG(before.st_mode):
        raise ManuscriptReadyError(f"artifact is not a regular file: {resolved}")
    digest = sha256_file(resolved)
    after = resolved.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise ManuscriptReadyError(f"artifact changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": after.st_size,
        "sha256": digest,
    }


def verify_binding(value: object, label: str) -> dict[str, object]:
    declared = require_dict(value, f"{label} binding")
    current = artifact_binding(Path(require_text(declared.get("path"), f"{label} path")))
    if (
        declared.get("sha256") != current["sha256"]
        or declared.get("size_bytes") != current["size_bytes"]
    ):
        raise ManuscriptReadyError(f"{label} differs from its declared binding")
    return current


def require_same_binding(
    value: object, expected: object, label: str
) -> dict[str, object]:
    observed = verify_binding(value, label)
    authority = verify_binding(expected, f"expected {label}")
    if observed != authority:
        raise ManuscriptReadyError(f"{label} differs from its authority")
    return observed


def load_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ManuscriptReadyError(f"cannot load {label}: {path}") from error
    return require_dict(value, label)


def load_bound_json(
    path: Path, label: str
) -> tuple[dict[str, object], dict[str, object]]:
    binding = artifact_binding(path)
    return load_json(Path(str(binding["path"])), label), binding


def stable_json(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def publication_text(value: object, *, csv_value: bool = False) -> str:
    if value is None:
        return "--"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "--"
        return repr(value) if csv_value else f"{value:.8g}"
    if isinstance(value, (list, tuple)):
        return "; ".join(
            publication_text(item, csv_value=csv_value) for item in value
        )
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return str(value)


def tex_unescape(value: str) -> str:
    placeholder = "\0"
    result = value.replace(r"\textbackslash{}", placeholder)
    for encoded, plain in (
        (r"\&", "&"),
        (r"\%", "%"),
        (r"\_", "_"),
        (r"\#", "#"),
        (r"\$", "$"),
        (r"\{", "{"),
        (r"\}", "}"),
    ):
        result = result.replace(encoded, plain)
    return result.replace(placeholder, "\\")


def load_csv_table(
    path: Path, columns: tuple[str, ...], label: str
) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != list(columns):
                raise ManuscriptReadyError(f"{label} columns differ")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeError, csv.Error) as error:
        raise ManuscriptReadyError(f"cannot load {label}: {path}") from error
    return rows


def load_tex_table(
    path: Path, columns: tuple[str, ...], label: str
) -> list[dict[str, str]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise ManuscriptReadyError(f"cannot load {label}: {path}") from error
    if len(lines) < 6 or lines[1] != r"\hline" or lines[3] != r"\hline":
        raise ManuscriptReadyError(f"{label} layout differs")

    def cells(line: str) -> list[str]:
        if not line.endswith(r" \\"):
            raise ManuscriptReadyError(f"{label} row terminator differs")
        return [tex_unescape(value) for value in line[:-3].split(" & ")]

    if cells(lines[2]) != list(columns):
        raise ManuscriptReadyError(f"{label} columns differ")
    if lines[-2:] != [r"\hline", r"\end{tabular}"]:
        raise ManuscriptReadyError(f"{label} footer differs")
    result = []
    for index, line in enumerate(lines[4:-2]):
        values = cells(line)
        if len(values) != len(columns):
            raise ManuscriptReadyError(f"{label} row {index} width differs")
        result.append(dict(zip(columns, values)))
    return result


def canonical_directory(path: Path, label: str) -> Path:
    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise ManuscriptReadyError(f"{label} does not resolve: {path}") from error
    if not resolved.is_dir():
        raise ManuscriptReadyError(f"{label} is not a directory: {resolved}")
    return resolved


def path_contains(parent: Path, child: Path) -> bool:
    parent = parent.resolve(strict=False)
    child = child.resolve(strict=False)
    return child == parent or parent in child.parents


def recursive_file_bindings(value: object) -> Iterable[dict[str, object]]:
    if isinstance(value, dict):
        if (
            isinstance(value.get("path"), str)
            and isinstance(value.get("sha256"), str)
            and isinstance(value.get("size_bytes"), int)
        ):
            yield value
        for child in value.values():
            yield from recursive_file_bindings(child)
    elif isinstance(value, list):
        for child in value:
            yield from recursive_file_bindings(child)


def expected_science_classification() -> dict[str, list[str]]:
    return {
        "corrected_active": list(ACTIVE_CASES),
        "authenticated_legacy_passive": list(PASSIVE_CASES),
        "selected": list(ALL_CASES),
    }


def expected_downstream_classification() -> dict[str, list[str]]:
    return {
        "corrected_active": list(ACTIVE_CASES),
        "reused_passive": list(PASSIVE_CASES),
        "selected": list(ALL_CASES),
    }


def marker_path(path: Path) -> Path:
    raw = path.expanduser().absolute()
    if raw.is_symlink():
        raise ManuscriptReadyError(f"manuscript-ready marker may not be a symlink: {raw}")
    return raw.resolve(strict=False)


def refuse_existing_marker(path: Path) -> Path:
    """Fail before any validation or output mutation when the marker exists."""

    resolved = marker_path(path)
    if resolved.exists():
        raise ManuscriptReadyError(
            f"refusing run because manuscript-ready marker already exists: {resolved}"
        )
    return resolved


def write_new_marker(path: Path, value: object) -> None:
    """Exclusively create the only output this adapter is permitted to mutate."""

    path = marker_path(path)
    parent = canonical_directory(path.parent, "manuscript-ready marker parent")
    if path.parent.resolve(strict=True) != parent:
        raise ManuscriptReadyError("manuscript-ready marker parent differs")
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            staged = Path(stream.name)
            stream.write(stable_json(value))
            stream.flush()
            os.fsync(stream.fileno())
        os.link(staged, path)
    except FileExistsError as error:
        raise ManuscriptReadyError(
            f"manuscript-ready marker appeared concurrently: {path}"
        ) from error
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def validate_corrected_context(
    args: argparse.Namespace,
) -> tuple[object, dict[str, object]]:
    """Authenticate the exact corrected-active/legacy-passive composite."""

    science = load_module("_cgl_manuscript_ready_science", CORRECTED_SCIENCE_TOOL)
    try:
        context = science.validate_corrected_context(
            args.identity, args.inventory, args.inventory_sha256
        )
    except Exception as error:
        raise ManuscriptReadyError(
            f"corrected composite inventory validation failed: {error}"
        ) from error
    if (
        context.get("inventory_kind") != "corrected-composite"
        or context.get("active_cases") != list(ACTIVE_CASES)
        or context.get("passive_cases") != list(PASSIVE_CASES)
        or context.get("selected_cases") != list(ALL_CASES)
    ):
        raise ManuscriptReadyError(
            "manuscript-ready gate requires exact corrected-active and "
            "authenticated-legacy-passive R02-R17"
        )
    inventory = require_dict(context.get("inventory"), "composite inventory")
    expected_classes = {
        ACTIVE_EVIDENCE_CLASS: list(ACTIVE_CASES),
        PASSIVE_EVIDENCE_CLASS: list(PASSIVE_CASES),
    }
    if (
        inventory.get("record_type") != "cgl_lf_stage_i_corrected_composite_report"
        or inventory.get("evidence_classes") != expected_classes
    ):
        raise ManuscriptReadyError("corrected composite inventory identity differs")
    cases = require_dict(inventory.get("cases"), "composite inventory cases")
    if set(cases) != set(ALL_CASES):
        raise ManuscriptReadyError("composite inventory does not cover exact R02-R17")
    for case_id in ALL_CASES:
        case = require_dict(cases[case_id], f"{case_id} composite case")
        authority = require_dict(
            case.get("execution_authority"), f"{case_id} execution authority"
        )
        expected = (
            ACTIVE_EVIDENCE_CLASS if case_id in ACTIVE_CASES else PASSIVE_EVIDENCE_CLASS
        )
        if (
            case.get("status") != "complete"
            or case.get("evidence_class") != expected
            or authority.get("evidence_class") != expected
        ):
            raise ManuscriptReadyError(
                f"{case_id} is not complete evidence from its exact authority"
            )
    return science, context


def reviewed_contrast_authority(
    direct: dict[str, object], intervention: dict[str, object]
) -> dict[tuple[str, str, str], dict[str, object]]:
    """Validate and flatten the reviewed contrast values used by claim products."""

    if (
        direct.get("authority") != SCIENCE_AUTHORITY
        or direct.get("release_authorizing") is not False
    ):
        raise ManuscriptReadyError(
            "direct reviewed science lacks its exact non-authorizing authority"
        )
    families = require_dict(direct.get("families"), "direct reviewed science families")
    records: dict[tuple[str, str, str], dict[str, object]] = {}
    for family, contrasts_value in sorted(families.items()):
        contrasts = require_dict(
            contrasts_value, f"direct reviewed science family {family}"
        )
        for contrast_name, contrast_value in sorted(contrasts.items()):
            contrast = require_dict(
                contrast_value,
                f"direct reviewed science contrast {family}.{contrast_name}",
            )
            left = contrast.get("active", contrast.get("left"))
            right = contrast.get("passive", contrast.get("right"))
            if left not in ALL_CASES or right not in ALL_CASES:
                raise ManuscriptReadyError(
                    f"direct reviewed science contrast {family}.{contrast_name} "
                    "has invalid case semantics"
                )
            if family == "active_passive":
                if [left, right] not in [list(pair) for pair in ACTIVE_PASSIVE_PAIRS]:
                    raise ManuscriptReadyError(
                        f"direct reviewed science contrast {contrast_name} "
                        "has invalid active/passive ordering"
                    )
                if contrast.get("intervention_scope") != intervention:
                    raise ManuscriptReadyError(
                        f"direct reviewed science contrast {contrast_name} "
                        "intervention scope differs"
                    )
            metrics = require_list(
                contrast.get("metrics"),
                f"direct reviewed science contrast {family}.{contrast_name} metrics",
            )
            if not metrics:
                raise ManuscriptReadyError(
                    f"direct reviewed science contrast {family}.{contrast_name} "
                    "has no reviewed metrics"
                )
            for index, metric_value in enumerate(metrics):
                label = (
                    f"direct reviewed science contrast {family}.{contrast_name} "
                    f"metric {index}"
                )
                metric = require_dict(metric_value, label)
                metric_name = require_text(metric.get("metric"), f"{label} name")
                available = require_bool(metric.get("available"), f"{label} availability")
                key = (str(family), str(contrast_name), metric_name)
                if key in records:
                    raise ManuscriptReadyError(
                        f"direct reviewed science contains duplicate metric {key}"
                    )
                record = {
                    "family": str(family),
                    "contrast": str(contrast_name),
                    "left": str(left),
                    "right": str(right),
                    "result": require_text(
                        contrast.get("result"), f"{label} contrast result"
                    ),
                    "claim_eligible": require_bool(
                        contrast.get("claim_eligible"), f"{label} claim eligibility"
                    ),
                    "metric": metric_name,
                    "available": available,
                    "reason": metric.get("reason", contrast.get("reason")),
                    "source_claim_scope": metric.get("claim_scope"),
                    "claim_scope": publication_claim_scope(left, right),
                    "inference_scope": (
                        intervention.get("claim_scope")
                        if family == "active_passive"
                        else "descriptive_within_realization"
                    ),
                }
                if available:
                    left_key, right_key = (
                        ("active_mean", "passive_mean")
                        if family == "active_passive"
                        else ("left_mean", "right_mean")
                    )
                    left_mean = require_finite(metric.get(left_key), f"{label} {left_key}")
                    right_mean = require_finite(
                        metric.get(right_key), f"{label} {right_key}"
                    )
                    difference = require_finite(
                        metric.get(
                            "difference",
                            metric.get("difference_left_minus_right"),
                        ),
                        f"{label} left-minus-right difference",
                    )
                    require_close(
                        difference,
                        left_mean - right_mean,
                        f"{label} left-minus-right sign convention",
                    )
                    combined_error = require_finite(
                        metric.get("combined_standard_error"),
                        f"{label} combined standard error",
                    )
                    if combined_error < 0.0:
                        raise ManuscriptReadyError(
                            f"{label} combined standard error is negative"
                        )
                    pooled = require_finite(
                        metric.get(
                            "pooled_within_realization_standard_deviation"
                        ),
                        f"{label} pooled standard deviation",
                    )
                    if pooled < 0.0:
                        raise ManuscriptReadyError(
                            f"{label} pooled standard deviation is negative"
                        )
                    effect = require_finite(
                        metric.get("standardized_effect"),
                        f"{label} standardized effect",
                    )
                    expected_effect = (
                        0.0
                        if difference == 0.0
                        else sys.float_info.max
                        if pooled == 0.0
                        else abs(difference / pooled)
                    )
                    require_close(
                        effect,
                        expected_effect,
                        f"{label} standardized-effect magnitude semantics",
                    )
                    if effect < 0.0:
                        raise ManuscriptReadyError(
                            f"{label} standardized effect is not a magnitude"
                        )
                    signed_effect = (
                        math.copysign(effect, difference)
                        if difference != 0.0 else 0.0
                    )
                    record.update(
                        {
                            "left_mean": left_mean,
                            "right_mean": right_mean,
                            "difference": difference,
                            "combined_standard_error": combined_error,
                            "pooled_within_realization_standard_deviation": pooled,
                            "standardized_effect": effect,
                            "standardized_effect_scope": require_text(
                                metric.get("standardized_effect_scope"),
                                f"{label} standardized-effect scope",
                            ),
                            "signed_standardized_effect": signed_effect,
                        }
                    )
                    if (
                        record["standardized_effect_scope"]
                        != "descriptive_within_realization"
                    ):
                        raise ManuscriptReadyError(
                            f"{label} standardized-effect scope differs"
                        )
                    if family == "active_passive":
                        expected_direction = require_text(
                            metric.get("expected_direction"),
                            f"{label} expected direction",
                        )
                        coherent = require_bool(
                            metric.get("direction_coherent"),
                            f"{label} direction coherence",
                        )
                        expected_coherent = (
                            difference < 0.0
                            if expected_direction == "active_lower"
                            else difference > 0.0
                            if expected_direction == "active_higher"
                            else None
                        )
                        if expected_coherent is None or coherent != expected_coherent:
                            raise ManuscriptReadyError(
                                f"{label} expected-direction semantics differ"
                            )
                        record.update(
                            {
                                "expected_direction": expected_direction,
                                "direction_coherent": coherent,
                                "large_direction_coherent_effect": require_bool(
                                    metric.get("large_direction_coherent_effect"),
                                    f"{label} large-effect coherence",
                                ),
                            }
                        )
                records[key] = record
    if not records:
        raise ManuscriptReadyError("direct reviewed science has no contrast authority")
    return records


def validate_reviewed_science(
    args: argparse.Namespace, science: object, context: dict[str, object]
) -> dict[str, object]:
    """Revalidate corrected/composite reviewed science without judging its result."""

    acceptance = canonical_directory(args.acceptance_output, "acceptance output")
    science_root = canonical_directory(args.science_output, "reviewed science output")
    criteria = args.criteria.expanduser().resolve(strict=True)
    review = args.criteria_review.expanduser().resolve(strict=True)
    try:
        expected = science.validate_products(
            context, acceptance, science_root, criteria, review
        )
        limitation, intervention = science.validated_scope_records(criteria, review)
    except Exception as error:
        raise ManuscriptReadyError(
            f"corrected/composite reviewed science validation failed: {error}"
        ) from error
    record_path = science_root / "corrected-composite-science.json"
    record, record_binding = load_bound_json(
        record_path, "corrected/composite reviewed-science record"
    )
    if record != expected:
        raise ManuscriptReadyError(
            "corrected/composite reviewed-science record differs from current evidence"
        )
    if (
        record.get("schema") != SCIENCE_SCHEMA
        or record.get("schema_version") != 1
        or record.get("record_type") != SCIENCE_RECORD_TYPE
        or record.get("status") != "complete"
        or record.get("campaign_kind") != "corrected-composite"
        or record.get("result") not in VALID_SCIENCE_RESULTS
        or record.get("case_classification") != expected_science_classification()
    ):
        raise ManuscriptReadyError(
            "corrected/composite reviewed-science identity or coverage differs"
        )
    require_exact_scope_records(
        record, limitation, intervention, "corrected/composite reviewed science"
    )
    require_same_binding(
        record.get("campaign_identity"), context["identity_binding"], "science identity"
    )
    require_same_binding(
        record.get("inventory"), context["inventory_binding"], "science inventory"
    )
    acceptance_record = require_dict(record.get("acceptance"), "science acceptance")
    products = require_dict(record.get("science"), "science products")
    if (
        Path(require_text(acceptance_record.get("root"), "acceptance root")).resolve()
        != acceptance
        or Path(require_text(products.get("root"), "science root")).resolve()
        != science_root
    ):
        raise ManuscriptReadyError("reviewed-science output roots differ")
    direct_binding = verify_binding(products.get("record"), "direct reviewed science")
    direct = load_json(Path(str(direct_binding["path"])), "direct reviewed science")
    if (
        Path(str(direct_binding["path"])) != science_root / "science.json"
        or direct.get("record_type") != DIRECT_SCIENCE_RECORD_TYPE
        or direct.get("selected_cases") != list(ALL_CASES)
        or direct.get("result") != record["result"]
    ):
        raise ManuscriptReadyError("direct reviewed-science record differs")
    require_exact_scope_records(
        direct, limitation, intervention, "direct reviewed science"
    )
    contrasts = reviewed_contrast_authority(direct, intervention)
    return {
        "root": str(science_root),
        "acceptance_root": str(acceptance),
        "active_passive_intervention_scope": intervention,
        "current_science_scope_limitation": limitation,
        "record": record,
        "record_binding": record_binding,
        "direct": direct,
        "direct_record": direct_binding,
        "contrasts": contrasts,
    }


def validate_downstream_completion(
    context: dict[str, object], workflow_root: Path
) -> dict[str, object]:
    """Revalidate and inspect the supervisor-created downstream completion."""

    downstream = load_module(
        "_cgl_manuscript_ready_downstream", CORRECTED_DOWNSTREAM_TOOL
    )
    root = canonical_directory(workflow_root, "corrected downstream workflow root")
    try:
        pointer_path = Path(downstream.validate_completion(root)).resolve(strict=True)
    except Exception as error:
        raise ManuscriptReadyError(
            f"corrected downstream completion validation failed: {error}"
        ) from error
    expected_pointer = root / "completion/corrected-production-complete.json"
    if pointer_path != expected_pointer:
        raise ManuscriptReadyError("corrected downstream completion pointer path differs")
    pointer, pointer_binding = load_bound_json(pointer_path, "downstream completion pointer")
    if (
        pointer.get("schema") != DOWNSTREAM_POINTER_SCHEMA
        or pointer.get("schema_version") != 1
        or pointer.get("status") != "complete"
        or pointer.get("campaign_kind") != "corrected-production"
    ):
        raise ManuscriptReadyError("corrected downstream completion pointer differs")
    require_same_binding(
        pointer.get("campaign_identity"), context["identity_binding"], "pointer identity"
    )
    require_same_binding(
        pointer.get("inventory"), context["inventory_binding"], "pointer inventory"
    )

    completion_binding = verify_binding(
        pointer.get("completion_record"), "downstream completion record"
    )
    completion = load_json(
        Path(str(completion_binding["path"])), "downstream completion record"
    )
    if (
        completion.get("schema") != DOWNSTREAM_COMPLETION_SCHEMA
        or completion.get("schema_version") != 1
        or completion.get("status") != "complete"
        or completion.get("campaign_kind") != "corrected-production"
        or completion.get("case_classification") != expected_downstream_classification()
    ):
        raise ManuscriptReadyError("corrected downstream completion record differs")
    require_same_binding(
        completion.get("campaign_identity"),
        context["identity_binding"],
        "completion identity",
    )
    require_same_binding(
        completion.get("inventory"), context["inventory_binding"], "completion inventory"
    )
    analysis = require_dict(completion.get("analysis"), "completion analysis")
    hyper = require_dict(completion.get("hyperbolicity"), "completion hyperbolicity")
    if set(analysis) != set(ALL_CASES) or set(hyper) != set(ACTIVE_CASES):
        raise ManuscriptReadyError("downstream completion coverage differs")

    ct = require_dict(completion.get("ct"), "completion CT")
    ct_binding = verify_binding(ct.get("audit"), "completion CT audit")
    ct_record = load_json(Path(str(ct_binding["path"])), "completion CT audit")
    selection = require_dict(ct_record.get("selection"), "CT audit selection")
    cases = require_dict(ct_record.get("cases"), "CT audit cases")
    if (
        ct_record.get("record_type") != CT_RECORD_TYPE
        or ct_record.get("result") not in VALID_CT_RESULTS
        or selection.get("cases") != list(ALL_CASES)
        or selection.get("snapshot_policy") != "all"
        or set(cases) != set(ALL_CASES)
    ):
        raise ManuscriptReadyError("CT audit identity or exact coverage differs")
    for case_id, value in cases.items():
        case = require_dict(value, f"{case_id} CT audit")
        native = require_dict(case.get("native_restart_ct"), f"{case_id} native CT")
        if (
            case.get("provenance_authenticated") is not True
            or native.get("coverage_complete") is not True
        ):
            raise ManuscriptReadyError(
                f"{case_id} CT authentication or native-restart coverage is incomplete"
            )

    publication = require_dict(completion.get("publication"), "completion publication")
    publication_binding = verify_binding(
        publication.get("manifest"), "completion publication manifest"
    )
    exact_manifest = root / "publication/manifest.json"
    if Path(str(publication_binding["path"])) != exact_manifest:
        raise ManuscriptReadyError(
            "completion publication manifest is not exact workflow_root/publication"
        )
    return {
        "root": str(root),
        "pointer": pointer_binding,
        "completion_record": completion_binding,
        "completion": completion,
        "analysis": analysis,
        "hyperbolicity": hyper,
        "ct_audit": ct_binding,
        "ct_record": ct_record,
        "publication_manifest": publication_binding,
    }


def publication_roots(
    context: dict[str, object],
    downstream: dict[str, object],
    science: dict[str, object],
) -> tuple[Path, Path, Path, Path]:
    analysis = canonical_directory(
        Path(str(context["inventory_output"])), "composite report"
    )
    workflow = canonical_directory(Path(str(downstream["root"])), "workflow root")
    return (
        canonical_directory(Path(str(science["acceptance_root"])), "acceptance root"),
        canonical_directory(Path(str(science["root"])), "science root"),
        canonical_directory(workflow / "ct", "workflow CT root"),
        canonical_directory(
            analysis / "corrected-downstream/hyperbolicity",
            "composite hyperbolicity root",
        ),
    )


def validated_manifest_bindings(
    values: object, label: str, *, root: Path | None = None
) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for index, value in enumerate(require_list(values, label)):
        current = verify_binding(value, f"{label} {index}")
        path = Path(str(current["path"]))
        if root is not None and not path_contains(root, path):
            raise ManuscriptReadyError(f"{label} contains a file outside {root}: {path}")
        if str(path) in records:
            raise ManuscriptReadyError(f"{label} contains duplicate path: {path}")
        records[str(path)] = current
    if not records:
        raise ManuscriptReadyError(f"{label} is empty")
    return records


def source_record_paths(
    sources: Iterable[dict[str, object]],
) -> dict[str, set[str]]:
    records: dict[str, set[str]] = {}
    for source in sources:
        path = Path(str(source["path"]))
        if path.suffix != ".json":
            continue
        record = load_json(path, f"publication JSON source {path}")
        record_type = record.get("record_type")
        if isinstance(record_type, str):
            records.setdefault(record_type, set()).add(str(path))
    return records


def table_index(
    rows: list[dict[str, str]], keys: tuple[str, ...], label: str
) -> dict[tuple[str, ...], dict[str, str]]:
    result: dict[tuple[str, ...], dict[str, str]] = {}
    for index, row in enumerate(rows):
        key = tuple(row.get(name, "") for name in keys)
        if not all(key) or key in result:
            raise ManuscriptReadyError(f"{label} row {index} has invalid key {key}")
        result[key] = row
    return result


def require_table_text(
    row: dict[str, str], field: str, expected: object, label: str, *, tex: bool = False
) -> None:
    current = row.get(field)
    wanted = publication_text(expected, csv_value=not tex)
    if current != wanted:
        raise ManuscriptReadyError(
            f"{label} {field} differs: observed={current!r}, expected={wanted!r}"
        )


def require_table_number(
    row: dict[str, str], field: str, expected: float, label: str
) -> None:
    value = row.get(field)
    try:
        current = float(value) if value not in (None, "", "--") else math.nan
    except ValueError as error:
        raise ManuscriptReadyError(f"{label} {field} is not numeric") from error
    require_close(current, expected, f"{label} {field}")


def require_table_value(
    row: dict[str, str],
    field: str,
    expected: object,
    label: str,
    *,
    tex: bool,
) -> None:
    """Compare one deterministic publication cell to its reconstructed value."""

    if (
        not tex
        and isinstance(expected, (int, float))
        and not isinstance(expected, bool)
    ):
        require_table_number(row, field, float(expected), label)
    else:
        require_table_text(row, field, expected, label, tex=tex)


def reviewed_publication_row(
    record: dict[str, object],
    intervention: dict[str, object],
    limitation: dict[str, object],
) -> dict[str, object]:
    """Reconstruct one reviewed-science publication row from its authority."""

    available = record["available"] is True
    active_passive = record["family"] == "active_passive"
    difference = float(record["difference"]) if available else None
    signed_effect = (
        float(record["signed_standardized_effect"]) if available else None
    )
    return {
        "family": record["family"],
        "contrast": record["contrast"],
        "left": record["left"],
        "right": record["right"],
        "result": record["result"],
        "claim_eligible": record["claim_eligible"],
        "metric": record["metric"],
        "available": record["available"],
        "left_mean": record.get("left_mean"),
        "right_mean": record.get("right_mean"),
        "left_minus_right": difference,
        "right_minus_left": -difference if difference is not None else None,
        "active_mean": record.get("left_mean") if active_passive else None,
        "passive_mean": record.get("right_mean") if active_passive else None,
        "active_minus_passive": difference if active_passive else None,
        "passive_minus_active": (
            -difference if active_passive and difference is not None else None
        ),
        "combined_standard_error": record.get("combined_standard_error"),
        "pooled_within_realization_standard_deviation": record.get(
            "pooled_within_realization_standard_deviation"
        ),
        "standardized_effect": record.get("standardized_effect"),
        "standardized_effect_scope": record.get("standardized_effect_scope"),
        "signed_standardized_left_minus_right_effect": signed_effect,
        "signed_standardized_active_minus_passive_effect": (
            signed_effect if active_passive else None
        ),
        "expected_direction": record.get("expected_direction"),
        "direction_coherent": record.get("direction_coherent"),
        "large_direction_coherent_effect": record.get(
            "large_direction_coherent_effect"
        ),
        "intervention_estimand": (
            intervention.get("estimand") if active_passive else None
        ),
        "intervention_enabled_components": (
            intervention.get("enabled_components") if active_passive else None
        ),
        "excluded_interpretation": (
            intervention.get("excluded_interpretation") if active_passive else None
        ),
        "intervention_declaration": (
            intervention.get("declaration") if active_passive else None
        ),
        "inference_scope": record["inference_scope"],
        "current_science_scope_disposition": limitation.get("disposition"),
        "full_scope_independent_review_complete": limitation.get(
            "full_scope_independent_review_complete"
        ),
        "reason": record.get("reason"),
        "authority": SCIENCE_AUTHORITY,
        "release_authorizing": False,
        "claim_scope": record["claim_scope"],
    }


def validate_reviewed_contrast_table(
    root: Path,
    science: dict[str, object],
) -> None:
    """Require CSV/TeX reviewed claims to preserve the authority's exact semantics."""

    csv_rows = load_csv_table(
        root / "tables/reviewed_science_contrasts.csv",
        REVIEWED_CONTRAST_COLUMNS,
        "publication reviewed-science CSV",
    )
    tex_rows = load_tex_table(
        root / "tables/reviewed_science_contrasts.tex",
        REVIEWED_CONTRAST_COLUMNS,
        "publication reviewed-science TeX",
    )
    keys = ("family", "contrast", "metric")
    csv_index = table_index(csv_rows, keys, "publication reviewed-science CSV")
    tex_index = table_index(tex_rows, keys, "publication reviewed-science TeX")
    authority = science["contrasts"]
    if set(csv_index) != set(authority) or set(tex_index) != set(authority):
        raise ManuscriptReadyError(
            "publication reviewed-science metric inventory differs from authority"
        )
    intervention = science["active_passive_intervention_scope"]
    limitation = science["current_science_scope_limitation"]
    for key, expected in authority.items():
        values = reviewed_publication_row(expected, intervention, limitation)
        for rows, tex, label in (
            (csv_index, False, "publication reviewed-science CSV"),
            (tex_index, True, "publication reviewed-science TeX"),
        ):
            row = rows[key]
            for field in REVIEWED_CONTRAST_COLUMNS:
                require_table_value(
                    row, field, values[field], f"{label} {key}", tex=tex
                )


def active_passive_publication_metrics(
    science: dict[str, object],
) -> dict[tuple[str, str], dict[str, object]]:
    result: dict[tuple[str, str], dict[str, object]] = {}
    for record in science["contrasts"].values():
        if record["family"] != "active_passive":
            continue
        key = (f"{record['left']}/{record['right']}", str(record["metric"]))
        result[key] = record
    return result


def validate_active_passive_summary_table(
    root: Path, science: dict[str, object]
) -> None:
    """Reject stale fast-report values and passive-minus-active sign drift."""

    expected = active_passive_publication_metrics(science)
    if not expected:
        raise ManuscriptReadyError(
            "reviewed science has no active/passive scalar claim authority"
        )
    csv_rows = load_csv_table(
        root / "tables/active_passive_summary.csv",
        ACTIVE_PASSIVE_SUMMARY_COLUMNS,
        "publication active/passive CSV",
    )
    tex_rows = load_tex_table(
        root / "tables/active_passive_summary.tex",
        ACTIVE_PASSIVE_SUMMARY_COLUMNS,
        "publication active/passive TeX",
    )
    keys = ("pair", "metric")
    csv_index = table_index(csv_rows, keys, "publication active/passive CSV")
    tex_index = table_index(tex_rows, keys, "publication active/passive TeX")
    if set(csv_index) != set(expected) or set(tex_index) != set(expected):
        raise ManuscriptReadyError(
            "publication active/passive summary inventory differs from reviewed science"
        )
    for key, record in expected.items():
        available = record["available"] is True
        difference = float(record["difference"]) if available else None
        values = {
            "pair": key[0],
            "metric": key[1],
            "active": record.get("left_mean"),
            "passive": record.get("right_mean"),
            "active_minus_passive": difference,
            "passive_minus_active": (
                -difference if difference is not None else None
            ),
            "combined_standard_error": record.get("combined_standard_error"),
            "pooled_within_realization_standard_deviation": record.get(
                "pooled_within_realization_standard_deviation"
            ),
            "standardized_effect": record.get("standardized_effect"),
            "standardized_effect_scope": record.get("standardized_effect_scope"),
            "signed_standardized_active_minus_passive_effect": record.get(
                "signed_standardized_effect"
            ),
            "expected_direction": record.get("expected_direction"),
            "direction_coherent": record.get("direction_coherent"),
            "large_direction_coherent_effect": record.get(
                "large_direction_coherent_effect"
            ),
            "result": record["result"],
            "claim_eligible": record["claim_eligible"],
            "claim_scope": record["claim_scope"],
            "inference_scope": record["inference_scope"],
            "authority": SCIENCE_AUTHORITY,
            "release_authorizing": False,
            "reason": record.get("reason"),
        }
        for rows, tex, label in (
            (csv_index, False, "publication active/passive CSV"),
            (tex_index, True, "publication active/passive TeX"),
        ):
            row = rows[key]
            for field in ACTIVE_PASSIVE_SUMMARY_COLUMNS:
                require_table_value(
                    row, field, values[field], f"{label} {key}", tex=tex
                )


def robustness_publication_metrics(
    science: dict[str, object],
) -> dict[tuple[str, str], dict[str, object]]:
    """Reconstruct the exact reviewed robustness rows rendered for publication."""

    authority = science["contrasts"]
    result: dict[tuple[str, str], dict[str, object]] = {}
    for label, family, reference, variant in ROBUSTNESS_CONTRASTS:
        contrast = f"{reference}_{variant}"
        for metric, source_metric in ROBUSTNESS_METRICS:
            source = authority.get((family, contrast, source_metric))
            if not isinstance(source, dict):
                raise ManuscriptReadyError(
                    "reviewed science omits required robustness metric "
                    f"{family}.{contrast}.{source_metric}"
                )
            available = source["available"] is True
            left = float(source["left_mean"]) if available else None
            right = float(source["right_mean"]) if available else None
            difference = float(source["difference"]) if available else None
            scale = (
                max(abs(left), abs(right), 1.0e-300)
                if left is not None and right is not None
                else None
            )
            result[(label, metric)] = {
                "contrast": label,
                "family": family,
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
                "reviewed_left_minus_right": difference,
                "combined_standard_error": source.get("combined_standard_error"),
                "pooled_within_realization_standard_deviation": source.get(
                    "pooled_within_realization_standard_deviation"
                ),
                "standardized_effect": source.get("standardized_effect"),
                "standardized_effect_scope": source.get(
                    "standardized_effect_scope"
                ),
                "signed_standardized_reference_minus_variant_effect": source.get(
                    "signed_standardized_effect"
                ),
                "result": source["result"],
                "claim_eligible": source["claim_eligible"],
                "claim_scope": source["claim_scope"],
                "inference_scope": source["inference_scope"],
                "authority": SCIENCE_AUTHORITY,
                "release_authorizing": False,
                "reason": source.get("reason"),
            }
    return result


def validate_robustness_summary_table(
    root: Path, science: dict[str, object]
) -> None:
    """Reject diagnostic fallback values in claim-facing sensitivity summaries."""

    expected = robustness_publication_metrics(science)
    csv_rows = load_csv_table(
        root / "tables/robustness_summary.csv",
        ROBUSTNESS_SUMMARY_COLUMNS,
        "publication robustness CSV",
    )
    tex_rows = load_tex_table(
        root / "tables/robustness_summary.tex",
        ROBUSTNESS_SUMMARY_COLUMNS,
        "publication robustness TeX",
    )
    keys = ("contrast", "metric")
    csv_index = table_index(csv_rows, keys, "publication robustness CSV")
    tex_index = table_index(tex_rows, keys, "publication robustness TeX")
    if set(csv_index) != set(expected) or set(tex_index) != set(expected):
        raise ManuscriptReadyError(
            "publication robustness summary inventory differs from reviewed science"
        )
    for key, values in expected.items():
        for rows, tex, label in (
            (csv_index, False, "publication robustness CSV"),
            (tex_index, True, "publication robustness TeX"),
        ):
            row = rows[key]
            for field in ROBUSTNESS_SUMMARY_COLUMNS:
                require_table_value(
                    row, field, values[field], f"{label} {key}", tex=tex
                )


def signed_work_expected_rows(
    downstream: dict[str, object],
) -> dict[str, dict[str, object]]:
    """Build the signed applied-ledger authority directly from bound diagnostics."""

    rows: dict[str, dict[str, object]] = {}
    for case_id in ALL_CASES:
        item = require_dict(
            downstream["analysis"].get(case_id), f"{case_id} completion analysis"
        )
        binding = verify_binding(item.get("diagnostics"), f"{case_id} diagnostics")
        diagnostics = load_json(Path(str(binding["path"])), f"{case_id} diagnostics")
        lf = nested(diagnostics, "windows.steady.lf_history")
        lf = lf if isinstance(lf, dict) else {}
        heat = lf.get("applied_heat_flux_work")
        pressure = lf.get("applied_pressure_work")
        heat = heat if isinstance(heat, dict) else None
        pressure = pressure if isinstance(pressure, dict) else None
        if heat is not None and heat.get("signed") is not True:
            raise ManuscriptReadyError(
                f"{case_id} applied heat-flux work lost signed semantics"
            )
        if pressure is not None and pressure.get("signed") is not True:
            raise ManuscriptReadyError(
                f"{case_id} applied pressure work lost signed semantics"
            )
        if heat is not None:
            parallel = require_finite(
                heat.get("parallel"), f"{case_id} applied heat-flux parallel work"
            )
            perpendicular = require_finite(
                heat.get("perpendicular"),
                f"{case_id} applied heat-flux perpendicular work",
            )
            total = require_finite(
                heat.get("total"), f"{case_id} applied heat-flux total work"
            )
            require_close(
                total,
                parallel + perpendicular,
                f"{case_id} applied heat-flux signed total",
            )
        else:
            parallel = perpendicular = total = None
        if pressure is not None:
            pressure_total = require_finite(
                pressure.get("total"), f"{case_id} applied pressure total work"
            )
            anisotropic = require_finite(
                pressure.get("anisotropic"),
                f"{case_id} applied pressure anisotropic work",
            )
        else:
            pressure_total = anisotropic = None
        caps = lf.get("heat_flux_cap_fractions")
        caps = caps if isinstance(caps, dict) else {}
        rows[case_id] = {
            "case_id": case_id,
            "availability": (
                "available" if heat is not None or pressure is not None else "inconclusive"
            ),
            "diagnostics_provenance": "authenticated",
            "applied_heat_flux_availability": (
                "available" if heat is not None else "inconclusive"
            ),
            "applied_pressure_work_availability": (
                "available" if pressure is not None else "inconclusive"
            ),
            "applied_ledgers_signed": (
                True if heat is not None and pressure is not None else None
            ),
            "applied_heat_flux_parallel": parallel,
            "applied_heat_flux_perpendicular": perpendicular,
            "applied_heat_flux_total": total,
            "applied_pressure_work_total": pressure_total,
            "applied_pressure_work_anisotropic": anisotropic,
            "cap_parallel_over_1": caps.get("parallel_over_1"),
            "cap_parallel_over_10": caps.get("parallel_over_10"),
            "cap_perpendicular_over_1": caps.get("perpendicular_over_1"),
            "cap_perpendicular_over_10": caps.get("perpendicular_over_10"),
            "semantics": SIGNED_WORK_SEMANTICS,
        }
    return rows


def validate_signed_work_table(root: Path, downstream: dict[str, object]) -> None:
    expected = signed_work_expected_rows(downstream)
    csv_rows = load_csv_table(
        root / "tables/signed_lf_cap_work_ledger.csv",
        SIGNED_WORK_COLUMNS,
        "publication signed-work CSV",
    )
    tex_rows = load_tex_table(
        root / "tables/signed_lf_cap_work_ledger.tex",
        SIGNED_WORK_COLUMNS,
        "publication signed-work TeX",
    )
    csv_index = table_index(csv_rows, ("case_id",), "publication signed-work CSV")
    tex_index = table_index(tex_rows, ("case_id",), "publication signed-work TeX")
    expected_keys = {(case_id,) for case_id in ALL_CASES}
    if set(csv_index) != expected_keys or set(tex_index) != expected_keys:
        raise ManuscriptReadyError(
            "publication signed-work ledger case inventory differs"
        )
    fields = (
        "availability",
        "diagnostics_provenance",
        "applied_heat_flux_availability",
        "applied_pressure_work_availability",
        "applied_ledgers_signed",
        "applied_heat_flux_parallel",
        "applied_heat_flux_perpendicular",
        "applied_heat_flux_total",
        "applied_pressure_work_total",
        "applied_pressure_work_anisotropic",
        "cap_parallel_over_1",
        "cap_parallel_over_10",
        "cap_perpendicular_over_1",
        "cap_perpendicular_over_10",
        "semantics",
    )
    for case_id, values in expected.items():
        key = (case_id,)
        for rows, tex, label in (
            (csv_index, False, "publication signed-work CSV"),
            (tex_index, True, "publication signed-work TeX"),
        ):
            row = rows[key]
            for field in fields:
                expected_value = values[field]
                if (
                    not tex
                    and isinstance(expected_value, (int, float))
                    and not isinstance(expected_value, bool)
                ):
                    require_table_number(
                        row,
                        field,
                        float(expected_value),
                        f"{label} {case_id}",
                    )
                else:
                    require_table_text(
                        row,
                        field,
                        expected_value,
                        f"{label} {case_id}",
                        tex=tex,
                    )
    try:
        report = (root / "report.md").read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise ManuscriptReadyError("cannot load publication report") from error
    if SIGNED_WORK_REPORT_STATEMENT not in report:
        raise ManuscriptReadyError(
            "publication report loses signed applied-work semantics"
        )


def validate_publication(
    context: dict[str, object],
    downstream: dict[str, object],
    science: dict[str, object],
) -> dict[str, object]:
    """Validate the exact completion-bound four-root publication."""

    workflow = Path(str(downstream["root"]))
    root = canonical_directory(workflow / "publication", "workflow publication")
    manifest, manifest_binding = load_bound_json(root / "manifest.json", "publication")
    if manifest_binding != downstream["publication_manifest"]:
        raise ManuscriptReadyError(
            "downstream completion publication binding differs from final publication"
        )
    reviewed = require_dict(manifest.get("reviewed_science"), "publication science")
    ct = require_dict(manifest.get("direct_ct_audit"), "publication CT audit")
    if (
        manifest.get("schema_version") != 2
        or manifest.get("record_type") != PUBLICATION_RECORD_TYPE
        or manifest.get("analysis_output") != str(context["inventory_output"])
        or manifest.get("evidence_state")
        not in {"complete", "complete integration / partial evidence"}
        or manifest.get("renderer_ingestion_warnings") != []
        or reviewed.get("record_type") != DIRECT_SCIENCE_RECORD_TYPE
        or reviewed.get("result") != science["record"]["result"]
        or reviewed.get("selected_cases") != list(ALL_CASES)
        or reviewed.get("release_authorizing") is not False
        or ct.get("record_type") != CT_RECORD_TYPE
        or ct.get("numerical_result") != downstream["ct_record"]["result"]
        or ct.get("selected_cases") != list(ALL_CASES)
        or ct.get("full_stage_i_coverage") is not True
        or ct.get("release_authorizing") is not False
    ):
        raise ManuscriptReadyError(
            "publication omits complete integration, reviewed science, or CT coverage"
        )

    roots = publication_roots(context, downstream, science)
    expected_tail = [
        str(PUBLICATION_TOOL.resolve()),
        str(Path(str(context["inventory_output"])).resolve()),
        "--output",
        str(root),
    ]
    for path in sorted(roots):
        expected_tail.extend(["--acceptance", str(path)])
    invocation = [
        require_text(value, "publication invocation argument")
        for value in require_list(
            manifest.get("normalized_invocation"), "publication normalized invocation"
        )
    ]
    if len(invocation) < 2 or invocation[1:] != expected_tail:
        raise ManuscriptReadyError(
            "publication was not rendered with the exact acceptance, science, CT, "
            "and composite hyperbolicity roots"
        )

    products = validated_manifest_bindings(
        manifest.get("products"), "publication products", root=root
    )
    relative_products = {
        str(Path(path).relative_to(root)): value for path, value in products.items()
    }
    if set(relative_products) != set(REQUIRED_PUBLICATION_PRODUCTS):
        missing = sorted(set(REQUIRED_PUBLICATION_PRODUCTS) - set(relative_products))
        extra = sorted(set(relative_products) - set(REQUIRED_PUBLICATION_PRODUCTS))
        raise ManuscriptReadyError(
            f"publication product inventory differs: missing={missing}, extra={extra}"
        )

    sources = validated_manifest_bindings(manifest.get("sources"), "publication sources")
    record_paths = source_record_paths(sources.values())
    exact_records = {
        DIRECT_SCIENCE_RECORD_TYPE: {str(science["direct_record"]["path"])},
        SCIENCE_RECORD_TYPE: {str(science["record_binding"]["path"])},
        CT_RECORD_TYPE: {str(downstream["ct_audit"]["path"])},
    }
    for record_type, expected in exact_records.items():
        if record_paths.get(record_type, set()) != expected:
            raise ManuscriptReadyError(
                f"publication {record_type} source inventory differs"
            )

    required: dict[str, dict[str, object]] = {}

    def add_required(value: object, label: str) -> None:
        current = verify_binding(value, label)
        required[str(current["path"])] = current

    add_required(context["inventory_binding"], "publication inventory")
    add_required(downstream["ct_audit"], "publication CT audit")
    add_required(science["record_binding"], "publication corrected science")
    add_required(science["direct_record"], "publication direct science")
    for case_id, record in downstream["hyperbolicity"].items():
        item = require_dict(record, f"{case_id} completion hyperbolicity")
        for key in ("manifest", "result"):
            add_required(item.get(key), f"{case_id} publication hyperbolicity {key}")
    for case_id, record in downstream["analysis"].items():
        item = require_dict(record, f"{case_id} completion analysis")
        add_required(item.get("diagnostics"), f"{case_id} publication diagnostics")
    for section in (
        require_dict(science["record"].get("acceptance"), "science acceptance"),
        require_dict(science["record"].get("science"), "science products"),
    ):
        for value in recursive_file_bindings(section):
            current = verify_binding(value, "reviewed-science publication source")
            if Path(str(current["path"])).suffix == ".json":
                required[str(current["path"])] = current
    missing = sorted(set(required) - set(sources))
    stale = sorted(
        path
        for path, value in required.items()
        if path in sources and sources[path] != value
    )
    if missing or stale:
        raise ManuscriptReadyError(
            f"publication omits or stale-binds required evidence: "
            f"missing={missing}, stale={stale}"
        )
    validate_reviewed_contrast_table(root, science)
    validate_active_passive_summary_table(root, science)
    validate_robustness_summary_table(root, science)
    validate_signed_work_table(root, downstream)
    return {
        "root": str(root),
        "manifest": manifest_binding,
        "acceptance_roots": [str(path) for path in roots],
        "products": dict(sorted(relative_products.items())),
    }


def validate_manuscript_products(context: dict[str, object]) -> dict[str, object]:
    """Validate the corrected report's manuscript-ready retained products."""

    analysis = canonical_directory(
        Path(str(context["inventory_output"])), "composite report"
    )
    products = {
        relative: artifact_binding(analysis / relative)
        for relative in MANUSCRIPT_PRODUCTS
    }
    verify = load_json(analysis / "verify.json", "composite report verification")
    if (
        verify.get("require_complete") is not True
        or verify.get("result") not in {"pass", "warnings"}
        or verify.get("errors") != []
        or set(require_dict(verify.get("cases"), "verify cases")) != set(ALL_CASES)
    ):
        raise ManuscriptReadyError("composite report verification is not complete")
    require_same_binding(
        verify.get("inventory"), context["inventory_binding"], "verify inventory"
    )

    results = load_json(analysis / "manuscript/results.json", "manuscript results")
    if set(require_dict(results.get("case_results"), "manuscript case results")) != set(
        ALL_CASES
    ):
        raise ManuscriptReadyError("manuscript results do not cover exact R02-R17")
    health = require_dict(results.get("campaign_health"), "manuscript campaign health")
    if health.get("partial_or_unavailable_cases") not in (None, []):
        raise ManuscriptReadyError("manuscript report retains partial cases")
    if health.get("structural_error_cases") not in (None, []):
        raise ManuscriptReadyError("manuscript report retains structural errors")

    figure_manifest = load_json(
        analysis / "manuscript/figure_manifest.json", "manuscript figure manifest"
    )
    figures = require_list(figure_manifest.get("figures"), "manuscript figures")
    if not figures:
        raise ManuscriptReadyError("manuscript figure manifest is empty")
    figure_bindings: dict[str, dict[str, object]] = {}
    for value in figures:
        relative = Path(require_text(value, "manuscript figure"))
        path = (analysis / relative).resolve(strict=True)
        if not path_contains(analysis, path):
            raise ManuscriptReadyError(f"manuscript figure escapes report root: {path}")
        figure_bindings[str(relative)] = artifact_binding(path)

    provenance = load_json(
        analysis / "manuscript/provenance_manifest.json", "manuscript provenance"
    )
    inputs = validated_manifest_bindings(provenance.get("inputs"), "manuscript inputs")
    inventory = verify_binding(context["inventory_binding"], "composite inventory")
    if str(inventory["path"]) not in inputs or inputs[str(inventory["path"])] != inventory:
        raise ManuscriptReadyError("manuscript provenance omits the composite inventory")
    if set(require_dict(provenance.get("case_inputs"), "manuscript case inputs")) != set(
        ALL_CASES
    ):
        raise ManuscriptReadyError("manuscript provenance case coverage differs")
    return {
        "root": str(analysis / "manuscript"),
        "products": dict(sorted(products.items())),
        "figures": dict(sorted(figure_bindings.items())),
    }


def validate_marker_layout(
    args: argparse.Namespace,
    context: dict[str, object],
    downstream: dict[str, object],
    science: dict[str, object],
) -> Path:
    path = marker_path(args.marker)
    expected = Path(str(downstream["root"])).resolve(strict=True) / RECORD_NAME
    if path != expected:
        raise ManuscriptReadyError(
            f"manuscript-ready marker must be exact workflow marker: {expected}"
        )
    for root, label in (
        (Path(str(context["inventory_output"])), "composite report"),
        (Path(str(science["acceptance_root"])), "acceptance"),
        (Path(str(science["root"])), "science"),
        (Path(str(downstream["root"])) / "publication", "publication"),
        (Path(str(downstream["root"])) / "ct", "CT"),
    ):
        if path_contains(root, path):
            raise ManuscriptReadyError(
                f"manuscript-ready marker may not be inside {label} inputs"
            )
    return path


def build_marker(args: argparse.Namespace) -> dict[str, object]:
    """Revalidate all retained inputs in supervisor order and build the marker."""

    science_module, context = validate_corrected_context(args)
    science = validate_reviewed_science(args, science_module, context)
    downstream = validate_downstream_completion(context, args.workflow_root)
    publication = validate_publication(context, downstream, science)
    manuscript = validate_manuscript_products(context)
    marker = validate_marker_layout(args, context, downstream, science)
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": "manuscript_ready",
        "marker_path": str(marker),
        "campaign_kind": "corrected-composite",
        "reviewed_science_result": science["record"]["result"],
        "active_passive_intervention_scope": science[
            "active_passive_intervention_scope"
        ],
        "current_science_scope_limitation": science[
            "current_science_scope_limitation"
        ],
        "ct_numerical_result": downstream["ct_record"]["result"],
        "ct_authentication_complete": True,
        "ct_coverage_complete": True,
        "case_classification": expected_science_classification(),
        "campaign_identity": verify_binding(
            context["identity_binding"], "marker campaign identity"
        ),
        "inventory": verify_binding(context["inventory_binding"], "marker inventory"),
        "tools": {
            "adapter": artifact_binding(SCRIPT_PATH),
            "corrected_report": artifact_binding(CORRECTED_REPORT_TOOL),
            "corrected_downstream": artifact_binding(CORRECTED_DOWNSTREAM_TOOL),
            "corrected_science": artifact_binding(CORRECTED_SCIENCE_TOOL),
            "publication": artifact_binding(PUBLICATION_TOOL),
        },
        "reviewed_science": {
            "record": science["record_binding"],
            "direct_record": science["direct_record"],
            "result": science["record"]["result"],
            "active_passive_intervention_scope": science[
                "active_passive_intervention_scope"
            ],
            "current_science_scope_limitation": science[
                "current_science_scope_limitation"
            ],
        },
        "downstream": {
            "root": downstream["root"],
            "pointer": downstream["pointer"],
            "completion_record": downstream["completion_record"],
            "ct_audit": downstream["ct_audit"],
            "publication_manifest": downstream["publication_manifest"],
        },
        "publication": publication,
        "manuscript_report": manuscript,
        "gate_policy": {
            "adapter_mode": "validation-only",
            "legacy_active_evidence_permitted": False,
            "authenticated_legacy_passive_controls": list(PASSIVE_CASES),
            "reviewed_science_results_permitted": sorted(VALID_SCIENCE_RESULTS),
            "ct_authentication_and_coverage_required": True,
            "marker_is_final_initiative_release": False,
            "required_publication_product_count": len(REQUIRED_PUBLICATION_PRODUCTS),
        },
    }


def run_marker(args: argparse.Namespace) -> Path:
    """Validate retained products and exclusively create the manuscript marker."""

    path = refuse_existing_marker(args.marker)
    record = build_marker(args)
    write_new_marker(path, record)
    if load_json(path, "manuscript-ready marker") != record:
        raise ManuscriptReadyError("written manuscript-ready marker differs")
    print(f"wrote corrected/composite manuscript-ready marker: {path}")
    return path


def validate_marker(args: argparse.Namespace) -> Path:
    """Revalidate one retained manuscript-ready marker and every bound product."""

    expected = build_marker(args)
    path = marker_path(args.marker).resolve(strict=True)
    if load_json(path, "manuscript-ready marker") != expected:
        raise ManuscriptReadyError(
            "manuscript-ready marker differs from current retained evidence"
        )
    print(f"validated corrected/composite manuscript-ready marker: {path}")
    return path


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--inventory-sha256", required=True)
    parser.add_argument("--workflow-root", type=Path, default=DEFAULT_WORKFLOW)
    parser.add_argument("--acceptance-output", type=Path, default=DEFAULT_ACCEPTANCE)
    parser.add_argument("--science-output", type=Path, default=DEFAULT_SCIENCE)
    parser.add_argument("--marker", type=Path, default=DEFAULT_MARKER)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--criteria-review", type=Path, default=DEFAULT_CRITERIA_REVIEW)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="validate retained products and write the immutable marker"
    )
    add_common_arguments(run)
    validate = commands.add_parser(
        "validate", help="revalidate the retained manuscript-ready marker"
    )
    add_common_arguments(validate)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "run":
            run_marker(args)
        else:
            validate_marker(args)
    except (
        ManuscriptReadyError,
        OSError,
        ValueError,
        TypeError,
        KeyError,
    ) as error:
        print(f"manuscript-ready gate error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
