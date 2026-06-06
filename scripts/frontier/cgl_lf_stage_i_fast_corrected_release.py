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
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import sys
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
ACTIVE_EVIDENCE_CLASS = "corrected_active_production"
PASSIVE_EVIDENCE_CLASS = "authenticated_legacy_passive_control"
VALID_SCIENCE_RESULTS = {"pass", "fail", "inconclusive"}
VALID_CT_RESULTS = {"pass", "fail", "inconclusive"}

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
    try:
        with path.open("xb") as stream:
            stream.write(stable_json(value))
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise ManuscriptReadyError(
            f"manuscript-ready marker appeared concurrently: {path}"
        ) from error


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
    return {
        "root": str(science_root),
        "acceptance_root": str(acceptance),
        "record": record,
        "record_binding": record_binding,
        "direct_record": direct_binding,
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
    validate_marker_layout(args, context, downstream, science)
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": "manuscript_ready",
        "campaign_kind": "corrected-composite",
        "reviewed_science_result": science["record"]["result"],
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
