#!/usr/bin/env python3
"""Fail-closed read-only preflight for final direct-fast Stage I analysis.

The reviewed scientific-acceptance policy already authenticates the exact
criteria, review, Stage I manifest, and top-level reference-archive manifests.
This preflight reuses that path, then verifies every hash-declared member of
both reference archives and every admitted Stage I reference product before
acceptance, CT, or publication adapters are run.

The command only reads inputs and emits deterministic JSON or Markdown to
stdout.  It never writes to simulations, analysis outputs, or reference trees.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
ACCEPTANCE_UTILITY = (
    REPOSITORY_ROOT / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
)
DEFAULT_CRITERIA = REPOSITORY_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
DEFAULT_CRITERIA_REVIEW = REPOSITORY_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)
ARCHIVE_SOURCE_KEYS = (
    "reference_archive_manifest",
    "verified_reference_archive_manifest",
)


class PreflightError(RuntimeError):
    """Raised when required final-analysis reference inputs are not ready."""


def load_acceptance_module() -> Any:
    """Import the existing reviewed policy and reference verification path."""

    name = "_cgl_lf_stage_i_scientific_acceptance_for_fast_preflight"
    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, ACCEPTANCE_UTILITY)
    if spec is None or spec.loader is None:
        raise PreflightError(
            f"cannot import reviewed acceptance utility: {ACCEPTANCE_UTILITY}"
        )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def archive_member_path(root: Path, value: object, label: str) -> Path:
    """Resolve one required relative archive member without permitting escape."""

    if not isinstance(value, str) or not value:
        raise PreflightError(f"{label} path must be a nonempty string")
    relative = Path(value)
    if relative.is_absolute():
        raise PreflightError(f"{label} path must be relative to its archive")
    try:
        resolved = (root / relative).resolve(strict=True)
    except OSError as error:
        raise PreflightError(f"{label} does not resolve: {root / relative}") from error
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise PreflightError(f"{label} escapes its bound archive: {relative}") from error
    return resolved


def verify_archive_member(
    acceptance: Any,
    root: Path,
    relative_path: object,
    expected_sha256: object,
    label: str,
) -> dict[str, object]:
    """Verify one hash-bound archive file with the existing stable-file reader."""

    path = archive_member_path(root, relative_path, label)
    expected = acceptance.require_sha256(expected_sha256, f"{label} sha256")
    observed = acceptance.regular_file_binding(path, label)
    if observed["sha256"] != expected:
        raise PreflightError(f"{label} SHA-256 differs from its declaration")
    return {
        "relative_path": str(relative_path),
        "path": observed["path"],
        "size_bytes": observed["size_bytes"],
        "sha256": observed["sha256"],
    }


def verify_archive_tree(
    acceptance: Any, source_key: str, binding: object
) -> dict[str, object]:
    """Verify one bound archive manifest and every file it hash-declares."""

    bound = acceptance.require_dict(binding, f"{source_key} binding")
    manifest_path = Path(str(bound.get("path")))
    manifest, observed_manifest = acceptance.load_json(
        manifest_path, f"{source_key} manifest"
    )
    if observed_manifest != bound:
        raise PreflightError(f"{source_key} changed after policy validation")
    if manifest.get("schema_version") != 1:
        raise PreflightError(f"{source_key} schema_version must be 1")

    root = manifest_path.parent.resolve(strict=True)
    if not root.is_dir():
        raise PreflightError(f"{source_key} archive root is not a directory: {root}")

    archive = acceptance.require_dict(
        manifest.get("archive"), f"{source_key} archive declaration"
    )
    files = [
        {
            "role": "source_archive",
            **verify_archive_member(
                acceptance,
                root,
                archive.get("path"),
                archive.get("sha256"),
                f"{source_key} source archive",
            ),
        }
    ]
    extracted = acceptance.require_dict(
        manifest.get("extracted_file_sha256"),
        f"{source_key} extracted-file declarations",
    )
    if not extracted:
        raise PreflightError(f"{source_key} extracted-file declarations are empty")
    for relative_path in sorted(extracted):
        files.append(
            {
                "role": "extracted_file",
                **verify_archive_member(
                    acceptance,
                    root,
                    relative_path,
                    extracted[relative_path],
                    f"{source_key} extracted file {relative_path}",
                ),
            }
        )
    source_tex = manifest.get("source_tex")
    if source_tex not in extracted:
        raise PreflightError(f"{source_key} source_tex is not hash-declared")
    return {
        "source_key": source_key,
        "root": str(root),
        "manifest": observed_manifest,
        "verified_file_count": len(files),
        "files": files,
    }


def required_reference_product_ids(
    acceptance: Any, policy: dict[str, object]
) -> list[str]:
    """Return the exact products required by criteria-admitted comparison panels."""

    criteria = acceptance.require_dict(policy.get("criteria"), "criteria")
    admitted = {
        str(acceptance.require_dict(value, "criteria comparison panel").get("id"))
        for value in acceptance.require_list(
            criteria.get("comparison_panels"), "criteria comparison panels"
        )
    }
    manifest = acceptance.require_dict(policy.get("manifest"), "Stage I manifest")
    panel_status = acceptance.require_dict(
        manifest.get("panel_status"), "Stage I manifest panel_status"
    )
    required: set[str] = set()
    seen: set[str] = set()
    for value in acceptance.require_list(
        panel_status.get("panels"), "Stage I manifest panels"
    ):
        panel = acceptance.require_dict(value, "Stage I manifest panel")
        panel_id = str(panel.get("id"))
        if panel_id not in admitted:
            continue
        if panel.get("disposition") != "comparison":
            raise PreflightError(
                f"criteria-admitted panel {panel_id} is not a comparison panel"
            )
        seen.add(panel_id)
        for product_id in acceptance.require_list(
            panel.get("reference_products"), f"panel {panel_id} reference products"
        ):
            if not isinstance(product_id, str) or not product_id:
                raise PreflightError(
                    f"panel {panel_id} has an invalid reference product id"
                )
            required.add(product_id)
    if seen != admitted:
        missing = ", ".join(sorted(admitted - seen))
        raise PreflightError(
            f"criteria-admitted panels are absent from the manifest: {missing}"
        )
    return sorted(required)


def verify_stage_i_references(
    acceptance: Any, policy: dict[str, object]
) -> dict[str, object]:
    """Verify Stage I reference manifests and every admitted reference CSV."""

    manifest = acceptance.require_dict(policy.get("manifest"), "Stage I manifest")
    panel_status = acceptance.require_dict(
        manifest.get("panel_status"), "Stage I manifest panel_status"
    )
    reference_manifests = acceptance.require_dict(
        panel_status.get("reference_manifests"), "reference_manifests"
    )
    product_bindings = acceptance.require_dict(
        panel_status.get("reference_product_bindings"), "reference_product_bindings"
    )
    verified_sources = acceptance.require_dict(
        policy.get("verified_sources"), "verified_sources"
    )
    archive_binding = acceptance.require_dict(
        verified_sources.get("reference_archive_manifest"),
        "reference archive manifest binding",
    )
    archive_root = Path(str(archive_binding.get("path"))).parent.resolve(strict=True)

    manifests: list[dict[str, object]] = []
    for manifest_id in sorted(reference_manifests):
        declaration = acceptance.require_dict(
            reference_manifests[manifest_id], f"reference manifest {manifest_id}"
        )
        manifests.append(
            {
                "manifest_id": manifest_id,
                **verify_archive_member(
                    acceptance,
                    archive_root,
                    declaration.get("path"),
                    declaration.get("sha256"),
                    f"reference manifest {manifest_id}",
                ),
            }
        )

    products: list[dict[str, object]] = []
    for product_id in required_reference_product_ids(acceptance, policy):
        binding = acceptance.require_dict(
            product_bindings.get(product_id), f"reference product {product_id} binding"
        )
        source = acceptance.verified_reference_product(policy, product_id, binding)
        data_path = Path(str(source["data_path"])).resolve(strict=True)
        try:
            data_path.relative_to(archive_root)
        except ValueError as error:
            raise PreflightError(
                f"reference product {product_id} escapes the bound archive"
            ) from error
        products.append(
            {
                "product_id": product_id,
                "reference_manifest": binding.get("reference_manifest"),
                "data_path": str(data_path),
                "data_sha256": source["data_sha256"],
                "manifest_sha256": source["manifest_sha256"],
            }
        )
    return {
        "archive_root": str(archive_root),
        "declared_manifest_count": len(reference_manifests),
        "verified_manifests": manifests,
        "declared_product_binding_count": len(product_bindings),
        "required_product_count": len(products),
        "verified_products": products,
    }


def build_preflight(
    criteria_path: Path, review_path: Path, *, acceptance: Any | None = None
) -> dict[str, object]:
    """Build deterministic pass evidence or raise before final analysis starts."""

    module = acceptance if acceptance is not None else load_acceptance_module()
    policy = module.load_validated_policy(criteria_path, review_path)
    verified_sources = module.require_dict(
        policy.get("verified_sources"), "verified_sources"
    )
    archives = [
        verify_archive_tree(module, key, verified_sources.get(key))
        for key in ARCHIVE_SOURCE_KEYS
    ]
    references = verify_stage_i_references(module, policy)
    return {
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-final-analysis-preflight",
        "result": "pass",
        "read_only": True,
        "fail_closed": True,
        "criteria": policy["criteria_binding"],
        "criteria_review": policy["review_binding"],
        "acceptance_utility": verified_sources["acceptance_utility"],
        "stage_i_manifest": verified_sources["stage_i_manifest"],
        "preflight_utility": module.regular_file_binding(
            Path(__file__).resolve(), "final-analysis preflight utility"
        ),
        "summary": {
            "verified_archive_count": len(archives),
            "verified_archive_file_count": sum(
                int(archive["verified_file_count"]) for archive in archives
            ),
            "verified_reference_manifest_count": len(
                references["verified_manifests"]
            ),
            "verified_reference_product_count": len(
                references["verified_products"]
            ),
        },
        "archives": archives,
        "stage_i_references": references,
    }


def failure_record(
    criteria_path: Path, review_path: Path, error: Exception
) -> dict[str, object]:
    """Return deterministic fail-closed command output."""

    return {
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-final-analysis-preflight",
        "result": "fail",
        "read_only": True,
        "fail_closed": True,
        "criteria_path": str(criteria_path.expanduser().absolute()),
        "criteria_review_path": str(review_path.expanduser().absolute()),
        "error": str(error),
    }


def render_markdown(record: dict[str, object]) -> str:
    """Render concise deterministic preflight Markdown."""

    lines = [
        "# Stage I Final Analysis Preflight",
        "",
        f"Result: **{record['result']}**.",
        "",
        "This fail-closed preflight is read-only and does not modify simulations.",
    ]
    if record["result"] != "pass":
        lines.extend(["", f"Error: `{record.get('error', 'unknown error')}`"])
        return "\n".join(lines) + "\n"
    summary = record["summary"]
    lines.extend(
        [
            "",
            "## Verified Inputs",
            "",
            f"- Reference archives: {summary['verified_archive_count']}",
            f"- Hash-bound archive files: {summary['verified_archive_file_count']}",
            "- Stage I reference manifests: "
            f"{summary['verified_reference_manifest_count']}",
            "- Required Stage I reference products: "
            f"{summary['verified_reference_product_count']}",
            "",
            "## Archives",
            "",
            "| Source key | Root | Verified files |",
            "|---|---|---:|",
        ]
    )
    for archive in record["archives"]:
        lines.append(
            f"| {archive['source_key']} | `{archive['root']}` | "
            f"{archive['verified_file_count']} |"
        )
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--criteria-review", type=Path, default=DEFAULT_CRITERIA_REVIEW)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        record = build_preflight(args.criteria, args.criteria_review)
    except Exception as error:
        record = failure_record(args.criteria, args.criteria_review, error)
    if args.format == "markdown":
        sys.stdout.write(render_markdown(record))
    else:
        json.dump(record, sys.stdout, indent=2, sort_keys=True, allow_nan=False)
        sys.stdout.write("\n")
    return 0 if record["result"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
