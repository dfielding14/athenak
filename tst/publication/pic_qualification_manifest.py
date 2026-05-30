#!/usr/bin/env python3
"""Validate and freeze a PIC qualification manifest.

This gate is intentionally separate from the engineering-proxy publication
helpers.  It does not run a campaign or promote a claim.  It verifies that an
already prepared manifest is structurally complete, points at a clean source
candidate, references known claims, and binds existing artifacts by checksum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_PATH = READINESS_DIR / "schemas" / "validation_manifest.schema.json"
CLAIMS_PATH = READINESS_DIR / "claims_registry.json"

QUALIFYING_EVIDENCE_CLASSES = {
    "physics_validation",
    "sun_bai_2023_reproduction",
    "athenak_production_mode",
    "cross_code_comparison",
    "scoped_state_of_the_art",
}


def _schema_matches(value: object, schema: dict[str, Any]) -> bool:
    try:
        validate_schema(value, schema)
    except ValueError:
        return False
    return True


def validate_schema(
    value: object, schema: dict[str, Any], path: str = "$"
) -> None:
    """Validate the JSON-schema keywords used by PIC readiness records."""
    if "const" in schema and value != schema["const"]:
        raise ValueError(f"{path} does not match const")
    if "enum" in schema and value not in schema["enum"]:
        raise ValueError(f"{path} is not in enum")

    expected_type = schema.get("type")
    type_matches = {
        "object": isinstance(value, dict),
        "array": isinstance(value, list),
        "string": isinstance(value, str),
        "number": isinstance(value, (int, float)) and not isinstance(value, bool),
    }
    if expected_type is not None and not type_matches[str(expected_type)]:
        raise ValueError(f"{path} is not a {expected_type}")

    if isinstance(value, dict):
        required = schema.get("required", [])
        for key in required:
            if key not in value:
                raise ValueError(f"{path} is missing {key}")
        min_properties = schema.get("minProperties")
        if min_properties is not None and len(value) < min_properties:
            raise ValueError(f"{path} has too few properties")
        properties = schema.get("properties", {})
        if schema.get("additionalProperties") is False:
            extra = set(value) - set(properties)
            if extra:
                raise ValueError(f"{path} has unexpected properties: {extra}")
        for key, child_schema in properties.items():
            if key in value:
                validate_schema(value[key], child_schema, f"{path}.{key}")

    if isinstance(value, list):
        min_items = schema.get("minItems")
        if min_items is not None and len(value) < min_items:
            raise ValueError(f"{path} has too few items")
        item_schema = schema.get("items")
        if item_schema is not None:
            for index, item in enumerate(value):
                validate_schema(item, item_schema, f"{path}[{index}]")

    if isinstance(value, str):
        min_length = schema.get("minLength")
        if min_length is not None and len(value) < min_length:
            raise ValueError(f"{path} is too short")
        pattern = schema.get("pattern")
        if pattern is not None and re.search(str(pattern), value) is None:
            raise ValueError(f"{path} does not match pattern")

    minimum = schema.get("minimum")
    if minimum is not None and value < minimum:
        raise ValueError(f"{path} is below minimum")

    rejected_schema = schema.get("not")
    if rejected_schema is not None and _schema_matches(value, rejected_schema):
        raise ValueError(f"{path} matches rejected schema")


def _load_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact_path(artifact_root: Path, raw_path: str, label: str) -> Path:
    path = Path(raw_path)
    if not path.is_absolute():
        path = artifact_root / path
    path = path.resolve()
    try:
        path.relative_to(artifact_root)
    except ValueError as exc:
        raise ValueError(f"{label} escapes artifact root: {raw_path}") from exc
    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")
    return path


def validate_qualification_manifest(
    manifest: dict[str, Any], *, verify_files: bool = True
) -> None:
    """Fail closed unless a manifest is ready for external scientific review."""
    validate_schema(manifest, _load_object(SCHEMA_PATH))

    evidence_class = manifest["evidence_class"]
    if evidence_class not in QUALIFYING_EVIDENCE_CLASSES:
        raise ValueError(
            f"evidence class is not qualifying evidence: {evidence_class}"
        )

    claim_registry = _load_object(CLAIMS_PATH)
    known_claim_ids = {
        claim["claim_id"] for claim in claim_registry["claims"]
    }
    unknown_claim_ids = set(manifest["claim_ids"]) - known_claim_ids
    if unknown_claim_ids:
        raise ValueError(f"manifest references unknown claims: {unknown_claim_ids}")

    if manifest["git"]["status"]:
        raise ValueError("qualification manifest requires a clean git status")

    if not verify_files:
        return

    artifact_root = Path(manifest["resources"]["artifact_root"]).resolve()
    if not artifact_root.is_dir():
        raise ValueError(f"artifact root is not a directory: {artifact_root}")

    executable = _artifact_path(
        artifact_root, manifest["executable"]["path"], "executable"
    )
    if _sha256(executable) != manifest["executable"]["sha256"]:
        raise ValueError("executable checksum mismatch")

    for field in ("cmake_cache", "modules", "environment_allowlist"):
        _artifact_path(
            artifact_root,
            manifest["executable"][field],
            f"executable provenance field {field}",
        )

    for index, artifact in enumerate(manifest["artifacts"]):
        path = _artifact_path(
            artifact_root, artifact["path"], f"artifact {index}"
        )
        if _sha256(path) != artifact["sha256"]:
            raise ValueError(f"artifact {index} checksum mismatch")


def freeze_qualification_manifest(input_path: Path, output_path: Path) -> None:
    """Validate a prepared manifest and write a canonical immutable candidate."""
    manifest = _load_object(input_path)
    validate_qualification_manifest(manifest)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--validate",
        type=Path,
        metavar="MANIFEST",
        help="validate an existing qualification manifest",
    )
    group.add_argument(
        "--input",
        type=Path,
        metavar="MANIFEST",
        help="validate and canonically freeze a prepared manifest",
    )
    parser.add_argument(
        "--output",
        type=Path,
        metavar="MANIFEST",
        help="exclusive output path required with --input",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.validate is not None:
        if args.output is not None:
            raise ValueError("--output is only valid with --input")
        validate_qualification_manifest(_load_object(args.validate))
    else:
        if args.output is None:
            raise ValueError("--output is required with --input")
        freeze_qualification_manifest(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
