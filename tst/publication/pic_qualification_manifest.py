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
import os
from pathlib import Path, PurePosixPath
import re
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_PATH = READINESS_DIR / "schemas" / "validation_manifest.schema.json"
CLAIMS_PATH = READINESS_DIR / "claims_registry.json"
CONTROL_PLANE_DIR = Path(__file__).with_name("frontier_control_plane")
sys.path.insert(0, str(CONTROL_PLANE_DIR))

from control_plane_common import read_stable_regular_file  # noqa: E402
from control_plane_common import require_canonical_path_below  # noqa: E402
from control_plane_common import sha256_bytes  # noqa: E402
from control_plane_common import source_bundle_sha256  # noqa: E402
from control_plane_common import validate_clean_candidate_bundle  # noqa: E402

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


def _source_bundle_sha256(
    source_archive_sha256: str, submodules: list[dict[str, Any]]
) -> str:
    value = {
        "source_archive_sha256": source_archive_sha256,
        "submodules": [
            {
                "path": record["path"],
                "archive_sha256": record["archive_sha256"],
                "git_commit": record["git_commit"],
                "git_tree": record["git_tree"],
            }
            for record in submodules
        ],
    }
    return source_bundle_sha256(
        str(value["source_archive_sha256"]), value["submodules"]
    )


def _artifact_path(artifact_root: Path, raw_path: str, label: str) -> Path:
    if not isinstance(raw_path, str):
        raise ValueError(f"{label} path must be a string")
    pure = PurePosixPath(raw_path)
    if (
        not raw_path
        or raw_path != pure.as_posix()
        or any(part in {"", ".", ".."} for part in pure.parts)
    ):
        raise ValueError(f"{label} path is not canonical: {raw_path!r}")
    path = Path(raw_path)
    if not path.is_absolute():
        path = artifact_root / path
    path = require_canonical_path_below(path, artifact_root)
    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")
    return path


def _artifact_bytes(artifact_root: Path, raw_path: str, label: str) -> bytes:
    return read_stable_regular_file(_artifact_path(artifact_root, raw_path, label))


def _load_object_bytes(data: bytes, *, label: str) -> dict[str, Any]:
    value = json.loads(data.decode("utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return value


def _project_submodules(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {
            "path": str(record["path"]),
            "archive_sha256": str(record["archive_sha256"]),
            "git_commit": str(record["git_commit"]),
            "git_tree": str(record["git_tree"]),
        }
        for record in records
    ]


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
    submodules = manifest["git"]["submodules"]
    paths = [record["path"] for record in submodules]
    if paths != sorted(set(paths)):
        raise ValueError("qualification manifest submodules must be unique and sorted")
    expected_status = "clean_pinned_archived" if submodules else "absent"
    if manifest["git"]["submodule_status"] != expected_status:
        raise ValueError("qualification manifest submodule status does not match records")
    if manifest["git"]["source_bundle_sha256"] != _source_bundle_sha256(
        manifest["git"]["source_archive"]["sha256"], submodules
    ):
        raise ValueError("qualification manifest source-bundle checksum mismatch")

    if not verify_files:
        return

    artifact_root = Path(manifest["resources"]["artifact_root"])
    if not artifact_root.is_absolute():
        raise ValueError("artifact root must be absolute")
    lexical_root = Path(os.path.abspath(artifact_root))
    if lexical_root != artifact_root or lexical_root.resolve() != lexical_root:
        raise ValueError("artifact root must use its canonical spelling")
    artifact_root = lexical_root
    if not artifact_root.is_dir():
        raise ValueError(f"artifact root is not a directory: {artifact_root}")

    source_archive = _artifact_bytes(
        artifact_root, manifest["git"]["source_archive"]["path"], "source archive"
    )
    if sha256_bytes(source_archive) != manifest["git"]["source_archive"]["sha256"]:
        raise ValueError("source archive checksum mismatch")
    submodule_archives = []
    for index, submodule in enumerate(submodules):
        archive = _artifact_bytes(
            artifact_root, submodule["archive_path"], f"submodule archive {index}"
        )
        if sha256_bytes(archive) != submodule["archive_sha256"]:
            raise ValueError(f"submodule archive {index} checksum mismatch")
        submodule_archives.append(archive)
    candidate_bytes = _artifact_bytes(
        artifact_root,
        manifest["git"]["clean_candidate_manifest"]["path"],
        "clean-candidate manifest",
    )
    if sha256_bytes(candidate_bytes) != manifest["git"]["clean_candidate_manifest"]["sha256"]:
        raise ValueError("clean-candidate manifest checksum mismatch")

    executable = _artifact_bytes(
        artifact_root, manifest["executable"]["path"], "executable"
    )
    executable_sha256 = sha256_bytes(executable)
    if executable_sha256 != manifest["executable"]["sha256"]:
        raise ValueError("executable checksum mismatch")
    candidate = _load_object_bytes(candidate_bytes, label="clean-candidate manifest")
    candidate_submodules = validate_clean_candidate_bundle(
        candidate,
        source_archive=source_archive,
        submodule_archives=submodule_archives,
        executable_sha256=executable_sha256,
    )
    candidate_source = candidate["source"]
    if not isinstance(candidate_source, dict):
        raise ValueError("clean-candidate source attestation is missing")
    if (
        manifest["git"]["commit"] != candidate_source["git_commit"]
        or manifest["git"]["tree"] != candidate_source["git_tree"]
        or manifest["git"]["source_archive"]["sha256"]
        != candidate_source["archive_sha256"]
        or manifest["git"]["source_bundle_sha256"]
        != candidate_source["source_bundle_sha256"]
        or manifest["git"]["submodule_status"]
        != candidate_source["submodule_status"]
        or _project_submodules(submodules) != candidate_submodules
    ):
        raise ValueError("qualification manifest source projection differs from clean candidate")

    for field in ("cmake_cache", "modules", "environment_allowlist"):
        _artifact_bytes(
            artifact_root,
            manifest["executable"][field],
            f"executable provenance field {field}",
        )

    for index, artifact in enumerate(manifest["artifacts"]):
        data = _artifact_bytes(
            artifact_root, artifact["path"], f"artifact {index}"
        )
        if sha256_bytes(data) != artifact["sha256"]:
            raise ValueError(f"artifact {index} checksum mismatch")


def freeze_qualification_manifest(input_path: Path, output_path: Path) -> None:
    """Validate a prepared manifest and write a canonical immutable candidate."""
    manifest = _load_object(input_path)
    validate_qualification_manifest(manifest)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    output_path.chmod(0o444)


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
