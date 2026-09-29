#!/usr/bin/env python3
"""Write and validate an isolated corrected-production campaign identity.

The retained ``campaign-identity.json`` binds the corrected campaign root to
the exact source revision and executable, EOS, matrix, and audit SHA-256
identities.  Selected run paths are stored relative to that root and are
rechecked on every validation so legacy and corrected outputs cannot be mixed
silently by downstream workflows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from typing import Sequence


SCHEMA = "athenak-cgl-corrected-campaign-identity"
SCHEMA_VERSION = 1
CAMPAIGN_KIND = "corrected-production"
IDENTITY_NAME = "campaign-identity.json"
ARTIFACT_NAMES = ("audit", "eos", "executable", "matrix")
IDENTITY_KEYS = {
    "artifacts",
    "campaign_id",
    "campaign_kind",
    "campaign_root",
    "schema",
    "schema_version",
    "selected_run_paths",
    "source_revision",
}
ARTIFACT_KEYS = {"path", "sha256"}
CAMPAIGN_ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
SOURCE_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


class IdentityError(RuntimeError):
    """Raised when corrected-campaign identity validation fails."""


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one regular file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_exact_keys(value: object, expected: set[str], label: str) -> dict:
    """Require a JSON object with exactly the declared keys."""

    if not isinstance(value, dict):
        raise IdentityError(f"{label} must be an object")
    observed = set(value)
    if observed != expected:
        raise IdentityError(
            f"{label} keys differ: expected {sorted(expected)}, "
            f"observed {sorted(observed)}"
        )
    return value


def require_campaign_id(value: object) -> str:
    """Require a corrected-production campaign identifier."""

    if not isinstance(value, str) or CAMPAIGN_ID_PATTERN.fullmatch(value) is None:
        raise IdentityError(
            "campaign_id must contain only letters, numbers, '.', '_', and '-'"
        )
    if "corrected" not in value.lower():
        raise IdentityError("campaign_id must identify the campaign as corrected")
    return value


def require_source_revision(value: object) -> str:
    """Require a full lowercase Git SHA-1 source revision."""

    if not isinstance(value, str) or SOURCE_REVISION_PATTERN.fullmatch(value) is None:
        raise IdentityError("source_revision must be a full lowercase 40-hex revision")
    return value


def require_sha256(value: object, label: str) -> str:
    """Require a lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise IdentityError(f"{label} must be a lowercase 64-hex SHA-256")
    return value


def canonical_campaign_root(path: Path, *, create: bool) -> Path:
    """Return an existing canonical campaign root, creating it only for write."""

    expanded = path.expanduser()
    if create:
        expanded.mkdir(parents=True, exist_ok=True)
    try:
        root = expanded.resolve(strict=True)
    except OSError as error:
        raise IdentityError(f"campaign root does not resolve: {expanded}") from error
    if not root.is_dir():
        raise IdentityError(f"campaign root is not a directory: {root}")
    return root


def normalize_run_path(root: Path, value: str | Path, label: str) -> str:
    """Return one selected run path relative to root, rejecting any escape."""

    if isinstance(value, str) and not value:
        raise IdentityError(f"{label} must be a nonempty path")
    path = Path(value).expanduser()
    candidate = path if path.is_absolute() else root / path
    try:
        resolved = candidate.resolve(strict=False)
    except OSError as error:
        raise IdentityError(f"{label} does not resolve: {candidate}") from error
    try:
        relative = resolved.relative_to(root)
    except ValueError as error:
        raise IdentityError(
            f"{label} escapes corrected campaign root {root}: {candidate}"
        ) from error
    if relative == Path("."):
        raise IdentityError(f"{label} must be a child of corrected campaign root {root}")
    return relative.as_posix()


def selected_run_paths(root: Path, values: Sequence[str | Path]) -> list[str]:
    """Normalize, deduplicate, and sort the selected corrected run paths."""

    if not values:
        raise IdentityError("at least one selected run path is required")
    normalized = {
        normalize_run_path(root, value, f"selected run path {index}")
        for index, value in enumerate(values, start=1)
    }
    return sorted(normalized)


def artifact_binding(path: Path, label: str) -> dict[str, str]:
    """Return the canonical path and SHA-256 for one required regular file."""

    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise IdentityError(f"{label} does not resolve: {path}") from error
    if not resolved.is_file():
        raise IdentityError(f"{label} is not a regular file: {resolved}")
    return {"path": str(resolved), "sha256": sha256(resolved)}


def stable_json(value: object) -> bytes:
    """Serialize retained identity JSON deterministically."""

    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def build_identity(
    campaign_id: str,
    campaign_root: Path,
    source_revision: str,
    artifacts: dict[str, Path],
    run_paths: Sequence[str | Path],
) -> dict[str, object]:
    """Build a validated corrected-production identity record."""

    root = canonical_campaign_root(campaign_root, create=True)
    if set(artifacts) != set(ARTIFACT_NAMES):
        raise IdentityError(
            f"artifacts must be exactly {list(ARTIFACT_NAMES)}"
        )
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "campaign_kind": CAMPAIGN_KIND,
        "campaign_id": require_campaign_id(campaign_id),
        "campaign_root": str(root),
        "source_revision": require_source_revision(source_revision),
        "artifacts": {
            name: artifact_binding(artifacts[name], f"{name} artifact")
            for name in ARTIFACT_NAMES
        },
        "selected_run_paths": selected_run_paths(root, run_paths),
    }


def write_identity(
    campaign_id: str,
    campaign_root: Path,
    source_revision: str,
    artifacts: dict[str, Path],
    run_paths: Sequence[str | Path],
) -> Path:
    """Create the canonical identity, allowing only an identical existing file."""

    identity = build_identity(
        campaign_id, campaign_root, source_revision, artifacts, run_paths
    )
    root = Path(str(identity["campaign_root"]))
    path = root / IDENTITY_NAME
    payload = stable_json(identity)
    if path.exists():
        if not path.is_file() or path.read_bytes() != payload:
            raise IdentityError(f"refusing to replace differing identity: {path}")
        validate_identity(path)
        return path
    try:
        with path.open("xb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise IdentityError(f"identity appeared concurrently: {path}") from error
    validate_identity(path)
    return path


def load_identity(path: Path) -> tuple[dict, Path]:
    """Load one canonical identity JSON object."""

    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise IdentityError(f"identity does not resolve: {path}") from error
    if not resolved.is_file():
        raise IdentityError(f"identity is not a regular file: {resolved}")
    try:
        value = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise IdentityError(f"identity is not valid UTF-8 JSON: {resolved}") from error
    return require_exact_keys(value, IDENTITY_KEYS, "identity"), resolved


def validate_artifacts(value: object) -> dict[str, dict[str, str]]:
    """Authenticate every artifact binding in an identity."""

    artifacts = require_exact_keys(value, set(ARTIFACT_NAMES), "identity artifacts")
    verified: dict[str, dict[str, str]] = {}
    for name in ARTIFACT_NAMES:
        binding = require_exact_keys(
            artifacts[name], ARTIFACT_KEYS, f"{name} artifact binding"
        )
        declared_path = binding["path"]
        if not isinstance(declared_path, str) or not Path(declared_path).is_absolute():
            raise IdentityError(f"{name} artifact path must be absolute")
        observed = artifact_binding(Path(declared_path), f"{name} artifact")
        if observed["path"] != declared_path:
            raise IdentityError(f"{name} artifact path is not canonical: {declared_path}")
        if observed["sha256"] != require_sha256(
            binding["sha256"], f"{name} artifact SHA-256"
        ):
            raise IdentityError(f"{name} artifact SHA-256 differs from identity")
        verified[name] = observed
    return verified


def validate_identity(
    identity_path: Path,
    requested_run_paths: Sequence[str | Path] = (),
) -> dict[str, object]:
    """Validate the identity and optional downstream-selected run paths."""

    identity, path = load_identity(identity_path)
    if (
        not isinstance(identity["schema"], str)
        or identity["schema"] != SCHEMA
        or type(identity["schema_version"]) is not int
        or identity["schema_version"] != SCHEMA_VERSION
    ):
        raise IdentityError("identity schema or schema_version is unsupported")
    if (
        not isinstance(identity["campaign_kind"], str)
        or identity["campaign_kind"] != CAMPAIGN_KIND
    ):
        raise IdentityError(f"identity campaign_kind must be {CAMPAIGN_KIND}")
    require_campaign_id(identity["campaign_id"])
    require_source_revision(identity["source_revision"])

    root_value = identity["campaign_root"]
    if not isinstance(root_value, str) or not Path(root_value).is_absolute():
        raise IdentityError("campaign_root must be an absolute path")
    root = canonical_campaign_root(Path(root_value), create=False)
    if str(root) != root_value:
        raise IdentityError(f"campaign_root is not canonical: {root_value}")
    if path != root / IDENTITY_NAME:
        raise IdentityError(f"identity must be located at {root / IDENTITY_NAME}")

    raw_selected = identity["selected_run_paths"]
    if not isinstance(raw_selected, list) or not raw_selected:
        raise IdentityError("selected_run_paths must be a nonempty list")
    if not all(isinstance(value, str) for value in raw_selected):
        raise IdentityError("selected_run_paths entries must be strings")
    normalized = selected_run_paths(root, raw_selected)
    if normalized != raw_selected:
        raise IdentityError("selected_run_paths must be sorted unique relative paths")

    requested = (
        selected_run_paths(root, requested_run_paths) if requested_run_paths else []
    )
    unselected = sorted(set(requested) - set(normalized))
    if unselected:
        raise IdentityError(
            f"requested run paths are not selected by identity: {unselected}"
        )

    artifacts = validate_artifacts(identity["artifacts"])
    return {
        "result": "pass",
        "identity": {
            "path": str(path),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        },
        "campaign_id": identity["campaign_id"],
        "campaign_kind": identity["campaign_kind"],
        "campaign_root": str(root),
        "source_revision": identity["source_revision"],
        "artifacts": artifacts,
        "selected_run_paths": normalized,
        "requested_run_paths": requested,
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    write = subparsers.add_parser("write", help="write an isolated campaign identity")
    write.add_argument("--campaign-id", required=True)
    write.add_argument("--campaign-root", type=Path, required=True)
    write.add_argument("--source-revision", required=True)
    write.add_argument("--executable", type=Path, required=True)
    write.add_argument("--eos", type=Path, required=True)
    write.add_argument("--matrix", type=Path, required=True)
    write.add_argument("--audit", type=Path, required=True)
    write.add_argument("--run-path", action="append", required=True)

    validate = subparsers.add_parser("validate", help="validate a campaign identity")
    validate.add_argument("--identity", type=Path, required=True)
    validate.add_argument("--run-path", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the writer or validator and emit a deterministic validation result."""

    args = build_parser().parse_args(argv)
    try:
        if args.command == "write":
            identity_path = write_identity(
                args.campaign_id,
                args.campaign_root,
                args.source_revision,
                {
                    "audit": args.audit,
                    "eos": args.eos,
                    "executable": args.executable,
                    "matrix": args.matrix,
                },
                args.run_path,
            )
            result = validate_identity(identity_path, args.run_path)
        else:
            result = validate_identity(args.identity, args.run_path)
    except IdentityError as error:
        print(f"corrected campaign identity error: {error}", file=sys.stderr)
        return 2
    json.dump(result, sys.stdout, indent=2, sort_keys=True, allow_nan=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
