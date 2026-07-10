#!/usr/bin/env python3
"""Create and verify content-addressed Frontier release artifacts."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


GIB = 1024**3
SOURCE_PATTERNS = (
    "src/srcterms/cooling_tables.hpp",
    "tools/gotham_pdf_rebuild/CMakeLists.txt",
    "tools/gotham_pdf_rebuild/gotham_pdf_rebuild.cpp",
    "tools/gotham_pdf_rebuild/plot_rebuilt_pdfs.py",
    "tools/gotham_pdf_rebuild/validate_real_8pc.py",
    "tools/gotham_pdf_rebuild/tests/*.py",
    "tools/gotham_pdf_rebuild/tests/pytest.ini",
    "tools/gotham_pdf_rebuild/jobs/*.py",
    "tools/gotham_pdf_rebuild/jobs/*.sh",
    "tools/gotham_pdf_rebuild/jobs/*.sbatch",
    "tools/gotham_pdf_rebuild/jobs/*.json",
    "tools/gotham_pdf_rebuild/jobs/manifests/*.tsv",
)


class GateError(RuntimeError):
    """Raised when a release artifact cannot be trusted."""


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GateError(f"unable to read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GateError(f"expected a JSON object: {path}")
    return value


def require_file(path: Path, description: str) -> Path:
    path = path.resolve()
    if not path.is_file():
        raise GateError(f"{description} not found: {path}")
    return path


def write_content_addressed(
    directory: Path, prefix: str, payload: dict[str, Any]
) -> Path:
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    encoded = canonical_json(payload)
    digest = sha256_bytes(encoded)
    safe_prefix = re.sub(r"[^A-Za-z0-9_.-]+", "-", prefix)
    path = directory / f"{safe_prefix}.{digest}.json"
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if path.exists() or sidecar.exists():
        if path.is_file() and sha256_file(path) == digest:
            return path
        raise GateError(f"refusing conflicting release artifact: {path}")
    partial = path.with_name(path.name + ".partial")
    partial.write_bytes(encoded)
    os.chmod(partial, 0o444)
    os.replace(partial, path)
    sidecar_partial = sidecar.with_name(sidecar.name + ".partial")
    sidecar_partial.write_text(f"{digest}  {path.name}\n", encoding="ascii")
    os.chmod(sidecar_partial, 0o444)
    os.replace(sidecar_partial, sidecar)
    return path


def verify_content_addressed(path: Path) -> str:
    path = require_file(path, "content-addressed artifact")
    sidecar = require_file(path.with_suffix(path.suffix + ".sha256"), "SHA-256 sidecar")
    fields = sidecar.read_text(encoding="ascii").strip().split()
    if len(fields) != 2 or fields[1] != path.name:
        raise GateError(f"malformed SHA-256 sidecar: {sidecar}")
    actual = sha256_file(path)
    if fields[0] != actual:
        raise GateError(f"artifact SHA-256 mismatch: {path}")
    if f".{actual}.json" not in path.name:
        raise GateError(f"artifact filename is not content-addressed: {path}")
    return actual


def source_files(repo_root: Path) -> list[Path]:
    files: set[Path] = set()
    for pattern in SOURCE_PATTERNS:
        files.update(path.resolve() for path in repo_root.glob(pattern) if path.is_file())
    if not files:
        raise GateError(f"no release source files found beneath {repo_root}")
    return sorted(files)


def source_manifest(repo_root: Path) -> tuple[list[dict[str, Any]], str]:
    """Return the exact release source inventory and its ordered bundle digest."""
    repo_root = repo_root.resolve()
    entries: list[dict[str, Any]] = []
    bundle = hashlib.sha256()
    for path in source_files(repo_root):
        relative = str(path.relative_to(repo_root))
        digest = sha256_file(path)
        entries.append({"path": relative, "sha256": digest, "bytes": path.stat().st_size})
        bundle.update(relative.encode("utf-8") + b"\0" + digest.encode("ascii") + b"\0")
    return entries, bundle.hexdigest()


def git_output(repo_root: Path, *args: str) -> str:
    result = subprocess.run(
        ("git", *args),
        cwd=repo_root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def build_identity(repo_root: Path, executable: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    executable = require_file(executable, "Frontier executable")
    entries, bundle_sha256 = source_manifest(repo_root)
    try:
        revision = git_output(repo_root, "rev-parse", "HEAD")
        branch = git_output(repo_root, "branch", "--show-current")
        status = git_output(repo_root, "status", "--short", "--untracked-files=all")
    except (OSError, subprocess.CalledProcessError) as exc:
        raise GateError(f"unable to record git identity: {exc}") from exc
    return {
        "schema_version": 1,
        "artifact_type": "frontier-release-identity",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "repo_root": str(repo_root),
        "git_revision": revision,
        "git_branch": branch,
        "git_status_at_freeze": status.splitlines(),
        "source_bundle_sha256": bundle_sha256,
        "source_files": entries,
        "executable": {
            "path": str(executable),
            "sha256": sha256_file(executable),
            "bytes": executable.stat().st_size,
        },
    }


def verify_identity(identity_path: Path, repo_root: Path, executable: Path) -> dict[str, Any]:
    identity_sha = verify_content_addressed(identity_path)
    identity = load_json(identity_path)
    if identity.get("artifact_type") != "frontier-release-identity":
        raise GateError(f"not a Frontier release identity: {identity_path}")
    repo_root = repo_root.resolve()
    executable = require_file(executable, "Frontier executable")
    frozen_entries = identity.get("source_files")
    if not isinstance(frozen_entries, list) or not frozen_entries:
        raise GateError("release identity has no frozen source inventory")
    frozen_by_path: dict[str, dict[str, Any]] = {}
    for entry in frozen_entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise GateError("release identity has a malformed source entry")
        relative = entry["path"]
        if relative in frozen_by_path:
            raise GateError(f"release identity repeats frozen source: {relative}")
        frozen_by_path[relative] = entry

    current_entries, current_bundle_sha256 = source_manifest(repo_root)
    current_by_path = {entry["path"]: entry for entry in current_entries}
    added = sorted(set(current_by_path) - set(frozen_by_path))
    removed = sorted(set(frozen_by_path) - set(current_by_path))
    if added or removed:
        details = []
        if added:
            details.append("added: " + ", ".join(added))
        if removed:
            details.append("removed: " + ", ".join(removed))
        raise GateError("frozen source inventory mismatch (" + "; ".join(details) + ")")
    for relative, current in current_by_path.items():
        frozen = frozen_by_path[relative]
        if current["sha256"] != frozen.get("sha256"):
            raise GateError(f"frozen source checksum mismatch: {repo_root / relative}")
        if current["bytes"] != frozen.get("bytes"):
            raise GateError(f"frozen source byte count mismatch: {repo_root / relative}")
    if current_bundle_sha256 != identity.get("source_bundle_sha256"):
        raise GateError("frozen source bundle checksum mismatch")
    expected_executable = identity.get("executable", {})
    if sha256_file(executable) != expected_executable.get("sha256"):
        raise GateError(f"Frontier executable checksum mismatch: {executable}")
    identity["_artifact_sha256"] = identity_sha
    return identity


def require_number(mapping: dict[str, Any], key: str) -> float:
    value = mapping.get(key)
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise GateError(f"{key} is not finite: {value!r}")
    return float(value)


def attest_run(args: argparse.Namespace) -> Path:
    identity = verify_identity(args.identity, args.repo_root, args.executable)
    output_dir = args.output_dir.resolve()
    manifest_path = require_file(output_dir / "rebuild_manifest.json", "rebuild manifest")
    validator_path = require_file(args.validator_report, "validator report")
    manifest = load_json(manifest_path)
    validator = load_json(validator_path)
    if validator.get("status") != "pass":
        raise GateError(f"validator report did not pass: {validator_path}")
    if Path(validator.get("inputs", {}).get("rebuilt_output_dir", "")).resolve() != output_dir:
        raise GateError("validator report does not identify the attested output directory")
    if manifest.get("geometry_to_logical_key_validated") is not True:
        raise GateError("rebuild manifest lacks passing geometry-to-logical-key validation")

    budgets = load_json(args.budget_file).get("budgets", {})
    budget = budgets.get(args.budget_key)
    if not isinstance(budget, dict):
        raise GateError(f"performance budget not found: {args.budget_key}")
    elapsed = require_number(manifest, "elapsed_seconds")
    cells = require_number(manifest, "cells_processed")
    payload_bytes = require_number(manifest, "payload_bytes_read")
    planned_peak = require_number(manifest, "planned_peak_buffer_bytes_per_rank") / GIB
    metrics = {
        "elapsed_seconds": elapsed,
        "cells_per_second": cells / elapsed,
        "payload_gib_per_second": payload_bytes / GIB / elapsed,
        "node_hours": args.nodes * elapsed / 3600.0,
        "planned_peak_buffer_gib_per_rank": planned_peak,
    }
    checks = {
        "required_product_set": manifest.get("product_set")
        == budget.get("required_product_set"),
        "required_nodes": args.nodes == budget.get("required_nodes"),
        "max_elapsed_seconds": elapsed <= require_number(budget, "max_elapsed_seconds"),
        "max_node_hours": metrics["node_hours"]
        <= require_number(budget, "max_node_hours"),
        "min_cells_per_second": metrics["cells_per_second"]
        >= require_number(budget, "min_cells_per_second"),
        "max_planned_peak_buffer_gib_per_rank": planned_peak
        <= require_number(budget, "max_planned_peak_buffer_gib_per_rank"),
        "minimum_timeout_margin_seconds": elapsed
        <= require_number(budget, "job_timeout_seconds")
        - require_number(budget, "minimum_timeout_margin_seconds"),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise GateError(f"performance/release checks failed: {', '.join(failed)}")

    payload = {
        "schema_version": 1,
        "artifact_type": "frontier-run-pass",
        "status": "pass",
        "kind": args.kind,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "release_identity": {
            "path": str(args.identity.resolve()),
            "artifact_sha256": identity["_artifact_sha256"],
            "source_bundle_sha256": identity["source_bundle_sha256"],
            "executable_sha256": identity["executable"]["sha256"],
        },
        "output_dir": str(output_dir),
        "rebuild_manifest_sha256": sha256_file(manifest_path),
        "validator_report_sha256": sha256_file(validator_path),
        "budget_key": args.budget_key,
        "budget": budget,
        "metrics": metrics,
        "checks": checks,
        "rebuild_manifest": manifest,
        "validator_report": validator,
    }
    return write_content_addressed(args.artifact_dir, args.kind, payload)


def verify_pass(args: argparse.Namespace) -> dict[str, Any]:
    identity = verify_identity(args.identity, args.repo_root, args.executable)
    report_sha = verify_content_addressed(args.pass_report)
    report = load_json(args.pass_report)
    if report.get("artifact_type") != "frontier-run-pass" or report.get("status") != "pass":
        raise GateError(f"not a passing Frontier run report: {args.pass_report}")
    if args.expected_kind and report.get("kind") != args.expected_kind:
        raise GateError(
            f"pass report kind mismatch: expected {args.expected_kind}, found {report.get('kind')}"
        )
    release_identity = report.get("release_identity", {})
    if release_identity.get("artifact_sha256") != identity["_artifact_sha256"]:
        raise GateError("pass report was produced under a different release identity")
    if release_identity.get("executable_sha256") != identity["executable"]["sha256"]:
        raise GateError("pass report was produced by a different executable")
    checks = report.get("checks")
    if not isinstance(checks, dict) or not checks or not all(checks.values()):
        raise GateError("pass report contains a failed or missing release check")
    validator = report.get("validator_report", {})
    if validator.get("status") != "pass":
        raise GateError("embedded validator report did not pass")
    report["_artifact_sha256"] = report_sha
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    freeze = subparsers.add_parser("freeze")
    freeze.add_argument("--repo-root", type=Path, required=True)
    freeze.add_argument("--executable", type=Path, required=True)
    freeze.add_argument("--artifact-dir", type=Path, required=True)

    verify_identity_parser = subparsers.add_parser("verify-identity")
    verify_identity_parser.add_argument("--identity", type=Path, required=True)
    verify_identity_parser.add_argument("--repo-root", type=Path, required=True)
    verify_identity_parser.add_argument("--executable", type=Path, required=True)

    attest = subparsers.add_parser("attest-run")
    attest.add_argument("--identity", type=Path, required=True)
    attest.add_argument("--repo-root", type=Path, required=True)
    attest.add_argument("--executable", type=Path, required=True)
    attest.add_argument("--output-dir", type=Path, required=True)
    attest.add_argument("--validator-report", type=Path, required=True)
    attest.add_argument("--budget-file", type=Path, required=True)
    attest.add_argument("--budget-key", required=True)
    attest.add_argument("--kind", required=True)
    attest.add_argument("--nodes", type=int, required=True)
    attest.add_argument("--artifact-dir", type=Path, required=True)

    verify = subparsers.add_parser("verify-pass")
    verify.add_argument("--identity", type=Path, required=True)
    verify.add_argument("--repo-root", type=Path, required=True)
    verify.add_argument("--executable", type=Path, required=True)
    verify.add_argument("--pass-report", type=Path, required=True)
    verify.add_argument("--expected-kind")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "freeze":
            identity = build_identity(args.repo_root, args.executable)
            path = write_content_addressed(
                args.artifact_dir, "frontier-release-identity", identity
            )
            print(path)
        elif args.command == "verify-identity":
            identity = verify_identity(args.identity, args.repo_root, args.executable)
            print(
                f"verified identity {identity['_artifact_sha256']} "
                f"executable {identity['executable']['sha256']}"
            )
        elif args.command == "attest-run":
            print(attest_run(args))
        elif args.command == "verify-pass":
            report = verify_pass(args)
            print(
                f"verified {report['kind']} pass {report['_artifact_sha256']} "
                f"for {report['output_dir']}"
            )
        else:
            raise GateError(f"unsupported command: {args.command}")
    except Exception as exc:
        print(f"release gate error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
