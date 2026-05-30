"""Checksummed lineage companions for retained exploratory PIC artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Iterable


PUBLICATION_DIR = Path(__file__).resolve().parent
REPO_ROOT = PUBLICATION_DIR.parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fobj:
        for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, str]:
    resolved = path.resolve()
    return {"path": str(resolved), "sha256": sha256(resolved)}


def file_records(paths: Iterable[Path]) -> list[dict[str, str]]:
    unique = {path.resolve() for path in paths}
    missing = [path for path in sorted(unique) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing requested lineage files: " + ", ".join(str(path) for path in missing)
        )
    return [file_record(path) for path in sorted(unique)]


def _git_output(*args: str) -> str:
    proc = subprocess.run(
        ["git", "-C", str(REPO_ROOT), *args],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to read Git revision provenance: " + " ".join(args))
    return proc.stdout.strip()


def write_companion_manifest(
    path: Path,
    *,
    generator: str,
    physical_model: str,
    inputs: Iterable[Path],
    outputs: Iterable[Path],
    executables: Iterable[Path] = (),
) -> None:
    generator_source = (PUBLICATION_DIR / generator).resolve()
    if not generator_source.is_file():
        raise FileNotFoundError(f"Missing lineage generator source: {generator_source}")
    git_commit = _git_output("rev-parse", "HEAD")
    if not git_commit:
        raise RuntimeError("Git revision provenance returned an empty commit")
    git_branch = _git_output("rev-parse", "--abbrev-ref", "HEAD")
    git_status = _git_output("status", "--short").splitlines()
    manifest = {
        "schema_version": 1,
        "generator": generator,
        "generator_source": str(generator_source),
        "generator_source_sha256": sha256(generator_source),
        "lineage_helper_source": str(Path(__file__).resolve()),
        "lineage_helper_source_sha256": sha256(Path(__file__).resolve()),
        "git_commit": git_commit,
        "git_branch": git_branch,
        "git_status": git_status or ["<clean>"],
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_engineering_artifact",
        "physical_model": physical_model,
        "inputs": file_records(inputs),
        "outputs": file_records(outputs),
        "executables": file_records(executables),
    }
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
