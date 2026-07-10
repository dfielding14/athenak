"""Tests for the content-addressed Frontier release gate."""

from __future__ import annotations

import copy
import importlib.util
import subprocess
from pathlib import Path

import pytest


TEST_DIR = Path(__file__).resolve().parent
RELEASE_GATE_PATH = TEST_DIR.parent / "jobs" / "release_gate.py"
SPEC = importlib.util.spec_from_file_location("gotham_pdf_release_gate", RELEASE_GATE_PATH)
assert SPEC is not None and SPEC.loader is not None
release_gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_gate)


def git(repo: Path, *args: str) -> None:
    subprocess.run(
        ("git", *args),
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


@pytest.fixture()
def frozen_release(tmp_path: Path) -> tuple[Path, Path, dict, Path]:
    repo = tmp_path / "repo"
    reducer = repo / "tools" / "gotham_pdf_rebuild"
    cooling_tables = repo / "src" / "srcterms" / "cooling_tables.hpp"
    reducer.mkdir(parents=True)
    cooling_tables.parent.mkdir(parents=True)
    (reducer / "CMakeLists.txt").write_text("project(test)\n", encoding="utf-8")
    (reducer / "gotham_pdf_rebuild.cpp").write_text(
        "int main() { return 0; }\n", encoding="utf-8"
    )
    cooling_tables.write_text("// frozen cooling tables\n", encoding="utf-8")

    git(repo, "init")
    git(repo, "add", ".")
    git(
        repo,
        "-c",
        "user.name=GOTHAM release test",
        "-c",
        "user.email=gotham-release-test@example.invalid",
        "commit",
        "-m",
        "fixture",
    )

    executable = tmp_path / "gotham_pdf_rebuild"
    executable.write_bytes(b"frozen executable\n")
    identity = release_gate.build_identity(repo, executable)
    identity_path = release_gate.write_content_addressed(
        tmp_path / "release", "frontier-release-identity", identity
    )
    return repo, executable, identity, identity_path


def test_verify_identity_requires_exact_current_source_inventory(
    frozen_release: tuple[Path, Path, dict, Path],
) -> None:
    repo, executable, _, identity_path = frozen_release
    release_gate.verify_identity(identity_path, repo, executable)

    added = repo / "tools" / "gotham_pdf_rebuild" / "jobs" / "late_guard.sh"
    added.parent.mkdir(parents=True)
    added.write_text("#!/bin/sh\n", encoding="utf-8")

    with pytest.raises(release_gate.GateError, match=r"inventory mismatch .*late_guard\.sh"):
        release_gate.verify_identity(identity_path, repo, executable)


@pytest.mark.parametrize("field", ["bytes", "source_bundle_sha256"])
def test_verify_identity_rejects_frozen_source_metadata_tampering(
    frozen_release: tuple[Path, Path, dict, Path], tmp_path: Path, field: str
) -> None:
    repo, executable, identity, _ = frozen_release
    tampered = copy.deepcopy(identity)
    if field == "bytes":
        tampered["source_files"][0]["bytes"] += 1
    else:
        tampered[field] = "0" * 64
    identity_path = release_gate.write_content_addressed(
        tmp_path / f"tampered-{field}", "frontier-release-identity", tampered
    )

    with pytest.raises(release_gate.GateError):
        release_gate.verify_identity(identity_path, repo, executable)
