"""Focused direct coverage for shared output-file helpers."""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = (
    ROOT / "tst" / "test_suite" / "io" / "rcp02_output_file_utils_harness.cpp"
)


@pytest.fixture(scope="session")
def rcp02_output_file_utils_harness(tmp_path_factory):
    compiler_parts = shlex.split(os.environ.get("CXX", "c++"))
    compiler = shutil.which(compiler_parts[0])
    assert compiler is not None
    output = (
        tmp_path_factory.mktemp("rcp02_output_file_utils_harness")
        / "rcp02_output_file_utils_harness"
    )
    subprocess.run(
        [
            compiler,
            *compiler_parts[1:],
            "-std=c++17",
            "-I",
            str(ROOT / "src"),
            str(HARNESS_SOURCE),
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return output


def _run_harness(rcp02_output_file_utils_harness: Path, tmp_path: Path, mode: str):
    root = tmp_path / mode
    root.mkdir()
    return subprocess.run(
        [str(rcp02_output_file_utils_harness), mode, str(root)],
        capture_output=True,
        text=True,
        timeout=30,
    )


@pytest.mark.parametrize(
    "mode",
    (
        "existing_directory",
        "conflicting_non_directory",
        "sequence_formatting",
        "sequence_rejections",
        "advance_boundary",
        "stale_temporary_replacement",
        "failed_publication_cleanup",
        "failed_owned_path_discard",
        "lexical_target_normalization",
        "adjacent_radius_tokens",
    ),
)
def test_output_file_utils_contract(
    rcp02_output_file_utils_harness, tmp_path, mode
):
    proc = _run_harness(rcp02_output_file_utils_harness, tmp_path, mode)
    assert proc.returncode == 0, proc.stdout + proc.stderr
