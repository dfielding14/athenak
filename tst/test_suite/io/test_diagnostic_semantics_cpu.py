"""Direct analytic coverage for generic output diagnostic semantics."""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = ROOT / "tst" / "test_suite" / "io" / "diagnostic_semantics_harness.cpp"


@pytest.fixture(scope="session")
def diagnostic_semantics_harness(tmp_path_factory):
    compiler_parts = shlex.split(os.environ.get("CXX", "c++"))
    compiler = shutil.which(compiler_parts[0])
    assert compiler is not None
    output = (
        tmp_path_factory.mktemp("diagnostic_semantics_harness")
        / "diagnostic_semantics_harness"
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


def test_diagnostic_semantics_contract(diagnostic_semantics_harness):
    proc = subprocess.run(
        [str(diagnostic_semantics_harness)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
