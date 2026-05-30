"""Direct CPU regression coverage for checked restart-layout arithmetic."""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = ROOT / "tst" / "test_suite" / "io" / "restart_layout_harness.cpp"


@pytest.fixture(scope="session")
def restart_layout_harness(tmp_path_factory):
    compiler_parts = shlex.split(os.environ.get("CXX", "c++"))
    compiler = shutil.which(compiler_parts[0])
    assert compiler is not None
    output = tmp_path_factory.mktemp("restart_layout_harness") / "restart_layout_harness"
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


@pytest.mark.parametrize(
    "mode",
    (
        "large_field",
        "multiply_overflow",
        "add_overflow",
        "offset_overflow",
        "payload_bytes_overflow",
        "subtract_underflow",
        "manifest_budget_overflow",
    ),
)
def test_restart_layout_checked_arithmetic(restart_layout_harness, mode):
    subprocess.run(
        [str(restart_layout_harness), mode],
        check=True,
        capture_output=True,
        text=True,
    )
