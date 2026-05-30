"""Focused checked-arithmetic coverage for coarsened-binary layouts."""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = ROOT / "tst" / "test_suite" / "io" / "cbin_layout_harness.cpp"


@pytest.fixture(scope="session")
def cbin_layout_harness(tmp_path_factory):
    compiler_parts = shlex.split(os.environ.get("CXX", "c++"))
    compiler = shutil.which(compiler_parts[0])
    assert compiler is not None
    output = tmp_path_factory.mktemp("cbin_layout_harness") / "cbin_layout_harness"
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
        "valid",
        "invalid_divisibility",
        "factor_cube_overflow",
        "coarsening_range_overflow",
        "normalization_range_overflow",
        "allocation_overflow",
        "allocation_byte_overflow",
        "wide_kernel_intermediates",
    ),
)
def test_cbin_layout_contract(cbin_layout_harness, mode):
    proc = subprocess.run(
        [str(cbin_layout_harness), mode],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
