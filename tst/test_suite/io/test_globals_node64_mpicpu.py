"""Direct MPI coverage for 64-bit node-local collective helpers."""

from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = ROOT / "tst" / "test_suite" / "io" / "globals_node64_harness.cpp"


@pytest.fixture(scope="session")
def globals_node64_harness(tmp_path_factory):
    compiler = shutil.which("mpicxx")
    launcher = shutil.which("mpirun")
    assert compiler is not None
    assert launcher is not None
    build_root = Path.cwd().parent
    config = build_root / "config.hpp"
    assert config.is_file(), (
        "run MPI IO tests from the configured MPI build tree's src directory"
    )
    output = tmp_path_factory.mktemp("globals_node64_harness") / "globals_node64_harness"
    subprocess.run(
        [
            compiler,
            "-std=c++17",
            "-I",
            str(build_root),
            "-I",
            str(ROOT / "src"),
            str(HARNESS_SOURCE),
            str(ROOT / "src" / "globals.cpp"),
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return launcher, output


def test_node_collective_helpers_preserve_values_above_int_max(globals_node64_harness):
    launcher, harness = globals_node64_harness
    proc = subprocess.run(
        [launcher, "-np", "2", str(harness)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
