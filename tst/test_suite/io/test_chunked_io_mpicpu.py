"""MPI regressions for forced-small chunked I/O and restart metadata broadcasts."""

import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402


INPUT_FILE = ROOT / "tst" / "inputs" / "io_finalization_timing.athinput"
MAX_MPI_BYTES_ENV = "ATHENAK_TEST_MAX_MPI_BYTES"


def _env(max_mpi_bytes: str):
    env = os.environ.copy()
    env[MAX_MPI_BYTES_ENV] = max_mpi_bytes
    return env


def _run_initial(run_dir: Path, max_mpi_bytes: str = "7"):
    run_dir.mkdir(exist_ok=True)
    return subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "mesh/nx1=48",
            "meshblock/nx1=16",
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
        env=_env(max_mpi_bytes),
        timeout=90,
    )


def test_forced_small_chunks_shared_restart_round_trip_and_truncation(tmp_path: Path):
    run_dir = tmp_path / "initial"
    _run_initial(run_dir)

    binary = run_dir / "bin" / "io_policy.diag.00000.bin"
    terminal_restart = run_dir / "rst" / "io_policy.00001.rst"
    assert bin_convert.read_binary(str(binary))["n_mbs"] == 3
    assert terminal_restart.exists()
    original_binary = binary.read_bytes()
    original_restart = terminal_restart.read_bytes()

    # Reusing the directory forces MPI output files through the truncation path.
    _run_initial(run_dir)
    assert binary.read_bytes() == original_binary
    assert terminal_restart.read_bytes() == original_restart

    resume_dir = tmp_path / "resume"
    resume_dir.mkdir()
    subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-r",
            str(terminal_restart.resolve()),
            "-d",
            str(resume_dir),
            "time/nlim=2",
            "time/tlim=0.02",
            "time/final_output_policy=none",
        ],
        check=True,
        capture_output=True,
        text=True,
        env=_env("7"),
        timeout=90,
    )


@pytest.mark.parametrize("max_mpi_bytes", ("", "bogus", "0", "-1", "2147483648"))
def test_invalid_forced_chunk_limit_fails_explicitly(tmp_path: Path, max_mpi_bytes: str):
    run_dir = tmp_path / "invalid"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "time/final_output_policy=none",
        ],
        capture_output=True,
        text=True,
        env=_env(max_mpi_bytes),
        timeout=30,
    )
    assert proc.returncode != 0
    assert (
        f"{MAX_MPI_BYTES_ENV} must be a decimal integer between 1 and INT_MAX."
        in (proc.stdout + proc.stderr)
    )
