"""Direct MPI regression coverage for IOWrapper defensive paths."""

from pathlib import Path
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = ROOT / "tst" / "test_suite" / "io" / "io_wrapper_mpi_harness.cpp"
WRAPPER_SOURCE = ROOT / "src" / "outputs" / "io_wrapper.cpp"


def _read_make_variable(flags_file: Path, name: str):
    prefix = f"{name} = "
    for line in flags_file.read_text().splitlines():
        if line.startswith(prefix):
            return shlex.split(line[len(prefix):])
    raise AssertionError(f"{name} not found in {flags_file}")


@pytest.fixture(scope="session")
def io_wrapper_harness(tmp_path_factory):
    compiler = shutil.which("mpicxx")
    assert compiler is not None
    flags_file = Path.cwd() / "CMakeFiles" / "athena.dir" / "flags.make"
    assert flags_file.exists()
    output = tmp_path_factory.mktemp("io_wrapper_harness") / "io_wrapper_mpi_harness"
    subprocess.run(
        [
            compiler,
            *_read_make_variable(flags_file, "CXX_DEFINES"),
            *_read_make_variable(flags_file, "CXX_INCLUDES"),
            *_read_make_variable(flags_file, "CXX_FLAGS"),
            str(HARNESS_SOURCE),
            str(WRAPPER_SOURCE),
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return output


def _run_harness(io_wrapper_harness: Path, tmp_path: Path, mode: str, ranks: int):
    return subprocess.run(
        [
            "mpirun",
            "-np",
            str(ranks),
            str(io_wrapper_harness),
            mode,
            str(tmp_path / f"{mode}.bin"),
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_collective_helpers_support_asymmetric_and_zero_byte_ranks(
    io_wrapper_harness, tmp_path
):
    proc = _run_harness(io_wrapper_harness, tmp_path, "collective", 3)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "collective participation ok" in proc.stdout


@pytest.mark.parametrize(
    ("mode", "expected"),
    (
        ("offset_range_write", "MPI_File_write_at_all exceeds MPI_Offset range."),
        ("offset_range_read", "MPI_File_read_at_all exceeds MPI_Offset range."),
        ("byte_count_overflow", "Write_any_type_at_all byte count overflow."),
        (
            "chunk_limit_disagreement",
            "ATHENAK_TEST_MAX_MPI_BYTES must have the same value",
        ),
    ),
)
def test_wrapper_rejects_invalid_mpi_io_before_progress(
    io_wrapper_harness, tmp_path, mode, expected
):
    proc = _run_harness(io_wrapper_harness, tmp_path, mode, 2)
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)
