"""Direct serial regression coverage for IOWrapper defensive paths."""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
HARNESS_SOURCE = ROOT / "tst" / "test_suite" / "io" / "io_wrapper_serial_harness.cpp"
WRAPPER_SOURCE = ROOT / "src" / "outputs" / "io_wrapper.cpp"


def _read_make_variable(flags_file: Path, name: str):
    prefix = f"{name} = "
    for line in flags_file.read_text().splitlines():
        if line.startswith(prefix):
            return shlex.split(line[len(prefix):])
    raise AssertionError(f"{name} not found in {flags_file}")


def _serial_flags_file():
    candidates = (
        Path.cwd() / "CMakeFiles" / "athena.dir" / "flags.make",
        Path.cwd() / "src" / "CMakeFiles" / "athena.dir" / "flags.make",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise AssertionError("Run from a configured serial CMake build directory.")


@pytest.fixture(scope="session")
def io_wrapper_serial_harness(tmp_path_factory):
    compiler_parts = shlex.split(os.environ.get("CXX", "c++"))
    compiler = shutil.which(compiler_parts[0])
    assert compiler is not None
    flags_file = _serial_flags_file()
    defines = _read_make_variable(flags_file, "CXX_DEFINES")
    config = flags_file.parents[3] / "config.hpp"
    assert config.exists()
    assert "#define MPI_PARALLEL_ENABLED 0" in config.read_text()
    output = (
        tmp_path_factory.mktemp("io_wrapper_serial_harness")
        / "io_wrapper_serial_harness"
    )
    subprocess.run(
        [
            compiler,
            *compiler_parts[1:],
            *defines,
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


def _run_harness(io_wrapper_serial_harness: Path, tmp_path: Path, mode: str):
    return subprocess.run(
        [
            str(io_wrapper_serial_harness),
            mode,
            str(tmp_path / f"{mode}.bin"),
        ],
        input="x",
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_seek_reports_failure_for_nonseekable_stream(io_wrapper_serial_harness, tmp_path):
    proc = _run_harness(io_wrapper_serial_harness, tmp_path, "seek_stdin")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "serial seek failure observed" in proc.stdout


@pytest.mark.parametrize(
    ("mode", "expected"),
    (
        ("maximum_offset", "Seek exceeds off_t range."),
        ("negative_position", "GetPosition returned a negative offset."),
        ("positioned_read_seek_failure", "Read_bytes_at serial seek failed."),
        ("terminal_offset_overflow", "Read_bytes_at exceeds off_t range."),
    ),
)
def test_serial_wrapper_rejects_invalid_positioned_io(
    io_wrapper_serial_harness, tmp_path, mode, expected
):
    proc = _run_harness(io_wrapper_serial_harness, tmp_path, mode)
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)
