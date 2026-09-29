"""Focused custom-pgen regression for generic CGL C2P pressure floors."""

from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]


def _run(command, cwd=None):
    result = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def _cmake_cache_value(build_dir, key):
    prefix = f"{key}:"
    for line in (build_dir / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1]
    raise AssertionError(f"{key} is absent from CMakeCache.txt")


@pytest.mark.parametrize("single_precision", [False, True])
def test_cgl_c2p_pressure_floor_energy_consistency(tmp_path, single_precision):
    build_dir = tmp_path / "build"
    configure = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(build_dir),
        "-G",
        "Unix Makefiles",
        "-DPROBLEM=unit_tests/cgl_c2p_pressure_floor_test",
        "-DAthena_ENABLE_MPI=OFF",
        "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
    ]
    if single_precision:
        configure.append("-DAthena_SINGLE_PRECISION=ON")
    _run(configure)
    # Link the exported checker directly to keep this a single-state unit test.
    test_object_target = (
        "src/CMakeFiles/athena.dir/pgen/unit_tests/"
        "cgl_c2p_pressure_floor_test.cpp.o"
    )
    production_object_target = "src/CMakeFiles/athena.dir/eos/cgl_mhd.cpp.o"
    for target in (test_object_target, production_object_target):
        _run(
            ["make", "-f", "src/CMakeFiles/athena.dir/build.make", target, "-j4"],
            cwd=build_dir,
        )

    main_source = tmp_path / "main.cpp"
    main_source.write_text(
        "void RunCglC2PPressureFloorChecks();\n"
        "int main() { RunCglC2PPressureFloorChecks(); }\n"
    )
    test_object = build_dir / test_object_target
    test_executable = tmp_path / "cgl_c2p_pressure_floor_test"
    compiler = _cmake_cache_value(build_dir, "CMAKE_CXX_COMPILER")
    _run([
        compiler,
        str(main_source),
        str(test_object),
        "-std=c++17",
        "-o",
        str(test_executable),
    ])
    result = _run([str(test_executable)])
    assert "CGL C2P pressure-floor checks passed" in result.stdout
