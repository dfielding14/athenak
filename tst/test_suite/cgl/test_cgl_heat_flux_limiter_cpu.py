"""Focused regression for overflow-safe CGL heat-flux limiting."""

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
def test_cgl_heat_flux_limiter_extreme_ranges(tmp_path, single_precision):
    build_dir = tmp_path / "build"
    configure = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(build_dir),
        "-G",
        "Unix Makefiles",
        "-DPROBLEM=unit_tests/cgl_heat_flux_limiter_test",
        "-DAthena_ENABLE_MPI=OFF",
        "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
    ]
    if single_precision:
        configure.append("-DAthena_SINGLE_PRECISION=ON")
    _run(configure)

    object_target = (
        "src/CMakeFiles/athena.dir/pgen/unit_tests/"
        "cgl_heat_flux_limiter_test.cpp.o"
    )
    _run(
        ["make", "-f", "src/CMakeFiles/athena.dir/build.make", object_target, "-j4"],
        cwd=build_dir,
    )

    main_source = tmp_path / "main.cpp"
    main_source.write_text(
        "void RunCglHeatFluxLimiterChecks();\n"
        "int main() { RunCglHeatFluxLimiterChecks(); }\n"
    )
    executable = tmp_path / "cgl_heat_flux_limiter_test"
    compiler = _cmake_cache_value(build_dir, "CMAKE_CXX_COMPILER")
    _run([
        compiler,
        str(main_source),
        str(build_dir / object_target),
        "-std=c++17",
        "-o",
        str(executable),
    ])
    result = _run([str(executable)])
    assert "CGL heat-flux limiter checks passed" in result.stdout
