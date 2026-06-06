"""Focused regression for the corrected active-CGL fast-speed and HLLE route."""

from pathlib import Path
import subprocess


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


def test_cgl_active_hlle_uses_literature_discriminant(tmp_path):
    build_dir = tmp_path / "build"
    _run([
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(build_dir),
        "-G",
        "Unix Makefiles",
        "-DPROBLEM=unit_tests/cgl_fast_speed_test",
        "-DAthena_ENABLE_MPI=OFF",
        "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
    ])
    object_target = (
        "src/CMakeFiles/athena.dir/pgen/unit_tests/cgl_fast_speed_test.cpp.o"
    )
    _run(
        ["make", "-f", "src/CMakeFiles/athena.dir/build.make", object_target, "-j4"],
        cwd=build_dir,
    )
    _run([
        "cmake",
        "--build",
        str(build_dir),
        "--target",
        "kokkoscontainers",
        "kokkosalgorithms",
        "kokkoscore",
        "kokkossimd",
        "-j4",
    ])

    main_source = tmp_path / "main.cpp"
    main_source.write_text(
        "void RunCglFastSpeedStandaloneChecks();\n"
        "int main() { RunCglFastSpeedStandaloneChecks(); }\n"
    )
    executable = tmp_path / "cgl_fast_speed_test"
    compiler = _cmake_cache_value(build_dir, "CMAKE_CXX_COMPILER")
    kokkos_libs = [
        build_dir / "kokkos" / "containers" / "src" / "libkokkoscontainers.a",
        build_dir / "kokkos" / "algorithms" / "src" / "libkokkosalgorithms.a",
        build_dir / "kokkos" / "core" / "src" / "libkokkoscore.a",
        build_dir / "kokkos" / "simd" / "src" / "libkokkossimd.a",
    ]
    _run([
        compiler,
        str(main_source),
        str(build_dir / object_target),
        *(str(path) for path in kokkos_libs),
        "-std=c++17",
        "-ldl",
        "-pthread",
        "-o",
        str(executable),
    ])
    result = _run([str(executable)])
    assert "CGL active reconstructed-HLLE route checks passed" in result.stdout
    assert "CGL fast-speed checks passed" in result.stdout
