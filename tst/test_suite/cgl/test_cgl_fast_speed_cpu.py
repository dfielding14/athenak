"""Focused custom-pgen regression for the CGL fast-magnetosonic speed."""

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


def test_cgl_fast_speed_uses_literature_discriminant(tmp_path):
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

    main_source = tmp_path / "main.cpp"
    main_source.write_text(
        "void RunCglFastSpeedChecks();\n"
        "int main() { RunCglFastSpeedChecks(); }\n"
    )
    executable = tmp_path / "cgl_fast_speed_test"
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
    assert "CGL fast-speed checks passed" in result.stdout
