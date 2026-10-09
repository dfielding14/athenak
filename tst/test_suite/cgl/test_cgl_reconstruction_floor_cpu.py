"""Call the production directional reconstruction wrappers in a focused Serial build."""

import pytest

from test_cgl_c2p_pressure_floor_cpu import REPO_ROOT, _cmake_cache_value, _run


@pytest.mark.parametrize("single_precision", [False, True])
def test_cgl_directional_reconstruction_pressure_floors(tmp_path, single_precision):
    build = tmp_path / "build"
    _run([
        "cmake", "-S", str(REPO_ROOT), "-B", str(build), "-G", "Unix Makefiles",
        "-DPROBLEM=unit_tests/cgl_reconstruction_floor_test",
        "-DAthena_ENABLE_MPI=OFF", "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
        f"-DAthena_SINGLE_PRECISION={'ON' if single_precision else 'OFF'}",
    ])
    object_target = (
        "src/CMakeFiles/athena.dir/pgen/unit_tests/cgl_reconstruction_floor_test.cpp.o"
    )
    _run([
        "make", "-f", "src/CMakeFiles/athena.dir/build.make", object_target, "-j4",
    ], cwd=build)
    _run([
        "cmake", "--build", str(build), "--target", "kokkoscontainers", "kokkoscore", "-j4",
    ])
    main = tmp_path / "main.cpp"
    main.write_text(
        "void RunCglReconstructionFloorStandaloneChecks();\n"
        "int main() { RunCglReconstructionFloorStandaloneChecks(); }\n"
    )
    executable = tmp_path / "cgl_reconstruction_floor_test"
    _run([
        _cmake_cache_value(build, "CMAKE_CXX_COMPILER"), str(main),
        str(build / object_target),
        str(build / "kokkos/containers/src/libkokkoscontainers.a"),
        str(build / "kokkos/core/src/libkokkoscore.a"),
        "-std=c++17", "-ldl", "-pthread", "-o", str(executable),
    ])
    result = _run([str(executable)])
    assert "CGL directional reconstruction floor checks passed" in result.stdout
