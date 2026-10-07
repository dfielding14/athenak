"""Focused custom-pgen regression for generic CGL C2P pressure floors."""

from pathlib import Path
import os
import shlex
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
    # The shared passive decoder uses Kokkos::abort for an invalid encoding.
    # Link the real implementation even though these checks exercise active CGL.
    _run([
        "cmake", "--build", str(build_dir), "--target", "kokkoscore", "-j4",
    ])

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
        str(build_dir / "kokkos/core/src/libkokkoscore.a"),
        "-std=c++17",
        "-ldl",
        "-pthread",
        "-o",
        str(test_executable),
    ])
    result = _run([str(test_executable)])
    assert "CGL C2P pressure-floor checks passed" in result.stdout


def test_cgl_production_wall_encoding_survives_independent_c2p(tmp_path):
    """Use actual CGLMHD::Collisions, including its compiled grid-kernel context.

    The optional custom-pgen executable/launcher allow this exact check on HIP.
    The default builds the ordinary serial custom pgen from current sources.
    """
    supplied = os.environ.get("ATHENAK_CGL_PRODUCTION_WALL_BINARY")
    if supplied:
        executable = Path(supplied).resolve()
        assert executable.is_file()
    else:
        build_dir = tmp_path / "production-build"
        _run(["cmake", "-S", str(REPO_ROOT), "-B", str(build_dir),
              "-DPROBLEM=unit_tests/cgl_c2p_pressure_floor_test",
              "-DAthena_ENABLE_MPI=OFF", "-DKokkos_ENABLE_HIP=OFF",
              "-DKokkos_ENABLE_SERIAL=ON", "-DCMAKE_BUILD_TYPE=Release"])
        _run(["cmake", "--build", str(build_dir), "--target", "athena", "-j4"])
        executable = build_dir / "src/athena"
    input_file = tmp_path / "production-wall.athinput"
    input_file.write_text("""<job>
basename = production_wall
<mesh>
nghost = 2
nx1 = 8
nx2 = 1
nx3 = 1
x1min = 0
x1max = 1
x2min = 0
x2max = 1
x3min = 0
x3max = 1
ix1_bc = periodic
ox1_bc = periodic
ix2_bc = periodic
ox2_bc = periodic
ix3_bc = periodic
ox3_bc = periodic
<meshblock>
nx1 = 8
nx2 = 1
nx3 = 1
<time>
evolution = dynamic
integrator = rk2
nlim = 0
tlim = 0.001
cfl_number = 0.3
<mhd>
eos = cgl
passive = false
reconstruct = plm
rsolver = hlle
fofc = false
mirror_limiter = true
firehose_limiter = true
mirror_threshold = 1
firehose_threshold = 2
backup_limiters = false
limiter_nu_coll = 1e10
nu_coll = 0
dfloor = 1e-12
pfloor = 1e-12
bfloor = 1e-10
<problem>
""")
    launcher = shlex.split(os.environ.get("ATHENAK_CGL_PRODUCTION_WALL_LAUNCHER", ""))
    command = launcher + [str(executable), "-i", str(input_file)]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True,
                            text=True, check=False)
    (tmp_path / "production-wall.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CGL production wall encoding/restart checks passed: 6 states" in result.stdout
    assert (tmp_path / "wall_state.bin").stat().st_size == 6 * 6 * 8
    old = os.environ.get("ATHENAK_CGL_PRODUCTION_WALL_OLD_BINARY")
    if old:
        control = tmp_path / "old-negative-control"
        control.mkdir()
        result = subprocess.run(launcher + [old, "-i", str(input_file)],
                                cwd=control, capture_output=True, text=True,
                                check=False)
        (control / "production-wall.log").write_text(result.stdout + result.stderr)
        assert result.returncode != 0, "old production wall encoding unexpectedly passed"
        assert "production wall survives independent C2P" in result.stdout
