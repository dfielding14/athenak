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

    # The custom pgen also contains the real grid wall check.
    _run(["cmake", "--build", str(build_dir), "--target", "kokkoscore", "-j4"])

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


@pytest.mark.parametrize("backup", (False, True), ids=("fluid_wall", "backup_walls"))
def test_cgl_production_wall_encoding_survives_independent_c2p(tmp_path, backup):
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
    source = """<job>
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
"""
    if backup:
        source += """
<mhd>
backup_limiters = true
firehose_threshold = 1.4
firehose_backup_factor = 1.1
mirror_threshold = 0.7
mirror_backup_factor = 1.3
nu_coll = 0.3
limiter_nu_coll = 7
"""
    input_file.write_text(source)
    launcher = shlex.split(os.environ.get("ATHENAK_CGL_PRODUCTION_WALL_LAUNCHER", ""))
    command = launcher + [str(executable), "-i", str(input_file)]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True,
                            text=True, check=False)
    (tmp_path / "production-wall.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CGL production wall encoding/restart checks passed: 6 states" in result.stdout
    assert "CGL LF conversion wall encoding checks passed: 6 states" in result.stdout
    assert "CGL production tiny-rate cancellation checks passed: 2 states" in result.stdout
    if backup:
        assert "CGL asymmetric backup wall/full-rate checks passed" in result.stdout
    assert (tmp_path / "wall_state.bin").stat().st_size == 6 * 6 * 8


def test_cgl_strict_lf_rejects_nonfinite_state_before_weak_field_refresh(tmp_path):
    """Exercise the actual lightweight EOS kernel, including HIP when supplied."""
    supplied = os.environ.get("ATHENAK_CGL_STORAGE_BINARY")
    if supplied:
        executable = Path(supplied).resolve()
        assert executable.is_file()
    else:
        build_dir = tmp_path / "storage-build"
        _run(["cmake", "-S", str(REPO_ROOT), "-B", str(build_dir),
              "-DPROBLEM=unit_tests/cgl_c2p_storage_test",
              "-DAthena_ENABLE_MPI=OFF", "-DKokkos_ENABLE_HIP=OFF",
              "-DKokkos_ENABLE_SERIAL=ON", "-DCMAKE_BUILD_TYPE=Release"])
        _run(["cmake", "--build", str(build_dir), "--target", "athena", "-j4"])
        executable = build_dir / "src/athena"
    launcher = shlex.split(os.environ.get("ATHENAK_CGL_STORAGE_LAUNCHER", ""))
    source = (REPO_ROOT / "inputs/unit_tests/cgl_lf_sts_grad_b.athinput").read_text()
    source = source.replace("<mhd>", "<mhd>\ncgl_lf_strict_admissibility = false", 1)
    source = source.replace("<problem>", "<problem>\nlf_nonfinite_case = 0", 1)
    input_file = tmp_path / "lf-nonfinite-refresh.athinput"
    input_file.write_text(source)
    failure = "strict admissibility failed before primitive refresh"
    for strict in (False, True):
        # Healthy control, then E and moment each poisoned with NaN, +Inf, and -Inf.
        for test_case in range(7):
            command = launcher + [str(executable), "-i", str(input_file),
                                  "time/nlim=0", f"problem/lf_nonfinite_case={test_case}",
                                  f"mhd/cgl_lf_strict_admissibility={str(strict).lower()}"]
            result = subprocess.run(command, cwd=tmp_path, capture_output=True,
                                    text=True, check=False)
            output = result.stdout + result.stderr
            (tmp_path / f"lf-refresh-{strict}-{test_case}.log").write_text(output)
            if strict and test_case != 0:
                assert result.returncode != 0, output
                assert failure in output, output
                assert "CGL LF nonfinite refresh returned" not in output
            else:
                assert result.returncode == 0, output
                assert f"CGL LF nonfinite refresh returned: case={test_case}" in output
                assert failure not in output
