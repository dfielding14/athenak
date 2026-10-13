"""Build the dedicated forcing fixture when the test runner uses built-in generators."""

from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]


@pytest.fixture(scope="module", autouse=True)
def forcing_executable(request, tmp_path_factory):
    """Reuse a matching executable, or build one without changing the runner's binary."""
    if request.module.__name__.split(".")[-1] not in (
            "test_turb_driving_cpu", "test_turb_driving_mpicpu"):
        yield
        return
    original = request.module.ATHENA
    binary = original
    mpi = "ON" if request.module.__name__.endswith("_mpicpu") else "OFF"
    configuration = {}
    if binary.is_file():
        result = subprocess.run([str(binary), "-c"], capture_output=True,
                                text=True, check=False, timeout=30)
        if result.returncode == 0:
            configuration = {
                key.strip(): value.strip()
                for key, value in (
                    line.split(":", 1)
                    for line in result.stdout.splitlines() if ":" in line
                )
            }
    if (configuration.get("Problem generator") != "tests/turb_forcing"
            or configuration.get("MPI parallelism") != mpi
            or configuration.get("Floating-point precision") != "double"):
        build_dir = tmp_path_factory.mktemp("turb_forcing_" + mpi.lower())
        commands = [
            ["cmake", "-S", str(REPO_ROOT), "-B", str(build_dir),
             "-DPROBLEM=tests/turb_forcing", "-DCMAKE_BUILD_TYPE=Release",
             "-DAthena_SINGLE_PRECISION=OFF", "-DAthena_ENABLE_MPI=" + mpi,
             "-DAthena_ENABLE_OPENMP=OFF", "-DKokkos_ENABLE_SERIAL=ON",
             "-DKokkos_ENABLE_OPENMP=OFF"],
            ["cmake", "--build", str(build_dir), "--target", "athena", "-j", "2"],
        ]
        for command in commands:
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            assert result.returncode == 0, result.stdout + result.stderr
        binary = build_dir / "src" / "athena"
    request.module.ATHENA = binary
    yield
    request.module.ATHENA = original
