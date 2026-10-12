"""MPI forcing regressions; build with -DPROBLEM=tests/turb_forcing."""

import subprocess
from pathlib import Path

import pytest

from test_suite.turb.test_turb_driving_cpu import latest_force, read_force_blocks
from test_suite.turb.test_turb_driving_cpu import rms_acceleration
from test_suite.turb.test_turb_driving_cpu import assert_step_power, forcing_input


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "turb_driving_amr_exclude.athinput"
ATHENA = Path.cwd() / "athena"


def test_amr_localized_driving_under_mpi(tmp_path):
    """AMR reduction and rendering preserve requested RMS forcing under MPI."""
    output_dir = tmp_path / "mpi_amr"
    output_dir.mkdir(parents=True)
    command = [
        "mpirun",
        "-np",
        "2",
        str(ATHENA),
        "-d",
        str(output_dir),
        "-i",
        str(INPUT),
    ]
    result = subprocess.run(command, capture_output=True, text=True,
                            check=False, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    output = read_force_blocks(latest_force(output_dir))
    assert max(block["level"] for block in output["blocks"]) > 0
    assert rms_acceleration(output) == pytest.approx(0.20, rel=2.0e-6)


def test_periodic_power_under_mpi(tmp_path):
    """MPI reductions supply the requested total power over a two-unit-volume box."""
    path = forcing_input(tmp_path, hydro_eos=None, mhd_eos="ideal")
    run_dir = tmp_path / "mpi_power"
    run_dir.mkdir()
    command = ["mpirun", "-np", "2", str(ATHENA), "-d", str(run_dir), "-i", str(path),
               "time/integrator=rk3", "mesh/x1max=1.5"]
    result = subprocess.run(command, capture_output=True, text=True,
                            check=False, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert_step_power(run_dir, "mhd", power=0.2)
