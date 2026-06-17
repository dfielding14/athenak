"""MPI regression test for passive-scalar forcing reductions."""

import subprocess
from pathlib import Path

import numpy as np

from test_suite.scalar.test_scalar_driving_cpu import read_history


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "scalar_driving_source_rms.athinput"
ATHENA = Path.cwd() / "athena"


def test_scalar_forcing_normalization_under_mpi(tmp_path):
    """Global projection preserves zero mean and requested RMS under MPI."""
    output_dir = tmp_path / "mpi"
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
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    history = read_history(output_dir / "scalar_force_test.user.hst")
    active = history["time"] > 0.0
    rho = history["rho"][active]
    np.testing.assert_allclose(history["rf_s0"][active] / rho, 0.0, atol=2.0e-15)
    np.testing.assert_allclose(
        np.sqrt(history["rf2_s0"][active] / rho), 0.2, rtol=2.0e-14
    )
