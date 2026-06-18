"""MPI regression for fully 2D viscous turbulence and SGS output."""

from pathlib import Path
import subprocess
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "turb_sgs_2d.athinput"
ATHENA = Path.cwd() / "athena"
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary, read_coarsened_binary  # noqa: E402


def latest(path, pattern):
    """Return the lexically latest numbered output file."""
    outputs = sorted(path.glob(pattern))
    assert outputs
    return outputs[-1]


def test_viscous_2d_sgs_output_under_mpi(tmp_path):
    """Four ranks preserve the 2D contract and write the viscous SGS state."""
    output_dir = tmp_path / "mpi_viscous_sgs"
    output_dir.mkdir()
    result = subprocess.run(
        [
            "mpirun",
            "-np",
            "4",
            str(ATHENA),
            "-d",
            str(output_dir),
            "-i",
            str(INPUT),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    state = read_binary(str(latest(output_dir / "bin", "*.state.*.bin")))
    force = read_binary(str(latest(output_dir / "bin", "*.force.*.bin")))
    sgs = read_coarsened_binary(
        str(latest(output_dir / "cbin_sgs_2", "*.sgs.*.cbin"))
    )
    assert state["n_mbs"] == 4
    assert sgs["var_names"] == [
        "dens",
        "velx",
        "vely",
        "tau_xx",
        "tau_xy",
        "tau_yy",
    ]
    assert np.all(np.asarray(force["mb_data"]["force3"]) == 0.0)
    assert np.all(np.asarray(state["mb_data"]["mom3"]) == 0.0)
