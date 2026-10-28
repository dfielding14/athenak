"""MPI regression test for AMR turbulence driving normalization."""

import subprocess
from pathlib import Path

import numpy as np
import pytest

from test_suite.turb.test_turb_driving_cpu import (
    assert_equal_force_outputs, latest_force, read_force_blocks,
)
from test_suite.turb.test_turb_driving_cpu import rms_acceleration, staged_turb_input


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
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    output = read_force_blocks(latest_force(output_dir))
    assert max(block["level"] for block in output["blocks"]) > 0
    assert rms_acceleration(output) == pytest.approx(0.20, rel=2.0e-6)


def test_full_step_power_agrees_across_mpi_decomposition(tmp_path):
    """Real RK3 cycles inject the prescribed power on one and four ranks."""
    input_file = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": "record_injected_work = true",
        "output1": "data_format = %24.16e",
    })
    histories = []
    for ranks in (1, 4):
        output_dir = tmp_path / f"ranks_{ranks}"
        output_dir.mkdir()
        result = subprocess.run(
            ["mpirun", "-np", str(ranks), str(ATHENA), "-d", str(output_dir),
             "-i", str(input_file), "time/integrator=rk3"],
            capture_output=True, text=True, check=False,
        )
        (output_dir / "run.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stdout + result.stderr
        history = next(output_dir.glob("*.hydro.hst"))
        names = [field.split("=", 1)[1] for field in history.read_text().splitlines()[1].split()
                 if "=" in field]
        data = np.loadtxt(history)
        times, energy = data[:, names.index("time")], data[:, names.index("tot-E")]
        dt, work = np.diff(times), np.diff(energy)
        assert np.all(dt >= 0.0) and np.count_nonzero(dt > 0.0) == 6
        np.testing.assert_array_equal(work[dt == 0.0], 0.0)
        np.testing.assert_allclose(work[dt > 0.0]/dt[dt > 0.0], 0.1,
                                   rtol=2e-7, atol=2e-8)
        histories.append((times, energy))
    np.testing.assert_allclose(histories[0], histories[1], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("single_file_per_rank", [False, True])
def test_perpendicular_forcing_restart_under_mpi(tmp_path, single_file_per_rank):
    """Shared and rank-local checkpoints retain the exact two-rank OU trajectory."""
    input_file = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": "projection_policy = mks24_alfvenic_perpendicular\n"
                        "spectrum = power_law\nrecord_injected_work = true\n"
                        "min_kx = -3\nmax_kx = 3\nmin_ky = -3\nmax_ky = 3\n"
                        "min_kz = 1\nmax_kz = 3",
        "output1": "data_format = %24.16e",
        "output2": "single_file_per_rank = false",
    })

    def run(directory, *arguments):
        directory.mkdir()
        result = subprocess.run(
            ["mpirun", "-np", "2", str(ATHENA), "-d", str(directory), *arguments],
            capture_output=True, text=True, check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    reference, split, resumed = [tmp_path / name for name in ("reference", "split", "resumed")]
    run(reference, "-i", str(input_file))
    run(split, "-i", str(input_file), "time/nlim=3", "output2/file_type=rst",
        "output2/dcycle=3", "output2/id=restart",
        f"output2/single_file_per_rank={str(single_file_per_rank).lower()}")
    restart_dir = split / "rst"
    if single_file_per_rank:
        restart_dir /= "rank_00000000"
        assert list((split / "rst" / "rank_00000001").glob("*.rst"))
    restart = max(restart_dir.glob("*.rst"), key=lambda path: path.stat().st_mtime)
    run(resumed, "-r", str(restart), "time/nlim=6", "output2/file_type=bin",
        "output2/variable=turb_force", "output2/id=resume_force", "output2/dcycle=6",
        "output2/single_file_per_rank=false")
    assert_equal_force_outputs(latest_force(reference), latest_force(resumed))
    reference_history = np.loadtxt(next(reference.glob("*.hydro.hst")))[-1]
    resumed_history = np.loadtxt(next(resumed.glob("*.hydro.hst")))[-1]
    np.testing.assert_array_equal(reference_history, resumed_history)
