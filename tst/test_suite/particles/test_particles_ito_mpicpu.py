"""MPI and AMR regression tests for second-moment Ito mass-flux tracers."""

from pathlib import Path
import shutil
import sys

import numpy as np

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


UNIFORM_INPUT = str(ROOT / "inputs/particles/ito_tracers.athinput")
AMR_INPUT = str(ROOT / "inputs/particles/ito_tracers_amr.athinput")
RUN_UNIFORM = Path("run_particles_ito_mpi")
RUN_AMR = Path("run_particles_ito_amr_mpi")


def test_ito2_mpi_and_amr():
    """Ito-2 migrates across ranks and preserves continuous positions through AMR."""
    shutil.rmtree(RUN_UNIFORM, ignore_errors=True)
    shutil.rmtree(RUN_AMR, ignore_errors=True)
    try:
        assert testutils.mpi_run(UNIFORM_INPUT, ["-d", str(RUN_UNIFORM)], threads=2)
        assert testutils.mpi_run(AMR_INPUT, ["-d", str(RUN_AMR)], threads=2)

        history = RUN_AMR / (
            "prtcl_thermo_history/ito_tracers_amr.prtcl_thermo_history.thp"
        )
        assert history.exists()
        assert history.stat().st_size > 0
        data = read_history(history)
        moved = data["cycle"] > 0
        assert np.any(moved)

        # A continuous Ito trajectory must not be snapped back to AMR cell centers.
        distance_to_center = np.full(np.count_nonzero(moved), np.inf)
        for dx in (1.0 / 32.0, 1.0 / 64.0):
            cell_coordinate = (data["x1"][moved] + 0.5) / dx - 0.5
            distance_to_center = np.minimum(
                distance_to_center, np.abs(cell_coordinate - np.rint(cell_coordinate))
            )
        assert np.any(distance_to_center > 1.0e-8)

        restart_files = sorted((RUN_AMR / "rst").glob("rank_*/*.rst"))
        assert restart_files
        sizes = np.array([path.stat().st_size for path in restart_files])
        assert np.all(sizes > 0)
        latest_rank0_restart = sorted(
            (RUN_AMR / "rst/rank_00000000").glob("*.rst")
        )[-1]
        assert testutils.run_command(
            [
                "mpirun",
                "-np",
                "2",
                "./athena",
                "-r",
                str(latest_rank0_restart),
                "-d",
                str(RUN_AMR / "restart"),
                "time/nlim=6",
                "time/tlim=0.035",
            ]
        )
    finally:
        shutil.rmtree(RUN_UNIFORM, ignore_errors=True)
        shutil.rmtree(RUN_AMR, ignore_errors=True)
