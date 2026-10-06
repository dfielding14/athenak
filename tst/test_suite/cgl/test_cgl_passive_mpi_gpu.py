"""Passive CGL physical and exact-flow acceptance on mpi_gpu."""
import os

import pytest
from test_suite.cgl import passive_acceptance

pytestmark = pytest.mark.skipif(
    os.environ.get("ATHENAK_RUN_MPI_GPU") != "1",
    reason="set ATHENAK_RUN_MPI_GPU=1 inside an MPI+GPU allocation",
)

@pytest.mark.parametrize("suite", ["identity", "heating", "linear", "advection", "restart", "fences"])
def test_cgl_passive_mpi_gpu(tmp_path, suite):
    launcher = os.environ["ATHENAK_MPI_GPU_LAUNCHER"] + " -N 1 -n 4 --ntasks-per-node 4"
    passive_acceptance.run_validation(tmp_path, suite, launcher)
