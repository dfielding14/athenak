"""Passive CGL physical and exact-flow acceptance on mpicpu."""
import os

import pytest
from test_suite.cgl import passive_acceptance

@pytest.mark.parametrize("suite", ["identity", "heating", "linear", "advection", "restart", "fences"])
def test_cgl_passive_mpicpu(tmp_path, suite):
    launcher = os.environ.get("ATHENAK_CGL_PASSIVE_MPI_LAUNCHER", "mpirun -n 4")
    passive_acceptance.run_validation(tmp_path, suite, launcher)
