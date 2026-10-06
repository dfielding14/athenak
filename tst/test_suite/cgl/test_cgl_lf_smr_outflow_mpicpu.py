"""The uniform boundary/refinement state must survive MPI decomposition."""
import numpy as np
import pytest

import test_suite.testutils as testutils
from test_suite.cgl import cgl_lf_outflow as uniform


@pytest.mark.parametrize("nlim", [0, 4])
def test_cgl_lf_smr_outflow_mpicpu(tmp_path, nlim):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        testutils.mpi_run(uniform.INPUT, uniform.flags(directory, nlim), threads=ranks)
        states.append(uniform.check(directory, nlim))
    for variable in states[0]["var_names"]:
        np.testing.assert_array_equal(states[0]["mb_data"][variable],
                                      states[1]["mb_data"][variable])
