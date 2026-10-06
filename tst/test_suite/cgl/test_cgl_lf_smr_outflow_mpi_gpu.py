"""Opt-in MPI+GPU coverage for physical/coarse-fine LF corner ghosts."""
import numpy as np
import pytest

from test_suite.cgl import cgl_lf_outflow as uniform
from test_suite.cgl import test_cgl_amr_mpi_gpu as amr

pytestmark = amr.pytestmark


@pytest.mark.parametrize("nlim", [0, 4])
def test_cgl_lf_smr_outflow_mpi_gpu(tmp_path, nlim):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        amr._launch(1, ranks, "cgl_lf_smr_outflow_uniform.athinput",
                    uniform.BASENAME, *uniform.flags(directory, nlim))
        states.append(uniform.check(directory, nlim))
    for variable in states[0]["var_names"]:
        np.testing.assert_array_equal(states[0]["mb_data"][variable],
                                      states[1]["mb_data"][variable])
