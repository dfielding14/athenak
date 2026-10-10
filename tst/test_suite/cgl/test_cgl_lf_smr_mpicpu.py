"""LF corner states and shared-face conservation on one and four MPI ranks."""

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl import cgl_lf_outflow as outflow
from test_suite.cgl import cgl_lf_smr_conservation as conservation
from test_suite.cgl import test_cgl_amr_gpu as amr


@pytest.mark.parametrize("nlim", [0, 4])
def test_cgl_lf_smr_outflow(tmp_path, nlim):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        testutils.mpi_run(outflow.INPUT, outflow.flags(directory, nlim), threads=ranks)
        states.append(outflow.check(directory, nlim))
    for variable in states[0]["var_names"]:
        np.testing.assert_array_equal(states[0]["mb_data"][variable],
                                      states[1]["mb_data"][variable])


def test_cgl_lf_smr_outflow_nonuniform(tmp_path):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        flags = outflow.flags(directory, 4) + ["problem/density_amp=0.2",
                 "problem/scalar_amp=0.1", "problem/field_amp=0.02"]
        testutils.mpi_run(outflow.INPUT, flags, threads=ranks)
        state = outflow.bin_convert.read_binary(str(sorted(
            (directory / "bin").glob(f"{outflow.BASENAME}.state.*.bin"))[-1]))
        assert state["cycle"] == 4 and state["n_mbs"] == 25
        for variable in state["var_names"]:
            assert np.all(np.isfinite(state["mb_data"][variable]))
        for variable in ("dens", "eint", "p_perp"):
            assert np.min(state["mb_data"][variable]) > 0
        amr._assert_clean_lf(testutils.athena_read.hst(
            str(directory / f"{outflow.BASENAME}.mhd.hst")))
        user = testutils.athena_read.hst(str(directory / f"{outflow.BASENAME}.user.hst"))
        amr._assert_clean_user(user, max_ndiv=1e-12)
        amr._assert_no_amr_repairs(user)
        states.append(state)
    amr._assert_same_state(*states)


@pytest.mark.parametrize("dimension,primitive,integrator", conservation.CASES)
def test_cgl_lf_smr_conservation(tmp_path, dimension, primitive, integrator):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        testutils.mpi_run(conservation.INPUT,
                         conservation.flags(directory, dimension, primitive, integrator),
                         threads=ranks)
        states.append(conservation.check(directory, dimension))
    conservation.assert_same(*states)
