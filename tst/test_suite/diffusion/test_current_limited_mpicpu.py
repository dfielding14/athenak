"""Nonlinear resistivity on changing meshes and across cycle-boundary restarts."""

import numpy as np
import pytest

from test_suite.diffusion.test_current_limited_cpu import energy_budget, run_model
from test_suite.diffusion.test_sts_rkl2_cpu import read_state
from test_suite.diffusion.test_sts_rkl2_mpicpu import AMR, compare_states
from athena_read import hst


@pytest.mark.parametrize("guide", [0, 1])
def test_harris_boundary_reference_history(tmp_path, guide):
    # With no perturbation the open sheet has the same nonzero resistive Ez at
    # its center and boundary, hence zero reference-subtracted flux-change rate.
    flags = ["problem/test=harris", "problem/user_hist=true",
             "problem/perturbation_flux=0", f"problem/guide_field={guide}",
             "mesh/nx1=32", "mesh/nx2=32", "meshblock/nx2=16", "time/nlim=0",
             "mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow",
             "mesh/ix2_bc=outflow", "mesh/ox2_bc=outflow",
             "mhd/d_i=0.01", "mhd/ohmic_resistivity=0.001", "mhd/eta_max=0.01"]
    extra = "\n<output2>\nfile_type=hst\ndata_format=%24.16e\ndt=1\n"
    dy, width = 1/32, 0.1
    jz = -2*width*np.log(np.cosh(dy/width))/dy**2
    rho_edge = 1 + 1/np.cosh(0.5*dy/width)**2
    q = abs(jz)*0.01/np.sqrt(rho_edge)
    expected = 0.001/(1-q)*jz
    histories = []
    for ranks in (1, 2):
        directory = tmp_path / str(ranks)
        run_model(directory, flags, ranks=ranks, extra=extra)
        data = hst(str(directory / "CurrentLimited.user.hst"), raw=True)
        for label in ("x_etaJz", "x_Ez", "ref_etaJz", "ref_Ez"):
            np.testing.assert_allclose(data[label], expected, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(data["x_Ez"]-data["ref_Ez"], 0, atol=2e-14)
        histories.append(data)
    for label in histories[0]:
        np.testing.assert_allclose(histories[0][label], histories[1][label],
                                   rtol=2e-13, atol=2e-14)


def test_current_limited_amr_rank_agreement(tmp_path):
    flags = ["mesh/nx1=32", "mesh/nx2=16", "meshblock/nx1=8", "meshblock/nx2=8",
             "problem/wave_n2=1", "time/nlim=4", "time/tlim=10"]
    outputs = []
    for ranks in (1, 2):
        directory = tmp_path / str(ranks)
        run_model(directory, flags, mode="sts", ranks=ranks, extra=AMR)
        initial = read_state(directory, first=True)
        final = read_state(directory)
        assert final["n_mbs"] > initial["n_mbs"]
        np.testing.assert_allclose(energy_budget(final)[0], energy_budget(initial)[0],
                                   rtol=2e-12, atol=2e-12)
        outputs.append(final)
    compare_states(*outputs)


def test_current_limited_cycle_boundary_restart(tmp_path):
    common = ["time/tlim=10", "output2/dcycle=2"]
    extra = "\n<output2>\nfile_type=rst\ndcycle=2\n"
    direct, split, resumed = (tmp_path / name for name in ("direct", "split", "resumed"))
    run_model(direct, common + ["time/nlim=4"], mode="sts", ranks=2, extra=extra)
    run_model(split, common + ["time/nlim=2"], mode="sts", ranks=2, extra=extra)
    checkpoints = sorted((split / "rst").rglob("*.rst"))
    checkpoints = [p for p in checkpoints if p.parent.name in ("rst", "rank_00000000")]
    assert checkpoints
    run_model(resumed, ["time/nlim=4", "output2/dcycle=0"], mode="sts", ranks=2,
              restart=checkpoints[-1])
    compare_states(read_state(direct), read_state(resumed))
