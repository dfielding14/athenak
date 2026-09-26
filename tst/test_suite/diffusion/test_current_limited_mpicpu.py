"""Nonlinear resistivity on changing meshes and across cycle-boundary restarts."""

import numpy as np

from test_suite.diffusion.test_current_limited_cpu import energy_budget, run_model
from test_suite.diffusion.test_sts_rkl2_cpu import read_state
from test_suite.diffusion.test_sts_rkl2_mpicpu import AMR, compare_states


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
