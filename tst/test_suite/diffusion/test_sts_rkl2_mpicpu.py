"""RKL2 MPI agreement, AMR, Ohmic heating, and full-state restart checks."""

import numpy as np
import pytest

from test_suite.diffusion.test_sts_rkl2_cpu import (
    cells, flags_for, read_state, run_case,
)

AMR = """
<mesh_refinement>
refinement = adaptive
num_levels = 2
ncycle_check = 1
refinement_interval = 1
max_nmb_per_rank = 128
<amr_criterion1>
method = location
location_x1 = 0
location_x2 = 0
location_rad = 0.4
"""


def compare_states(left, right):
    assert left["cycle"] == right["cycle"]
    assert left["time"] == right["time"]
    xyz1, volume1, fields1 = cells(left)
    xyz2, volume2, fields2 = cells(right)
    np.testing.assert_array_equal(xyz1, xyz2)
    np.testing.assert_array_equal(volume1, volume2)
    for field in fields1:
        np.testing.assert_allclose(fields1[field], fields2[field], rtol=2e-13, atol=2e-14)


@pytest.mark.parametrize("process", ["conduction", "ohmic"])
def test_sts_amr_agrees_on_one_and_two_ranks(tmp_path, process):
    fluid, flags = flags_for(process)
    flags += ["mesh/nx1=32", "mesh/nx2=16", "meshblock/nx1=8", "meshblock/nx2=8",
              "mesh/x2min=-2", "mesh/x2max=2", "time/nlim=4", "time/tlim=10"]
    if process == "conduction":
        flags += ["problem/shock_dir=2"]
    outputs = []
    for ranks in (1, 2):
        directory = tmp_path / str(ranks)
        run_case(directory, fluid, flags, ranks=ranks, extra=AMR)
        initial = read_state(directory, first=True)
        final = read_state(directory)
        assert final["n_mbs"] > initial["n_mbs"]
        assert len(np.unique(final["mb_logical"][:, 3])) > 1
        _, ivol, istate = cells(initial)
        _, fvol, fstate = cells(final)
        for field in ("dens", "ener"):
            np.testing.assert_allclose(np.dot(fvol, fstate[field]),
                                       np.dot(ivol, istate[field]),
                                       rtol=2e-12, atol=2e-12)
        if process == "ohmic":
            emag0 = np.dot(ivol, sum(istate[f"bcc{a}"]**2 for a in (1, 2, 3))) / 2
            emag1 = np.dot(fvol, sum(fstate[f"bcc{a}"]**2 for a in (1, 2, 3))) / 2
            assert emag1 < 0.99 * emag0
            ekin0 = np.dot(ivol, sum(istate[f"mom{a}"]**2 for a in (1, 2, 3))
                           / istate["dens"]) / 2
            ekin1 = np.dot(fvol, sum(fstate[f"mom{a}"]**2 for a in (1, 2, 3))
                           / fstate["dens"]) / 2
            eint0 = np.dot(ivol, istate["ener"]) - emag0 - ekin0
            eint1 = np.dot(fvol, fstate["ener"]) - emag1 - ekin1
            assert eint1 > eint0
        outputs.append(final)
    compare_states(*outputs)


def test_sts_cycle_boundary_restart(tmp_path):
    fluid, flags = flags_for("ohmic")
    common = flags + ["time/tlim=10", "output2/dcycle=2"]
    restart_output = "\n<output2>\nfile_type=rst\ndcycle=2\n"
    direct = tmp_path / "direct"
    split = tmp_path / "split"
    resumed = tmp_path / "resumed"
    run_case(direct, fluid, common + ["time/nlim=4"], ranks=2, extra=restart_output)
    run_case(split, fluid, common + ["time/nlim=2"], ranks=2, extra=restart_output)
    checkpoints = sorted((split / "rst").rglob("*.rst"))
    checkpoints = [p for p in checkpoints if p.parent.name in ("rst", "rank_00000000")]
    assert checkpoints
    run_case(resumed, fluid, ["time/nlim=4", "output2/dcycle=0"], ranks=2,
             restart=checkpoints[-1])
    compare_states(read_state(direct), read_state(resumed))
