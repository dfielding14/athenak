"""3D hydro/MHD SGS products survive restarts with a different MPI rank count."""

import numpy as np
import pytest

from .test_turb_sgs_3d_cpu import (INPUT, MHD_INPUT, check_sgs, latest,
                                   require_success, run_athena)


@pytest.mark.parametrize("mhd", [False, True])
def test_3d_sgs_restart_with_changed_rank_count(tmp_path, mhd):
    input_file = tmp_path / "restart.athinput"
    input_file.write_text((MHD_INPUT if mhd else INPUT).read_text() +
                         "\n<output3>\nfile_type = rst\ndcycle = 4\n")
    overrides = ["time/tlim=1", "turb_driving/sol_fraction=0.3",
                 "turb_driving/accel_rms=2", "output2/coarsen_factor=12",
                 "output1/dcycle=8", "output2/dcycle=8"]
    overrides += [f"mesh/nx{axis}=24" for axis in (1, 2, 3)]
    overrides += [f"meshblock/nx{axis}=12" for axis in (1, 2, 3)]
    reference, split, resumed = [tmp_path / name
                                 for name in ("reference", "split", "resumed")]
    for directory, cycles, ranks, restart in (
        (reference, 8, 4, False), (split, 4, 4, False), (resumed, 8, 2, True)
    ):
        require_success(run_athena(
            directory, *overrides, f"time/nlim={cycles}", input_file=input_file,
            launcher=("mpirun", "-np", str(ranks)),
            restart=latest(split / "rst", "*.rst") if restart else None,
        ))
    expected = check_sgs(reference, 12, mhd)
    actual = check_sgs(resumed, 12, mhd)
    for before, after in zip(expected, actual):
        assert before["n_mbs"] == after["n_mbs"] == 8
        assert before["time"] == pytest.approx(after["time"], rel=1.0e-14)
        np.testing.assert_array_equal(before["mb_logical"], after["mb_logical"])
        for name in before["var_names"]:
            np.testing.assert_array_equal(before["mb_data"][name], after["mb_data"][name])
