"""AMR refresh preserves finite anisotropy without an extra collision projection."""

import os
from pathlib import Path
import subprocess

import numpy as np

from test_suite.cgl import test_cgl_amr_gpu as amr


def test_cgl_amr_refresh_preserves_finite_anisotropy(tmp_path):
    source = Path(__file__).resolve().parents[3] / "inputs/tests"
    # End at the time-triggered derefinement, independently of the CFL timestep.
    switch_time = 2.7e-4
    result = subprocess.run(
        [str(Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()), "-i",
         str(source / "cgl_lf_amr_3d_current.athinput"), "-d", str(tmp_path),
         "time/nlim=-1", f"time/tlim={switch_time:.17g}",
         f"problem/current_threshold_switch_time={switch_time:.17g}",
         "output2/variable=mhd_w_bcc", "output2/id=state"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    states = [amr.bin_convert.read_binary(str(path))
              for path in sorted((tmp_path / "bin").glob("*.state.*.bin"))]
    assert states[-2]["cycle"] == states[-1]["cycle"]
    assert states[-2]["time"] == states[-1]["time"] == switch_time
    assert states[-2]["n_mbs"] == 216
    assert states[-1]["n_mbs"] == 27

    def margin(state):
        fields = state["mb_data"]
        return (np.asarray(fields["p_perp"], dtype=float)
                - np.asarray(fields["eint"], dtype=float)
                + sum(np.asarray(fields[name], dtype=float)**2
                      for name in ("bcc1", "bcc2", "bcc3")))

    assert np.min(margin(states[0])) > 0.03
    assert np.min(margin(states[-2])) > 0.05
    # Derefinement leaves positive pressures outside the firehose boundary.
    # A zero-time AMR refresh must not add an unrequested infinite-rate collision.
    tol = 2.0e-6
    final_margin = margin(states[-1])
    assert np.min(final_margin) < -tol
    for state in states:
        for name in ("dens", "eint", "p_perp"):
            values = np.asarray(state["mb_data"][name], dtype=float)
            assert np.all(np.isfinite(values)) and np.all(values > 0.0)
    history = amr.testutils.athena_read.hst(
        str(tmp_path / "cgl_lf_amr_3d_current.mhd.hst"))
    for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos"):
        assert np.all(np.asarray(history[name]) == 0.0)
    amr._assert_conserved(history, ("mass", "1-mom", "2-mom", "3-mom", "tot-E"))
