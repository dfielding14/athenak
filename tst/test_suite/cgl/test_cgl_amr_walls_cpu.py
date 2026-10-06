"""Post-AMR CGL walls must precede the next LF sweep and preserve energy."""

from pathlib import Path
import subprocess

import numpy as np

from test_suite.cgl import test_cgl_amr_gpu as amr


def test_cgl_amr_walls_after_final_field_refresh(tmp_path):
    source = Path(__file__).resolve().parents[3] / "inputs/tests"
    result = subprocess.run(
        [str(Path("athena").resolve()), "-i",
         str(source / "cgl_lf_amr_3d_current.athinput"), "-d", str(tmp_path),
         # Restore the original cycle-three event using only a timestep cap.
         "time/sts_max_dt_ratio=0.13", "time/nlim=3",
         "output2/variable=mhd_w_bcc", "output2/id=state"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    states = [amr.bin_convert.read_binary(str(path))
              for path in sorted((tmp_path / "bin").glob("*.state.*.bin"))]
    assert states[-2]["cycle"] == states[-1]["cycle"] == 3
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
    # Binary field output is float32; these eight cells cross the wall by 0.0046
    # without the post-AMR projection. Check the state before any LF stage runs.
    tol = 2.0e-6
    final_margin = margin(states[-1])
    assert np.min(final_margin) >= -tol
    assert np.count_nonzero(np.abs(final_margin) <= tol) == 8
    history = amr.testutils.athena_read.hst(
        str(tmp_path / "cgl_lf_amr_3d_current.mhd.hst"))
    amr._assert_clean_lf(history)
    amr._assert_conserved(history, ("mass", "1-mom", "2-mom", "3-mom", "tot-E"))
