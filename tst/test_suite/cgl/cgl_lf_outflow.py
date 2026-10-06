"""Analytic LF state at a physical boundary meeting a coarse/fine interface."""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis/python"))
import bin_convert  # noqa: E402

INPUT = str(ROOT / "inputs/tests/cgl_lf_smr_outflow_uniform.athinput")
BASENAME = "cgl_lf_smr_outflow_uniform"


def flags(directory, nlim):
    return ["-d", str(directory), f"time/nlim={nlim}"]


def check(directory, nlim):
    """Check initialization and final outputs, including the LF stencil halo."""
    expected = {
        "dens": 1.0, "velx": 0.08, "vely": 0.05, "velz": 0.0,
        "eint": 1.0, "p_perp": 1.01, "s_00": 0.5,
        "bcc1": 0.7, "bcc2": 0.2, "bcc3": -0.15,
    }
    paths = sorted((directory / "bin").glob(f"{BASENAME}.state.*.bin"))
    assert paths, "uniform outflow regression produced no primitive states"
    states = [bin_convert.read_binary(str(path)) for path in paths]
    assert states[0]["cycle"] == 0
    assert states[-1]["cycle"] == nlim
    for state in states:
        assert state["n_mbs"] == 25, "fixture must contain the boundary/refinement corners"
        assert set(state["var_names"]) == set(expected)
        # Two output ghost layers are present. The nearest layer in each active
        # dimension is read by the transverse LF stencil, including its corners.
        for block in range(state["n_mbs"]):
            slices = tuple(slice(1, size + 3) if size > 1 else slice(None)
                           for size in (state["nx3_mb"], state["nx2_mb"],
                                        state["nx1_mb"]))
            for variable, value in expected.items():
                values = np.asarray(state["mb_data"][variable][block])[slices]
                assert np.all(np.isfinite(values)), variable
                # Field output is float32. Round the analytic target once and
                # allow a little output rounding, well below the baseline defect.
                np.testing.assert_allclose(
                    values, np.float32(value), rtol=0.0, atol=2.0e-7,
                    err_msg=(f"{variable}: cycle {state['cycle']}, block "
                             f"{state['mb_logical'][block]} active/ghost state"),
                )
    return states[-1]
