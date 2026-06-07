"""Timestep safety regressions for Ito-2 flux tracers."""

from pathlib import Path
import shutil
import sys

import numpy as np

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "inputs/particles/ito_tracers.athinput")
RUN_DIR = Path("run_particles_ito_timestep_cpu")


def test_first_step_is_limited_before_ito_probabilities_exist():
    """A large configured CFL is conservatively limited on cycle zero."""
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    try:
        assert testutils.run(
            INPUT,
            [
                "-d",
                str(RUN_DIR),
                "job/basename=ito_timestep",
                "mesh/nx3=8",
                "meshblock/nx3=8",
                "time/cfl_number=0.9",
                "time/nlim=1",
                "time/tlim=1.0",
                "hydro/iso_sound_speed=0.01",
                "problem/vx0=1.0",
                "problem/vy0=0.15",
                "problem/vz0=-0.05",
                "particles/ito_probability_target=0.99",
                "tracer_seed1/count_per_event=64",
                "output2/dt=-1.0",
                "output3/dt=-1.0",
            ],
        )
        history = read_history(
            RUN_DIR / "prtcl_thermo_history/ito_timestep.prtcl_thermo_history.thp"
        )
        final_time = float(np.max(history["time"]))
        fluid_directional_dt = (1.0 / 32.0) / (1.0 + 0.01)
        expected_limit = 0.99 * fluid_directional_dt / 3.0
        # The history payload uses Athena's Real precision.  The limiter should
        # be exact to a few single-precision ulps in float builds.
        np.testing.assert_allclose(final_time, expected_limit, rtol=5.0e-7)
        assert final_time < 0.9 * fluid_directional_dt
    finally:
        shutil.rmtree(RUN_DIR, ignore_errors=True)
