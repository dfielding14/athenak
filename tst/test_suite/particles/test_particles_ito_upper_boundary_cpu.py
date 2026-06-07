"""CPU regression for Ito tracers exactly on periodic upper boundaries."""

from pathlib import Path
import shutil

import test_suite.testutils as testutils

from test_suite.particles.ito_upper_boundary_test_utils import (
    INPUT,
    assert_wrapped_state,
    final_history,
    latest_rank_zero_restart,
    patch_exact_upper_coordinates,
)


RUN_ROOT = Path("run_particles_ito_upper_boundary_cpu")


def test_ito_exact_upper_faces_corner_restart_and_second_push():
    """Exact upper faces/corner wrap, reassign GIDs, restart, and push again."""
    initial = RUN_ROOT / "initial"
    first = RUN_ROOT / "first"
    second = RUN_ROOT / "second"
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        RUN_ROOT.mkdir(parents=True)
        assert testutils.run(INPUT, ["-d", str(initial)])
        initial_restart, initial_rank_files = latest_rank_zero_restart(initial)
        _, initial_gid = patch_exact_upper_coordinates(initial_rank_files)

        assert testutils.run_command(
            [
                "./athena",
                "-r",
                str(initial_restart),
                "-d",
                str(first),
                "time/nlim=1",
                "time/tlim=0.05",
            ]
        )
        cycle, first_state = final_history(first)
        assert cycle == 1
        assert_wrapped_state(first_state, initial_gid)

        first_restart, _ = latest_rank_zero_restart(first)
        assert testutils.run_command(
            [
                "./athena",
                "-r",
                str(first_restart),
                "-d",
                str(second),
                "time/nlim=2",
                "time/tlim=0.1",
            ]
        )
        cycle, second_state = final_history(second)
        assert cycle == 2
        assert_wrapped_state(second_state, initial_gid)
        assert second_state == first_state
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)
