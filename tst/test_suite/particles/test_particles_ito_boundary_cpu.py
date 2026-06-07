"""Exact periodic-boundary regression tests for Ito flux tracers."""

from pathlib import Path
import shutil
import struct
import sys

import numpy as np

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "tst/inputs/particles_ito_tags.athinput")
RUN_ROOT = Path("run_particles_ito_boundary_cpu")
PARTICLE_HEADER = struct.Struct("@16s6iQ")


def _set_first_particle_x1(restart_path, x1):
    payload = bytearray(restart_path.read_bytes())
    particle_offset = payload.index(b"ATHKPRTCLMC")
    fields = PARTICLE_HEADER.unpack_from(payload, particle_offset)
    _, version, enabled, nrdata, _, nlocal, nschedules, _ = fields
    assert version == 3
    assert enabled == 2
    assert nrdata >= 1
    assert nlocal == 1

    real_size = 8
    data_offset = (
        particle_offset
        + PARTICLE_HEADER.size
        + nschedules * real_size
        + 3 * nschedules * struct.calcsize("@i")
    )
    struct.pack_into("@d", payload, data_offset, x1)
    restart_path.write_bytes(payload)


def _last_position(run_dir):
    path = (
        run_dir
        / "prtcl_thermo_history"
        / "ito_tags.prtcl_thermo_history.thp"
    )
    data = read_history(path)
    final = data["cycle"] == np.max(data["cycle"])
    assert np.count_nonzero(final) >= 1
    x1 = data["x1"][final]
    x2 = data["x2"][final]
    np.testing.assert_array_equal(x1, np.full_like(x1, x1[-1]))
    np.testing.assert_array_equal(x2, np.full_like(x2, x2[-1]))
    return np.array([x1[-1], x2[-1]])


def test_exact_periodic_upper_boundary_wraps_before_following_step():
    """A tracer at xmax uses periodic CIC ghosts, wraps to xmin, and remains valid."""
    initial = RUN_ROOT / "initial"
    first = RUN_ROOT / "first"
    second = RUN_ROOT / "second"
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        RUN_ROOT.mkdir(parents=True)
        assert testutils.run(
            INPUT,
            [
                "-d",
                str(initial),
                "time/nlim=0",
                "problem/vx0=0.0",
                "tracer_seed1/count_per_event=1",
                "output3/dt=-1.0",
            ],
        )
        restart = sorted(
            (initial / "rst/rank_00000000").glob("ito_tags.*.rst")
        )[-1]
        initial_y = _last_position(initial)[1]
        _set_first_particle_x1(restart, 1.0)

        assert testutils.run_command(
            [
                "./athena",
                "-r",
                str(restart),
                "-d",
                str(first),
                "time/nlim=1",
                "output3/dt=-1.0",
            ]
        )
        np.testing.assert_array_equal(_last_position(first), np.array([0.0, initial_y]))

        continued_restart = sorted(
            (first / "rst/rank_00000000").glob("ito_tags.*.rst")
        )[-1]
        assert testutils.run_command(
            [
                "./athena",
                "-r",
                str(continued_restart),
                "-d",
                str(second),
                "time/nlim=2",
                "output3/dt=-1.0",
            ]
        )
        np.testing.assert_array_equal(_last_position(second), np.array([0.0, initial_y]))
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)
