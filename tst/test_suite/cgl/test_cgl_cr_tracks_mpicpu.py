"""Rich tracks follow their selected particles across MPI ownership changes."""

import numpy as np
import pytest

from test_suite.cgl.test_cgl_cr_cpu import OVERLAY, background, frozen, particles
from test_suite.cgl.test_cgl_cr_tracks_cpu import frames


def test_rich_tracks_follow_migrating_particles(tmp_path, monkeypatch):
    restart = background(tmp_path / "source", uniform=True, along_x=True, ranks=2)
    overlay = tmp_path / "rich.athinput"
    overlay.write_text(OVERLAY.read_text() +
                       "\n<output98>\nparticles_per_species = 1\nlast_time = 0\n")
    monkeypatch.setattr("test_suite.cgl.test_cgl_cr_cpu.OVERLAY", overlay)
    flags = ("particles/ppc=0.001953125", "particles/nspecies=2",
             "particles/mass_log_spacing=2", "problem/particle_position=fixed",
             "problem/particle_x1=0.49", "problem/particle_x2=0.5",
             "problem/particle_x3=0.5", "problem/v0x=1", "problem/v0z=0",
             "output98/dt=0.01", "output98/nparticles=1")
    initial = tmp_path / "initial"
    tracked, control = tmp_path / "tracked", tmp_path / "control"
    frozen(initial, restart, 0.08, *flags, "time/nlim=0", ranks=2)
    frozen(tracked, restart, 0.08, *flags, ranks=2)
    frozen(control, restart, 0.08, *flags, "output98/particles_per_species=0", ranks=2)
    time, rows = particles(tracked)
    assert time == pytest.approx(0.08)
    np.testing.assert_array_equal(rows, particles(control)[1])
    np.testing.assert_array_equal(rows[:, 0], 1)
    shards = [frames(path) for path in sorted((tracked / "trk").rglob("*.trk"))]
    assert len(shards) == 2
    assert [(t, c) for t, c, _ in shards[0]] == [(t, c) for t, c, _ in shards[1]]
    loaded = [frames(path) for path in sorted((initial / "trk").rglob("*.trk"))]
    assert loaded[0][0][0] == 0
    assert len(loaded[0][0][2]) == len(shards[1][-1][2]) == 2
    assert len(loaded[1][0][2]) == len(shards[0][-1][2]) == 0
    for left, right in zip(*shards):
        data = np.concatenate([left[2], right[2]])
        data = data[np.argsort(data[:, 0])]
        np.testing.assert_array_equal(data[:, 0], [0, 1])
        np.testing.assert_array_equal(data[:, 1], [0, 2])
        np.testing.assert_array_equal(data[:, 2], [0, 1])
        np.testing.assert_allclose(data[:, 10:13], [[2, 0, 0]]*2, atol=1e-7)
    np.testing.assert_array_equal(data[:, 4:10], rows[[0, 2], 3:9].astype(np.float32))
    np.testing.assert_array_equal(data[:, 13:17], rows[[0, 2], 13:17].astype(np.float32))
