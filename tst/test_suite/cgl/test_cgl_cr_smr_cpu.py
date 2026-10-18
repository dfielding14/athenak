"""Frozen CGL particles crossing an actual static coarse/fine interface."""

import numpy as np
import pytest

from test_suite.cgl.test_cgl_cr_cpu import (
    TEMPLATE, checkpoint, frozen, particles, run, state,
)


def test_frozen_cgl_boris_crosses_static_refinement(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    text = TEMPLATE.read_text().replace("nx1 = 24", "nx1 = 16")
    text = text.replace("nx2 = 24", "nx2 = 8").replace("nx3 = 1", "nx3 = 8")
    text = text.replace("refinement = adaptive", "refinement = static")
    text = text.replace("variable = turb_force", "variable = mhd_w_bcc")
    text = text.replace("id = force", "id = state")
    text += ("\n<refined_region1>\nlevel = 1\nx1min = 0.51\nx1max = 0.99\n"
             "x2min = 0.1\nx2max = 0.9\nx3min = 0.1\nx3max = 0.9\n")
    deck = source / "background.athinput"
    deck.write_text(text)
    run(source, "-i", deck, "time/nlim=0", "time/tlim=0.1",
        "problem/field_amp=0", "problem/guide_b1=2", "problem/guide_b2=0",
        "problem/guide_b3=0")
    original = state(source)
    assert original["n_mbs"] == 9
    levels = original["mb_logical"][:, 3]
    assert len(np.unique(levels)) == 2

    flags = ("particles/ppc=0.001953125", "problem/particle_position=fixed",
             "problem/particle_x1=0.49", "problem/particle_x2=0.3",
             "problem/particle_x3=0.3", "problem/v0x=1", "problem/v0y=0",
             "problem/v0z=0", "output98/nparticles=1")
    loaded = tmp_path / "loaded"
    advanced = tmp_path / "advanced"
    frozen(loaded, checkpoint(source), 0.08, *flags, "time/nlim=0")
    frozen(advanced, checkpoint(source), 0.08, *flags)
    before = state(loaded)
    after = state(advanced)
    np.testing.assert_array_equal(after["mb_logical"], original["mb_logical"])
    np.testing.assert_array_equal(after["mb_geometry"], original["mb_geometry"])
    for name, values in before["mb_data"].items():
        np.testing.assert_array_equal(after["mb_data"][name], values)
    # Native restart payload: six CGL variables and three face fields, including ghosts.
    block_values = 6 * 12**3 + 3 * 13 * 12**2
    payload_bytes = 9 * block_values * np.dtype(np.float64).itemsize
    assert (checkpoint(advanced).read_bytes()[-payload_bytes:] ==
            checkpoint(loaded).read_bytes()[-payload_bytes:])

    _, initial_particles = particles(loaded)
    time, final_particles = particles(advanced)
    assert len(initial_particles) == len(final_particles) == 9
    assert time == pytest.approx(0.08)
    np.testing.assert_array_equal(levels[initial_particles[:, 0].astype(int)], levels.min())
    np.testing.assert_array_equal(levels[final_particles[:, 0].astype(int)], levels.max())
    np.testing.assert_array_equal(final_particles[:, 1:3], initial_particles[:, 1:3])
    np.testing.assert_allclose(final_particles[:, 3:6], [[0.57, 0.3, 0.3]] * 9,
                               rtol=0, atol=2e-13)
    np.testing.assert_allclose(final_particles[:, 6:9], [[1, 0, 0]] * 9, atol=2e-13)
    np.testing.assert_allclose(final_particles[:, 10:13], [[2, 0, 0]] * 9, atol=2e-13)
