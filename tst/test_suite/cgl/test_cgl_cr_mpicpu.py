"""Two-rank frozen-field particle ownership and crossing regression."""

import re

import numpy as np
import pytest

from test_suite.cgl.test_cgl_cr_cpu import (
    OVERLAY, background, checkpoint, check_tracks, frozen, latest, particles,
)


def test_frozen_cgl_boris_crosses_mpi_boundary(tmp_path, monkeypatch):
    overlay = tmp_path / "particles_vtk.athinput"
    overlay.write_text(OVERLAY.read_text() + "\n<output97>\nfile_type = pvtk\n"
                       "variable = prtcl_all\ndt = 0.01\n")
    monkeypatch.setattr("test_suite.cgl.test_cgl_cr_cpu.OVERLAY", overlay)
    restart = background(tmp_path / "source", uniform=True, along_x=True, ranks=2)
    initial = tmp_path / "initial"
    advanced = tmp_path / "advanced"
    flags = ("particles/ppc=0.001953125", "particles/check_consistency_mode=none",
             "particles/update_global_counts_each_exchange=false",
             "problem/particle_position=fixed",
             "problem/particle_x1=0.49", "problem/particle_x2=0.5",
             "problem/particle_x3=0.5", "problem/v0x=1", "problem/v0z=0",
             "output98/dt=0.01", "output98/nparticles=1")
    frozen(initial, restart, 0.08, *flags, "time/nlim=0", ranks=2)
    frozen(advanced, restart, 0.08, *flags, ranks=2)
    _, before = particles(initial)
    time, after = particles(advanced)
    assert len(before) == len(after) == 2
    assert time == pytest.approx(0.08)
    np.testing.assert_array_equal(before[:, 0], 0)
    np.testing.assert_array_equal(after[:, 0], 1)
    np.testing.assert_array_equal(before[:, 1:3], after[:, 1:3])
    np.testing.assert_allclose(after[:, 3], 0.57, atol=2e-13)
    np.testing.assert_allclose(after[:, 6:9], [[1, 0, 0], [1, 0, 0]], atol=2e-13)
    np.testing.assert_allclose(after[:, 10:13], [[2, 0, 0], [2, 0, 0]], atol=2e-13)
    continued = tmp_path / "continued"
    prst = latest(initial / "prst", "rank_00000000", "*.prst")
    frozen(continued, checkpoint(initial), 0.08, *flags,
           "particles/seed_on_restart=false", "problem/prtcl_rst_flag=1",
           f"problem/prtcl_res_file={prst}", ranks=2)
    _, result = particles(continued)
    np.testing.assert_allclose(result, after, rtol=2e-13, atol=2e-13)
    check_tracks(advanced, after, 1)
    check_tracks(continued, result, 1)
    for directory in (advanced, continued):
        data = latest(directory, "pvtk", "*.part.vtk").read_bytes()
        points = re.search(rb"\nPOINTS (\d+) float\n", data)
        assert points and int(points[1]) == 2
        xyz = np.frombuffer(data, dtype=">f4", count=6,
                            offset=points.end()).reshape(2, 3)
        np.testing.assert_allclose(xyz, [[0.57, 0.5, 0.5]] * 2, atol=1e-7)
        assert b"\nPOINT_DATA 2\n" in data
        species = re.search(rb"\nSCALARS species float\nLOOKUP_TABLE default\n", data)
        assert species
        np.testing.assert_array_equal(np.frombuffer(data, dtype=">f4", count=2,
                                                    offset=species.end()), [0, 0])
