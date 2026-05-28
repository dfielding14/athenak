"""Regression tests for fixed-grid and refined projection output."""

from pathlib import Path
import shutil

import numpy as np
import projection
import pytest
import test_suite.testutils as testutils


def _remove_projection_outputs():
    path = Path("proj")
    if path.exists():
        shutil.rmtree(path)


def test_fixed_grid_projection_weighting_bounds_and_derived_variable():
    try:
        assert testutils.run("inputs/projection.athinput")

        bounded = projection.read_projection(
            "proj/projection.bounded_integral.00000.proj"
        )
        volume = projection.read_projection("proj/projection.volume_mean.00000.proj")
        mass = projection.read_projection("proj/projection.mass_mean.00000.proj")
        vorticity = projection.read_projection(
            "proj/projection.derived_vorticity.00000.proj"
        )
        axis_x1 = projection.read_projection("proj/projection.axis_x1.00000.proj")
        narrow_x1 = projection.read_projection(
            "proj/projection.narrow_x1_slab.00000.proj"
        )
        profile = projection.read_projection("proj/projection.bounded_profile.00000.proj")
        stddev = projection.read_projection("proj/projection.volume_stddev.00000.proj")
        mass_stddev = projection.read_projection("proj/projection.mass_stddev.00000.proj")
        uniform = projection.read_projection(
            "proj/projection.uniform_reference.00000.proj"
        )

        assert bounded["metadata"]["layout"] == "native_amr"
        assert bounded["patches"]
        assert np.allclose(bounded["fields"]["dens"], 0.5, atol=1.0e-13)
        assert np.allclose(uniform["fields"]["dens"], bounded["fields"]["dens"],
                           atol=1.0e-13)
        assert np.allclose(volume["fields"]["dens"], 1.0, atol=1.0e-13)
        assert np.allclose(mass["fields"]["dens"], 1.125, atol=1.0e-13)
        assert np.allclose(vorticity["fields"]["vor2"], 0.0, atol=1.0e-13)

        z = axis_x1["y"]
        expected_z = 1.0 + 0.5 * np.sin(2.0 * np.pi * (z + 0.5))
        assert np.allclose(axis_x1["fields"]["dens"], expected_z[:, None], atol=1.0e-13)
        assert np.allclose(narrow_x1["fields"]["dens"], 0.1 * expected_z[:, None],
                           atol=1.0e-13)
        assert profile["fields"]["dens"].shape == (1, 8)
        assert profile["metadata"]["projection_axes"] == "x1,x2"
        assert np.allclose(profile["fields"]["dens"][0], 0.05 * expected_z, atol=1.0e-13)
        assert np.allclose(stddev["fields"]["dens"], np.sqrt(0.125), atol=1.0e-13)
        assert np.allclose(mass_stddev["fields"]["dens"], np.sqrt(0.109375),
                           atol=1.0e-13)
    finally:
        _remove_projection_outputs()
        testutils.cleanup()


def test_uniform_projection_rejects_working_set_over_memory_guard():
    try:
        with pytest.raises(RuntimeError):
            testutils.run(
                "inputs/projection.athinput",
                flags=["output10/projection_max_megabytes=0"],
            )
    finally:
        _remove_projection_outputs()
        testutils.cleanup()


def test_uniform_projection_rejects_full_domain_request():
    try:
        with pytest.raises(RuntimeError):
            testutils.run(
                "inputs/projection.athinput",
                flags=[
                    "output10/projection_min=-0.5",
                    "output10/projection_max=0.5",
                ],
            )
    finally:
        _remove_projection_outputs()
        testutils.cleanup()


def test_projection_preserves_uniform_field_on_mixed_static_refinement():
    try:
        assert testutils.run("inputs/projection_smr.athinput")
        output = projection.read_projection("proj/projection_smr.mixed_levels.00000.proj")
        assert output["fields"]["dens"].shape == (32, 32)
        assert {patch["level"] for patch in output["patches"]} == {0, 1}
        assert np.allclose(output["fields"]["dens"], 0.5, atol=1.0e-13)
    finally:
        _remove_projection_outputs()
        testutils.cleanup()


def test_projection_preserves_uniform_field_after_adaptive_refinement():
    try:
        assert testutils.run("inputs/projection_amr.athinput")
        output = projection.read_projection("proj/projection_amr.adaptive.00001.proj")
        profile = projection.read_projection(
            "proj/projection_amr.adaptive_profile.00001.proj"
        )
        coarse = projection.read_projection(
            "proj/projection_amr.adaptive.00001.proj", level=0
        )
        assert output["fields"]["dens"].shape == (16, 16)
        assert np.allclose(output["fields"]["dens"], 0.5, atol=1.0e-13)
        assert coarse["fields"]["dens"].shape == (8, 8)
        assert np.allclose(coarse["fields"]["dens"], 0.5, atol=1.0e-13)
        assert profile["fields"]["dens"].shape == (1, 16)
        assert np.allclose(profile["fields"]["dens"], 0.05, atol=1.0e-13)
    finally:
        _remove_projection_outputs()
        testutils.cleanup()
