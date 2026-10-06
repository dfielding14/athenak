"""Representation-aware CGL-LF boundaries: analytic state and independent callback."""
import pytest
import test_suite.testutils as testutils
from test_suite.cgl import cgl_lf_boundary as boundary


@pytest.mark.parametrize("integrator", ["sts", "explicit"])
@pytest.mark.parametrize("geometry", ["1d", "mixed2d", "smr", "3d"])
def test_cgl_lf_uniform_boundary_cpu(tmp_path, integrator, geometry):
    for kind in ["inflow", "user"]:
        directory = tmp_path / kind
        directory.mkdir()
        testutils.run(boundary.INPUT,
                      boundary.flags(directory, integrator, geometry, kind))
        boundary.check_uniform(directory)


@pytest.mark.parametrize("integrator", ["sts", "explicit"])
def test_cgl_lf_user_matches_inflow_cpu(tmp_path, integrator):
    final = []
    for kind in ["inflow", "user"]:
        directory = tmp_path / kind
        directory.mkdir()
        testutils.run(boundary.INPUT,
                      boundary.flags(directory, integrator, "1d", kind, 0.01)
                      + ["problem/pperp0=1.0"])
        final.append(boundary.states(directory)[-1])
    boundary.same_state(*final)
