"""MPI coverage of current A/mu boundary representation and corner ownership."""
import pytest
from test_suite.cgl import cgl_lf_boundary as boundary
from test_suite.cgl import test_cgl_amr_mpi_gpu as amr

pytestmark = amr.pytestmark


@pytest.mark.parametrize("integrator", ["sts", "explicit"])
@pytest.mark.parametrize("geometry", ["smr", "3d"])
def test_cgl_lf_uniform_boundary_mpi_gpu(tmp_path, integrator, geometry):
    for kind in ["inflow", "user"]:
        final = []
        for ranks in [1, 4]:
            directory = tmp_path / f"{kind}-{ranks}"
            directory.mkdir()
            flags = boundary.flags(directory, integrator, geometry, kind)
            amr._launch(1, ranks, boundary.INPUT_NAME, boundary.BASENAME, *flags)
            final.append(boundary.check_uniform(directory))
        boundary.same_state(*final)
