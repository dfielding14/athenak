"""LF shared-face conservation with smooth density and magnetic field."""
import pytest

import test_suite.testutils as testutils
from test_suite.cgl import cgl_lf_smr_conservation as conservation
from test_suite.cgl import test_cgl_amr_mpi_gpu as amr

pytestmark = amr.pytestmark


@pytest.mark.parametrize("dimension,primitive,integrator", conservation.CASES)
def test_cgl_lf_smr_conservation_mpi_gpu(tmp_path, dimension, primitive, integrator):
    states = []
    for ranks in (1, 4):
        directory = tmp_path / f"ranks{ranks}"
        directory.mkdir()
        flags = conservation.flags(directory, dimension, primitive, integrator)
        amr._launch(1, ranks, conservation.INPUT_NAME, conservation.BASENAME, *flags)
        states.append(conservation.check(directory, dimension))
    conservation.assert_same(*states)
