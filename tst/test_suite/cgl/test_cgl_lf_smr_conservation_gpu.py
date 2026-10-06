"""LF shared-face conservation with smooth density and magnetic field."""
import pytest

import test_suite.testutils as testutils
from test_suite.cgl import cgl_lf_smr_conservation as conservation


@pytest.mark.parametrize("dimension,primitive,integrator", conservation.CASES)
def test_cgl_lf_smr_conservation_gpu(tmp_path, dimension, primitive, integrator):
    testutils.run(conservation.INPUT,
                  conservation.flags(tmp_path, dimension, primitive, integrator))
    conservation.check(tmp_path, dimension)
