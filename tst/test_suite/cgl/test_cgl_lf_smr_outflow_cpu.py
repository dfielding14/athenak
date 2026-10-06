"""Physical/coarse-fine corner ghosts must be valid before the first LF stage."""
import pytest

import test_suite.testutils as testutils
from test_suite.cgl import cgl_lf_outflow as uniform


@pytest.mark.parametrize("nlim", [0, 4])
def test_cgl_lf_smr_outflow_cpu(tmp_path, nlim):
    testutils.run(uniform.INPUT, uniform.flags(tmp_path, nlim))
    uniform.check(tmp_path, nlim)
