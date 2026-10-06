"""Passive CGL physical and exact-flow acceptance on gpu."""
import os

import pytest
from test_suite.cgl import passive_acceptance

@pytest.mark.parametrize("suite", ["identity", "heating", "linear", "advection", "restart", "fences"])
def test_cgl_passive_gpu(tmp_path, suite):
    launcher = os.environ.get("ATHENAK_CGL_PASSIVE_LAUNCHER", "")
    passive_acceptance.run_validation(tmp_path, suite, launcher)
