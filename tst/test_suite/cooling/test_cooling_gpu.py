"""GPU-targeted cooling smoke tests selected by run_test_suite.py --gpu."""

from pathlib import Path

import numpy as np
import pytest

import test_suite.testutils as testutils

import sys

sys.path.insert(0, "../vis/python")
import athena_read  # noqa: E402


HYDRO_INPUT = "inputs/cooling_test_hydro.athinput"


def _remove_outputs(basename):
    Path(f"{basename}.hydro.hst").unlink(missing_ok=True)


def _run_case(basename, flags):
    _remove_outputs(basename)
    testutils.run(HYDRO_INPUT, [f"job/basename={basename}", *flags])
    return athena_read.hst(f"{basename}.hydro.hst", raw=True)


def _first_evolved_row(data):
    for row in range(1, len(data["time"])):
        if data["time"][row] > data["time"][row - 1]:
            return row
    pytest.fail("history file did not contain an evolved row")


def test_gpu_table_interpolation_with_scalar_axis():
    data = _run_case(
        "cooling_gpu_table_3d",
        ["cooling/cooling_model=table",
         "cooling/cooling_table=inputs/tables/cooling_lambda_3d.tbl",
         "problem/scalar0=1.0"],
    )
    row = _first_evolved_row(data)
    assert np.isclose(data["cool_gross"][row], 1.1e-2, rtol=1.0e-10,
                      atol=1.0e-12)


def test_gpu_user_hook_path():
    data = _run_case(
        "cooling_gpu_user_hook",
        ["cooling/cooling_model=user",
         "cooling/heating_model=user",
         "problem/register_user_cooling=true",
         "problem/user_lambda=1.0e-2",
         "problem/user_gamma=2.0e-3"],
    )
    row = _first_evolved_row(data)
    assert np.isclose(data["cool_gross"][row], 1.0e-2, rtol=1.0e-10,
                      atol=1.0e-12)
    assert np.isclose(data["cool_net"][row], 8.0e-3, rtol=1.0e-10,
                      atol=1.0e-12)
