"""MPI coverage for domain-integrated cooling history reductions."""

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


def _first_evolved_row(data):
    for row in range(1, len(data["time"])):
        if data["time"][row] > data["time"][row - 1]:
            return row
    pytest.fail("history file did not contain an evolved row")


def test_cooling_history_is_reduced_across_mpi_ranks():
    basename = "cooling_mpi_history"
    _remove_outputs(basename)
    testutils.mpi_run(
        HYDRO_INPUT,
        [f"job/basename={basename}",
         "mesh/nx1=16",
         "meshblock/nx1=8"],
        threads=2,
    )
    data = athena_read.hst(f"{basename}.hydro.hst", raw=True)
    row = _first_evolved_row(data)
    assert np.isclose(data["cool_gross"][row], 1.0e-2, rtol=2.0e-12,
                      atol=1.0e-13)
    assert np.isclose(data["cool_net"][row], 1.0e-2, rtol=2.0e-12,
                      atol=1.0e-13)
