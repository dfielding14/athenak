"""
External-gravity source-term smoke test.

This checks that a standalone <external_gravity> block enables source-term parsing and
adds gravitational potential energy to the default hydro history output.
"""

import os

import athena_read
import numpy as np
import test_suite.testutils as testutils


input_file = "inputs/external_gravity.athinput"
history_file = "external_gravity.hydro.hst"
force_history_file = "external_gravity_force.hydro.hst"


def test_external_gravity_history():
    try:
        results = testutils.run(input_file)
        assert results, "External-gravity test run failed."

        data = athena_read.hst(history_file)
        assert "grav-PE" in data, "grav-PE was not added to hydro history output."

        nx1 = 16
        dx1 = 1.0/nx1
        x1 = (np.arange(nx1) + 0.5)*dx1
        expected = np.sum(0.5*x1*x1*dx1)
        assert np.isclose(data["grav-PE"][0], expected, rtol=1.0e-12, atol=1.0e-12)

        os.remove(history_file)

        results = testutils.run(
            input_file,
            [
                "job/basename=external_gravity_force",
                "time/nlim=1",
                "time/tlim=1.0",
                "external_gravity/model=uniform",
                "external_gravity/g1=0.25",
                "external_gravity/g2=0.0",
                "external_gravity/g3=0.0",
            ],
        )
        assert results, "External-gravity force test run failed."

        data = athena_read.hst(force_history_file)
        expected_momentum = 0.25*data["time"][-1]
        assert np.isclose(data["1-mom"][-1], expected_momentum, rtol=1.0e-12,
                          atol=1.0e-12)
    finally:
        if os.path.exists(history_file):
            os.remove(history_file)
        if os.path.exists(force_history_file):
            os.remove(force_history_file)
        testutils.cleanup()
