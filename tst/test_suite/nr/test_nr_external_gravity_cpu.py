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
hydrostatic_input_file = "inputs/external_gravity_hydrostatic.athinput"
mhd_input_file = "inputs/external_gravity_mhd.athinput"
history_file = "external_gravity.hydro.hst"
force_history_file = "external_gravity_force.hydro.hst"
hydrostatic_tab_file = "tab/external_gravity_hydrostatic.hydro_w.00001.tab"
mhd_history_file = "external_gravity_mhd.mhd.hst"


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


def test_external_gravity_hydrostatic_boundary():
    try:
        results = testutils.run(hydrostatic_input_file)
        assert results, "External-gravity hydrostatic boundary run failed."

        data = athena_read.tab(hydrostatic_tab_file)
        expected_density = np.exp(-0.1*data["x1v"])
        assert np.max(np.abs(data["dens"] - expected_density)) < 2.0e-4
        assert np.max(np.abs(data["velx"])) < 2.0e-3
    finally:
        testutils.cleanup()


def test_external_gravity_mhd_history():
    try:
        results = testutils.run(mhd_input_file)
        assert results, "External-gravity MHD history run failed."

        data = athena_read.hst(mhd_history_file)
        assert "grav-PE" in data, "grav-PE was not added to MHD history output."
    finally:
        if os.path.exists(mhd_history_file):
            os.remove(mhd_history_file)
        testutils.cleanup()
