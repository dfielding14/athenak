"""Uniform-grid streaming-instability linear-growth test."""

from subprocess import Popen, PIPE

import numpy as np
import pytest

import athena_read
import test_suite.testutils as testutils


EXPECTED_GROWTH = 0.1424527496587504


def test_streaming_instability_uniform_growth():
    """Run the streaming eigenmode and compare against the linear growth rate."""
    try:
        results = testutils.run("inputs/particles/streaming_instability.athinput")
        assert results, "streaming-instability run failed"

        data = athena_read.hst("streaming_instability.user.hst")
        gas_mode = np.hypot(data["gmodec"], data["gmodes"])
        particle_mode = np.hypot(data["pmodec"], data["pmodes"])
        mode = gas_mode + particle_mode

        if not np.all(np.isfinite(mode)):
            pytest.fail("streaming mode history contains non-finite values")

        mass_drift = np.max(np.abs(data["pmass"] - data["pmass"][0])) / data["pmass"][0]
        if mass_drift > 1.0e-12:
            pytest.fail(f"particle mass changed during streaming run: {mass_drift:g}")

        mask = (data["time"] >= 0.0) & (data["time"] <= 8.0) & (mode > 0.0)
        gamma = np.polyfit(data["time"][mask], np.log(mode[mask]), 1)[0]
        rel_err = abs(gamma - EXPECTED_GROWTH) / EXPECTED_GROWTH
        if rel_err > 0.10:
            pytest.fail(
                f"streaming growth {gamma:g} differs from linear value "
                f"{EXPECTED_GROWTH:g} by {rel_err:g}"
            )
    finally:
        Popen(["rm -f streaming_instability.user.hst"], shell=True,
              stdout=PIPE).communicate()
