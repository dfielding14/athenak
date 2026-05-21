"""Streaming-instability linear-growth test with AMR and MPI."""

from subprocess import Popen, PIPE

import numpy as np
import pytest

import athena_read
import test_suite.testutils as testutils


EXPECTED_GROWTH = 0.1424527496587504


def test_streaming_instability_amr_mpi_growth():
    """Run the AMR/MPI streaming eigenmode and compare the linear growth rate."""
    try:
        args = [
            "time/tlim=4.0",
            "output1/dt=0.05",
        ]
        results = testutils.mpi_run(
            "inputs/particles/streaming_instability_amr.athinput", args, threads=2
        )
        assert results, "AMR/MPI streaming-instability run failed"

        data = athena_read.hst("streaming_instability_amr.user.hst")
        pmass = data["pmass"]
        mass_drift = np.max(np.abs(pmass - pmass[0])) / pmass[0]
        if mass_drift > 1.0e-12:
            pytest.fail(f"particle mass changed during AMR/MPI run: {mass_drift:g}")

        gas_mode = np.hypot(data["gmodec"], data["gmodes"])
        particle_mode = np.hypot(data["pmodec"], data["pmodes"])
        if not np.all(np.isfinite(gas_mode + particle_mode)):
            pytest.fail("streaming mode history contains non-finite values")

        mode = gas_mode + particle_mode
        mask = (data["time"] >= 0.0) & (data["time"] <= 4.0) & (mode > 0.0)
        gamma = np.polyfit(data["time"][mask], np.log(mode[mask]), 1)[0]
        rel_err = abs(gamma - EXPECTED_GROWTH) / EXPECTED_GROWTH
        if rel_err > 0.15:
            pytest.fail(
                f"AMR/MPI streaming growth {gamma:g} differs from "
                f"linear value {EXPECTED_GROWTH:g} by {rel_err:g}"
            )
    finally:
        Popen(["rm -f streaming_instability_amr.user.hst"], shell=True,
              stdout=PIPE).communicate()


def test_streaming_instability_dynamic_amr_mpi_remap():
    """Trigger AMR creation under MPI and verify particle remapping conserves mass."""
    try:
        results = testutils.mpi_run(
            "inputs/particles/streaming_instability_amr_dynamic.athinput", threads=2
        )
        assert results, "dynamic AMR/MPI streaming-instability run failed"

        data = athena_read.hst("streaming_instability_amr_dynamic.user.hst")
        pmass = data["pmass"]
        mass_drift = np.max(np.abs(pmass - pmass[0])) / pmass[0]
        if mass_drift > 1.0e-12:
            pytest.fail(f"particle mass changed during AMR remap: {mass_drift:g}")

        gas_mode = np.hypot(data["gmodec"], data["gmodes"])
        particle_mode = np.hypot(data["pmodec"], data["pmodes"])
        if not np.all(np.isfinite(gas_mode + particle_mode)):
            pytest.fail("dynamic AMR mode history contains non-finite values")
    finally:
        Popen(["rm -f streaming_instability_amr_dynamic.user.hst"], shell=True,
              stdout=PIPE).communicate()


def test_streaming_instability_dynamic_amr_mpi_many_particles():
    """Stress AMR/MPI remapping with more particles and smaller MeshBlocks."""
    history = "streaming_instability_amr_dynamic_many_particles.user.hst"
    try:
        args = [
            "job/basename=streaming_instability_amr_dynamic_many_particles",
            "meshblock/nx1=8",
            "meshblock/nx3=8",
            "particles/ppc=2.0",
            "time/tlim=0.2",
            "output1/dt=0.1",
        ]
        results = testutils.mpi_run(
            "inputs/particles/streaming_instability_amr_dynamic.athinput",
            args,
            threads=4,
        )
        assert results, "higher-count dynamic AMR/MPI particle run failed"

        data = athena_read.hst(history)
        pmass = data["pmass"]
        mass_drift = np.max(np.abs(pmass - pmass[0])) / pmass[0]
        if mass_drift > 1.0e-12:
            pytest.fail(
                f"particle mass changed during higher-count AMR remap: {mass_drift:g}"
            )

        mode = np.hypot(data["gmodec"], data["gmodes"]) + np.hypot(
            data["pmodec"], data["pmodes"]
        )
        if not np.all(np.isfinite(mode)):
            pytest.fail("higher-count dynamic AMR history contains non-finite values")
    finally:
        Popen(["rm", "-f", history], stdout=PIPE).communicate()
