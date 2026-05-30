"""Streaming-instability linear-growth test with AMR and MPI."""

from subprocess import Popen, PIPE
from pathlib import Path

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


def test_drag_fast_particles_migrate_across_ranks():
    """Limit fast particles to neighboring MeshBlock migrations under MPI."""
    history = "particle_drag_fast_mpi.user.hst"
    try:
        args = [
            "job/basename=particle_drag_fast_mpi",
            "problem/particle_vx=100.0",
            "drag_particles/stopping_time=100.0",
            "drag_particles/cfl_drag=1.0",
            "time/tlim=0.03",
            "output1/dt=0.01",
        ]
        results = testutils.mpi_run(
            "inputs/tests/particle_drag_relaxation.athinput", args, threads=2
        )
        assert results, "fast MPI particle migration run failed"

        data = athena_read.hst(history)
        pmass = data["pmass"]
        mass_drift = np.max(np.abs(pmass - pmass[0])) / pmass[0]
        if mass_drift > 1.0e-12:
            pytest.fail(
                f"particle mass changed during fast MPI migration: {mass_drift:g}"
            )
        if not np.all(np.isfinite(data["gmom1"] + data["pmom1"])):
            pytest.fail("fast MPI particle migration history contains non-finite values")
    finally:
        Popen(["rm", "-f", history], stdout=PIPE).communicate()


def test_drag_particle_outputs_with_sparse_mpi_rank():
    """Write particle outputs when one MPI rank owns no particles."""
    try:
        Popen(["rm", "-rf", "pvtk", "trk"], stdout=PIPE).communicate()
        args = [
            "mesh/nx1=24",
            "meshblock/nx1=8",
            "particles/ppc=0.005",
            "time/tlim=0.001",
            "output1/dt=0.001",
            "output2/dt=0.001",
        ]
        results = testutils.mpi_run(
            "inputs/tests/particle_drag_outputs.athinput", args, threads=2
        )
        assert results, "sparse-rank MPI particle output run failed"

        pvtk = list(Path("pvtk").glob("particle_drag_outputs.particles.*.part.vtk"))
        if not pvtk or any(path.stat().st_size == 0 for path in pvtk):
            pytest.fail("sparse-rank MPI particle VTK output was not written")
        tracked = Path("trk/particle_drag_outputs.trk")
        if not tracked.exists() or tracked.stat().st_size == 0:
            pytest.fail("sparse-rank MPI tracked-particle output was not written")
    finally:
        Popen(["rm", "-rf", "pvtk", "trk"], stdout=PIPE).communicate()


def test_drag_restart_user_hook_reenrollment_mpi():
    """Restart user-hook drag particles across two MPI ranks."""
    history = "particle_drag_restart.user.hst"
    try:
        Popen(["rm", "-rf", "rst", history], stdout=PIPE).communicate()
        results = testutils.mpi_run(
            "inputs/tests/particle_drag_restart.athinput", threads=2
        )
        assert results, "pre-restart MPI drag run failed"
        rst = Path("rst/particle_drag_restart.00000.rst")
        if not rst.exists():
            pytest.fail(f"MPI restart file was not written: {rst}")

        Popen(["rm", "-f", history], stdout=PIPE).communicate()
        restarted = testutils.run_command(
            ["mpirun", "-np", "2", "./athena", "-r", str(rst), "time/tlim=0.1"]
        )
        assert restarted, "MPI drag restart run failed"

        data = athena_read.hst(history)
        pmass = data["pmass"]
        mass_drift = np.max(np.abs(pmass - pmass[0])) / pmass[0]
        if mass_drift > 1.0e-12:
            pytest.fail(f"particle mass changed across MPI restart: {mass_drift:g}")
    finally:
        Popen(["rm", "-rf", "rst", history], stdout=PIPE).communicate()
