"""MPI regression coverage for uniform-grid fourth-derivative hyperviscosity."""

import shutil

import numpy as np

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"


def test_uniform_mpi_hyperviscosity_matches_single_rank():
    """Uniform multi-rank STS exchange preserves the analytic shear damping result."""
    try:
        flags = [
            "mesh/nx1=128",
            "meshblock/nx1=32",
            "time/tlim=0.02",
            "output1/dt=0.02",
            "output2/dt=0.02",
        ]
        testutils.run(
            f"{INPUT_ROOT}/sts_hyperviscous_shear.athinput",
            ["job/basename=hypervisc_single", *flags],
        )
        testutils.mpi_run(
            f"{INPUT_ROOT}/sts_hyperviscous_shear.athinput",
            ["job/basename=hypervisc_mpi", *flags],
            threads=2,
        )
        single = np.loadtxt("hypervisc_single-hypervisc-errors.dat")
        parallel = np.loadtxt("hypervisc_mpi-hypervisc-errors.dat")
        assert np.max(np.abs(single[2:] - parallel[2:])) < 1.0e-12
    finally:
        testutils.cleanup()
        shutil.rmtree("bin", ignore_errors=True)
