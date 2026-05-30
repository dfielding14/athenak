"""GPU smoke coverage for fourth-derivative hyperviscosity."""

import shutil

import numpy as np

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"


def test_gpu_hyperviscosity_shear_mode_is_finite_and_damped():
    """Execute the Kokkos device kernels and check the expected shear damping."""
    try:
        testutils.run(
            f"{INPUT_ROOT}/sts_hyperviscous_shear.athinput",
            [
                "job/basename=hypervisc_gpu",
                "mesh/nx1=64",
                "meshblock/nx1=64",
                "time/tlim=0.01",
                "output1/dt=0.01",
                "output2/dt=0.01",
            ],
        )
        errors = np.loadtxt("hypervisc_gpu-hypervisc-errors.dat")
        assert np.all(np.isfinite(errors))
        assert errors[4] < 0.1
    finally:
        testutils.cleanup()
        shutil.rmtree("bin", ignore_errors=True)
