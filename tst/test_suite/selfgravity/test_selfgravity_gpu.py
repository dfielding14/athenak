"""
GPU self-gravity smoke tests.

These are run only by `run_test_suite.py --gpu` on machines with a CUDA-capable
AthenaK build. They mirror the CPU correctness checks without assuming a GPU is
available on developer laptops.
"""

import numpy as np

import test_suite.testutils as testutils
from test_suite.selfgravity import selfgravity_utils as sg


def test_periodic_jeans_potential_gpu():
    basename = "selfgravity_jeans_gpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            basename,
            extra_args=["gravity/profile=true"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        output = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
        k_wave, four_pi_g, _ = sg.jeans_geometry(n_jeans=0.5)
        analytic = -four_pi_g * 1.0e-3 / (k_wave*k_wave)
        analytic *= np.sin(k_wave*sg.jeans_phase(x1, x2, x3))
        assert sg.relative_l2((phi - np.mean(phi)) - analytic, analytic) < 4.0e-2
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_root_host_and_device_paths_match_gpu():
    device_basename = "selfgravity_root_device_gpu"
    host_basename = "selfgravity_root_host_gpu"
    try:
        sg.run_athena("inputs/tests/selfgravity.athinput", device_basename,
                      expect_jeans_banner=True, max_final_defect=2.0e-8)
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            host_basename,
            extra_args=["gravity/root_on_host=true"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        device = sg.parse_binary_output(sg.latest_binary_output(device_basename, "grav_phi"))
        host = sg.parse_binary_output(sg.latest_binary_output(host_basename, "grav_phi"))
        device_phi, _, _, _ = sg.concatenate_field(device, "grav_phi")
        host_phi, _, _, _ = sg.concatenate_field(host, "grav_phi")
        assert np.max(np.abs(device_phi - host_phi)) < 1.0e-10
    finally:
        sg.cleanup_outputs(device_basename, host_basename)
        testutils.cleanup()
