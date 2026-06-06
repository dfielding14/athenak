"""
MPI self-gravity regression tests.

These tests are intentionally small but exercise distributed source loading,
root transfers, and multipole coefficient reductions.
"""

import math

import numpy as np
import pytest

import test_suite.testutils as testutils
from test_suite.selfgravity import selfgravity_utils as sg


@pytest.mark.parametrize("ranks", [2, 4])
def test_mpi_periodic_jeans_converges_mpicpu(ranks):
    basename = f"selfgravity_jeans_{ranks}rank_mpicpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            basename,
            mpi_ranks=ranks,
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        output = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
        k_wave, four_pi_g, _ = sg.jeans_geometry(n_jeans=0.5)
        analytic = -four_pi_g * 1.0e-3 / (k_wave * k_wave)
        analytic *= np.sin(k_wave * sg.jeans_phase(x1, x2, x3))
        assert sg.relative_l2((phi - np.mean(phi)) - analytic, analytic) < 4.0e-2
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_mpi_root_host_and_device_paths_match_mpicpu():
    device_basename = "selfgravity_root_device_mpicpu"
    host_basename = "selfgravity_root_host_mpicpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            device_basename,
            mpi_ranks=2,
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            host_basename,
            extra_args=["gravity/root_on_host=true"],
            mpi_ranks=2,
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        device = sg.parse_binary_output(
            sg.latest_binary_output(device_basename, "grav_phi")
        )
        host = sg.parse_binary_output(sg.latest_binary_output(host_basename, "grav_phi"))
        device_phi, _, _, _ = sg.concatenate_field(device, "grav_phi")
        host_phi, _, _, _ = sg.concatenate_field(host, "grav_phi")
        assert np.max(np.abs(device_phi - host_phi)) < 1.0e-10
    finally:
        sg.cleanup_outputs(device_basename, host_basename)
        testutils.cleanup()


def test_mpi_multipole_binary_reduction_mpicpu():
    basename = "selfgravity_binary_multipole_mpicpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity_binary.athinput",
            basename,
            mpi_ranks=2,
            max_final_defect=2.0e-6,
        )
        output = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
        gconst = 1.0 / (4.0 * math.pi)
        radius = 0.10
        r1 = np.sqrt((x1 - 0.15) ** 2 + x2 * x2 + x3 * x3)
        r2 = np.sqrt((x1 + 0.15) ** 2 + x2 * x2 + x3 * x3)
        analytic = -gconst * 1.0 / np.maximum(r1, 1.0e-300)
        analytic += -gconst * 0.5 / np.maximum(r2, 1.0e-300)
        mask = (r1 > 2.0 * radius) & (r2 > 2.0 * radius)
        assert np.all(phi[mask] < 0.0)
        assert sg.relative_l2(phi[mask] - analytic[mask], analytic[mask]) < 4.0e-1
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_mpi_partial_smr_matches_serial_mpicpu():
    # Fixed work isolates rank-dependent fine/coarse communication from convergence.
    serial_basename = "selfgravity_binary_partial_smr_serial"
    mpi_basename = "selfgravity_binary_partial_smr_mpicpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity_binary_smr.athinput",
            serial_basename,
            max_final_defect=None,
        )
        sg.run_athena(
            "inputs/tests/selfgravity_binary_smr.athinput",
            mpi_basename,
            mpi_ranks=4,
            max_final_defect=None,
        )

        serial = sg.parse_binary_output(
            sg.latest_binary_output(serial_basename, "grav_phi")
        )
        mpi = sg.parse_binary_output(sg.latest_binary_output(mpi_basename, "grav_phi"))
        assert {block["level"] for block in serial["blocks"]} == {0, 1}
        assert {block["level"] for block in mpi["blocks"]} == {0, 1}

        serial_phi, serial_x1, serial_x2, serial_x3 = sg.concatenate_field(
            serial, "grav_phi"
        )
        mpi_phi, mpi_x1, mpi_x2, mpi_x3 = sg.concatenate_field(mpi, "grav_phi")
        serial_order = np.lexsort((serial_x1, serial_x2, serial_x3))
        mpi_order = np.lexsort((mpi_x1, mpi_x2, mpi_x3))
        serial_coords = np.column_stack(
            (serial_x1[serial_order], serial_x2[serial_order], serial_x3[serial_order])
        )
        mpi_coords = np.column_stack(
            (mpi_x1[mpi_order], mpi_x2[mpi_order], mpi_x3[mpi_order])
        )
        assert np.max(np.abs(serial_coords - mpi_coords)) < 1.0e-14
        assert np.max(np.abs(serial_phi[serial_order] - mpi_phi[mpi_order])) < 1.0e-10
    finally:
        sg.cleanup_outputs(serial_basename, mpi_basename)
        testutils.cleanup()
