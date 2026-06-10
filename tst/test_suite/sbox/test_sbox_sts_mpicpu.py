"""STS regression coverage with shearing-periodic boundaries."""

from pathlib import Path
import shutil
import sys

import numpy as np

VIS_PYTHON = Path(__file__).resolve().parents[3] / "vis" / "python"
sys.path.insert(0, str(VIS_PYTHON))

import test_suite.testutils as testutils  # noqa: E402
import bin_convert  # noqa: E402


HYDRO_INPUT = "inputs/hydro_sts_sbox.athinput"
SCALAR_INPUT = "inputs/hydro_scalar_sts_sbox.athinput"
MHD_INPUT = "inputs/mhd_sts_sbox.athinput"


def _clean_outputs():
    testutils.cleanup()
    shutil.rmtree("bin", ignore_errors=True)
    shutil.rmtree("rst", ignore_errors=True)


def _run(input_file, basename, *flags, mpi=False):
    arguments = [f"job/basename={basename}", *flags]
    if mpi:
        return testutils.mpi_run(input_file, arguments, threads=2)
    return testutils.run(input_file, arguments)


def _read(basename, variable, index=1):
    return bin_convert.read_binary(
        f"bin/{basename}.{variable}.{index:05d}.bin"
    )


def _sorted_field(data, variable):
    order = sorted(
        range(data["n_mbs"]),
        key=lambda block: tuple(data["mb_logical"][block]),
    )
    return np.concatenate(
        [np.asarray(data["mb_data"][variable][block]).ravel() for block in order]
    )


def _compare_outputs(left, right, *, atol):
    assert left["var_names"] == right["var_names"]
    for variable in left["var_names"]:
        np.testing.assert_allclose(
            _sorted_field(left, variable),
            _sorted_field(right, variable),
            rtol=2.0e-6,
            atol=atol,
            err_msg=f"mismatch in {variable}",
        )


def _assert_magnetic_fluxes(output):
    time = output["time"]
    np.testing.assert_allclose(
        np.mean(_sorted_field(output, "bcc1")), 0.1, rtol=0.0, atol=2.0e-6
    )
    np.testing.assert_allclose(
        np.mean(_sorted_field(output, "bcc2")),
        0.2 - 1.5*0.1*time,
        rtol=0.0,
        atol=2.0e-6,
    )
    np.testing.assert_allclose(
        np.mean(_sorted_field(output, "bcc3")), 0.0, rtol=0.0, atol=2.0e-6
    )


def _restart(restart_file, basename, *, mpi=False):
    command = []
    if mpi:
        command.extend(["mpirun", "-np", "2"])
    command.extend([
        "./athena",
        "-r",
        restart_file,
        f"job/basename={basename}",
    ])
    if not testutils.run_command(command):
        raise RuntimeError(f"Failed to restart from {restart_file}")


def test_hydro_sts_shearing_box_serial_mpi_and_explicit():
    """Hydro cell-centered STS is decomposition-independent and matches explicit."""
    try:
        _clean_outputs()
        _run(HYDRO_INPUT, "hydro_sts")
        _run(HYDRO_INPUT, "hydro_mpi", mpi=True)
        oracle_flags = ("time/tlim=0.08", "output1/dt=0.08")
        _run(
            HYDRO_INPUT,
            "hydro_capped",
            *oracle_flags,
            "time/sts_max_dt_ratio=1.0",
        )
        _run(
            HYDRO_INPUT,
            "hydro_explicit",
            *oracle_flags,
            "hydro/viscosity_integrator=explicit",
            "time/sts_integrator=none",
        )

        sts = _read("hydro_sts", "hydro_w")
        mpi = _read("hydro_mpi", "hydro_w")
        capped = _read("hydro_capped", "hydro_w")
        explicit = _read("hydro_explicit", "hydro_w")

        _compare_outputs(sts, mpi, atol=2.0e-5)
        _compare_outputs(capped, explicit, atol=2.0e-5)
        assert sts["cycle"]/sts["time"] < explicit["cycle"]/explicit["time"]
    finally:
        _clean_outputs()


def test_scalar_sts_shearing_box_conserves_mass_across_mpi_seam():
    """Scalar STS conserves its integral with remote shearing partners."""
    try:
        _clean_outputs()
        _run(SCALAR_INPUT, "scalar_sts")
        _run(SCALAR_INPUT, "scalar_mpi", mpi=True)

        sts = _read("scalar_sts", "hydro_u_s")
        mpi = _read("scalar_mpi", "hydro_u_s")

        _compare_outputs(sts, mpi, atol=3.0e-6)
        assert sts["var_names"] == ["r_00"]
        np.testing.assert_allclose(
            np.mean(_sorted_field(sts, "r_00")), 1.0, rtol=0.0, atol=2.0e-12
        )
        np.testing.assert_allclose(
            np.mean(_sorted_field(mpi, "r_00")), 1.0, rtol=0.0, atol=2.0e-12
        )
    finally:
        _clean_outputs()


def test_mhd_resistive_sts_shearing_box_serial_mpi_divb_and_restart():
    """Resistive MHD STS preserves CT, decomposition, explicit, and restart results."""
    try:
        _clean_outputs()
        _run(MHD_INPUT, "mhd_sts")
        _run(MHD_INPUT, "mhd_mpi", mpi=True)
        midpoint = Path("rst/mhd_mpi.00001.rst")
        assert midpoint.exists()
        oracle_flags = (
            "time/tlim=0.08",
            "output1/dt=0.08",
            "output2/dt=0.08",
            "output3/dt=-1",
        )
        _run(
            MHD_INPUT,
            "mhd_capped",
            *oracle_flags,
            "time/sts_max_dt_ratio=1.0",
        )
        _run(
            MHD_INPUT,
            "mhd_explicit",
            *oracle_flags,
            "mhd/viscosity_integrator=explicit",
            "mhd/ohmic_resistivity_integrator=explicit",
            "time/sts_integrator=none",
        )
        _restart(str(midpoint), "mhd_mpi_restart", mpi=True)

        sts = _read("mhd_sts", "mhd_w_bcc")
        mpi = _read("mhd_mpi", "mhd_w_bcc")
        capped = _read("mhd_capped", "mhd_w_bcc")
        explicit = _read("mhd_explicit", "mhd_w_bcc")
        restarted = _read("mhd_mpi_restart", "mhd_w_bcc")
        divb = _read("mhd_sts", "mhd_divb")
        mpi_divb = _read("mhd_mpi", "mhd_divb")
        restarted_divb = _read("mhd_mpi_restart", "mhd_divb")

        _compare_outputs(sts, mpi, atol=3.0e-6)
        _compare_outputs(capped, explicit, atol=3.0e-5)
        _compare_outputs(mpi, restarted, atol=3.0e-6)
        assert restarted["cycle"] == mpi["cycle"]
        assert restarted["time"] == mpi["time"]
        _assert_magnetic_fluxes(sts)
        _assert_magnetic_fluxes(mpi)
        _assert_magnetic_fluxes(restarted)
        assert sts["cycle"]/sts["time"] < explicit["cycle"]/explicit["time"]
        assert np.max(np.abs(_sorted_field(divb, "divb"))) < 5.0e-6
        assert np.max(np.abs(_sorted_field(mpi_divb, "divb"))) < 5.0e-6
        assert np.max(np.abs(_sorted_field(restarted_divb, "divb"))) < 5.0e-6
    finally:
        _clean_outputs()
