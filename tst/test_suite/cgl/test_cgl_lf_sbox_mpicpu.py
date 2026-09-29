"""CGL Landau-fluid STS regression with shearing-periodic boundaries."""

from pathlib import Path
import shutil
import sys

import numpy as np

VIS_PYTHON = Path(__file__).resolve().parents[3] / "vis" / "python"
sys.path.insert(0, str(VIS_PYTHON))

import bin_convert  # noqa: E402
import test_suite.testutils as testutils  # noqa: E402


INPUT_FILE = "inputs/cgl_lf_sts_sbox.athinput"


def _clean_outputs():
    testutils.cleanup()
    shutil.rmtree("bin", ignore_errors=True)
    shutil.rmtree("rst", ignore_errors=True)
    for path in Path(".").glob("cgl_sbox_*.hst"):
        path.unlink()


def _run(basename, *flags, mpi=False):
    arguments = [f"job/basename={basename}", *flags]
    if mpi:
        return testutils.mpi_run(INPUT_FILE, arguments, threads=2)
    return testutils.run(INPUT_FILE, arguments)


def _read(basename, variable):
    return bin_convert.read_binary(f"bin/{basename}.{variable}.00001.bin")


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
            rtol=3.0e-6,
            atol=atol,
            err_msg=f"mismatch in {variable}",
        )


def _assert_admissible(output, history):
    for variable in output["var_names"]:
        assert np.all(np.isfinite(_sorted_field(output, variable)))
    assert np.min(_sorted_field(output, "dens")) > 0.0
    assert np.min(_sorted_field(output, "eint")) > 0.0
    assert np.min(_sorted_field(output, "p_perp")) > 0.0
    assert history["lf_nstage"][-1] > 0.0
    assert history["lf_qface"][-1] > 0.0
    assert abs(history["lf_qprwrk"][-1]) > 0.0
    assert abs(history["lf_qpewrk"][-1]) > 0.0
    for column in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd"):
        assert history[column][-1] == 0.0


def _assert_magnetic_state(output, divb):
    time = output["time"]
    np.testing.assert_allclose(
        np.mean(_sorted_field(output, "bcc1")), 0.3, rtol=0.0, atol=3.0e-6
    )
    np.testing.assert_allclose(
        np.mean(_sorted_field(output, "bcc2")),
        0.2 - 1.5*0.3*time,
        rtol=0.0,
        atol=3.0e-6,
    )
    np.testing.assert_allclose(
        np.mean(_sorted_field(output, "bcc3")), 0.1, rtol=0.0, atol=3.0e-6
    )
    assert np.max(np.abs(_sorted_field(divb, "divb"))) < 5.0e-6


def _assert_restarted_diagnostics(reference, resumed):
    for column in (
        "lf_nstage",
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_hardbd",
        "lf_qface",
        "lf_qprcap",
        "lf_qpr10",
        "lf_qpecap",
        "lf_qpe10",
        "lf_mirror",
        "lf_firehs",
        "lf_hwproj",
    ):
        assert resumed[column][-1] == reference[column][-1]
    for column in ("lf_qprwrk", "lf_qpewrk"):
        np.testing.assert_allclose(
            resumed[column][-1],
            reference[column][-1],
            rtol=0.0,
            atol=2.0e-11,
            err_msg=f"restart mismatch in {column}",
        )


def _assert_decomposition_independent_diagnostics(serial, mpi):
    for column in (
        "lf_nstage",
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_hardbd",
        "lf_qface",
        "lf_qprcap",
        "lf_qpr10",
        "lf_qpecap",
        "lf_qpe10",
        "lf_mirror",
        "lf_firehs",
        "lf_hwproj",
    ):
        assert mpi[column][-1] == serial[column][-1]
    for column in ("lf_qprwrk", "lf_qpewrk"):
        np.testing.assert_allclose(
            mpi[column][-1],
            serial[column][-1],
            rtol=0.0,
            atol=2.0e-11,
            err_msg=f"decomposition mismatch in {column}",
        )


def _restart(restart_file, basename):
    command = [
        "mpirun",
        "-np",
        "2",
        "./athena",
        "-r",
        restart_file,
        f"job/basename={basename}",
    ]
    if not testutils.run_command(command):
        raise RuntimeError(f"Failed to restart from {restart_file}")


def test_cgl_lf_sts_shearing_box_serial_mpi_explicit_and_restart():
    try:
        _clean_outputs()
        _run("cgl_sbox_sts")
        _run("cgl_sbox_mpi", mpi=True)
        midpoint = Path("rst/cgl_sbox_mpi.00001.rst")
        assert midpoint.exists()

        oracle_flags = (
            "time/tlim=0.04",
            "output1/dt=0.04",
            "output2/dt=0.04",
            "output3/dt=-1",
        )
        _run(
            "cgl_sbox_capped",
            *oracle_flags,
            "time/sts_max_dt_ratio=1.0",
        )
        _run(
            "cgl_sbox_explicit",
            *oracle_flags,
            "mhd/cgl_heat_flux_integrator=explicit",
            "time/sts_integrator=none",
        )
        _restart(str(midpoint), "cgl_sbox_restart")

        sts = _read("cgl_sbox_sts", "mhd_w_bcc")
        mpi = _read("cgl_sbox_mpi", "mhd_w_bcc")
        capped = _read("cgl_sbox_capped", "mhd_w_bcc")
        explicit = _read("cgl_sbox_explicit", "mhd_w_bcc")
        restarted = _read("cgl_sbox_restart", "mhd_w_bcc")
        sts_divb = _read("cgl_sbox_sts", "mhd_divb")
        mpi_divb = _read("cgl_sbox_mpi", "mhd_divb")
        restarted_divb = _read("cgl_sbox_restart", "mhd_divb")
        sts_history = testutils.athena_read.hst("cgl_sbox_sts.mhd.hst")
        mpi_history = testutils.athena_read.hst("cgl_sbox_mpi.mhd.hst")
        restarted_history = testutils.athena_read.hst("cgl_sbox_restart.mhd.hst")

        _compare_outputs(sts, mpi, atol=5.0e-6)
        _compare_outputs(capped, explicit, atol=5.0e-5)
        _compare_outputs(mpi, restarted, atol=5.0e-6)
        _assert_admissible(sts, sts_history)
        _assert_admissible(mpi, mpi_history)
        _assert_decomposition_independent_diagnostics(sts_history, mpi_history)
        _assert_magnetic_state(sts, sts_divb)
        _assert_magnetic_state(mpi, mpi_divb)
        _assert_magnetic_state(restarted, restarted_divb)
        _assert_restarted_diagnostics(mpi_history, restarted_history)
        assert restarted["cycle"] == mpi["cycle"]
        assert restarted["time"] == mpi["time"]
        assert sts["cycle"]/sts["time"] < explicit["cycle"]/explicit["time"]
    finally:
        _clean_outputs()
