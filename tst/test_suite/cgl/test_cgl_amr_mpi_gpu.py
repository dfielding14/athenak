"""Opt-in MPI+GPU regressions for CGL-aware AMR transfers."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import numpy as np
import pytest

VIS_PYTHON = Path(__file__).resolve().parents[3] / "vis" / "python"
sys.path.insert(0, str(VIS_PYTHON))

import bin_convert  # noqa: E402
import athena_read  # noqa: E402

athena_read.check_nan_flag = True


INPUT_ROOT = "../../../inputs/tests"
RUN_MPI_GPU = os.environ.get("ATHENAK_RUN_MPI_GPU") == "1"
MPI_GPU_LAUNCHER = os.environ.get("ATHENAK_MPI_GPU_LAUNCHER", "")
AMR_REPAIR_COLUMNS = (
    "amr_cell",
    "amr_nan",
    "amr_rho",
    "amr_U",
    "amr_ppar",
    "amr_pperp",
    "amr_lowB",
    "amr_fh",
    "amr_mirr",
    "amr_dlt",
    "amr_int",
    "amr_slp",
)

pytestmark = pytest.mark.skipif(
    not RUN_MPI_GPU,
    reason="set ATHENAK_RUN_MPI_GPU=1 inside an MPI+GPU allocation",
)


def _launch(nodes, ranks, input_name, basename, *flags, restart_file=None):
    assert MPI_GPU_LAUNCHER, "ATHENAK_MPI_GPU_LAUNCHER is required"
    tasks_per_node = (ranks + nodes - 1) // nodes
    command = [
        *shlex.split(MPI_GPU_LAUNCHER),
        "-N",
        str(nodes),
        "-n",
        str(ranks),
        "--ntasks-per-node",
        str(tasks_per_node),
        "./athena",
    ]
    if restart_file is None:
        command.extend(["-i", f"{INPUT_ROOT}/{input_name}"])
    else:
        command.extend(["-r", str(restart_file)])
    command.extend([f"job/basename={basename}", *flags])
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


def _cleanup():
    for path in Path(".").glob("cgl_mpigpu*.hst"):
        path.unlink()
    for path in Path(".").glob("*.dat"):
        path.unlink()
    shutil.rmtree("bin", ignore_errors=True)
    shutil.rmtree("rst", ignore_errors=True)
    shutil.rmtree("tab", ignore_errors=True)


def _history(basename, kind):
    return athena_read.hst(f"{basename}.{kind}.hst")


def _latest_state(basename):
    paths = sorted(Path("bin").glob(f"{basename}.state.*.bin"))
    assert paths, f"no state output found for {basename}"
    return bin_convert.read_binary(str(paths[-1]))


def _assert_clean(user, mhd):
    assert np.max(user["bad_state"]) == 0.0
    assert np.max(user["max_ndiv"]) < 1.0e-12
    for column in AMR_REPAIR_COLUMNS:
        assert user[column][-1] == 0.0
    for column in (
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_hardbd",
        "lf_mirror",
        "lf_firehs",
        "lf_hwproj",
    ):
        if column in mhd:
            assert mhd[column][-1] == 0.0


def _assert_conserved(mhd):
    for column in ("mass", "1-mom", "2-mom", "3-mom", "tot-E"):
        np.testing.assert_allclose(
            mhd[column],
            mhd[column][0],
            rtol=0.0,
            atol=5.0e-12,
            err_msg=f"MPI+GPU AMR changed {column}",
        )


def _assert_same_state(reference, candidate):
    assert candidate["time"] == reference["time"]
    assert candidate["cycle"] == reference["cycle"]
    assert candidate["var_names"] == reference["var_names"]
    assert candidate["n_mbs"] == reference["n_mbs"]
    reference_order = sorted(
        range(reference["n_mbs"]),
        key=lambda block: tuple(reference["mb_logical"][block]),
    )
    candidate_order = sorted(
        range(candidate["n_mbs"]),
        key=lambda block: tuple(candidate["mb_logical"][block]),
    )
    for reference_block, candidate_block in zip(reference_order, candidate_order):
        np.testing.assert_array_equal(
            reference["mb_logical"][reference_block],
            candidate["mb_logical"][candidate_block],
        )
        for variable in reference["var_names"]:
            np.testing.assert_array_equal(
                reference["mb_data"][variable][reference_block],
                candidate["mb_data"][variable][candidate_block],
            )


def test_cgl_amr_restart_and_conservation_mpi_gpu():
    try:
        reference = "cgl_mpigpu_restart_reference"
        partial = "cgl_mpigpu_restart_partial"
        resumed = "cgl_mpigpu_restart_resumed"
        _launch(
            1,
            4,
            "cgl_amr_primitive_current_churn_restart.athinput",
            reference,
            "time/tlim=0.002",
            "time/nlim=-1",
        )
        _launch(
            1,
            4,
            "cgl_amr_primitive_current_churn_restart.athinput",
            partial,
            "time/tlim=0.002",
            "time/nlim=3",
        )
        restarts = sorted(Path("rst").rglob(f"{partial}*.rst"))
        assert restarts
        _launch(
            1,
            4,
            "",
            resumed,
            "time/nlim=-1",
            restart_file=restarts[-1],
        )

        reference_user = _history(reference, "user")
        reference_mhd = _history(reference, "mhd")
        resumed_user = _history(resumed, "user")
        resumed_mhd = _history(resumed, "mhd")
        _assert_clean(reference_user, reference_mhd)
        _assert_clean(resumed_user, resumed_mhd)
        _assert_conserved(reference_mhd)
        transitions = np.diff(reference_user["ncell"])
        assert np.any(transitions > 0.0)
        assert np.any(transitions < 0.0)
        _assert_same_state(_latest_state(reference), _latest_state(resumed))
    finally:
        _cleanup()


def test_cgl_lf_amr_3d_churn_mpi_gpu():
    try:
        basename = "cgl_mpigpu_lf_3d_churn"
        _launch(2, 16, "cgl_lf_amr_3d_current.athinput", basename)
        user = _history(basename, "user")
        mhd = _history(basename, "mhd")
        _assert_clean(user, mhd)
        _assert_conserved(mhd)
        assert user["ncell"][0] == 13824.0
        assert np.max(user["ncell"]) == 110592.0
        assert user["ncell"][-1] == user["ncell"][0]
        assert mhd["lf_nstage"][-1] > 0.0
    finally:
        _cleanup()
