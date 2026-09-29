"""GPU smoke coverage for CGL-aware AMR primitive transfers."""

from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

VIS_PYTHON = Path(__file__).resolve().parents[3] / "vis" / "python"
sys.path.insert(0, str(VIS_PYTHON))

import bin_convert  # noqa: E402
import test_suite.testutils as testutils  # noqa: E402


INPUT_ROOT = "../../../inputs/tests"
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


def _run(input_name, basename, *flags):
    testutils.run(
        f"{INPUT_ROOT}/{input_name}",
        [f"job/basename={basename}", *flags],
    )


def _cleanup(prefix="cgl_amr_gpu"):
    for path in Path(".").glob(f"{prefix}*.hst"):
        path.unlink()
    shutil.rmtree("bin", ignore_errors=True)
    shutil.rmtree("rst", ignore_errors=True)
    testutils.cleanup()


def _user_history(basename):
    return testutils.athena_read.hst(f"{basename}.user.hst")


def _mhd_history(basename):
    return testutils.athena_read.hst(f"{basename}.mhd.hst")


def _restart(restart_file, basename, *flags):
    command = ["./athena", "-r", str(restart_file), f"job/basename={basename}", *flags]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stdout + result.stderr)


def _latest_state(basename):
    paths = sorted(Path("bin").glob(f"{basename}.state.*.bin"))
    assert paths, f"no state output found for {basename}"
    return bin_convert.read_binary(str(paths[-1]))


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
                err_msg=f"restart changed {variable}",
            )


def _assert_clean_user(history, max_ndiv=1.0e-10):
    assert np.max(history["bad_state"]) == 0.0
    assert np.all(np.isfinite(history["abs_anis"]))
    assert np.max(history["max_ndiv"]) < max_ndiv


def _assert_no_amr_repairs(history):
    for column in AMR_REPAIR_COLUMNS:
        assert history[column][-1] == 0.0


def _assert_clean_lf(history):
    assert history["lf_nstage"][-1] > 0.0
    for column in (
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_hardbd",
    ):
        assert history[column][-1] == 0.0


def _assert_conserved(history, columns, atol=5.0e-12):
    for column in columns:
        np.testing.assert_allclose(
            history[column],
            history[column][0],
            rtol=0.0,
            atol=atol,
            err_msg=f"AMR changed conserved history column {column}",
        )


def test_cgl_amr_primitive_smooth_gpu():
    try:
        _run("cgl_amr_primitive_uniform.athinput", "cgl_amr_gpu_uniform")
        uniform = _user_history("cgl_amr_gpu_uniform")
        _assert_clean_user(uniform)
        _assert_no_amr_repairs(uniform)
        assert uniform["ncell"][-1] > uniform["ncell"][0]

        _run(
            "cgl_amr_primitive_strong_anisotropy.athinput",
            "cgl_amr_gpu_aniso",
        )
        aniso = _user_history("cgl_amr_gpu_aniso")
        _assert_clean_user(aniso)
        _assert_no_amr_repairs(aniso)
        assert np.max(aniso["ncell"]) > aniso["ncell"][0]
    finally:
        _cleanup()


def test_cgl_amr_primitive_stage_order_gpu():
    try:
        common_flags = (
            "time/nlim=1",
            "problem/current_refine_threshold=1.0e30",
        )
        _run(
            "cgl_amr_primitive_current.athinput",
            "cgl_amr_gpu_stage_primitive",
            *common_flags,
        )
        _run(
            "cgl_amr_primitive_current.athinput",
            "cgl_amr_gpu_stage_conserved",
            *common_flags,
            "mesh_refinement/prolong_primitives=false",
        )

        primitive = _mhd_history("cgl_amr_gpu_stage_primitive")
        conserved = _mhd_history("cgl_amr_gpu_stage_conserved")
        primitive_user = _user_history("cgl_amr_gpu_stage_primitive")
        conserved_user = _user_history("cgl_amr_gpu_stage_conserved")
        assert np.all(primitive_user["ncell"] == primitive_user["ncell"][0])
        np.testing.assert_array_equal(
            primitive_user["ncell"], conserved_user["ncell"]
        )
        _assert_no_amr_repairs(primitive_user)
        _assert_no_amr_repairs(conserved_user)
        assert np.max(primitive_user["crs_d_err"]) < 1.0e-12

        # With no coarse/fine transfer, selecting primitive AMR must not alter the
        # hyperbolic stage state before CornerE and CT.
        for column in (
            "time",
            "dt",
            "mass",
            "1-mom",
            "2-mom",
            "3-mom",
            "tot-E",
            "aam-D",
            "1-KE",
            "2-KE",
            "3-KE",
            "1-ME",
            "2-ME",
            "3-ME",
        ):
            np.testing.assert_allclose(
                primitive[column],
                conserved[column],
                rtol=0.0,
                atol=1.0e-12,
                err_msg=f"primitive restriction changed {column} without regridding",
            )
    finally:
        _cleanup()


def test_cgl_amr_conservative_refinement_gpu():
    try:
        _run(
            "cgl_amr_primitive_current_churn.athinput",
            "cgl_amr_gpu_conservative_refine",
            "time/nlim=1",
        )
        user = _user_history("cgl_amr_gpu_conservative_refine")
        mhd = _mhd_history("cgl_amr_gpu_conservative_refine")
        _assert_clean_user(user)
        _assert_no_amr_repairs(user)
        assert user["ncell"][-1] > user["ncell"][0]
        _assert_conserved(
            mhd,
            ("mass", "1-mom", "2-mom", "3-mom", "tot-E"),
        )
    finally:
        _cleanup()


def test_cgl_amr_current_churn_gpu():
    try:
        _run("cgl_amr_primitive_current_churn.athinput", "cgl_amr_gpu_churn")
        history = _user_history("cgl_amr_gpu_churn")
        mhd = _mhd_history("cgl_amr_gpu_churn")
        _assert_clean_user(history)
        _assert_no_amr_repairs(history)
        assert np.max(history["ncell"]) > history["ncell"][0]
        assert history["ncell"][-1] < np.max(history["ncell"])
        _assert_conserved(
            mhd,
            ("mass", "1-mom", "2-mom", "3-mom", "tot-E"),
        )
    finally:
        _cleanup()


def test_cgl_amr_restart_through_regrid_gpu():
    try:
        _run(
            "cgl_amr_primitive_current_churn_restart.athinput",
            "cgl_amr_gpu_restart_ref",
            "time/tlim=0.002",
            "time/nlim=-1",
        )
        reference = _user_history("cgl_amr_gpu_restart_ref")
        reference_state = _latest_state("cgl_amr_gpu_restart_ref")
        _assert_clean_user(reference)
        _assert_no_amr_repairs(reference)
        assert np.max(reference["ncell"]) > reference["ncell"][0]
        assert reference["ncell"][-1] < np.max(reference["ncell"])

        for split_cycle in (2, 3):
            split_basename = f"cgl_amr_gpu_restart_split_{split_cycle}"
            resumed_basename = f"cgl_amr_gpu_restart_resume_{split_cycle}"
            _run(
                "cgl_amr_primitive_current_churn_restart.athinput",
                split_basename,
                "time/tlim=0.002",
                f"time/nlim={split_cycle}",
            )
            restarts = sorted(Path("rst").rglob(f"{split_basename}*.rst"))
            assert restarts
            _restart(restarts[-1], resumed_basename, "time/nlim=-1")

            resumed = _user_history(resumed_basename)
            _assert_clean_user(resumed)
            _assert_no_amr_repairs(resumed)
            assert resumed["ncell"][-1] == reference["ncell"][-1]
            _assert_same_state(reference_state, _latest_state(resumed_basename))
    finally:
        _cleanup()


def test_cgl_amr_projection_stress_gpu():
    try:
        _run("cgl_amr_primitive_low_b.athinput", "cgl_amr_gpu_lowb")
        lowb = _user_history("cgl_amr_gpu_lowb")
        _assert_clean_user(lowb)
        assert lowb["amr_lowB"][-1] > 0.0

        _run("cgl_amr_primitive_firehose.athinput", "cgl_amr_gpu_firehose")
        firehose = _user_history("cgl_amr_gpu_firehose")
        _assert_clean_user(firehose)
        assert firehose["amr_fh"][-1] > 0.0

        _run("cgl_amr_primitive_mirror.athinput", "cgl_amr_gpu_mirror")
        mirror = _user_history("cgl_amr_gpu_mirror")
        _assert_clean_user(mirror)

        _run(
            "cgl_amr_primitive_mirror_slope_stress.athinput",
            "cgl_amr_gpu_mirror_slope",
        )
        mirror_slope = _user_history("cgl_amr_gpu_mirror_slope")
        _assert_clean_user(mirror_slope)
        assert mirror_slope["amr_mirr"][-1] > 0.0
        assert mirror_slope["amr_slp"][-1] > 0.0
    finally:
        _cleanup()


def test_cgl_amr_repair_accounting_tracks_only_transferred_states_gpu():
    try:
        _run(
            "cgl_amr_primitive_low_b.athinput",
            "cgl_amr_gpu_lowb_no_transfer",
            "problem/refine_levels=0",
            "time/nlim=1",
        )
        no_transfer = _user_history("cgl_amr_gpu_lowb_no_transfer")
        _assert_no_amr_repairs(no_transfer)

        _run(
            "cgl_amr_primitive_low_b.athinput",
            "cgl_amr_gpu_lowb_transfer",
            "time/nlim=1",
        )
        transferred = _user_history("cgl_amr_gpu_lowb_transfer")
        assert transferred["amr_cell"][-1] > 0.0
        assert transferred["amr_lowB"][-1] > 0.0

        _run(
            "cgl_amr_primitive_low_b_static.athinput",
            "cgl_amr_gpu_lowb_boundary",
        )
        boundary = _user_history("cgl_amr_gpu_lowb_boundary")
        assert np.all(boundary["ncell"] == boundary["ncell"][0])
        assert boundary["amr_cell"][-1] > boundary["amr_cell"][0]
        assert boundary["amr_lowB"][-1] > boundary["amr_lowB"][0]
    finally:
        _cleanup()


def test_cgl_amr_repair_counters_survive_restart_gpu():
    try:
        for label, per_rank in (("shared", False), ("rank_local", True)):
            partial_basename = f"cgl_amr_gpu_counter_{label}_partial"
            resumed_basename = f"cgl_amr_gpu_counter_{label}_resumed"
            flags = ["time/nlim=2", "output2/dcycle=2"]
            if per_rank:
                flags.append("output2/single_file_per_rank=true")
            _run(
                "cgl_amr_primitive_low_b.athinput",
                partial_basename,
                *flags,
            )
            partial = _user_history(partial_basename)
            assert partial["amr_lowB"][-1] > 0.0
            restarts = sorted(Path("rst").rglob(f"{partial_basename}*.rst"))
            assert restarts
            assert b"cgl_amr_repair_restart_version" in (
                restarts[-1].read_bytes()[:40000]
            )

            _restart(
                restarts[-1],
                resumed_basename,
                "time/nlim=2",
            )
            resumed = _user_history(resumed_basename)
            for column in AMR_REPAIR_COLUMNS:
                assert resumed[column][-1] == partial[column][-1]
    finally:
        _cleanup()


def test_cgl_lf_amr_primitive_churn_gpu():
    try:
        _run("cgl_lf_amr_primitive_churn.athinput", "cgl_amr_gpu_lf_churn")
        user = _user_history("cgl_amr_gpu_lf_churn")
        mhd = _mhd_history("cgl_amr_gpu_lf_churn")
        _assert_clean_user(user)
        _assert_clean_lf(mhd)
        _assert_no_amr_repairs(user)
        transitions = np.diff(user["ncell"])
        assert np.any(transitions > 0.0)
        assert np.any(transitions < 0.0)
        _assert_conserved(
            mhd,
            ("mass", "1-mom", "2-mom", "3-mom", "tot-E"),
        )
    finally:
        _cleanup()


def test_cgl_lf_amr_3d_churn_gpu():
    try:
        _run("cgl_lf_amr_3d_current.athinput", "cgl_amr_gpu_lf_3d_churn")
        user = _user_history("cgl_amr_gpu_lf_3d_churn")
        mhd = _mhd_history("cgl_amr_gpu_lf_3d_churn")
        _assert_clean_user(user, max_ndiv=1.0e-12)
        _assert_clean_lf(mhd)
        _assert_no_amr_repairs(user)
        assert user["ncell"][0] == 13824.0
        assert np.max(user["ncell"]) == 110592.0
        assert user["ncell"][-1] == user["ncell"][0]
        assert np.count_nonzero(user["ncell"] == np.max(user["ncell"])) >= 2
        assert np.count_nonzero(user["ncell"] == user["ncell"][0]) >= 4
        divb_paths = sorted(Path("bin").glob("*divb_resize*.bin"))
        assert len(divb_paths) >= 2
        divb_states = [bin_convert.read_binary(str(path)) for path in divb_paths]
        assert divb_states[0]["n_mbs"] == 27
        assert max(state["n_mbs"] for state in divb_states) > 27
        for state in divb_states:
            assert np.all(np.isfinite(state["mb_data"]["divb"]))
        for column in ("lf_mirror", "lf_firehs", "lf_hwproj"):
            assert mhd[column][-1] == 0.0
        _assert_conserved(
            mhd,
            ("mass", "1-mom", "2-mom", "3-mom", "tot-E"),
        )
    finally:
        _cleanup()


def test_cgl_amr_3d_and_passive_gpu():
    try:
        _run(
            "cgl_amr_paper_oblique_wave_static.athinput",
            "cgl_amr_gpu_wave_static",
        )

        _run("cgl_amr_primitive_3d_current.athinput", "cgl_amr_gpu_3d")
        three_d = _user_history("cgl_amr_gpu_3d")
        _assert_clean_user(three_d, max_ndiv=1.0e-9)
        assert np.max(three_d["crs_d_err"]) < 1.0e-12
        assert np.max(three_d["ncell"]) > three_d["ncell"][0]

        _run(
            "cgl_amr_passive_primitive_current.athinput",
            "cgl_amr_gpu_passive",
            "time/nlim=1",
        )
        passive = _user_history("cgl_amr_gpu_passive")
        passive_mhd = _mhd_history("cgl_amr_gpu_passive")
        _assert_clean_user(passive)
        _assert_no_amr_repairs(passive)
        assert np.max(passive["ncell"]) > passive["ncell"][0]
        _assert_conserved(
            passive_mhd,
            ("mass", "1-mom", "2-mom", "3-mom", "tot-E", "scal-0"),
        )
    finally:
        _cleanup()
