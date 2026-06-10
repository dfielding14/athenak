"""CPU regressions for the CGL Landau-fluid closure and CGL FOFC path."""

import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import struct
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

import test_suite.testutils as testutils
from test_suite.turb.test_turb_driving_cpu import read_force_blocks


INPUT_ROOT = "../../../inputs/tests"
PAPER_INPUT = "../../../inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"
PAPER_PASSIVE_INPUT = (
    "../../../inputs/cgl_lf_paper/cgl_lf_paper_smoke_passive_beta10.athinput"
)
PAPER_PRODUCTION_INPUT_ROOT = Path("../../../inputs/cgl_lf_paper")
PAPER_STAGE_I_MANIFEST = PAPER_PRODUCTION_INPUT_ROOT / "mks24_stage_i_manifest.json"
PAPER_WORKFLOW_PATH = Path("../../../scripts/cgl_lf_workflow.py")
PAPER_ANALYZER_PATH = Path("../../../scripts/analyze_cgl_lf_paper.py")
PAPER_STAGE_I_TOOL = Path("../../../scripts/frontier/cgl_lf_stage_i.py")


def _run(input_name, basename, *flags):
    testutils.run(
        f"{INPUT_ROOT}/{input_name}",
        [f"job/basename={basename}", *flags],
    )


def _run_paper(basename, *flags):
    testutils.run(PAPER_INPUT, [f"job/basename={basename}", *flags])


def _run_paper_passive(basename, *flags):
    testutils.run(PAPER_PASSIVE_INPUT, [f"job/basename={basename}", *flags])


def _cleanup():
    testutils.cleanup()
    for path in Path(".").glob("cgl_*.hst"):
        path.unlink()
    shutil.rmtree("bin", ignore_errors=True)


def _tab(basename):
    return testutils.athena_read.tab(f"tab/{basename}.mhd_w.00001.tab")


def _final_tab(basename):
    paths = sorted(Path("tab").glob(f"{basename}.mhd_w.*.tab"))
    assert paths, f"no table output found for {basename}"
    return testutils.athena_read.tab(str(paths[-1]))


def _final_variable_tab(basename, variable):
    paths = sorted(Path("tab").glob(f"{basename}.{variable}.*.tab"))
    assert paths, f"no {variable} table output found for {basename}"
    return testutils.athena_read.tab(str(paths[-1]))


def _final_force_binary(basename):
    paths = sorted(Path("bin").glob(f"{basename}.force_bin.*.bin"))
    assert paths, f"no full-field force output found for {basename}"
    return read_force_blocks(paths[-1])


def _force_fourier_components(output):
    """Return one fixed-grid force snapshot and its physical Fourier wavevectors."""

    assert len(output["blocks"]) == 1
    block = output["blocks"][0]
    force = np.asarray(block["force"], dtype=float)
    _, nz, ny, nx = force.shape
    limits = block["limits"]
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=(limits[1] - limits[0]) / nx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=(limits[3] - limits[2]) / ny)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=(limits[5] - limits[4]) / nz)
    kz_grid, ky_grid, kx_grid = np.meshgrid(kz, ky, kx, indexing="ij")
    fourier = np.fft.fftn(force, axes=(1, 2, 3))
    return force, fourier, kx_grid, ky_grid, kz_grid


def _assert_clean_lf_history(history):
    assert history["lf_nstage"][-1] > 0.0
    assert history["lf_dfloor"][-1] == 0.0
    assert history["lf_pfloor"][-1] == 0.0
    assert history["lf_nonfin"][-1] == 0.0
    assert history["lf_nonpos"][-1] == 0.0
    assert history["lf_hardbd"][-1] == 0.0
    assert history["lf_qface"][-1] > 0.0
    for name in ("lf_qprcap", "lf_qpr10", "lf_qpecap", "lf_qpe10"):
        assert 0.0 <= history[name][-1] <= history["lf_qface"][-1]
    for name in ("lf_qprwrk", "lf_qpewrk", "lf_cpwrk", "lf_cawrk"):
        if name in history:
            assert np.all(np.isfinite(history[name]))


def _assert_restarted_lf_diagnostics(reference, resumed):
    for column in (
        "lf_nstage",
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_mirror",
        "lf_firehs",
        "lf_hardbd",
        "lf_hwproj",
        "lf_qface",
        "lf_qprcap",
        "lf_qpr10",
        "lf_qpecap",
        "lf_qpe10",
        "lf_qprwrk",
        "lf_qpewrk",
        "lf_cpwrk",
        "lf_cawrk",
    ):
        assert np.isclose(
            reference[column][-1],
            resumed[column][-1],
            rtol=1.0e-12,
            atol=1.0e-14,
        )


def test_cgl_lf_quantitative_decay_and_diagnostics():
    try:
        _run("cgl_lf_decay.athinput", "cgl_ci_decay")
        history = testutils.athena_read.hst("cgl_ci_decay.mhd.hst")
        _assert_clean_lf_history(history)
        assert history["lf_qprwrk"][-1] > history["lf_qprwrk"][0]
        assert abs(history["lf_qpewrk"][-1]) < 1.0e-6 * history["lf_qprwrk"][-1]
    finally:
        _cleanup()


def test_cgl_lf_limiter_occupancy_remains_admissible():
    try:
        _run("cgl_lf_limiter.athinput", "cgl_ci_limiter")
        history = testutils.athena_read.hst("cgl_ci_limiter.mhd.hst")
        assert history["lf_mirror"][-1] > 0.0
        assert history["lf_nonfin"][-1] == 0.0
        assert history["lf_nonpos"][-1] == 0.0
        assert history["lf_hardbd"][-1] == 0.0
        assert np.all(np.diff(history["lf_nstage"]) >= 0.0)
    finally:
        _cleanup()


def test_cgl_lf_heat_flux_cap_face_activity_is_reported():
    try:
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_flux_cap",
            "time/evolution=kinematic",
            "time/nlim=1",
            "time/tlim=1.0e-5",
            "mhd/rsolver=advect",
            "mhd/lf_k_parallel=1.0e-2",
            "problem/test_mode=flux_limiter",
            "problem/amp=0.5",
        )
        history = testutils.athena_read.hst("cgl_ci_flux_cap.mhd.hst")
        _assert_clean_lf_history(history)
        assert history["lf_qprcap"][-1] > 0.0
        assert history["lf_qpecap"][-1] > 0.0
        assert history["lf_qpr10"][-1] > 0.0
        assert history["lf_qpe10"][-1] > 0.0
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                "--lf-history",
                "cgl_ci_flux_cap.mhd.hst",
                "--time-start",
                "0.0",
                "--time-end",
                "1.0",
                "--output-dir",
                "cgl_ci_flux_cap_analysis",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        diagnostics = json.loads(
            Path("cgl_ci_flux_cap_analysis/diagnostics.json").read_text()
        )
        cap = diagnostics["lf_histories"][0]["heat_flux_cap_fractions"]
        assert cap["parallel_over_10"] > 0.0
        assert cap["perpendicular_over_10"] > 0.0
        work = diagnostics["lf_histories"][0]["applied_heat_flux_work"]
        assert np.isfinite(work["total"])
        assert work["total"] > 0.0
    finally:
        shutil.rmtree("cgl_ci_flux_cap_analysis", ignore_errors=True)
        _cleanup()


def test_cgl_lf_heat_flux_cap_fractions_are_meshblock_layout_independent():
    common_flags = (
        "time/nlim=1",
        "time/tlim=2.0e-5",
        "mhd/lf_k_parallel=1.0e-8",
        "output2/dt=1.0",
        "output3/dt=1.0",
    )
    try:
        _run_paper("cgl_ci_cap_single_block", *common_flags)
        _run_paper("cgl_ci_cap_split_blocks", *common_flags, "meshblock/nx1=4")
        single = testutils.athena_read.hst("cgl_ci_cap_single_block.mhd.hst")
        split = testutils.athena_read.hst("cgl_ci_cap_split_blocks.mhd.hst")
        for column in (
            "lf_qface",
            "lf_qprcap",
            "lf_qpr10",
            "lf_qpecap",
            "lf_qpe10",
        ):
            assert single[column][-1] == split[column][-1]
        for column in ("lf_qprwrk", "lf_qpewrk"):
            assert np.isclose(
                single[column][-1], split[column][-1], rtol=5.0e-12, atol=1.0e-20
            )
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_low_field_faces_disable_transport_cleanly():
    try:
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_low_field",
            "time/evolution=kinematic",
            "time/nlim=1",
            "time/tlim=1.0e-5",
            "mhd/rsolver=advect",
            "problem/test_mode=low_field",
            "problem/b0=5.0e-11",
            "problem/amp=0.5",
        )
        history = testutils.athena_read.hst("cgl_ci_low_field.mhd.hst")
        assert history["lf_nstage"][-1] > 0.0
        assert history["lf_qface"][-1] == 0.0
        assert history["lf_qprwrk"][-1] == 0.0
        assert history["lf_qpewrk"][-1] == 0.0
        for column in (
            "lf_dfloor",
            "lf_pfloor",
            "lf_nonfin",
            "lf_nonpos",
            "lf_hardbd",
        ):
            assert history[column][-1] == 0.0
    finally:
        _cleanup()


def test_cgl_lf_background_collision_advances_one_physical_timestep():
    try:
        _run("cgl_lf_collision_relaxation.athinput", "cgl_ci_collision_relaxation")
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_collision_relaxation.mhd.hst")
        )
    finally:
        _cleanup()


def test_cgl_lf_firehose_threshold_policies_are_distinct():
    try:
        _run("cgl_lf_firehose_policy.athinput", "cgl_ci_firehose_oblique")
        _run(
            "cgl_lf_firehose_policy.athinput",
            "cgl_ci_firehose_parallel",
            "mhd/cgl_firehose_threshold=parallel",
        )
        oblique = testutils.athena_read.hst("cgl_ci_firehose_oblique.mhd.hst")
        parallel = testutils.athena_read.hst("cgl_ci_firehose_parallel.mhd.hst")
        assert oblique["lf_firehs"][-1] > 0.0
        assert parallel["lf_firehs"][-1] == 0.0
        _assert_clean_lf_history(oblique)
        _assert_clean_lf_history(parallel)
    finally:
        _cleanup()


def test_cgl_lf_hardwall_projects_to_selected_firehose_threshold():
    try:
        _run(
            "cgl_lf_firehose_policy.athinput",
            "cgl_ci_firehose_hardwall",
            "mhd/limiter_hardwall=true",
        )
        history = testutils.athena_read.hst("cgl_ci_firehose_hardwall.mhd.hst")
        assert history["lf_hwproj"][-1] > 0.0
        assert history["lf_hardbd"][-1] == 0.0
        assert history["lf_nonfin"][-1] == 0.0
        assert history["lf_nonpos"][-1] == 0.0
    finally:
        _cleanup()


def test_cgl_lf_strict_hard_bound_is_reported_without_backup_correction():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_firehose_policy.athinput",
        "mhd/cgl_firehose_threshold=parallel",
        "problem/ppar0=3.0",
        "problem/pperp0=1.0",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "strict admissibility failed" in result.stdout
    assert "hard_bound=" in result.stdout


def test_cgl_lf_explicit_reference_agrees_with_capped_sts():
    try:
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_sts_reference",
            "time/sts_max_dt_ratio=1.0",
            "time/nlim=-1",
        )
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_explicit_reference",
            "mhd/cgl_heat_flux_integrator=explicit",
            "time/sts_integrator=none",
            "time/nlim=-1",
        )
        sts = _tab("cgl_ci_sts_reference")
        explicit = _tab("cgl_ci_explicit_reference")
        fields = set(sts).intersection(explicit) - {"time", "cycle"}
        max_difference = max(np.max(np.abs(sts[field] - explicit[field]))
                             for field in fields)
        assert max_difference < 1.0e-3
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_sts_reference.mhd.hst")
        )
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_explicit_reference.mhd.hst")
        )
    finally:
        _cleanup()


def test_cgl_lf_explicit_reference_executes_both_half_sweeps():
    try:
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_explicit_two_sweeps",
            "mhd/cgl_heat_flux_integrator=explicit",
            "time/sts_integrator=none",
            "time/nlim=1",
            "time/tlim=1.0e-5",
        )
        history = testutils.athena_read.hst(
            "cgl_ci_explicit_two_sweeps.mhd.hst"
        )
        _assert_clean_lf_history(history)
        assert history["lf_nstage"][-1] == 128.0
        assert history["lf_qface"][-1] == 128.0
    finally:
        _cleanup()


def test_cgl_lf_explicit_reference_finite_collision_split():
    try:
        common = ("time/nlim=-1", "mhd/nu_coll=1.0")
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_sts_collision",
            "time/sts_max_dt_ratio=1.0",
            *common,
        )
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_explicit_collision",
            "mhd/cgl_heat_flux_integrator=explicit",
            "time/sts_integrator=none",
            *common,
        )
        sts = _tab("cgl_ci_sts_collision")
        explicit = _tab("cgl_ci_explicit_collision")
        assert np.max(np.abs(sts["eint"] - explicit["eint"])) < 1.0e-3
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_sts_collision.mhd.hst")
        )
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_explicit_collision.mhd.hst")
        )
    finally:
        _cleanup()


def test_cgl_fofc_live_flux_mutation():
    try:
        _run("cgl_fofc.athinput", "cgl_ci_fofc")
    finally:
        _cleanup()


def test_cgl_lf_amr_conserved_prolongation_stays_admissible():
    try:
        _run("cgl_lf_amr_2d.athinput", "cgl_ci_amr")
        user = testutils.athena_read.hst("cgl_ci_amr.user.hst")
        mhd = testutils.athena_read.hst("cgl_ci_amr.mhd.hst")
        assert user["ncell"][-1] > user["ncell"][0]
        assert np.max(user["max_ndiv"]) < 1.0e-12
        assert np.max(user["bad_state"]) == 0.0
        assert np.all(np.isfinite(user["abs_anis"]))
        _assert_clean_lf_history(mhd)
        assert abs(mhd["lf_qprwrk"][-1]) + abs(mhd["lf_qpewrk"][-1]) > 0.0
        assert abs(mhd["lf_cpwrk"][-1]) + abs(mhd["lf_cawrk"][-1]) > 0.0
        energy_scale = max(abs(mhd["tot-E"][0]), 1.0e-30)
        energy_residual = abs(mhd["tot-E"][-1] - mhd["tot-E"][0]) / energy_scale
        assert energy_residual < 5.0e-3
    finally:
        _cleanup()


def test_cgl_lf_restart_marker_round_trips_terminal_time():
    try:
        basename = "cgl_ci_restart_precision"
        _run(
            "cgl_lf_restart.athinput",
            basename,
            "time/nlim=1",
            "time/cfl_number=0.371234567890123",
        )
        restart_paths = sorted(Path("rst").glob(f"{basename}.*.rst"))
        assert restart_paths, "partial run did not write a restart checkpoint"
        marker = re.search(
            rb"(?m)^restart_time\s*=\s*(\S+)",
            restart_paths[-1].read_bytes()[:40000],
        )
        assert marker is not None, "restart checkpoint did not contain time/restart_time"
        terminal_time = testutils.athena_read.hst(f"{basename}.mhd.hst")["time"][-1]
        assert float(marker.group(1)) == terminal_time
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_restart_preserves_final_state_and_admissibility():
    try:
        _run("cgl_lf_restart.athinput", "cgl_ci_restart_reference")
        # Stop on a shared cycle so the checkpoint does not change timesteps.
        _run(
            "cgl_lf_restart.athinput",
            "cgl_ci_restart_partial",
            "time/nlim=1",
        )
        restart_paths = sorted(Path("rst").glob("cgl_ci_restart_partial.*.rst"))
        assert restart_paths, "partial run did not write a restart checkpoint"
        assert b"restart_time" in restart_paths[-1].read_bytes()[:40000]
        command = [
            "./athena",
            "-r",
            str(restart_paths[-1]),
            "job/basename=cgl_ci_restart_resumed",
            "time/nlim=-1",
        ]
        assert testutils.run_command(command)

        reference = _final_tab("cgl_ci_restart_reference")
        resumed = _final_tab("cgl_ci_restart_resumed")
        fields = set(reference).intersection(resumed) - {"time", "cycle"}
        maximum = max(
            np.max(np.abs(reference[field] - resumed[field]))
            for field in fields
        )
        assert maximum < 1.0e-12
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_restart_resumed.mhd.hst")
        )
        reference_history = testutils.athena_read.hst(
            "cgl_ci_restart_reference.mhd.hst"
        )
        resumed_history = testutils.athena_read.hst("cgl_ci_restart_resumed.mhd.hst")
        _assert_restarted_lf_diagnostics(reference_history, resumed_history)
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_restart_with_finite_collision_preserves_corrected_split():
    try:
        collision = "mhd/nu_coll=1.0"
        _run(
            "cgl_lf_restart.athinput",
            "cgl_ci_restart_collision_reference",
            collision,
        )
        _run(
            "cgl_lf_restart.athinput",
            "cgl_ci_restart_collision_partial",
            collision,
            "time/nlim=1",
        )
        restart_paths = sorted(
            Path("rst").glob("cgl_ci_restart_collision_partial.*.rst")
        )
        assert restart_paths, "finite-collision partial run did not write a checkpoint"
        command = [
            "./athena",
            "-r",
            str(restart_paths[-1]),
            "job/basename=cgl_ci_restart_collision_resumed",
            "time/nlim=-1",
        ]
        assert testutils.run_command(command)

        reference = _final_tab("cgl_ci_restart_collision_reference")
        resumed = _final_tab("cgl_ci_restart_collision_resumed")
        fields = set(reference).intersection(resumed) - {"time", "cycle"}
        maximum = max(
            np.max(np.abs(reference[field] - resumed[field]))
            for field in fields
        )
        assert maximum < 1.0e-12
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_restart_collision_resumed.mhd.hst")
        )
        reference_history = testutils.athena_read.hst(
            "cgl_ci_restart_collision_reference.mhd.hst"
        )
        resumed_history = testutils.athena_read.hst(
            "cgl_ci_restart_collision_resumed.mhd.hst"
        )
        _assert_restarted_lf_diagnostics(reference_history, resumed_history)
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_invalid_integrator_is_rejected():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_decay.athinput",
        "mhd/cgl_heat_flux_integrator=invalid",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "cgl_heat_flux_integrator" in result.stdout


def test_cgl_lf_invalid_firehose_threshold_is_rejected():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_firehose_policy.athinput",
        "mhd/cgl_firehose_threshold=invalid",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "cgl_firehose_threshold" in result.stdout


def test_cgl_lf_hardwall_requires_instability_limiter():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_decay.athinput",
        "mhd/limiter_hardwall=true",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "limiter_hardwall requires" in result.stdout


def test_cgl_lf_paper_active_alfvenic_smoke_injects_energy_without_parallel_force():
    try:
        _run_paper("cgl_ci_paper_alfvenic", "time/nlim=1")
        mhd = testutils.athena_read.hst("cgl_ci_paper_alfvenic.mhd.hst")
        user = testutils.athena_read.hst("cgl_ci_paper_alfvenic.user.hst")
        primitive = _final_variable_tab("cgl_ci_paper_alfvenic", "mhd_w_bcc")
        _assert_clean_lf_history(mhd)
        assert user["p_parallel"][0] == user["p_perp"][0]
        assert user["force_prp2"][-1] > 0.0
        assert user["force_prl2"][-1] == 0.0
        assert np.isclose(user["mass"][0], user["volume"][0])
        assert np.isclose(user["b2"][0], user["volume"][0])
        assert np.isclose(user["b4"][0], user["volume"][0])
        assert np.isclose(user["beta"][0] / user["volume"][0], 10.0)
        assert user["delta_p"][0] == 0.0
        assert user["mirror_vol"][-1] == 0.0
        assert user["fire_vol"][-1] == 0.0
        assert user["hard_vol"][-1] == 0.0
        assert user["force_pwr"][-1] > 0.0
        assert user["force_work"][-1] > user["force_work"][0]
        assert "p_perp" in primitive
        assert np.all(np.isfinite(primitive["p_perp"]))
        assert mhd["tot-E"][-1] > mhd["tot-E"][0]
        measured_work = mhd["tot-E"][-1] - mhd["tot-E"][0]
        applied_work = user["force_work"][-1] - user["force_work"][0]
        assert np.isclose(applied_work, measured_work, rtol=1.0e-10, atol=1.0e-12)
        force, fourier, kx, ky, kz = _force_fourier_components(
            _final_force_binary("cgl_ci_paper_alfvenic")
        )
        assert np.max(np.abs(force[2])) == 0.0
        perpendicular_divergence = kx * fourier[0] + ky * fourier[1]
        perpendicular_scale = np.abs(kx * fourier[0]) + np.abs(ky * fourier[1])
        retained_kz = (
            (np.abs(kz) > 0.0)
            & (perpendicular_scale > 1.0e-10 * np.max(perpendicular_scale))
        )
        assert np.any(retained_kz)
        assert np.max(np.abs(perpendicular_divergence[retained_kz])) < (
            2.0e-6 * np.max(perpendicular_scale[retained_kz])
        )
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_random_forcing_leaves_cartesian_amplitudes_unprojected():
    try:
        _run_paper(
            "cgl_ci_paper_random",
            "time/nlim=1",
            "turb_driving/driving_type=0",
            "turb_driving/projection_policy=mks24_random_unprojected",
            "turb_driving/rseed=314159",
        )
        force, fourier, kx, ky, kz = _force_fourier_components(
            _final_force_binary("cgl_ci_paper_random")
        )
        assert np.max(np.abs(force[2])) > 0.0
        divergence = kx * fourier[0] + ky * fourier[1] + kz * fourier[2]
        divergence_scale = (
            np.abs(kx * fourier[0])
            + np.abs(ky * fourier[1])
            + np.abs(kz * fourier[2])
        )
        active = divergence_scale > 1.0e-10 * np.max(divergence_scale)
        assert np.any(active)
        assert np.max(np.abs(divergence[active])) > (
            1.0e-3 * np.max(divergence_scale[active])
        )
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_fixed_edot_normalization_sets_first_cycle_amplitude():
    try:
        _run_paper("cgl_ci_ou_correlated", "time/nlim=1")
        _run_paper("cgl_ci_ou_white", "time/nlim=1", "turb_driving/tcorr=0.0")
        correlated = _final_variable_tab("cgl_ci_ou_correlated", "turb_force")
        white = _final_variable_tab("cgl_ci_ou_white", "turb_force")
        for field in ("force1", "force2", "force3"):
            assert np.allclose(
                correlated[field], white[field], rtol=5.0e-12, atol=5.0e-14
            )
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_multicycle_forcing_work_follows_rk_state_recurrence():
    try:
        _run_paper("cgl_ci_paper_multicycle_work", "time/nlim=4")
        mhd = testutils.athena_read.hst("cgl_ci_paper_multicycle_work.mhd.hst")
        user = testutils.athena_read.hst("cgl_ci_paper_multicycle_work.user.hst")
        energy_delta = mhd["tot-E"][-1] - mhd["tot-E"][0]
        applied_work = user["force_work"][-1] - user["force_work"][0]
        assert applied_work > 0.0
        assert np.isclose(applied_work, energy_delta, rtol=1.0e-10, atol=1.0e-12)
        assert np.all(np.isfinite(mhd["lf_cpwrk"]))
        assert np.all(np.isfinite(mhd["lf_cawrk"]))
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_forcing_seed_is_deterministic_and_selectable():
    try:
        _run_paper("cgl_ci_seed_a", "time/nlim=1")
        _run_paper("cgl_ci_seed_b", "time/nlim=1")
        _run_paper("cgl_ci_seed_c", "time/nlim=1", "turb_driving/rseed=42")
        force_a = _final_variable_tab("cgl_ci_seed_a", "turb_force")
        force_b = _final_variable_tab("cgl_ci_seed_b", "turb_force")
        force_c = _final_variable_tab("cgl_ci_seed_c", "turb_force")
        for field in ("force1", "force2", "force3"):
            assert np.array_equal(force_a[field], force_b[field])
        changed = max(
            np.max(np.abs(force_a[field] - force_c[field]))
            for field in ("force1", "force2", "force3")
        )
        assert changed > 1.0e-12
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_forcing_restart_preserves_rng_and_force_state():
    try:
        _run_paper("cgl_ci_paper_restart_reference")
        _run_paper(
            "cgl_ci_paper_restart_partial",
            "time/nlim=1",
            "output4/single_file_per_rank=true",
        )
        restart_paths = sorted(
            Path("rst/rank_00000000").glob("cgl_ci_paper_restart_partial.*.rst")
        )
        assert restart_paths, "paper smoke partial run did not write a checkpoint"
        incompatible = subprocess.run(
            [
                "./athena",
                "-r",
                str(restart_paths[-1]),
                "turb_driving/projection_policy=solenoidal_compressive",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert incompatible.returncode != 0
        assert "configuration differs for 'projection_policy'" in (
            incompatible.stdout + incompatible.stderr
        )
        command = [
            "./athena",
            "-r",
            str(restart_paths[-1]),
            "job/basename=cgl_ci_paper_restart_resumed",
            "time/nlim=-1",
        ]
        assert testutils.run_command(command)

        for variable in ("mhd_w_bcc", "turb_force"):
            reference = _final_variable_tab("cgl_ci_paper_restart_reference", variable)
            resumed = _final_variable_tab("cgl_ci_paper_restart_resumed", variable)
            fields = set(reference).intersection(resumed) - {"time", "cycle"}
            maximum = max(
                np.max(np.abs(reference[field] - resumed[field]))
                for field in fields
            )
            assert maximum < 1.0e-12
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_paper_restart_resumed.mhd.hst")
        )
        reference_user = testutils.athena_read.hst(
            "cgl_ci_paper_restart_reference.user.hst"
        )
        resumed_user = testutils.athena_read.hst("cgl_ci_paper_restart_resumed.user.hst")
        assert np.isclose(
            reference_user["force_work"][-1],
            resumed_user["force_work"][-1],
            rtol=1.0e-12,
            atol=1.0e-14,
        )
        reference_mhd = testutils.athena_read.hst(
            "cgl_ci_paper_restart_reference.mhd.hst"
        )
        resumed_mhd = testutils.athena_read.hst("cgl_ci_paper_restart_resumed.mhd.hst")
        _assert_restarted_lf_diagnostics(reference_mhd, resumed_mhd)
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_passive_delta_has_no_anisotropic_flow_feedback():
    try:
        _run_paper_passive("cgl_ci_passive_iso", "time/nlim=4")
        _run_paper_passive(
            "cgl_ci_passive_aniso",
            "time/nlim=4",
            "problem/p_parallel0=5.2",
            "problem/p_perp0=4.9",
        )
        isotropic = _final_variable_tab("cgl_ci_passive_iso", "mhd_w_bcc")
        anisotropic = _final_variable_tab("cgl_ci_passive_aniso", "mhd_w_bcc")
        for field in ("dens", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3"):
            assert np.max(np.abs(isotropic[field] - anisotropic[field])) < 1.0e-12
        assert np.max(np.abs(isotropic["eint"] - anisotropic["eint"])) > 1.0e-4
        passive_iso = testutils.athena_read.hst("cgl_ci_passive_iso.mhd.hst")
        passive_aniso = testutils.athena_read.hst("cgl_ci_passive_aniso.mhd.hst")
        _assert_clean_lf_history(passive_iso)
        _assert_clean_lf_history(passive_aniso)
        assert np.all(passive_iso["lf_cpwrk"] == 0.0)
        assert np.all(passive_iso["lf_cawrk"] == 0.0)
        assert np.all(passive_aniso["lf_cpwrk"] == 0.0)
        assert np.all(passive_aniso["lf_cawrk"] == 0.0)
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_passive_delta_must_match_eos_mode():
    command = [
        "./athena",
        "-i",
        PAPER_INPUT,
        "problem/passive_delta=true",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "passive_delta must match" in result.stdout


def test_cgl_lf_paper_rejects_unsupported_forcing_mode():
    command = [
        "./athena",
        "-i",
        PAPER_INPUT,
        "turb_driving/driving_type=2",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "driving_type must be 0" in result.stdout


def test_cgl_lf_paper_alfvenic_policy_rejects_compressive_blend():
    command = [
        "./athena",
        "-i",
        PAPER_INPUT,
        "turb_driving/sol_fraction=0.5",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "mks24_alfvenic_perpendicular requires sol_fraction = 1" in result.stdout


def test_cgl_lf_paper_physical_forcing_shell_requires_positive_unit():
    command = [
        "./athena",
        "-i",
        PAPER_INPUT,
        "turb_driving/k_shell_unit=0.0",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "k_shell_unit must be positive" in result.stdout


def test_cgl_lf_paper_production_inputs_explicitly_use_rank_local_io():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_workflow_test", PAPER_WORKFLOW_PATH
    )
    assert spec is not None and spec.loader is not None
    workflow = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = workflow
    spec.loader.exec_module(workflow)
    args = workflow.parser().parse_args([
        "paper-analyze",
        "--reference-curves",
        "fig2.json",
        "--reference-curves",
        "fig13.json",
        "--allow-partial-reference-cases",
    ])
    assert args.reference_curves == ["fig2.json", "fig13.json"]
    assert args.allow_partial_reference_cases
    all_paper_inputs = sorted(PAPER_PRODUCTION_INPUT_ROOT.glob("*.athinput"))
    assert len(all_paper_inputs) == 20
    for input_path in all_paper_inputs:
        source = input_path.read_text()
        turbulence = source.split("<turb_driving>", 1)[1].split("<", 1)[0]
        assert "spectrum = power_law" in turbulence
        assert "physical_k_shell = true" in turbulence
        assert "k_shell_unit = 3.141592653589793" in turbulence
        assert "nlow = 1" in turbulence
        assert "nhigh = 3" in turbulence
        assert "isotropic_power_spectrum = true" in turbulence
        assert "expo = 2.0" in turbulence
        expected_policy = (
            "mks24_random_unprojected"
            if "driving_type = 0" in turbulence
            else "mks24_alfvenic_perpendicular"
        )
        assert f"projection_policy = {expected_policy}" in turbulence
        if expected_policy == "mks24_alfvenic_perpendicular":
            assert "sol_fraction = 1.0" in turbulence
    input_paths = sorted(
        PAPER_PRODUCTION_INPUT_ROOT.glob("cgl_lf_paper_standard_*.athinput")
    )
    input_paths += sorted(
        PAPER_PRODUCTION_INPUT_ROOT.glob("cgl_lf_paper_nulim_*.athinput")
    )
    input_paths += sorted(
        PAPER_PRODUCTION_INPUT_ROOT.glob("cgl_lf_paper_heat_flux_*.athinput")
    )
    input_paths += sorted(
        PAPER_PRODUCTION_INPUT_ROOT.glob("cgl_lf_paper_compressive_*.athinput")
    )
    input_paths += sorted(
        PAPER_PRODUCTION_INPUT_ROOT.glob(
            "cgl_lf_paper_scale_separation_*.athinput"
        )
    )
    assert len(input_paths) == 18
    standard_input_names = {
        path.name
        for path in PAPER_PRODUCTION_INPUT_ROOT.glob("cgl_lf_paper_standard_*.athinput")
    }
    workflow_standard_names = {
        Path(case.input_path).name
        for case in workflow.workflow_cases("paper-standard")
    }
    assert workflow_standard_names == standard_input_names - {
        "cgl_lf_paper_standard_active_alfvenic_beta1.athinput"
    }
    heat_flux_input_names = {
        path.name
        for path in PAPER_PRODUCTION_INPUT_ROOT.glob(
            "cgl_lf_paper_heat_flux_*.athinput"
        )
    }
    workflow_heat_flux_names = {
        Path(case.input_path).name
        for case in workflow.workflow_cases("paper-heat-flux")
    }
    assert workflow_heat_flux_names == heat_flux_input_names
    compressive_input_names = {
        path.name
        for path in PAPER_PRODUCTION_INPUT_ROOT.glob(
            "cgl_lf_paper_compressive_*.athinput"
        )
    }
    workflow_compressive_names = {
        Path(case.input_path).name
        for case in workflow.workflow_cases("paper-compressive")
    }
    assert workflow_compressive_names == compressive_input_names | {
        "cgl_lf_paper_standard_active_alfvenic_beta100.athinput",
        "cgl_lf_paper_standard_active_random_beta100.athinput",
        "cgl_lf_paper_standard_active_random_beta10.athinput",
        "cgl_lf_paper_standard_passive_random_beta100.athinput",
    }
    scale_separation_input_names = {
        path.name
        for path in PAPER_PRODUCTION_INPUT_ROOT.glob(
            "cgl_lf_paper_scale_separation_*.athinput"
        )
    }
    workflow_scale_separation_names = {
        Path(case.input_path).name
        for case in workflow.workflow_cases("paper-scale-separation")
    }
    assert workflow_scale_separation_names == scale_separation_input_names | {
        "cgl_lf_paper_standard_active_alfvenic_beta10.athinput"
    }
    stage_i_names = [
        Path(case.input_path).name
        for case in workflow.workflow_cases("paper-mks24-stage-i")
    ]
    stage_i_manifest = json.loads(PAPER_STAGE_I_MANIFEST.read_text())
    manifest_names = [
        Path(case["input"]).name for case in stage_i_manifest["cases"]
    ]
    assert len(stage_i_names) == 16
    assert len(set(stage_i_names)) == 16
    assert stage_i_manifest["authorization"]["mapped_unique_runs"] == 16
    assert stage_i_names == manifest_names
    assert "cgl_lf_paper_standard_active_alfvenic_beta1.athinput" not in stage_i_names
    assert "cgl_lf_paper_nulim_beta100_hardwall.athinput" not in stage_i_names
    assert set(stage_i_names) == {
        "cgl_lf_paper_standard_active_alfvenic_beta10.athinput",
        "cgl_lf_paper_standard_active_alfvenic_beta100.athinput",
        "cgl_lf_paper_standard_active_random_beta10.athinput",
        "cgl_lf_paper_standard_active_random_beta100.athinput",
        "cgl_lf_paper_standard_passive_alfvenic_beta10.athinput",
        "cgl_lf_paper_standard_passive_alfvenic_beta100.athinput",
        "cgl_lf_paper_standard_passive_random_beta10.athinput",
        "cgl_lf_paper_standard_passive_random_beta100.athinput",
        "cgl_lf_paper_compressive_active_random_beta1.athinput",
        "cgl_lf_paper_compressive_active_random_beta100_sonic.athinput",
        "cgl_lf_paper_heat_flux_beta10_strong.athinput",
        "cgl_lf_paper_heat_flux_beta10_weak.athinput",
        "cgl_lf_paper_nulim_beta100_20.athinput",
        "cgl_lf_paper_nulim_beta100_200.athinput",
        "cgl_lf_paper_scale_separation_beta10_nperp96.athinput",
        "cgl_lf_paper_scale_separation_beta10_nperp384.athinput",
    }
    hardwall_paths = {
        "cgl_lf_paper_standard_active_alfvenic_beta1.athinput",
        "cgl_lf_paper_standard_active_alfvenic_beta10.athinput",
        "cgl_lf_paper_standard_active_alfvenic_beta100.athinput",
        "cgl_lf_paper_standard_active_random_beta10.athinput",
        "cgl_lf_paper_standard_active_random_beta100.athinput",
        "cgl_lf_paper_standard_passive_alfvenic_beta10.athinput",
        "cgl_lf_paper_standard_passive_alfvenic_beta100.athinput",
        "cgl_lf_paper_standard_passive_random_beta10.athinput",
        "cgl_lf_paper_standard_passive_random_beta100.athinput",
        "cgl_lf_paper_nulim_beta100_hardwall.athinput",
        "cgl_lf_paper_heat_flux_beta10_strong.athinput",
        "cgl_lf_paper_heat_flux_beta10_weak.athinput",
        "cgl_lf_paper_compressive_active_random_beta1.athinput",
        "cgl_lf_paper_compressive_active_random_beta100_sonic.athinput",
        "cgl_lf_paper_scale_separation_beta10_nperp96.athinput",
        "cgl_lf_paper_scale_separation_beta10_nperp384.athinput",
    }
    for input_path in input_paths:
        source = input_path.read_text()
        for block in ("output2", "output3"):
            body = source.split(f"<{block}>", 1)[1].split("<", 1)[0]
            assert "single_file_per_rank = true" in body
        choices = workflow.model_choices(source, [])
        assert choices["time_integrator"] == "rk2"
        assert choices["time_sts_integrator"] == "rkl2"
        assert choices["time_sts_max_dt_ratio"] == "-1.0"
        assert choices["time_cfl_number"] == "0.3"
        assert choices["time_tlim"] == "10.0"
        assert choices["output2_file_type"] == "bin"
        assert choices["output2_dt"] == "0.25"
        assert choices["output2_single_file_per_rank"] == "true"
        assert choices["output3_file_type"] == "rst"
        assert choices["output3_dt"] == "1.0"
        assert choices["output3_single_file_per_rank"] == "true"
        expected_hardwall = "true" if input_path.name in hardwall_paths else "false"
        assert choices["limiter_hardwall"] == expected_hardwall
    scale_resolutions = {
        "cgl_lf_paper_scale_separation_beta10_nperp96.athinput": (
            "96", "96", "192"
        ),
        "cgl_lf_paper_scale_separation_beta10_nperp384.athinput": (
            "384", "384", "768"
        ),
    }
    for name, expected in scale_resolutions.items():
        source = (PAPER_PRODUCTION_INPUT_ROOT / name).read_text()
        choices = workflow.model_choices(source, [])
        assert (
            choices["mesh_nx1"], choices["mesh_nx2"], choices["mesh_nx3"]
        ) == expected
    compressive_parameters = {
        "cgl_lf_paper_compressive_active_random_beta1.athinput": ("1.0", "2.0"),
        "cgl_lf_paper_compressive_active_random_beta100_sonic.athinput": (
            "100.0", "0.2"
        ),
    }
    for name, (beta0, tcorr) in compressive_parameters.items():
        source = (PAPER_PRODUCTION_INPUT_ROOT / name).read_text()
        choices = workflow.model_choices(source, [])
        assert choices["beta0"] == beta0
        assert choices["forcing_mode"] == "isotropic_random"
        assert choices["forcing_tcorr"] == tcorr
    for workflow_name in (
        "paper-compressive", "paper-scale-separation", "paper-mks24-stage-i"
    ):
        output_dir = f"cgl_ci_denied_{workflow_name.replace('-', '_')}"
        denied = subprocess.run(
            [
                sys.executable,
                str(PAPER_WORKFLOW_PATH),
                workflow_name,
                "--output-dir",
                output_dir,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert denied.returncode != 0
        assert "--authorize-paper-execution" in denied.stderr
        assert not Path(output_dir).exists()


def test_cgl_lf_stage_i_acceptance_requires_clean_complete_segment(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    root = tmp_path / "root"
    output_dir = (
        root / "runs" / "mks24-stage-i" / stage_i.EXECUTION_EPOCH
        / "R16" / "s00" / "output"
    )
    manifest_dir = output_dir.parent / "manifest"
    (output_dir / "bin").mkdir(parents=True)
    (output_dir / "rst").mkdir()
    manifest_dir.mkdir()
    mhd_history = output_dir / "case.mhd.hst"
    (output_dir / "case.user.hst").write_text("# retained user history\n")
    snapshot = output_dir / "bin" / "case.00000.bin"
    restart = output_dir / "rst" / "case.00000.rst"
    restart.write_bytes(b"restart")

    def write_snapshot_time(time):
        snapshot.write_bytes(
            (
                "Athena binary output version=1.1\n"
                "  size of preheader=5\n"
                f"  time={time}\n"
                "  cycle=0\n"
                "  size of location=8\n"
                "  size of variable=4\n"
            ).encode()
        )

    def write_history(rows):
        mhd_history.write_text(
            "# [0]=time [1]=lf_dfloor [2]=lf_pfloor [3]=lf_nonfin "
            "[4]=lf_nonpos [5]=lf_hardbd [6]=lf_hwproj\n"
            + "\n".join(" ".join(str(value) for value in row) for row in rows)
            + "\n"
        )
        restart.write_text(
            f"<time>\nrestart_time = {rows[-1][0]}\n<par_end>\n"
        )

    manifest_path = manifest_dir / "prepared_run.json"
    manifest = {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "state": "submitted",
        "project_root": str(root),
        "job_id": "12345",
        "run": {
            "case_id": "R16",
            "case_name": "case",
            "segment": "s00",
        },
        "allocation": {
            "nodes": 1,
            "requested_walltime": "00:10:00",
            "reserved_node_hours": 1.0 / 6.0,
        },
        "command": {
            "executable_revision": "a" * 40,
            "executable_sha256": "b" * 64,
            "input_revision": "a" * 40,
            "input_sha256": "c" * 64,
            "input_file": "submitted_input.athinput",
            "overrides": ["time/tlim=2.0"],
            "time_tlim_target": 2.0,
        },
        "paths": {"output_dir": str(output_dir)},
    }
    stage_i.write_json(manifest_path, manifest)
    paths = stage_i.initialize(root)
    stage_i.write_json(paths["reservations"], [{
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "manifest": str(manifest_path),
        "case_id": "R16",
        "case_name": "case",
        "segment": "s00",
        "nodes": 1,
        "requested_walltime": "00:10:00",
        "reserved_node_hours": 1.0 / 6.0,
        "state": "submitted",
        "prepared_utc": stage_i.utc_now(),
        "job_id": "12345",
    }])
    inspect_args = SimpleNamespace(
        manifest=str(manifest_path), required_time=2.0, allow_local_root=True
    )
    write_snapshot_time(0.5)
    write_history([(0.0, 0, 0, 0, 0, 0, 0), (1.0, 0, 0, 0, 0, 0, 1)])
    with pytest.raises(ValueError, match="prepared time/tlim target"):
        stage_i.inspect_segment(SimpleNamespace(
            manifest=str(manifest_path), required_time=2.5,
            allow_local_root=True,
        ))
    assert stage_i.inspect_segment(inspect_args) == 1
    incomplete = json.loads((manifest_dir / "segment_inspection.json").read_text())
    assert not incomplete["checks"]["terminal_snapshot_retained"]
    assert not incomplete["clean_for_continuation"]
    write_snapshot_time(2.0)
    assert stage_i.inspect_segment(inspect_args) == 1
    partial = json.loads((manifest_dir / "segment_inspection.json").read_text())
    assert partial["clean_for_continuation"]

    write_history([
        (0.0, 0, 0, 0, 0, 0, 0),
        (1.0, 0, 0, 1, 0, 0, 1),
        (2.0, 0, 0, 0, 0, 0, 2),
    ])
    assert stage_i.inspect_segment(inspect_args) == 1
    inspection_path = manifest_dir / "segment_inspection.json"
    inspection = json.loads(inspection_path.read_text())
    assert not inspection["checks"]["strict_lf_failure_counters_zero"]

    restart.write_text("<time>\nrestart_time = 1.5\n<par_end>\n")
    with pytest.raises(ValueError, match="explicit physical time"):
        stage_i.inspect_segment(inspect_args)
    write_history([(0.0, 0, 0, 0, 0, 0, 0), (2.0, 0, 0, 0, 0, 0, 2)])
    assert stage_i.inspect_segment(inspect_args) == 0
    inspection_path.unlink()
    sacct_path = tmp_path / "job.sacct"
    sacct_path.write_text(
        f"12345|cgl_mks24_{stage_i.EXECUTION_EPOCH_SLUG}_R16_s00|"
        "COMPLETED|0:0|1|60|"
        "submit-time|end-time|\n"
    )
    record_args = SimpleNamespace(
        manifest=str(manifest_path),
        allow_local_root=True,
        job_id="12345",
        result="accepted",
        notes="test",
        sacct_file=str(sacct_path),
    )
    with pytest.raises(ValueError, match="inspect-segment evidence"):
        stage_i.record(record_args)
    assert stage_i.inspect_segment(inspect_args) == 0
    sacct_path.write_text(
        "12345|wrong_name|COMPLETED|0:0|1|60|submit-time|end-time|\n"
    )
    with pytest.raises(ValueError, match="sacct job name"):
        stage_i.record(record_args)
    sacct_path.write_text(
        f"12345|cgl_mks24_{stage_i.EXECUTION_EPOCH_SLUG}_R16_s00|"
        "COMPLETED|0:0|1|60|"
        "submit-time|end-time|\n"
    )
    snapshot.write_bytes(snapshot.read_bytes().replace(b"variable=4", b"variable=5"))
    with pytest.raises(ValueError, match="inspection-retained file checksum"):
        stage_i.record(record_args)
    write_snapshot_time(2.0)
    assert stage_i.inspect_segment(inspect_args) == 0
    reservations = json.loads(paths["reservations"].read_text())
    unrelated = {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "manifest": str(
            paths["runs"] / "R04" / "s00"
            / "manifest" / "prepared_run.json"
        ),
        "case_id": "R04",
        "case_name": "unrelated",
        "segment": "s00",
        "nodes": 2,
        "requested_walltime": "00:10:00",
        "reserved_node_hours": 1.0 / 3.0,
        "state": "submitted",
        "prepared_utc": stage_i.utc_now(),
        "job_id": "67890",
    }
    reservations.append(unrelated)
    stage_i.write_json(paths["reservations"], reservations)
    assert stage_i.record(record_args) == 0
    reservations = json.loads(paths["reservations"].read_text())
    assert reservations[1] == unrelated
    assert stage_i.active_reservations(reservations) == [unrelated]
    stage_i.write_json(paths["reservations"], [reservations[0]])
    accounted = json.loads(manifest_path.read_text())
    assert accounted["state"] == "recorded"
    assert accounted["scientific_inspection"]["accepted"]
    parent = stage_i.verify_continuation_restart(
        output_dir / "rst" / "case.00000.rst"
    )
    assert parent["result"] == "accepted"
    assert parent["execution_epoch"] == stage_i.EXECUTION_EPOCH
    accounted["execution_epoch"] = "E01-pre-modal-driver"
    stage_i.write_json(manifest_path, accounted)
    with pytest.raises(ValueError, match="execution epoch"):
        stage_i.verify_continuation_restart(
            output_dir / "rst" / "case.00000.rst"
        )
    accounted["execution_epoch"] = stage_i.EXECUTION_EPOCH
    stage_i.write_json(manifest_path, accounted)
    continuation = tmp_path / "continuation.mhd.hst"
    continuation.write_text(
        "# [0]=time [1]=lf_dfloor [2]=lf_pfloor [3]=lf_nonfin "
        "[4]=lf_nonpos [5]=lf_hardbd [6]=lf_hwproj\n"
        "2.0 0 0 0 0 0 2\n"
        "3.0 0 0 0 0 0 3\n"
    )
    merged = tmp_path / "merged.mhd.hst"
    stage_i.merge_history_files([mhd_history, continuation], merged)
    merged_history = stage_i.parse_history(merged)
    assert merged_history["time"] == [0.0, 2.0, 3.0]
    before_reconcile = paths["reservations"].read_text()
    report = stage_i.reconcile_report(root)
    assert report["consistent"], report["issues"]
    assert paths["reservations"].read_text() == before_reconcile
    detached_manifest = json.loads(manifest_path.read_text())
    detached_manifest["allocation"]["nodes"] = 8
    detached_manifest["allocation"]["reserved_node_hours"] = stage_i.node_hours(
        8, stage_i.parse_walltime(
            detached_manifest["allocation"]["requested_walltime"]
        )
    )
    detached_reservations = json.loads(before_reconcile)
    detached_reservations[0]["nodes"] = 8
    detached_reservations[0]["reserved_node_hours"] = (
        detached_manifest["allocation"]["reserved_node_hours"]
    )
    if "execution_intent_sha256" in detached_reservations[0]:
        detached_reservations[0]["execution_intent_sha256"] = (
            stage_i.execution_intent_sha256(detached_manifest)
        )
    stage_i.write_json(manifest_path, detached_manifest)
    stage_i.write_json(paths["reservations"], detached_reservations)
    detached = stage_i.reconcile_report(root)
    assert not detached["consistent"]
    assert any(
        "recorded ledger provenance differs" in issue
        and "nodes differs from manifest" in issue
        for issue in detached["issues"]
    )
    stage_i.write_json(manifest_path, accounted)
    stage_i.write_json(paths["reservations"], json.loads(before_reconcile))
    reservations = json.loads(before_reconcile)
    reservations[0]["segment"] = "wrong"
    stage_i.write_json(paths["reservations"], reservations)
    inconsistent = stage_i.reconcile_report(root)
    assert not inconsistent["consistent"]
    assert any(
        "reservation record is invalid" in issue
        and "identity differs from manifest path" in issue
        for issue in inconsistent["issues"]
    )


def test_cgl_lf_stage_i_hardens_identifiers_overrides_json_and_locking(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_hardening_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    assert stage_i.require_safe_segment("s06_rankio_t6p438280_t6p5")
    assert stage_i.require_numeric_job_id("4745305") == "4745305"
    assert stage_i.expected_job_name({
        "run": {"case_id": "R16", "segment": "s-01"},
    }) == f"cgl_mks24_{stage_i.EXECUTION_EPOCH_SLUG}_R16_s_01"
    for value in ("../escape", "bad/name", "analysis", "x" * 30):
        with pytest.raises(ValueError, match="--segment"):
            stage_i.require_safe_segment(value)
    for value in ("0", "-1", "123.batch", "not-a-job"):
        with pytest.raises(ValueError, match="--job-id"):
            stage_i.require_numeric_job_id(value)

    input_path = tmp_path / "case.athinput"
    input_path.write_text("<time>\ntlim = 10.0\n<job>\nbasename = case\n")
    assert stage_i.validate_prepare_overrides(
        input_path, ["time/tlim=2.5"]
    ) == 2.5
    assert stage_i.validate_prepare_overrides(
        input_path, [], allow_missing_time_target=True
    ) is None
    with pytest.raises(ValueError, match="absent from input deck"):
        stage_i.validate_prepare_overrides(
            input_path, ["time/tlim=2.5", "mhd/missing=true"]
        )
    with pytest.raises(ValueError, match="exactly one"):
        stage_i.validate_prepare_overrides(input_path, [])
    with pytest.raises(ValueError, match="numeric"):
        stage_i.validate_prepare_overrides(input_path, ["time/tlim=not-a-time"])
    for target in ("0", "-1"):
        with pytest.raises(ValueError, match="positive"):
            stage_i.validate_prepare_overrides(
                input_path, [f"time/tlim={target}"]
            )

    metadata = tmp_path / "metadata.json"
    metadata.write_text('{"old": true}\n')
    original_replace = stage_i.os.replace
    replacements = []

    def capture_replace(source, destination):
        replacements.append((Path(source), Path(destination)))
        original_replace(source, destination)

    monkeypatch.setattr(stage_i.os, "replace", capture_replace)
    stage_i.write_json(metadata, {"new": True})
    assert json.loads(metadata.read_text()) == {"new": True}
    assert len(replacements) == 1
    assert replacements[0][0] != metadata
    assert replacements[0][1] == metadata
    assert not list(tmp_path.glob(".metadata.json.*.tmp"))

    local_root = tmp_path / "offline"
    with stage_i.canonical_root_lock(local_root):
        pass
    assert not local_root.exists()

    paths = stage_i.initialize(tmp_path / "stores")
    original_ledger = paths["ledger"].read_text()
    paths["ledger"].write_text("")
    with pytest.raises(ValueError, match="ledger is empty"):
        stage_i.read_ledger(paths)
    paths["ledger"].write_text("wrong,header\n")
    with pytest.raises(ValueError, match="ledger header is invalid"):
        stage_i.read_ledger(paths)
    invalid_ledger_row = {column: "" for column in stage_i.LEDGER_COLUMNS}
    invalid_ledger_row.update({
        "nodes": "1",
        "requested_walltime": "01:00:00",
        "elapsed_seconds": "1",
        "reserved_node_hours": "1.0",
        "actual_node_hours": "nan",
        "cumulative_stage_i_node_hours": "1.0",
    })
    paths["ledger"].write_text(
        ",".join(stage_i.LEDGER_COLUMNS)
        + "\n"
        + ",".join(
            invalid_ledger_row[column] for column in stage_i.LEDGER_COLUMNS
        )
        + "\n"
    )
    with pytest.raises(ValueError, match="invalid numeric fields"):
        stage_i.read_ledger(paths)
    paths["ledger"].write_text(original_ledger)
    with pytest.raises(ValueError, match="invalid numeric fields"):
        stage_i.append_ledger_row(paths["ledger"], invalid_ledger_row, [])
    assert paths["ledger"].read_text() == original_ledger
    ceiling_ledger_row = dict(invalid_ledger_row)
    ceiling_ledger_row["elapsed_seconds"] = "3600"
    ceiling_ledger_row["actual_node_hours"] = "1.000000"
    ceiling_ledger_row["cumulative_stage_i_node_hours"] = str(
        stage_i.CURRENT_STAGE_I_RESERVED_NODE_HOURS + 1.0
    )
    with pytest.raises(ValueError, match="accounting ceiling"):
        stage_i.validate_ledger_cumulative_fields(
            ceiling_ledger_row,
            stage_i.CURRENT_STAGE_I_RESERVED_NODE_HOURS,
            "prospective ledger row",
        )
    scheduler_fixture = tmp_path / "fixture.sacct"
    scheduler_fixture.write_text(
        "999|fixture|COMPLETED|0:0|1|1|submit|end\n"
    )
    scheduler_fixture_args = SimpleNamespace(
        sacct_file=str(scheduler_fixture), job_id="999"
    )
    scheduler_archive = paths["accounting"] / "999.stage_i.sacct.txt"
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", paths["root"])
        with pytest.raises(ValueError, match="only for offline local roots"):
            stage_i.sacct_output(
                scheduler_fixture_args, paths,
                allow_fixture=stage_i.is_offline_local_root(
                    paths["root"], allow_local_root=True
                ),
            )
    assert not scheduler_archive.exists()
    assert stage_i.sacct_output(
        scheduler_fixture_args, paths, allow_fixture=True
    ) == scheduler_fixture.read_text()
    assert scheduler_archive.read_text() == scheduler_fixture.read_text()
    scheduler_archive.unlink()

    analysis = paths["runs"] / "R02" / "analysis"
    analysis.mkdir(parents=True)
    (analysis / "recost.json").write_text("{}\n")
    assert stage_i.orphaned_segment_run_directories(paths) == []
    orphan = paths["runs"] / "R02" / "s_orphan"
    orphan.mkdir()
    assert stage_i.orphaned_segment_run_directories(paths) == [orphan]
    with pytest.raises(ValueError, match="interrupted-prepare cleanup"):
        stage_i.require_no_orphaned_segment_runs(paths)
    orphan.rmdir()

    resources = {
        "allocation": {
            "nodes": 1,
            "requested_walltime": "00:20:00",
            "requested_seconds": 1200,
            "reserved_node_hours": 1.0 / 3.0,
            "ranks_per_node": 8,
            "cpus_per_task": 7,
        },
        "command": {"athena_walltime": "00:10:00"},
        "run": {"case_id": "R02"},
    }
    stage_i.validate_prepared_resources(resources, canonical_production=True)
    resources["allocation"]["nodes"] = 8
    resources["allocation"]["reserved_node_hours"] = 8.0 / 3.0
    with pytest.raises(ValueError, match="R02 canonical Stage I preparation"):
        stage_i.validate_prepared_resources(resources, canonical_production=True)
    stage_i.validate_prepared_resources(resources, canonical_production=False)
    resources["allocation"]["nodes"] = 1
    resources["allocation"]["reserved_node_hours"] = 1.0 / 3.0
    resources["allocation"]["reserved_node_hours"] = float("nan")
    with pytest.raises(ValueError, match="resource values must be positive"):
        stage_i.validate_prepared_resources(resources, canonical_production=True)
    resources["allocation"]["reserved_node_hours"] = 1.0 / 3.0
    resources["allocation"]["ranks_per_node"] = 0
    with pytest.raises(ValueError, match="resource shape"):
        stage_i.validate_prepared_resources(resources, canonical_production=True)
    resources["allocation"]["ranks_per_node"] = 8
    resources["command"]["athena_walltime"] = "00:20:00"
    with pytest.raises(ValueError, match="not shorter"):
        stage_i.validate_prepared_resources(resources, canonical_production=True)

    restart_a = tmp_path / "a.rst"
    restart_b = tmp_path / "b.rst"
    restart_a.write_text("<time>\nrestart_time = 1.0\n<par_end>\n")
    restart_b.write_text("<time>\nrestart_time = 1.5\n<par_end>\n")
    with pytest.raises(ValueError, match="markers disagree"):
        stage_i.restart_product_time([restart_a, restart_b])

    transaction = paths["transactions"] / "outside.json"
    stage_i.write_json(transaction, {
        "schema_version": 1,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "transaction_id": transaction.stem,
        "kind": "submit_pending",
        "created_utc": stage_i.utc_now(),
        "manifest_path": str(tmp_path / "outside.json"),
        "prior_reservations": [],
        "prior_reservations_sha256": stage_i.stable_json_sha256([]),
        "prepared_manifest_sha256": "a" * 64,
        "submission_audit": {
            "created_utc": stage_i.utc_now(),
            "offline_local_root": True,
            "skip_slurm_test": True,
            "slurm_test_only": "offline local-root fixture",
            "acknowledged_shared_root_campaigns": [],
        },
    })
    with pytest.raises(ValueError, match="outside the E03 run store"):
        stage_i.read_transaction(paths, transaction)

    canonical_root = tmp_path / "canonical"
    canonical_root.mkdir()
    monkeypatch.setattr(stage_i, "DEFAULT_ROOT", canonical_root)
    lock_path = stage_i.canonical_root_lock_path(canonical_root)
    prior_umask = os.umask(0o077)
    try:
        with stage_i.canonical_root_lock(canonical_root):
            assert lock_path.is_file()
            with stage_i.canonical_root_lock(canonical_root):
                pass
    finally:
        os.umask(prior_umask)
    lock_profile = lock_path.stat()
    assert stat.S_IMODE(lock_profile.st_mode) == 0o644
    assert lock_profile.st_nlink == 1
    assert lock_profile.st_uid == os.geteuid()
    with lock_path.open("a+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(ValueError, match="another Stage I mutation"):
            with stage_i.canonical_root_lock(canonical_root):
                pass
    lock_path.chmod(0o600)
    with stage_i.canonical_root_lock(canonical_root):
        pass
    assert stat.S_IMODE(lock_path.stat().st_mode) == 0o644

    lock_path.unlink()
    target = canonical_root / "lock-target"
    target.write_text("")
    lock_path.symlink_to(target.name)
    with pytest.raises(ValueError, match="must not be a symlink"):
        with stage_i.canonical_root_lock(canonical_root):
            pass
    lock_path.unlink()
    os.link(target, lock_path)
    with pytest.raises(ValueError, match="link count differs"):
        with stage_i.canonical_root_lock(canonical_root):
            pass
    lock_path.unlink()
    target.unlink()

    lock_path.write_text("")
    original_flock = stage_i.fcntl.flock

    def replace_lock_path(descriptor, operation):
        original_flock(descriptor, operation)
        if operation & fcntl.LOCK_EX:
            lock_path.unlink()
            lock_path.write_text("")

    with monkeypatch.context() as policy:
        policy.setattr(stage_i.fcntl, "flock", replace_lock_path)
        with pytest.raises(ValueError, match="path changed while locking"):
            with stage_i.canonical_root_lock(canonical_root):
                pass
    lock_path.chmod(0o666)
    with pytest.raises(ValueError, match="mode is too permissive"):
        with stage_i.canonical_root_lock(canonical_root):
            pass

    queue_fixture = tmp_path / "fixture.squeue"
    queue_fixture.write_text("")
    fixture_args = SimpleNamespace(
        squeue_file=str(queue_fixture), skip_slurm_test=False
    )
    assert stage_i.production_queue_output(fixture_args, True) == ""
    with pytest.raises(ValueError, match="restricted to offline local roots"):
        stage_i.validate_submission_fixture_options(fixture_args, False)
    with pytest.raises(ValueError, match="requires --squeue-file"):
        stage_i.production_queue_output(SimpleNamespace(squeue_file=None), True)
    with pytest.raises(ValueError, match="restricted to offline local roots"):
        stage_i.validate_submission_fixture_options(
            SimpleNamespace(squeue_file=None, skip_slurm_test=True), False
        )
    canonical_manifest = {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "project_root": str(canonical_root),
    }
    with pytest.raises(ValueError, match="--squeue-file"):
        stage_i.submission_preflight(
            SimpleNamespace(
                allow_local_root=False,
                squeue_file=str(queue_fixture),
                skip_slurm_test=False,
            ),
            tmp_path / "canonical.json",
            canonical_manifest,
            run_slurm_test=True,
        )
    with pytest.raises(ValueError, match="--skip-slurm-test"):
        stage_i.submission_preflight(
            SimpleNamespace(
                allow_local_root=False,
                squeue_file=None,
                skip_slurm_test=True,
            ),
            tmp_path / "canonical.json",
            canonical_manifest,
            run_slurm_test=True,
        )

    captured = {}

    def capture_squeue(command, **kwargs):
        captured["command"] = command
        captured["environment"] = kwargs["env"]
        return SimpleNamespace(stdout="")

    with monkeypatch.context() as policy:
        policy.setenv("USER", "forged-user")
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_squeue)
        assert stage_i.production_queue_output(
            SimpleNamespace(squeue_file=None), False
        ) == ""
    assert captured["command"] == [
        "/usr/bin/squeue", "-h", "-u", stage_i.pwd.getpwuid(os.geteuid()).pw_name,
        "-o", "%i|%P|%T|%j",
    ]
    assert "SLURM_CONF" not in captured["environment"]
    expected_scheduler_user = stage_i.pwd.getpwuid(os.geteuid()).pw_name
    assert captured["command"][3] == expected_scheduler_user
    assert captured["command"][3] != "forged-user"
    assert stage_i.SACCT == Path("/usr/bin/sacct")
    assert stage_i.SBATCH == Path("/usr/bin/sbatch")
    assert stage_i.SCONTROL == Path("/usr/bin/scontrol")
    assert stage_i.SQUEUE == Path("/usr/bin/squeue")

    scheduler_calls = []
    submitted_script = tmp_path / "submitted.sbatch"
    submitted_script.write_text("#!/bin/bash\n")
    submitted_manifest_path = tmp_path / "submitted.json"
    submitted_manifest = {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "project_root": str(tmp_path / "production-root"),
        "state": "prepared",
        "run": {"case_id": "R16", "segment": "s00"},
        "paths": {"batch_script": str(submitted_script)},
    }

    def capture_submit(command, **kwargs):
        scheduler_calls.append((command, kwargs["env"]))
        return SimpleNamespace(stdout="12345\n")

    with monkeypatch.context() as policy:
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_submit)
        assert stage_i.scheduler_test_only_output(submitted_script) == "12345\n"
    assert scheduler_calls == [(
        ["/usr/bin/sbatch", "--test-only", str(submitted_script)],
        scheduler_calls[0][1],
    )]
    assert "SLURM_CONF" not in scheduler_calls[0][1]

    scheduler_calls.clear()

    with monkeypatch.context() as policy:
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_submit)
        policy.setattr(
            stage_i, "read_manifest", lambda _path: submitted_manifest
        )
        policy.setattr(stage_i, "require_root", lambda root, _allowed: root)
        policy.setattr(stage_i, "is_offline_local_root", lambda *_args: False)
        policy.setattr(
            stage_i,
            "submission_preflight",
            lambda *_args, **_kwargs: (paths, submitted_script, {}),
        )
        policy.setattr(
            stage_i,
            "authenticated_production_queue_evidence",
            lambda *_args, **_kwargs: {
                "checked_utc": stage_i.utc_now(),
                "rows": [],
                "rows_sha256": stage_i.stable_json_sha256([]),
            },
        )
        policy.setattr(
            stage_i,
            "write_submit_pending_transaction",
            lambda *_args: tmp_path / "submit-pending.json",
        )
        policy.setattr(stage_i, "read_reservations", lambda _paths: [])
        policy.setattr(stage_i, "finish_submit_transaction", lambda *_args: None)
        assert stage_i.submit.__wrapped__(SimpleNamespace(
            manifest=str(submitted_manifest_path),
            allow_local_root=False,
            sbatch_output_file=None,
        )) == 0
    assert scheduler_calls == [(
        ["/usr/bin/sbatch", "--parsable", str(submitted_script)],
        scheduler_calls[0][1],
    )]
    assert "SLURM_CONF" not in scheduler_calls[0][1]

    scheduler_calls.clear()
    transaction = {"created_utc": stage_i.utc_now()}
    expected_name = stage_i.expected_job_name(submitted_manifest)

    def capture_recovery(command, **kwargs):
        scheduler_calls.append((command, kwargs["env"]))
        assert command[0] == "/usr/bin/scontrol"
        return SimpleNamespace(
            returncode=0,
            stdout=(
                f"JobId=12345 JobName={expected_name} Account={stage_i.ACCOUNT} "
                f"Partition={stage_i.PARTITION} SubmitTime={transaction['created_utc']} "
                f"Command={submitted_script}\n"
            ),
        )

    with monkeypatch.context() as policy:
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_recovery)
        recovered = stage_i.verify_recovered_scheduler_job(
            submitted_manifest, "12345", False, transaction
        )
    assert recovered["mode"] == "scontrol"
    assert scheduler_calls[0][0] == [
        "/usr/bin/scontrol", "show", "job", "-o", "12345",
    ]
    assert "SLURM_CONF" not in scheduler_calls[0][1]

    scheduler_calls.clear()

    def capture_recovery_fallback(command, **kwargs):
        scheduler_calls.append((command, kwargs["env"]))
        if command[0] == "/usr/bin/scontrol":
            return SimpleNamespace(returncode=1, stdout="")
        assert command[0] == "/usr/bin/sacct"
        return SimpleNamespace(
            stdout=(
                f"12345|{expected_name}|{stage_i.ACCOUNT}|{stage_i.PARTITION}|"
                f"{transaction['created_utc']}\n"
            ),
        )

    with monkeypatch.context() as policy:
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_recovery_fallback)
        recovered = stage_i.verify_recovered_scheduler_job(
            submitted_manifest, "12345", False, transaction
        )
    assert recovered["mode"] == "sacct"
    assert [call[0][0] for call in scheduler_calls] == [
        "/usr/bin/scontrol", "/usr/bin/sacct",
    ]
    assert all("SLURM_CONF" not in environment for _, environment in scheduler_calls)

    scheduler_calls.clear()

    def capture_absence(command, **kwargs):
        scheduler_calls.append((command, kwargs["env"]))
        return SimpleNamespace(stdout="")

    with monkeypatch.context() as policy:
        policy.setenv("USER", "forged-user")
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_absence)
        absence = stage_i.scheduler_absence_evidence(
            SimpleNamespace(scheduler_absence_evidence_file=None),
            submitted_manifest,
            transaction,
            False,
        )
    assert absence["mode"] == "live scheduler absence query"
    assert [call[0][0] for call in scheduler_calls] == [
        "/usr/bin/squeue", "/usr/bin/sacct",
    ]
    assert scheduler_calls[0][0][3] == expected_scheduler_user
    assert scheduler_calls[0][0][3] != "forged-user"
    assert all("SLURM_CONF" not in environment for _, environment in scheduler_calls)

    scheduler_calls.clear()

    def capture_sacct(command, **kwargs):
        scheduler_calls.append((command, kwargs["env"]))
        return SimpleNamespace(stdout="")

    with monkeypatch.context() as policy:
        policy.setenv("SLURM_CONF", "/tmp/forged-slurm.conf")
        policy.setattr(stage_i.subprocess, "run", capture_sacct)
        assert stage_i.sacct_output(
            SimpleNamespace(sacct_file=None, job_id="12345"),
            paths,
            allow_fixture=False,
        ) == ""
    assert scheduler_calls[0][0][0] == "/usr/bin/sacct"
    assert "SLURM_CONF" not in scheduler_calls[0][1]


def test_cgl_lf_stage_i_groups_rank_local_output_products(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_rank_output_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    restart_dir = tmp_path / "rst"
    for rank in range(2):
        path = restart_dir / f"rank_{rank:08d}" / "case.00000.rst"
        path.parent.mkdir(parents=True)
        path.write_bytes(f"rank-{rank}".encode())
    groups = stage_i.output_product_groups(restart_dir, "*.rst", expected_ranks=2)
    assert len(groups) == 1
    product = stage_i.retained_product(groups[0])
    assert product["storage"] == "per_rank"
    assert len(product["rank_files"]) == 2
    assert stage_i.retained_product_paths(product) == groups[0]


def test_cgl_lf_stage_i_requires_retained_source_bundle_provenance(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_bundle_provenance_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    repository = tmp_path / "source"
    repository.mkdir()
    subprocess.run(["git", "init", "--quiet"], cwd=repository, check=True)
    subprocess.run(
        ["git", "config", "user.email", "test@example.invalid"],
        cwd=repository, check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "CGL-LF test"],
        cwd=repository, check=True,
    )
    tracked = repository / "tracked.txt"
    tracked.write_text("retained source provenance\n")
    subprocess.run(["git", "add", "tracked.txt"], cwd=repository, check=True)
    subprocess.run(
        ["git", "commit", "--quiet", "-m", "test"], cwd=repository, check=True
    )
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository,
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    root = tmp_path / "root"
    bundle = root / "source-archives" / "test.bundle"
    bundle.parent.mkdir(parents=True)
    subprocess.run(
        ["git", "bundle", "create", str(bundle), "--all"],
        cwd=repository, check=True,
    )

    provenance = stage_i.source_bundle_provenance(
        str(bundle), [revision], root, allow_local_root=False
    )
    assert provenance is not None
    assert provenance["path"] == str(bundle)
    assert provenance["sha256"] == hashlib.sha256(bundle.read_bytes()).hexdigest()
    assert provenance["verified_revisions"] == [revision]
    with pytest.raises(ValueError, match="--source-bundle is required"):
        stage_i.source_bundle_provenance(
            None, [revision], root, allow_local_root=False
        )
    with pytest.raises(ValueError, match="does not contain revision"):
        stage_i.source_bundle_provenance(
            str(bundle), ["0" * 40], root, allow_local_root=False
        )


def test_cgl_lf_stage_i_authenticates_legacy_restart_marker_with_binary_time(
    tmp_path,
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_restart_binary_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    executable_revision, executable_sha256 = next(
        iter(stage_i.QUALIFIED_RESTART_BINARY_ABIS)
    )
    command = {
        "executable_revision": executable_revision,
        "executable_sha256": executable_sha256,
    }
    binary_time = 0.31282347945569927
    offset = stage_i.QUALIFIED_RESTART_BINARY_ABIS[
        executable_revision, executable_sha256
    ]["mesh_time_offset_after_parameter_dump"]
    assert offset == 232

    def write_restart(path, marker, value, prefix=b""):
        path.write_bytes(
            prefix
            + f"<time>\nrestart_time = {marker}\n<par_end>\n".encode()
            + b"\0" * offset
            + struct.pack("<d", value)
        )

    restart_a = tmp_path / "rank_00000000.rst"
    restart_b = tmp_path / "rank_00000001.rst"
    write_restart(restart_a, "0.312823", binary_time)
    write_restart(restart_b, "0.312823", binary_time)
    authenticated = stage_i.authenticated_restart_product_time(
        [restart_a, restart_b], command
    )
    assert authenticated["binary_time"] == binary_time
    assert authenticated["marker_modes"] == [
        "legacy_default_precision", "legacy_default_precision",
    ]

    write_restart(restart_b, "0.31282347945569927", binary_time)
    with pytest.raises(ValueError, match="marker modes disagree"):
        stage_i.authenticated_restart_product_time([restart_a, restart_b], command)
    write_restart(restart_a, "0.31282347945569927", binary_time)
    assert stage_i.authenticated_restart_product_time(
        [restart_a, restart_b], command
    )["marker_modes"] == ["full_precision", "full_precision"]
    write_restart(restart_b, "0.31282347945569927", binary_time + 5.0e-13)
    with pytest.raises(ValueError, match="binary physical times disagree"):
        stage_i.authenticated_restart_product_time([restart_a, restart_b], command)
    write_restart(restart_b, "0.312824", binary_time)
    with pytest.raises(ValueError, match="does not authenticate"):
        stage_i.authenticated_restart_product_time([restart_a, restart_b], command)
    with pytest.raises(ValueError, match="no qualified binary-restart ABI"):
        stage_i.authenticated_restart_product_time(
            [restart_a], {**command, "executable_sha256": "0" * 64}
        )
    prefix = b"<time>\nrestart_time = 1\n"
    padding = b"#" * (
        stage_i.MAX_RESTART_PARAMETER_DUMP_BYTES
        - len(prefix)
        - len(b"<par_end>\n")
    )
    boundary = tmp_path / "boundary.rst"
    boundary.write_bytes(
        prefix + padding + b"<par_end>\n" + b"\0" * offset + struct.pack("<d", 1.0)
    )
    assert stage_i.restart_binary_time(boundary, command) == 1.0
    oversized = tmp_path / "oversized.rst"
    oversized.write_bytes(
        prefix + padding + b"#<par_end>\n"
        + b"\0" * offset + struct.pack("<d", 1.0)
    )
    with pytest.raises(ValueError, match="lacks loadable"):
        stage_i.restart_binary_time(oversized, command)


def test_cgl_lf_stage_i_authorizes_only_reviewed_historical_r03_submission():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_historical_r03_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    expected = stage_i.HISTORICAL_SUBMITTED_R03_UTILITY_TRANSITION
    manifest = {
        "project_root": expected["project_root"],
        "state": expected["state"],
        "job_id": expected["job_id"],
        "run": {
            "case_id": expected["case_id"],
            "segment": expected["segment"],
        },
        "command": {
            "production_utility": {
                "revision": expected["production_utility_revision"],
                "sha256": expected["production_utility_sha256"],
            },
            "source_bundle": {"sha256": expected["source_bundle_sha256"]},
            "executable_revision": expected["executable_revision"],
            "executable_sha256": expected["executable_sha256"],
        },
    }
    assert stage_i.historical_submitted_utility_transition_authorized(manifest)
    manifest["job_id"] = "4762473"
    assert not stage_i.historical_submitted_utility_transition_authorized(manifest)


def test_cgl_lf_stage_i_rejects_prospective_budget_overruns_and_counts_rows_once(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_budget_replay_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    reservation = {"state": "prepared", "reserved_node_hours": 1.0 / 3.0}
    with pytest.raises(ValueError, match="Stage I reservation ceiling"):
        stage_i.require_reservation_budget(
            [{"actual_node_hours": "1399.75"}],
            [reservation],
            "prospective fixture",
        )
    monkeypatch.setattr(stage_i, "CURRENT_STAGE_I_RESERVED_NODE_HOURS", 5000.0)
    with pytest.raises(ValueError, match="incremental project ceiling"):
        stage_i.require_reservation_budget(
            [{"actual_node_hours": "3999.75"}],
            [reservation],
            "prospective fixture",
        )

    paths = stage_i.initialize(tmp_path / "root")
    row = {"job_id": "123", "actual_node_hours": "0.25"}
    monkeypatch.setattr(stage_i, "read_ledger", lambda _paths: [row])
    baseline, prospective = stage_i.transaction_ledger_views(paths, row)
    assert baseline == []
    assert prospective == [row]


def test_cgl_lf_stage_i_replay_budget_gate_precedes_controlled_writes(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_replay_apply_budget_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    paths = stage_i.initialize(tmp_path / "root")
    manifest_path = (
        paths["runs"] / "R03" / "s_replay"
        / "manifest" / "prepared_run.json"
    )
    journal = paths["transactions"] / "replay.json"
    journal.write_text("{}\n")
    row = {"job_id": "123", "actual_node_hours": "0.25"}
    transaction = {
        "kind": "recorded",
        "manifest_path": str(manifest_path),
        "manifest": {},
        "reservations": [{
            "state": "submitted",
            "reserved_node_hours": 0.1,
        }],
        "ledger_row": row,
    }
    monkeypatch.setattr(stage_i, "read_transaction", lambda *_args: transaction)
    monkeypatch.setattr(
        stage_i, "validate_transaction_reservation_baseline", lambda *_args: None
    )
    monkeypatch.setattr(
        stage_i,
        "read_ledger",
        lambda _paths: [{"job_id": "seed", "actual_node_hours": "1399.75"}],
    )
    writes = []
    appends = []
    monkeypatch.setattr(stage_i, "write_json", lambda *args: writes.append(args))
    monkeypatch.setattr(
        stage_i, "append_ledger_row", lambda *args: appends.append(args)
    )
    with pytest.raises(ValueError, match="Stage I reservation ceiling"):
        stage_i.apply_transaction(paths, journal)
    assert writes == []
    assert appends == []
    assert journal.is_file()

    transaction["reservations"] = []
    monkeypatch.setattr(stage_i, "read_ledger", lambda _paths: [row])
    unlinked = []
    monkeypatch.setattr(stage_i, "refresh_summary", lambda *_args: None)
    monkeypatch.setattr(stage_i, "unlink_durable", lambda path: unlinked.append(path))
    stage_i.apply_transaction(paths, journal)
    assert appends == []
    assert writes == [
        (paths["reservations"], []),
        (manifest_path, {}),
    ]
    assert unlinked == [journal]


def test_cgl_lf_stage_i_retains_r17_last_across_lifecycle(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_retained_r17_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    paths = stage_i.initialize(tmp_path / "root")
    monkeypatch.setattr(stage_i, "DEFAULT_ROOT", paths["root"])
    monkeypatch.setattr(stage_i, "accepted_case_lineage", lambda *_args: [])
    with pytest.raises(ValueError, match="R17 must remain last"):
        stage_i.require_retained_r17_policy(
            paths, [{"case_id": "R17", "state": "submitted"}]
        )

    r17 = paths["runs"] / "R17" / "s00" / "manifest" / "prepared_run.json"
    r17.parent.mkdir(parents=True)
    r17.write_text("{}\n")
    monkeypatch.setattr(
        stage_i,
        "accepted_case_lineage",
        lambda _paths, _case_id: [{"scientific_inspection": {"final_time": 10.0}}],
    )
    with pytest.raises(ValueError, match="active lower-resolution lanes"):
        stage_i.require_retained_r17_policy(
            paths, [{"case_id": "R03", "state": "submitted"}]
        )


def test_cgl_lf_stage_i_rejects_late_r17_replay_and_reconcile(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_r17_replay_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    root = tmp_path / "root"
    paths = stage_i.initialize(root)
    monkeypatch.setattr(stage_i, "DEFAULT_ROOT", root)
    monkeypatch.setattr(stage_i, "accepted_case_lineage", lambda *_args: [])
    manifest_path = (
        paths["runs"] / "R17" / "s00"
        / "manifest" / "prepared_run.json"
    )

    def reservation(state):
        record = {
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "manifest": str(manifest_path),
            "case_id": "R17",
            "case_name": "case",
            "segment": "s00",
            "nodes": 8,
            "requested_walltime": "00:10:00",
            "reserved_node_hours": 8.0 / 6.0,
            "state": state,
            "prepared_utc": stage_i.utc_now(),
            "execution_intent_sha256": "a" * 64,
        }
        if state in {"submitted", "recorded"}:
            record["job_id"] = "321"
        if state == "recorded":
            record["actual_node_hours"] = 0.25
            record["result"] = "clean_partial"
        return record

    def write_transaction(kind, prior, payload):
        journal = paths["transactions"] / f"r17-{kind}.json"
        value = {
            "schema_version": 1,
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "transaction_id": journal.stem,
            "kind": kind,
            "created_utc": stage_i.utc_now(),
            "manifest_path": str(manifest_path),
            "prior_reservations": [prior],
            "prior_reservations_sha256": stage_i.stable_json_sha256([prior]),
            "manifest": {
                "execution_epoch": stage_i.EXECUTION_EPOCH,
                "project_root": str(root),
                "state": kind,
                "run": {
                    "case_id": "R17",
                    "case_name": "case",
                    "segment": "s00",
                },
            },
            "reservations": [payload],
            "ledger_row": None,
        }
        if kind == "submitted":
            value.update({
                "prepared_manifest_sha256": "b" * 64,
                "submission_audit": {},
                "job_id": "321",
                "submitted_recorded_utc": stage_i.utc_now(),
            })
        stage_i.write_json(journal, value)
        return journal

    submitted = reservation("submitted")
    submitted_journal = write_transaction(
        "submitted", reservation("prepared"), submitted
    )
    with pytest.raises(ValueError, match="R17 must remain last"):
        stage_i.read_transaction(paths, submitted_journal)
    submitted_journal.unlink()

    recorded_journal = write_transaction(
        "recorded", submitted, reservation("recorded")
    )
    with pytest.raises(ValueError, match="R17 must remain last"):
        stage_i.read_transaction(paths, recorded_journal)
    recorded_journal.unlink()

    stage_i.write_json(paths["reservations"], [submitted])
    report = stage_i.reconcile_report(root)
    assert any(
        issue.startswith("retained R17 policy is invalid: R17 must remain last")
        for issue in report["issues"]
    )


def test_cgl_lf_stage_i_authenticates_historical_production_utility(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_historical_utility_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    repository = tmp_path / "source"
    script = repository / "scripts" / "frontier" / "stage_i.py"
    script.parent.mkdir(parents=True)
    subprocess.run(["git", "init", "--quiet"], cwd=repository, check=True)
    subprocess.run(
        ["git", "config", "user.email", "test@example.invalid"],
        cwd=repository, check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "CGL-LF test"],
        cwd=repository, check=True,
    )
    script.write_text("historical helper\n")
    subprocess.run(["git", "add", "."], cwd=repository, check=True)
    subprocess.run(
        ["git", "commit", "--quiet", "-m", "historical"], cwd=repository, check=True
    )
    historical_revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository,
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    historical_sha = hashlib.sha256(script.read_bytes()).hexdigest()
    bundle = tmp_path / "historical.bundle"
    subprocess.run(
        ["git", "bundle", "create", str(bundle), "--all"],
        cwd=repository, check=True,
    )
    bundle_record = {
        "path": str(bundle),
        "sha256": hashlib.sha256(bundle.read_bytes()).hexdigest(),
        "verified_revisions": [historical_revision],
    }

    script.write_text("current helper\n")
    subprocess.run(["git", "add", "."], cwd=repository, check=True)
    subprocess.run(
        ["git", "commit", "--quiet", "-m", "current"], cwd=repository, check=True
    )
    monkeypatch.setattr(stage_i, "ROOT_DIR", repository)
    monkeypatch.setattr(stage_i, "__file__", str(script))
    monkeypatch.setattr(
        stage_i, "PRODUCTION_UTILITY_RELATIVE", Path("scripts/frontier/stage_i.py")
    )
    record = {
        "path": str(script),
        "revision": historical_revision,
        "sha256": historical_sha,
        "committed": True,
    }
    with pytest.raises(ValueError, match="checksum has changed"):
        stage_i.authenticate_production_utility(record)
    stage_i.authenticate_production_utility(
        record, source_bundle=bundle_record, allow_historical=True
    )
    script.write_text("dirty current helper\n")
    with pytest.raises(ValueError, match="must be committed"):
        stage_i.authenticate_production_utility(
            record, source_bundle=bundle_record, allow_historical=True
        )
    with pytest.raises(ValueError, match="must be committed"):
        stage_i.authenticate_production_utility(record)
    subprocess.run(
        ["git", "checkout", "--", "scripts/frontier/stage_i.py"],
        cwd=repository, check=True,
    )
    monkeypatch.setattr(
        stage_i, "PRODUCTION_UTILITY_RELATIVE", Path("scripts/frontier/wrong.py")
    )
    with pytest.raises(ValueError, match="path is inconsistent"):
        stage_i.authenticate_production_utility(
            record, source_bundle=bundle_record, allow_historical=True
        )
    monkeypatch.setattr(
        stage_i, "PRODUCTION_UTILITY_RELATIVE", Path("scripts/frontier/stage_i.py")
    )
    with pytest.raises(ValueError, match="checksum has changed"):
        stage_i.authenticate_production_utility(
            {**record, "sha256": "0" * 64},
            source_bundle=bundle_record,
            allow_historical=True,
        )
    with pytest.raises(ValueError, match="lacks the historical"):
        stage_i.authenticate_production_utility(
            record,
            source_bundle={**bundle_record, "verified_revisions": []},
            allow_historical=True,
        )

    manifest_dir = tmp_path / "run" / "manifest"
    manifest_dir.mkdir(parents=True)
    batch_script = manifest_dir / "cgl_lf_stage_i.sbatch"
    batch_digest = "a" * 64
    batch_script.write_text(f"BATCH_SCRIPT_SHA256={batch_digest}\n")
    input_file = manifest_dir / "submitted_input.athinput"
    matrix_file = manifest_dir / "mks24_stage_i_manifest.json"
    executable = tmp_path / "athena"
    input_file.write_text("input\n")
    matrix_file.write_text("{}\n")
    executable.write_text("#!/bin/sh\n")
    executable.chmod(0o755)
    manifest_path = manifest_dir / "prepared_run.json"
    manifest = {
        "state": "recorded",
        "command": {
            "batch_script_sha256": batch_digest,
            "input_file": str(input_file),
            "input_sha256": stage_i.sha256(input_file),
            "matrix_file": str(matrix_file),
            "matrix_sha256": stage_i.sha256(matrix_file),
            "executable": str(executable),
            "executable_sha256": stage_i.sha256(executable),
            "restart_files": [],
            "production_utility": record,
            "source_bundle": bundle_record,
            "input_revision": historical_revision,
            "executable_revision": historical_revision,
        },
        "paths": {"batch_script": str(batch_script)},
    }
    monkeypatch.setattr(stage_i, "validate_prepared_continuation_target", lambda _: None)
    monkeypatch.setattr(stage_i, "validate_prepared_resources", lambda *_args, **_kw: None)
    monkeypatch.setattr(stage_i, "normalized_batch_script_sha256", lambda _: batch_digest)
    monkeypatch.setattr(
        stage_i, "generated_batch_script",
        lambda *_: pytest.fail("recorded manifest regenerated a live batch script"),
    )
    stage_i.authenticate_prepared_execution(
        manifest, manifest_path, allow_legacy_local=True
    )
    with pytest.raises(ValueError, match="checksum has changed"):
        stage_i.authenticate_prepared_execution(
            {**manifest, "state": "prepared"},
            manifest_path,
            allow_legacy_local=True,
        )
    submitted = {**manifest, "state": "submitted"}
    with pytest.raises(ValueError, match="checksum has changed"):
        stage_i.authenticate_prepared_execution(
            submitted, manifest_path, allow_legacy_local=True,
        )
    monkeypatch.setattr(
        stage_i,
        "historical_submitted_utility_transition_authorized",
        lambda candidate: candidate.get("state") == "submitted",
    )
    stage_i.authenticate_prepared_execution(
        submitted, manifest_path, allow_legacy_local=True,
    )

    bundle.write_bytes(bundle.read_bytes() + b"tampered\n")
    with pytest.raises(ValueError, match="source bundle checksum has changed"):
        stage_i.authenticate_production_utility(
            record, source_bundle=bundle_record, allow_historical=True
        )
    with pytest.raises(ValueError, match="source bundle checksum has changed"):
        stage_i.authenticate_prepared_execution(
            manifest, manifest_path, allow_legacy_local=True
        )


def test_cgl_lf_stage_i_isolates_epoch_and_checks_all_shared_root_jobs(
    tmp_path, capsys, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_epoch_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    root = tmp_path / "root"
    paths = stage_i.initialize(root)
    assert paths["runs"] == (
        root / "runs" / "mks24-stage-i" / stage_i.EXECUTION_EPOCH
    )
    assert paths["ledger"].name == (
        f"mks24_stage_i_{stage_i.EXECUTION_EPOCH_SLUG}_node_hours.csv"
    )
    assert paths["reservations"].name == (
        f"mks24_stage_i_{stage_i.EXECUTION_EPOCH_SLUG}_reservations.json"
    )
    assert stage_i.COMPLETED_R16_NODE_HOURS == 6.145556
    assert stage_i.COMPLETED_R02_STANDARD_LAYOUT_PILOT_NODE_HOURS == 0.473333
    assert stage_i.COMPLETED_R17_HIGH_RESOLUTION_PILOT_NODE_HOURS == 4.235556
    assert stage_i.MEASURED_STAGE_I_RESERVED_NODE_HOURS == 900.0
    assert stage_i.PROMOTED_STAGE_I_RESERVED_NODE_HOURS == 1400.0
    assert stage_i.CURRENT_STAGE_I_RESERVED_NODE_HOURS == 1400.0
    assert (
        "- Current E03 mapped-matrix planning envelope: `1400.000000` node-hours"
        in paths["summary"].read_text()
    )
    assert stage_i.MAX_SEGMENT_SECONDS == 2 * 60 * 60
    assert stage_i.R17_CASE_ID == "R17"
    assert stage_i.MAX_ACTIVE_STAGE_I_SEGMENTS == 4
    assert stage_i.MAX_ACTIVE_STAGE_I_NODES == 10
    assert stage_i.R17_PREDECESSOR_CASE_IDS == tuple(
        f"R{number:02d}" for number in range(2, 17)
    )
    assert stage_i.REQUIRED_CASE_FINAL_TIME == 10.0
    stage_i.require_authorized_case("R02")
    stage_i.require_authorized_case("R17")
    with pytest.raises(
        ValueError, match="authorized only for mapped matrix cases R02-R17"
    ):
        stage_i.require_authorized_case("R18")
    stage_i.require_case_node_count("R02", 1)
    stage_i.require_case_node_count("R03", 1)
    stage_i.require_case_node_count("R04", 1)
    stage_i.require_case_node_count("R04", 2)
    stage_i.require_case_node_count("R04", 4)
    stage_i.require_case_node_count("R16", 1)
    stage_i.require_case_node_count("R16", 2)
    stage_i.require_case_node_count("R17", 8)
    with pytest.raises(ValueError, match="R02 canonical Stage I preparation"):
        stage_i.require_case_node_count("R02", 8)
    with pytest.raises(ValueError, match="R03 canonical Stage I preparation"):
        stage_i.require_case_node_count("R03", 2)
    with pytest.raises(ValueError, match="R16 canonical Stage I preparation"):
        stage_i.require_case_node_count("R16", 4)
    with pytest.raises(ValueError, match="R17 canonical Stage I preparation"):
        stage_i.require_case_node_count("R17", 1)
    active = lambda case_id, state="submitted", nodes=1: {
        "case_id": case_id, "state": state, "nodes": nodes,
    }
    stage_i.require_active_reservation_policy([
        active("R03"), active("R04", nodes=4), active("R05", nodes=4),
    ], "R06", 1)
    stage_i.require_active_reservation_policy([
        active("R04", nodes=4), active("R05", nodes=4), active("R06", nodes=2),
    ])
    with pytest.raises(ValueError, match="10-node concurrency limit"):
        stage_i.require_active_reservation_policy([
            active("R04", nodes=4), active("R05", nodes=4),
            active("R06", nodes=4),
        ])
    with pytest.raises(ValueError, match="would exceed the 10-node"):
        stage_i.require_active_reservation_policy([
            active("R04", nodes=4), active("R05", nodes=4),
        ], "R06", 4)
    with pytest.raises(ValueError, match="invalid nodes"):
        stage_i.require_active_reservation_policy([
            {"case_id": "R03", "state": "submitted"},
        ])
    with pytest.raises(ValueError, match="invalid nodes"):
        stage_i.require_active_reservation_policy([], "R03", 0)
    with pytest.raises(ValueError, match="candidate nodes require"):
        stage_i.require_active_reservation_policy([], candidate_nodes=1)
    with pytest.raises(ValueError, match="duplicate case R03"):
        stage_i.require_active_reservation_policy([
            active("R03"), active("R03"),
        ])
    with pytest.raises(ValueError, match="only one Stage I segment may be prepared"):
        stage_i.require_active_reservation_policy([
            active("R03", "prepared"), active("R04", "prepared"),
        ])
    with pytest.raises(ValueError, match="4-segment concurrency limit"):
        stage_i.require_active_reservation_policy([
            active("R03"), active("R04"), active("R05"), active("R06"),
        ], "R07", 1)
    with pytest.raises(ValueError, match="R17 requires exclusive"):
        stage_i.require_active_reservation_policy([active("R03")], "R17", 8)
    with pytest.raises(ValueError, match="R17 requires exclusive"):
        stage_i.require_active_reservation_policy([
            active("R17"), active("R03"),
        ])
    with pytest.raises(ValueError, match="distinct R03-R16"):
        stage_i.require_active_reservation_policy([active("R02")], "R03", 1)
    with pytest.raises(ValueError, match="another Stage I segment is prepared"):
        stage_i.require_active_reservation_policy([
            active("R03", "prepared"),
        ], "R04", 1)
    assert not stage_i.retained_case_has_started(paths, "R17")
    stage_i.require_r17_last(paths, "R16")
    stage_i.require_prepare_case_policy(paths, "R17", 1, offline_local_root=True)
    lineages = {
        case_id: [{"scientific_inspection": {"final_time": 10.0}}]
        for case_id in stage_i.R17_PREDECESSOR_CASE_IDS
    }
    monkeypatch.setattr(
        stage_i, "accepted_case_lineage",
        lambda _paths, case_id: lineages.get(case_id, []),
    )
    stage_i.require_r17_last(paths, "R17")
    stage_i.require_prepare_case_policy(paths, "R17", 8, offline_local_root=False)
    lineages.pop("R08")
    with pytest.raises(ValueError, match="R17 must remain last.*R08"):
        stage_i.require_r17_last(paths, "R17")
    lineages["R08"] = [{"scientific_inspection": {"final_time": 9.5}}]
    with pytest.raises(ValueError, match="R17 must remain last.*R08"):
        stage_i.require_r17_last(paths, "R17")
    lineages["R08"] = [{"scientific_inspection": {"final_time": float("nan")}}]
    with pytest.raises(ValueError, match="R17 must remain last.*R08"):
        stage_i.require_r17_last(paths, "R17")
    lineages["R08"] = [{"scientific_inspection": {"final_time": 10.0}}]
    failed_r17_manifest = (
        paths["runs"] / "R17" / "s_failed"
        / "manifest" / "prepared_run.json"
    )
    failed_r17_manifest.parent.mkdir(parents=True)
    stage_i.write_json(failed_r17_manifest, {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "state": "recorded",
        "accounting": {"result": "failed"},
    })
    assert stage_i.retained_case_has_started(paths, "R17")
    with pytest.raises(ValueError, match="R17 has started.*R16"):
        stage_i.require_r17_last(paths, "R16")
    failed_r17_manifest.unlink()
    failed_r17_manifest.parent.rmdir()
    failed_r17_manifest.parent.parent.rmdir()
    assert stage_i.parser().parse_args([
        "bundle-case", "--case-id", "R02",
    ]).required_final_time == stage_i.REQUIRED_CASE_FINAL_TIME
    assert stage_i.parser().parse_args([
        "bundle-campaign",
    ]).required_final_time == stage_i.REQUIRED_CASE_FINAL_TIME
    for value in ("nan", "inf", "0", "-1"):
        with pytest.raises(ValueError, match="positive and finite"):
            stage_i.require_positive_finite_float(value, "--required-final-time")
        with pytest.raises(SystemExit):
            stage_i.parser().parse_args([
                "bundle-case", "--case-id", "R02",
                "--required-final-time", value,
            ])
        with pytest.raises(SystemExit):
            stage_i.parser().parse_args([
                "bundle-campaign", "--required-final-time", value,
            ])
        direct_args = SimpleNamespace(
            root=str(root), allow_local_root=True, required_final_time=value
        )
        with pytest.raises(ValueError, match="positive and finite"):
            stage_i.bundle_case(direct_args)
        with pytest.raises(ValueError, match="positive and finite"):
            stage_i.bundle_campaign(direct_args)

    def reservation(case_id, segment, nodes):
        return {
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "manifest": str(
                paths["runs"] / case_id / segment
                / "manifest" / "prepared_run.json"
            ),
            "case_id": case_id,
            "case_name": "case",
            "segment": segment,
            "nodes": nodes,
            "requested_walltime": "00:10:00",
            "reserved_node_hours": nodes / 6.0,
            "state": "prepared",
            "prepared_utc": stage_i.utc_now(),
            "execution_intent_sha256": "a" * 64,
        }

    local_r17 = reservation("R17", "s_local_shape", 1)
    stage_i.validate_reservation_record(paths, local_r17)
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        stage_i.validate_reservation_record(
            paths, reservation("R02", "s_canonical_shape", 1)
        )
        stage_i.validate_reservation_record(
            paths, reservation("R17", "s_canonical_shape", 8)
        )
        with pytest.raises(ValueError, match="R02 canonical Stage I preparation"):
            stage_i.validate_reservation_record(
                paths, reservation("R02", "s_wrong_shape", 8)
            )
        with pytest.raises(ValueError, match="R17 canonical Stage I preparation"):
            stage_i.validate_reservation_record(paths, local_r17)
    replay_manifest_path = (
        paths["runs"] / "R02" / "s_replay_mismatch"
        / "manifest" / "prepared_run.json"
    )
    replay_manifest = {
        "schema_version": 3,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "project_root": str(root),
        "state": "prepared",
        "policy": {},
        "run": {
            "case_id": "R02",
            "case_name": "case",
            "segment": "s_replay_mismatch",
        },
        "allocation": {
            "nodes": 1,
            "requested_walltime": "00:20:00",
            "requested_seconds": 1200,
            "reserved_node_hours": 1.0 / 3.0,
            "ranks_per_node": 8,
            "cpus_per_task": 7,
        },
        "command": {"athena_walltime": "00:10:00"},
        "paths": {},
    }
    replay_reservation = reservation("R02", "s_replay_mismatch", 1)
    replay_reservation["execution_intent_sha256"] = (
        stage_i.execution_intent_sha256(replay_manifest)
    )
    replay_transaction = paths["transactions"] / "replay-mismatch.json"
    stage_i.write_json(replay_transaction, {
        "schema_version": 1,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "transaction_id": replay_transaction.stem,
        "kind": "prepared",
        "created_utc": stage_i.utc_now(),
        "manifest_path": str(replay_manifest_path),
        "prior_reservations": [],
        "prior_reservations_sha256": stage_i.stable_json_sha256([]),
        "manifest": replay_manifest,
        "reservations": [replay_reservation],
        "ledger_row": None,
    }, mode=0o644)
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        with pytest.raises(
            ValueError, match="reservation allocation differs from manifest"
        ):
            stage_i.read_transaction(paths, replay_transaction)
    assert not replay_manifest_path.exists()
    assert json.loads(paths["reservations"].read_text()) == []
    replay_transaction.unlink()
    r17_replay_manifest_path = (
        paths["runs"] / "R17" / "s_replay_before_predecessors"
        / "manifest" / "prepared_run.json"
    )
    r17_replay_manifest = json.loads(json.dumps(replay_manifest))
    r17_replay_manifest["run"] = {
        "case_id": "R17",
        "case_name": "case",
        "segment": "s_replay_before_predecessors",
    }
    r17_replay_manifest["allocation"]["nodes"] = 8
    r17_replay_manifest["allocation"]["reserved_node_hours"] = 8.0 / 3.0
    r17_replay_reservation = reservation(
        "R17", "s_replay_before_predecessors", 8
    )
    r17_replay_reservation["requested_walltime"] = "00:20:00"
    r17_replay_reservation["reserved_node_hours"] = 8.0 / 3.0
    r17_replay_reservation["execution_intent_sha256"] = (
        stage_i.execution_intent_sha256(r17_replay_manifest)
    )
    r17_replay_transaction = paths["transactions"] / "r17-replay-policy.json"
    stage_i.write_json(r17_replay_transaction, {
        "schema_version": 1,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "transaction_id": r17_replay_transaction.stem,
        "kind": "prepared",
        "created_utc": stage_i.utc_now(),
        "manifest_path": str(r17_replay_manifest_path),
        "prior_reservations": [],
        "prior_reservations_sha256": stage_i.stable_json_sha256([]),
        "manifest": r17_replay_manifest,
        "reservations": [r17_replay_reservation],
        "ledger_row": None,
    }, mode=0o644)
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        policy.setattr(
            stage_i, "accepted_case_lineage", lambda _paths, _case_id: []
        )
        with pytest.raises(ValueError, match="R17 must remain last"):
            stage_i.read_transaction(paths, r17_replay_transaction)
    assert not r17_replay_manifest_path.exists()
    assert json.loads(paths["reservations"].read_text()) == []
    r17_replay_transaction.unlink()
    recorded_replay_manifest_path = (
        paths["runs"] / "R02" / "s_recorded_nan_replay"
        / "manifest" / "prepared_run.json"
    )
    recorded_replay_manifest = json.loads(json.dumps(replay_manifest))
    recorded_replay_manifest["state"] = "recorded"
    recorded_replay_manifest["job_id"] = "123"
    recorded_replay_manifest["run"]["segment"] = "s_recorded_nan_replay"
    recorded_ledger_row = {column: "" for column in stage_i.LEDGER_COLUMNS}
    recorded_ledger_row.update({
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "job_id": "123",
        "case_id": "R02",
        "case_name": "case",
        "segment": "s_recorded_nan_replay",
        "nodes": "1",
        "requested_walltime": "00:20:00",
        "elapsed_seconds": "1",
        "reserved_node_hours": "0.333333",
        "actual_node_hours": "nan",
        "cumulative_stage_i_node_hours": "0.1",
        "executable_revision": None,
        "executable_sha256": None,
        "input_revision": None,
        "input_file": None,
        "output_dir": None,
        "result": "failed",
    })
    recorded_replay_manifest["accounting"] = recorded_ledger_row
    recorded_replay_reservation = reservation(
        "R02", "s_recorded_nan_replay", 1
    )
    recorded_replay_reservation["requested_walltime"] = "00:20:00"
    recorded_replay_reservation["reserved_node_hours"] = 1.0 / 3.0
    recorded_replay_reservation["state"] = "recorded"
    recorded_replay_reservation["job_id"] = "123"
    recorded_replay_reservation["actual_node_hours"] = 0.1
    recorded_replay_reservation["result"] = "failed"
    recorded_replay_reservation["execution_intent_sha256"] = (
        stage_i.execution_intent_sha256(recorded_replay_manifest)
    )
    submitted_replay_reservation = dict(recorded_replay_reservation)
    submitted_replay_reservation["state"] = "submitted"
    submitted_replay_reservation.pop("actual_node_hours")
    submitted_replay_reservation.pop("result")
    recorded_replay_transaction = (
        paths["transactions"] / "recorded-nan-replay.json"
    )
    stage_i.write_json(recorded_replay_transaction, {
        "schema_version": 1,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "transaction_id": recorded_replay_transaction.stem,
        "kind": "recorded",
        "created_utc": stage_i.utc_now(),
        "manifest_path": str(recorded_replay_manifest_path),
        "prior_reservations": [submitted_replay_reservation],
        "prior_reservations_sha256": stage_i.stable_json_sha256([
            submitted_replay_reservation
        ]),
        "manifest": recorded_replay_manifest,
        "reservations": [recorded_replay_reservation],
        "ledger_row": recorded_ledger_row,
    }, mode=0o644)
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        with pytest.raises(ValueError, match="invalid numeric fields"):
            stage_i.read_transaction(paths, recorded_replay_transaction)
    assert not recorded_replay_manifest_path.exists()
    assert json.loads(paths["reservations"].read_text()) == []
    recorded_replay_transaction.unlink()
    recorded_ledger_row["elapsed_seconds"] = "3600"
    recorded_ledger_row["actual_node_hours"] = "0.000000"
    recorded_ledger_row["cumulative_stage_i_node_hours"] = "0.000000"
    recorded_replay_manifest["accounting"] = recorded_ledger_row
    recorded_replay_reservation["actual_node_hours"] = 0.0
    recorded_replay_reservation["execution_intent_sha256"] = (
        stage_i.execution_intent_sha256(recorded_replay_manifest)
    )
    submitted_replay_reservation = dict(recorded_replay_reservation)
    submitted_replay_reservation["state"] = "submitted"
    submitted_replay_reservation.pop("actual_node_hours")
    submitted_replay_reservation.pop("result")
    stage_i.write_json(recorded_replay_transaction, {
        "schema_version": 1,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "transaction_id": recorded_replay_transaction.stem,
        "kind": "recorded",
        "created_utc": stage_i.utc_now(),
        "manifest_path": str(recorded_replay_manifest_path),
        "prior_reservations": [submitted_replay_reservation],
        "prior_reservations_sha256": stage_i.stable_json_sha256([
            submitted_replay_reservation
        ]),
        "manifest": recorded_replay_manifest,
        "reservations": [recorded_replay_reservation],
        "ledger_row": recorded_ledger_row,
    }, mode=0o644)
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        with pytest.raises(ValueError, match="invalid numeric fields"):
            stage_i.read_transaction(paths, recorded_replay_transaction)
    assert not recorded_replay_manifest_path.exists()
    assert json.loads(paths["reservations"].read_text()) == []
    recorded_replay_transaction.unlink()
    recorded_ledger_row["elapsed_seconds"] = "0"
    recorded_ledger_row["actual_node_hours"] = "0.000000"
    recorded_ledger_row["cumulative_stage_i_node_hours"] = "0.000000"
    recorded_ledger_row["state"] = "COMPLETED"
    recorded_ledger_row["exit_code"] = "0:0"
    recorded_ledger_row["submitted_utc"] = "2026-05-30T00:00:00"
    recorded_ledger_row["completed_utc"] = "2026-05-30T01:00:00"
    recorded_replay_manifest["accounting"] = recorded_ledger_row
    recorded_replay_reservation["execution_intent_sha256"] = (
        stage_i.execution_intent_sha256(recorded_replay_manifest)
    )
    submitted_replay_reservation = dict(recorded_replay_reservation)
    submitted_replay_reservation["state"] = "submitted"
    submitted_replay_reservation.pop("actual_node_hours")
    submitted_replay_reservation.pop("result")
    scheduler_evidence = paths["accounting"] / "123.stage_i.sacct.txt"
    scheduler_evidence.write_text(
        "123|"
        + stage_i.expected_job_name(recorded_replay_manifest)
        + "|COMPLETED|0:0|1|3600|2026-05-30T00:00:00|2026-05-30T01:00:00\n"
    )
    stage_i.write_json(recorded_replay_transaction, {
        "schema_version": 1,
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "transaction_id": recorded_replay_transaction.stem,
        "kind": "recorded",
        "created_utc": stage_i.utc_now(),
        "manifest_path": str(recorded_replay_manifest_path),
        "prior_reservations": [submitted_replay_reservation],
        "prior_reservations_sha256": stage_i.stable_json_sha256([
            submitted_replay_reservation
        ]),
        "manifest": recorded_replay_manifest,
        "reservations": [recorded_replay_reservation],
        "ledger_row": recorded_ledger_row,
    }, mode=0o644)
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        with pytest.raises(ValueError, match="scheduler evidence elapsed_seconds"):
            stage_i.read_transaction(paths, recorded_replay_transaction)
    assert not recorded_replay_manifest_path.exists()
    assert json.loads(paths["reservations"].read_text()) == []
    recorded_replay_transaction.unlink()
    scheduler_evidence.unlink()
    reached_policy = []
    with monkeypatch.context() as policy:
        policy.setattr(stage_i, "DEFAULT_ROOT", root)
        policy.setattr(
            stage_i, "require_reconciled_store_consistency", lambda *_args: None
        )

        def stop_at_policy(*args):
            reached_policy.append(args)
            raise RuntimeError("prepare policy sentinel")

        policy.setattr(stage_i, "require_prepare_case_policy", stop_at_policy)
        with pytest.raises(RuntimeError, match="prepare policy sentinel"):
            stage_i.prepare(SimpleNamespace(
                root=str(root),
                allow_local_root=False,
                case_id="R02",
                segment="s_prepare_policy_route",
                nodes=1,
            ))
    assert reached_policy == [(paths, "R02", 1, False)]

    manifest_path = (
        paths["runs"] / "R16" / "s00" / "manifest" / "prepared_run.json"
    )
    manifest_path.parent.mkdir(parents=True)
    batch_script = manifest_path.parent / "cgl_lf_stage_i.sbatch"
    batch_template = (
        "#!/bin/bash\n"
        f"BATCH_SCRIPT_SHA256={stage_i.BATCH_SCRIPT_DIGEST_PLACEHOLDER}\n"
    )
    batch_text, batch_digest = stage_i.finalize_batch_script(batch_template)
    batch_script.write_text(batch_text)
    batch_script.chmod(0o750)
    monkeypatch.setattr(
        stage_i, "generated_batch_script", lambda *_args: batch_template
    )
    stage_i.write_json(manifest_path, {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "project_root": str(root),
        "state": "prepared",
        "command": {"batch_script_sha256": batch_digest},
        "paths": {"batch_script": str(batch_script)},
    })
    stage_i.write_json(paths["reservations"], [{
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "manifest": str(manifest_path),
        "case_id": "R16",
        "case_name": "case",
        "segment": "s00",
        "nodes": 1,
        "requested_walltime": "00:10:00",
        "reserved_node_hours": 1.0 / 6.0,
        "state": "prepared",
    }])
    queue = tmp_path / "squeue.txt"
    queue.write_text("")
    args = SimpleNamespace(
        manifest=str(manifest_path),
        allow_local_root=True,
        squeue_file=str(queue),
        skip_slurm_test=True,
        allow_shared_root_campaign=[],
    )
    assert stage_i.check_submit(args) == 0
    assert "--allow-shared-root-campaign" not in capsys.readouterr().out

    queue.write_text("123|batch|RUNNING|unrelated_job\n")
    assert stage_i.check_submit(args) == 0
    capsys.readouterr()
    overlap_manifest_path = (
        paths["runs"] / "R04" / "s00" / "manifest" / "prepared_run.json"
    )
    overlap_manifest_path.parent.mkdir(parents=True)
    overlap_manifest = {
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "project_root": str(root),
        "state": "submitted",
        "job_id": "321",
        "run": {
            "case_id": "R04",
            "case_name": "overlap",
            "segment": "s00",
        },
        "allocation": {
            "nodes": 2,
            "requested_walltime": "00:10:00",
            "reserved_node_hours": 1.0 / 3.0,
        },
    }
    stage_i.write_json(overlap_manifest_path, overlap_manifest)
    reservations = json.loads(paths["reservations"].read_text())
    reservations.append({
        "execution_epoch": stage_i.EXECUTION_EPOCH,
        "manifest": str(overlap_manifest_path),
        "case_id": "R04",
        "case_name": "overlap",
        "segment": "s00",
        "nodes": 2,
        "requested_walltime": "00:10:00",
        "reserved_node_hours": 1.0 / 3.0,
        "state": "submitted",
        "job_id": "321",
        "prepared_utc": stage_i.utc_now(),
    })
    stage_i.write_json(paths["reservations"], reservations)
    overlap_queue_row = (
        f"321|batch|RUNNING|{stage_i.expected_job_name(overlap_manifest)}\n"
    )
    queue.write_text(overlap_queue_row)
    assert stage_i.check_submit(args) == 0
    capsys.readouterr()
    queue.write_text(overlap_queue_row + "123|batch|RUNNING|unrelated_job\n")
    assert stage_i.check_submit(args) == 0
    capsys.readouterr()
    queue.write_text(overlap_queue_row)
    for invalid_row in (
        f"321|debug|RUNNING|{stage_i.expected_job_name(overlap_manifest)}",
        "321|batch|RUNNING|wrong_name",
        overlap_queue_row.strip() + "\n" + overlap_queue_row.strip(),
    ):
        with pytest.raises(ValueError, match="queued CGL job"):
            stage_i.authenticate_production_queue(
                paths, reservations, invalid_row.splitlines()
            )

    shared_manifest = (
        root / "runs" / "exploratory" / "manifest" / "prepared_run.json"
    )
    shared_manifest.parent.mkdir(parents=True)
    stage_i.write_json(shared_manifest, {
        "campaign_id": "exploratory",
        "state": "running",
    })
    with pytest.raises(ValueError, match="shared-root campaign records"):
        stage_i.check_submit(args)
    args.allow_shared_root_campaign = [
        "exploratory", "beta 25", "-reviewed", "exploratory",
    ]
    assert stage_i.check_submit(args) == 0
    python = stage_i.authenticated_python_binary()
    helper = str(PAPER_STAGE_I_TOOL.resolve())
    submit_line = next(
        line.strip() for line in capsys.readouterr().out.splitlines()
        if line.strip().startswith(
            f"{shlex.quote(python)} -I -S -B {shlex.quote(helper)} submit "
        )
    )
    tokens = shlex.split(submit_line)
    assert tokens[:5] == [python, "-I", "-S", "-B", helper]
    parsed = stage_i.parser().parse_args(tokens[5:])
    assert parsed.allow_shared_root_campaign == [
        "-reviewed", "beta 25", "exploratory",
    ]

    manifest = json.loads(manifest_path.read_text())
    manifest["execution_epoch"] = "E01-pre-modal-driver"
    stage_i.write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="execution epoch"):
        stage_i.check_submit(args)


def test_cgl_lf_stage_i_bundle_selects_terminal_restart_lineage(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_lineage_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    root = tmp_path / "root"
    paths = stage_i.initialize(root)

    def record_segment(segment, times, result, input_sha, parent=None):
        manifest_path = (
            root / "runs" / "mks24-stage-i" / stage_i.EXECUTION_EPOCH
            / "R16" / segment
            / "manifest" / "prepared_run.json"
        )
        manifest_path.parent.mkdir(parents=True)
        output_dir = manifest_path.parent.parent / "output"
        (output_dir / "bin").mkdir(parents=True)
        input_file = manifest_path.parent / "submitted_input.athinput"
        input_file.write_text("retained input\n")
        rows = "".join(f"{time} 0.001 0\n" for time in times)
        for name in ("case.mhd.hst", "case.user.hst"):
            (output_dir / name).write_text(
                "# [0]=time [1]=dt [2]=clean\n" + rows
            )
        (output_dir / "bin" / f"{segment}.00000.bin").write_bytes(
            (
                "Athena binary output version=1.1\n"
                "  size of preheader=5\n"
                f"  time={times[-1]}\n"
                "  cycle=0\n"
                "  size of location=8\n"
                "  size of variable=4\n"
            ).encode()
        )
        command = {
            "input_sha256": input_sha,
            "input_revision": "r" * 40,
            "input_file": str(input_file),
            "executable_sha256": "x" * 64,
            "executable": "athena",
        }
        if parent is not None:
            command["parent_segment"] = {"manifest": str(parent)}
        stage_i.write_json(manifest_path, {
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "state": "recorded",
            "accounting": {"result": result},
            "command": command,
            "allocation": {"nodes": 1, "ranks_per_node": 1},
            "paths": {"output_dir": str(output_dir)},
            "scientific_inspection": {
                "accepted": result == "accepted",
                "clean_for_continuation": True,
                "final_time": times[-1],
            },
        })
        return manifest_path

    diagnostic = record_segment(
        "s00_legacy", [0.0, 1.6], "clean_partial", "a" * 64
    )
    ranked = record_segment(
        "s01_rankio", [0.0, 1.8], "clean_partial", "b" * 64
    )
    terminal = record_segment(
        "s02_rankio", [1.823, 2.0], "accepted", "b" * 64, parent=ranked
    )
    stage_i.write_json(paths["reservations"], [
        {
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "manifest": str(manifest),
        }
        for manifest in (diagnostic, ranked, terminal)
    ])

    lineage = stage_i.accepted_case_lineage(paths, "R16")
    assert [Path(segment["_manifest_path"]) for segment in lineage] == [
        ranked, terminal
    ]
    assert diagnostic not in [
        Path(segment["_manifest_path"]) for segment in lineage
    ]
    monkeypatch.setattr(stage_i, "validate_matrix", lambda *_: {
        "cases": [{"id": "R16", "name": "case", "input": "ignored"}]
    })
    monkeypatch.setattr(stage_i, "analysis_model_choices", lambda _: {
        "output1_file_type": "hst",
        "output1_dt": "0.02",
    })
    bundle = root / "runs" / "bundles" / "test"
    args = SimpleNamespace(
        root=str(root),
        allow_local_root=True,
        source_dir=str(tmp_path),
        matrix=str(tmp_path / "matrix.json"),
        case_id="R16",
        required_final_time=2.0,
        output_dir=str(bundle),
        replace=False,
    )
    assert stage_i.bundle_case(args) == 0
    manifest = json.loads((bundle / "manifest.json").read_text())
    assert manifest["production_segment_manifests"] == [
        str(ranked), str(terminal)
    ]
    merged = stage_i.parse_history(bundle / "history" / "case.mhd.hst")
    assert merged["time"] == [0.0, 1.8, 1.823, 2.0]
    (terminal.parent.parent / "output" / "case.mhd.hst").write_text(
        "# [0]=time [1]=dt [2]=clean\n1.85 0.001 0\n2.0 0.001 0\n"
    )
    args.output_dir = str(root / "runs" / "bundles" / "gap")
    with pytest.raises(ValueError, match="configured sampling cadence"):
        stage_i.bundle_case(args)


def test_cgl_lf_stage_i_qualification_token_binds_corrected_build(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_qualification_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    root = tmp_path / "root"
    paths = stage_i.initialize(root)
    executable = root / "build" / "athena"
    executable.parent.mkdir()
    executable.write_bytes(b"corrected executable\n")
    executable.chmod(0o755)
    executable_sha = hashlib.sha256(executable.read_bytes()).hexdigest()
    revision = "a" * 40
    build_manifest = root / "build-manifest"
    build_manifest.mkdir()
    (build_manifest / "athena.sha256").write_text(f"{executable_sha}  athena\n")
    (build_manifest / "environment.txt").write_text(
        f"git_revision={revision}\n"
    )
    assert stage_i.require_qualification_approval(
        paths, executable_sha, revision, offline_local_root=True
    ) is None
    with pytest.raises(ValueError, match="token is absent"):
        stage_i.require_qualification_approval(
            paths, executable_sha, revision, offline_local_root=False
        )

    approval_args = SimpleNamespace(
        root=str(root),
        allow_local_root=True,
        executable=str(executable),
        build_manifest=str(build_manifest),
        approved_by="test-reviewer",
        review_notes="reviewed corrected-build Frontier qualification",
        confirm_corrected_build_frontier_qualified=True,
        replace_existing_approval=False,
    )
    assert stage_i.approve_qualification(approval_args) == 0
    qualification = stage_i.require_qualification_approval(
        paths, executable_sha, revision, offline_local_root=False
    )
    assert qualification is not None
    assert qualification["sha256"] == stage_i.sha256(paths["qualification"])
    assert "approved by `test-reviewer`" in paths["summary"].read_text()
    with pytest.raises(ValueError, match="does not approve this executable"):
        stage_i.require_qualification_approval(
            paths, "b" * 64, revision, offline_local_root=False
        )
    with pytest.raises(ValueError, match="does not approve this git revision"):
        stage_i.require_qualification_approval(
            paths, executable_sha, "b" * 40, offline_local_root=False
        )

    manifest_path = (
        paths["runs"] / "R02" / "s00" / "manifest" / "prepared_run.json"
    )
    manifest_path.parent.mkdir(parents=True)
    input_path = manifest_path.parent / "submitted_input.athinput"
    input_path.write_text("<time>\ntlim = 10.0\n")
    matrix_path = manifest_path.parent / "mks24_stage_i_manifest.json"
    matrix_path.write_text("{}\n")
    manifest = {
        "run": {
            "case_id": "R02",
            "segment": "s00",
            "run_basename": "test",
        },
        "allocation": {
            "requested_walltime": "00:10:00",
            "nodes": 1,
            "ranks_per_node": 8,
            "cpus_per_task": 7,
        },
        "command": {
            "overrides": ["time/tlim=1.0"],
            "athena_walltime": "00:09:00",
            "executable": str(executable),
            "executable_sha256": executable_sha,
            "input_file": str(input_path),
            "input_sha256": stage_i.sha256(input_path),
            "matrix_file": str(matrix_path),
            "matrix_sha256": stage_i.sha256(matrix_path),
            "production_utility": {
                "path": str(PAPER_STAGE_I_TOOL.resolve()),
                "sha256": stage_i.sha256(PAPER_STAGE_I_TOOL.resolve()),
            },
            "qualification_approval": qualification,
            "restart_files": [],
        },
        "paths": {
            "slurm_log": str(root / "slurm.log"),
            "output_dir": str(root / "output"),
            "environment_log": str(root / "environment.log"),
        },
    }
    script = stage_i.generated_batch_script(manifest, manifest_path)
    assert "qualification_approval" in script
    assert qualification["sha256"] in script
    assert script.index("qualification_approval") < script.index("srun -N")

    stage_i.write_json(manifest_path, {"state": "prepared"})
    approval_args.replace_existing_approval = True
    with pytest.raises(ValueError, match="cannot change after segment preparation"):
        stage_i.approve_qualification(approval_args)


def test_cgl_lf_stage_i_recovers_ambiguous_atomic_submit(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_submit_recovery_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)
    batch_template = (
        "#!/bin/bash\n"
        f"BATCH_SCRIPT_SHA256={stage_i.BATCH_SCRIPT_DIGEST_PLACEHOLDER}\n"
    )
    monkeypatch.setattr(
        stage_i, "generated_batch_script", lambda *_args: batch_template
    )

    def prepared_fixture(root, segment):
        paths = stage_i.initialize(root)
        manifest_path = (
            paths["runs"] / "R16" / segment / "manifest" / "prepared_run.json"
        )
        manifest_path.parent.mkdir(parents=True)
        batch_script = manifest_path.parent / "cgl_lf_stage_i.sbatch"
        batch_text, batch_digest = stage_i.finalize_batch_script(batch_template)
        batch_script.write_text(batch_text)
        batch_script.chmod(0o750)
        manifest = {
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "project_root": str(root),
            "state": "prepared",
            "run": {
                "case_id": "R16",
                "case_name": "case",
                "segment": segment,
            },
            "allocation": {
                "nodes": 1,
                "requested_walltime": "00:10:00",
                "reserved_node_hours": 1.0 / 6.0,
            },
            "command": {"batch_script_sha256": batch_digest},
            "paths": {"batch_script": str(batch_script)},
        }
        stage_i.write_json(manifest_path, manifest)
        stage_i.write_json(paths["reservations"], [{
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "manifest": str(manifest_path),
            "case_id": "R16",
            "case_name": "case",
            "segment": segment,
            "nodes": 1,
            "requested_walltime": "00:10:00",
            "reserved_node_hours": 1.0 / 6.0,
            "state": "prepared",
            "prepared_utc": stage_i.utc_now(),
        }])
        queue = root / "squeue.txt"
        queue.write_text("")
        malformed = root / "sbatch.txt"
        malformed.write_text("scheduler output unavailable\n")
        args = SimpleNamespace(
            manifest=str(manifest_path),
            allow_local_root=True,
            squeue_file=str(queue),
            skip_slurm_test=True,
            allow_shared_root_campaign=[],
            sbatch_output_file=str(malformed),
        )
        return paths, manifest_path, args

    paths, manifest_path, args = prepared_fixture(tmp_path / "recover", "s00")
    before_manifest = manifest_path.read_text()
    before_reservations = paths["reservations"].read_text()
    assert stage_i.check_submit(args) == 0
    assert manifest_path.read_text() == before_manifest
    assert paths["reservations"].read_text() == before_reservations
    with pytest.raises(ValueError, match="did not return one numeric job ID"):
        stage_i.submit(args)
    assert len(stage_i.pending_transaction_paths(paths)) == 1
    pending = stage_i.read_transaction(
        paths, stage_i.pending_transaction_paths(paths)[0]
    )
    assert "initial_queue_authentication" in pending["submission_audit"]
    assert "final_queue_authentication" in pending["submission_audit"]
    assert json.loads(manifest_path.read_text())["state"] == "prepared"
    assert stage_i.recover_submit(SimpleNamespace(
        manifest=str(manifest_path),
        allow_local_root=True,
        job_id="12345",
    )) == 0
    assert not stage_i.pending_transaction_paths(paths)
    assert json.loads(manifest_path.read_text())["state"] == "submitted"

    paths, manifest_path, args = prepared_fixture(tmp_path / "clear", "s01")
    with pytest.raises(ValueError, match="did not return one numeric job ID"):
        stage_i.submit(args)
    assert stage_i.clear_submit_pending(SimpleNamespace(
        manifest=str(manifest_path),
        allow_local_root=True,
        confirm_no_job_submitted=True,
        notes="reviewed scheduler and confirmed no submitted job",
    )) == 0
    cleared = json.loads(manifest_path.read_text())
    assert cleared["state"] == "prepared"
    assert cleared["submission_recovery_notes"]
    assert not stage_i.pending_transaction_paths(paths)

    paths, manifest_path, args = prepared_fixture(tmp_path / "race", "s02")
    Path(args.sbatch_output_file).write_text("12345\n")
    before_manifest = manifest_path.read_text()
    before_reservations = paths["reservations"].read_text()
    queue_results = iter(["", "999|batch|RUNNING|cgl_untracked_stage_i_job\n"])
    monkeypatch.setattr(
        stage_i,
        "production_queue_output",
        lambda *_args, **_kwargs: next(queue_results),
    )
    with pytest.raises(ValueError, match="queued CGL job"):
        stage_i.submit(args)
    assert not stage_i.pending_transaction_paths(paths)
    assert manifest_path.read_text() == before_manifest
    assert paths["reservations"].read_text() == before_reservations


def test_cgl_lf_stage_i_cancels_only_terminal_never_started_submitted_job(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_submitted_cancel_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    def fixture(root, scheduler_state="CANCELLED by 1234", elapsed="0",
                allocated_nodes="0"):
        paths = stage_i.initialize(root)
        manifest_path = (
            paths["runs"] / "R16" / "s_cancel"
            / "manifest" / "prepared_run.json"
        )
        manifest_path.parent.mkdir(parents=True)
        batch_script = manifest_path.parent / "cgl_lf_stage_i.sbatch"
        batch_script.write_text(
            "#!/bin/bash\nBATCH_SCRIPT_SHA256=" + "0" * 64 + "\n"
        )
        source_bundle = root / "source-archives" / "corrupt.bundle"
        source_bundle.parent.mkdir()
        source_bundle.write_text("corrupt\n")
        manifest = {
            "schema_version": 3,
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "project_root": str(root),
            "state": "submitted",
            "job_id": "12345",
            "policy": {},
            "run": {
                "case_id": "R16",
                "case_name": "case",
                "segment": "s_cancel",
            },
            "allocation": {
                "nodes": 1,
                "requested_walltime": "00:10:00",
                "reserved_node_hours": 1.0 / 6.0,
            },
            "command": {
                "batch_script_sha256": stage_i.normalized_batch_script_sha256(
                    batch_script
                ),
                "source_bundle": {
                    "path": str(source_bundle),
                    "sha256": "a" * 64,
                    "verified_revisions": ["b" * 40],
                },
            },
            "paths": {
                "batch_script": str(batch_script),
                "output_dir": str(manifest_path.parents[1] / "output"),
            },
        }
        stage_i.write_json(manifest_path, manifest)
        reservation = {
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "manifest": str(manifest_path),
            "case_id": "R16",
            "case_name": "case",
            "segment": "s_cancel",
            "nodes": 1,
            "requested_walltime": "00:10:00",
            "reserved_node_hours": 1.0 / 6.0,
            "state": "submitted",
            "prepared_utc": stage_i.utc_now(),
            "job_id": "12345",
            "execution_intent_sha256": stage_i.execution_intent_sha256(manifest),
        }
        stage_i.write_json(paths["reservations"], [reservation])
        incident = paths["accounting"] / "incidents" / "incident.json"
        incident.parent.mkdir()
        stage_i.write_json(incident, {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-corruption-incident",
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "expected": {"sha256": "a" * 64},
            "observed": {
                "path": str(source_bundle),
                "sha256": stage_i.sha256(source_bundle),
            },
            "scheduler_containment": {
                "job_id": "12345",
                "reason": "JobHeldUser",
                "released": False,
            },
        })
        incident.chmod(0o444)
        evidence = stage_i.submitted_cancellation_evidence_paths(paths, "12345")
        submitted_snapshot = evidence["submitted_manifest"]
        stage_i.write_json(submitted_snapshot, manifest)
        submitted_snapshot.chmod(0o444)
        scheduler_batch_script = evidence["scheduler_batch_script"]
        scheduler_batch_script.write_bytes(batch_script.read_bytes())
        scheduler_batch_script.chmod(0o444)
        hold = evidence["pre_cancel_hold"]
        hold.write_text(
            f"12345|{stage_i.expected_job_name(manifest)}|PENDING|0:00|1|"
            "(JobHeldUser)\n"
        )
        hold.chmod(0o444)
        authorization = evidence["authorization"]
        stage_i.write_json(authorization, {
            "schema_version": 1,
            "record_type":
                "stage-i-source-bundle-recovery-cancellation-authorization",
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "job_id": "12345",
            "decision": "cancel-held-job-without-start",
            "reason": "the immutable launch source bundle is corrupt",
            "manifest": {
                "live_path": str(manifest_path),
                "snapshot_path": str(submitted_snapshot),
                "sha256": stage_i.sha256(submitted_snapshot),
            },
            "batch_script": {
                "path": str(batch_script),
                "sha256": stage_i.sha256(batch_script),
                "normalized_sha256": stage_i.normalized_batch_script_sha256(
                    batch_script
                ),
                "scheduler_snapshot_path": str(scheduler_batch_script),
                "scheduler_snapshot_sha256": stage_i.sha256(
                    scheduler_batch_script
                ),
            },
            "pre_cancel_hold": {
                "path": str(hold),
                "sha256": stage_i.sha256(hold),
            },
            "source_bundle": {
                "path": str(source_bundle),
                "expected_sha256": "a" * 64,
                "observed_sha256": stage_i.sha256(source_bundle),
                "incident": {
                    "path": str(incident),
                    "sha256": stage_i.sha256(incident),
                },
            },
        })
        authorization.chmod(0o444)
        audit = evidence["publication_audit"]
        stage_i.write_json(audit, {
            "schema_version": 1,
            "record_type":
                "stage-i-source-bundle-recovery-cancellation-publication-audit",
            "execution_epoch": stage_i.EXECUTION_EPOCH,
            "published_utc": stage_i.utc_now(),
            "authorization": {
                "path": str(authorization),
                "sha256": stage_i.sha256(authorization),
                "mode": "0444",
            },
            "review": {
                "status": "approved",
                "authority": "independent-source-bundle-recovery-review",
                "reviewed_by": "test reviewer",
            },
        })
        audit.chmod(0o444)
        sacct = evidence["post_cancel_sacct"]
        sacct.write_text(
            f"12345|{stage_i.expected_job_name(manifest)}|{scheduler_state}|0:0|"
            f"{allocated_nodes}|{elapsed}|2026-06-05T00:00:00|"
            "2026-06-05T00:01:00|\n"
        )
        sacct.chmod(0o444)
        squeue = evidence["post_cancel_queue"]
        squeue.write_text("")
        squeue.chmod(0o444)
        args = SimpleNamespace(
            manifest=str(manifest_path),
            allow_local_root=True,
            job_id="12345",
            notes="cancelled before start for formal source-bundle recovery",
            sacct_file=str(sacct),
            squeue_file=str(squeue),
        )
        return paths, manifest_path, args

    paths, manifest_path, args = fixture(tmp_path / "accepted")
    assert stage_i.cancel_submitted(args) == 0
    cancelled = json.loads(manifest_path.read_text())
    assert cancelled["state"] == "cancelled"
    assert cancelled["cancellation"]["mode"] == "authenticated-submitted-no-start"
    assert cancelled["cancellation"]["scheduler"]["state"] == "CANCELLED"
    assert json.loads(paths["reservations"].read_text())[0]["state"] == "cancelled"
    assert not stage_i.pending_transaction_paths(paths)
    with monkeypatch.context() as context:
        context.setattr(
            stage_i,
            "authenticate_prepared_execution",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("cancelled-submitted reconcile used launch authentication")
            ),
        )
        report = stage_i.reconcile_report(tmp_path / "accepted")
    assert report["consistent"]
    assert report["counts"]["active_reservations"] == 0

    root = tmp_path / "accepted"
    evidence = stage_i.submitted_cancellation_evidence_paths(paths, "12345")
    source_bundle = root / "source-archives" / "corrupt.bundle"
    incident = paths["accounting"] / "incidents" / "incident.json"
    with monkeypatch.context() as context:
        context.setattr(stage_i, "DEFAULT_ROOT", root.resolve())
        for key, value in {
            "job_id": "12345",
            "manifest_relative": str(manifest_path.relative_to(root)),
            "source_bundle_relative": str(source_bundle.relative_to(root)),
            "source_bundle_expected_sha256": "a" * 64,
            "source_bundle_observed_sha256": stage_i.sha256(source_bundle),
            "incident_relative": str(incident.relative_to(root)),
            "incident_sha256": stage_i.sha256(incident),
            "authorization_sha256": stage_i.sha256(evidence["authorization"]),
            "publication_audit_sha256": stage_i.sha256(
                evidence["publication_audit"]
            ),
        }.items():
            context.setitem(stage_i.CANONICAL_SOURCE_BUNDLE_RECOVERY, key, value)
        context.setattr(
            stage_i,
            "live_cancellation_sacct_output",
            lambda *_args: (_ for _ in ()).throw(
                AssertionError("reconcile queried live historical sacct")
            ),
        )
        context.setattr(
            stage_i,
            "live_cancelled_job_queue_output",
            lambda *_args: (_ for _ in ()).throw(
                AssertionError("reconcile queried live historical squeue")
            ),
        )
        context.setattr(
            stage_i,
            "qualification_approval_status",
            lambda *_args: {"state": "approved"},
        )
        report = stage_i.reconcile_report(root)
    assert report["consistent"], report["issues"]

    audit = stage_i.submitted_cancellation_evidence_paths(
        paths, "12345"
    )["publication_audit"]
    audit.chmod(0o644)
    report = stage_i.reconcile_report(tmp_path / "accepted")
    assert not report["consistent"]
    assert any(
        "submitted cancellation evidence drift" in issue
        for issue in report["issues"]
    )

    paths, manifest_path, args = fixture(
        tmp_path / "running", scheduler_state="CANCELLED by 1234", elapsed="1"
    )
    with pytest.raises(ValueError, match="terminal CANCELLED no-start"):
        stage_i.cancel_submitted(args)
    assert json.loads(manifest_path.read_text())["state"] == "submitted"
    assert json.loads(paths["reservations"].read_text())[0]["state"] == "submitted"

    paths, manifest_path, args = fixture(
        tmp_path / "allocated", scheduler_state="CANCELLED by 1234",
        allocated_nodes="1"
    )
    with pytest.raises(ValueError, match="terminal CANCELLED no-start"):
        stage_i.cancel_submitted(args)
    assert json.loads(manifest_path.read_text())["state"] == "submitted"

    paths, manifest_path, args = fixture(tmp_path / "queued")
    Path(args.squeue_file).chmod(0o644)
    Path(args.squeue_file).write_text(
        "12345|batch|PENDING|"
        + stage_i.expected_job_name(json.loads(manifest_path.read_text()))
        + "\n"
    )
    Path(args.squeue_file).chmod(0o444)
    with pytest.raises(ValueError, match="remains present in squeue"):
        stage_i.cancel_submitted(args)
    assert json.loads(manifest_path.read_text())["state"] == "submitted"

    paths, manifest_path, args = fixture(tmp_path / "scheduler-script-drift")
    evidence = stage_i.submitted_cancellation_evidence_paths(paths, "12345")
    authorization = json.loads(evidence["authorization"].read_text())
    Path(authorization["batch_script"]["scheduler_snapshot_path"]).chmod(0o644)
    Path(authorization["batch_script"]["scheduler_snapshot_path"]).write_text(
        "#!/bin/bash\nchanged\n"
    )
    with pytest.raises(
        ValueError, match="scheduler-retained submitted batch script"
    ):
        stage_i.cancel_submitted(args)
    assert json.loads(manifest_path.read_text())["state"] == "submitted"

    paths, manifest_path, args = fixture(tmp_path / "submitted-snapshot-drift")
    evidence = stage_i.submitted_cancellation_evidence_paths(paths, "12345")
    authorization = json.loads(evidence["authorization"].read_text())
    Path(authorization["manifest"]["snapshot_path"]).chmod(0o644)
    Path(authorization["manifest"]["snapshot_path"]).write_text("{}\n")
    with pytest.raises(ValueError, match="submitted manifest snapshot"):
        stage_i.cancel_submitted(args)
    assert json.loads(manifest_path.read_text())["state"] == "submitted"

    paths, manifest_path, args = fixture(tmp_path / "job-id-mismatch")
    assert stage_i.cancel_submitted(args) == 0
    reservations = json.loads(paths["reservations"].read_text())
    reservations[0]["job_id"] = "99999"
    stage_i.write_json(paths["reservations"], reservations)
    report = stage_i.reconcile_report(tmp_path / "job-id-mismatch")
    assert not report["consistent"]
    assert any("reservation job ID differs from manifest" in issue
               for issue in report["issues"])

    paths, manifest_path, args = fixture(tmp_path / "job-id-erasure")
    assert stage_i.cancel_submitted(args) == 0
    manifest = json.loads(manifest_path.read_text())
    manifest.pop("job_id")
    stage_i.write_json(manifest_path, manifest)
    reservations = json.loads(paths["reservations"].read_text())
    reservations[0].pop("job_id")
    stage_i.write_json(paths["reservations"], reservations)
    report = stage_i.reconcile_report(tmp_path / "job-id-erasure")
    assert not report["consistent"]
    assert any(
        "submitted cancellation evidence drift" in issue
        for issue in report["issues"]
    )

    root = tmp_path / "prepared-cancel"
    paths, manifest_path, _ = fixture(root)
    manifest = json.loads(manifest_path.read_text())
    manifest["state"] = "prepared"
    manifest.pop("job_id")
    stage_i.write_json(manifest_path, manifest)
    reservations = json.loads(paths["reservations"].read_text())
    reservations[0]["state"] = "prepared"
    reservations[0].pop("job_id")
    reservations[0]["execution_intent_sha256"] = (
        stage_i.execution_intent_sha256(manifest)
    )
    stage_i.write_json(paths["reservations"], reservations)
    with monkeypatch.context() as context:
        context.setattr(
            stage_i, "require_reconciled_store_consistency", lambda *_args: None
        )
        context.setattr(
            stage_i, "authenticate_prepared_execution", lambda *_args, **_kwargs: None
        )
        assert stage_i.cancel(SimpleNamespace(
            manifest=str(manifest_path),
            allow_local_root=True,
            notes="cancel ordinary prepared packet",
        )) == 0
        context.setattr(
            stage_i,
            "validate_submitted_cancellation_metadata",
            lambda *_args: (_ for _ in ()).throw(
                AssertionError("prepared cancellation treated as submitted cancellation")
            ),
        )
        report = stage_i.reconcile_report(root)
    assert report["consistent"], report["issues"]
    assert report["counts"]["active_reservations"] == 0

    root = tmp_path / "canonical-erasure"
    paths, manifest_path, args = fixture(root)
    assert stage_i.cancel_submitted(args) == 0
    manifest = json.loads(manifest_path.read_text())
    manifest.pop("job_id")
    manifest.pop("cancellation")
    stage_i.write_json(manifest_path, manifest)
    reservations = json.loads(paths["reservations"].read_text())
    reservations[0].pop("job_id")
    stage_i.write_json(paths["reservations"], reservations)
    with monkeypatch.context() as context:
        context.setattr(stage_i, "DEFAULT_ROOT", root.resolve())
        context.setitem(
            stage_i.CANONICAL_SOURCE_BUNDLE_RECOVERY,
            "manifest_relative",
            str(manifest_path.relative_to(root)),
        )
        report = stage_i.reconcile_report(root)
    assert not report["consistent"]
    assert any(
        "submitted cancellation evidence drift" in issue
        for issue in report["issues"]
    )


def test_cgl_lf_stage_i_reconcile_holds_root_lock(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_locked_reconcile_test", PAPER_STAGE_I_TOOL
    )
    assert spec is not None and spec.loader is not None
    stage_i = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage_i
    spec.loader.exec_module(stage_i)

    events = []

    class Lock:
        def __enter__(self):
            events.append("lock-enter")

        def __exit__(self, *_args):
            events.append("lock-exit")

    monkeypatch.setattr(
        stage_i, "canonical_root_lock", lambda _root: Lock()
    )
    monkeypatch.setattr(
        stage_i,
        "reconcile_report",
        lambda _root: events.append("report") or {"consistent": True},
    )
    assert stage_i.reconcile(
        SimpleNamespace(root=str(tmp_path), allow_local_root=True)
    ) == 0
    assert events == ["lock-enter", "report", "lock-exit"]


def test_cgl_lf_stage_i_panel_schema_pins_reference_inventory_and_admission(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_panel_schema_test", PAPER_ANALYZER_PATH
    )
    assert spec is not None and spec.loader is not None
    analyzer = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = analyzer
    spec.loader.exec_module(analyzer)

    configuration = analyzer.stage_i_panels_configuration(PAPER_STAGE_I_MANIFEST)
    panels = configuration["panels"]
    bindings = configuration["reference_product_bindings"]
    assert len(configuration["reference_manifests"]) == 10
    assert len(bindings) == 58
    assert len(panels) == 23
    assert sum(panel["disposition"] == "comparison" for panel in panels) == 11
    assert sum(panel["disposition"] == "blocked_reference" for panel in panels) == 11
    assert sum(panel["disposition"] == "external_model" for panel in panels) == 1
    assert all(
        panel["criterion_state"] == "pending_review"
        for panel in panels if panel["disposition"] == "comparison"
    )
    assert configuration["analysis_case_aliases"][
        "paper_nulim_beta100_hardwall"
    ] == "paper_standard_active_alfvenic_beta100"

    product_id = "fig2a_parallel_athenak_pressure_units"
    binding = bindings[product_id]
    assert analyzer.validate_stage_i_reference_binding(
        bindings, product_id, binding["kind"], binding["case"],
        binding["product"], binding["data_file"], binding["data_sha256"],
        binding["reference_manifest_sha256"],
    )
    with pytest.raises(ValueError, match="Stage I reference binding mismatch"):
        analyzer.validate_stage_i_reference_binding(
            bindings, product_id, "curve", binding["case"],
            binding["product"], binding["data_file"], binding["data_sha256"],
            binding["reference_manifest_sha256"],
        )

    monkeypatch.setattr(analyzer, "bundle_cases", lambda *_: [{
        "case_id": "R02",
        "name": "paper_standard_active_alfvenic_beta10",
    }])
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    with pytest.raises(ValueError, match="accepted_for_analysis"):
        analyzer.stage_i_bundle_case_ids(bundle, {
            "workflow": analyzer.STAGE_I_PRODUCTION_WORKFLOW,
            "status": "pending",
        }, configuration)
    assert analyzer.stage_i_bundle_case_ids(bundle, {
        "workflow": analyzer.STAGE_I_PRODUCTION_WORKFLOW,
        "status": "accepted_for_analysis",
    }, configuration) == ["R02"]


def test_rank_local_binary_reader_keeps_unequal_rank_files(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "bin_convert_rank_output_test", "../../../vis/python/bin_convert.py"
    )
    assert spec is not None and spec.loader is not None
    bin_convert = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = bin_convert
    spec.loader.exec_module(bin_convert)

    rank0 = tmp_path / "rank_00000000" / "case.00000.bin"
    rank1 = tmp_path / "rank_00000001" / "case.00000.bin"
    rank0.parent.mkdir()
    rank1.parent.mkdir()
    rank0.write_bytes(b"short")
    rank1.write_bytes(b"deliberately-longer-rank-output")

    def fake_read(path):
        value = 0.0 if Path(path).parent.name == "rank_00000000" else 1.0
        return {
            "header": [],
            "time": 0.0,
            "cycle": 0,
            "var_names": ["dens"],
            "Nx1": 2, "Nx2": 1, "Nx3": 1, "nvars": 1,
            "x1min": 0.0, "x1max": 1.0,
            "x2min": 0.0, "x2max": 1.0,
            "x3min": 0.0, "x3max": 1.0,
            "n_mbs": 1,
            "nx1_mb": 1, "nx2_mb": 1, "nx3_mb": 1,
            "nx1_out_mb": 1, "nx2_out_mb": 1, "nx3_out_mb": 1,
            "mb_index": np.array([[0, 0, 0, 0, 0, 0]]),
            "mb_logical": np.array([[int(value), 0, 0, 0]]),
            "mb_geometry": np.zeros((1, 6)),
            "mb_data": {"dens": [np.asarray([[[value]]])]},
        }

    monkeypatch.setattr(bin_convert, "read_binary", fake_read)
    combined = bin_convert.read_all_ranks_binary(str(rank0))
    assert combined["n_mbs"] == 2
    assert [values[0, 0, 0] for values in combined["mb_data"]["dens"]] == [0.0, 1.0]


def test_cgl_lf_paper_snapshot_window_skips_unselected_field_reads(
    tmp_path, monkeypatch
):
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_paper_snapshot_selection_test", PAPER_ANALYZER_PATH
    )
    assert spec is not None and spec.loader is not None
    analyzer = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = analyzer
    spec.loader.exec_module(analyzer)

    def write_header(path, time):
        path.write_bytes(
            (
                "Athena binary output version=1.1\n"
                "  size of preheader=5\n"
                f"  time={time}\n"
                "  cycle=0\n"
                "  size of location=8\n"
                "  size of variable=4\n"
            ).encode()
        )

    excluded = tmp_path / "excluded.bin"
    selected = tmp_path / "selected.bin"
    write_header(excluded, 0.0)
    write_header(selected, 1.0)
    read_paths = []

    def fake_read(path):
        read_paths.append(path)
        return {}, (1.0, 1.0, 1.0), 1.0

    monkeypatch.setattr(analyzer, "read_snapshot", fake_read)
    monkeypatch.setattr(analyzer, "pdf_fields", lambda *_: {})
    monkeypatch.setattr(analyzer, "pressure_density_fields", lambda *_: {})
    monkeypatch.setattr(analyzer, "analyze_fields", lambda *_, **__: {})
    monkeypatch.setattr(
        analyzer, "average_snapshot_records",
        lambda records: {"snapshot_count": len(records)},
    )
    records, ensemble = analyzer.analyze_snapshot_paths(
        [excluded, selected], bins=4, alignment_shells=[],
        time_start=1.0, time_end=1.0,
    )
    assert list(records) == [str(selected)]
    assert read_paths == [selected, selected]
    assert ensemble["snapshot_count"] == 1


def test_cgl_lf_paper_snapshot_analysis_uses_both_pressures():
    try:
        _run_paper(
            "cgl_ci_paper_analysis",
            "time/nlim=1",
            "output2/file_type=bin",
            "output2/dt=0.01",
            "output2/single_file_per_rank=true",
        )
        snapshots = sorted(
            Path("bin/rank_00000000").glob(
                "cgl_ci_paper_analysis.mhd_w_bcc.*.bin"
            )
        )
        assert snapshots
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--history",
                "cgl_ci_paper_analysis.user.hst",
                "--time-start",
                "0.0",
                "--time-end",
                "1.0",
                "--synthetic-test",
                "--alignment-shells",
                "2,4,6",
                "--eddy-samples",
                "200000",
                "--eddy-bins",
                "10",
                "--eddy-seed",
                "731",
                "--output-dir",
                "cgl_ci_paper_analysis_products",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        diagnostics = json.loads(
            Path("cgl_ci_paper_analysis_products/diagnostics.json").read_text()
        )
        assert diagnostics["synthetic_test"]["passed"]
        assert (
            diagnostics["synthetic_test"]["positive_perpendicular_heat_flux_proxy"]
            > 0.0
        )
        assert diagnostics["synthetic_test"]["zero_parallel_heat_flux_proxy"] == 0.0
        assert diagnostics["synthetic_test"]["zero_anisotropic_pressure_work"] == 0.0
        assert (
            diagnostics["synthetic_test"]["negative_correlated_anisotropic_pressure_work"]
            < 0.0
        )
        assert abs(
            diagnostics["synthetic_test"]["correlated_transfer_work_difference"]
        ) < 1.0e-14
        assert diagnostics["synthetic_test"]["finite_normalized_transfer"]
        assert (
            diagnostics["synthetic_test"]["positive_correlated_transfer_normalization"]
            > 0.0
        )
        assert diagnostics["synthetic_test"]["passive_pressure_work_is_diagnostic_only"]
        assert diagnostics["synthetic_test"]["finite_eddy_anisotropy"]
        assert diagnostics["synthetic_test"]["finite_pressure_density_joint_pdf"]
        assert (
            diagnostics["synthetic_test"][
                "pressure_density_five_thirds_coordinate_error"
            ] < 1.0e-14
        )
        assert (
            diagnostics["synthetic_test"]["positive_longitudinal_compressive_spectrum"]
            > 0.0
        )
        assert (
            diagnostics["synthetic_test"]["zero_transverse_compressive_spectrum"]
            < 1.0e-28
        )
        assert abs(
            diagnostics["synthetic_test"]["constant_power_quadrature_error"]
        ) < 1.0e-14
        snapshot = next(iter(diagnostics["snapshots"].values()))
        assert "beta_delta" in snapshot["pdf"]
        assert "pressure_density_joint" in snapshot
        assert "parallel" in snapshot["pressure_density_joint"]
        assert "perpendicular" in snapshot["pressure_density_joint"]
        assert "delta_p" in snapshot["spectra"]
        for product in (
            "compressive_velocity",
            "density_fluctuation",
            "p_parallel",
            "p_perp",
            "magnetic_pressure",
        ):
            assert product in snapshot["spectra"]
            assert "field_definition" in snapshot["spectra"][product]
        assert abs(snapshot["pressure_transfer"]["closure_error"]) < 1.0e-12
        assert snapshot["pressure_transfer"]["normalization_available"]
        assert snapshot["pressure_work_decomposition"]["available"]
        assert np.isfinite(
            snapshot["pressure_work_decomposition"]["anisotropic_stress_power"]
        )
        assert diagnostics["snapshot_ensemble"]["snapshot_count"] == 1
        assert "beta_delta" in diagnostics["snapshot_ensemble"]["pdf"]
        assert "pressure_density_joint" in diagnostics["snapshot_ensemble"]
        assert "compressive_velocity" in diagnostics["snapshot_ensemble"]["spectra"]
        assert "magnetic_pressure" in diagnostics["snapshot_ensemble"]["spectra"]
        assert "field_definition" in diagnostics["snapshot_ensemble"]["spectra"][
            "compressive_velocity"
        ]
        assert diagnostics["snapshot_ensemble"]["pressure_transfer"][
            "normalization_available"
        ]
        assert diagnostics["snapshot_ensemble"]["pressure_work_decomposition"][
            "available"
        ]
        assert not diagnostics["snapshot_ensemble"]["pressure_work_decomposition"][
            "time_integral_estimate"
        ]["available"]
        assert diagnostics["histories"][0]["analysis_window"]["rows_selected"] > 0
        assert "eddy_anisotropy" in diagnostics["snapshot_ensemble"]
        module_spec = importlib.util.spec_from_file_location(
            "cgl_lf_paper_analysis", "../../../scripts/analyze_cgl_lf_paper.py"
        )
        assert module_spec is not None and module_spec.loader is not None
        module = importlib.util.module_from_spec(module_spec)
        module_spec.loader.exec_module(module)
        eddy_x, eddy_y = module.analyzed_product_curve({
            "eddy_anisotropy": {
                "velocity_perp": {
                    "available": True,
                    "ell_perp_over_lperp": [0.05, 0.1],
                    "ell_parallel_over_lperp": [0.15, 0.25],
                }
            }
        }, "eddy_anisotropy.velocity_perp")
        assert np.allclose(eddy_x, [0.05, 0.1])
        assert np.allclose(eddy_y, [0.15, 0.25])
        history_series = diagnostics["histories"][0]["time_series"]
        assert "unstable_fraction" in history_series
        reference_dir = Path("cgl_ci_paper_reference_inputs")
        reference_dir.mkdir()
        spectrum = diagnostics["snapshot_ensemble"]["spectra"]["delta_p"]
        sampled = [
            (x, y) for x, y in zip(spectrum["k"], spectrum["power_per_dk"])
            if x > 0.0 and y > 0.0
        ][:3]
        assert len(sampled) >= 2
        curve_path = reference_dir / "delta_p.csv"
        curve_path.write_text(
            "x,y,y_uncertainty\n" + "".join(
                f"{x:.17g},{y:.17g},{max(y * 0.01, 1.0e-30):.17g}\n"
                for x, y in sampled
            )
        )
        curve_digest = hashlib.sha256(curve_path.read_bytes()).hexdigest()
        alignment = diagnostics["snapshot_ensemble"]["alignment"]
        alignment_dk = diagnostics["snapshot_ensemble"]["spectra"]["velocity"]["dk"]
        alignment_sampled = []
        for shell, distribution in sorted(
            alignment.items(), key=lambda item: int(item[0])
        ):
            edges = np.asarray(distribution["edges"])
            density = np.asarray(distribution["density"])
            centers = 0.5 * (edges[1:] + edges[:-1])
            alignment_sampled.append(
                (float(shell) * alignment_dk, centers[np.argmax(density)])
            )
        assert len(alignment_sampled) >= 2
        alignment_curve_path = reference_dir / "alignment_peak.csv"
        alignment_curve_path.write_text(
            "x,y,y_uncertainty\n" + "".join(
                f"{x:.17g},{y:.17g},1.0e-12\n"
                for x, y in alignment_sampled
            )
        )
        alignment_curve_digest = hashlib.sha256(
            alignment_curve_path.read_bytes()
        ).hexdigest()
        normalized_transfer = diagnostics["snapshot_ensemble"]["pressure_transfer"]
        transfer_sampled = list(zip(
            normalized_transfer["k_perp"],
            normalized_transfer["transfer_normalized_by_total"],
        ))[1:4]
        assert len(transfer_sampled) >= 2
        transfer_curve_path = reference_dir / "normalized_transfer.csv"
        transfer_curve_path.write_text(
            "x,y,y_uncertainty\n" + "".join(
                f"{x:.17g},{y:.17g},{max(abs(y) * 0.01, 1.0e-12):.17g}\n"
                for x, y in transfer_sampled
            )
        )
        transfer_curve_digest = hashlib.sha256(
            transfer_curve_path.read_bytes()
        ).hexdigest()
        history_curve_path = reference_dir / "unstable_fraction.csv"
        history_sampled = list(
            zip(history_series["time"], history_series["unstable_fraction"])
        )[:2]
        assert len(history_sampled) == 2
        history_curve_path.write_text(
            "x,y,y_uncertainty\n" + "".join(
                f"{x:.17g},{y:.17g},1.0e-12\n" for x, y in history_sampled
            )
        )
        history_curve_digest = hashlib.sha256(
            history_curve_path.read_bytes()
        ).hexdigest()
        joint = diagnostics["snapshot_ensemble"]["pressure_density_joint"]["parallel"]
        joint_x_edges = np.asarray(joint["x_edges"])
        joint_y_edges = np.asarray(joint["y_edges"])
        joint_density = np.asarray(joint["density"])
        joint_x = 0.5 * (joint_x_edges[1:] + joint_x_edges[:-1])
        joint_y = 0.5 * (joint_y_edges[1:] + joint_y_edges[:-1])
        surface_sampled = [
            (joint_x[index], joint_y[index], joint_density[index, index])
            for index in (16, 24, 32)
        ]
        surface_path = reference_dir / "pressure_density_joint.csv"
        surface_path.write_text(
            "x,y,z,z_uncertainty\n" + "".join(
                f"{x:.17g},{y:.17g},{z:.17g},{max(abs(z) * 0.01, 1.0e-12):.17g}\n"
                for x, y, z in surface_sampled
            )
        )
        surface_digest = hashlib.sha256(surface_path.read_bytes()).hexdigest()
        source_figures = [
            reference_dir / "synthetic_a.pdf",
            reference_dir / "synthetic_b.pdf",
        ]
        source_figures[0].write_bytes(b"synthetic reference panel a\n")
        source_figures[1].write_bytes(b"synthetic reference panel b\n")
        manifest_path = reference_dir / "curves.json"
        manifest_path.write_text(json.dumps({
            "schema_version": 1,
            "provenance": {
                "method": "digitized",
                "source_description": "synthetic regression reference",
                "source_figures": [{
                    "source_figure": source_figure.name,
                    "source_figure_sha256": hashlib.sha256(
                        source_figure.read_bytes()
                    ).hexdigest(),
                } for source_figure in source_figures],
                "digitization_tool": "test fixture",
                "uncertainty_description": "one-percent fixture uncertainty",
            },
            "curves": [{
                "id": "delta_p_exact",
                "case": "direct",
                "product": "spectra.delta_p",
                "data_file": curve_path.name,
                "data_sha256": curve_digest,
                "interpolation": "linear",
            }, {
                "id": "alignment_peak_exact",
                "case": "direct",
                "product": "alignment_peak.cos_theta",
                "data_file": alignment_curve_path.name,
                "data_sha256": alignment_curve_digest,
                "interpolation": "linear",
            }, {
                "id": "normalized_transfer_exact",
                "case": "direct",
                "product": "pressure_transfer.transfer_normalized_by_total",
                "data_file": transfer_curve_path.name,
                "data_sha256": transfer_curve_digest,
                "interpolation": "linear",
            }, {
                "id": "unstable_fraction_exact",
                "case": "direct",
                "product": "history.unstable_fraction",
                "data_file": history_curve_path.name,
                "data_sha256": history_curve_digest,
                "interpolation": "linear",
            }],
            "surfaces": [{
                "id": "pressure_density_joint_exact",
                "case": "direct",
                "product": "pressure_density_joint.parallel",
                "data_file": surface_path.name,
                "data_sha256": surface_digest,
                "interpolation": "bilinear",
            }],
        }))
        second_manifest = json.loads(manifest_path.read_text())
        second_manifest["curves"] = [{
            **second_manifest["curves"][0],
            "id": "delta_p_exact_second",
        }]
        second_manifest["surfaces"] = []
        second_manifest_path = reference_dir / "second_curves.json"
        second_manifest_path.write_text(json.dumps(second_manifest))
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--history",
                "cgl_ci_paper_analysis.user.hst",
                "--reference-curves",
                str(manifest_path),
                "--reference-curves",
                str(second_manifest_path),
                "--alignment-shells",
                "2,4,6",
                "--eddy-samples",
                "200000",
                "--eddy-bins",
                "10",
                "--eddy-seed",
                "731",
                "--output-dir",
                "cgl_ci_paper_reference_products",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        reference = json.loads(
            Path("cgl_ci_paper_reference_products/diagnostics.json").read_text()
        )["reference_curve_comparisons"]
        assert reference["available"]
        assert len(reference["manifests"]) == 2
        comparison = reference["comparisons"]["delta_p_exact"]
        assert comparison["sample_count"] == len(sampled)
        assert comparison["maximum_absolute_residual"] < 1.0e-14
        assert comparison["rms_normalized_by_reported_uncertainty"] < 1.0e-12
        assert (
            reference["comparisons"]["delta_p_exact_second"][
                "maximum_absolute_residual"
            ] < 1.0e-14
        )
        alignment_comparison = reference["comparisons"]["alignment_peak_exact"]
        assert alignment_comparison["sample_count"] == len(alignment_sampled)
        assert alignment_comparison["maximum_absolute_residual"] < 1.0e-14
        transfer_comparison = reference["comparisons"]["normalized_transfer_exact"]
        assert transfer_comparison["sample_count"] == len(transfer_sampled)
        assert transfer_comparison["maximum_absolute_residual"] < 1.0e-14
        history_comparison = reference["comparisons"]["unstable_fraction_exact"]
        assert history_comparison["sample_count"] == len(history_sampled)
        assert history_comparison["maximum_absolute_residual"] < 1.0e-14
        surface_comparison = reference["surface_comparisons"][
            "pressure_density_joint_exact"
        ]
        assert surface_comparison["sample_count"] == len(surface_sampled)
        assert surface_comparison["maximum_absolute_residual"] < 1.0e-14
        rendered = subprocess.run(
            [
                sys.executable,
                "../../../scripts/plot_cgl_lf_paper.py",
                "--diagnostics",
                "cgl_ci_paper_reference_products/diagnostics.json",
                "--figure-dir",
                "cgl_ci_paper_reference_figures",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert rendered.returncode == 0, rendered.stdout + rendered.stderr
        assert Path(
            "cgl_ci_paper_reference_figures/"
            "paper_reference_surface_pressure_density_joint_exact.pdf"
        ).is_file()
        invalid_manifest = json.loads(manifest_path.read_text())
        invalid_manifest["provenance"]["source_figures"][1][
            "source_figure_sha256"
        ] = "0" * 64
        invalid_path = reference_dir / "invalid_curves.json"
        invalid_path.write_text(json.dumps(invalid_manifest))
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--reference-curves",
                str(invalid_path),
                "--output-dir",
                "cgl_ci_paper_invalid_reference_products",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode != 0
        assert "source figure checksum does not match" in result.stderr
        duplicate_manifest = json.loads(second_manifest_path.read_text())
        duplicate_path = reference_dir / "duplicate_curves.json"
        duplicate_path.write_text(json.dumps(duplicate_manifest))
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--reference-curves",
                str(second_manifest_path),
                "--reference-curves",
                str(duplicate_path),
                "--output-dir",
                "cgl_ci_paper_duplicate_reference_products",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode != 0
        assert "duplicated across manifests" in result.stderr
        partial_manifest = json.loads(manifest_path.read_text())
        partial_manifest["curves"].append({
            **partial_manifest["curves"][0],
            "id": "absent_case_reference",
            "case": "case_not_in_bundle",
        })
        partial_path = reference_dir / "partial_curves.json"
        partial_path.write_text(json.dumps(partial_manifest))
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--history",
                "cgl_ci_paper_analysis.user.hst",
                "--reference-curves",
                str(partial_path),
                "--alignment-shells",
                "2,4,6",
                "--output-dir",
                "cgl_ci_paper_unscoped_reference_products",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode != 0
        assert "selects missing case case_not_in_bundle" in result.stderr
        result = subprocess.run(
            [
                sys.executable,
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--history",
                "cgl_ci_paper_analysis.user.hst",
                "--reference-curves",
                str(partial_path),
                "--allow-partial-reference-cases",
                "--alignment-shells",
                "2,4,6",
                "--output-dir",
                "cgl_ci_paper_partial_reference_products",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        partial = json.loads(
            Path("cgl_ci_paper_partial_reference_products/diagnostics.json").read_text()
        )["reference_curve_comparisons"]
        assert partial["allow_missing_cases"]
        assert partial["comparisons"]["delta_p_exact"]["available"]
        assert partial["omitted_products"][0]["id"] == "absent_case_reference"
        assert partial["omitted_products"][0]["case"] == "case_not_in_bundle"
        assert "absent_case_reference" not in partial["comparisons"]
    finally:
        shutil.rmtree("bin", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_analysis_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_reference_inputs", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_reference_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_reference_figures", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_invalid_reference_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_duplicate_reference_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_unscoped_reference_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_partial_reference_products", ignore_errors=True)
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_explicit_reference_rejects_sts_configuration():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_decay.athinput",
        "mhd/cgl_heat_flux_integrator=explicit",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "explicit reference integration requires" in result.stdout


def test_cgl_lf_amr_rejects_primitive_prolongation():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_amr_2d.athinput",
        "mesh_refinement/prolong_primitives=true",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "use conserved prolongation for LF AMR runs" in result.stdout
