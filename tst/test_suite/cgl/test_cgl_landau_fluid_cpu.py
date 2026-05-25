"""CPU regressions for the CGL Landau-fluid closure and CGL FOFC path."""

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"
PAPER_INPUT = "../../../inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"
PAPER_PASSIVE_INPUT = (
    "../../../inputs/cgl_lf_paper/cgl_lf_paper_smoke_passive_beta10.athinput"
)
PAPER_PRODUCTION_INPUT_ROOT = Path("../../../inputs/cgl_lf_paper")
PAPER_STAGE_I_MANIFEST = PAPER_PRODUCTION_INPUT_ROOT / "mks24_stage_i_manifest.json"
PAPER_WORKFLOW_PATH = Path("../../../scripts/cgl_lf_workflow.py")
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
                "python3",
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
        dt = mhd["time"][-1] - mhd["time"][0]
        expected_work = 0.5 * user["force_prp2"][-1] * dt * dt
        measured_work = mhd["tot-E"][-1] - mhd["tot-E"][0]
        assert np.isclose(measured_work, expected_work, rtol=1.0e-10, atol=1.0e-12)
        applied_work = user["force_work"][-1] - user["force_work"][0]
        assert np.isclose(applied_work, measured_work, rtol=1.0e-10, atol=1.0e-12)
    finally:
        shutil.rmtree("rst", ignore_errors=True)
        _cleanup()


def test_cgl_lf_paper_forcing_ou_state_advances_once_per_cycle():
    try:
        _run_paper("cgl_ci_ou_correlated", "time/nlim=1")
        _run_paper("cgl_ci_ou_white", "time/nlim=1", "turb_driving/tcorr=0.0")
        correlated = _final_variable_tab("cgl_ci_ou_correlated", "turb_force")
        white = _final_variable_tab("cgl_ci_ou_white", "turb_force")
        history = testutils.athena_read.hst("cgl_ci_ou_correlated.mhd.hst")
        dt = history["time"][-1] - history["time"][0]
        fcorr = np.exp(-dt / 2.0)
        gcorr = np.sqrt(1.0 - fcorr * fcorr)
        for field in ("force1", "force2", "force3"):
            assert np.allclose(
                correlated[field], gcorr * white[field], rtol=5.0e-12, atol=5.0e-14
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


def test_cgl_lf_paper_production_inputs_explicitly_use_shared_mpiio():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_workflow_test", PAPER_WORKFLOW_PATH
    )
    assert spec is not None and spec.loader is not None
    workflow = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = workflow
    spec.loader.exec_module(workflow)
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
                "python3",
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
    output_dir = root / "runs" / "mks24-stage-i" / "R16" / "s00" / "output"
    manifest_dir = output_dir.parent / "manifest"
    (output_dir / "bin").mkdir(parents=True)
    (output_dir / "rst").mkdir()
    manifest_dir.mkdir()
    mhd_history = output_dir / "case.mhd.hst"
    (output_dir / "case.user.hst").write_text("# retained user history\n")
    (output_dir / "bin" / "case.00000.bin").write_bytes(b"snapshot")
    (output_dir / "rst" / "case.00000.rst").write_bytes(b"restart")

    def write_history(rows):
        mhd_history.write_text(
            "# [0]=time [1]=lf_dfloor [2]=lf_pfloor [3]=lf_nonfin "
            "[4]=lf_nonpos [5]=lf_hardbd [6]=lf_hwproj\n"
            + "\n".join(" ".join(str(value) for value in row) for row in rows)
            + "\n"
        )

    manifest_path = manifest_dir / "prepared_run.json"
    manifest = {
        "state": "submitted",
        "project_root": str(root),
        "job_id": "test-job",
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
        },
        "paths": {"output_dir": str(output_dir)},
    }
    stage_i.write_json(manifest_path, manifest)
    paths = stage_i.initialize(root)
    stage_i.write_json(paths["reservations"], [{
        "manifest": str(manifest_path),
        "case_id": "R16",
        "case_name": "case",
        "segment": "s00",
        "nodes": 1,
        "requested_walltime": "00:10:00",
        "reserved_node_hours": 1.0 / 6.0,
        "state": "submitted",
    }])
    inspect_args = SimpleNamespace(
        manifest=str(manifest_path), required_time=2.0
    )
    write_history([(0.0, 0, 0, 0, 0, 0, 0), (1.0, 0, 0, 0, 0, 0, 1)])
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

    write_history([(0.0, 0, 0, 0, 0, 0, 0), (2.0, 0, 0, 0, 0, 0, 2)])
    assert stage_i.inspect_segment(inspect_args) == 0
    inspection_path.unlink()
    sacct_path = tmp_path / "job.sacct"
    sacct_path.write_text(
        "test-job|case|COMPLETED|0:0|1|60|submit-time|end-time|\n"
    )
    record_args = SimpleNamespace(
        manifest=str(manifest_path),
        allow_local_root=True,
        job_id="test-job",
        result="accepted",
        notes="test",
        sacct_file=str(sacct_path),
    )
    with pytest.raises(ValueError, match="inspect-segment evidence"):
        stage_i.record(record_args)
    assert stage_i.inspect_segment(inspect_args) == 0
    assert stage_i.record(record_args) == 0
    accounted = json.loads(manifest_path.read_text())
    assert accounted["state"] == "recorded"
    assert accounted["scientific_inspection"]["accepted"]
    parent = stage_i.verify_continuation_restart(
        output_dir / "rst" / "case.00000.rst"
    )
    assert parent["result"] == "accepted"
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
                "python3",
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
        result = subprocess.run(
            [
                "python3",
                "../../../scripts/analyze_cgl_lf_paper.py",
                str(snapshots[-1]),
                "--history",
                "cgl_ci_paper_analysis.user.hst",
                "--reference-curves",
                str(manifest_path),
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
        comparison = reference["comparisons"]["delta_p_exact"]
        assert comparison["sample_count"] == len(sampled)
        assert comparison["maximum_absolute_residual"] < 1.0e-14
        assert comparison["rms_normalized_by_reported_uncertainty"] < 1.0e-12
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
                "python3",
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
                "python3",
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
    finally:
        shutil.rmtree("bin", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_analysis_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_reference_inputs", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_reference_products", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_reference_figures", ignore_errors=True)
        shutil.rmtree("cgl_ci_paper_invalid_reference_products", ignore_errors=True)
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
