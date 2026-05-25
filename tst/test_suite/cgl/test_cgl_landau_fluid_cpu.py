"""CPU regressions for the CGL Landau-fluid closure and CGL FOFC path."""

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"
PAPER_INPUT = "../../../inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"
PAPER_PASSIVE_INPUT = (
    "../../../inputs/cgl_lf_paper/cgl_lf_paper_smoke_passive_beta10.athinput"
)
PAPER_PRODUCTION_INPUT_ROOT = Path("../../../inputs/cgl_lf_paper")
PAPER_WORKFLOW_PATH = Path("../../../scripts/cgl_lf_workflow.py")


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
        _run_paper("cgl_ci_paper_restart_partial", "time/nlim=1")
        restart_paths = sorted(Path("rst").glob("cgl_ci_paper_restart_partial.*.rst"))
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
    assert len(input_paths) == 8
    hardwall_paths = {
        "cgl_lf_paper_standard_active_alfvenic_beta1.athinput",
        "cgl_lf_paper_standard_active_alfvenic_beta10.athinput",
        "cgl_lf_paper_standard_active_alfvenic_beta100.athinput",
        "cgl_lf_paper_standard_active_random_beta10.athinput",
        "cgl_lf_paper_standard_passive_alfvenic_beta10.athinput",
        "cgl_lf_paper_nulim_beta100_hardwall.athinput",
    }
    for input_path in input_paths:
        source = input_path.read_text()
        for block in ("output2", "output3"):
            body = source.split(f"<{block}>", 1)[1].split("<", 1)[0]
            assert "single_file_per_rank = false" in body
        choices = workflow.model_choices(source, [])
        assert choices["time_integrator"] == "rk2"
        assert choices["time_sts_integrator"] == "rkl2"
        assert choices["time_sts_max_dt_ratio"] == "-1.0"
        assert choices["time_cfl_number"] == "0.3"
        assert choices["time_tlim"] == "10.0"
        assert choices["output2_file_type"] == "bin"
        assert choices["output2_dt"] == "0.25"
        assert choices["output2_single_file_per_rank"] == "false"
        assert choices["output3_file_type"] == "rst"
        assert choices["output3_dt"] == "1.0"
        assert choices["output3_single_file_per_rank"] == "false"
        expected_hardwall = "true" if input_path.name in hardwall_paths else "false"
        assert choices["limiter_hardwall"] == expected_hardwall


def test_cgl_lf_paper_snapshot_analysis_uses_both_pressures():
    try:
        _run_paper(
            "cgl_ci_paper_analysis",
            "time/nlim=1",
            "output2/file_type=bin",
            "output2/dt=0.01",
        )
        snapshots = sorted(Path("bin").glob("cgl_ci_paper_analysis.mhd_w_bcc.*.bin"))
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
        assert diagnostics["synthetic_test"]["passive_pressure_work_is_diagnostic_only"]
        assert abs(
            diagnostics["synthetic_test"]["constant_power_quadrature_error"]
        ) < 1.0e-14
        snapshot = next(iter(diagnostics["snapshots"].values()))
        assert "beta_delta" in snapshot["pdf"]
        assert "delta_p" in snapshot["spectra"]
        assert abs(snapshot["pressure_transfer"]["closure_error"]) < 1.0e-12
        assert snapshot["pressure_work_decomposition"]["available"]
        assert np.isfinite(
            snapshot["pressure_work_decomposition"]["anisotropic_stress_power"]
        )
        assert diagnostics["snapshot_ensemble"]["snapshot_count"] == 1
        assert "beta_delta" in diagnostics["snapshot_ensemble"]["pdf"]
        assert diagnostics["snapshot_ensemble"]["pressure_work_decomposition"][
            "available"
        ]
        assert not diagnostics["snapshot_ensemble"]["pressure_work_decomposition"][
            "time_integral_estimate"
        ]["available"]
        assert diagnostics["histories"][0]["analysis_window"]["rows_selected"] > 0
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
        source_figure = reference_dir / "synthetic.pdf"
        source_figure.write_bytes(b"synthetic reference panel\n")
        manifest_path = reference_dir / "curves.json"
        manifest_path.write_text(json.dumps({
            "schema_version": 1,
            "provenance": {
                "method": "digitized",
                "source_description": "synthetic regression reference",
                "source_figure": source_figure.name,
                "source_figure_sha256": hashlib.sha256(
                    source_figure.read_bytes()
                ).hexdigest(),
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
                "id": "unstable_fraction_exact",
                "case": "direct",
                "product": "history.unstable_fraction",
                "data_file": history_curve_path.name,
                "data_sha256": history_curve_digest,
                "interpolation": "linear",
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
        history_comparison = reference["comparisons"]["unstable_fraction_exact"]
        assert history_comparison["sample_count"] == len(history_sampled)
        assert history_comparison["maximum_absolute_residual"] < 1.0e-14
        invalid_manifest = json.loads(manifest_path.read_text())
        invalid_manifest["provenance"]["source_figure_sha256"] = "0" * 64
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
