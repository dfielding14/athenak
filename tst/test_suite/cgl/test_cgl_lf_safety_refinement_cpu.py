"""Full-model temporal refinement through a finite-rate LF soft threshold.

This checks an unforced smooth flow, not turbulent switching or split second order.
"""

import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np

import test_suite.testutils as testutils


REPO_ROOT = Path(__file__).resolve().parents[3]
END_TIME = 0.2
OUTER_BUDGETS = (4.8, 2.4, 1.2, 0.6)
SAFETIES = (0.3, 0.9)


def _threshold_run(directory, safety, outer_budget):
    directory.mkdir()
    source = (REPO_ROOT / "inputs/unit_tests/cgl_lf_limiter_mirror.athinput").read_text()
    for old, new in (
        ("evolution = kinematic", "evolution = dynamic"),
        ("rsolver = advect", "rsolver = hlle"),
        ("nlim = 1", "nlim = -1"),
        ("tlim = 0.005", f"tlim = {END_TIME}"),
        ("sts_max_dt_ratio = 1000.0", f"sts_max_dt_ratio = {outer_budget/safety:.17g}"),
        ("limiter_nu_coll = 100.0", "limiter_nu_coll = 10.0"),
        ("lf_coefficient_mode = background", "lf_coefficient_mode = local"),
    ):
        assert source.count(old) == 1, old
        source = source.replace(old, new)
    # Background isotropization carries faces across the soft threshold; soft
    # relaxation alone approaches its threshold asymptotically in a uniform state.
    source += f"""
<time>
sts_safety = {safety:.17g}
<mhd>
nu_coll = 2.0
mirror_threshold = 1.0
mirror_backup_factor = 2.0
cgl_lf_arithmetic = safe
cgl_lf_sts_flux = weighted
cgl_lf_diagnostics = full
cgl_lf_profile = false
cgl_lf_profile_detail = false
cgl_lf_record_pressure_work = true
<output1>
file_type = hst
data_format = %24.16e
dcycle = 1
<output2>
file_type = tab
variable = mhd_w_bcc
data_format = %24.16e
dcycle = 1
"""
    deck = directory / "crossing.athinput"
    deck.write_text(source)
    environment = os.environ.copy()
    for key in ("PROFILE", "PROFILE_DETAIL", "DIAGNOSTICS", "ARITHMETIC", "STS_FLUX"):
        environment.pop(f"ATHENAK_CGL_LF_{key}", None)
    command = [str(Path("athena").resolve()), "-i", str(deck), "-d", str(directory)]
    result = subprocess.run(command, capture_output=True, text=True, env=environment)
    (directory / "run.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CGL LF limiter_stress passed" in result.stdout

    # Finalize writes a second terminal snapshot. Keep one state per cycle.
    snapshots = {}
    for path in sorted((directory / "tab").glob("*.mhd_w_bcc.*.tab")):
        # This 1D writer labels unused j/k coordinates but omits their data.
        with path.open() as stream:
            header, fields = stream.readline(), stream.readline().split()[-9:]
        assert fields == ["dens", "velx", "vely", "velz", "eint", "p_perp",
                          "bcc1", "bcc2", "bcc3"]
        data = np.loadtxt(path)
        assert data.shape == (64, 12) and np.all(np.isfinite(data))
        snapshot = dict(zip(fields, data[:, 3:].T))
        snapshot["time"] = float(re.search(r"time=(\S+)", header)[1])
        snapshot["cycle"] = int(re.search(r"cycle=(\d+)", header)[1])
        snapshot["x1v"] = data[:, 2]
        if snapshot["cycle"] in snapshots:
            previous = snapshots[snapshot["cycle"]]
            for field in snapshot.keys() - {"time", "cycle"}:
                np.testing.assert_array_equal(snapshot[field], previous[field])
        snapshots[snapshot["cycle"]] = snapshot
    assert snapshots
    cycles = sorted(snapshots)
    assert cycles == list(range(cycles[-1] + 1))
    initial, final = snapshots[0], snapshots[cycles[-1]]
    assert len(initial["eint"]) == 64
    assert initial["time"] == 0 and final["time"] == END_TIME
    face_distances = []
    for cycle in cycles:
        snapshot = snapshots[cycle]
        for field in ("dens", "eint", "p_perp"):
            assert np.all(np.isfinite(snapshot[field])) and np.min(snapshot[field]) > 0
        magnitude = np.sqrt(sum(snapshot[f"bcc{i}"]**2 for i in (1, 2, 3)))
        face_b = 0.5*magnitude + 0.5*np.roll(magnitude, 1)
        anisotropy = snapshot["p_perp"] - snapshot["eint"]
        face_anisotropy = 0.5*anisotropy + 0.5*np.roll(anisotropy, 1)
        face_distances.append(face_anisotropy - 0.5*face_b**2)
    distances = np.asarray(face_distances)
    crossed = np.flatnonzero((distances[0] > 1e-4) & (np.min(distances[1:], axis=0) < -1e-4))
    assert crossed.size, "no initially active face crossed the enabled mirror threshold"
    face = int(crossed[0])
    crossing = int(np.flatnonzero(distances[:, face] < -1e-4)[0])

    history_path, = directory.glob("*.mhd.hst")
    history = testutils.athena_read.hst(str(history_path))
    steps = np.diff(history["time"])
    assert len(steps) == cycles[-1] and np.all(steps > 0)
    assert history["time"][0] == 0 and history["time"][-1] == END_TIME
    events = {key: float(history[key][-1]) for key in (
        "lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd", "lf_hwproj")}
    assert all(value == 0 for value in events.values()), events
    energy_error = float(np.max(abs(history["tot-E"]/history["tot-E"][0] - 1)))
    assert energy_error < 5e-12
    q_work = float(np.ptp(history["lf_qprwrk"]) + np.ptp(history["lf_qpewrk"]))
    assert q_work > 1e-8, "threshold fixture has no resolved LF transport"
    rhs = np.diff(history["lf_nstage"])/64
    assert np.all(rhs >= 6) and np.all(rhs == np.rint(rhs))
    report = {
        "safety": safety, "outer_budget": outer_budget,
        "sts_max_dt_ratio": outer_budget/safety, "cycles": cycles[-1],
        "dt_min": float(np.min(steps)), "dt_max": float(np.max(steps)),
        "rhs_total": int(np.sum(rhs)), "rhs_per_cycle_min": float(np.min(rhs)),
        "rhs_per_cycle_max": float(np.max(rhs)), "events": events,
        "relative_total_energy_drift": energy_error, "q_work_variation": q_work,
        "crossed_faces": int(crossed.size), "witness_face": face,
        "witness_initial_threshold_distance": float(distances[0, face]),
        "witness_crossing_threshold_distance": float(distances[crossing, face]),
        "witness_crossing_cycle": cycles[crossing],
        "witness_crossing_time": float(history["time"][crossing]),
        "initial_pressure_means": [float(np.mean(initial[key])) for key in ("eint", "p_perp")],
        "final_pressure_means": [float(np.mean(final[key])) for key in ("eint", "p_perp")],
        "final_pressure_minima": [float(np.min(final[key])) for key in ("eint", "p_perp")],
        "final_thermal_energy_mean": float(np.mean(0.5*final["eint"] + final["p_perp"])),
    }
    (directory / "measurements.json").write_text(json.dumps(report, indent=2) + "\n")
    return final, report


def test_cgl_lf_sts_safety_threshold_crossing_refines(tmp_path):
    """Compare both pressure fields at one endpoint; do not assume second order."""
    runs = {}
    for safety in SAFETIES:
        for level, outer_budget in enumerate(OUTER_BUDGETS):
            runs[safety, level] = _threshold_run(
                tmp_path / f"safety{int(10*safety)}_level{level}", safety, outer_budget)
    runs[0.3, 4] = _threshold_run(tmp_path / "safety3_level4", 0.3, 0.3)
    reference = runs[0.3, 4][0]
    thermal_scale = float(np.mean(0.5*reference["eint"] + reference["p_perp"]))
    errors = {}
    for (safety, level), (state, report) in runs.items():
        report["pressure_errors"] = {}
        for field in ("eint", "p_perp"):
            difference = state[field] - reference[field]
            error = {
                "l1": float(np.mean(abs(difference))),
                "l2": float(np.sqrt(np.mean(difference**2))),
                "linf": float(np.max(abs(difference))),
                "max_relative": float(np.max(abs(difference)/reference[field])),
                "thermal_scaled_linf": float(np.max(abs(difference))/thermal_scale),
                "mean_difference": float(np.mean(difference)),
            }
            report["pressure_errors"][field] = error
            errors[safety, level, field] = error
    report = {
        "scope": "unforced full-model finite-rate mirror crossing; no second-order claim",
        "reference": "safety3_level4", "endpoint": END_TIME,
        "grid": [64, 1, 1], "background_nu": 2.0, "soft_nu": 10.0,
        "reference_refinement_difference": errors[0.3, 3, "eint"],
        "reference_refinement_difference_perp": errors[0.3, 3, "p_perp"],
        "runs": [report for _, report in runs.values()],
    }
    (tmp_path / "threshold-refinement.json").write_text(json.dumps(report, indent=2) + "\n")
    for safety in SAFETIES:
        for level in range(3):
            ratio = runs[safety, level+1][1]["dt_max"]/runs[safety, level][1]["dt_max"]
            assert 0.4 < ratio < 0.6, "outer-step ladder did not halve actual steps"
        for field in ("eint", "p_perp"):
            for norm in ("l1", "l2", "linf"):
                ladder = [errors[safety, level, field][norm] for level in range(3)]
                assert 0 < ladder[2] < ladder[1] < ladder[0], (safety, field, norm, ladder)
                uncertainty = errors[0.3, 3, field][norm]
                assert uncertainty < 0.5*ladder[2], (safety, field, norm, uncertainty, ladder)
    # The safety comparison must actually select different RKL work at the top step.
    assert runs[0.9, 0][1]["rhs_per_cycle_max"] < runs[0.3, 0][1]["rhs_per_cycle_max"]
