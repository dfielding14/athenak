"""Transverse LF slope checks with independent flux and diffusion oracles."""

import csv
import json

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl.test_cgl_lf_geometry_cpu import INPUTS, _run, _work


@pytest.mark.parametrize("arithmetic", ["safe", "fast"])
@pytest.mark.parametrize("diagnostics", ["full", "none"])
def test_cgl_lf_transverse_face_flux(tmp_path, arithmetic, diagnostics):
    output = _run(
        tmp_path / "faces", "cgl_lf_geometry.athinput",
        "problem/test_mode=transverse_flux", "mesh/x2max=2", "mesh/x3max=3",
        f"mhd/cgl_lf_arithmetic={arithmetic}",
        f"mhd/cgl_lf_diagnostics={diagnostics}", "mhd/cgl_lf_profile=true",
    )
    assert "transverse_flux passed" in output


@pytest.mark.parametrize("angle", [30, 45])
@pytest.mark.parametrize("arithmetic", ["safe", "fast"])
@pytest.mark.parametrize("diagnostics", ["full", "none"])
def test_cgl_lf_hotspot(tmp_path, angle, arithmetic, diagnostics):
    case = tmp_path / "hotspot"
    theta = np.deg2rad(angle)
    output = _run(
        case, "cgl_lf_hotspot.athinput", f"problem/b0={10*np.cos(theta):.17g}",
        f"problem/by0={10*np.sin(theta):.17g}", f"mhd/cgl_lf_arithmetic={arithmetic}",
        f"mhd/cgl_lf_diagnostics={diagnostics}", "mhd/cgl_lf_profile=true",
    )
    assert "CGL LF hotspot:" in output
    _hotspot_checks(case, output)


def _hotspot_checks(case, output):
    history, _ = _work(case, output, 64**2, False)
    np.testing.assert_allclose(history["time"][-1], 0.201, rtol=0, atol=1e-15)
    monitor = testutils.athena_read.hst(str(next(case.glob("*.user.hst"))))
    for name in ("min_tpar", "min_tperp"):
        assert np.all(np.isfinite(monitor[name]))
        assert np.min(monitor[name]) >= monitor[name][0]*(1-1e-12)
    for name in ("err_energy", "err_moment"):
        assert np.all(np.isfinite(monitor[name]))
        assert np.max(monitor[name]) <= 5e-13
    # Pressure recovery subtracts magnetic energy; use the total-E absolute budget.
    e0 = history["tot-E"][0]
    nonthermal = sum(history[f"{axis}-{kind}"][0]
                     for axis in (1, 2, 3) for kind in ("KE", "ME"))
    thermal_fraction = (e0-nonthermal)/e0
    assert 0 < thermal_fraction <= 1
    thermal = monitor["err_therma"]  # History labels are limited to ten characters.
    assert np.all(np.isfinite(thermal))
    assert np.max(thermal)*thermal_fraction <= 5e-13
    return history, monitor


def _hotspot_state(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    assert data.shape == (64**2,) and np.all(np.isfinite(data.tolist()))
    np.testing.assert_array_equal(data["i"], np.tile(np.arange(64), 64))
    np.testing.assert_array_equal(data["j"], np.repeat(np.arange(64), 64))
    for field in ("dens", "ppar", "pperp"):
        assert np.min(data[field]) > 0
    magnetic = sum(data[f"bcc{axis}"]**2 for axis in (1, 2, 3))
    kinetic = 0.5*data["dens"]*sum(data[f"vel{axis}"]**2 for axis in "xyz")
    np.testing.assert_allclose(data["energy"],
                               0.5*data["ppar"] + data["pperp"] +
                               0.5*magnetic + kinetic, rtol=2e-15, atol=0)
    np.testing.assert_allclose(data["A"], data["dens"]*np.log(
        data["pperp"]*data["dens"]**2/(data["ppar"]*magnetic**1.5)),
                               rtol=2e-15, atol=0)
    return data


def _hotspot_error(coarse, fine, thermal_scale):
    errors = {}
    for field in ("ppar", "pperp"):
        difference = coarse[field] - fine[field]
        errors[field] = {
            "l1": float(np.mean(abs(difference))),
            "l2": float(np.sqrt(np.mean(difference**2))),
            "linf": float(np.max(abs(difference))),
            "max_relative": float(np.max(abs(difference)/fine[field])),
            "thermal_scaled_linf": float(np.max(abs(difference))/thermal_scale),
            "mean_difference": float(np.mean(difference)),
        }
    return errors


def test_cgl_lf_local_hotspot_temporal_refinement(tmp_path):
    """Refine nonlinear LF at fixed geometry; the finest state is not an exact oracle."""
    source = (INPUTS / "cgl_lf_hotspot.athinput").read_text()
    source += ("\n<time>\nsts_safety = 0.3\n"
               "<problem>\nvalidation_output = true\nvalidation_output_dir = validation\n"
               "<output2>\nfile_type = log\ndcycle = 1\n")
    # Use the prescribed truncation window; finer sweeps separately stress roundoff.
    ratios = (8.0, 4.0, 2.0, 1.0)
    states, measurements = [], []
    for level, ratio in enumerate(ratios):
        case = tmp_path / f"level{level}"
        output = _run(
            case, "cgl_lf_hotspot.athinput", f"time/sts_max_dt_ratio={ratio}",
            "mhd/lf_coefficient_mode=local", "mhd/cgl_lf_diagnostics=none",
            source=source,
        )
        assert "CGL LF hotspot:" in output
        history, monitor = _hotspot_checks(case, output)
        initial = _hotspot_state(next((case / "validation").glob("*.hotspot_initial.csv")))
        final = _hotspot_state(next((case / "validation").glob("*.hotspot_final.csv")))
        np.testing.assert_array_equal(initial["time"], 0)
        np.testing.assert_array_equal(final["time"], 0.201)
        steps = np.diff(history["time"])
        np.testing.assert_array_equal(initial["cycle"], 0)
        np.testing.assert_array_equal(final["cycle"], len(steps))
        gaussian = 1 + 99*np.exp(-(((initial["i"]-32)/64)**2 +
                                    ((initial["j"]-32)/64)**2)/(2*0.04**2))
        for field in ("ppar", "pperp"):
            np.testing.assert_allclose(initial[field], gaussian, rtol=3e-15, atol=0)
            assert np.max(abs(final[field]-initial[field])) > 1
            assert final[field].min() >= initial[field].min()*(1-1e-12)
        for field in ("dens", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3"):
            np.testing.assert_array_equal(final[field], initial[field])
        if states:
            np.testing.assert_array_equal(initial, initial_reference)
        else:
            initial_reference = initial
        for field in ("energy", "pperp"):
            np.testing.assert_allclose(final[field].mean(), initial[field].mean(),
                                       rtol=5e-13, atol=0)
        event_log = (case / "cgl_lf_hotspot.log").read_text()
        assert "# Athena event counter data" in event_log
        assert not [line for line in event_log.splitlines()
                    if line.strip() and not line.startswith("#")], event_log
        rhs = np.diff(history["lf_nstage"])/64**2
        if ratio == 8:
            assert np.all(np.isin(rhs, (6, 8, 10))) and rhs.max() == 10, rhs
        else:
            np.testing.assert_array_equal(rhs, 6)
        rhs_values, rhs_counts = np.unique(rhs, return_counts=True)
        states.append(final)
        measurements.append({
            "ratio": ratio, "cycles": len(steps), "rhs": int(np.sum(rhs)),
            "rhs_per_cycle_distribution": {int(value): int(count)
                                           for value, count in zip(rhs_values, rhs_counts)},
            "dt_min": float(steps.min()), "dt_max": float(steps.max()),
            "pressure_minima": [float(final[key].min()) for key in ("ppar", "pperp")],
            "pressure_maxima": [float(final[key].max()) for key in ("ppar", "pperp")],
            "pressure_means": [float(final[key].mean()) for key in ("ppar", "pperp")],
            "monitor_maxima": {key: float(monitor[key].max()) for key in
                               ("err_energy", "err_moment", "err_therma")},
            "minimum_temperatures": [float(monitor[key].min()) for key in
                                     ("min_tpar", "min_tperp")],
            "lf_events": {key: float(history[key][-1]) for key in
                          ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos",
                           "lf_hardbd", "lf_hwproj")},
            "eos_event_rows": 0,
        })
    thermal_scale = float(np.mean(0.5*initial["ppar"]+initial["pperp"]))
    adjacent = [_hotspot_error(coarse, fine, thermal_scale)
                for coarse, fine in zip(states[:-1], states[1:])]
    errors = [_hotspot_error(state, states[-1], thermal_scale) for state in states[:-1]]
    report = {
        "scope": "nonlinear local LF, frozen uniform field, zero velocity, no collisions",
        "grid": [64, 64, 1], "endpoint": 0.201, "safety": 0.3,
        "nominal_diffusion_times": 0.201*np.sqrt(2/np.pi)/10/0.04**2,
        "thermal_energy_density_scale": thermal_scale,
        "runs": measurements, "adjacent_refinement_differences": adjacent,
        "errors_against_ratio_1": errors,
        "fixed_stage_4_2_1_log2_refinement_indicators": {
            field: {norm: (float(np.log2(adjacent[1][field][norm]/adjacent[2][field][norm]))
                           if adjacent[1][field][norm] > 0 and adjacent[2][field][norm] > 0
                           else None) for norm in ("l1", "l2", "linf")}
            for field in ("ppar", "pperp")},
        "stage_transition_limit": "Ratio 8 may use five stages per half sweep, while "
                                  "4, 2, 1 use three; their leading error constants differ.",
        "reference_limit": "The last adjacent difference is a refinement indicator, "
                           "not a bound on the finest state's error.",
    }
    (tmp_path / "hotspot-refinement.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    assert report["nominal_diffusion_times"] >= 5
    steps = np.array([run["dt_max"] for run in measurements])
    assert np.all((steps[1:]/steps[:-1] > 0.4) & (steps[1:]/steps[:-1] < 0.6)), steps
    for field in ("ppar", "pperp"):
        for norm in ("l1", "l2", "linf"):
            for differences in (adjacent, errors):
                ladder = [entry[field][norm] for entry in differences]
                assert 0 < ladder[2] < ladder[1] < ladder[0], (field, norm, ladder)


def _oblique(case, nx, cfl):
    # Both the field and wave have three nonzero components. At one e-fold,
    # missing transverse transport leaves a large, resolution-independent error.
    field = np.array([1., 2., 3.])
    wave = 2*np.pi*np.ones(3)
    chi = np.sqrt(8/np.pi)/(2*np.pi)
    end_time = 1/(chi*(np.dot(field, wave)/np.linalg.norm(field))**2)
    source = (INPUTS / "cgl_lf_rotated_decay.athinput").read_text()
    source += ("\n<mhd>\nnu_coll = 0\nmirror_limiter = false\n"
               "firehose_limiter = false\nbackup_limiters = false\n")
    output = _run(
        case, "cgl_lf_rotated_decay.athinput", "problem/rotated_axis=oblique",
        "problem/b0=1", "problem/by0=2", "problem/bz0=3", "problem/amp=1e-5",
        # The coarse rung is intentionally unresolved; convergence and the finest
        # 1% gate below are the acceptance test, not this single-run smoke bound.
        "problem/decay_rel_tol=0.2", "problem/validation_output=true",
        "problem/validation_output_dir=validation",
        *(f"{block}/nx{direction}={nx}" for block in ("mesh", "meshblock")
          for direction in (1, 2, 3)),
        "mhd/cgl_heat_flux_integrator=explicit", "time/sts_integrator=none",
        f"time/cfl_number={cfl}", "time/nlim=100000", f"time/tlim={end_time:.17g}",
        source=source,
    )
    history, steps = _work(case, output, nx**3, True)
    np.testing.assert_allclose(history["time"][-1], end_time, rtol=0, atol=1e-14)
    with next((case / "validation").glob("*.rotated_decay.csv")).open() as stream:
        values = {row["metric"]: float(row["value"]) for row in csv.DictReader(stream)}
    return values["measured_amp"], values["expected_amp"], steps.max()


def test_cgl_lf_oblique_convergence(tmp_path):
    runs = [_oblique(tmp_path / f"nx{nx}", nx, 0.05) for nx in (16, 32, 64)]
    amplitudes, exact, steps = np.asarray(runs).T
    errors = np.abs(amplitudes-exact)/exact
    print(f"oblique amplitude relative errors: {errors}; actual dt: {steps}")
    assert np.all(errors > 1e-7)
    assert np.all(errors[1:] < 0.45*errors[:-1]), errors
    assert errors[-1] < 0.01
    # Spatial refinement also changes the slope-dependent row envelope.
    assert np.all(steps[1:] < 0.5*steps[:-1]), steps
    # Independently check that time truncation is small compared with spatial error.
    finer, expected, finer_step = _oblique(tmp_path / "half_dt", 32, 0.025)
    # The evolving VL4 derivative changes the bound slightly between trajectories.
    assert 0.45 < finer_step/steps[1] < 0.55
    assert abs(finer-amplitudes[1]) < 0.1*abs(amplitudes[1]-expected)
