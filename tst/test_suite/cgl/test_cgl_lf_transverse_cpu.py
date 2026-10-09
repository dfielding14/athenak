"""Transverse LF slope checks with independent flux and diffusion oracles."""

import csv

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
        "mhd/cgl_lf_profile_detail=true",
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
        "mhd/cgl_lf_profile_detail=true",
    )
    assert "CGL LF hotspot:" in output
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
    np.testing.assert_allclose(steps[1:], 0.25*steps[:-1], rtol=1e-11, atol=1e-16)
    # Independently check that time truncation is small compared with spatial error.
    finer, expected, finer_step = _oblique(tmp_path / "half_dt", 32, 0.025)
    np.testing.assert_allclose(finer_step, 0.5*steps[1], rtol=1e-11, atol=1e-16)
    assert abs(finer-amplitudes[1]) < 0.1*abs(amplitudes[1]-expected)
