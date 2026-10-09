"""Independent face geometry and frozen-field LF evolution checks."""

import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

from test_suite import testutils


ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "inputs/unit_tests"


def _run(case, input_name, *flags, source=None):
    case.mkdir()
    input_path = INPUTS / input_name
    if source is not None:
        input_path = case / "input.athinput"
        input_path.write_text(source)
    env = os.environ.copy()
    for name in ("PROFILE", "PROFILE_DETAIL", "DIAGNOSTICS", "ARITHMETIC", "STS_FLUX"):
        env.pop(f"ATHENAK_CGL_LF_{name}", None)
    result = subprocess.run(
        [str(Path("athena").resolve()), "-i", str(input_path), *flags],
        cwd=case, env=env, capture_output=True, text=True, check=False,
    )
    output = result.stdout + result.stderr
    (case / "run.log").write_text(output)
    assert result.returncode == 0, output
    return output


def _work(case, output, cells, explicit):
    history = testutils.athena_read.hst(str(next(case.glob("*.mhd.hst"))))
    cycles = int(re.findall(r"cycle=(\d+)", output)[-1])
    steps = np.diff(history["time"])
    assert len(steps) == cycles and np.all(steps > 0.0)
    rhs = history["lf_nstage"][-1] / cells
    assert rhs == int(rhs) and rhs >= 2*cycles
    if explicit:
        assert rhs == 2*cycles
    # Fixed density and background coefficients keep the pre-piece-9 budget fixed.
    if len(steps) > 2:
        np.testing.assert_allclose(steps[:-1], steps[0], rtol=1e-11, atol=1e-16)
    for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd",
                 "lf_hwproj"):
        assert np.all(history[name] == 0.0), (name, history[name])
    print(f"{case.name}: time={history['time'][-1]:.17g} cycles={cycles} "
          f"dt=[{steps.min():.17g},{steps.max():.17g}] RHS={rhs:g}")
    return history, steps


@pytest.mark.parametrize("arithmetic", ["safe", "fast"])
@pytest.mark.parametrize("diagnostics", ["full", "none"])
def test_cgl_lf_face_geometry(tmp_path, arithmetic, diagnostics):
    output = _run(
        tmp_path / "faces", "cgl_lf_geometry.athinput",
        f"mhd/cgl_lf_arithmetic={arithmetic}",
        f"mhd/cgl_lf_diagnostics={diagnostics}",
        "mhd/cgl_lf_profile=true", "mhd/cgl_lf_profile_detail=true",
    )
    assert "face_geometry passed" in output


@pytest.mark.parametrize("geometry", ["reversal", "ct"])
def test_cgl_lf_geometry_negative_control(tmp_path, geometry):
    """Keep the two defects independently selectable before applying the correction."""
    output = _run(tmp_path / geometry, "cgl_lf_geometry.athinput",
                  f"problem/geometry_case={geometry}")
    assert "face_geometry passed" in output


@pytest.mark.parametrize("dimension", [1, 2])
@pytest.mark.parametrize("integrator", ["explicit", "sts"])
def test_cgl_lf_reversal_stability(tmp_path, dimension, integrator):
    case = tmp_path / f"reversal{dimension}d_{integrator}"
    flags = [f"mhd/cgl_heat_flux_integrator={integrator}"]
    if dimension == 2:
        flags += ["mesh/nx2=16", "meshblock/nx2=16"]
    if integrator == "sts":
        flags += ["time/sts_integrator=rkl2", "time/cfl_number=0.4"]
    output = _run(case, "cgl_lf_field_reversal.athinput", *flags)
    assert "field_reversal" in output
    cells = 64 if dimension == 1 else 64*16
    history, _ = _work(case, output, cells, integrator == "explicit")
    assert len(history["time"]) == 21
    report = next(line for line in output.splitlines()
                  if line.startswith("CGL LF field_reversal:"))
    totals = {name: float(value) for name, value in re.findall(r"(\w+)=(\S+)", report)}
    np.testing.assert_allclose([totals["energy"], totals["moment"]],
                               [totals["energy0"], totals["moment0"]], rtol=2e-12, atol=0)
    # The C++ monitor retains the original positive-pressure and <10 seed-growth gates.
    if integrator == "sts":
        assert history["lf_nstage"][-1] > cells*20*2


def _resolved_reversal(case, nx, cfl):
    source = (INPUTS / "cgl_lf_field_reversal.athinput").read_text()
    for number, variable in ((2, "mhd_w_bcc"), (3, "mhd_u")):
        source += (f"\n<output{number}>\nfile_type = tab\nvariable = {variable}\n"
                   "data_format = %24.16e\ndcycle = 100000\n")
    output = _run(
        case, "cgl_lf_field_reversal.athinput", f"mesh/nx1={nx}",
        f"meshblock/nx1={nx}", "mesh/ix1_bc=periodic", "mesh/ox1_bc=periodic",
        "problem/reversal_width=0.1", "problem/b0=0.5", "problem/amp=1.0e-4",
        "time/nlim=100000", "time/tlim=0.002", f"time/cfl_number={cfl}", source=source,
    )
    history, steps = _work(case, output, nx, True)
    assert abs(history["time"][-1] - 0.002) < 1e-15
    paths = sorted((case / "tab").glob("*.mhd_w_bcc.*.tab"))
    assert len(paths) == 2
    initial, final = np.loadtxt(paths[0]), np.loadtxt(paths[-1])
    assert final.shape == initial.shape == (nx, 12)
    assert np.all(np.isfinite(final)) and np.min(final[:, [3, 7, 8]]) > 0.0
    np.testing.assert_array_equal(final[:, 3:7], initial[:, 3:7])
    np.testing.assert_array_equal(final[:, 9:12], initial[:, 9:12])
    initial_mu = initial[:, 8] / np.linalg.norm(initial[:, 9:12], axis=1)
    final_mu = final[:, 8] / np.linalg.norm(final[:, 9:12], axis=1)
    np.testing.assert_allclose(final_mu.mean(), initial_mu.mean(), rtol=0, atol=1e-11)
    conserved = sorted((case / "tab").glob("*.mhd_u.*.tab"))
    assert len(conserved) == 2
    initial_u, final_u = np.loadtxt(conserved[0]), np.loadtxt(conserved[-1])
    np.testing.assert_allclose(final_u[:, 7].mean(), initial_u[:, 7].mean(),
                               rtol=0, atol=1e-11)
    return final[:, 7:9], steps.max()


def test_cgl_lf_resolved_reversal_converges(tmp_path):
    """Refine fixed physical sheets and a sine seed, independently of the time step."""
    spatial = []
    for nx in (32, 64, 128):
        state, base_step = _resolved_reversal(tmp_path / f"nx{nx}", nx, 0.2)
        spatial.append(state)
    errors = [np.mean(np.abs(coarse - fine.reshape(len(coarse), 2, 2).mean(axis=1)),
                      axis=0)
              for coarse, fine in zip(spatial[:-1], spatial[1:])]
    errors = np.asarray(errors)
    print(f"resolved reversal spatial ppar/pperp errors: {errors}")
    assert np.all(errors > 1e-14)
    assert np.all(errors[1] < 0.6*errors[0]), errors

    temporal = [spatial[-1]]
    steps = [base_step]
    for cfl in (0.1, 0.05):
        state, step = _resolved_reversal(tmp_path / f"dt{cfl}", 128, cfl)
        temporal.append(state)
        steps.append(step)
    np.testing.assert_allclose(steps[1:], 0.5*np.asarray(steps[:-1]),
                               rtol=1e-11, atol=1e-16)
    errors = np.array([np.mean(np.abs(coarse - fine), axis=0)
                       for coarse, fine in zip(temporal[:-1], temporal[1:])])
    print(f"resolved reversal temporal ppar/pperp errors: {errors}")
    assert np.all(errors > 1e-14)
    # The explicit reference uses two forward-Euler half-sweeps per fluid cycle.
    assert np.all(errors[1] < 0.7*errors[0]), errors


@pytest.mark.parametrize("axis,field", [
    ("x", (1.0, 0.0, 0.0)), ("y", (0.0, 1.0, 0.0)),
    ("z", (0.0, 0.0, 1.0)), ("oblique", (1.0, 2.0, 3.0)),
])
def test_cgl_lf_uniform_oblique_decay(tmp_path, axis, field):
    case = tmp_path / axis
    # Resolve the limited transverse slopes to the existing 3% e-fold bound.
    nx = 32 if axis == "oblique" else 16
    wave = 2*np.pi*np.array([1.0, 1.0, 1.0] if axis == "oblique" else
                            [float(axis == direction) for direction in ("x", "y", "z")])
    k_parallel = np.dot(field, wave) / np.linalg.norm(field)
    chi = np.sqrt(8/np.pi) / (2*np.pi)
    end_time = 1/(chi*k_parallel**2)
    source = (INPUTS / "cgl_lf_rotated_decay.athinput").read_text()
    source += ("\n<mhd>\nnu_coll = 0\nmirror_limiter = false\n"
               "firehose_limiter = false\nbackup_limiters = false\n")
    output = _run(
        case, "cgl_lf_rotated_decay.athinput", f"problem/rotated_axis={axis}",
        f"problem/b0={field[0]}", f"problem/by0={field[1]}", f"problem/bz0={field[2]}",
        *(f"{block}/nx{direction}={nx}" for block in ("mesh", "meshblock")
          for direction in (1, 2, 3)),
        "mhd/cgl_heat_flux_integrator=explicit", "time/sts_integrator=none",
        "time/cfl_number=0.1", "time/nlim=100000", f"time/tlim={end_time:.17g}",
        source=source,
    )
    assert "rotated_decay passed" in output
    history, _ = _work(case, output, nx**3, True)
    np.testing.assert_allclose(history["time"][-1], end_time, rtol=0, atol=1e-14)
