"""Perpendicular BGK response and independently coupled temperature decay."""

import csv
import re

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl.test_cgl_lf_geometry_cpu import INPUTS, _run


END_TIME = 0.2972165436900063  # One e-fold of the slow coupled mode.
AMPLITUDE = 1.0e-4


def test_cgl_lf_closure_device_endpoints(tmp_path):
    output = _run(tmp_path / "device", "cgl_lf_geometry.athinput",
                  "problem/test_mode=closure_response")
    assert "closure_response passed" in output


def _continuum(seed):
    # BGK preserves Tpar/2+Tperp and damps Tperp-Tpar at nu. A similarity
    # transform gives a symmetric matrix, independent of the C++ closed form.
    nu, c, k = 10.0, 1.0, 2*np.pi
    chi_par = 8*c*c/(np.sqrt(8*np.pi)*c*k+(3*np.pi-8)*nu)
    chi_perp = c*c/(np.sqrt(np.pi/2)*c*k+nu)
    symmetric = np.array([
        [-chi_par*k*k-2*nu/3, np.sqrt(2)*nu/3],
        [np.sqrt(2)*nu/3, -chi_perp*k*k-nu/3],
    ])
    rates, vectors = np.linalg.eigh(symmetric)
    weights = np.array([1.0, np.sqrt(2)])
    initial = AMPLITUDE*np.eye(2)[seed]
    return (vectors @ (np.exp(rates*END_TIME)*(vectors.T @ (weights*initial))))/weights


def _coupled(case, seed, nx, cfl):
    input_name = f"cgl_lf_quant_{('parallel', 'perp')[seed]}_collisional.athinput"
    source = (INPUTS / input_name).read_text()
    source += ("\n<mhd>\nmirror_limiter = false\nfirehose_limiter = false\n"
               "backup_limiters = false\ncgl_lf_diagnostics = full\n"
               "cgl_lf_arithmetic = safe\n"
               "\n<output1>\nfile_type = hst\ndcycle = 100\ndata_format = %24.16e\n"
               "\n<output2>\nfile_type = tab\nvariable = mhd_w_bcc\n"
               "dcycle = 100000\ndata_format = %24.16e\n")
    output = _run(
        case, input_name, f"mesh/nx1={nx}", f"meshblock/nx1={nx}",
        "mhd/cgl_heat_flux_integrator=explicit", "time/sts_integrator=none",
        f"time/cfl_number={cfl}", "time/nlim=100000", f"time/tlim={END_TIME:.17g}",
        # This smoke bound lets every rung finish on the old closure. The
        # independent two-channel 1% gate below is the acceptance criterion.
        "problem/decay_rel_tol=0.25", "problem/validation_output=true",
        "problem/validation_output_dir=validation", source=source,
    )
    with next((case / "validation").glob("*.decay.csv")).open() as stream:
        values = {row["metric"]: float(row["value"]) for row in csv.DictReader(stream)}
    np.testing.assert_allclose(values["time"], END_TIME, rtol=0, atol=1e-14)
    measured = np.array([values[f"measured_{channel}_amp"] for channel in ("tpar", "tperp")])
    for channel in ("tpar", "tperp"):
        assert abs(values[f"measured_{channel}_cos"]) < 1e-6*AMPLITUDE
        assert abs(values[f"mean_{channel}"]-1.0) < 1e-10

    history = testutils.athena_read.hst(str(next(case.glob("*.mhd.hst"))))
    diagnostics = np.array(re.findall(
        r"elapsed=\S+ cycle=(\d+) time=(\S+) dt=(\S+)", output), dtype=float)
    cycles = int(diagnostics[-1, 0])
    np.testing.assert_array_equal(diagnostics[:, 0], np.arange(cycles+1))
    assert cycles > 100 and history["lf_nstage"][-1] == 2*cycles*nx
    step = history["dt"][0]
    # Background coefficients fix the actual budget; the final update is clipped.
    np.testing.assert_allclose(history["dt"][:-1], step, rtol=1e-12, atol=0)
    np.testing.assert_allclose(diagnostics[:-2, 2], step, rtol=1e-6, atol=0)
    assert 0 < diagnostics[-2, 2] <= step*(1+1e-6)
    np.testing.assert_allclose(diagnostics[:-1, 2].sum(), END_TIME, rtol=1e-6, atol=0)
    np.testing.assert_allclose(history["time"][-1], END_TIME, rtol=0, atol=1e-14)
    for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd", "lf_hwproj"):
        assert np.all(history[name] == 0.0), (name, history[name])
    for name in ("lf_qprwrk", "lf_qpewrk"):
        assert history[name][-1] > 0.0
    np.testing.assert_allclose(history["tot-E"], history["tot-E"][0], rtol=1e-11, atol=0)
    paths = sorted((case / "tab").glob("*.mhd_w_bcc.*.tab"))
    assert len(paths) == 2
    initial, final = (np.loadtxt(path) for path in paths)
    assert initial.shape == final.shape == (nx, 12)
    assert np.all(np.isfinite(final)) and np.min(final[:, [3, 7, 8]]) > 0
    np.testing.assert_array_equal(final[:, 3:7], initial[:, 3:7])
    np.testing.assert_array_equal(final[:, 9:12], initial[:, 9:12])
    print(f"{case.name}: t={values['time']:.17g}, cycles={cycles}, dt={step:.17g}, "
          f"RHS={2*cycles}, signed amplitudes={measured}")
    return measured, step


@pytest.mark.parametrize("seed", [0, 1], ids=["parallel_seed", "perpendicular_seed"])
def test_cgl_lf_coupled_decay(tmp_path, seed):
    expected = _continuum(seed)
    coarse, coarse_dt = _coupled(tmp_path / "nx16", seed, 16, 0.1)
    fine, fine_dt = _coupled(tmp_path / "nx32", seed, 32, 0.1)
    half, half_dt = _coupled(tmp_path / "half_dt", seed, 32, 0.05)
    errors = np.abs(np.array([coarse, fine, half])-expected)/np.abs(expected)
    print(f"coupled continuum amplitudes={expected}; relative channel errors={errors}")
    assert np.all(np.isfinite(errors)) and np.all(errors[1:] < 0.01), errors
    assert np.all(errors[1] < 0.4*errors[0]), errors
    np.testing.assert_allclose([fine_dt, half_dt], [0.25*coarse_dt, 0.5*fine_dt],
                               rtol=1e-12, atol=0)
    # This bounds temporal contamination, not temporal convergence order: the
    # once-per-cycle collision composition can partially cancel spatial error.
    assert np.all(np.abs(half-fine) < 0.2*np.abs(fine-expected))
