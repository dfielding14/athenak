"""Resolved CGL/LF waves over noninteger periods, with retained Fourier oracles."""

import csv
import re

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl.test_cgl_lf_geometry_cpu import INPUTS, _run


WAVE_CASES = [
    f"cgl_{closure}_paper_eigen_{branch}"
    for closure in ("pure", "lf") for branch in ("alfven", "slow", "fast")
] + ["cgl_lf_field_wave", "cgl_lf_paper_oblique_wave", "cgl_pure_paper_oblique_wave"]


@pytest.mark.parametrize("name", WAVE_CASES)
def test_wave_convergence_over_noninteger_periods(tmp_path, name):
    source = (INPUTS / f"{name}.athinput").read_text()
    eigen = "eigen_lambda_im" in source
    if eigen:
        omega = abs(float(re.search(r"eigen_lambda_im\s*=\s*(\S+)", source)[1]))
    elif "field_wave" in name:
        # Linearized continuity, parallel momentum, and CGL pressure equations.
        k = 2*np.pi
        chi = np.sqrt(8/np.pi)/k
        matrix = np.array([[0, -1j*k, 0], [0, 0, -1j*k],
                           [chi*k*k, -3j*k, -chi*k*k]])
        omega = max(abs(np.linalg.eigvals(matrix).imag))
    else:
        # The oblique initial-value problem contains several modes.
        omega = 2*np.pi
    final_time = 1.25*2*np.pi/omega
    tolerance_key = "eigen_wave_rel_tol" if eigen else "wave_rel_tol"
    source += ("\n<time>\nsts_safety = 0.3\n"
               "\n<output1>\nfile_type = hst\ndcycle = 1\ndata_format = %24.16e\n"
               "\n<output2>\nfile_type = tab\nvariable = mhd_w_bcc\n"
               "dcycle = 1000000\ndata_format = %24.16e\n")
    errors = []
    for nx, tolerance in ((64, 0.02), (128, 0.005), (256, 0.001)):
        case = tmp_path / str(nx)
        output = _run(
            case, f"{name}.athinput", f"mesh/nx1={nx}", f"meshblock/nx1={nx}",
            "time/nlim=-1", f"time/tlim={final_time:.17g}",
            f"problem/{tolerance_key}={tolerance}", "problem/validation_output=true",
            "problem/validation_output_dir=validation", source=source,
        )
        with next((case / "validation").glob("*.csv")).open() as stream:
            rows = list(csv.DictReader(stream))
        errors.append(max(float(row["rel_err"]) for row in rows
                          if row.get("required", "1") == "1"))
        history = testutils.athena_read.hst(str(next(case.glob("*.mhd.hst"))))
        np.testing.assert_allclose(history["time"][-1], final_time, rtol=0, atol=1e-14)
        cycles = int(re.findall(r"cycle=(\d+)", output)[-1])
        assert cycles > 100
        if "lf_nstage" in history:
            assert history["lf_nstage"][-1] >= 6*cycles*nx
            for key in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos",
                        "lf_hardbd", "lf_hwproj"):
                assert np.all(history[key] == 0), key
        snapshots = sorted((case / "tab").glob("*.mhd_w_bcc.*.tab"))
        assert len(snapshots) == 2
        final = np.loadtxt(snapshots[-1])
        assert final.shape == (nx, 12) and np.all(np.isfinite(final))
        assert np.min(final[:, [3, 7, 8]]) > 0
    orders = np.log2(np.array(errors[:-1])/errors[1:])
    print(f"{name}: endpoint={final_time:.17g}, errors={errors}, orders={orders}")
    # Retain the archived physical thresholds; LF-dominated waves reach first order.
    minimum_order = 1.8 if "pure" in name or "alfven" in name else 1.0
    assert np.all(orders >= minimum_order), (name, errors, orders)
