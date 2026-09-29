"""Long, quantitative CGL/LF acceptance checks with independent Fourier references."""

import csv
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
WAVE_CASES = [
    f"cgl_{closure}_paper_eigen_{branch}"
    for closure in ("pure", "lf") for branch in ("alfven", "slow", "fast")
] + ["cgl_lf_field_wave", "cgl_lf_paper_oblique_wave", "cgl_pure_paper_oblique_wave"]


def run_case(tmp_path, name, *flags, source=None):
    """Keep every run's output isolated, including deliberate failures."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    input_file = ROOT / "inputs/unit_tests" / f"{name}.athinput"
    if source is not None:
        input_file = tmp_path / "input.athinput"
        input_file.write_text(source)
    return subprocess.run(
        [str(Path("athena").resolve()), "-i", str(input_file),
         "time/ndiag=100000", *flags], cwd=tmp_path,
        capture_output=True, text=True, check=False,
    )


@pytest.mark.parametrize("name", WAVE_CASES)
def test_short_wave_is_resolved(tmp_path, name):
    result = run_case(tmp_path, name)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("name", WAVE_CASES)
def test_wave_convergence_over_noninteger_periods(tmp_path, name):
    source = (ROOT / "inputs/unit_tests" / f"{name}.athinput").read_text()
    eigen = "eigen_lambda_im" in source
    if eigen:
        omega = abs(float(re.search(r"eigen_lambda_im\s*=\s*(\S+)", source)[1]))
    elif "field_wave" in name:
        # Linearized continuity, parallel momentum, and CGL pressure equations.
        k = 2.0*np.pi
        chi = np.sqrt(8.0/np.pi)/k
        matrix = np.array([[0, -1j*k, 0], [0, 0, -1j*k],
                           [chi*k*k, -3j*k, -chi*k*k]])
        omega = max(abs(np.linalg.eigvals(matrix).imag))
    else:
        # The oblique IVP contains several modes; use its Alfven period.
        omega = 2.0*np.pi
    final_time = 1.25*2.0*np.pi/omega
    tolerance_key = "eigen_wave_rel_tol" if eigen else "wave_rel_tol"
    errors = []
    for nx, tolerance in ((64, 0.02), (128, 0.005), (256, 0.001)):
        run_dir = tmp_path / str(nx)
        result = run_case(
            run_dir, name, f"mesh/nx1={nx}", f"meshblock/nx1={nx}",
            f"time/tlim={final_time:.17g}", f"problem/{tolerance_key}={tolerance}",
            "problem/validation_output=true",
            f"problem/validation_output_dir={run_dir}",
        )
        assert result.returncode == 0, result.stdout + result.stderr
        with next(run_dir.glob("*.csv")).open() as stream:
            rows = list(csv.DictReader(stream))
        errors.append(max(float(row["rel_err"]) for row in rows
                          if row.get("required", "1") == "1"))
    orders = np.log2(np.array(errors[:-1])/errors[1:])
    # Weakly damped waves are second order; LF-dominated waves must reach first.
    minimum_order = 1.8 if "pure" in name or "alfven" in name else 1.0
    assert np.all(orders >= minimum_order), (name, errors, orders)


@pytest.mark.parametrize("name", ["cgl_lf_field_wave", "cgl_lf_paper_oblique_wave",
                                  "cgl_lf_paper_eigen_slow"])
def test_wave_rejects_ineffective_reference(tmp_path, name):
    result = run_case(tmp_path, name, "time/tlim=1e-8")
    assert result.returncode != 0
    assert "a frozen state could pass" in result.stdout


@pytest.mark.parametrize("branch", ("slow", "fast"))
def test_eigen_wave_rejects_disabled_lf(tmp_path, branch):
    name = f"cgl_lf_paper_eigen_{branch}"
    source = (ROOT / "inputs/unit_tests" / f"{name}.athinput").read_text()
    source = re.sub(r"^(?:cgl_heat_flux|sts_integrator)\s*=.*\n", "", source,
                    flags=re.MULTILINE)
    result = run_case(tmp_path, name, source=source)
    assert result.returncode != 0
    assert "paper_eigen_wave" in result.stdout and "rel_err=" in result.stdout
