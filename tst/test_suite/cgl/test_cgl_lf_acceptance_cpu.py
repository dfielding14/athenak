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


def decay_reference(chi_parallel, chi_perp, nu, time, component):
    # Fourier amplitudes of the two pressure moments: LF conduction plus
    # isotropization that conserves (T_parallel + 2 T_perp)/3.
    k2 = (2.0*np.pi)**2
    matrix = np.array([[-chi_parallel*k2 - 2.0*nu/3.0, 2.0*nu/3.0],
                       [nu/3.0, -chi_perp*k2 - nu/3.0]])
    values, vectors = np.linalg.eig(matrix)
    initial = np.eye(2)[component]
    return (vectors @ (np.exp(values*time)*np.linalg.solve(vectors, initial)))[component]


@pytest.mark.parametrize("component", (0, 1), ids=("parallel", "perp"))
@pytest.mark.parametrize("nu", (0.0, 10.0), ids=("collisionless", "collisional"))
def test_decay_resolves_closure_coefficients(tmp_path, component, nu):
    mode = ("parallel", "perp")[component]
    name = f"cgl_lf_quant_{mode}" + ("_collisional" if nu else "")
    # For this aligned, fixed-coefficient closure, dt_FE = dx**2/(2*chi_parallel).
    # With sts_safety=0.9, a cap of 8 gives 2*128**2/(0.9*8*(2*pi)**2) > 115
    # full-step equivalents per parallel damping time; retain >=100 resolved
    # cycles without changing the physical duration or either closure coefficient.
    result = run_case(tmp_path, name, "time/sts_max_dt_ratio=8",
                      "problem/validation_output=true",
                      f"problem/validation_output_dir={tmp_path}")
    assert result.returncode == 0, result.stdout + result.stderr
    cycles = int(re.findall(r"cycle=(\d+)", result.stdout)[-1])
    assert cycles >= 100
    with next(tmp_path.glob("*.csv")).open() as stream:
        data = {row[0]: float(row[1]) for row in list(csv.reader(stream))[1:]}
    k = 2.0*np.pi
    chi_parallel = 8.0/(np.sqrt(8.0*np.pi)*k + (3.0*np.pi - 8.0)*nu)
    # SHD97/BGK perpendicular moment denominator is +2 nu.
    chi_perp = 2.0/(np.sqrt(2.0*np.pi)*k + 2.0*nu)
    assert abs((chi_parallel, chi_perp)[component]*k*k*data["time"] - 1) < 1e-14
    reference = decay_reference(chi_parallel, chi_perp, nu, data["time"], component)
    measured = data["measured_sin_amp"]/data["initial_amp"]
    tolerance = 3.0e-3
    assert abs(measured/reference - 1) < tolerance
    assert abs(reference - 1) >= 10.0*tolerance
    # These mistakes must lie far outside this test's acceptance window.
    swapped = decay_reference(chi_perp, chi_parallel, nu, data["time"], component)
    assert abs(swapped/reference - 1) > 10.0*tolerance
    if nu and component == 1:
        for wrong_perp in (0.5*chi_perp, 2.0*chi_perp,
                           2.0/(np.sqrt(2.0*np.pi)*k + nu)):
            wrong = decay_reference(chi_parallel, wrong_perp, nu, data["time"], component)
            assert abs(wrong/reference - 1) > 10.0*tolerance


@pytest.mark.parametrize("kind", ("mirror", "firehose"))
@pytest.mark.parametrize("lf", (False, True), ids=("varying_pure", "uniform_lf"))
@pytest.mark.parametrize("nudt", (0.0, 1.0, 1.0e10))
@pytest.mark.parametrize("background_nu", (0.0, 3.0))
def test_limiter_stress_matches_cellwise_relaxation(
    tmp_path, kind, lf, nudt, background_nu
):
    name = f"cgl_lf_limiter_{kind}"
    source = (ROOT / "inputs/unit_tests" / f"{name}.athinput").read_text()
    source = source.replace("<mhd>", f"<mhd>\nnu_coll = {background_nu}", 1)
    if not lf:
        source = re.sub(r"^(?:cgl_heat_flux|sts_integrator)\s*=.*\n", "", source,
                        flags=re.MULTILINE)
    result = run_case(
        tmp_path, name, f"problem/amp={0.0 if lf else 0.1}",
        f"mhd/limiter_nu_coll={nudt/0.005:.17g}", source=source,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "analytic_relaxation=1" in result.stdout


@pytest.mark.parametrize("kind", ("mirror", "firehose"))
@pytest.mark.parametrize("lf", (False, True), ids=("pure", "lf"))
def test_limiter_stress_matches_wall_ordering(tmp_path, kind, lf):
    name = f"cgl_lf_limiter_{kind}"
    source = (ROOT / "inputs/unit_tests" / f"{name}.athinput").read_text()
    if not lf:
        source = re.sub(r"^(?:cgl_heat_flux|sts_integrator)\s*=.*\n", "", source,
                        flags=re.MULTILINE)
    ppar, pperp = (0.5, 2.0) if kind == "mirror" else (3.0, 1.0)
    result = run_case(
        tmp_path, name, "problem/amp=0", f"problem/ppar0={ppar}",
        f"problem/pperp0={pperp}", "mhd/limiter_nu_coll=200",
        "mhd/cgl_lf_strict_admissibility=false", source=source,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "analytic_relaxation=1" in result.stdout
