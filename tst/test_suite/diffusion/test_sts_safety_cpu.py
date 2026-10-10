"""Runtime checks of separate STS safety, ordinary diffusion, and restart loading."""

import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl.test_cgl_lf_geometry_cpu import INPUTS, _run


UNIFORM = "cgl_lf_uniform_timestep.athinput"
CFL = 0.31234567890123457
CHI = np.sqrt(8/np.pi)/(2*np.pi)


def _history(case, fluid="mhd"):
    return testutils.athena_read.hst(str(next(case.glob(f"*.{fluid}.hst"))))


def _command(case, arguments):
    case.mkdir()
    env = os.environ.copy()
    for name in ("PROFILE", "PROFILE_DETAIL", "DIAGNOSTICS", "ARITHMETIC", "STS_FLUX"):
        env.pop(f"ATHENAK_CGL_LF_{name}", None)
    result = subprocess.run([str(Path("athena").resolve()), *arguments], cwd=case,
                            env=env, capture_output=True, text=True, check=False)
    (case / "run.log").write_text(result.stdout + result.stderr)
    return result


def _odd_stages(sweep_dt, budget):
    # RKL2's Legendre stability interval, independent of the controller's sqrt rule.
    stages = 3
    while (stages*stages + stages - 2)*budget < 4*sweep_dt:
        stages += 2
    return stages


@pytest.mark.parametrize("value", ["0", "-0.1", "1.01", "nan", "inf", "-inf"])
def test_invalid_sts_safety_is_rejected(tmp_path, value):
    staged = tmp_path / "invalid.athinput"
    staged.write_text((INPUTS / UNIFORM).read_text() + f"\n<time>\nsts_safety = {value}\n")
    result = _command(tmp_path / "invalid", ["-i", str(staged)])
    assert result.returncode != 0
    assert "<time>/sts_safety must be finite and in (0, 1]" in result.stdout + result.stderr


def test_absent_sts_safety_preserves_current_cfl(tmp_path):
    histories = []
    for label, safety in (("absent", ""), ("explicit", f"sts_safety = {CFL:.17g}\n")):
        case = tmp_path / label
        source = (INPUTS / UNIFORM).read_text() + f"\n<time>\n{safety}"
        # Nonzero anisotropy keeps the pgen's relative relaxation oracle well scaled.
        _run(case, UNIFORM, f"time/cfl_number={CFL:.17g}", "time/nlim=2",
             "mhd/nu_coll=0", "problem/pperp0=1.5", source=source)
        histories.append(_history(case))
    for name in histories[0]:
        np.testing.assert_array_equal(histories[0][name], histories[1][name])
    expected = 20*CFL/(2*16**2*CHI)
    np.testing.assert_allclose(histories[0]["dt"], expected, rtol=2e-12, atol=0)
    assert histories[0]["lf_nstage"][-1] == 2*2*7*16


@pytest.mark.parametrize("integrator", ["sts", "explicit"])
def test_sts_safety_changes_only_selected_diffusion_budget(tmp_path, integrator):
    source = (INPUTS / UNIFORM).read_text() + "\n<time>\nsts_safety = 0.9\n"
    for safety in (0.45, 0.9, 1.0):
        case = tmp_path / f"safety_{safety:g}"
        _run(case, UNIFORM, f"time/sts_safety={safety}", "time/cfl_number=0.3",
             f"mhd/cgl_heat_flux_integrator={integrator}", "mhd/nu_coll=0",
             f"time/sts_integrator={'rkl2' if integrator == 'sts' else 'none'}",
             source=source)
        expected = (20*safety if integrator == "sts" else 0.3)/(2*16**2*CHI)
        np.testing.assert_allclose(_history(case)["dt"][0], expected, rtol=2e-12, atol=0)


def test_sts_safety_leaves_advective_cfl_unchanged(tmp_path):
    source = (INPUTS / UNIFORM).read_text() + "\n<time>\nsts_safety = 0.9\n"
    for cfl, safety in ((0.3, 0.45), (0.3, 0.9), (0.15, 0.9)):
        case = tmp_path / f"cfl_{cfl:g}_safety_{safety:g}"
        _run(case, UNIFORM, f"time/sts_safety={safety}", f"time/cfl_number={cfl}",
             "time/evolution=dynamic", "mhd/rsolver=hlle", "mhd/nu_coll=0",
             "mhd/lf_k_parallel=1e8", source=source)
        # Isotropic unit pressure/density and Bx=1 give parallel CGL speed sqrt(3).
        np.testing.assert_allclose(_history(case)["dt"][0], cfl/(16*np.sqrt(3)),
                                   rtol=2e-12, atol=0)


@pytest.mark.parametrize("safety", [0.3, 0.9])
def test_sts_safety_controls_initial_and_refreshed_stages(tmp_path, safety):
    name = "cgl_lf_timestep_refresh.athinput"
    source = (INPUTS / name).read_text() + f"\n<time>\nsts_safety = {safety}\n"
    case = tmp_path / "heated"
    output = _run(case, name, source=source)
    history = _history(case)
    raw_budget = 1/(2*64**2*CHI)
    cycle_dt = 20*safety*raw_budget
    pressure = 1 + (2/3)*3000*cycle_dt
    expected = [_odd_stages(cycle_dt/2, safety*raw_budget),
                _odd_stages(cycle_dt/2, safety*raw_budget/np.sqrt(pressure))]
    match = re.search(r"timestep_refresh: pre_stages=(\d+) post_stages=(\d+)", output)
    assert match is not None, output
    assert [int(value) for value in match.groups()] == expected
    assert expected[1] > expected[0]
    np.testing.assert_allclose(history["dt"][0], cycle_dt, rtol=2e-12, atol=0)
    np.testing.assert_allclose(history["time"][-1], cycle_dt, rtol=2e-12, atol=0)
    assert history["lf_nstage"][-1] == sum(expected)*64
    np.testing.assert_allclose(history["tot-E"][-1] - history["tot-E"][0],
                               3000*cycle_dt, rtol=2e-12, atol=0)
    for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hwproj"):
        assert np.all(history[name] == 0), name


@pytest.mark.parametrize("safety", [None, 0.9], ids=["default", "explicit"])
def test_sts_safety_restart_loading(tmp_path, safety):
    source = (INPUTS / UNIFORM).read_text()
    source += ("\n<output2>\nfile_type = rst\ndcycle = 1\n"
               f"\n<time>\ncfl_number = {CFL:.17g}\n")
    if safety is not None:
        source += f"sts_safety = {safety}\n"
    original = tmp_path / "original"
    # Uniform background coefficients keep the analytic budget independent of anisotropy.
    _run(original, UNIFORM, "time/nlim=1", "mhd/nu_coll=0",
         "problem/pperp0=1.5", source=source)
    checkpoint = sorted((original / "rst").glob("*.rst"))[-1]
    record = checkpoint.read_bytes().split(b"<par_end>", 1)[0].decode("ascii")
    marker = re.search(r"(?m)^sts_safety\s*=\s*(\S+)", record)
    if safety is not None:
        assert marker is not None and float(marker.group(1)) == safety
    elif marker is not None:
        # An implementation may persist the default, but must not round it to six digits.
        assert float(marker.group(1)) == CFL
    original_dt = 20*(CFL if safety is None else safety)/(2*16**2*CHI)
    np.testing.assert_allclose(_history(original)["dt"], original_dt, rtol=2e-12, atol=0)
    for label, cfl in (("unchanged", CFL), ("changed_cfl", 0.2)):
        case = tmp_path / label
        result = _command(case, ["-r", str(checkpoint), "time/nlim=3",
                                f"time/cfl_number={cfl:.17g}"])
        assert result.returncode == 0, result.stdout + result.stderr
        history = _history(case)
        expected_safety = cfl if safety is None else safety
        expected_dt = 20*expected_safety/(2*16**2*CHI)
        np.testing.assert_allclose(history["dt"], expected_dt, rtol=2e-12, atol=0)
        np.testing.assert_allclose(history["time"][-1], original_dt + 2*expected_dt,
                                   rtol=2e-12, atol=0)
        assert history["lf_nstage"][-1] == 3*2*7*16
        for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hwproj"):
            assert np.all(history[name] == 0), name


@pytest.mark.parametrize("integrator", ["sts", "explicit"])
def test_sts_safety_applies_to_ordinary_viscosity(tmp_path, integrator):
    source = (INPUTS.parent / "tests/sts_viscous_shear.athinput").read_text()
    source += ("\n<time>\nsts_safety = 0.9\nnlim = 0\ntlim = 1\n"
               "sts_max_dt_ratio = 2\n<hydro>\ntdep_viscosity = false\n"
               "viscosity = 1\n<problem>\namp = 0\n"
               "<output2>\nfile_type = hst\ndcycle = 1\ndata_format = %24.16e\n")
    for safety in (0.45, 0.9):
        case = tmp_path / f"viscosity_{safety:g}"
        _run(case, UNIFORM, "mesh/nx1=32", "meshblock/nx1=32",
             f"time/sts_safety={safety}", "time/cfl_number=0.3",
             f"hydro/viscosity_integrator={integrator}",
             f"time/sts_integrator={'rkl2' if integrator == 'sts' else 'none'}",
             source=source)
        expected = (2*safety if integrator == "sts" else 0.3)/(2*32**2)
        np.testing.assert_allclose(_history(case, "hydro")["dt"][0], expected,
                                   rtol=2e-12, atol=0)


def test_absent_sts_safety_preserves_non_sts_rk4_cfl(tmp_path):
    source = (INPUTS.parent / "tests/sts_viscous_shear.athinput").read_text()
    source += ("\n<time>\nintegrator = rk4\ncfl_number = 1.2\nnlim = 0\ntlim = 1\n"
               "sts_integrator = none\n<hydro>\ntdep_viscosity = false\n"
               "viscosity = 1\nviscosity_integrator = explicit\n<problem>\namp = 0\n"
               "<output2>\nfile_type = hst\ndcycle = 1\ndata_format = %24.16e\n")
    case = tmp_path / "rk4"
    _run(case, UNIFORM, "mesh/nx1=32", "meshblock/nx1=32", source=source)
    np.testing.assert_allclose(_history(case, "hydro")["dt"][0], 1.2/(2*32**2),
                               rtol=2e-12, atol=0)


@pytest.mark.parametrize("safety", [0.3, 1.0])
def test_sts_safety_keeps_zero_operator_inactive(tmp_path, safety):
    source = (INPUTS.parent / "tests/cgl_lf_decay.athinput").read_text()
    source += f"\n<time>\nsts_safety = {safety}\n"
    case = tmp_path / "inactive"
    _run(case, UNIFORM, "time/evolution=kinematic", "time/nlim=1",
         "time/tlim=1e-5", "mhd/rsolver=advect", "problem/test_mode=low_field",
         "problem/b0=5e-11", "problem/amp=0.5", source=source)
    history = _history(case)
    np.testing.assert_allclose(history["time"][-1], 1e-5, rtol=0, atol=1e-19)
    for name in ("lf_nstage", "lf_qface", "lf_qprwrk", "lf_qpewrk", "lf_dfloor",
                 "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd", "lf_hwproj"):
        assert np.all(history[name] == 0), name
