"""Static contract for the matched 256^3 CGL closure matrix."""

from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
INPUT_ROOT = ROOT / "inputs" / "cgl_256_closure_matrix"

CASES = {
    "cgl256_b10_nolf_nu0": {"nu_coll": 0.0, "k_parallel": None},
    "cgl256_b10_nolf_nu10": {"nu_coll": 10.0, "k_parallel": None},
    "cgl256_b10_lf_k2pi_nu0": {
        "nu_coll": 0.0,
        "k_parallel": 6.283185307179586,
    },
    "cgl256_b10_lf_kpi_nu0": {
        "nu_coll": 0.0,
        "k_parallel": 3.141592653589793,
    },
}

LF_ONLY = {
    "time/sts_integrator",
    "time/sts_max_dt_ratio",
    "mhd/cgl_heat_flux",
    "mhd/cgl_heat_flux_integrator",
    "mhd/cgl_lf_strict_admissibility",
    "mhd/cgl_lf_record_pressure_work",
    "mhd/lf_k_parallel",
    "mhd/lf_coefficient_mode",
}


def _parse(path):
    block = None
    values = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1]
            continue
        assert block is not None and "=" in line, f"invalid line in {path}: {raw_line}"
        key, value = (token.strip() for token in line.split("=", 1))
        flat_key = f"{block}/{key}"
        assert flat_key not in values, f"duplicate parameter in {path}: {flat_key}"
        values[flat_key] = value
    return values


def test_cgl_256_closure_matrix_is_complete_and_matched():
    inputs = {
        path.stem: _parse(path)
        for path in INPUT_ROOT.glob("*.athinput")
    }
    assert set(inputs) == set(CASES)

    invariant = None
    for case, physics in CASES.items():
        values = inputs[case]
        assert values["job/basename"] == case
        assert float(values["mhd/nu_coll"]) == physics["nu_coll"]

        k_parallel = physics["k_parallel"]
        if k_parallel is None:
            assert LF_ONLY.isdisjoint(values)
        else:
            assert LF_ONLY.issubset(values)
            assert values["mhd/cgl_heat_flux"] == "landau_fluid"
            assert values["mhd/cgl_heat_flux_integrator"] == "sts"
            assert values["mhd/cgl_lf_strict_admissibility"] == "true"
            assert float(values["mhd/lf_k_parallel"]) == pytest.approx(k_parallel)

        common = {
            key: value
            for key, value in values.items()
            if key not in LF_ONLY | {"job/basename", "mhd/nu_coll"}
        }
        if invariant is None:
            invariant = common
        else:
            assert common == invariant

    assert invariant is not None
    expected = {
        "mesh/nx1": "256",
        "mesh/nx2": "256",
        "mesh/nx3": "256",
        "mesh/x1max": "1.0",
        "mesh/x2max": "1.0",
        "mesh/x3max": "1.0",
        "time/tlim": "10.0",
        "mhd/eos": "cgl",
        "mhd/passive": "false",
        "mhd/limiter_nu_coll": "1.0e10",
        "mhd/limiter_hardwall": "true",
        "problem/beta0": "10.0",
        "problem/b0": "1.0",
        "turb_driving/projection_policy": "mks24_alfvenic_perpendicular",
        "turb_driving/rseed": "271828",
        "turb_driving/k_shell_unit": "6.283185307179586",
        "turb_driving/dedt": "0.16",
        "turb_driving/tcorr": "1.0",
        "output1/dt": "0.02",
        "output2/variable": "mhd_w_bcc",
        "output2/dt": "0.25",
        "output3/dt": "0.5",
    }
    for key, value in expected.items():
        assert invariant[key] == value
