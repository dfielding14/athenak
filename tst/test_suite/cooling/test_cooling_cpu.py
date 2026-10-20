"""
Regression tests for the general <cooling> source term.

The tests use the built-in cooling_test pgen, which initializes a uniform 1D
domain.  With no fluid gradients, the history cooling columns have direct
analytic expectations and can be checked without contamination from hydro fluxes.
"""

from pathlib import Path
import subprocess

import numpy as np
import pytest

import test_suite.testutils as testutils

import sys

sys.path.insert(0, "../vis/python")
import athena_read  # noqa: E402


HYDRO_INPUT = "inputs/cooling_test_hydro.athinput"
MHD_INPUT = "inputs/cooling_test_mhd.athinput"
ROOT = Path(__file__).resolve().parents[3]
VALIDATOR = ROOT / "tools" / "validate_cooling_table.py"
AMU_CGS = 1.660538921e-24


def _remove_outputs(basename, physics):
    Path(f"{basename}.{physics}.hst").unlink(missing_ok=True)


def _run_case(basename, flags=None, input_file=HYDRO_INPUT, physics="hydro"):
    flags = [] if flags is None else list(flags)
    _remove_outputs(basename, physics)
    testutils.run(input_file, [f"job/basename={basename}"] + flags)
    return athena_read.hst(f"{basename}.{physics}.hst", raw=True)


def _first_evolved_row(data):
    for row in range(1, len(data["time"])):
        if data["time"][row] > data["time"][row - 1]:
            return row
    pytest.fail("history file did not contain an evolved row")


def _assert_history_rate(data, gross, net, rtol=2.0e-12, atol=1.0e-13):
    row = _first_evolved_row(data)
    assert np.isclose(data["cool_gross"][row], gross, rtol=rtol, atol=atol)
    assert np.isclose(data["cool_net"][row], net, rtol=rtol, atol=atol)
    dt = data["time"][row] - data["time"][row - 1]
    energy_loss_rate = (data["tot-E"][row - 1] - data["tot-E"][row])/dt
    assert np.isclose(data["cool_net"][row], energy_loss_rate, rtol=1.0e-10,
                      atol=1.0e-12)


@pytest.mark.parametrize("integrator", ["rk1", "rk2", "rk3", "rk4"])
def test_constant_powerlaw_history_weights(integrator):
    """Constant cooling should give the same interval rate for all RK schemes."""
    data = _run_case(
        f"cooling_powerlaw_{integrator}",
        [f"time/integrator={integrator}"],
    )
    _assert_history_rate(data, gross=1.0e-2, net=1.0e-2)


@pytest.mark.parametrize("physics", ["hydro", "mhd"])
def test_rk4_temperature_dependent_history(physics):
    """History retains the actual stage cooling over several timesteps."""
    data = _run_case(
        f"cooling_rk4_interval_{physics}",
        ["time/integrator=rk4", "time/nlim=8", "time/tlim=1.0",
         "cooling/cooling_slope=2.0", "cooling/cooling_reference_value=1.0",
         "output1/dcycle=4"],
        input_file=HYDRO_INPUT if physics == "hydro" else MHD_INPUT,
        physics=physics,
    )
    elapsed = np.diff(data["time"])
    valid = elapsed > 0.0
    assert np.count_nonzero(valid) >= 2
    energy_loss = -np.diff(data["tot-E"])[valid]
    np.testing.assert_allclose(data["cool_net"][1:][valid]*elapsed[valid],
                               energy_loss, rtol=1.0e-10, atol=1.0e-12)
    np.testing.assert_allclose(data["cool_gross"], data["cool_net"],
                               rtol=1.0e-12, atol=1.0e-13)


@pytest.mark.parametrize(
    "basename,flags,expected_gross,expected_net",
    [
        (
            "cooling_table_1d",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_1d.tbl"],
            1.0e-2,
            1.0e-2,
        ),
        (
            "cooling_table_2d",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_2d.tbl"],
            1.0e-2,
            1.0e-2,
        ),
        (
            "cooling_table_3d",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_3d.tbl"],
            1.0e-2,
            1.0e-2,
        ),
        (
            "cooling_table_default_log_axis",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_default_log.tbl",
             "problem/pressure=10.0"],
            3.0e-2,
            3.0e-2,
        ),
        (
            "cooling_table_mixed_axes",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_mixed_axes.tbl",
             "problem/pressure=10.0"],
            3.0e-2,
            3.0e-2,
        ),
        (
            "cooling_piecewise_powerlaw",
            ["cooling/cooling_model=piecewise_powerlaw",
             "cooling/cooling_breaks=0.5,1.5",
             "cooling/cooling_slopes=0.0,0.0,0.0"],
            1.0e-2,
            1.0e-2,
        ),
        (
            "cooling_heating_constant",
            ["cooling/heating_model=constant", "cooling/heating_rate=2.0e-3"],
            1.0e-2,
            8.0e-3,
        ),
        (
            "cooling_heating_table",
            ["cooling/heating_model=table",
             "cooling/heating_table=inputs/tables/cooling_gamma_1d.tbl"],
            1.0e-2,
            8.0e-3,
        ),
        (
            "cooling_modifier_constant",
            ["cooling/cooling_modifier_model=constant",
             "cooling/cooling_modifier=2.0"],
            2.0e-2,
            2.0e-2,
        ),
        (
            "cooling_modifier_table",
            ["cooling/cooling_modifier_model=table",
             "cooling/cooling_modifier_table=inputs/tables/cooling_modifier_1d.tbl"],
            2.0e-2,
            2.0e-2,
        ),
        (
            "cooling_user_hook",
            ["cooling/cooling_model=user",
             "cooling/heating_model=user",
             "problem/register_user_cooling=true",
             "problem/user_lambda=1.0e-2",
             "problem/user_gamma=2.0e-3"],
            1.0e-2,
            8.0e-3,
        ),
        (
            "cooling_cgm_no_shield",
            ["cooling/cooling_model=cgm",
             "cooling/cgm_pie_table=inputs/tables/cooling_cgm_pie.tbl",
             "cooling/cgm_cie_table=inputs/tables/cooling_cgm_cie.tbl",
             "cooling/cooling_density=mass_density",
             "cooling/shielding_model=none"],
            1.0e-2,
            1.0e-2,
        ),
    ],
)
def test_cooling_models_and_history_columns(basename, flags, expected_gross,
                                            expected_net):
    data = _run_case(basename, flags)
    _assert_history_rate(data, gross=expected_gross, net=expected_net,
                         rtol=1.0e-5, atol=1.0e-12)


def test_cgm_shielding_uses_composition_metadata():
    data = _run_case(
        "cooling_cgm_shielded",
        input_file="inputs/cooling_test_hydro_composition.athinput",
    )
    _assert_history_rate(data, gross=3.0e-2, net=3.0e-2)


def test_mhd_source_path_and_history():
    data = _run_case("cooling_mhd_powerlaw", input_file=MHD_INPUT, physics="mhd")
    _assert_history_rate(data, gross=1.0e-2, net=1.0e-2)


def test_scalar_axis_uses_selected_passive_scalar():
    data = _run_case(
        "cooling_scalar_axis_valid",
        ["cooling/cooling_model=table",
         "cooling/cooling_table=inputs/tables/cooling_lambda_3d.tbl",
         "problem/scalar0=1.0",
         "time/tlim=1.0e-8",
         "time/nlim=1"],
    )
    row = _first_evolved_row(data)
    assert np.isclose(data["cool_gross"][row], 1.1e-2, rtol=1.0e-8,
                      atol=1.0e-12)
    assert np.isclose(data["cool_net"][row], 1.1e-2, rtol=1.0e-8,
                      atol=1.0e-12)


@pytest.mark.parametrize(
    "basename,density_kind,composition_flags,expected_rate",
    [
        (
            "cooling_cgs_number_density_units",
            "number_density",
            ["cooling/mu=0.5"],
            2.0e-2,
        ),
        (
            "cooling_cgs_hydrogen_density_units",
            "hydrogen_number_density",
            ["cooling/hydrogen_mass_fraction=0.4"],
            3.0e-2,
        ),
    ],
)
def test_cgs_density_conversion_with_nontrivial_units(
        basename, density_kind, composition_flags, expected_rate):
    length_cgs = 2.0
    mass_cgs = 16.0
    time_cgs = 4.0
    density_code = 0.75
    density_cgs = density_code*mass_cgs/length_cgs**3
    pressure_cgs = mass_cgs/(length_cgs*time_cgs**2)
    edot_cgs_to_code = time_cgs/pressure_cgs
    if density_kind == "number_density":
        mu = 0.5
        q = density_cgs/(mu*AMU_CGS)
    else:
        hydrogen_mass_fraction = 0.4
        q = density_cgs*hydrogen_mass_fraction/AMU_CGS
    lambda_cgs = expected_rate/(q*q*edot_cgs_to_code)
    data = _run_case(
        basename,
        ["problem/density=0.75",
         "problem/pressure=1.2",
         "units/length_cgs=2.0",
         "units/mass_cgs=16.0",
         "units/time_cgs=4.0",
         "cooling/cooling_model=powerlaw",
         "cooling/shielding_model=none",
         f"cooling/cooling_density={density_kind}",
         "cooling/cooling_powerlaw_value_units=cgs",
         f"cooling/cooling_reference_value={lambda_cgs:.17e}",
         *composition_flags],
        input_file="inputs/cooling_test_hydro_composition.athinput",
    )
    _assert_history_rate(data, gross=expected_rate, net=expected_rate,
                         rtol=2.0e-12, atol=1.0e-13)


def test_history_can_be_disabled():
    data = _run_case("cooling_history_off", ["cooling/history=false"])
    assert "cool_gross" not in data
    assert "cool_net" not in data


@pytest.mark.parametrize(
    "basename,flags,expected_gross",
    [
        (
            "cooling_bounds_zero",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_1d.tbl",
             "cooling/cooling_bounds=zero",
             "problem/pressure=3.0"],
            0.0,
        ),
        (
            "cooling_bounds_clamp",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_1d.tbl",
             "cooling/cooling_bounds=clamp",
             "problem/pressure=3.0"],
            2.0e-2,
        ),
    ],
)
def test_table_bounds_zero_and_clamp(basename, flags, expected_gross):
    data = _run_case(basename, flags)
    _assert_history_rate(data, gross=expected_gross, net=expected_gross)


def test_cooling_timestep_constraint_and_bounds():
    included = _run_case(
        "cooling_dt_included",
        ["time/tlim=1.0",
         "cooling/timestep=true",
         "cooling/timestep_factor=0.1",
         "cooling/cooling_reference_value=1.0e2",
         "cooling/timestep_use_temperature_bounds=true",
         "cooling/timestep_temperature_units=code",
         "cooling/timestep_temperature_min=0.5",
         "cooling/timestep_temperature_max=1.5"],
    )
    row = _first_evolved_row(included)
    assert np.isclose(included["time"][row], 6.0e-4, rtol=1.0e-12,
                      atol=1.0e-14)

    excluded = _run_case(
        "cooling_dt_excluded",
        ["time/tlim=1.0",
         "cooling/timestep=true",
         "cooling/timestep_factor=0.1",
         "cooling/cooling_reference_value=1.0e2",
         "cooling/timestep_use_temperature_bounds=true",
         "cooling/timestep_temperature_units=code",
         "cooling/timestep_temperature_min=2.0",
         "cooling/timestep_temperature_max=3.0"],
    )
    row = _first_evolved_row(excluded)
    assert excluded["time"][row] > 1.0e-2


GM1 = 2.0/3.0
DT = 1.0e-3  # the uniform test takes one step to tlim


def _step_energies(data):
    row = _first_evolved_row(data)
    return data["tot-E"][row - 1], data["tot-E"][row]


def test_exact_integrator_solves_stiff_cooling():
    """Lambda proportional to T: e decays exponentially, exact at any dt/t_cool."""
    lam0 = 1.0e4  # dt/t_cool ~ 6.7, where explicit RK2 is unstable
    data = _run_case(
        "cooling_exact_linear",
        ["cooling/integrator=exact",
         "cooling/cooling_slope=1.0",
         f"cooling/cooling_reference_value={lam0}"],
    )
    e0, e1 = _step_energies(data)
    assert np.isclose(e1, e0*np.exp(-lam0*GM1*DT), rtol=1.0e-12, atol=0.0)
    _assert_history_rate(data, gross=(e0 - e1)/DT, net=(e0 - e1)/DT, rtol=1.0e-12)


def test_exact_integrator_relaxes_to_heating_equilibrium():
    """Linear cooling against constant heating relaxes to T_eq = heating/Lambda_0."""
    lam0, heat = 1.0e4, 5.0e3
    data = _run_case(
        "cooling_exact_equilibrium",
        ["cooling/integrator=exact",
         "cooling/cooling_slope=1.0",
         f"cooling/cooling_reference_value={lam0}",
         "cooling/heating_model=constant",
         f"cooling/heating_rate={heat}"],
    )
    e0, e1 = _step_energies(data)
    t0, t_eq, x = GM1*e0, heat/lam0, lam0*GM1*DT
    assert np.isclose(e1, (t_eq + (t0 - t_eq)*np.exp(-x))/GM1, rtol=1.0e-12, atol=0.0)
    gross = lam0*(t_eq + (t0 - t_eq)*(1.0 - np.exp(-x))/x)  # time-averaged cooling
    _assert_history_rate(data, gross=gross, net=gross - heat, rtol=1.0e-12)


def test_exact_integrator_rests_on_table_edge():
    """Lambda = 1e4*T on 1 <= T <= 2 and zero below (bounds=zero), against heating:
    gas cooling onto the T = 1 edge stays there, with cooling balancing heating."""
    heat = 5.0e3
    data = _run_case(
        "cooling_exact_table_edge",
        ["cooling/integrator=exact",
         "cooling/cooling_model=table",
         "cooling/cooling_table=inputs/tables/cooling_lambda_1d_edge.tbl",
         "cooling/cooling_bounds=zero",
         "cooling/heating_model=constant",
         f"cooling/heating_rate={heat}",
         "problem/pressure=1.8"],
    )
    e0, e1 = _step_energies(data)
    assert np.isclose(e1, 1.0/GM1, rtol=1.0e-12, atol=0.0)
    net = (e0 - e1)/DT
    _assert_history_rate(data, gross=net + heat, net=net, rtol=1.0e-12)


def test_exact_integrator_mhd_path():
    lam0 = 1.0e4
    data = _run_case(
        "cooling_exact_mhd",
        ["cooling/integrator=exact",
         "cooling/cooling_slope=1.0",
         f"cooling/cooling_reference_value={lam0}"],
        input_file=MHD_INPUT, physics="mhd")
    e0, e1 = _step_energies(data)
    emag = 0.5*(0.1**2 + 0.2**2 + 0.3**2)
    assert np.isclose(e1 - emag, (e0 - emag)*np.exp(-lam0*GM1*DT), rtol=1.0e-11,
                      atol=0.0)


def _run_failure(basename, flags, input_file=HYDRO_INPUT):
    _remove_outputs(basename, "hydro")
    command = ["./athena", "-i", input_file, f"job/basename={basename}"] + flags
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, check=False)
    assert result.returncode != 0
    return result.stdout + result.stderr


@pytest.mark.parametrize(
    "basename,flags,expected_text",
    [
        (
            "cooling_legacy_rejected",
            ["hydro/ism_cooling=true"],
            "Use the standalone <cooling> block instead",
        ),
        (
            "cooling_integrator_unknown",
            ["cooling/integrator=implicit"],
            "integrator must be 'explicit' or 'exact'",
        ),
        (
            "cooling_exact_rejects_user_model",
            ["cooling/integrator=exact", "cooling/cooling_model=user"],
            "integrator=exact requires built-in cooling and heating models",
        ),
        (
            "cooling_table_bounds_fatal",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/cooling_lambda_1d.tbl",
             "cooling/cooling_bounds=fatal",
             "problem/pressure=3.0"],
            "fatal cooling/heating bounds were exceeded",
        ),
        (
            "cooling_cgs_number_density_requires_mu",
            ["cooling/cooling_density=number_density",
             "cooling/cooling_powerlaw_axis=density",
             "cooling/cooling_powerlaw_axis_units=code",
             "cooling/cooling_powerlaw_density_kind=mass_density",
             "cooling/cooling_powerlaw_value_units=cgs"],
            "requires <cooling>/mu",
        ),
        (
            "cooling_cgs_hydrogen_density_requires_x",
            ["cooling/cooling_density=hydrogen_number_density",
             "cooling/cooling_powerlaw_axis=density",
             "cooling/cooling_powerlaw_axis_units=code",
             "cooling/cooling_powerlaw_density_kind=mass_density",
             "cooling/cooling_powerlaw_value_units=cgs"],
            "requires <cooling>/hydrogen_mass_fraction or <cooling>/mu_H",
        ),
        (
            "cooling_cgs_dt_temperature_requires_mu",
            ["cooling/timestep=true",
             "cooling/timestep_use_temperature_bounds=true",
             "cooling/timestep_temperature_units=cgs",
             "cooling/timestep_temperature_min=1.0e-20",
             "cooling/timestep_temperature_max=1.0e20"],
            "cgs temperature evaluation requires <cooling>/temperature_mu",
        ),
        (
            "cooling_missing_axis_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/cooling_missing_axis.tbl"],
            "missing axis1",
        ),
        (
            "cooling_wrong_value_count_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/"
             "cooling_wrong_value_count.tbl"],
            "values but expected",
        ),
        (
            "cooling_bad_value_kind_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/"
             "cooling_bad_value_kind.tbl"],
            "this use requires value_kind lambda",
        ),
        (
            "cooling_nonfinite_value_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/"
             "cooling_nonfinite_value.tbl"],
            "contains a non-finite value",
        ),
        (
            "cooling_bad_log_axis_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/cooling_bad_log_axis.tbl"],
            "log10 axis0 extent",
        ),
        (
            "cooling_bad_axis_scale_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/"
             "cooling_bad_axis_scale.tbl"],
            "optional axis scale must be 'linear' or 'log10'",
        ),
        (
            "cooling_scalar_index_rejected",
            ["cooling/cooling_model=table",
             "cooling/cooling_table=inputs/tables/malformed/"
             "cooling_scalar_bad_index.tbl"],
            "scalar_index=1",
        ),
    ],
)
def test_invalid_cooling_inputs_fail_cleanly(basename, flags, expected_text):
    input_file = "inputs/cooling_test_legacy_hydro.athinput" if (
        basename == "cooling_legacy_rejected") else HYDRO_INPUT
    output = _run_failure(basename, flags, input_file=input_file)
    assert expected_text in output


def test_table_validator_accepts_runtime_table_and_plots(tmp_path):
    table = ROOT / "tst" / "inputs" / "tables" / "cooling_lambda_3d.tbl"
    plot = tmp_path / "table.png"
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(table),
         "--expect-value-kind", "lambda", "--plot", str(plot)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "OK:" in result.stdout
    assert plot.exists()


def test_table_validator_rejects_malformed_table():
    table = ROOT / "tst" / "inputs" / "tables" / "malformed" / (
        "cooling_wrong_value_count.tbl")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(table),
         "--expect-value-kind", "lambda"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "values but expected" in result.stderr


def test_table_validator_can_generate_example(tmp_path):
    table = tmp_path / "generated.tbl"
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "--write-example", str(table),
         "--ndim", "2", "--value-kind", "lambda"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert table.exists()
    assert "ndim=2" in result.stdout
