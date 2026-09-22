"""Turbulence history regressions; build with PROBLEM=built_in_pgens."""

import numpy as np
import pytest

from .test_linear_drag_cpu import run_athena
from .test_turb_sgs_2d_cpu import (
    INPUT, assemble_2d_blocks, latest, require_success,
    run_athena as run_turbulence,
)

from athena_read import hst
from bin_convert import read_binary


BUDGET_NAMES = {
    "enstrophy", "drag-KE", "drag-enst", "visc-KE", "visc-enst",
    "force-KE", "force-enst", "p-dilat", "enst-comp",
}


@pytest.mark.parametrize("eos, nz, axis", [
    ("ideal", 1, 0), ("isothermal", 1, 1), ("ideal", 4, 2),
])
def test_shear_history_rates_and_block_decomposition(tmp_path, eos, nz, axis):
    """A transverse sine wave has known drag and discrete viscous loss rates."""
    # Non-unit density and volume distinguish integrals from volume/mass averages.
    rho, volume, amplitude, speed, nu, alpha = 2.0, 6.0, 0.2, 0.125, 0.03, 2.0
    length = (2.0, 3.0, 1.0)[axis]
    dx = length / (32, 8, nz)[axis]
    k = 2 * np.pi / length
    enstrophy = volume * amplitude**2 * (np.sin(k * dx) / dx)**2 / 4
    kinetic = 0.5 * rho * volume * (speed**2 + amplitude**2 / 2)
    visc_eigenvalue = -4 * nu * np.sin(k * dx / 2)**2 / dx**2
    expected = {
        "mass": rho * volume,
        "enstrophy": enstrophy,
        "drag-KE": 2 * alpha * kinetic,
        "drag-enst": 2 * alpha * enstrophy,
        "visc-KE": -visc_eigenvalue * rho * volume * amplitude**2 / 2,
        "visc-enst": 2 * visc_eigenvalue * enstrophy,
        "force-KE": 0.0,
        "force-enst": 0.0,
        "p-dilat": 0.0,
        "enst-comp": 0.0,
    }
    evolve = eos == "ideal" and nz == 1
    for bx, by in [(32, 8), (8, 4)]:
        output_dir = tmp_path / f"blocks_{bx}_{by}"
        require_success(run_athena(
            output_dir, f"time/nlim={4 if evolve else 0}", "mesh/nx1=32", "mesh/nx2=8",
            "mesh/x1min=0", "mesh/x1max=2", "mesh/x2min=0", "mesh/x2max=3",
            f"mesh/nx3={nz}", f"meshblock/nx3={nz}",
            f"meshblock/nx1={bx}", f"meshblock/nx2={by}",
            f"hydro/eos={eos}", "hydro/iso_sound_speed=1.0", "hydro/rsolver=llf",
            f"hydro/viscosity={nu}", "problem/pgen_name=linear_wave",
            f"problem/wave_flag={2 if eos == 'ideal' else 1}",
            f"problem/along_x{axis + 1}=true", f"problem/dens={rho}", "problem/pgas=1.0",
            f"problem/amp={amplitude}", f"problem/vx0={speed}",
            "output2/file_type=hst", "output2/dcycle=1",
            "output2/turbulence=true", "output2/data_format=%23.16e",
            "output3/file_type=bin", "output3/variable=hydro_wz", "output3/id=wz",
            "output3/dcycle=1", "output4/file_type=bin", "output4/variable=hydro_w2",
            "output4/id=w2", "output4/dcycle=1",
        ))
        history = hst(output_dir / "linear_drag_test.hydro.hst")
        assert len(history["time"]) == (5 if evolve else 1)
        for name, value in expected.items():
            np.testing.assert_allclose(history[name][0], value, rtol=2e-13, atol=1e-13,
                                       err_msg=name)
        if evolve:
            decay = np.exp(-2 * (alpha - visc_eigenvalue) * history["time"])
            np.testing.assert_allclose(history["enstrophy"], enstrophy * decay, rtol=0.01)
            np.testing.assert_allclose(history["drag-enst"], 2 * alpha * history["enstrophy"])
        for name, field, power in (("wz", "vorz", 1), ("w2", "vor2", 2)):
            data = read_binary(str(sorted((output_dir / "bin").glob(f"*.{name}.*.bin"))[0]))
            for values, geometry in zip(data["mb_data"][field], data["mb_geometry"]):
                lo, hi = geometry[2 * axis:2 * axis + 2]
                count = values.shape[2 - axis]
                position = lo + (np.arange(count) + 0.5) * (hi - lo) / count
                omega = amplitude * np.sin(k * dx) / dx * np.cos(k * position)
                expected_field = omega**power if name == "w2" or axis < 2 else 0 * omega
                shape = [1, 1, 1]
                shape[2 - axis] = count
                np.testing.assert_allclose(values, np.broadcast_to(
                    expected_field.reshape(shape), values.shape), rtol=1e-6, atol=1e-8)


def test_driven_variable_density_history_matches_instantaneous_fields(tmp_path):
    """Forcing and face-stress rates agree with the resolved variable-density state."""
    input_file = tmp_path / "history.athinput"
    input_file.write_text(INPUT.read_text() +
                         "\n<output1>\nturbulence=true\ndata_format=%23.16e\n"
                         "\n<output2>\nvariable=hydro_w\n"
                         "\n<hydro>\neos=ideal\ngamma=1.4\nrsolver=hllc\n"
                         "\n<problem>\npgen_name=linear_wave\nwave_flag=1\n"
                         "amp=0.3\ndens=1\npgas=1\nalong_x1=true\nvx0=0.125\n")
    require_success(run_turbulence(tmp_path, input_file=input_file))
    state = read_binary(str(latest(tmp_path / "bin", "*.state.*.bin")))
    force = read_binary(str(latest(tmp_path / "bin", "*.force.*.bin")))
    fields = assemble_2d_blocks(state, ("dens", "velx", "vely", "eint"))
    acceleration = assemble_2d_blocks(force, ("force1", "force2"))
    rho, ux, uy = [fields[name] for name in ("dens", "velx", "vely")]
    assert np.ptp(rho) > 0.25
    ax, ay = [acceleration[name] for name in ("force1", "force2")]
    dx, dy = 1 / state["Nx1"], 1 / state["Nx2"]

    def derivative(field, axis, spacing):
        return (np.roll(field, -1, axis) - np.roll(field, 1, axis)) / (2 * spacing)

    omega = derivative(uy, 1, dx) - derivative(ux, 0, dy)
    divu = derivative(ux, 1, dx) + derivative(uy, 0, dy)
    velocity = (ux, uy)
    grad = [[derivative(u, axis, h) for axis, h in ((1, dx), (0, dy))]
            for u in velocity]
    visc_rhs = np.zeros((2, *rho.shape))
    for d, (axis, h) in enumerate(((1, dx), (0, dy))):
        face_grad = [[(np.roll(u, -1, axis) - u) / h if j == d else
                      0.5 * (grad[n][j] + np.roll(grad[n][j], -1, axis))
                      for j in range(2)] for n, u in enumerate(velocity)]
        div_face = face_grad[0][0] + face_grad[1][1]
        mu_face = 0.01 * 0.5 * (rho + np.roll(rho, -1, axis))
        for n in range(2):
            stress = mu_face * (face_grad[n][d] + face_grad[d][n]
                                - (2 / 3 * div_face if n == d else 0))
            visc_rhs[n] += (stress - np.roll(stress, 1, axis)) / h
    expected = {
        "enstrophy": 0.5 * np.sum(omega**2) * dx * dy,
        "visc-KE": -np.sum(ux * visc_rhs[0] + uy * visc_rhs[1]) * dx * dy,
        "visc-enst": np.sum(visc_rhs[0] / rho * derivative(omega, 0, dy)
                             - visc_rhs[1] / rho * derivative(omega, 1, dx)) * dx * dy,
        "force-KE": np.sum(rho * (ux * ax + uy * ay)) * dx * dy,
        "force-enst": np.sum(omega * (derivative(ay, 1, dx)
                                      - derivative(ax, 0, dy))) * dx * dy,
        "p-dilat": np.sum(0.4 * fields["eint"] * divu) * dx * dy,
        "enst-comp": -0.5 * np.sum(omega**2 * divu) * dx * dy,
    }
    history = hst(tmp_path / "turb_sgs_2d_test.hydro.hst")
    assert history["time"][-1] == pytest.approx(state["time"])
    assert history["force-KE"][-1] > 0
    for name, value in expected.items():
        # Fine output is float32, while history reductions use simulation precision.
        np.testing.assert_allclose(history[name][-1], value, rtol=2e-5, atol=2e-10,
                                   err_msg=name)
    np.testing.assert_array_equal(history["drag-KE"], 0.0)
    np.testing.assert_array_equal(history["drag-enst"], 0.0)


def test_turbulence_history_is_opt_in(tmp_path):
    require_success(run_athena(tmp_path, "time/nlim=0", "output2/file_type=hst",
                              "output2/dcycle=1"))
    history = hst(tmp_path / "linear_drag_test.hydro.hst")
    assert BUDGET_NAMES.isdisjoint(history)
    assert len(history) == 10  # time, dt, and the eight standard ideal-hydro columns.


@pytest.mark.parametrize("overrides, message", [
    (("mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow"),
     "requires a periodic, uniform Newtonian fluid grid"),
    (("hydro/nscalars=4",), "exceeds NHISTORY_VARIABLES"),
])
def test_turbulence_history_rejects_unsupported_configuration(tmp_path, overrides, message):
    result = run_athena(tmp_path, "time/nlim=0", "output2/file_type=hst",
                       "output2/dcycle=1", "output2/turbulence=true", *overrides)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
