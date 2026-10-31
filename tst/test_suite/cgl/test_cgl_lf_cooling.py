"""History-sampled CGL cooling, floors, and native checkpoint continuation.

Run with ATHENAK_TEST_EXE pointing to a double-precision built_in_pgens executable.
"""

import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


ROOT = Path(__file__).resolve().parents[3]
TEMPLATE = ROOT / "inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"


def input_text(*, ppar=2.0, pperp=2.2, target=0.5, tau=0.1, floor=1e-12,
               integrator="rk2", history="dcycle = 100", nlim=6):
    text = TEMPLATE.read_text().split("<output1>", 1)[0]
    text = text.replace("cgl_heat_flux = landau_fluid\n", "")
    text = text.replace("user_hist = true\n", "")
    return text + f"""
<job>
basename = cooling
<mesh>
nghost = 4
nx1 = 8
nx2 = 4
nx3 = 4
x3max = 1
<meshblock>
nx1 = 4
nx2 = 4
nx3 = 4
<time>
integrator = {integrator}
sts_integrator = none
nlim = {nlim}
tlim = 10
<mhd>
mirror_limiter = false
firehose_limiter = false
pfloor = {floor}
tfloor = 0
sfloor = 0
<problem>
pgen_name = cgl_lf_cooling
p_parallel0 = {ppar}
p_perp0 = {pperp}
cooling_piso = {target}
cooling_tau = {tau}
<turb_driving>
dedt = 0
<output1>
file_type = hst
data_format = %25.17e
{history}
<output2>
file_type = rst
dcycle = 1
"""


def run(directory, source=None, *, restart=None, flags=(), ranks=1, adiabatic=False):
    directory.mkdir()
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    command = [str(executable), "-d", str(directory)]
    if source is not None:
        deck = directory / "input.athinput"
        deck.write_text(source)
        command += ["-i", str(deck)]
    if restart is not None:
        command += ["-r", str(restart)]
    if ranks > 1:
        command = ["mpirun", "-np", str(ranks), *command]
    result = subprocess.run(command + list(flags), cwd=directory, capture_output=True,
                            text=True, check=False, timeout=120)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output)
    assert result.returncode == 0, output
    # Final output can repeat the last cycle; retain the final version of each state.
    states = {}
    for path in sorted((directory / "rst").glob("*.rst")):
        state = passive_acceptance.decode_restart(path, adiabatic,
                                                  active_cgl=not adiabatic,
                                                  adiabatic=adiabatic)
        header = path.read_bytes().split(b"<par_end>\n", 1)[0].decode()
        for name in ("cooling_alpha", "cooling_history_last_time",
                     "cooling_history_last_cycle", "cooling_target_fraction",
                     "cooling_sample_time"):
            match = re.search(rf"^{name}\s*=\s*(\S+)", header, re.MULTILINE)
            assert match, name
            state[name] = float(match[1])
        state["path"] = path
        states[state["cycle"]] = state
    assert states
    return states


def pressures(state):
    if "pp" in state:
        return state["pp"], state["pt"]
    u, b = state["active"], state["bcc"]
    rho = u[:, 0]
    bsqr = np.sum(b*b, axis=1)
    thermal = u[:, 4] - 0.5*np.sum(u[:, 1:4]**2, axis=1)/rho - 0.5*bsqr
    ratio = np.exp(u[:, 5]/rho - 2*np.log(rho) + 1.5*np.log(bsqr))
    ppar = thermal/(0.5 + ratio)
    return ppar, ratio*ppar


@pytest.mark.parametrize("adiabatic", [False, True])
def test_frozen_cooling_checkpoint_preserves_controller(tmp_path, adiabatic):
    source = input_text(history="dcycle = 1", nlim=1)
    if adiabatic:
        source += adiabatic_settings()
    seed = run(tmp_path / "seed", source, adiabatic=adiabatic)[1]
    overlay = (ROOT / "inputs/particles/frozen_cgl.athinput").read_text()
    final_time = seed["time"] + 0.02
    overlay += ("\n<particles>\nppc = 0.015625\n"
                f"<time>\ntlim = {final_time:.17g}\n")
    frozen = run(tmp_path / "frozen", overlay, restart=seed["path"],
                 adiabatic=adiabatic)
    final = frozen[max(frozen)]
    assert final["time"] == final_time
    assert final["cycle"] > seed["cycle"]
    np.testing.assert_array_equal(final["active"], seed["active"])
    for first, second in zip(final["faces"], seed["faces"]):
        np.testing.assert_array_equal(first, second)
    for name in ("cooling_alpha", "cooling_history_last_time",
                 "cooling_history_last_cycle", "cooling_target_fraction",
                 "cooling_sample_time"):
        assert final[name] == seed[name]
    hist = passive_acceptance.history(tmp_path / "frozen/cooling.user.hst")
    np.testing.assert_array_equal(hist["cool_alpha"], seed["cooling_alpha"])
    checkpoint = sorted((tmp_path / "frozen/prst/rank_00000000").glob("*.prst"))[-1]
    particles = np.fromfile(checkpoint, dtype=np.float64)
    assert particles[0] == final_time and particles[2] == 2
    resumed_overlay = overlay + (
        "\n<particles>\nseed_on_restart = false\n"
        "<problem>\nprtcl_rst_flag = 1\n"
        f"prtcl_res_file = {checkpoint}\n"
    )
    resumed = run(tmp_path / "resumed", resumed_overlay, restart=final["path"],
                  flags=(f"time/tlim={final_time + 0.01:.17g}",), adiabatic=adiabatic)
    last = resumed[max(resumed)]
    np.testing.assert_array_equal(last["active"], final["active"])
    assert last["cooling_alpha"] == seed["cooling_alpha"]


def assert_background_unchanged(initial, final):
    np.testing.assert_allclose(final["active"][:, :4], initial["active"][:, :4],
                               rtol=0, atol=3e-13)
    for before, after in zip(initial["faces"], final["faces"]):
        np.testing.assert_allclose(after, before, rtol=0, atol=3e-13)


def adiabatic_settings(*, pgas=2.3, gamma=5/3):
    return f"""
<mhd>
passive = true
passive_eos = adiabatic
gamma = {gamma:.17g}
<problem>
passive_delta = true
p_gas0 = {pgas}
"""


def gas_pressure(state, gamma=5/3):
    u, b = state["active"], state["bcc"]
    thermal = u[:, 6] - 0.5*np.sum(u[:, 1:4]**2, axis=1)/u[:, 0] - 0.5*np.sum(b*b, axis=1)
    return (gamma - 1)*thermal


@pytest.mark.parametrize("gamma", [5/3, 1.4])
@pytest.mark.parametrize("integrator", ["rk2", "rk4"])
def test_adiabatic_cooling_uses_gas_controller(tmp_path, gamma, integrator):
    source = input_text(integrator=integrator) + adiabatic_settings(gamma=gamma)
    states = run(tmp_path / "gas", source, adiabatic=True)
    alpha = (1 - 0.5/2.3)/0.1
    piso0 = (2.0 + 2*2.2)/3
    for state in states.values():
        factor = 0.5/2.3 + (1 - 0.5/2.3)*np.exp(-state["time"]/0.1)
        pp, pt = pressures(state)
        np.testing.assert_allclose(gas_pressure(state, gamma), 2.3*factor,
                                   rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose((pp + 2*pt)/3, piso0*factor, rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose(pt - pp, 0.2, rtol=0, atol=3e-12)
        assert state["cooling_alpha"] == pytest.approx(alpha, rel=3e-12)
        assert_background_unchanged(states[0], state)


@pytest.mark.parametrize("floor_kind", ["pfloor", "tfloor", "sfloor", "roundoff"])
def test_adiabatic_cooling_gas_floors(tmp_path, floor_kind):
    gamma, rho = 5/3, 2.0
    floor = 1e-20 if floor_kind == "roundoff" else 0.5
    source = (input_text(target=1e-25, tau=1e-5, floor=1e-20, nlim=3)
              + adiabatic_settings() + f"\n<problem>\nrho0 = {rho}\n")
    key = "pfloor" if floor_kind == "roundoff" else floor_kind
    source += f"\n<mhd>\n{key} = {floor}\n"
    states = run(tmp_path / floor_kind, source, adiabatic=True)
    expected = floor * {"pfloor": 1, "tfloor": rho, "sfloor": rho**gamma,
                        "roundoff": 1}[floor_kind]
    for cycle, state in states.items():
        if cycle == 0:
            continue
        pgas = gas_pressure(state)
        assert np.min(pgas) >= expected*(1 - 1e-14)
        np.testing.assert_allclose(pgas, expected, rtol=3e-12, atol=3e-15)
        assert np.all(state["active"][:, 6] <= states[0]["active"][:, 6])
        assert_background_unchanged(states[0], state)


def test_passive_cooling_cap_does_not_limit_gas_cooling(tmp_path):
    source = (input_text(ppar=0.8, pperp=1.0, target=0.05, tau=1e-5,
                         floor=1e-6, nlim=3) + adiabatic_settings())
    states = run(tmp_path / "cap", source, adiabatic=True)
    for cycle, state in states.items():
        if cycle == 0:
            continue
        pp, pt = pressures(state)
        np.testing.assert_allclose(pt/pp, 10, rtol=3e-12, atol=0)
        np.testing.assert_allclose(pt - pp, 0.2, rtol=0, atol=3e-13)
        np.testing.assert_allclose(gas_pressure(state), 0.05, rtol=1e-9, atol=0)
        assert_background_unchanged(states[0], state)


@pytest.mark.parametrize("integrator", ["rk2", "rk4"])
@pytest.mark.parametrize("history", ["dcycle = 100", "dt = 0.005"])
def test_cooling_exact_sampled_interval(tmp_path, integrator, history):
    source = input_text(integrator=integrator, history=history)
    states = run(tmp_path / integrator, source)
    initial = states[0]
    piso0, delta = (2.0 + 2*2.2)/3, 0.2
    alpha = (1 - 0.5/piso0)/0.1
    for state in states.values():
        ppar, pperp = pressures(state)
        piso = 0.5 + (piso0 - 0.5)*np.exp(-state["time"]/0.1)
        np.testing.assert_allclose(ppar, piso - 2*delta/3, rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose(pperp, piso + delta/3, rtol=3e-12, atol=3e-13)
        if history.startswith("dcycle"):
            assert state["cooling_alpha"] == pytest.approx(alpha, rel=2e-13)
            assert state["cooling_history_last_cycle"] == 0
        else:
            assert state["cooling_sample_time"] == state["time"]
        assert_background_unchanged(initial, state)
    assert states[max(states)]["time"] > 0


@pytest.mark.parametrize("landau_fluid", [False, True])
def test_cooling_and_collisions_keep_their_independent_rates(tmp_path, landau_fluid):
    source = input_text() + "\n<mhd>\nnu_coll = 0.7\n"
    if landau_fluid:
        source += ("cgl_heat_flux = landau_fluid\ncgl_heat_flux_integrator = sts\n"
                   "<time>\nsts_integrator = rkl2\n")
    states = run(tmp_path / "collisions", source)
    piso0 = (2.0 + 2*2.2)/3
    for state in states.values():
        parallel, perpendicular = pressures(state)
        piso = 0.5 + (piso0 - 0.5)*np.exp(-state["time"]/0.1)
        delta = 0.2*np.exp(-0.7*state["time"])
        np.testing.assert_allclose((parallel + 2*perpendicular)/3, piso,
                                   rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose(perpendicular - parallel, delta,
                                   rtol=3e-12, atol=3e-13)
        assert_background_unchanged(states[0], state)


def test_cooling_below_target_is_off(tmp_path):
    states = run(tmp_path / "off", input_text(target=3, history="dcycle = 2"))
    for state in states.values():
        assert state["cooling_alpha"] == 0
        np.testing.assert_allclose(state["payload"], states[0]["payload"],
                                   rtol=0, atol=3e-13)


def test_cooling_coefficient_changes_only_at_history(tmp_path):
    states = run(tmp_path / "adaptive", input_text(history="dcycle = 2"))
    previous = states[0]
    coefficients = []
    for cycle, state in states.items():
        ppar, pperp = pressures(state)
        piso = (ppar + 2*pperp)/3
        if cycle:
            old_par, old_perp = pressures(previous)
            expected = 0.5 + ((old_par + 2*old_perp)/3 - 0.5)*np.exp(
                -(state["time"] - previous["time"])/0.1)
            np.testing.assert_allclose(piso, expected, rtol=3e-12, atol=3e-13)
        expected_alpha = ((1 - 0.5/np.mean(piso))/0.1 if cycle % 2 == 0
                          else previous["cooling_alpha"])
        assert state["cooling_alpha"] == pytest.approx(expected_alpha, rel=3e-12)
        assert state["cooling_history_last_cycle"] == cycle - cycle % 2
        np.testing.assert_allclose(pperp - ppar, 0.2, rtol=0, atol=3e-12)
        coefficients.append(state["cooling_alpha"])
        previous = state
    assert len(set(coefficients)) >= 3
    history = passive_acceptance.history(tmp_path / "adaptive/cooling.user.hst")
    assert "cool_alpha" in history
    # Scheduled history and finalization at cycle 6 must not update twice.
    assert history["cool_alpha"][-1] == history["cool_alpha"][-2]


@pytest.mark.parametrize("ppar,pperp", [(0.8, 1.0), (1.0, 0.8)])
def test_cooling_floor_preserves_anisotropy(tmp_path, ppar, pperp):
    states = run(tmp_path / "floor", input_text(ppar=ppar, pperp=pperp,
                 target=0.05, tau=1e-5, floor=0.5, nlim=3))
    for cycle, state in states.items():
        if cycle == 0:
            continue
        parallel, perpendicular = pressures(state)
        np.testing.assert_allclose(np.minimum(parallel, perpendicular), 0.5,
                                   rtol=0, atol=3e-13)
        np.testing.assert_allclose(perpendicular - parallel, pperp - ppar,
                                   rtol=0, atol=3e-13)
        assert_background_unchanged(states[0], state)


@pytest.mark.parametrize("landau_fluid", [False, True])
@pytest.mark.parametrize("ppar,pperp", [(0.8, 1.0), (1.0, 0.8)])
def test_cooling_stops_at_pressure_ratio_cap(tmp_path, ppar, pperp, landau_fluid):
    source = input_text(ppar=ppar, pperp=pperp, target=0.05, tau=1e-5,
                        floor=1e-6, nlim=3)
    if landau_fluid:
        source += ("\n<mhd>\ncgl_heat_flux = landau_fluid\n"
                   "cgl_heat_flux_integrator = sts\n<time>\nsts_integrator = rkl2\n")
    states = run(tmp_path / "guarded", source)
    for cycle, state in states.items():
        if cycle == 0:
            continue
        parallel, perpendicular = pressures(state)
        smaller = np.minimum(parallel, perpendicular)
        larger = np.maximum(parallel, perpendicular)
        np.testing.assert_allclose(larger/smaller, 10, rtol=3e-12, atol=0)
        np.testing.assert_allclose(perpendicular - parallel, pperp - ppar,
                                   rtol=0, atol=3e-13)
        removed = states[0]["active"][:, 4] - state["active"][:, 4]
        np.testing.assert_allclose(removed, 1.5*(min(ppar, pperp) - smaller),
                                   rtol=0, atol=3e-13)
        assert_background_unchanged(states[0], state)
    resumed = run(tmp_path / "resumed", restart=states[1]["path"],
                  flags=("time/nlim=3",))
    np.testing.assert_array_equal(resumed[3]["payload"], states[3]["payload"])


@pytest.mark.parametrize("ppar,pperp", [(1e-6, 0.392), (0.392, 1e-6)])
def test_cooling_skips_existing_extreme_ratio(tmp_path, ppar, pperp):
    states = run(tmp_path / "extreme", input_text(ppar=ppar, pperp=pperp,
                 target=0.05, tau=1e-5, floor=1e-8, nlim=3))
    for state in states.values():
        assert state["cooling_alpha"] > 0
        np.testing.assert_allclose(state["payload"], states[0]["payload"],
                                   rtol=0, atol=3e-13)


def test_cooling_cold_floor_retains_representable_pressure(tmp_path):
    floor = 1e-20
    states = run(tmp_path / "cold", input_text(ppar=2, pperp=2, target=1e-25,
                 tau=1e-5, floor=floor, nlim=3))
    for cycle, state in states.items():
        if cycle == 0:
            continue
        parallel, perpendicular = pressures(state)
        # Magnetic energy can set a representable thermal minimum above pfloor.
        assert np.all(parallel >= floor) and np.all(perpendicular >= floor)
        assert np.max(parallel) < 1e-10
        np.testing.assert_allclose(perpendicular - parallel, 0, rtol=0, atol=5e-15)
        assert_background_unchanged(states[0], state)


@pytest.mark.parametrize("ranks", [1, 2])
@pytest.mark.parametrize("case,expected", [
    ("ratio", 1), ("blocked", 1), ("floor", 0), ("free", 0),
])
def test_cooling_guard_counts_owned_cell_updates(tmp_path, case, expected, ranks):
    if ranks == 2:
        executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
        config = subprocess.run([str(executable), "-c"], capture_output=True,
                                text=True, check=True).stdout
        if not re.search(r"MPI parallelism:\s+ON", config) or not shutil.which("mpirun"):
            pytest.skip("requires an MPI executable and mpirun")
    pp, pt = (1e-6, 0.392) if case == "blocked" else (0.8, 1.0)
    source = input_text(ppar=pp, pperp=pt, target=0.05,
                        tau=100 if case == "free" else 1e-5,
                        floor=0.5 if case == "floor" else 1e-8,
                        history="dcycle = 1", nlim=3)
    states = run(tmp_path / case, source, ranks=ranks)
    history = passive_acceptance.history(tmp_path / case / "cooling.user.hst")
    np.testing.assert_array_equal(history["cool_hits"][1:], 128*expected)
    np.testing.assert_array_equal(history["cool_frac"][1:], expected)
    assert history["cool_hits"][0] == history["cool_frac"][0] == 0
    resumed = run(tmp_path / "resumed", restart=states[2]["path"],
                  flags=("time/nlim=3",), ranks=ranks)
    hist2 = passive_acceptance.history(tmp_path / "resumed/cooling.user.hst")
    np.testing.assert_array_equal(hist2["cool_hits"], 128*expected)
    np.testing.assert_array_equal(hist2["cool_frac"], expected)
    np.testing.assert_array_equal(resumed[3]["payload"], states[3]["payload"])


@pytest.mark.parametrize("checkpoint_kind", ["interior", "off_cadence_final"])
@pytest.mark.parametrize("adiabatic", [False, True])
def test_cooling_restart_preserves_history_schedule(tmp_path, checkpoint_kind, adiabatic):
    source = input_text(history="dt = 0.04", nlim=12)
    if adiabatic:
        source += adiabatic_settings()
    direct = run(tmp_path / "direct", source, adiabatic=adiabatic)
    # Cycle 1 is before the first scheduled update, with a different current pressure.
    assert direct[1]["time"] < 0.04 and direct[1]["cooling_history_last_cycle"] == 0
    if checkpoint_kind == "interior":
        checkpoint = direct[1]
    else:
        checkpoint = run(tmp_path / "partial", source, flags=("time/nlim=1",),
                         adiabatic=adiabatic)[1]
    resumed = run(tmp_path / "resumed", restart=checkpoint["path"],
                  flags=("time/nlim=12",), adiabatic=adiabatic)
    assert direct[12]["cooling_history_last_cycle"] > 1
    for cycle, state in resumed.items():
        reference = direct[cycle]
        assert state["time"] == reference["time"]
        np.testing.assert_array_equal(state["payload"], reference["payload"])
        for name in ("cooling_alpha", "cooling_history_last_time",
                     "cooling_history_last_cycle", "cooling_target_fraction",
                     "cooling_sample_time"):
            assert state[name] == reference[name], (cycle, name)


@pytest.mark.parametrize("history", ["dt = 0.1", "dcycle = 100"])
@pytest.mark.parametrize("adiabatic", [False, True])
def test_cooling_long_interval_relaxes_to_target(tmp_path, history, adiabatic):
    source = input_text(ppar=2, pperp=2, target=1, tau=0.001,
                        history=history, nlim=-1)
    if adiabatic:
        source += adiabatic_settings(pgas=2)
    final_time = 0.004
    single = run(tmp_path / "single", source, flags=(f"time/tlim={final_time}",),
                 adiabatic=adiabatic)
    subdivided = run(tmp_path / "subdivided", source,
                     flags=(f"time/tlim={final_time}", "time/cfl_number=0.005"),
                     adiabatic=adiabatic)
    assert max(single) == 1 and max(subdivided) > 4
    for states in (single, subdivided):
        final = states[max(states)]
        pp, pt = pressures(final)
        np.testing.assert_allclose(pp, 1 + np.exp(-4), rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose(pt, pp, rtol=3e-12, atol=3e-13)
        if adiabatic:
            np.testing.assert_allclose(gas_pressure(final), pp, rtol=3e-12, atol=3e-13)
        assert np.min(pp) >= 1
        assert_background_unchanged(states[0], final)


@pytest.mark.parametrize("adiabatic", [False, True])
def test_cooling_large_contrast_retains_target(tmp_path, adiabatic):
    target = 1e-16
    source = input_text(ppar=2, pperp=2, target=target, tau=1e-5,
                        floor=1e-27, nlim=1)
    source += "\n<problem>\nb0 = 1e-4\n"
    if adiabatic:
        source += adiabatic_settings(pgas=2)
    states = run(tmp_path / "contrast", source, flags=("time/tlim=0.001",),
                 adiabatic=adiabatic)
    pp, pt = pressures(states[1])
    # The magnetic energy sets an absolute storage error of about one ULP.
    np.testing.assert_allclose(pp, target, rtol=0, atol=2e-24)
    np.testing.assert_allclose(pt, target, rtol=0, atol=2e-24)
    if adiabatic:
        np.testing.assert_allclose(gas_pressure(states[1]), target, rtol=0, atol=2e-24)
    assert_background_unchanged(states[0], states[1])


def test_legacy_cooling_restart_waits_for_thermal_sample(tmp_path):
    source = input_text(history="dcycle = 3", nlim=1)
    seed = run(tmp_path / "seed", source)[1]
    header, payload = seed["path"].read_bytes().split(b"<par_end>\n", 1)
    header = re.sub(rb"^cooling_(target_fraction|sample_time)\s*=.*\n", b"",
                    header, flags=re.MULTILINE)
    checkpoint = tmp_path / "legacy.rst"
    checkpoint.write_bytes(header + b"<par_end>\n" + payload)
    resumed = run(tmp_path / "resumed", restart=checkpoint, flags=("time/nlim=4",))
    np.testing.assert_array_equal(resumed[2]["payload"], seed["payload"])
    np.testing.assert_array_equal(resumed[3]["payload"], seed["payload"])
    assert resumed[2]["cooling_alpha"] == 0
    assert resumed[3]["cooling_alpha"] > 0
    assert resumed[3]["cooling_sample_time"] == resumed[3]["time"]
    assert np.all(resumed[4]["active"][:, 4] < resumed[3]["active"][:, 4])


@pytest.mark.parametrize("adiabatic", [False, True])
def test_cooling_global_coefficient_matches_mpi(tmp_path, adiabatic):
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    config = subprocess.run([str(executable), "-c"], capture_output=True,
                            text=True, check=True).stdout
    if not re.search(r"MPI parallelism:\s+ON", config) or not shutil.which("mpirun"):
        pytest.skip("requires an MPI executable and mpirun")
    source = input_text(history="dcycle = 2", nlim=0)
    if adiabatic:
        source += adiabatic_settings()
    seed = run(tmp_path / "seed", source, adiabatic=adiabatic)[0]
    # Give the two blocks different thermal states without adding a production IC knob.
    raw = bytearray(seed["path"].read_bytes())
    offset = len(raw) - seed["payload"].nbytes
    payload = np.frombuffer(raw, dtype="=f8", offset=offset).reshape(seed["payload"].shape)
    u = payload[:, :seed["u"][0].size].reshape(seed["u"].shape)
    for block, (ppar, pperp) in enumerate(((2.0, 2.2), (3.0, 3.2))):
        if adiabatic:
            u[block, 6] = 0.5*ppar + pperp + 0.5
            u[block, 4] = np.log(ppar)
        else:
            u[block, 4] = 0.5*ppar + pperp + 0.5
        u[block, 5] = np.log(pperp/ppar)  # rho = |B| = 1.
    checkpoint = tmp_path / "nonuniform.rst"
    checkpoint.write_bytes(raw)
    one = run(tmp_path / "one", restart=checkpoint, flags=("time/nlim=4",),
              adiabatic=adiabatic)
    two = run(tmp_path / "two", restart=checkpoint, flags=("time/nlim=4",), ranks=2,
              adiabatic=adiabatic)
    for cycle, state in two.items():
        np.testing.assert_allclose(state["payload"], one[cycle]["payload"],
                                   rtol=2e-12, atol=2e-13)
        assert state["cooling_alpha"] == pytest.approx(one[cycle]["cooling_alpha"],
                                                       rel=2e-12)
    for directory in ("one", "two"):
        history = passive_acceptance.history(tmp_path / directory / "cooling.user.hst")
        thermal = "therm_mhd" if adiabatic else "therm_cgl"
        expected = np.maximum(1 - 1.5*0.5*history["volume"]/history[thermal], 0)/0.1
        np.testing.assert_allclose(history["cool_alpha"], expected, rtol=2e-12, atol=2e-13)


def test_forced_cooled_gas_is_independent_of_passive_physics(tmp_path):
    source = input_text(history="dcycle = 2", nlim=8) + adiabatic_settings()
    source += "\n<turb_driving>\ndedt = 0.2\n"
    reference = run(tmp_path / "reference", source, adiabatic=True)
    modified = source + """
<problem>
p_parallel0 = 0.0001
p_perp0 = 1.0
<mhd>
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
nu_coll = 0.7
<time>
sts_integrator = rkl2
"""
    passive = run(tmp_path / "passive", modified, adiabatic=True)
    for cycle, state in passive.items():
        expected = reference[cycle]
        assert state["time"] == expected["time"]
        assert state["dt"] == expected["dt"]
        np.testing.assert_array_equal(state["active"][:, [0, 1, 2, 3, 6]],
                                      expected["active"][:, [0, 1, 2, 3, 6]])
        for first, second in zip(state["faces"], expected["faces"]):
            np.testing.assert_array_equal(first, second)
        assert state["cooling_alpha"] == expected["cooling_alpha"]
    assert np.ptp(gas_pressure(passive[8])) > 0


@pytest.mark.parametrize("landau_fluid", [False, True])
def test_passive_ratio_cap_after_cooling(tmp_path, landau_fluid):
    source = input_text(ppar=1.0, pperp=2.0, target=0.05, tau=1e-5,
                        nlim=3) + adiabatic_settings()
    if landau_fluid:
        source += ("\n<mhd>\ncgl_heat_flux = landau_fluid\n"
                   "cgl_heat_flux_integrator = sts\n<time>\nsts_integrator = rkl2\n")
    control = run(tmp_path / "control", source, adiabatic=True)
    capped = run(tmp_path / "capped", source, adiabatic=True,
                 flags=("mhd/passive_max_pressure_ratio=3",))
    for cycle, state in capped.items():
        pp, pt = pressures(state)
        assert np.max(np.maximum(pp/pt, pt/pp)) <= 3*(1 + 2e-12)
        if cycle <= 1:
            pp0, pt0 = pressures(control[cycle])
            np.testing.assert_allclose(pt + 0.5*pp, pt0 + 0.5*pp0,
                                       rtol=3e-12, atol=3e-13)
        for index in (0, 1, 2, 3, 6):
            np.testing.assert_allclose(state["active"][:, index],
                                       control[cycle]["active"][:, index],
                                       rtol=0, atol=3e-13)
    resumed = run(tmp_path / "resumed", restart=capped[1]["path"],
                  flags=("time/nlim=3",), adiabatic=True)
    np.testing.assert_array_equal(resumed[3]["payload"], capped[3]["payload"])
