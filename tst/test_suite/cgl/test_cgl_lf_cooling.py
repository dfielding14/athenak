"""History-held CGL cooling, floors, and native checkpoint continuation.

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


def run(directory, source=None, *, restart=None, flags=(), ranks=1):
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
        state = passive_acceptance.decode_restart(path, False, active_cgl=True)
        header = path.read_bytes().split(b"<par_end>\n", 1)[0].decode()
        for name in ("cooling_alpha", "cooling_history_last_time",
                     "cooling_history_last_cycle"):
            match = re.search(rf"^{name}\s*=\s*(\S+)", header, re.MULTILINE)
            assert match, name
            state[name] = float(match[1])
        state["path"] = path
        states[state["cycle"]] = state
    assert states
    return states


def pressures(state):
    u, b = state["active"], state["bcc"]
    rho = u[:, 0]
    bsqr = np.sum(b*b, axis=1)
    thermal = u[:, 4] - 0.5*np.sum(u[:, 1:4]**2, axis=1)/rho - 0.5*bsqr
    ratio = np.exp(u[:, 5]/rho - 2*np.log(rho) + 1.5*np.log(bsqr))
    ppar = thermal/(0.5 + ratio)
    return ppar, ratio*ppar


def assert_background_unchanged(initial, final):
    np.testing.assert_allclose(final["active"][:, :4], initial["active"][:, :4],
                               rtol=0, atol=3e-13)
    for before, after in zip(initial["faces"], final["faces"]):
        np.testing.assert_allclose(after, before, rtol=0, atol=3e-13)


@pytest.mark.parametrize("integrator", ["rk2", "rk4"])
def test_cooling_exact_held_interval(tmp_path, integrator):
    states = run(tmp_path / integrator, input_text(integrator=integrator))
    initial = states[0]
    piso0, delta = (2.0 + 2*2.2)/3, 0.2
    alpha = (1 - 0.5/piso0)/0.1
    for state in states.values():
        ppar, pperp = pressures(state)
        piso = piso0*np.exp(-alpha*state["time"])
        np.testing.assert_allclose(ppar, piso - 2*delta/3, rtol=3e-12, atol=3e-13)
        np.testing.assert_allclose(pperp, piso + delta/3, rtol=3e-12, atol=3e-13)
        assert state["cooling_alpha"] == pytest.approx(alpha, rel=2e-13)
        assert state["cooling_history_last_cycle"] == 0
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
    alpha = (1 - 0.5/piso0)/0.1
    for state in states.values():
        parallel, perpendicular = pressures(state)
        piso = piso0*np.exp(-alpha*state["time"])
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
            expected = (old_par + 2*old_perp)/3 * np.exp(
                -previous["cooling_alpha"]*(state["time"] - previous["time"]))
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


@pytest.mark.parametrize("checkpoint_kind", ["interior", "off_cadence_final"])
def test_cooling_restart_preserves_history_schedule(tmp_path, checkpoint_kind):
    source = input_text(history="dt = 0.04", nlim=12)
    direct = run(tmp_path / "direct", source)
    # Cycle 1 is before the first scheduled update, with a different current pressure.
    assert direct[1]["time"] < 0.04 and direct[1]["cooling_history_last_cycle"] == 0
    if checkpoint_kind == "interior":
        checkpoint = direct[1]
    else:
        checkpoint = run(tmp_path / "partial", source, flags=("time/nlim=1",))[1]
    resumed = run(tmp_path / "resumed", restart=checkpoint["path"],
                  flags=("time/nlim=12",))
    assert direct[12]["cooling_history_last_cycle"] > 1
    for cycle, state in resumed.items():
        reference = direct[cycle]
        assert state["time"] == reference["time"]
        np.testing.assert_array_equal(state["payload"], reference["payload"])
        for name in ("cooling_alpha", "cooling_history_last_time",
                     "cooling_history_last_cycle"):
            assert state[name] == reference[name], (cycle, name)


def test_cooling_global_coefficient_matches_mpi(tmp_path):
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    config = subprocess.run([str(executable), "-c"], capture_output=True,
                            text=True, check=True).stdout
    if not re.search(r"MPI parallelism:\s+ON", config) or not shutil.which("mpirun"):
        pytest.skip("requires an MPI executable and mpirun")
    source = input_text(history="dcycle = 2", nlim=0)
    seed = run(tmp_path / "seed", source)[0]
    # Give the two blocks different thermal states without adding a production IC knob.
    raw = bytearray(seed["path"].read_bytes())
    offset = len(raw) - seed["payload"].nbytes
    payload = np.frombuffer(raw, dtype="=f8", offset=offset).reshape(seed["payload"].shape)
    u = payload[:, :seed["u"][0].size].reshape(seed["u"].shape)
    for block, (ppar, pperp) in enumerate(((2.0, 2.2), (3.0, 3.2))):
        u[block, 4] = 0.5*ppar + pperp + 0.5
        u[block, 5] = np.log(pperp/ppar)  # rho = |B| = 1.
    checkpoint = tmp_path / "nonuniform.rst"
    checkpoint.write_bytes(raw)
    one = run(tmp_path / "one", restart=checkpoint, flags=("time/nlim=4",))
    two = run(tmp_path / "two", restart=checkpoint, flags=("time/nlim=4",), ranks=2)
    for cycle, state in two.items():
        np.testing.assert_allclose(state["payload"], one[cycle]["payload"],
                                   rtol=2e-12, atol=2e-13)
        assert state["cooling_alpha"] == pytest.approx(one[cycle]["cooling_alpha"],
                                                       rel=2e-12)
    for directory in ("one", "two"):
        history = passive_acceptance.history(tmp_path / directory / "cooling.user.hst")
        expected = np.maximum(1 - 1.5*0.5*history["volume"]/history["therm_cgl"], 0)/0.1
        np.testing.assert_allclose(history["cool_alpha"], expected, rtol=2e-12, atol=2e-13)
