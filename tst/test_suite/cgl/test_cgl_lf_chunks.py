"""LF interval refinement preserves outer cadence, storage and restart semantics.

Use the suite's prebuilt athena executable (or its MPI launch adapter).
"""
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

from test_suite import testutils
from test_suite.cgl import passive_acceptance as passive
from test_suite.cgl import cgl_lf_smr_conservation as conservation


def run(directory, extra="", *, passive_model=True, restart=None, nlim=3,
        tlim=1, thermal=.2, expect_failure=False, source=None):
    directory.mkdir()
    if source is None:
        source = passive.input_text(nx=64, mode=3, lf="sts", nlim=nlim,
                                    tlim=tlim, thermal=thermal)
        source = source.replace("passive = true", f"passive = {str(passive_model).lower()}")
        source += "\n<problem>\nadvect_velocity = 0\n<time>\nsts_max_dt_ratio = 16\n"
    deck = directory / "input.athinput"
    deck.write_text(source + extra)
    command = [str(Path("athena").resolve()), "-i", str(deck), "-d", str(directory)]
    if restart:
        command += ["-r", str(restart)]
    result = subprocess.run(command, text=True, capture_output=True)
    log = result.stdout + result.stderr
    (directory / "stdout.log").write_text(log)
    if expect_failure:
        assert result.returncode != 0 and "cgl_lf_max_chunk_ratio" in log, log[-3000:]
        assert not list(directory.glob("*.hst"))
        return
    assert result.returncode == 0, log[-4000:]
    path = sorted((directory / "rst").glob("*.rst"))[-1]
    state = passive.decode_restart(path, passive_model, active_cgl=not passive_model)
    state.update(restart=path, hist=passive.history(next(directory.glob("*.mhd.hst"))),
                 log=log)
    for key in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos"):
        assert np.all(state["hist"][key] == 0), key
    return state


def counts(state):
    match = re.search(r"CGL LF chunks: max_ratio=(\S+) logical_sweeps=(\d+) "
                      r"chunks=(\d+) rhs_evaluations=(\d+)", state["log"])
    assert match, state["log"][-3000:]
    _, sweeps, chunks, rhs = match.groups()
    return tuple(map(int, (sweeps, chunks, rhs)))


@pytest.mark.parametrize("passive_model", [False, True])
def test_default_absent_and_zero_are_identical(tmp_path, passive_model):
    a = run(tmp_path / "absent", passive_model=passive_model)
    b = run(tmp_path / "zero", "\n<time>\ncgl_lf_max_chunk_ratio = 0\n",
            passive_model=passive_model)
    assert np.array_equal(a["payload"].view("=u8"), b["payload"].view("=u8"))
    assert a["dt"] == b["dt"]
    for name in a["hist"]:
        np.testing.assert_array_equal(a["hist"][name], b["hist"][name])
    assert "CGL LF chunks:" not in a["log"] + b["log"]


@pytest.mark.parametrize("zero_operator", [False, True])
def test_uniform_collision_and_zero_operator_cadence(tmp_path, zero_operator):
    # Use the existing uniform limiter fixture; a transverse field gives an
    # exactly zero 1D LF operator while retaining a nonzero mirror threshold.
    repo = Path(__file__).resolve().parents[3]
    source = (repo / "inputs/unit_tests/cgl_lf_limiter_mirror.athinput").read_text()
    source += """
<mesh>
nghost = 4
<time>
evolution = dynamic
nlim = 3
tlim = 1
sts_max_dt_ratio = 16
<mhd>
rsolver = hlle
nu_coll = 2
limiter_nu_coll = 10
backup_limiters = false
<problem>
amp = 0
ppar0 = 1
pperp0 = 1.8
<output1>
file_type = hst
dcycle = 1
data_format = %25.17e
<output2>
file_type = rst
dcycle = 100000
"""
    if zero_operator:
        source += "\n<problem>\nb0 = 0\nby0 = 1\n"
    a = run(tmp_path / "off", source=source, passive_model=False)
    b = run(tmp_path / "on", "\n<time>\ncgl_lf_max_chunk_ratio = .5\n",
            source=source, passive_model=False)
    np.testing.assert_array_equal(a["hist"]["time"], b["hist"]["time"])
    np.testing.assert_allclose(a["active"], b["active"], rtol=0, atol=3e-12)
    h = b["hist"]
    time, first = np.unique(h["time"], return_index=True)
    ratio = np.exp(h["aam-D"][first] / h["mass"][first])  # rho=B=1
    measured = 4.6 * (ratio - 1) / (1 + 2 * ratio)
    expected = [.8]
    for dt in np.diff(time):
        delta = expected[-1] * np.exp(-2 * dt)
        expected.append(.5 + (delta - .5)/(1 + 10*dt) if delta > .5 else delta)
    np.testing.assert_allclose(measured, expected, rtol=2e-11, atol=2e-13)
    sweeps, chunks, rhs = counts(b)
    if zero_operator:
        assert rhs == chunks == sweeps == 0 and h["lf_nstage"][-1] == 0
    else:
        assert sweeps == 6 and chunks > sweeps and rhs >= 3*chunks


def test_smooth_chunk_refinement_and_conservation(tmp_path):
    states, errors = {}, []
    for ratio in (2, 1, .5, .125, .0625):
        # The native isothermal limit fixes the outer step independently of LF pressures.
        state = run(tmp_path / f"ratio-{ratio}",
                    f"\n<time>\ncgl_lf_max_chunk_ratio = {ratio}\n"
                    "sts_max_dt_ratio = -1\n"
                    "<mhd>\nlf_coefficient_mode = local\n", nlim=-1, tlim=.002)
        states[ratio] = state
        assert np.all(state["active"][:, 0] == 1) and np.all(state["active"][:, 1:4] == 0)
        assert abs(float(state["U"].mean() - 3)) < 5e-12
        sweeps, chunks, rhs = counts(state)
        assert sweeps == 2*state["cycle"] and chunks > sweeps and rhs == 3*chunks
        np.testing.assert_array_equal(state["hist"]["time"], states[2]["hist"]["time"])
    pressure = lambda s: np.stack((s["pp"], s["pt"]))
    reference = pressure(states[.0625])
    uncertainty = np.mean(abs(reference - pressure(states[.125])), axis=(1, 2, 3, 4))
    for ratio in (2, 1, .5):
        errors.append(np.mean(abs(pressure(states[ratio]) - reference), axis=(1, 2, 3, 4)))
    errors = np.array(errors)
    order = np.log2(errors[:-1]/errors[1:])
    (tmp_path / "refinement.json").write_text(json.dumps(dict(
        pressure_L1_errors=errors.tolist(), orders=order.tolist(),
        reference_difference=uncertainty.tolist()), indent=2) + "\n")
    assert np.all((order > 1.7) & (order < 2.4)), order
    assert np.all(uncertainty < .1*errors[-1]), (uncertainty, errors)


@pytest.mark.parametrize("passive_model", [False, True])
def test_chunked_restart_matches_uninterrupted(tmp_path, passive_model):
    extra = "\n<time>\ncgl_lf_max_chunk_ratio = 1\n<mhd>\nnu_coll = .3\n"
    first = run(tmp_path / "first", extra, passive_model=passive_model, nlim=2)
    resumed = run(tmp_path / "resumed", "", passive_model=passive_model,
                  source="<time>\nnlim = 4\n", restart=first["restart"])
    direct = run(tmp_path / "direct", extra, passive_model=passive_model, nlim=4)
    assert np.array_equal(resumed["payload"].view("=u8"), direct["payload"].view("=u8"))
    assert resumed["time"] == direct["time"] and resumed["dt"] == direct["dt"]
    for name in ("time", "dt", "lf_nstage"):
        assert resumed["hist"][name][-1] == direct["hist"][name][-1]
    assert counts(resumed)[0] == 4


def test_chunked_mpi_state_and_schedule(tmp_path):
    common = ["mesh_refinement/refinement=none", "time/evolution=dynamic",
              "mhd/rsolver=hlle", "time/nlim=2", "time/tlim=1",
              "mhd/cgl_lf_strict_admissibility=true"]
    common += [f"mesh/nx{n}=16" for n in (1, 2, 3)]
    common += [f"meshblock/nx{n}=8" for n in (1, 2, 3)]
    # An overlay supplies the new key; all ordinary controls exist in this deck.
    overlay = tmp_path / "chunk.athinput"
    overlay.write_text(Path(conservation.INPUT).read_text() +
                       "\n<time>\ncgl_lf_max_chunk_ratio = .5\n")
    reference = None
    for ranks in (1, 4):
        directory = tmp_path / f"ranks-{ranks}"
        directory.mkdir()
        testutils.mpi_run(str(overlay), ["-d", str(directory), *common], threads=ranks)
        path = sorted((directory / "rst").glob("*.rst"))[-1]
        state = conservation.read_state(path)
        assert state["nmb"] == 8 and state["cycle"] == 2
        h = passive.history(next(directory.glob("*.mhd.hst")))
        assert h["lf_nstage"][-1]/16**3 > 6*state["cycle"]
        payload_bytes = 8*(7*12**3 + 3*13*12**2)*np.dtype("=f8").itemsize
        payload = path.read_bytes()[-payload_bytes:]
        if reference is None:
            reference = (state, payload, h)
        conservation.assert_same(reference[0], state)
        assert payload == reference[1]
        for name in ("time", "dt", "lf_nstage"):
            np.testing.assert_array_equal(h[name], reference[2][name])
        for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos"):
            assert np.all(h[name] == 0)


@pytest.mark.parametrize("extra", [
    "<time>\ncgl_lf_max_chunk_ratio = -1\n",
    "<time>\ncgl_lf_max_chunk_ratio = nan\n",
    "<time>\ncgl_lf_max_chunk_ratio = inf\n",
    "<time>\ncgl_lf_max_chunk_ratio = 1\nsts_integrator = none\n"
    "<mhd>\ncgl_heat_flux_integrator = explicit\n",
    "<time>\ncgl_lf_max_chunk_ratio = 1\n<mesh>\nix1_bc = outflow\nox1_bc = outflow\n",
    "<time>\ncgl_lf_max_chunk_ratio = 1\n<mhd>\nviscosity = .001\n",
])
def test_invalid_chunk_configurations(tmp_path, extra):
    run(tmp_path / "invalid", "\n" + extra, passive_model=False, expect_failure=True)
