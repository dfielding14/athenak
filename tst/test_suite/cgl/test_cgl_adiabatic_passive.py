"""Adiabatic MHD flow with independently passive CGL pressures and heat fluxes."""

import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


def run(directory, *, passive=True, restart=None, flags=(), ranks=1,
        error=None, **kwargs):
    directory.mkdir()
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    command = [str(executable), "-d", str(directory)]
    if restart is None:
        source = passive_acceptance.input_text(
            passive=passive, adiabatic=passive, ideal_mhd=not passive,
            nx=32, nlim=6, tlim=1, **kwargs)
        source += ("\n<problem>\np_gas0 = 2\n<output2>\nsingle_file_per_rank = false\n"
                   "<output3>\nfile_type = log\ndcycle = 1\n")
        deck = directory / "input.athinput"
        deck.write_text(source)
        command += ["-i", str(deck)]
    else:
        command += ["-r", str(restart)]
    if ranks > 1:
        command = ["mpirun", "-np", str(ranks), *command]
    result = subprocess.run(command + list(flags), cwd=directory, capture_output=True,
                            text=True, check=False, timeout=120)
    log = result.stdout + result.stderr
    (directory / "run.log").write_text(log)
    if error is not None:
        assert result.returncode != 0 and error in log, log
        return None
    assert result.returncode == 0, log
    checkpoint = sorted((directory / "rst").glob("*.rst"))[-1]
    state = passive_acceptance.decode_restart(
        checkpoint, passive, adiabatic=passive, ideal_mhd=not passive)
    state["path"] = checkpoint
    state["history"] = passive_acceptance.history(directory / "passive_validation.mhd.hst")
    event_rows = (directory / "passive_validation.log").read_text().splitlines()
    state["fofc"] = sum(int(row.split()[-1]) for row in event_rows
                        if row.strip() and not row.startswith("#"))
    return state


def assert_flow_equal(passive, ideal):
    assert passive["time"] == ideal["time"]
    assert passive["dt"] == ideal["dt"]
    np.testing.assert_array_equal(passive["active"][:, :4], ideal["active"][:, :4])
    np.testing.assert_array_equal(passive["active"][:, 6], ideal["active"][:, 4])
    for first, second in zip(passive["faces"], ideal["faces"]):
        np.testing.assert_array_equal(first, second)


@pytest.mark.parametrize("recon", ["dc", "plm", "ppm4", "ppmx", "wenoz"])
def test_adiabatic_passive_matches_native_mhd(tmp_path, recon):
    reference = run(tmp_path / "ideal", passive=False, recon=recon)
    passive = run(tmp_path / "passive", recon=recon)
    assert_flow_equal(passive, reference)
    assert np.ptp(passive["pp"]) > 0.01
    assert np.ptp(passive["pt"]) > 0.01


@pytest.mark.parametrize("source", ["forcing", "acceleration"])
def test_adiabatic_passive_work_updates_native_energy(tmp_path, source):
    reference = run(tmp_path / "ideal", passive=False, **{source: True})
    passive = run(tmp_path / "passive", **{source: True})
    assert_flow_equal(passive, reference)
    assert reference["history"]["tot-E"][-1] > reference["history"]["tot-E"][0]


def test_passive_collisions_and_lf_leave_adiabatic_flow_unchanged(tmp_path):
    reference = run(tmp_path / "ideal", passive=False)
    passive = run(tmp_path / "passive", lf="sts", nu=0.7,
                  flags=("problem/p0=0.2",))
    assert_flow_equal(passive, reference)
    history = passive["history"]
    assert history["lf_nstage"][-1] > 0
    assert history["lf_qface"][-1] > 0


def test_nonuniform_passive_lf_diffuses_without_heating_gas(tmp_path):
    # Case 1 with zero shear has uniform rho, gas pressure, B, and passive p_iso.
    initial = run(tmp_path / "initial", mode=1, amp=0,
                  flags=("time/nlim=0",))
    passive = run(tmp_path / "lf", mode=1, amp=0, lf="sts")
    np.testing.assert_allclose(passive["active"][:, :4], initial["active"][:, :4],
                               rtol=0, atol=3e-15)
    np.testing.assert_array_equal(passive["active"][:, 6], initial["active"][:, 6])
    for first, second in zip(passive["faces"], initial["faces"]):
        np.testing.assert_array_equal(first, second)
    assert np.std(passive["pp"]) < np.std(initial["pp"])
    assert np.std(passive["pt"]) < np.std(initial["pt"])
    np.testing.assert_allclose(np.mean(passive["U"]), np.mean(initial["U"]),
                               rtol=2e-13, atol=2e-13)


def test_adiabatic_passive_density_floor_matches_native_flow(tmp_path):
    reference = run(tmp_path / "ideal", passive=False, amp=0.95, dfloor=0.2)
    passive = run(tmp_path / "passive", amp=0.95, dfloor=0.2,
                  flags=("problem/p0=0.02", "mhd/nu_coll=3"))
    assert_flow_equal(passive, reference)


def test_passive_pressure_repairs_leave_native_gas_unchanged(tmp_path):
    reference = run(tmp_path / "ideal", passive=False)
    passive = run(tmp_path / "passive", flags=("problem/p0=1e-14",))
    assert_flow_equal(passive, reference)
    # Independent log decoding can differ from the device's FMA ordering by a few ulps.
    assert np.min(passive["pp"]) >= 1e-12*(1 - 1e-14)
    assert np.min(passive["pt"]) >= 1e-12*(1 - 1e-14)


def test_adiabatic_passive_fofc_matches_native_mhd_in_3d(tmp_path):
    # Cold, compressive flow activates fallback on faces in all three directions.
    options = dict(forcing=True, amp=0.95, dfloor=0.2, recon="ppmx",
                   flags=("problem/p_gas0=1e-6", "time/nlim=12"))
    reference = run(tmp_path / "ideal", passive=False, **options)
    passive = run(tmp_path / "passive", **options)
    assert reference["fofc"] > 0
    assert passive["fofc"] == reference["fofc"]
    assert_flow_equal(passive, reference)


@pytest.mark.parametrize("flag,error", [
    ("mhd/passive_eos=unknown", "passive_eos must be"),
    ("mhd/passive=false", "passive_eos must be"),
    ("mhd/gamma=1", "gamma"),
])
def test_adiabatic_passive_rejects_invalid_eos(tmp_path, flag, error):
    run(tmp_path / "invalid", flags=(flag,), error=error)


def test_adiabatic_passive_restart_rejects_eos_change(tmp_path):
    partial = run(tmp_path / "partial", flags=("time/nlim=0",))
    run(tmp_path / "wrong_eos", restart=partial["path"],
        flags=("mhd/passive_eos=isothermal", "mhd/passive_restart_encoding=1"),
        error="passive restarts require J/A encoding")


def mpi_available():
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    config = subprocess.run([str(executable), "-c"], capture_output=True,
                            text=True, check=True).stdout
    return re.search(r"MPI parallelism:\s+ON", config) and shutil.which("mpirun")


@pytest.mark.parametrize("ranks", [1, 2])
def test_adiabatic_passive_restart_preserves_all_seven_slots(tmp_path, ranks):
    if ranks > 1 and not mpi_available():
        pytest.skip("requires an MPI executable and mpirun")
    reference = run(tmp_path / "reference", lf="sts", nu=0.7)
    partial = run(tmp_path / "partial", lf="sts", nu=0.7,
                  flags=("time/nlim=2",))
    resumed = run(tmp_path / "resumed", restart=partial["path"], ranks=ranks,
                  flags=("time/nlim=6",))
    assert resumed["time"] == reference["time"]
    assert resumed["dt"] == reference["dt"]
    np.testing.assert_array_equal(resumed["payload"], reference["payload"])
