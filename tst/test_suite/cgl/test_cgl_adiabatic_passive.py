"""Adiabatic MHD flow with independently passive CGL pressures and heat fluxes."""

import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


def run(directory, *, passive=True, adiabatic=True, restart=None, flags=(), ranks=1,
        error=None, **kwargs):
    directory.mkdir()
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    command = [str(executable), "-d", str(directory)]
    if restart is None:
        source = passive_acceptance.input_text(
            passive=passive, adiabatic=passive and adiabatic, ideal_mhd=not passive,
            nx=32, nlim=6, tlim=1, **kwargs)
        source += ("\n<problem>\np_gas0 = 2\nadvect_velocity = 0.4\n"
                   "<output2>\nsingle_file_per_rank = false\n"
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
        checkpoint, passive, adiabatic=passive and adiabatic, ideal_mhd=not passive)
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


@pytest.mark.parametrize("arithmetic,flux", [
    ("safe", "weighted"), ("fast", "weighted"), ("fast", "physical"),
])
def test_nonuniform_passive_lf_diffuses_without_heating_gas(tmp_path, arithmetic, flux):
    # Case 1 with zero shear has uniform rho, gas pressure, B, and passive p_iso.
    initial = run(tmp_path / "initial", mode=1, amp=0, lf="sts",
                  flags=("time/nlim=0",))
    passive = run(tmp_path / "lf", restart=initial["path"],
                  flags=("time/nlim=6", f"mhd/cgl_lf_arithmetic={arithmetic}",
                         f"mhd/cgl_lf_sts_flux={flux}", "mhd/cgl_lf_diagnostics=none"))
    assert passive["history"]["lf_nstage"][-1] > 0
    assert passive["history"]["lf_qface"][-1] == 0
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


@pytest.mark.parametrize("max_ratio", [0, 1])
def test_adiabatic_passive_fofc_matches_native_mhd_in_3d(tmp_path, max_ratio):
    # Cold, compressive flow activates fallback on faces in all three directions.
    options = dict(forcing=True, amp=0.95, dfloor=0.2, recon="ppmx",
                   flags=("problem/p_gas0=1e-6", "time/nlim=12"))
    reference = run(tmp_path / "ideal", passive=False, **options)
    passive_options = dict(options)
    passive_options["flags"] += (f"mhd/passive_max_pressure_ratio={max_ratio}",)
    passive = run(tmp_path / "passive", **passive_options)
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


@pytest.mark.parametrize("adiabatic", [False, True])
@pytest.mark.parametrize("ratio", [1, 2, 100])
def test_passive_ratio_cap_conserves_thermal_energy(tmp_path, adiabatic, ratio):
    options = dict(adiabatic=adiabatic, mode=1, amp=0, thermal=2.5)
    initial = run(tmp_path / "initial", **options, flags=("time/nlim=0",))
    capped = run(tmp_path / "capped", **options, flags=(
        "time/nlim=0", f"mhd/passive_max_pressure_ratio={ratio}"))
    observed = initial["pt"] / initial["pp"]
    assert observed.min() < 0.5 and observed.max() > 2
    result = capped["pt"] / capped["pp"]
    np.testing.assert_allclose(result, np.clip(observed, 1 / ratio, ratio),
                               rtol=2e-13, atol=0)
    np.testing.assert_allclose(capped["U"], initial["U"], rtol=3e-14, atol=3e-14)
    np.testing.assert_array_equal(capped["active"][:, :4], initial["active"][:, :4])
    if adiabatic:
        np.testing.assert_array_equal(capped["active"][:, 6], initial["active"][:, 6])
    for first, second in zip(capped["faces"], initial["faces"]):
        np.testing.assert_array_equal(first, second)
    assert capped["fofc"] == initial["fofc"]
    if ratio == 100:
        np.testing.assert_array_equal(capped["payload"], initial["payload"])
    expected = np.count_nonzero((observed < 1 / ratio) | (observed > ratio))
    assert capped["history"]["ratio_hits"][-1] == expected
    maximum = np.maximum(observed, 1 / observed).max() if expected else 0
    assert capped["history"]["ratio_max"][-1] == pytest.approx(maximum, rel=3e-13)


@pytest.mark.parametrize("adiabatic", [False, True])
def test_passive_ratio_cap_survives_lf_and_restart(tmp_path, adiabatic):
    options = dict(adiabatic=adiabatic, mode=1, amp=0, thermal=2.5, lf="sts")
    flag = "mhd/passive_max_pressure_ratio=2"
    initial = run(tmp_path / "initial", **options, flags=(flag, "time/nlim=0"))
    reference = run(tmp_path / "reference", **options, flags=(flag,))
    partial = run(tmp_path / "partial", **options, flags=(flag, "time/nlim=2"))
    resumed = run(tmp_path / "resumed", adiabatic=adiabatic,
                  restart=partial["path"], flags=("time/nlim=6",))
    np.testing.assert_array_equal(resumed["payload"], reference["payload"])
    result = reference["pt"] / reference["pp"]
    assert result.min() >= 0.5*(1 - 2e-13)
    assert result.max() <= 2*(1 + 2e-13)
    np.testing.assert_allclose(np.mean(reference["U"]), np.mean(initial["U"]),
                               rtol=2e-13, atol=2e-13)
    assert reference["history"]["lf_nstage"][-1] > 0
    np.testing.assert_array_equal(reference["active"][:, :4], initial["active"][:, :4])
    if adiabatic:
        np.testing.assert_array_equal(reference["active"][:, 6], initial["active"][:, 6])


@pytest.mark.parametrize("adiabatic", [False, True])
def test_passive_ratio_cap_disabled_is_unchanged(tmp_path, adiabatic):
    options = dict(adiabatic=adiabatic, mode=1, amp=0, thermal=2.5, lf="sts")
    reference = run(tmp_path / "default", **options)
    disabled = run(tmp_path / "disabled", **options,
                   flags=("mhd/passive_max_pressure_ratio=0",))
    np.testing.assert_array_equal(disabled["payload"], reference["payload"])


def test_enable_ratio_cap_on_older_passive_checkpoint(tmp_path):
    original = run(tmp_path / "original", mode=1, amp=0, thermal=2.5,
                   flags=("time/nlim=0",))
    header, end, payload = original["path"].read_bytes().partition(b"<par_end>\n")
    header, removed = re.subn(rb"^passive_max_pressure_ratio\s*=.*\n", b"", header,
                              flags=re.MULTILINE)
    assert end and removed == 1
    legacy = tmp_path / "uncapped.rst"
    legacy.write_bytes(header + end + payload)
    resumed = run(tmp_path / "resumed", restart=legacy,
                  flags=("mhd/passive_max_pressure_ratio=2", "time/nlim=0"))
    ratio = original["pt"] / original["pp"]
    np.testing.assert_allclose(resumed["pt"] / resumed["pp"], np.clip(ratio, 0.5, 2),
                               rtol=2e-13)
    np.testing.assert_allclose(resumed["U"], original["U"], rtol=3e-14)


@pytest.mark.parametrize("adiabatic", [False, True])
def test_passive_ratio_cap_respects_pressure_floors(tmp_path, adiabatic):
    options = dict(adiabatic=adiabatic, mode=1, amp=0, thermal=2.5, pfloor=1.5)
    floored = run(tmp_path / "floored", **options, flags=("time/nlim=0",))
    capped = run(tmp_path / "capped", **options,
                 flags=("mhd/passive_max_pressure_ratio=1", "time/nlim=0"))
    assert floored["pp"].min() == pytest.approx(1.5)
    assert capped["history"]["ratio_hits"][-1] > 0
    assert np.minimum(capped["pp"], capped["pt"]).min() >= 1.5*(1 - 2e-13)
    np.testing.assert_allclose(capped["U"], floored["U"], rtol=3e-14)
    np.testing.assert_allclose(capped["pp"], capped["pt"], rtol=2e-13, atol=0)


@pytest.mark.parametrize("adiabatic", [False, True])
@pytest.mark.parametrize("ranks", [1, 2])
def test_passive_ratio_cap_counts_owned_cells_once(tmp_path, adiabatic, ranks):
    if ranks > 1 and not mpi_available():
        pytest.skip("requires an MPI executable and mpirun")
    initial = run(tmp_path / "initial", adiabatic=adiabatic, mode=1, amp=0,
                  thermal=2.5, flags=("mhd/passive_max_pressure_ratio=2", "time/nlim=0"))
    evolved = run(tmp_path / "evolved", adiabatic=adiabatic, mode=1, amp=0,
                  thermal=2.5, ranks=ranks, flags=("mhd/passive_max_pressure_ratio=2",))
    # There is no flow or LF here. Re-encoding must not repeatedly count boundary states.
    np.testing.assert_array_equal(evolved["history"]["ratio_hits"],
                                  initial["history"]["ratio_hits"][-1])
    np.testing.assert_array_equal(evolved["history"]["ratio_max"],
                                  initial["history"]["ratio_max"][-1])
    np.testing.assert_array_equal(evolved["payload"], initial["payload"])


@pytest.mark.parametrize("adiabatic", [False, True])
def test_passive_ratio_cap_repairs_new_lf_stage_anisotropy(tmp_path, adiabatic):
    options = dict(adiabatic=adiabatic, mode=3, amp=0, thermal=0.2)
    flags = ("mhd/passive_max_pressure_ratio=1", "problem/advect_velocity=0")
    initial = run(tmp_path / "initial", **options, flags=(*flags, "time/nlim=0"))
    control = run(tmp_path / "control", **options, flags=flags)
    conducted = run(tmp_path / "conducted", **options, lf="sts", flags=flags)
    assert control["history"]["ratio_hits"][-1] == initial["history"]["ratio_hits"][-1]
    assert conducted["history"]["lf_nstage"][-1] > 0
    # Each owned cell develops fresh anisotropy at every LF stage in this fixture.
    np.testing.assert_array_equal(
        conducted["history"]["ratio_hits"] - initial["history"]["ratio_hits"][-1],
        conducted["history"]["lf_nstage"])
    np.testing.assert_allclose(conducted["pt"], conducted["pp"], rtol=2e-13, atol=0)
    np.testing.assert_allclose(np.mean(conducted["U"]), np.mean(initial["U"]),
                               rtol=3e-13, atol=3e-13)


@pytest.mark.parametrize("value", ["-1", "0.5", "nan", "inf", "10oops", "1e-999"])
def test_passive_ratio_cap_rejects_invalid_values(tmp_path, value):
    run(tmp_path / "invalid", flags=(f"mhd/passive_max_pressure_ratio={value}",),
        error="passive_max_pressure_ratio")
