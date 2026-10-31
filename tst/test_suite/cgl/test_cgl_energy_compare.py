"""Thermal-energy comparisons use the same magnetic representation as the EOS."""

import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/plot_cgl_energy_compare.py"


@pytest.fixture(scope="module")
def plotter():
    spec = importlib.util.spec_from_file_location("cgl_energy_compare", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_history(path, fields):
    header = "# Athena++ history data\n#  " + " ".join(
        f"[{index}]={name}" for index, name in enumerate(fields, 1)) + "\n"
    with path.open("w") as stream:
        stream.write(header)
        np.savetxt(stream, np.column_stack(list(fields.values())), fmt="%.17e")


def energy_history():
    # Face checkerboards Bx and By have Bcc=0. A unit guide field contributes 0.5.
    # The face-square magnetic history includes another unit of magnetic energy.
    return {"time": np.array([0., 1.]), "mass": np.ones(2),
            "tot-E": np.array([2., 2.5]), "1-KE": np.zeros(2),
            "2-KE": np.zeros(2), "3-KE": np.zeros(2),
            "1-ME": np.full(2, .5), "2-ME": np.full(2, .5),
            "3-ME": np.full(2, .5)}


@pytest.mark.parametrize("physical", ["B^2", "therm_cgl", "magnetic"])
def test_checkerboard_thermal_energy(plotter, physical):
    value = {"B^2": [1., 1.], "therm_cgl": [1.5, 2.], "magnetic": [.5, .5]}
    history = energy_history()
    user = {"time": history["time"], physical: np.array(value[physical])}
    result = plotter.diagnostics(history, user)
    np.testing.assert_array_equal(result["thermal_energy"], [1.5, 2.])
    np.testing.assert_array_equal(result["thermal_energy_change"], [0., .5])
    np.testing.assert_array_equal(result["magnetic_energy"], [.5, .5])
    np.testing.assert_array_equal(result["alfven_rms"], [1., 1.])


def test_passive_thermal_energy_uses_independent_gas(plotter):
    history = energy_history()
    history["cgl-J"] = np.zeros(2)
    history["thermal-U"] = np.array([4., 5.])
    user = {"time": history["time"], "therm_cgl": np.array([4., 5.]),
            "therm_mhd": np.array([1.5, 2.])}
    result = plotter.diagnostics(history, user)
    np.testing.assert_array_equal(result["thermal_energy"], [1.5, 2.])


@pytest.mark.parametrize("magnetic_key,magnetic_value", [("B^2", 1.), ("magnetic", .5)])
@pytest.mark.parametrize("passive", [False, True])
def test_direct_magnetic_energy_survives_large_thermal_energy(
    plotter, magnetic_key, magnetic_value, passive,
):
    history = energy_history()
    history["tot-E"] = np.full(2, 1e30)
    if passive:
        history["cgl-J"] = np.zeros(2)
    thermal_key = "therm_mhd" if passive else "therm_cgl"
    user = {"time": history["time"], thermal_key: np.full(2, 1e30),
            magnetic_key: np.full(2, magnetic_value)}
    result = plotter.diagnostics(history, user)
    np.testing.assert_array_equal(result["thermal_energy"], [1e30, 1e30])
    np.testing.assert_array_equal(result["magnetic_energy"], [.5, .5])
    np.testing.assert_array_equal(result["alfven_rms"], [1., 1.])


@pytest.mark.parametrize("times", [[0., 1.001], [0.], [.1, 1.],
                                   [0., np.nextafter(1., 0.)],
                                   [0., np.nextafter(1., 2.)]])
def test_mismatched_sample_times_fail(plotter, times):
    user = {"time": np.array(times), "B^2": np.ones(len(times))}
    with pytest.raises(ValueError, match="matching sample times"):
        plotter.diagnostics(energy_history(), user)


def test_missing_physical_energy_reports_snapshot_requirement(plotter):
    with pytest.raises(ValueError, match="from snapshots"):
        plotter.diagnostics(energy_history())


def test_resumed_history_replaces_abandoned_branch(plotter, tmp_path):
    path = tmp_path / "run.mhd.hst"
    write_history(path, {"time": [0., .5, 1., .5, .75, 1., 1.],
                         "tot-E": [1., 2., 3., 4., 5., 6., 7.]})
    data = plotter.read_history(path)
    np.testing.assert_array_equal(data["time"], [0., .5, .75, 1.])
    np.testing.assert_array_equal(data["tot-E"], [1., 4., 5., 7.])


def test_restart_history_schema_change_is_reported(plotter, tmp_path):
    path = tmp_path / "run.mhd.hst"
    write_history(path, {"time": [0.], "tot-E": [1.]})
    with path.open("a") as stream:
        stream.write("# [1]=time [2]=different\n1 2\n")
    with pytest.raises(ValueError, match="columns change"):
        plotter.read_history(path)


def test_cli_uses_companion_user_histories(tmp_path):
    for name in ("mhd", "cgl"):
        write_history(tmp_path / f"{name}.mhd.hst", energy_history())
        # An extra user sample is allowed. Only matching MHD times are used.
        write_history(tmp_path / f"{name}.user.hst",
                      {"time": [0., .5, 1.], "B^2": [1., 100., 1.]})
    output = tmp_path / "plots"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--mhd-history", str(tmp_path / "mhd.mhd.hst"),
         "--cgl-history", str(tmp_path / "cgl.mhd.hst"), "--output-dir", str(output)],
        capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    values = np.genfromtxt(output / "energy_velocity_comparison.csv", delimiter=",",
                           names=True, dtype=None, encoding="utf-8")
    np.testing.assert_array_equal(values["thermal_energy"], [1.5, 2., 1.5, 2.])
    assert len(list(output.glob("*.png"))) == 5
