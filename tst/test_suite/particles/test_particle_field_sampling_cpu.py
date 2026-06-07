"""CPU regression tests for particle thermodynamic-history field sampling."""

from pathlib import Path
import shutil
import sys

import numpy as np

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


HYDRO_INPUT = str(ROOT / "tst/inputs/particle_field_sampling_hydro.athinput")
MHD_INPUT = str(ROOT / "tst/inputs/particle_field_sampling_mhd.athinput")
RUN_UNIFORM = Path("run_particle_field_sampling_uniform")
RUN_HYDRO = Path("run_particle_field_sampling_hydro")
RUN_MHD = Path("run_particle_field_sampling_mhd")
RUN_INVALID = Path("run_particle_field_sampling_invalid")
GAMMA = 5.0 / 3.0


def _history_path(run_dir, basename, sampling):
    return (
        run_dir
        / "prtcl_thermo_history"
        / f"{basename}.{sampling}.thp"
    )


def _final_particle_data(run_dir, basename):
    cell = read_history(_history_path(run_dir, basename, "cell"))
    cic = read_history(_history_path(run_dir, basename, "cic"))
    final_cycle = int(np.max(cell["cycle"]))
    cell_mask = cell["cycle"] == final_cycle
    cic_mask = cic["cycle"] == final_cycle
    cell_order = np.argsort(cell["tag"][cell_mask])
    cic_order = np.argsort(cic["tag"][cic_mask])
    cell_final = {name: values[cell_mask][cell_order] for name, values in cell.items()}
    cic_final = {name: values[cic_mask][cic_order] for name, values in cic.items()}
    np.testing.assert_array_equal(cell_final["tag"], cic_final["tag"])
    for coordinate in ("x1", "x2", "x3"):
        np.testing.assert_array_equal(cell_final[coordinate], cic_final[coordinate])
    return cell, cic, cell_final, cic_final


def _mesh_table(run_dir, basename):
    mesh = testutils.athena_read.tab(
        run_dir / "tab" / f"{basename}.mesh.00001.tab"
    )
    order = np.argsort(mesh["x1v"])
    return {name: values[order] if isinstance(values, np.ndarray) else values
            for name, values in mesh.items()}


def _periodic_indices(x, centers):
    dx = centers[1] - centers[0]
    xmin = centers[0] - 0.5 * dx
    grid_coordinate = (x - xmin) / dx - 0.5
    lower_unwrapped = np.floor(grid_coordinate).astype(int)
    upper_weight = grid_coordinate - lower_unwrapped
    lower = lower_unwrapped % centers.size
    upper = (lower + 1) % centers.size
    containing = np.floor((x - xmin) / dx).astype(int) % centers.size
    return containing, lower, upper, upper_weight


def _assert_reconstructed_sampling(cell, cic, centers, cell_fields):
    containing, lower, upper, upper_weight = _periodic_indices(cell["x1"], centers)
    for name, grid_values in cell_fields.items():
        expected_cell = grid_values[containing]
        expected_cic = (
            (1.0 - upper_weight) * grid_values[lower]
            + upper_weight * grid_values[upper]
        )
        np.testing.assert_allclose(cell[name], expected_cell, rtol=0.0, atol=2.0e-13)
        np.testing.assert_allclose(cic[name], expected_cic, rtol=0.0, atol=2.0e-13)

    periodic_stencil = (cell["x1"] < centers[0]) | (cell["x1"] >= centers[-1])
    assert np.count_nonzero(periodic_stencil) > 0
    assert np.max(np.abs(cic["density"] - cell["density"])) > 1.0e-8


def _hydro_cell_fields(mesh):
    pressure = (GAMMA - 1.0) * mesh["eint"]
    temperature = pressure / mesh["dens"]
    speed = np.sqrt(mesh["velx"] ** 2 + mesh["vely"] ** 2 + mesh["velz"] ** 2)
    sound_speed = np.sqrt(GAMMA * pressure / mesh["dens"])
    return {
        "density": mesh["dens"],
        "pressure": pressure,
        "temperature": temperature,
        "v1": mesh["velx"],
        "mach": speed / sound_speed,
    }


def _mhd_cell_fields(mesh):
    pressure = (GAMMA - 1.0) * mesh["eint"]
    bmag = np.sqrt(mesh["bcc1"] ** 2 + mesh["bcc2"] ** 2 + mesh["bcc3"] ** 2)
    magnetic_pressure = 0.5 * bmag ** 2
    return {
        "density": mesh["dens"],
        "v1": mesh["velx"],
        "b1": mesh["bcc1"],
        "b2": mesh["bcc2"],
        "b3": mesh["bcc3"],
        "bmag": bmag,
        "beta": pressure / magnetic_pressure,
        "alfven_speed": bmag / np.sqrt(mesh["dens"]),
    }


def test_particle_field_sampling_uniform_and_hydro_cic():
    """Default cell sampling is compatible; CIC matches periodic linear interpolation."""
    shutil.rmtree(RUN_UNIFORM, ignore_errors=True)
    shutil.rmtree(RUN_HYDRO, ignore_errors=True)
    try:
        assert testutils.run(
            HYDRO_INPUT,
            ["-d", str(RUN_UNIFORM), "problem/amp=0.0"],
        )
        cell_all, cic_all, _, _ = _final_particle_data(
            RUN_UNIFORM, "particle_field_sampling_hydro"
        )
        for field in ("density", "pressure", "temperature", "v1", "mach"):
            np.testing.assert_allclose(
                cell_all[field], cic_all[field], rtol=0.0, atol=2.0e-14
            )

        assert testutils.run(HYDRO_INPUT, ["-d", str(RUN_HYDRO)])
        _, _, cell, cic = _final_particle_data(
            RUN_HYDRO, "particle_field_sampling_hydro"
        )
        mesh = _mesh_table(RUN_HYDRO, "particle_field_sampling_hydro")
        _assert_reconstructed_sampling(
            cell, cic, mesh["x1v"], _hydro_cell_fields(mesh)
        )
    finally:
        shutil.rmtree(RUN_UNIFORM, ignore_errors=True)
        shutil.rmtree(RUN_HYDRO, ignore_errors=True)


def test_particle_field_sampling_mhd_and_validation():
    """MHD diagnostics use the same CIC stencil, and invalid modes are rejected."""
    shutil.rmtree(RUN_MHD, ignore_errors=True)
    shutil.rmtree(RUN_INVALID, ignore_errors=True)
    try:
        assert testutils.run(MHD_INPUT, ["-d", str(RUN_MHD)])
        _, _, cell, cic = _final_particle_data(
            RUN_MHD, "particle_field_sampling_mhd"
        )
        mesh = _mesh_table(RUN_MHD, "particle_field_sampling_mhd")
        _assert_reconstructed_sampling(
            cell, cic, mesh["x1v"], _mhd_cell_fields(mesh)
        )

        assert not testutils.run_command(
            [
                "./athena",
                "-i",
                HYDRO_INPUT,
                "-d",
                str(RUN_INVALID),
                "output1/particle_field_sampling=nearest",
            ]
        )
    finally:
        shutil.rmtree(RUN_MHD, ignore_errors=True)
        shutil.rmtree(RUN_INVALID, ignore_errors=True)
