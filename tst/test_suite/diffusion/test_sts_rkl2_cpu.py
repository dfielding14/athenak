"""Exact diffusion profiles, explicit defaults, mixed scheduling, and validation."""

import math
import os
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis/python"))
from bin_convert import read_binary  # noqa: E402

AMP = 0.01
# Diffusion acts before sound waves cross one cell, making the error-function
# solution accurate without the isothermal-only kinematic advection solver.
DIFFUSIVITY = 100.0
GM1 = 2.0 / 3.0
COEFFICIENTS = {"conduction": "conductivity", "viscosity": "viscosity",
                "ohmic": "ohmic_resistivity"}
SELECTORS = {"conduction": "conductivity", "viscosity": "viscosity",
             "ohmic": "resistivity"}


def flags_for(process, integrator="sts", fluid=None):
    fluid = fluid or ("mhd" if process == "ohmic" else "hydro")
    coefficient = DIFFUSIVITY / GM1 if process == "conduction" else DIFFUSIVITY
    args = [f"{fluid}/{COEFFICIENTS[process]}={coefficient}"]
    if integrator != "default":
        args += [f"{fluid}/{SELECTORS[process]}_integrator={integrator}"]
        args += ["time/sts_integrator=" + ("rkl2" if integrator == "sts" else "none")]
    if integrator == "sts":
        args += ["time/sts_max_dt_ratio=8"]
    if process == "conduction":
        args += [f"problem/pl={1 + AMP}", f"problem/pr={1 - AMP}"]
    elif process == "viscosity":
        args += [f"problem/vl={AMP}", f"problem/vr={-AMP}"]
    else:
        args += [f"problem/bzl={AMP}", f"problem/bzr={-AMP}"]
    return fluid, args


def run_case(directory, fluid, flags=(), ranks=1, extra="", restart=None, check=True):
    directory.mkdir(parents=True, exist_ok=True)
    binary = Path(os.environ.get("ATHENA", "./athena")).resolve()
    text = (ROOT / "tst/inputs/sts_diffusion.athinput").read_text()
    text = text.replace("<hydro>", f"<{fluid}>").replace(
        "variable = hydro_u", "variable = mhd_u_bcc" if fluid == "mhd"
        else "variable = hydro_u")
    if fluid == "mhd":
        text = text.replace("rsolver = hllc", "rsolver = hlld")
    source = directory / "input.athinput"
    text += extra
    if not restart:
        # This fork only permits command-line overrides for existing input keys.
        for setting in flags:
            block, assignment = setting.split("/", 1)
            text += f"\n<{block}>\n{assignment}\n"
    source.write_text(text)
    command = shlex.split(os.environ.get("MPIEXEC", "mpirun")) + ["-np", str(ranks)]
    command += ([str(binary), "-r", str(restart)] if restart
                else [str(binary), "-i", str(source)])
    result = subprocess.run(command + (list(flags) if restart else []), cwd=directory,
                            capture_output=True, text=True, timeout=90)
    if check:
        assert result.returncode == 0, result.stdout + result.stderr
    return result


def read_state(directory, first=False):
    files = sorted((directory / "bin").glob("*.state.*.bin"))
    assert files
    return read_binary(files[0] if first else files[-1])


def cells(data):
    """Return coordinates, cell volumes and all fields in stable physical order."""
    coords, volumes, fields = [], [], {name: [] for name in data["var_names"]}
    for m, extent in enumerate(data["mb_geometry"]):
        nz, ny, nx = np.asarray(data["mb_data"][data["var_names"][0]][m]).shape
        axes = [np.linspace(extent[2*d], extent[2*d+1], n, endpoint=False)
                + (extent[2*d+1] - extent[2*d]) / (2 * n)
                for d, n in enumerate((nx, ny, nz))]
        z, y, x = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
        coords.extend(zip(x.ravel(), y.ravel(), z.ravel()))
        volume = np.prod(extent[1::2] - extent[::2]) / (nx * ny * nz)
        volumes.extend([volume] * (nx * ny * nz))
        for name in fields:
            fields[name].extend(np.ravel(data["mb_data"][name][m]))
    coords = np.array(coords)
    order = np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))
    return coords[order], np.array(volumes)[order], {
        key: np.asarray(value)[order] for key, value in fields.items()}


@pytest.mark.parametrize("fluid, process", [
    ("hydro", "conduction"), ("hydro", "viscosity"), ("mhd", "conduction"),
    ("mhd", "viscosity"), ("mhd", "ohmic"),
])
@pytest.mark.parametrize("integrator", ["explicit", "sts"])
def test_exact_diffusion_profile_converges(tmp_path, fluid, process, integrator):
    errors = []
    for nx in (64, 128):
        fluid, flags = flags_for(process, integrator, fluid)
        directory = tmp_path / str(nx)
        run_case(directory, fluid, flags + [f"mesh/nx1={nx}"])
        data = read_state(directory)
        coords, volumes, state = cells(data)
        width = math.sqrt(4 * DIFFUSIVITY * data["time"])
        reference = -AMP * np.array([math.erf(x / width) for x in coords[:, 0]])
        if process == "conduction":
            measured = GM1 * state["ener"] - 1
        elif process == "viscosity":
            measured = state["mom2"]
        else:
            measured = state["bcc3"]
        errors.append(np.sum(np.abs(measured - reference) * volumes) / volumes.sum())
        assert np.max(np.abs(measured)) <= AMP * 1.00001
    assert errors[1] < 2e-5, errors
    assert errors[1] < 0.4 * errors[0], errors


@pytest.mark.parametrize("process", ["conduction", "viscosity", "ohmic"])
def test_explicit_defaults_are_unchanged(tmp_path, process):
    results = []
    for mode in ("default", "explicit"):
        fluid, flags = flags_for(process, mode)
        directory = tmp_path / mode
        run_case(directory, fluid, flags + ["time/nlim=3", "time/tlim=10"])
        results.append(read_state(directory))
    assert results[0]["cycle"] == results[1]["cycle"]
    assert results[0]["time"] == results[1]["time"]
    for field in results[0]["var_names"]:
        np.testing.assert_array_equal(results[0]["mb_data"][field],
                                      results[1]["mb_data"][field])


@pytest.mark.parametrize("second_mode", ["explicit", "sts"])
def test_mixed_diffusion_processes(tmp_path, second_mode):
    fluid, flags = flags_for("conduction")
    flags += [f"hydro/viscosity={DIFFUSIVITY}",
              f"hydro/viscosity_integrator={second_mode}",
              f"problem/vl={AMP}", f"problem/vr={-AMP}"]
    run_case(tmp_path, fluid, flags + ["mesh/nx1=128"])
    data = read_state(tmp_path)
    coords, _, state = cells(data)
    reference = -AMP * np.array([math.erf(x / math.sqrt(4 * DIFFUSIVITY * data["time"]))
                                for x in coords[:, 0]])
    np.testing.assert_allclose(state["mom2"], reference, atol=3e-5, rtol=0)
    # Viscous heat is O(AMP^2), so allow its physical contribution to temperature.
    np.testing.assert_allclose(GM1 * state["ener"] - 1, reference, atol=1e-4, rtol=0)


@pytest.mark.parametrize("fluid, process", [
    ("hydro", "conduction"), ("mhd", "ohmic"),
])
def test_dynamic_sts_matches_explicit_reference(tmp_path, fluid, process):
    """Resolve sound propagation as well as diffusion to exercise Strang coupling."""
    results = []
    for mode in ("explicit", "sts"):
        _, flags = flags_for(process, mode, fluid)
        coefficient = 0.25 / GM1 if process == "conduction" else 0.25
        flags += [f"{fluid}/{COEFFICIENTS[process]}={coefficient}", "time/tlim=0.2",
                  "time/cfl_number=0.1"]
        directory = tmp_path / mode
        run_case(directory, fluid, flags)
        results.append(read_state(directory))
    assert results[0]["time"] == results[1]["time"]
    assert results[1]["cycle"] < results[0]["cycle"]
    _, _, reference = cells(results[0])
    _, _, actual = cells(results[1])
    for field in reference:
        np.testing.assert_allclose(actual[field], reference[field], atol=1e-5, rtol=0)


@pytest.mark.parametrize("flags, message", [
    (["time/sts_integrator=rkl2"], "at least one active"),
    (["hydro/viscosity=1", "hydro/viscosity_integrator=sts"], "sts_integrator = none"),
    (["hydro/viscosity_integrator=bogus"], "must be 'explicit' or 'sts'"),
    (["time/sts_integrator=rkl2", "hydro/conductivity=1",
      "hydro/conductivity_integrator=sts", "hydro/sat_hflux=true"], "constant isotropic"),
    (["time/sts_integrator=rkl2", "hydro/viscosity_integrator=sts"], "positive"),
])
def test_invalid_sts_configuration_fails(tmp_path, flags, message):
    result = run_case(tmp_path, "hydro", flags, check=False)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
