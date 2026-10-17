"""Regression tests for the turbulence driving source term on CPU."""

import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUTS = REPO_ROOT / "tst" / "inputs"
ATHENA = Path.cwd() / "athena"
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert  # noqa: E402


def run_athena(output_dir, input_name, *overrides, restart=None):
    """Run AthenaK in an isolated output directory and return the process."""
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [str(ATHENA), "-d", str(output_dir)]
    if restart is None:
        command.extend(["-i", str(INPUTS / input_name)])
    else:
        command.extend(["-r", str(restart)])
    command.extend(overrides)
    return subprocess.run(command, capture_output=True, text=True, check=False)


def require_success(result):
    """Report executable output when a regression run fails."""
    assert result.returncode == 0, result.stdout + result.stderr


def read_force_blocks(filename):
    """Read force arrays and MeshBlock geometry from a binary output."""
    with open(filename, "rb") as file_obj:
        assert file_obj.readline().decode("ascii") == "Athena binary output version=1.1\n"
        file_obj.readline()
        time = float(file_obj.readline().decode("ascii").split("=", 1)[1])
        cycle = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        location_size = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        variable_size = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        file_obj.readline()
        variables = file_obj.readline().decode("ascii").split(":", 1)[1].split()
        header_offset = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        file_obj.read(header_offset)
        location_dtype = np.dtype("=f4" if location_size == 4 else "=f8")
        variable_dtype = np.dtype("=f4" if variable_size == 4 else "=f8")
        force_indices = [
            variables.index(name) for name in ("force1", "force2", "force3")
        ]
        blocks = []
        while True:
            raw_indices = file_obj.read(6 * 4)
            if not raw_indices:
                break
            indices = struct.unpack("=6i", raw_indices)
            logical = struct.unpack("=4i", file_obj.read(4 * 4))
            limits = np.frombuffer(file_obj.read(6 * location_size), dtype=location_dtype)
            nx = indices[1] - indices[0] + 1
            ny = indices[3] - indices[2] + 1
            nz = indices[5] - indices[4] + 1
            count = nx * ny * nz
            data = np.frombuffer(
                file_obj.read(len(variables) * count * variable_size),
                dtype=variable_dtype,
            ).reshape((len(variables), nz, ny, nx))
            blocks.append({
                "level": logical[3],
                "limits": limits,
                "force": data[force_indices, :, :, :],
            })
    return {"time": time, "cycle": cycle, "blocks": blocks}


def latest_force(run_dir):
    """Return the latest turbulence force output in a run directory."""
    outputs = []
    for path in (run_dir / "bin").glob("*.bin"):
        try:
            read_force_blocks(path)
            outputs.append(path)
        except ValueError:
            pass
    assert outputs, f"no turbulence force output found under {run_dir}"
    return max(outputs, key=lambda path: read_force_blocks(path)["cycle"])


def rms_acceleration(output):
    """Return the volume-weighted RMS acceleration over active blocks."""
    weighted_sum = 0.0
    volume_sum = 0.0
    for block in output["blocks"]:
        limits = block["limits"]
        volume = (limits[1] - limits[0]) * (limits[3] - limits[2]) * (
            limits[5] - limits[4]
        )
        magnitude_sq = np.sum(block["force"] * block["force"], axis=0)
        weighted_sum += volume * np.mean(magnitude_sq)
        volume_sum += volume
    return np.sqrt(weighted_sum / volume_sum)


def radial_magnitude_means(output):
    """Compute mean forcing inside and outside central radial regions."""
    inside = []
    outside = []
    for block in output["blocks"]:
        force = block["force"]
        limits = block["limits"]
        nx = force.shape[3]
        ny = force.shape[2]
        nz = force.shape[1]
        x = np.linspace(limits[0], limits[1], nx, endpoint=False)
        y = np.linspace(limits[2], limits[3], ny, endpoint=False)
        z = np.linspace(limits[4], limits[5], nz, endpoint=False)
        x += 0.5 * (limits[1] - limits[0]) / nx
        y += 0.5 * (limits[3] - limits[2]) / ny
        z += 0.5 * (limits[5] - limits[4]) / nz
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        radius = np.sqrt(xx * xx + yy * yy + zz * zz)
        magnitude = np.sqrt(np.sum(force * force, axis=0))
        inside.extend(magnitude[radius < 0.17])
        outside.extend(magnitude[radius > 0.38])
    return np.mean(inside), np.mean(outside)


def tile_periodicity_error(output):
    """Return the largest force mismatch between corresponding x/y tiles."""
    samples = {}
    for block in output["blocks"]:
        force = block["force"]
        limits = block["limits"]
        nx = force.shape[3]
        ny = force.shape[2]
        nz = force.shape[1]
        x = np.linspace(limits[0], limits[1], nx, endpoint=False)
        y = np.linspace(limits[2], limits[3], ny, endpoint=False)
        z = np.linspace(limits[4], limits[5], nz, endpoint=False)
        x += 0.5 * (limits[1] - limits[0]) / nx
        y += 0.5 * (limits[3] - limits[2]) / ny
        z += 0.5 * (limits[5] - limits[4]) / nz
        for k, x3 in enumerate(z):
            for j, x2 in enumerate(y):
                for i, x1 in enumerate(x):
                    key = tuple(np.round((x1, x2, x3), decimals=8))
                    samples[key] = force[:, k, j, i]
    max_error = 0.0
    for key, value in samples.items():
        if key[0] < 0.0:
            target = tuple(np.round((key[0] + 0.5, key[1], key[2]), decimals=8))
            max_error = max(max_error, np.max(np.abs(value - samples[target])))
        if key[1] < 0.0:
            target = tuple(np.round((key[0], key[1] + 0.5, key[2]), decimals=8))
            max_error = max(max_error, np.max(np.abs(value - samples[target])))
    return max_error


def assert_equal_force_outputs(left_path, right_path):
    """Check that a restart produces the identical rendered forcing field."""
    left = read_force_blocks(left_path)
    right = read_force_blocks(right_path)
    assert left["cycle"] == right["cycle"]
    assert left["time"] == right["time"]
    assert len(left["blocks"]) == len(right["blocks"])
    for left_block, right_block in zip(left["blocks"], right["blocks"]):
        assert left_block["level"] == right_block["level"]
        np.testing.assert_array_equal(left_block["limits"], right_block["limits"])
        np.testing.assert_array_equal(left_block["force"], right_block["force"])


def test_acceleration_normalization_and_localization(tmp_path):
    """Constant RMS normalization survives tiled localized fixed-grid and AMR runs."""
    tiled_dir = tmp_path / "tiled"
    require_success(run_athena(tiled_dir, "turb_driving_tiled_include.athinput"))
    tiled = read_force_blocks(latest_force(tiled_dir))
    assert rms_acceleration(tiled) == pytest.approx(0.25, rel=2.0e-6)
    inner, outer = radial_magnitude_means(tiled)
    assert inner > outer
    periodic_dir = tmp_path / "tiled_periodic"
    require_success(
        run_athena(
            periodic_dir,
            "turb_driving_tiled_include.athinput",
            "turb_driving/localization=none",
            "turb_driving/sigma_x1=-1.0",
            "turb_driving/sigma_x2=-1.0",
            "turb_driving/sigma_x3=-1.0",
        )
    )
    assert tile_periodicity_error(read_force_blocks(latest_force(periodic_dir))) < 1.0e-12

    amr_dir = tmp_path / "amr"
    require_success(run_athena(amr_dir, "turb_driving_amr_exclude.athinput"))
    amr = read_force_blocks(latest_force(amr_dir))
    assert max(block["level"] for block in amr["blocks"]) > 0
    assert rms_acceleration(amr) == pytest.approx(0.20, rel=2.0e-6)
    inner, outer = radial_magnitude_means(amr)
    assert inner < outer


def test_restart_preserves_modal_state_on_fixed_grid_and_amr(tmp_path):
    """Restarting retains the modal OU trajectory independently of mesh layout."""
    for name, split_cycle, final_cycle in (
        ("turb_driving_edot.athinput", 3, 6),
        ("turb_driving_amr_exclude.athinput", 2, 4),
    ):
        stem = Path(name).stem
        reference_dir = tmp_path / f"{stem}_reference"
        split_dir = tmp_path / f"{stem}_split"
        resume_dir = tmp_path / f"{stem}_resume"
        require_success(run_athena(reference_dir, name))
        if name == "turb_driving_edot.athinput":
            history = np.loadtxt(reference_dir / "turb_test_edot.hydro.hst")
            slope = np.polyfit(history[2:, 0], history[2:, 6], 1)[0]
            assert slope == pytest.approx(0.1, abs=0.015)
        require_success(
            run_athena(
                split_dir,
                name,
                f"time/nlim={split_cycle}",
                "output2/file_type=rst",
                f"output2/dcycle={split_cycle}",
                "output2/id=restart",
            )
        )
        restart_files = list((split_dir / "rst").glob("*.rst"))
        assert restart_files
        restart = max(restart_files, key=lambda path: path.stat().st_mtime)
        require_success(
            run_athena(
                resume_dir,
                name,
                f"time/nlim={final_cycle}",
                "output2/file_type=bin",
                "output2/variable=turb_force",
                "output2/id=resume_force",
                f"output2/dcycle={final_cycle}",
                restart=restart,
            )
        )
        assert_equal_force_outputs(latest_force(reference_dir), latest_force(resume_dir))
        if name == "turb_driving_edot.athinput":
            reference = np.loadtxt(next(reference_dir.glob("*.hydro.hst")))[-1]
            resumed = np.loadtxt(next(resume_dir.glob("*.hydro.hst")))[-1]
            assert np.all(np.isfinite(reference)) and np.all(np.isfinite(resumed))
            np.testing.assert_array_equal(reference, resumed)
            incompatible = run_athena(
                tmp_path / "incompatible_restart",
                name,
                "turb_driving/npeak=1.5",
                restart=restart,
            )
            assert incompatible.returncode != 0
            assert "configuration differs for 'npeak'" in (
                incompatible.stdout + incompatible.stderr
            )


def test_removed_parameter_and_degenerate_spectrum_are_rejected(tmp_path):
    """Removed aliases and singular spectra fail during input validation."""
    removed = run_athena(tmp_path / "removed", "turb_driving_invalid_removed.athinput")
    assert removed.returncode != 0
    assert "removed parameter 'constant_edot'" in removed.stdout + removed.stderr

    spectrum = run_athena(tmp_path / "spectrum", "turb_driving_invalid_spectrum.athinput")
    assert spectrum.returncode != 0
    assert "parabolic spectrum requires nhigh greater than nlow" in (
        spectrum.stdout + spectrum.stderr
    )


def test_zero_force_normalization_reports_replay_context(tmp_path):
    """A finite zero force fails with enough state to diagnose an exact replay."""
    result = run_athena(
        tmp_path / "zero_force",
        "turb_driving_tiled_include.athinput",
        "turb_driving/sigma_x1=1.0e-100",
        "turb_driving/sigma_x2=1.0e-100",
        "turb_driving/sigma_x3=1.0e-100",
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "cannot impose non-zero accel_rms with a zero forcing field" in output
    for field in (
        "time=",
        "cycle=",
        "update=",
        "mode_count=",
        "dt=",
        "t0=",
        "t1=",
        "totvol=",
        "m0=",
        "m1=",
    ):
        assert field in output


def test_edot_zero_force_is_not_reported_as_zero_timestep(tmp_path):
    """A zero edot field is diagnosed independently of its valid timestep."""
    result = run_athena(
        tmp_path / "zero_edot_force",
        "turb_driving_edot.athinput",
        "turb_driving/localization=include",
        "turb_driving/sigma_x1=1.0e-100",
        "turb_driving/sigma_x2=1.0e-100",
        "turb_driving/sigma_x3=1.0e-100",
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "cannot inject non-zero dedt with a zero forcing field" in output
    assert "edot forcing requires a finite positive timestep" not in output
    assert "dt=" in output


def test_nonpositive_final_cycle_timestep_is_rejected(tmp_path):
    """The mesh rejects a final zero timestep with cycle-budget diagnostics."""
    result = run_athena(
        tmp_path / "zero_dt",
        "turb_driving_edot.athinput",
        "time/cfl_number=0.0",
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "Final cycle timestep must be finite, positive, and advance time" in output
    for field in (
        "rank=",
        "cycle=",
        "time=",
        "tlim=",
        "dt=",
        "time_plus_dt=",
        "dt_before_tlim=",
        "dtold=",
        "dt_legacy=",
        "dt_cycle_candidate=",
        "dt_parabolic_sts=",
        "cfl_number=",
        "sts_integrator=",
        "sts_max_dt_ratio=",
    ):
        assert field in output


def staged_turb_input(tmp_path, name, additions):
    """Add canonical keys to a local input; Athena CLI cannot add missing keys."""
    text = (INPUTS / name).read_text()
    for block, keys in additions.items():
        text += f"\n<{block}>\n{keys}\n"
    path = tmp_path / "forcing.athinput"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def test_type_two_has_nonzero_force(tmp_path):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": "driving_type = 2",
    })
    result = run_athena(tmp_path / "run", path)
    require_success(result)
    assert rms_acceleration(read_force_blocks(latest_force(tmp_path / "run"))) > 0.0


def test_type_two_rejects_projected_policy(tmp_path):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": "driving_type = 2\nprojection_policy = solenoidal_compressive",
    })
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "driving_type = 2 requires mks24_random_unprojected" in (
        result.stdout + result.stderr
    )


@pytest.mark.parametrize("integrator", ["rk1", "rk2", "rk3"])
@pytest.mark.parametrize("driving_type", [0, 2])
def test_once_per_step_power(tmp_path, integrator, driving_type):
    """Periodic total energy gains dedt*dt on every step, including the first."""
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": f"driving_type = {driving_type}\nrecord_injected_work = true",
        "output1": "data_format = %24.16e",
    })
    result = run_athena(tmp_path / "run", path, f"time/integrator={integrator}")
    require_success(result)
    history = next((tmp_path / "run").glob("*.hydro.hst"))
    header = history.read_text().splitlines()[1].split()
    names = [field.split("=", 1)[1] for field in header if "=" in field]
    data = np.loadtxt(history)
    dt = np.diff(data[:, names.index("time")])
    work = np.diff(data[:, names.index("tot-E")])
    assert np.all(dt >= 0.0)
    assert np.count_nonzero(dt > 0.0) == 6
    np.testing.assert_array_equal(work[dt == 0.0], 0.0)
    np.testing.assert_allclose(work[dt > 0.0] / dt[dt > 0.0], 0.1,
                               rtol=2e-7, atol=2e-8)


@pytest.mark.parametrize("hydro_eos,mhd_eos", [
    (None, "ideal"), (None, "cgl"), (None, "isothermal"),
    ("ideal", None), ("isothermal", None),
    ("ideal", "ideal"), ("isothermal", "isothermal"),
    ("ideal", "isothermal"), ("isothermal", "ideal"),
])
def test_conservative_forcing_kick(tmp_path, hydro_eos, mhd_eos):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {})
    text = path.read_text()
    start = text.index("<hydro>")
    end = text.index("<problem>", start)
    blocks = []
    for fluid, eos in (("hydro", hydro_eos), ("mhd", mhd_eos)):
        if eos is not None:
            blocks.append(f"<{fluid}>\neos = {eos}\ngamma = 1.666666666666667\n"
                          "reconstruct = plm\nrsolver = hlle\nnscalars = 1\n"
                          "iso_sound_speed = 1.0\n")
    if hydro_eos and mhd_eos:
        blocks.append("<ion-neutral>\ndrag_coeff = 0.0\n")
    text = text[:start] + "\n".join(blocks) + "\n" + text[end:]
    if hydro_eos and mhd_eos:
        text = text.replace("integrator = rk2", "integrator = imex2")
    text = text.replace("pgen_name = turb", "pgen_name = turb_forcing")
    text = text[:text.index("<output1>")]
    text += "\n<turb_driving>\nspectrum = power_law\n"
    path.write_text(text)
    result = run_athena(tmp_path / "run", path, "time/nlim=0")
    require_success(result)
    assert "max kick error=" in result.stdout


@pytest.mark.parametrize("model", ["hydro", "mhd", "isothermal", "active", "passive"])
@pytest.mark.parametrize("normalization", ["edot", "accel_rms"])
def test_final_stage_kick_from_uniform_state(tmp_path, model, normalization):
    """The first cycle is an exact impulse after the uniform RK3 update."""
    fluid = "hydro" if model == "hydro" else "mhd"
    eos = {"active": "cgl", "passive": "cgl", "isothermal": "isothermal"}.get(model, "ideal")
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        fluid: f"eos = {eos}\n"
               f"passive = {str(model == 'passive').lower()}\n"
               "gamma = 1.666666666666667\nrsolver = hlle\niso_sound_speed = 1.0",
        "turb_driving": f"normalization = {normalization}",
    })
    text = path.read_text().replace("<hydro>", f"<{fluid}>")
    if normalization == "accel_rms":
        text = text.replace("dedt           = 0.1", "accel_rms = 0.25")
    if fluid == "mhd":
        text = text.replace("variable  = hydro_w", "variable  = mhd_w_bcc")
    path.write_text(text)
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/nlim=1", "time/integrator=rk3"))
    initial, final = [bin_convert.read_binary(str(path)) for path in sorted(
        (run_dir / "bin").glob("*.prim.*.bin"))]
    force = bin_convert.read_binary(str(latest_force(run_dir)))
    assert initial["cycle"] == 0 and final["cycle"] == force["cycle"] == 1
    dt = final["time"] - initial["time"]
    assert dt > 0.0
    if model == "isothermal":
        assert dt == pytest.approx(0.3 / (16 * np.sqrt(3.0)), rel=1e-12)
    for name, value in final["mb_data"].items():
        if name in ("velx", "vely", "velz"):
            axis = ("velx", "vely", "velz").index(name) + 1
            expected = dt * np.asarray(force["mb_data"][f"force{axis}"])
        else:
            expected = initial["mb_data"][name]
        # Binary fields are stored as float32; the direct kick test uses Real.
        np.testing.assert_allclose(value, expected, rtol=3e-7, atol=1e-9,
                                   err_msg=name)
    if normalization == "accel_rms":
        assert rms_acceleration(read_force_blocks(latest_force(run_dir))) == (
            pytest.approx(0.25, rel=2e-6))


@pytest.mark.parametrize("driving_type", [1, 2])
def test_every_selected_mode_is_nonzero(tmp_path, driving_type):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {})
    text = path.read_text().replace("pgen_name = turb", "pgen_name = turb_forcing")
    # Power-law modes have no legitimate zero-amplitude shell endpoints.
    text = text[:text.index("<output1>")]
    text += f"\n<turb_driving>\ndriving_type = {driving_type}\nspectrum = power_law\n"
    path.write_text(text)
    result = run_athena(tmp_path / "run", path, "time/nlim=0")
    require_success(result)
    assert "nonzero modes" in result.stdout


def test_type_two_matches_unprojected_type_zero(tmp_path):
    """Equal fixed-seed isotropic policies render the same OU force trajectory."""
    outputs = []
    for driving_type in (0, 2):
        directory = tmp_path / f"type_{driving_type}"
        path = staged_turb_input(directory, "turb_driving_edot.athinput", {
            "turb_driving": f"driving_type = {driving_type}\n"
                            "projection_policy = mks24_random_unprojected\n"
                            "spectrum = power_law\nrseed = 2718",
        })
        output_dir = directory / "run"
        require_success(run_athena(output_dir, path))
        output = latest_force(output_dir)
        assert rms_acceleration(read_force_blocks(output)) > 0.0
        outputs.append(output)
    assert_equal_force_outputs(*outputs)


@pytest.mark.parametrize("integrator", ["imex2", "imex3"])
@pytest.mark.parametrize("ifield", [1, 2])
def test_two_fluid_once_per_step_power(tmp_path, integrator, ifield):
    """Both periodic fluid energies gain dedt*dt*volume once per full step."""
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "hydro": "gamma = 1.666666666666667",
        "mhd": "eos = ideal\ngamma = 1.666666666666667\n"
               "reconstruct = plm\nrsolver = hlle",
        "ion-neutral": "drag_coeff = 0.0",
        "problem": f"d_n = 1.3\nd_i = 0.7\nifield = {ifield}",
        "turb_driving": "driving_type = 0\nrseed = 2718\ndedt = 0.1\n"
                        "record_injected_work = false",
        "output1": "data_format = %24.16e",
    })
    output_dir = tmp_path / "run"
    require_success(run_athena(output_dir, path, f"time/integrator={integrator}"))
    times, energies = [], []
    for fluid, density in (("hydro", 1.3), ("mhd", 0.7)):
        history = next(output_dir.glob(f"*.{fluid}.hst"))
        header = history.read_text().splitlines()[1].split()
        names = [field.split("=", 1)[1] for field in header if "=" in field]
        data = np.loadtxt(history)
        # The inherited periodic cube has unit volume.
        assert data[0, names.index("mass")] == pytest.approx(density, rel=1e-12)
        magnetic_energy = (2.0 if ifield == 2 else 1.0) if fluid == "mhd" else 0.0
        assert data[0, names.index("tot-E")] == pytest.approx(
            0.9*density + magnetic_energy, rel=1e-12)
        for axis in (1, 2, 3):
            assert data[0, names.index(f"{axis}-KE")] == 0.0
        times.append(data[:, names.index("time")])
        energies.append(data[:, names.index("tot-E")])
    np.testing.assert_array_equal(*times)
    dt = np.diff(times[0])
    work = np.diff(energies[0] + energies[1])
    assert np.all(dt >= 0.0)
    assert np.count_nonzero(dt > 0.0) == 6
    np.testing.assert_array_equal(work[dt == 0.0], 0.0)
    np.testing.assert_allclose(work[dt > 0.0] / dt[dt > 0.0], 0.1,
                               rtol=2e-7, atol=2e-8)


def test_zero_power_with_opposing_flow(tmp_path):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {})
    text = path.read_text().replace("pgen_name = turb", "pgen_name = turb_forcing")
    path.write_text(text[:text.index("<output1>")] +
                    "\n<turb_driving>\nspectrum = power_law\n")
    result = run_athena(tmp_path / "run", path, "time/nlim=0", "turb_driving/dedt=0")
    require_success(result)
    assert "zero power with opposing flow passed" in result.stdout
