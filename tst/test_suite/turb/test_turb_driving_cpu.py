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


@pytest.mark.parametrize("projection_policy", [
    "solenoidal_compressive", "mks24_alfvenic_perpendicular",
])
def test_restart_preserves_modal_state_on_fixed_grid_and_amr(tmp_path, projection_policy):
    """Restarting retains the modal OU trajectory independently of mesh layout."""
    for name, split_cycle, final_cycle in (
        ("turb_driving_edot.athinput", 3, 6),
        ("turb_driving_amr_exclude.athinput", 2, 4),
    ):
        stem = Path(name).stem
        policy_input = f"projection_policy = {projection_policy}"
        if projection_policy == "mks24_alfvenic_perpendicular":
            policy_input += "\nsol_fraction = 1.0"
        path = staged_turb_input(tmp_path / stem, name, {
            "turb_driving": policy_input,
        })
        reference_dir = tmp_path / f"{stem}_reference"
        split_dir = tmp_path / f"{stem}_split"
        resume_dir = tmp_path / f"{stem}_resume"
        require_success(run_athena(reference_dir, path))
        if name == "turb_driving_edot.athinput":
            history = np.loadtxt(reference_dir / "turb_test_edot.hydro.hst")
            slope = np.polyfit(history[2:, 0], history[2:, 6], 1)[0]
            assert slope == pytest.approx(0.1, abs=0.015)
        require_success(
            run_athena(
                split_dir,
                path,
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
                path,
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
            other_policy = ("mks24_random_unprojected"
                            if projection_policy == "solenoidal_compressive"
                            else "solenoidal_compressive")
            incompatible = run_athena(
                tmp_path / "incompatible_projection", name,
                f"turb_driving/projection_policy={other_policy}", restart=restart,
            )
            assert incompatible.returncode != 0
            assert "configuration differs for 'projection_policy'" in (
                incompatible.stdout + incompatible.stderr
            )
            if projection_policy == "solenoidal_compressive":
                data = bytearray(restart.read_bytes())
                # This fixture has 61 independent modes in the full signed shell.
                marker = struct.pack("=2i", 5, 61)
                assert data.count(marker) == 1
                struct.pack_into("=i", data, data.index(marker), 3)
                old_restart = tmp_path / "old_format.rst"
                old_restart.write_bytes(data)
                incompatible = run_athena(
                    tmp_path / "old_format", name, restart=old_restart,
                )
                assert incompatible.returncode != 0
                assert "version" in incompatible.stdout + incompatible.stderr


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


@pytest.mark.parametrize("key,value", [
    ("driving_type", "0"), ("isotropic_power_spectrum", "true"),
    ("exp_prp", "1.666666666666667"), ("exp_prl", "0.0"),
])
def test_removed_driving_controls_are_rejected(tmp_path, key, value):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": f"{key} = {value}",
    })
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert f"removed parameter '{key}'" in result.stdout + result.stderr


@pytest.mark.parametrize("integrator", ["rk1", "rk2", "rk3"])
@pytest.mark.parametrize("projection_policy", [
    "solenoidal_compressive", "mks24_random_unprojected",
    "mks24_alfvenic_perpendicular",
])
def test_once_per_step_power(tmp_path, integrator, projection_policy):
    """Periodic total energy gains dedt*dt on every step, including the first."""
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": f"projection_policy = {projection_policy}\n"
                        "record_injected_work = true",
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


@pytest.mark.parametrize("projection_policy", [
    "solenoidal_compressive", "mks24_random_unprojected",
    "mks24_alfvenic_perpendicular",
])
def test_every_selected_mode_is_nonzero(tmp_path, projection_policy):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {})
    text = path.read_text().replace("pgen_name = turb", "pgen_name = turb_forcing")
    # Power-law modes have no legitimate zero-amplitude shell endpoints.
    text = text[:text.index("<output1>")]
    text += (f"\n<turb_driving>\nprojection_policy = {projection_policy}\n"
             "spectrum = power_law\n")
    path.write_text(text)
    result = run_athena(tmp_path / "run", path, "time/nlim=0")
    require_success(result)
    assert "nonzero modes" in result.stdout


@pytest.mark.parametrize("projection_policy", [
    "solenoidal_compressive", "mks24_alfvenic_perpendicular",
])
def test_signed_modes_projection_and_radial_spectrum(tmp_path, projection_policy):
    """Signed modes, force polarization, and radial weights are independent."""
    powers = {}
    modes = np.fft.fftfreq(16) * 16
    kz, ky, kx = np.meshgrid(modes, modes, modes, indexing="ij")
    radius = np.sqrt(kx*kx + ky*ky + kz*kz)
    shell = (radius >= 1) & (radius <= 3) & (kz != 0)
    interior = shell & (radius > 1) & (radius < 3)
    assert np.count_nonzero(shell) == 94
    for spectrum in ("power_law", "parabolic"):
        directory = tmp_path / spectrum
        path = staged_turb_input(directory, "turb_driving_edot.athinput", {
            "meshblock": "nx1 = 16\nnx2 = 16\nnx3 = 16",
            "turb_driving": f"projection_policy = {projection_policy}\n"
                            f"spectrum = {spectrum}\nexpo = 2.0\n"
                            "min_kx = -3\nmax_kx = 3\n"
                            "min_ky = -3\nmax_ky = 3\n"
                            "min_kz = 1\nmax_kz = 3",
            "output2": "dcycle = 1",
        })
        run_dir = directory / "run"
        require_success(run_athena(run_dir, path, "time/nlim=1"))
        force = read_force_blocks(latest_force(run_dir))["blocks"][0]["force"]
        if projection_policy == "mks24_alfvenic_perpendicular":
            np.testing.assert_array_equal(force[2], 0.0)
        else:
            assert np.max(np.abs(force[2])) > 0.0
        amplitude = np.fft.fftn(force, axes=(1, 2, 3))
        power = np.sum(np.abs(amplitude)**2, axis=0)
        expected = shell if spectrum == "power_law" else interior
        np.testing.assert_array_equal(power > power.max()*1e-12, expected)
        divergence = kx*amplitude[0] + ky*amplitude[1] + kz*amplitude[2]
        assert np.max(np.abs(divergence)) < 2e-6*np.max(np.abs(amplitude))
        powers[spectrum] = power
    # Fixed-seed innovations are identical; only the prescribed radial weight changes.
    # For expo=2 and npeak=2, squared weights give P_par/P_pow proportional to
    # [1-(|n|-2)^2] |n|^2, independently of the force's overall normalization.
    ratio = powers["parabolic"][interior] / powers["power_law"][interior]
    expected = (1.0 - (radius[interior] - 2.0)**2) * radius[interior]**2
    np.testing.assert_allclose(ratio / ratio[0], expected / expected[0], rtol=2e-5)


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
        "turb_driving": "rseed = 2718\ndedt = 0.1\n"
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


@pytest.mark.parametrize("direction", [(0, 0, 1), (1, 0, 0), (1, 1, 1)])
@pytest.mark.parametrize("exclude_parallel_zero", [False, True])
def test_arbitrary_field_direction_and_individual_mode_exclusions(
    tmp_path, direction, exclude_parallel_zero
):
    """The physical signed spectrum obeys both filters and both projections."""
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "meshblock": "nx1 = 16\nnx2 = 16\nnx3 = 16",
        "turb_driving": "projection_policy = mks24_alfvenic_perpendicular\n"
                        "spectrum = power_law\n"
                        f"mean_field_direction = {' '.join(map(str, direction))}\n"
                        f"exclude_k_parallel_zero = {str(exclude_parallel_zero).lower()}\n"
                        "exclude_modes = 1 0 1; -1 -1 1",
        "output2": "dcycle = 1",
    })
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/nlim=1"))
    force = read_force_blocks(latest_force(run_dir))["blocks"][0]["force"]
    amplitude = np.fft.fftn(force, axes=(1, 2, 3))
    modes = np.fft.fftfreq(16) * 16
    kz, ky, kx = np.meshgrid(modes, modes, modes, indexing="ij")
    radius2 = kx*kx + ky*ky + kz*kz
    expected = (radius2 >= 1) & (radius2 <= 9)
    if exclude_parallel_zero:
        expected &= direction[0]*kx + direction[1]*ky + direction[2]*kz != 0
    for excluded in ((1, 0, 1), (-1, -1, 1)):
        for sign in (-1, 1):
            expected &= ~((kx == sign*excluded[0]) & (ky == sign*excluded[1])
                          & (kz == sign*excluded[2]))
    power = np.sum(np.abs(amplitude)**2, axis=0)
    np.testing.assert_array_equal(power > power.max()*1e-12, expected)
    scale = np.max(np.abs(amplitude))
    assert np.max(np.abs(kx*amplitude[0] + ky*amplitude[1] + kz*amplitude[2])) < 2e-6*scale
    assert np.max(np.abs(sum(direction[i]*amplitude[i] for i in range(3)))) < 2e-6*scale


@pytest.mark.parametrize("tile_nx", [1, 2])
def test_field_direction_filter_uses_physical_wavevectors(tmp_path, tile_nx):
    """Rectangular and tiled boxes use physical k, including the full shell."""
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "mesh": "x1min = -1.0\nx1max = 1.0",
        "meshblock": "nx1 = 16\nnx2 = 16\nnx3 = 16",
        "turb_driving": "projection_policy = mks24_alfvenic_perpendicular\n"
                        "spectrum = power_law\nmean_field_direction = 1 1 1\n"
                        "exclude_k_parallel_zero = true\nexclude_modes = 1 -1 0\n"
                        "physical_k_shell = true\nk_shell_unit = 6.283185307179586\n"
                        f"tile_nx = {tile_nx}",
        "output2": "dcycle = 1",
    })
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/nlim=1"))
    force = read_force_blocks(latest_force(run_dir))["blocks"][0]["force"]
    amplitude = np.fft.fftn(force, axes=(1, 2, 3))
    modes = np.fft.fftfreq(16) * 16
    nz, ny, nx = np.meshgrid(modes, modes, modes, indexing="ij")
    kx, ky, kz = nx/2.0, ny, nz
    radius2 = kx*kx + ky*ky + kz*kz
    expected = ((radius2 >= 1) & (radius2 <= 9) & (kx+ky+kz != 0)
                & (nx % tile_nx == 0))
    for sign in (-1, 1):
        expected &= ~((nx == sign*tile_nx) & (ny == -sign) & (nz == 0))
    power = np.sum(np.abs(amplitude)**2, axis=0)
    np.testing.assert_array_equal(power > power.max()*1e-12, expected)
    scale = np.max(np.abs(amplitude))
    assert np.max(np.abs(kx*amplitude[0] + ky*amplitude[1] + kz*amplitude[2])) < 3e-6*scale
    assert np.max(np.abs(amplitude[0] + amplitude[1] + amplitude[2])) < 2e-6*scale


@pytest.mark.parametrize("direction", [(0, 0, 1), (1, 0, 0), (1, 1, 1)])
def test_wavevector_parallel_to_field_keeps_both_transverse_polarizations(tmp_path, direction):
    """At k_perp=0 both independent force polarizations are admissible."""
    bounds = "\n".join(
        f"min_k{axis} = {value}\nmax_k{axis} = {value}"
        for axis, value in zip("xyz", direction)
    )
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "meshblock": "nx1 = 16\nnx2 = 16\nnx3 = 16",
        "turb_driving": "projection_policy = mks24_alfvenic_perpendicular\n"
                        "spectrum = power_law\n"
                        f"mean_field_direction = {' '.join(map(str, direction))}\n{bounds}",
        "output2": "dcycle = 1",
    })
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/nlim=1"))
    force = read_force_blocks(latest_force(run_dir))["blocks"][0]["force"]
    amplitude = np.fft.fftn(force, axes=(1, 2, 3))
    value = amplitude[:, direction[2], direction[1], direction[0]]
    assert np.linalg.norm(value) > 0.0
    assert abs(np.dot(direction, value)) < 2e-6*np.linalg.norm(value)
    assert np.linalg.norm(np.cross(value.real, value.imag)) > (
        1e-6*np.linalg.norm(value.real)*np.linalg.norm(value.imag))


@pytest.mark.parametrize("key,value", [
    ("mean_field_direction", "0 0 0"),
    ("mean_field_direction", "nan 0 1"),
    ("mean_field_direction", "inf 0 1"),
    ("mean_field_direction", "1 2"),
    ("mean_field_direction", "1 2 3 4"),
    ("exclude_modes", "1 2"),
    ("exclude_modes", "1 2 3 4"),
    ("exclude_modes", "1 0.5 0"),
    ("exclude_modes", "1 x 0"),
])
def test_invalid_direction_and_mode_exclusions_are_rejected(tmp_path, key, value):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": f"{key} = {value}",
    })
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert key in result.stdout + result.stderr


def test_excluding_the_last_conjugate_pair_is_rejected(tmp_path):
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": "spectrum = power_law\n"
                        "min_kx = 1\nmax_kx = 1\nmin_ky = 0\nmax_ky = 0\n"
                        "min_kz = 0\nmax_kz = 0\nexclude_modes = -1 0 0",
    })
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "mode_count" in result.stdout + result.stderr


def test_filtered_oblique_forcing_restart_and_configuration_validation(tmp_path):
    """Restart validates the direction and actual mode identities, not just count."""
    path = staged_turb_input(tmp_path, "turb_driving_edot.athinput", {
        "turb_driving": "projection_policy = mks24_alfvenic_perpendicular\n"
                        "spectrum = power_law\nmean_field_direction = 1 1 1\n"
                        "exclude_k_parallel_zero = true\n"
                        "exclude_modes = 1 0 1; 1 -1 1\nrecord_injected_work = true",
    })
    reference, split, resumed = [tmp_path/name for name in ("reference", "split", "resumed")]
    require_success(run_athena(reference, path))
    require_success(run_athena(split, path, "time/nlim=3", "output2/file_type=rst",
                               "output2/dcycle=3", "output2/id=restart"))
    restart = max((split/"rst").glob("*.rst"), key=lambda file: file.stat().st_mtime)
    require_success(run_athena(
        resumed, path, "time/nlim=6", "output2/file_type=bin", "output2/variable=turb_force",
        "output2/id=resume_force", "output2/dcycle=6", restart=restart,
    ))
    assert_equal_force_outputs(latest_force(reference), latest_force(resumed))
    np.testing.assert_array_equal(np.loadtxt(next(reference.glob("*.hydro.hst")))[-1],
                                  np.loadtxt(next(resumed.glob("*.hydro.hst")))[-1])
    for index, override in enumerate((
        "turb_driving/mean_field_direction=1 0 0",
        "turb_driving/exclude_k_parallel_zero=false",
        "turb_driving/exclude_modes=1 1 0; 1 -1 1",
    )):
        result = run_athena(tmp_path/f"changed_{index}", path, override, restart=restart)
        assert result.returncode != 0
        assert "configuration differs" in result.stdout + result.stderr
    equivalent = tmp_path/"equivalent"
    require_success(run_athena(
        equivalent, path, "time/nlim=6", "output2/file_type=bin", "output2/variable=turb_force",
        "output2/id=resume_force", "output2/dcycle=6",
        "turb_driving/exclude_modes=-1 1 -1; -1 0 -1", restart=restart,
    ))
    assert_equal_force_outputs(latest_force(reference), latest_force(equivalent))
