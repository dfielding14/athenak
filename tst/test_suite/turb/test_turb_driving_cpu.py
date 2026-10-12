"""Forcing regressions; build with -DPROBLEM=tests/turb_forcing."""

import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUTS = REPO_ROOT / "tst" / "inputs"
ATHENA = Path.cwd() / "athena"


def run_athena(output_dir, input_name, *overrides, restart=None):
    """Run AthenaK in an isolated output directory and return the process."""
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [str(ATHENA), "-d", str(output_dir)]
    if restart is None:
        command.extend(["-i", str(INPUTS / input_name)])
    else:
        command.extend(["-r", str(restart)])
    command.extend(overrides)
    return subprocess.run(command, capture_output=True, text=True, check=False, timeout=60)


def require_success(result):
    """Report executable output when a regression run fails."""
    assert result.returncode == 0, result.stdout + result.stderr


def forcing_input(tmp_path, hydro_eos="ideal", mhd_eos=None,
                  normalization="edot", target=0.1, check_kick=False):
    """Prepare a uniform periodic fixture using the existing forcing parameters."""
    text = (INPUTS / "turb_driving_edot.athinput").read_text()
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
    text = text.replace("<problem>", f"<problem>\ncheck_kick = {str(check_kick).lower()}")
    text = text.replace("<output1>", "<output1>\ndata_format = %24.16e")
    text = text.replace("normalization  = edot", f"normalization  = {normalization}")
    target_key = "dedt" if normalization == "edot" else "accel_rms"
    text = text.replace("dedt           = 0.1", f"{target_key} = {target}")
    if hydro_eos is None:
        text = text.replace("variable  = hydro_w", "variable  = mhd_w")
    path = tmp_path / "forcing.athinput"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def assert_step_power(run_dir, fluid="hydro", power=0.1, cycles=6):
    """Check independent periodic-box energy increments, including cycle zero."""
    data = np.loadtxt(next(run_dir.glob(f"*.{fluid}.hst")))
    np.testing.assert_allclose(data[:, 3:6], 0.0, atol=3.0e-12)
    dt = np.diff(data[:, 0])
    work = np.diff(data[:, 6])
    assert np.all(dt >= 0.0)
    assert np.count_nonzero(dt > 0.0) == cycles
    np.testing.assert_array_equal(work[dt == 0.0], 0.0)
    np.testing.assert_allclose(work[dt > 0.0], power * dt[dt > 0.0],
                               rtol=2.0e-10, atol=3.0e-12)


def read_binary_blocks(filename, names):
    """Read selected arrays and MeshBlock geometry from a binary output."""
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
        variable_indices = [variables.index(name) for name in names]
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
                "data": data[variable_indices, :, :, :],
            })
    return {"time": time, "cycle": cycle, "blocks": blocks}


def read_force_blocks(filename):
    """Read force arrays and MeshBlock geometry from a binary output."""
    output = read_binary_blocks(filename, ("force1", "force2", "force3"))
    for block in output["blocks"]:
        block["force"] = block.pop("data")
    return output


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
        path = tmp_path / name
        text = (INPUTS / name).read_text()
        if name == "turb_driving_amr_exclude.athinput":
            text += "\n<output3>\nfile_type = hst\ndcycle = 1\ndata_format = %24.16e\n"
        else:
            text = text.replace("<output1>", "<output1>\ndata_format = %24.16e")
        path.write_text(text)
        reference_dir = tmp_path / f"{stem}_reference"
        split_dir = tmp_path / f"{stem}_split"
        resume_dir = tmp_path / f"{stem}_resume"
        require_success(run_athena(reference_dir, path))
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
        reference = np.loadtxt(next(reference_dir.glob("*.hydro.hst")))[-1]
        resumed = np.loadtxt(next(resume_dir.glob("*.hydro.hst")))[-1]
        assert np.all(np.isfinite(reference)) and np.all(np.isfinite(resumed))
        np.testing.assert_array_equal(reference, resumed)
        if name == "turb_driving_edot.athinput":
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


@pytest.mark.parametrize("integrator", ["rk1", "rk2", "rk3", "rk4"])
@pytest.mark.parametrize("fluid", ["hydro", "mhd"])
def test_once_per_step_power(tmp_path, integrator, fluid):
    """The two half kicks supply dedt times box volume on every timestep."""
    path = forcing_input(tmp_path, "ideal" if fluid == "hydro" else None,
                         "ideal" if fluid == "mhd" else None)
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, f"time/integrator={integrator}",
                               "mesh/x1max=1.5"))
    assert_step_power(run_dir, fluid, power=0.2)


@pytest.mark.parametrize("hydro_eos,mhd_eos", [
    ("ideal", None), ("isothermal", None),
    (None, "ideal"), (None, "isothermal"),
    ("ideal", "ideal"), ("isothermal", "isothermal"),
    ("ideal", "isothermal"), ("isothermal", "ideal"),
])
def test_conservative_forcing_kick(tmp_path, hydro_eos, mhd_eos):
    """Stale primitives and unequal fluids cannot alter kick work or EOS guards."""
    path = forcing_input(tmp_path, hydro_eos, mhd_eos, check_kick=True)
    result = run_athena(tmp_path / "run", path, "time/nlim=0")
    require_success(result)
    assert "max kick error=" in result.stdout


def test_acceleration_default_matches_explicit(tmp_path):
    """Omitting normalization selects the same fixed RMS acceleration."""
    explicit = forcing_input(tmp_path, normalization="accel_rms", target=0.25)
    implicit = tmp_path / "default.athinput"
    implicit.write_text(explicit.read_text().replace("normalization  = accel_rms\n", ""))
    for label, path in (("explicit", explicit), ("default", implicit)):
        require_success(run_athena(tmp_path / label, path))
        assert rms_acceleration(read_force_blocks(latest_force(tmp_path / label))) == (
            pytest.approx(0.25, rel=2.0e-6))
    assert_equal_force_outputs(latest_force(tmp_path / "explicit"),
                               latest_force(tmp_path / "default"))


@pytest.mark.parametrize("normalization", ["edot", "accel_rms"])
def test_zero_target_applies_no_force(tmp_path, normalization):
    """A zero target preserves the uniform state and produces zero acceleration."""
    path = forcing_input(tmp_path, normalization=normalization, target=0.0)
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path))
    assert_step_power(run_dir, power=0.0)
    assert rms_acceleration(read_force_blocks(latest_force(run_dir))) == 0.0


@pytest.mark.parametrize("velocity,target", [(-1.0, 0.0), (1.0e12, 0.1)])
def test_edot_normalization_root(tmp_path, velocity, target):
    """Opposing flow at zero power stays unforced; positive work avoids cancellation."""
    path = forcing_input(tmp_path, target=target)
    path.write_text(path.read_text().replace(
        "<problem>", f"<problem>\nnormalization_velocity = {velocity}"))
    result = run_athena(tmp_path / "run", path, "time/nlim=0")
    require_success(result)
    assert "normalization check passed" in result.stdout


@pytest.mark.parametrize("integrator", ["imex2", "imex2+"])
def test_two_fluid_step_power(tmp_path, integrator):
    """Both fluids share the prescribed power across the two half kicks."""
    path = forcing_input(tmp_path, "ideal", "ideal")
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/nlim=3",
                               f"time/integrator={integrator}"))
    hydro = np.loadtxt(next(run_dir.glob("*.hydro.hst")))
    mhd = np.loadtxt(next(run_dir.glob("*.mhd.hst")))
    np.testing.assert_array_equal(hydro[:, 0], mhd[:, 0])
    np.testing.assert_allclose(hydro[:, 3:6] + mhd[:, 3:6], 0.0, atol=3.0e-12)
    dt = np.diff(hydro[:, 0])
    assert np.count_nonzero(dt > 0.0) == 3
    np.testing.assert_allclose(np.diff(hydro[:, 6] + mhd[:, 6]), 0.1*dt,
                               rtol=2.0e-10, atol=3.0e-12)
    assert np.sum(hydro[-1, 7:10]) > 0.0
    assert np.sum(mhd[-1, 7:10]) > 0.0


@pytest.mark.parametrize("integrator", ["rk2", "rk3"])
@pytest.mark.parametrize("fluid", ["hydro", "mhd"])
def test_frozen_forcing_second_order(tmp_path, integrator, fluid):
    """At fixed grid, smooth frozen forcing converges quadratically in time."""
    path = forcing_input(tmp_path, "ideal" if fluid == "hydro" else None,
                         "ideal" if fluid == "mhd" else None,
                         normalization="accel_rms", target=0.01)
    path.write_text(path.read_text().replace(
        "<turb_driving>", "<turb_driving>\nmin_kx = 1\nmax_kx = 1\nmin_ky = 0\n"
        "max_ky = 0\nmin_kz = 0\nmax_kz = 0").replace(
        "<output3>", "<output3>\ndata_precision = real"))
    states = []
    for cfl in (0.4, 0.2, 0.1, 0.05):
        run_dir = tmp_path / f"cfl_{cfl}"
        require_success(run_athena(
            run_dir, path, f"time/integrator={integrator}", f"time/cfl_number={cfl}",
            "time/nlim=-1", "time/tlim=0.2", "mesh/nx1=128", "mesh/nx2=1",
            "mesh/nx3=1", "meshblock/nx1=128", "meshblock/nx2=1", "meshblock/nx3=1",
            f"{fluid}/rsolver={'hllc' if fluid == 'hydro' else 'hlld'}",
            "turb_driving/dt_update=1", "turb_driving/tcorr=1e6",
            "turb_driving/nlow=1", "turb_driving/nhigh=1", "turb_driving/spectrum=power_law",
            "turb_driving/sol_fraction=0", "output1/dcycle=1000000",
            "output2/dcycle=1000000", "output3/dcycle=1000000", f"output3/variable={fluid}_u",
            "output3/id=cons"))
        output = read_binary_blocks(max((run_dir / "bin").glob("*.cons.*.bin")),
                                    ("dens", "mom1"))
        assert output["time"] == pytest.approx(0.2, abs=1.0e-14)
        assert len(output["blocks"]) == 1
        state = output["blocks"][0]["data"]
        assert np.all(np.isfinite(state))
        states.append(state)
    differences = [np.sqrt(np.sum((left - right)**2, axis=(1, 2, 3)))
                   for left, right in zip(states[:-1], states[1:])]
    orders = np.log2(np.array(differences[:-1]) / np.array(differences[1:]))
    assert np.all(orders > 1.8), f"density/momentum temporal orders: {orders}"


def test_forcing_start_stop_and_ou_boundaries(tmp_path):
    """Both halves use the same forcing interval, including its end boundary."""
    start, stop, update = 0.013, 0.047, 0.011
    path = forcing_input(tmp_path)
    path.write_text(path.read_text().replace(
        "<turb_driving>", "<turb_driving>\ntdriv_start = 0.013\ntdriv_duration = 0.034"
    ).replace("<output2>", "<output2>\ndata_precision = real"))
    run_dir = tmp_path / "run"
    require_success(run_athena(
        run_dir, path, "time/integrator=rk3", "time/cfl_number=0.15", "time/nlim=-1",
        "time/tlim=0.065", "turb_driving/turb_flag=1", f"turb_driving/dt_update={update}",
        "output2/dcycle=1", "output3/dcycle=1000000"))
    history = np.loadtxt(next(run_dir.glob("*.hydro.hst")))
    times = np.unique(history[:, 0])
    for event in (start, start + update, start + 2*update, start + 3*update, stop):
        assert np.any(np.abs(times - event) < 1.0e-14), f"missed forcing event at {event}"
    expected_work = 0.1*np.clip(history[:, 0] - start, 0.0, stop - start)
    np.testing.assert_allclose(history[:, 6] - history[0, 6], expected_work, atol=3.0e-12)
    np.testing.assert_allclose(history[:, 3:6], 0.0, atol=3.0e-12)

    shapes = {}
    repeats = 0
    for filename in sorted((run_dir / "bin").glob("*.force.*.bin")):
        output = read_force_blocks(filename)
        if output["cycle"] == 0:
            continue
        left_time = times[output["cycle"] - 1]
        force = np.concatenate([block["force"].reshape(3, -1)
                                for block in output["blocks"]], axis=1)
        if left_time < start - 1.0e-14 or left_time >= stop - 1.0e-14:
            np.testing.assert_array_equal(force, 0.0)
            continue
        interval = int(np.floor((left_time - start + 1.0e-14) / update))
        # Remove normalization and the density-dependent constant offset.
        shape = force - np.mean(force, axis=1, keepdims=True)
        shape /= np.linalg.norm(shape)
        if interval in shapes:
            repeats += 1
            np.testing.assert_allclose(shape, shapes[interval], rtol=3.0e-12, atol=3.0e-13)
        else:
            shapes[interval] = shape
    assert set(shapes) == {0, 1, 2, 3}
    assert repeats > 0
    for interval in range(1, 4):
        assert np.linalg.norm(shapes[interval] - shapes[interval - 1]) > 1.0e-4


def test_final_time_at_ou_boundary_has_no_tiny_step_tail(tmp_path):
    """A clipped step retains unused fluid CFL margin and reaches the final event."""
    path = forcing_input(tmp_path, normalization="accel_rms", target=0.25)
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/nlim=-1", "time/tlim=0.1",
                               "turb_driving/dt_update=0.01"))
    history = np.loadtxt(next(run_dir.glob("*.hydro.hst")))
    times = np.unique(history[:, 0])
    assert times[-1] == pytest.approx(0.1, abs=1.0e-14)
    # The weak forcing needs at most two fluid steps per 0.01 forcing interval.
    assert len(times) - 1 <= 20, f"excess steps near the final forcing boundary: {times}"


def test_first_half_kick_respects_source_cfl(tmp_path):
    """A strong first kick cannot use the unforced initial fluid timestep."""
    acceleration, cfl, dx = 1.0e5, 0.3, 1.0/16.0
    path = forcing_input(tmp_path, normalization="accel_rms", target=acceleration)
    path.write_text(path.read_text().replace(
        "<turb_driving>", "<turb_driving>\nmin_kx = 1\nmax_kx = 1\nmin_ky = 0\n"
        "max_ky = 0\nmin_kz = 0\nmax_kz = 0").replace(
        "<output3>", "<output3>\ndata_precision = real"))
    run_dir = tmp_path / "run"
    require_success(run_athena(
        run_dir, path, "time/nlim=1", f"time/cfl_number={cfl}", "mesh/nx2=1", "mesh/nx3=1",
        "meshblock/nx1=16", "meshblock/nx2=1", "meshblock/nx3=1", "hydro/rsolver=hllc",
        "turb_driving/nlow=1", "turb_driving/nhigh=1", "turb_driving/spectrum=power_law",
        "turb_driving/sol_fraction=0", "turb_driving/dt_update=1", "turb_driving/tcorr=1e6",
        "output3/variable=hydro_u", "output3/id=cons"))
    history = np.loadtxt(next(run_dir.glob("*.hydro.hst")))
    step = history[-1, 0] - history[0, 0]
    sound_speed = np.sqrt(5.0/3.0)
    assert 0.0 < step < cfl*dx/sound_speed
    # The peak acceleration of this cosine is bounded by sqrt(2) times its RMS.
    assert step*(sound_speed + 0.5*np.sqrt(2.0)*acceleration*step) <= cfl*dx
    output = read_binary_blocks(max((run_dir / "bin").glob("*.cons.*.bin")),
                                ("dens", "mom1", "mom2", "mom3", "ener"))
    state = output["blocks"][0]["data"]
    assert np.all(np.isfinite(state))
    assert np.all(state[0] > 0.0)
    assert np.all(state[4] - 0.5*np.sum(state[1:4]**2, axis=0)/state[0] > 0.0)


@pytest.mark.parametrize("fluid", ["hydro", "mhd"])
def test_sts_preserves_forcing_power(tmp_path, fluid):
    """Viscous half sweeps preserve the forcing work budget around RK."""
    path = forcing_input(tmp_path, "ideal" if fluid == "hydro" else None,
                         "ideal" if fluid == "mhd" else None)
    path.write_text(path.read_text().replace(
        "<time>", "<time>\nsts_integrator = rkl2").replace(
        f"<{fluid}>", f"<{fluid}>\nviscosity = 1.0\nviscosity_integrator = sts"))
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "time/integrator=rk3"))
    assert_step_power(run_dir, fluid)


@pytest.mark.parametrize("normalization", ["edot", "accel_rms"])
@pytest.mark.parametrize("target", [-1.0, float("nan"), float("inf")])
def test_invalid_forcing_strength(tmp_path, normalization, target):
    path = forcing_input(tmp_path, normalization=normalization, target=target)
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0


def test_dedt_requires_explicit_normalization(tmp_path):
    path = forcing_input(tmp_path)
    path.write_text(path.read_text().replace("normalization  = edot\n", ""))
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "dedt requires normalization = edot" in result.stdout + result.stderr


def test_relativistic_two_fluid_forcing_is_rejected(tmp_path):
    path = forcing_input(tmp_path, "ideal", "ideal")
    path.write_text(path.read_text() + "\n<coord>\nspecial_rel = true\n")
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "requires Newtonian fluids" in result.stdout + result.stderr


def test_general_relativistic_forcing_is_rejected(tmp_path):
    path = forcing_input(tmp_path)
    path.write_text(path.read_text() + "\n<coord>\ngeneral_rel = true\nminkowski = true\n")
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "does not support general" in result.stdout + result.stderr
