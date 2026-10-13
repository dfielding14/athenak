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
    return subprocess.run(command, capture_output=True, text=True,
                          check=False, timeout=60)


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
    """One complete kick supplies dedt times box volume on every timestep."""
    if fluid == "mhd" and integrator == "rk4":
        pytest.skip("Existing MHD RK4 register updates do not preserve a uniform state")
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


def test_edot_default_matches_explicit(tmp_path):
    """This branch retains fixed power as the default normalization."""
    explicit = forcing_input(tmp_path)
    implicit = tmp_path / "default.athinput"
    implicit.write_text(explicit.read_text().replace("normalization  = edot\n", ""))
    for label, path in (("explicit", explicit), ("default", implicit)):
        require_success(run_athena(tmp_path / label, path))
        assert_step_power(tmp_path / label)
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


@pytest.mark.parametrize("integrator", ["imex2", "imex2+", "imex3"])
def test_two_fluid_step_power(tmp_path, integrator):
    """Both fluids finish their implicit stages before receiving one common kick."""
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


@pytest.mark.parametrize("normalization", ["edot", "accel_rms"])
@pytest.mark.parametrize("target", [-1.0, float("nan"), float("inf")])
def test_invalid_forcing_strength(tmp_path, normalization, target):
    path = forcing_input(tmp_path, normalization=normalization, target=target)
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0


@pytest.mark.parametrize("normalization", ["edot", "accel_rms"])
def test_nonfinite_density_is_not_reported_as_zero_force(tmp_path, normalization):
    """A broken fluid state must be distinguished from a vanishing forcing field."""
    path = forcing_input(tmp_path, normalization=normalization)
    path.write_text(path.read_text().replace(
        "<problem>", "<problem>\nnonfinite_density = true"))
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "non-finite mass or acceleration integral" in result.stdout + result.stderr


def test_nonfinite_work_is_rejected(tmp_path):
    """Fixed-power normalization cannot propagate a non-finite momentum integral."""
    path = forcing_input(tmp_path)
    path.write_text(path.read_text().replace(
        "<problem>", "<problem>\nnormalization_velocity = nan"))
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "non-finite forcing normalization integral" in result.stdout + result.stderr


def test_driving_window_uses_cycle_start_time(tmp_path):
    """Start/stop times select full kicks; they do not clip the fluid timestep."""
    path = forcing_input(tmp_path)
    path.write_text(path.read_text().replace(
        "<turb_driving>", "<turb_driving>\ntdriv_start = 0.02\ntdriv_duration = 0.025"))
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir, path, "turb_driving/turb_flag=1"))
    data = np.loadtxt(next(run_dir.glob("*.hydro.hst")))
    dt = np.diff(data[:, 0])
    driven = (data[:-1, 0] >= 0.02) & (data[:-1, 0] < 0.045)
    assert np.any(driven) and np.any(~driven)
    assert np.any(driven & (data[1:, 0] > 0.045))
    np.testing.assert_allclose(np.diff(data[:, 6]), 0.1*dt*driven,
                               rtol=2.0e-10, atol=3.0e-12)


def test_accel_rms_requires_explicit_normalization(tmp_path):
    path = forcing_input(tmp_path, normalization="accel_rms")
    path.write_text(path.read_text().replace("normalization  = accel_rms\n", ""))
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "normalization = edot requires dedt" in result.stdout + result.stderr


def test_relativistic_two_fluid_forcing_is_rejected(tmp_path):
    path = forcing_input(tmp_path, "ideal", "ideal")
    path.write_text(path.read_text() + "\n<coord>\nspecial_rel = true\n")
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "requires Newtonian fluids" in result.stdout + result.stderr


def test_general_relativistic_forcing_is_rejected(tmp_path):
    path = forcing_input(tmp_path)
    path.write_text(path.read_text() +
                    "\n<coord>\ngeneral_rel = true\nminkowski = true\n")
    result = run_athena(tmp_path / "run", path)
    assert result.returncode != 0
    assert "does not support general" in result.stdout + result.stderr


@pytest.mark.parametrize("integrator", ["rk2", "rk3"])
def test_sparse_annulus_step_power(tmp_path, integrator):
    """The SGS branch's 2D sparse modes retain exact per-cycle power."""
    path = forcing_input(tmp_path)
    path.write_text(path.read_text().replace(
        "<turb_driving>", "<turb_driving>\nmin_kz = 0\nmax_kz = 0\n"
        "mode_sampling = sparse_annulus\nsparse_mode_count = 8"))
    run_dir = tmp_path / "run"
    require_success(run_athena(
        run_dir, path, f"time/integrator={integrator}",
        "mesh/nx3=1", "meshblock/nx3=1"))
    assert_step_power(run_dir)
