"""Regression tests for the turbulence driving source term on CPU."""

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
