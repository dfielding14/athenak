"""Regression tests for stochastic passive-scalar forcing on CPU."""

import re
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


def read_history(path):
    """Read one AthenaK history file into a label-to-column dictionary."""
    with open(path, encoding="ascii") as history_file:
        history_file.readline()
        labels = re.findall(r"\[\d+\]=(\S+)", history_file.readline())
    values = np.loadtxt(path)
    if values.ndim == 1:
        values = values[np.newaxis, :]
    return {label: values[:, index] for index, label in enumerate(labels)}


def read_binary(path):
    """Read cell-centered variables and MeshBlock geometry from binary output."""
    with open(path, "rb") as file_obj:
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
        blocks = []
        while True:
            raw_indices = file_obj.read(6 * 4)
            if not raw_indices:
                break
            indices = struct.unpack("=6i", raw_indices)
            logical = struct.unpack("=4i", file_obj.read(4 * 4))
            limits = np.frombuffer(
                file_obj.read(6 * location_size), dtype=location_dtype
            ).copy()
            nx = indices[1] - indices[0] + 1
            ny = indices[3] - indices[2] + 1
            nz = indices[5] - indices[4] + 1
            count = nx * ny * nz
            data = np.frombuffer(
                file_obj.read(len(variables) * count * variable_size),
                dtype=variable_dtype,
            ).reshape((len(variables), nz, ny, nx)).copy()
            blocks.append({"logical": logical, "limits": limits, "data": data})
    return {"time": time, "cycle": cycle, "variables": variables, "blocks": blocks}


def latest_binary(run_dir, required_variable):
    """Return the latest binary output containing one requested variable."""
    candidates = []
    for path in (run_dir / "bin").glob("*.bin"):
        output = read_binary(path)
        if required_variable in output["variables"]:
            candidates.append((output["cycle"], path.name, path))
    assert candidates, f"no output containing {required_variable} under {run_dir}"
    return max(candidates)[2]


def assert_equal_outputs(left_path, right_path):
    """Check that two cell-centered outputs are bitwise identical."""
    left = read_binary(left_path)
    right = read_binary(right_path)
    assert left["time"] == right["time"]
    assert left["cycle"] == right["cycle"]
    assert left["variables"] == right["variables"]
    assert len(left["blocks"]) == len(right["blocks"])
    for left_block, right_block in zip(left["blocks"], right["blocks"]):
        assert left_block["logical"] == right_block["logical"]
        np.testing.assert_array_equal(left_block["limits"], right_block["limits"])
        np.testing.assert_array_equal(left_block["data"], right_block["data"])


def test_mean_projection_and_normalizations(tmp_path):
    """Both normalizations preserve the scalar mean and meet their target."""
    rms_dir = tmp_path / "source_rms"
    require_success(run_athena(rms_dir, "scalar_driving_source_rms.athinput"))
    history = read_history(rms_dir / "scalar_force_test.user.hst")
    active = history["time"] > 0.0
    rho = history["rho"][active]
    np.testing.assert_allclose(history["rth_s0"][active] / rho, 0.5, atol=2.0e-15)
    np.testing.assert_allclose(history["rf_s0"][active] / rho, 0.0, atol=2.0e-15)
    np.testing.assert_allclose(
        np.sqrt(history["rf2_s0"][active] / rho), 0.2, rtol=2.0e-14
    )

    variance_dir = tmp_path / "variance_rate"
    require_success(
        run_athena(variance_dir, "scalar_driving_variance_rate.athinput")
    )
    history = read_history(variance_dir / "scalar_variance_test.user.hst")
    rho = history["rho"]
    mean = history["rth_s0"] / rho
    variance = 0.5 * (history["rth2_s0"] / rho - mean * mean)
    final = np.flatnonzero(history["time"] > 0.0)[0]
    measured_rate = (variance[final] - variance[0]) / (
        history["time"][final] - history["time"][0]
    )
    assert measured_rate == pytest.approx(0.04, rel=2.0e-13)
    assert mean[final] == pytest.approx(mean[0], abs=2.0e-15)


def test_restart_preserves_scalar_and_modal_state(tmp_path):
    """A split run reproduces both the scalar field and rendered OU source."""
    input_name = "scalar_driving_source_rms.athinput"
    reference_dir = tmp_path / "reference"
    split_dir = tmp_path / "split"
    resume_dir = tmp_path / "resume"
    require_success(run_athena(reference_dir, input_name))
    require_success(
        run_athena(
            split_dir,
            input_name,
            "time/nlim=3",
            "output2/file_type=rst",
            "output2/dcycle=3",
            "output2/id=restart",
        )
    )
    restart = max((split_dir / "rst").glob("*.rst"))
    require_success(
        run_athena(
            resume_dir,
            input_name,
            "time/nlim=6",
            "output2/file_type=bin",
            "output2/variable=scalar_force",
            "output2/id=resume_force",
            "output2/dcycle=6",
            "output3/file_type=bin",
            "output3/variable=hydro_w",
            "output3/id=resume_scalar",
            "output3/dcycle=6",
            restart=restart,
        )
    )
    assert_equal_outputs(
        latest_binary(reference_dir, "scalar_force"),
        latest_binary(resume_dir, "scalar_force"),
    )
    assert_equal_outputs(
        latest_binary(reference_dir, "s_00"),
        latest_binary(resume_dir, "s_00"),
    )


def test_invalid_normalization_is_rejected(tmp_path):
    """A normalization cannot retain the target belonging to the other mode."""
    result = run_athena(
        tmp_path / "invalid",
        "scalar_driving_source_rms.athinput",
        "scalar_driving/normalization=variance_rate",
    )
    assert result.returncode != 0
    assert "normalization = variance_rate requires variance_rate" in (
        result.stdout + result.stderr
    )
