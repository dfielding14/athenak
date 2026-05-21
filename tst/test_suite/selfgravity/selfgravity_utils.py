"""
Shared helpers for self-gravity regression tests.

The helpers intentionally inspect AthenaK binary outputs directly so the tests
can assert physical fields, not just process success.
"""

from pathlib import Path
import math
import re
import shutil
import struct
import subprocess

import numpy as np


DEFECT_RE = re.compile(r"Final defect norm = ([0-9.eE+-]+)")
REPO_ROOT = Path(__file__).resolve().parents[3]


def run_athena(
    input_file,
    basename,
    extra_args=None,
    max_final_defect=2.0e-8,
    expect_jeans_banner=False,
    mpi_ranks=None,
):
    command = []
    if mpi_ranks is not None:
        command = ["mpirun", "-np", str(mpi_ranks)]
    command.extend([
        "./athena",
        "-i",
        str(REPO_ROOT / input_file),
        f"job/basename={basename}",
        "gravity/show_defect=true",
    ])
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, check=False)
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "nan" not in combined.lower()
    matches = DEFECT_RE.findall(combined)
    assert matches, combined
    final_defect = float(matches[-1])
    if max_final_defect is not None:
        assert final_defect < max_final_defect
    if expect_jeans_banner:
        assert "k/k_J = 2" in combined
    return combined


def run_restart(restart_file, extra_args=None):
    command = ["./athena", "-r", str(restart_file)]
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, check=False)
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "nan" not in combined.lower()
    return combined


def expect_failure(input_file, extra_args=None, expected_text=None):
    command = ["./athena", "-i", str(REPO_ROOT / input_file)]
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, check=False)
    combined = result.stdout + result.stderr
    assert result.returncode != 0, combined
    if expected_text is not None:
        assert expected_text in combined
    return combined


def cleanup_outputs(*basenames):
    for directory in (Path("bin"), Path("rst"), Path("tab")):
        if not directory.exists():
            continue
        for basename in basenames:
            for path in directory.glob(f"{basename}*"):
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink(missing_ok=True)
        try:
            directory.rmdir()
        except OSError:
            pass
    for pattern in ("*.dat", "*.hst"):
        for path in Path(".").glob(pattern):
            path.unlink(missing_ok=True)


def latest_binary_output(basename, output_id=None):
    paths = sorted(Path("bin").glob(f"{basename}*.bin"))
    if output_id is not None:
        paths = [path for path in paths if output_id in path.name]
    assert paths, f"No binary outputs found for {basename}"
    return paths[-1]


def latest_restart_output(basename):
    paths = sorted(Path("rst").glob(f"{basename}*.rst"))
    assert paths, f"No restart outputs found for {basename}"
    return paths[-1]


def parse_binary_output(path):
    metadata = {}
    labels = []
    with open(path, "rb") as stream:
        while True:
            raw_line = stream.readline()
            assert raw_line, f"{path} ended before header offset"
            line = raw_line.decode("ascii")
            stripped = line.strip()
            if stripped.startswith("header offset="):
                header_offset = int(stripped.split("=", 1)[1])
                break
            if stripped.startswith("variables:"):
                labels = stripped.split(":", 1)[1].split()
                continue
            if "=" in stripped:
                key, value = stripped.split("=", 1)
                metadata[key.strip()] = value.strip()

        stream.seek(header_offset, 1)
        nvars = int(metadata["number of variables"])
        loc_size = int(metadata["size of location"])
        var_size = int(metadata["size of variable"])
        real_dtype = np.dtype("<f8" if loc_size == 8 else "<f4")
        var_dtype = np.dtype("<f4" if var_size == 4 else "<f8")

        blocks = []
        while True:
            index_bytes = stream.read(10 * 4)
            if not index_bytes:
                break
            assert len(index_bytes) == 10 * 4, f"Truncated block index in {path}"
            indices = struct.unpack("<10i", index_bytes)
            bounds_raw = stream.read(6 * loc_size)
            assert len(bounds_raw) == 6 * loc_size, f"Truncated block bounds in {path}"
            bounds = np.frombuffer(bounds_raw, dtype=real_dtype, count=6).astype(np.float64)
            ois, oie, ojs, oje, oks, oke, lx1, lx2, lx3, level = indices
            nx = oie - ois + 1
            ny = oje - ojs + 1
            nz = oke - oks + 1
            count = nvars * nx * ny * nz
            data_raw = stream.read(count * var_size)
            assert len(data_raw) == count * var_size, f"Truncated block data in {path}"
            data = np.frombuffer(data_raw, dtype=var_dtype, count=count)
            data = data.astype(np.float64).reshape((nvars, nz, ny, nx))
            blocks.append({
                "bounds": bounds,
                "indices": indices,
                "logical_location": (lx1, lx2, lx3),
                "level": level,
                "labels": labels,
                "data": data,
            })

    metadata["time"] = float(metadata.get("time", 0.0))
    metadata["cycle"] = int(metadata.get("cycle", 0))
    return {"metadata": metadata, "labels": labels, "blocks": blocks}


def block_centers(block):
    x1min, x1max, x2min, x2max, x3min, x3max = block["bounds"]
    nz, ny, nx = block["data"].shape[1:]
    x1 = np.linspace(x1min, x1max, nx, endpoint=False) + 0.5*(x1max - x1min)/nx
    x2 = np.linspace(x2min, x2max, ny, endpoint=False) + 0.5*(x2max - x2min)/ny
    x3 = np.linspace(x3min, x3max, nz, endpoint=False) + 0.5*(x3max - x3min)/nz
    x3g, x2g, x1g = np.meshgrid(x3, x2, x1, indexing="ij")
    return x1g, x2g, x3g


def variable_index(labels, name):
    if name in labels:
        return labels.index(name)
    matches = [i for i, label in enumerate(labels) if label.endswith(name)]
    assert len(matches) == 1, f"Could not resolve variable {name} in labels {labels}"
    return matches[0]


def concatenate_field(output, name):
    values = []
    x1_values = []
    x2_values = []
    x3_values = []
    for block in output["blocks"]:
        ivar = variable_index(block["labels"], name)
        x1, x2, x3 = block_centers(block)
        values.append(block["data"][ivar].ravel())
        x1_values.append(x1.ravel())
        x2_values.append(x2.ravel())
        x3_values.append(x3.ravel())
    return (
        np.concatenate(values),
        np.concatenate(x1_values),
        np.concatenate(x2_values),
        np.concatenate(x3_values),
    )


def uniform_grid_field(output, name):
    blocks = output["blocks"]
    x1_all = []
    x2_all = []
    x3_all = []
    for block in blocks:
        x1, x2, x3 = block_centers(block)
        x1_all.append(x1[0, 0, :])
        x2_all.append(x2[0, :, 0])
        x3_all.append(x3[:, 0, 0])
    x1_unique = np.unique(np.round(np.concatenate(x1_all), 14))
    x2_unique = np.unique(np.round(np.concatenate(x2_all), 14))
    x3_unique = np.unique(np.round(np.concatenate(x3_all), 14))
    grid = np.empty((len(x3_unique), len(x2_unique), len(x1_unique)))
    grid.fill(np.nan)
    for block in blocks:
        ivar = variable_index(block["labels"], name)
        x1, x2, x3 = block_centers(block)
        i1 = np.searchsorted(x1_unique, np.round(x1[0, 0, :], 14))
        i2 = np.searchsorted(x2_unique, np.round(x2[0, :, 0], 14))
        i3 = np.searchsorted(x3_unique, np.round(x3[:, 0, 0], 14))
        grid[np.ix_(i3, i2, i1)] = block["data"][ivar]
    assert not np.isnan(grid).any()
    return grid, x1_unique, x2_unique, x3_unique


def jeans_geometry(lengths=(1.0, 1.0, 1.0), n_jeans=0.5, rho0=1.0, cs=1.0):
    lx1, lx2, lx3 = lengths
    ang3 = math.atan(lx1/lx2)
    sin_a3 = math.sin(ang3)
    cos_a3 = math.cos(ang3)
    ang2 = math.atan(0.5*(lx1*cos_a3 + lx2*sin_a3)/lx3)
    sin_a2 = math.sin(ang2)
    cos_a2 = math.cos(ang2)
    wavelength = min(lx1*cos_a2*cos_a3, lx2*cos_a2*sin_a3, lx3*sin_a2)
    k_wave = 2.0*math.pi/wavelength
    four_pi_g = (n_jeans*k_wave*cs)**2/rho0
    direction = np.array([cos_a2*cos_a3, cos_a2*sin_a3, sin_a2])
    return k_wave, four_pi_g, direction


def jeans_phase(x1, x2, x3):
    _, _, direction = jeans_geometry()
    return direction[0]*x1 + direction[1]*x2 + direction[2]*x3


def relative_l2(error, reference):
    return np.linalg.norm(error.ravel()) / max(
        np.linalg.norm(reference.ravel()), 1.0e-300)
