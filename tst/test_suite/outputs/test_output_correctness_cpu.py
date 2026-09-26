"""Focused checks for derived output, PDF axes/weights, and coarsened headers."""

import os
from pathlib import Path
import struct
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis/python"))
from bin_convert import read_binary, read_coarsened_binary  # noqa: E402
from read_pdf import read_pdf  # noqa: E402


def run_case(tmp_path, outputs, flags=(), source="tst/inputs/lwave_mhd.athinput"):
    binary = Path(os.environ.get("ATHENA", "./athena")).resolve()
    input_file = tmp_path / "test.athinput"
    input_file.write_text((ROOT / source).read_text() + outputs)
    result = subprocess.run(
        [str(binary), "-i", str(input_file), *flags], cwd=tmp_path,
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def output(block, variable, file_type="bin"):
    return (f"\n<output{block}>\nfile_type={file_type}\nvariable={variable}\n"
            f"id={variable}\ndcycle=1\ndata_precision=real\n")


def grid(data, field):
    result = np.empty((data["Nx3"], data["Nx2"], data["Nx1"]))
    for location, values in zip(data["mb_logical"], data["mb_data"][field]):
        k, j, i = np.asarray(location[:3][::-1]) * np.asarray(values.shape)
        nz, ny, nx = values.shape
        result[k:k + nz, j:j + ny, i:i + nx] = values
    return result


def curl(fields, spacing):
    def derivative(field, axis):
        return (np.roll(field, -1, axis) - np.roll(field, 1, axis)) / (2 * spacing[axis])
    x, y, z = fields
    return np.array([derivative(z, 1) - derivative(y, 0),
                     derivative(x, 0) - derivative(z, 2),
                     derivative(y, 2) - derivative(x, 1)])


def test_curls_and_repeated_dumps(tmp_path):
    variables = ["mhd_w_bcc", "mhd_wz", "mhd_w2", "mhd_jz", "mhd_j2",
                 "mhd_k_jxb", "mhd_dynamo_ks"]
    run_case(tmp_path, "".join(output(i + 1, v) for i, v in enumerate(variables)), [
        "mesh/nx1=16", "mesh/nx2=8", "mesh/nx3=8", "meshblock/nx1=8",
        "meshblock/nx2=4", "meshblock/nx3=4", "mesh_refinement/refinement=none",
        "time/nlim=2", "problem/amp=0.001",
    ])
    snapshots = sorted(tmp_path.glob("bin/*.mhd_w_bcc.*.bin"))
    assert len(snapshots) >= 3
    for path in snapshots:
        primitives = read_binary(path)
        suffix = path.name.split(".")[-2]
        velocity = np.array([grid(primitives, f"vel{a}") for a in "xyz"])
        magnetic = np.array([grid(primitives, f"bcc{a}") for a in (1, 2, 3)])
        spacing = [1.5 / 8, 1.5 / 8, 3.0 / 16]
        omega = curl(velocity, spacing)
        current = curl(magnetic, spacing)
        expected = {"mhd_wz": ("vorz", omega[2]),
                    "mhd_w2": ("vor2", np.sum(omega**2, axis=0)),
                    "mhd_jz": ("jz", current[2]),
                    "mhd_j2": ("j2", np.sum(current**2, axis=0)),
                    "mhd_k_jxb": ("k_jxb", np.linalg.norm(
                        np.cross(current, magnetic, axisa=0, axisb=0, axisc=0), axis=0)
                        / np.sum(magnetic**2, axis=0))}
        assert np.max(np.abs(omega)) > 1e-5
        assert np.max(np.abs(current)) > 1e-5
        for variable, (field, reference) in expected.items():
            filename = next(tmp_path.glob(f"bin/*.{variable}.{suffix}.bin"))
            np.testing.assert_allclose(grid(read_binary(filename), field), reference,
                                       rtol=2e-10, atol=2e-12)
        # The rectangular blocks exercise the old k/i typo under bounds checking.
        dynamo = read_binary(next(tmp_path.glob(f"bin/*.mhd_dynamo_ks.{suffix}.bin")))
        assert all(np.all(np.isfinite(v)) for v in dynamo["mb_data"].values())


def test_derived_output_after_amr_growth(tmp_path):
    run_case(tmp_path, output(2, "mhd_bmag") + output(3, "mhd_w_bcc"), [
        "mesh/nx1=24", "mesh/nx2=24", "time/nlim=4",
        "mesh_refinement/max_nmb_per_rank=512",
    ], source="inputs/tests/divb_amr_2d.athinput")
    block_counts = []
    for path in sorted(tmp_path.glob("bin/*.mhd_bmag.*.bin")):
        data = read_binary(path)
        block_counts.append(data["n_mbs"])
        suffix = path.name.split(".")[-2]
        primitives = read_binary(next(tmp_path.glob(f"bin/*.mhd_w_bcc.{suffix}.bin")))
        reference = np.sqrt(sum(np.array(primitives["mb_data"][f"bcc{a}"])**2
                                for a in (1, 2, 3)))
        np.testing.assert_allclose(data["mb_data"]["bmag"], reference, rtol=2e-12)
    assert len(block_counts) >= 3
    assert max(block_counts) > block_counts[0], block_counts


@pytest.mark.parametrize("ndim", [1, 2, 3, 4])
def test_pdf_dimensions_and_variable_weight(tmp_path, ndim):
    variables = ["mhd_k_jxb", "mhd_bmag", "mhd_w_d", "temperature"]
    text = "\n<output1>\nfile_type=pdf\nid=pdf\ndcycle=1\n"
    for d in range(ndim):
        text += (f"variable_{d + 1}={variables[d]}\nnbin{d + 1}=1\n"
                 f"bin{d + 1}_min=-1\nbin{d + 1}_max=20\n")
    text += "weight=variable\nweight_variable=temperature\n"
    run_case(tmp_path, text, [
        "mesh/nx1=8", "mesh/nx2=4", "mesh/nx3=4", "meshblock/nx1=8",
        "meshblock/nx2=4", "meshblock/nx3=4", "mesh_refinement/refinement=none",
        "time/nlim=0", "problem/amp=0", "problem/dens=2",
    ])
    path = next(p for p in tmp_path.glob("pdf_*/*.pdf") if ".header." not in p.name)
    data = read_pdf(path)
    assert data["header"]["ndim"] == ndim
    expected_weight = (0.6 / (5 / 3 - 1) / 2) * (3 * 1.5 * 1.5)
    np.testing.assert_allclose(data["pdf"].sum(), expected_weight, rtol=1e-10)
    assert np.count_nonzero(data["pdf"]) == 1


def test_pdf_second_axis_requires_matching_physics(tmp_path):
    binary = Path(os.environ.get("ATHENA", "./athena")).resolve()
    text = (ROOT / "tst/inputs/lwave_hydro.athinput").read_text()
    text += """
<output1>
file_type=pdf
variable_1=hydro_w_d
variable_2=mhd_bmag
nbin1=1
nbin2=1
bin1_min=0
bin1_max=2
bin2_min=0
bin2_max=2
dcycle=1
"""
    source = tmp_path / "invalid.athinput"
    source.write_text(text)
    result = subprocess.run([str(binary), "-i", str(source),
                             "mesh_refinement/refinement=none", "time/nlim=0"],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0
    assert "no MHD object" in result.stdout + result.stderr


def test_legacy_pdf_axes_and_mixed_syntax(tmp_path):
    text = """
<output1>
file_type=pdf
variable=mhd_w_d
variable_2=mhd_bmag
nbin=1
nbin2=1
bin_min=0
bin_max=10
bin2_min=0
bin2_max=10
scale=linear
scale2=linear
dcycle=1
"""
    run_case(tmp_path, text, ["mesh/nx1=8", "mesh/nx2=8", "mesh/nx3=8",
                              "mesh_refinement/refinement=none", "time/nlim=0"])
    path = next(p for p in tmp_path.glob("pdf_*/*.pdf") if ".header." not in p.name)
    assert read_pdf(path)["header"]["ndim"] == 2
    # A third axis requires the N-D variable_1 syntax; never silently discard it.
    source = tmp_path / "test.athinput"
    source.write_text(source.read_text() + "variable_3=mhd_jz\n")
    binary = Path(os.environ.get("ATHENA", "./athena")).resolve()
    result = subprocess.run([str(binary), "-i", str(source), "time/nlim=0"],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0
    assert "variable_1" in result.stdout + result.stderr


def test_coarsened_header_preserves_equals(tmp_path):
    header = "<comment>\nnote = a=b\n<mesh>\nnghost=0\n"
    for axis in (1, 2, 3):
        header += f"nx{axis}=2\nx{axis}min=0\nx{axis}max=1\n"
    header += "<meshblock>\nnx1=2\nnx2=2\nnx3=2\n"
    data = ("Athena binary output version=1.1\nsize of preheader=8\ntime=0\n"
            "cycle=0\nsize of location=8\nsize of variable=8\n"
            "coarsening factor=2\nnumber of moments=1\nnote=a=b\n"
            f"nvars=1\nvariables: dens\nheader offset={len(header)}\n{header}").encode()
    data += struct.pack("=10i6dd", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0., 1., 0., 1., 0., 1., 3.25)
    path = tmp_path / "equals.cbin"
    path.write_bytes(data)
    result = read_coarsened_binary(path)
    assert result["n_mbs"] == 1
    np.testing.assert_array_equal(result["mb_data"]["dens"], [[[[3.25]]]])
