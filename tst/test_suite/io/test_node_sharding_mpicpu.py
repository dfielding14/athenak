"""MPI regression coverage for per-node diagnostic and restart output."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402
from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


INPUT_FILE = "inputs/io_node_sharding.athinput"
NODE_OVERRIDES = tuple(
    f"output{number}/single_file_per_node=true" for number in range(1, 7)
)


def _run(tmp_path: Path, name: str, *overrides: str):
    run_dir = tmp_path / name
    run_dir.mkdir()
    proc = subprocess.run(
        ["mpirun", "-np", "2", "./athena", "-i", INPUT_FILE, "-d", str(run_dir), *overrides],
        check=True,
        capture_output=True,
        text=True,
    )
    return run_dir, proc.stdout


def _compare_binary(shared: Path, node: Path, coarsened=False):
    reader = bin_convert.read_coarsened_binary if coarsened else bin_convert.read_binary
    reference = reader(str(shared))
    reconstructed = reader(str(node), assemble_shards=True)
    assert reference["var_names"] == reconstructed["var_names"]
    assert reference["n_mbs"] == reconstructed["n_mbs"]
    np.testing.assert_array_equal(reference["mb_logical"], reconstructed["mb_logical"])
    for variable in reference["var_names"]:
        np.testing.assert_allclose(
            reference["mb_data"][variable], reconstructed["mb_data"][variable]
        )


def test_node_sharded_diagnostics_reconstruct_to_shared_output(tmp_path):
    shared, _ = _run(tmp_path, "shared")
    node, _ = _run(tmp_path, "node", *NODE_OVERRIDES)
    node_dir = node / "bin" / "node_00000000"
    assert node_dir.is_dir()
    assert len(list((node / "bin").glob("node_*"))) == 1

    _compare_binary(
        shared / "bin" / "io_node.full.00000.bin",
        node_dir / "io_node.full.00000.bin",
    )
    _compare_binary(
        shared / "bin" / "io_node.slice.00000.bin",
        node_dir / "io_node.slice.00000.bin",
    )
    _compare_binary(
        shared / "cbin_coarse_2" / "io_node.coarse.00000.cbin",
        node / "cbin_coarse_2" / "node_00000000" / "io_node.coarse.00000.cbin",
        coarsened=True,
    )
    shared_pdf = read_pdf(
        str(shared / "pdf_node_pdf_hydro_w_d" / "io_node.00000.pdf")
    )
    node_pdf = read_pdf(
        str(
            node
            / "pdf_node_pdf_hydro_w_d"
            / "node_00000000"
            / "io_node.00000.pdf"
        )
    )
    np.testing.assert_allclose(node_pdf["pdf"], shared_pdf["pdf"])
    shared_surface = read_sphslice(
        str(shared / "bin" / "io_node.density.r_0.25.00000.sph.bin")
    )
    node_surface = read_sphslice(
        str(node_dir / "io_node.density.r_0.25.00000.sph.bin")
    )
    np.testing.assert_allclose(node_surface["data"], shared_surface["data"])

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "bin_convert.py"),
            "--assemble-shards",
            str(node_dir / "io_node.full.00000.bin"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert (node_dir / "io_node.full.00000.athdf").exists()


def test_node_restart_manifest_resumes_without_overwriting_terminal_checkpoint(tmp_path):
    run_dir, stdout = _run(
        tmp_path,
        "resume",
        *NODE_OVERRIDES,
        "time/output_timing=true",
        "time/final_output_policy=restart_only",
    )
    assert "event=initial block=output1 type=bin distribution=node elapsed_max_s=" in stdout
    assert "event=final block=output6 type=rst distribution=node elapsed_max_s=" in stdout
    terminal = run_dir / "rst" / "io_node.00001.rst"
    original_manifest = terminal.read_bytes()
    assert b"AthenaK node restart manifest version=1" in original_manifest
    assert list((run_dir / "rst" / "node_00000000").glob("*.payload.rst"))
    assert not list((run_dir / "rst").rglob("*.tmp"))

    subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-r",
            str(terminal),
            "-d",
            str(run_dir),
            "time/tlim=0.02",
            "time/nlim=2",
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert terminal.read_bytes() == original_manifest
    assert (run_dir / "rst" / "io_node.00002.rst").exists()


@pytest.mark.parametrize(
    ("corruption", "expected"),
    (
        ("traversal", "payload path"),
        ("absolute", "payload path"),
        ("incomplete", "completion record"),
        ("byte_count", "payload block count or byte count"),
    ),
)
def test_node_restart_rejects_corrupted_manifest(tmp_path, corruption, expected):
    run_dir, _ = _run(tmp_path, "corrupt", *NODE_OVERRIDES)
    manifest = run_dir / "rst" / "io_node.00000.rst"
    text = manifest.read_text()
    if corruption == "traversal":
        text = text.replace("node_00000000/", "../node_00000000/", 1)
    elif corruption == "absolute":
        text = text.replace("node_00000000/", "/tmp/node_00000000/", 1)
    elif corruption == "incomplete":
        text = text.replace("complete=1", "complete=0", 1)
    else:
        lines = text.splitlines()
        for index, line in enumerate(lines):
            if line.startswith("payload "):
                fields = line.split()
                fields[3] = str(int(fields[3]) + 1)
                lines[index] = " ".join(fields)
                break
        text = "\n".join(lines) + "\n"
    manifest.write_text(text)
    proc = subprocess.run(
        ["mpirun", "-np", "2", "./athena", "-r", str(manifest), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)


def test_existing_per_rank_restart_resumes_and_timing_identifies_layout(tmp_path):
    run_dir, stdout = _run(
        tmp_path,
        "per_rank",
        "output1/single_file_per_rank=true",
        "output6/single_file_per_rank=true",
        "time/output_timing=true",
        "time/final_output_policy=restart_only",
    )
    assert "block=output1 type=bin distribution=rank elapsed_max_s=" in stdout
    assert "block=output6 type=rst distribution=rank elapsed_max_s=" in stdout
    rank0 = run_dir / "rst" / "rank_00000000" / "io_node.00001.rst"
    rank1 = run_dir / "rst" / "rank_00000001" / "io_node.00001.rst"
    assert rank0.exists()
    assert rank1.exists()
    subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-r",
            str(rank0),
            "-d",
            str(run_dir),
            "time/tlim=0.02",
            "time/nlim=2",
            "time/final_output_policy=none",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def test_conflicting_rank_and_node_modes_are_rejected(tmp_path):
    run_dir = tmp_path / "conflict"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "output1/single_file_per_rank=true",
            "output1/single_file_per_node=true",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "cannot set both single_file_per_rank=true and single_file_per_node=true" in (
        proc.stdout + proc.stderr
    )


def test_promoted_node_example_generates_readable_outputs_and_manifest(tmp_path):
    run_dir = tmp_path / "node_example"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            str(ROOT / "inputs" / "io" / "node_sharded_outputs.athinput"),
            "-d",
            str(run_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    node_bin = run_dir / "bin" / "node_00000000" / "io_node_example.density.00000.bin"
    node_cbin = (
        run_dir
        / "cbin_coarse_density_2"
        / "node_00000000"
        / "io_node_example.coarse_density.00000.cbin"
    )
    node_pdf = (
        run_dir
        / "pdf_radius_density_hydro_w_d"
        / "node_00000000"
        / "io_node_example.00000.pdf"
    )
    node_surface = (
        run_dir / "bin" / "node_00000000" / "io_node_example.density.r_0.25.00000.sph.bin"
    )
    assert bin_convert.read_binary(str(node_bin), assemble_shards=True)["n_mbs"] == 4
    assert bin_convert.read_coarsened_binary(str(node_cbin), assemble_shards=True)["n_mbs"] == 4
    assert read_pdf(str(node_pdf))["header"]["distribution"] == "node"
    assert read_sphslice(str(node_surface))["data"].shape == (16, 32, 1)
    assert (run_dir / "rst" / "io_node_example.00001.rst").exists()
    assert "type=bin distribution=node elapsed_max_s=" in proc.stdout
    assert "type=rst distribution=node elapsed_max_s=" in proc.stdout

    summary = subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"),
            "bin",
            str(node_bin),
            "--assemble-shards",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "binary meshblocks=4" in summary.stdout
