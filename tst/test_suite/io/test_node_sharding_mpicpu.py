"""MPI regression coverage for per-node diagnostic and restart output."""

import os
from pathlib import Path
import shutil
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
MAX_MPI_BYTES_ENV = "ATHENAK_TEST_MAX_MPI_BYTES"
NODE_OVERRIDES = tuple(
    f"output{number}/single_file_per_node=true" for number in range(1, 7)
)


def _env(max_mpi_bytes: str):
    env = os.environ.copy()
    env[MAX_MPI_BYTES_ENV] = max_mpi_bytes
    return env


def _run(tmp_path: Path, name: str, *overrides: str, env=None, nranks=2):
    run_dir = tmp_path / name
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "mpirun",
            "-np",
            str(nranks),
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            *overrides,
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
        timeout=90,
    )
    return run_dir, proc.stdout


def _resume(run_dir: Path, restart: Path, *, env=None, nranks=2, check=True):
    run_dir.mkdir()
    return subprocess.run(
        [
            "mpirun",
            "-np",
            str(nranks),
            "./athena",
            "-r",
            str(restart),
            "-d",
            str(run_dir),
            "time/tlim=0.02",
            "time/nlim=2",
            "time/final_output_policy=none",
        ],
        check=check,
        capture_output=True,
        text=True,
        env=env,
        timeout=90,
    )


def _copy_node_checkpoint(tmp_path: Path, template: Path, name: str):
    run_dir = tmp_path / name
    shutil.copytree(template, run_dir)
    return run_dir, run_dir / "rst" / "io_node.00000.rst"


def _payload_path(manifest: Path):
    payload_record = next(
        line for line in manifest.read_text().splitlines() if line.startswith("payload ")
    )
    return manifest.parent / payload_record.split()[4]


@pytest.fixture(scope="module")
def node_restart_template(tmp_path_factory):
    template_root = tmp_path_factory.mktemp("node_restart_template")
    run_dir, _ = _run(template_root, "checkpoint", *NODE_OVERRIDES)
    return run_dir


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
    assert (
        "event=initial block=output1 type=bin distribution=node elapsed_max_s=" in stdout
    )
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
    assert not list((run_dir / "rst").rglob("*.assembled"))


@pytest.mark.parametrize(
    ("corruption", "expected"),
    (
        ("traversal", "payload path"),
        ("absolute", "payload path"),
        ("incomplete", "completion record"),
        ("byte_count", "payload block count or byte count"),
        ("missing_record", "expected 'data_size' record"),
        ("duplicate_record", "expected 'header_size' record"),
        ("reordered_record", "expected 'header_size' record"),
        ("unknown_record", "unrecognized or malformed inventory record"),
        ("trailing_record", "records found after the manifest terminator"),
        ("unsupported_version", "unsupported manifest signature"),
        ("mixed_generation", "mixed generations"),
        ("header_mismatch", "does not match canonical payload 0"),
        ("bad_node_directory", "declared node directory"),
        ("segment_gap", "leave a gap"),
        ("segment_overlap", "overlap"),
        ("segment_invalid_node", "invalid node"),
        ("segment_local_offset", "node-local payload offset"),
        ("segment_count", "node-local payload count"),
    ),
)
def test_node_restart_rejects_corrupted_manifest(
    tmp_path, node_restart_template, corruption, expected
):
    run_dir, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "corrupt")
    text = manifest.read_text()
    if corruption == "traversal":
        text = text.replace("node_00000000/", "../node_00000000/", 1)
    elif corruption == "absolute":
        text = text.replace("node_00000000/", "/tmp/node_00000000/", 1)
    elif corruption == "incomplete":
        text = text.replace("complete=1", "complete=0", 1)
    elif corruption == "byte_count":
        lines = text.splitlines()
        for index, line in enumerate(lines):
            if line.startswith("payload "):
                fields = line.split()
                fields[3] = str(int(fields[3]) + 1)
                lines[index] = " ".join(fields)
                break
        text = "\n".join(lines) + "\n"
    elif corruption == "missing_record":
        text = text.replace(next(line for line in text.splitlines()
                                 if line.startswith("data_size=")) + "\n", "", 1)
    elif corruption == "duplicate_record":
        record = next(line for line in text.splitlines() if line.startswith("nmb_total="))
        text = text.replace(record + "\n", record + "\n" + record + "\n", 1)
    elif corruption == "reordered_record":
        lines = text.splitlines()
        header = next(index for index, line in enumerate(lines)
                      if line.startswith("header_size="))
        data = next(index for index, line in enumerate(lines)
                    if line.startswith("data_size="))
        lines[header], lines[data] = lines[data], lines[header]
        text = "\n".join(lines) + "\n"
    elif corruption == "unknown_record":
        text = text.replace("end\n", "unknown=1\nend\n", 1)
    elif corruption == "trailing_record":
        text += "unknown=1\n"
    elif corruption == "unsupported_version":
        text = text.replace("version=1", "version=2", 1)
    elif corruption == "mixed_generation":
        lines = text.splitlines()
        payload_index = next(index for index, line in enumerate(lines)
                             if line.startswith("payload "))
        fields = lines[payload_index].split()
        header_size = next(line.split("=")[1] for line in lines
                           if line.startswith("header_size="))
        fields[1] = "1"
        fields[2] = "0"
        fields[3] = header_size
        fields[4] = fields[4].replace("node_00000000/", "node_00000001/", 1)
        fields[4] = fields[4].replace(".payload.rst", "1.payload.rst", 1)
        lines.insert(payload_index + 1, " ".join(fields))
        text = "\n".join(lines).replace("payload_count=1", "payload_count=2", 1) + "\n"
    elif corruption == "header_mismatch":
        lines = text.splitlines()
        payload_index = next(index for index, line in enumerate(lines)
                             if line.startswith("payload "))
        fields = lines[payload_index].split()
        header_size = next(int(line.split("=")[1]) for line in lines
                           if line.startswith("header_size="))
        fields[1] = "1"
        fields[2] = "0"
        fields[3] = str(header_size)
        fields[4] = fields[4].replace("node_00000000/", "node_00000001/", 1)
        duplicate = manifest.parent / fields[4]
        duplicate.parent.mkdir()
        header = bytearray(_payload_path(manifest).read_bytes()[:header_size])
        header[0] ^= 1
        duplicate.write_bytes(header)
        lines.insert(payload_index + 1, " ".join(fields))
        text = "\n".join(lines).replace("payload_count=1", "payload_count=2", 1) + "\n"
    elif corruption == "bad_node_directory":
        text = text.replace("node_00000000/", "node_00000001/", 1)
    else:
        lines = text.splitlines()
        segment_indexes = [
            index for index, line in enumerate(lines) if line.startswith("segment ")
        ]
        fields = lines[segment_indexes[0]].split()
        if corruption == "segment_gap":
            fields[2] = str(int(fields[2]) + 1)
        elif corruption == "segment_invalid_node":
            fields[1] = "1"
        elif corruption == "segment_local_offset":
            fields[4] = str(int(fields[4]) + 1)
        elif corruption == "segment_count":
            payload = next(line.split() for line in lines if line.startswith("payload "))
            fields[3] = str(int(payload[2]) + 1)
        else:
            assert corruption == "segment_overlap"
            fields = lines[segment_indexes[1]].split()
            fields[2] = str(int(fields[2]) - 1)
            lines[segment_indexes[1]] = " ".join(fields)
            text = "\n".join(lines) + "\n"
            manifest.write_text(text)
            fields = None
        if fields is not None:
            lines[segment_indexes[0]] = " ".join(fields)
            text = "\n".join(lines) + "\n"
    manifest.write_text(text)
    proc = subprocess.run(
        ["mpirun", "-np", "2", "./athena", "-r", str(manifest), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)


@pytest.mark.parametrize("damage", ("missing", "truncated"))
def test_node_restart_rejects_missing_or_truncated_payload(
    tmp_path, node_restart_template, damage
):
    run_dir, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, damage)
    payload = _payload_path(manifest)
    if damage == "missing":
        payload.unlink()
    else:
        payload.write_bytes(payload.read_bytes()[:-1])
    proc = _resume(tmp_path / f"{damage}_resume", manifest, check=False)
    assert proc.returncode != 0
    assert "is absent or incomplete" in (proc.stdout + proc.stderr)


@pytest.mark.parametrize(
    "alias",
    ("exact", "dot", "repeated_slash", "symlink", "temporary"),
)
def test_node_restart_rejects_payload_path_entry(tmp_path, node_restart_template, alias):
    run_dir, manifest = _copy_node_checkpoint(
        tmp_path, node_restart_template, f"payload_entry_{alias}"
    )
    payload = _payload_path(manifest)
    if alias == "exact":
        restart = payload
    elif alias == "dot":
        restart = str(payload.parent) + "/./" + payload.name
    elif alias == "repeated_slash":
        restart = str(payload.parent) + "//" + payload.name
    elif alias == "symlink":
        restart = tmp_path / "payload_alias.rst"
        restart.symlink_to(payload)
    else:
        restart = Path(str(payload) + ".tmp")
        shutil.copyfile(payload, restart)
    proc = _resume(
        tmp_path / f"payload_entry_resume_{alias}", restart, check=False
    )
    assert proc.returncode != 0
    assert "use the public manifest path" in (proc.stdout + proc.stderr)


def test_unrelated_shared_restart_payload_suffix_remains_compatible(tmp_path):
    run_dir, _ = _run(tmp_path, "shared_payload_suffix")
    shared_restart = run_dir / "rst" / "io_node.00000.rst"
    renamed = run_dir / "rst" / "ordinary.payload.rst"
    shutil.copyfile(shared_restart, renamed)
    _resume(tmp_path / "shared_payload_suffix_resume", renamed)


def test_node_restart_rejects_payload_symlink_escape(tmp_path, node_restart_template):
    _, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "symlink_escape")
    payload = _payload_path(manifest)
    escaped = tmp_path / "escaped.payload.rst"
    shutil.copyfile(payload, escaped)
    payload.unlink()
    payload.symlink_to(escaped)
    proc = _resume(tmp_path / "symlink_escape_resume", manifest, check=False)
    assert proc.returncode != 0
    assert "escapes the manifest directory through a symlink" in (
        proc.stdout + proc.stderr
    )


def test_node_restart_never_creates_assembled_files(tmp_path, node_restart_template):
    _, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "native")
    assembled = Path(str(manifest) + ".assembled")
    temporary = Path(str(manifest) + ".assembled.tmp")
    assembled.mkdir()
    temporary.mkdir()
    _resume(tmp_path / "native_resume", manifest)
    assert assembled.is_dir()
    assert temporary.is_dir()


def test_stale_node_files_do_not_change_shared_restart_interpretation(tmp_path):
    run_dir, _ = _run(tmp_path, "shared")
    shared_restart = run_dir / "rst" / "io_node.00000.rst"
    stale = run_dir / "rst" / "node_00000000" / "stale.payload.rst"
    stale.parent.mkdir()
    stale.write_bytes(b"not a restart")
    _resume(tmp_path / "shared_resume", shared_restart)


def test_forced_small_chunk_node_restart_write_and_read(tmp_path):
    run_dir, _ = _run(tmp_path, "forced_chunks", *NODE_OVERRIDES, env=_env("7"))
    manifest = run_dir / "rst" / "io_node.00000.rst"
    _resume(tmp_path / "forced_chunks_resume", manifest, env=_env("7"))
    assert not list((run_dir / "rst").rglob("*.assembled"))


def test_node_restart_supports_rank_count_change_on_local_node(
    tmp_path, node_restart_template
):
    _, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "rank_count")
    # This workstation exercises redistribution across ranks on one node. True
    # multi-node payload routing and non-owning-node qualification remain external gates.
    _resume(tmp_path / "rank_count_resume", manifest, nranks=1)


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
        run_dir
        / "bin"
        / "node_00000000"
        / "io_node_example.density.r_0.25.00000.sph.bin"
    )
    assert bin_convert.read_binary(str(node_bin), assemble_shards=True)["n_mbs"] == 4
    assert (
        bin_convert.read_coarsened_binary(str(node_cbin), assemble_shards=True)["n_mbs"]
        == 4
    )
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
