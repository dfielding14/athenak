"""MPI regression coverage for per-node diagnostic and restart output."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402
from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


INPUT_FILE = "inputs/io_node_sharding.athinput"
MAX_MPI_BYTES_ENV = "ATHENAK_TEST_MAX_MPI_BYTES"
NODE_PAYLOAD_MARKER = b"AthenaK node restart payload version=1\n"
MAX_NODE_MANIFEST_BYTES = 64 * 1024 * 1024
MAX_GENERATED_PAYLOAD_PATH_BYTES = 1024
OVERSIZED_MANIFEST_LINE = "x" * 4096
NODE_OVERRIDES = tuple(
    f"output{number}/single_file_per_node=true" for number in range(1, 7)
)
RESTART_ONLY_NODE_OVERRIDES = tuple(
    f"output{number}/dt=-1" for number in range(1, 6)
) + (
    "output6/single_file_per_node=true",
    "time/final_output_policy=none",
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


def _split_manifest_into_two_payloads(manifest: Path):
    lines = manifest.read_text().splitlines()
    header_size = next(
        int(line.split("=", 1)[1]) for line in lines
        if line.startswith("header_size=")
    )
    data_size = next(
        int(line.split("=", 1)[1]) for line in lines
        if line.startswith("data_size=")
    )
    nmb_total = next(
        int(line.split("=", 1)[1]) for line in lines
        if line.startswith("nmb_total=")
    )
    first_blocks = nmb_total // 2
    second_blocks = nmb_total - first_blocks
    assert first_blocks > 0
    assert second_blocks > 0

    payload_index = next(
        index for index, line in enumerate(lines) if line.startswith("payload ")
    )
    original_fields = lines[payload_index].split()
    original_payload = manifest.parent / original_fields[4]
    original_data = original_payload.read_bytes()
    assert len(original_data) == header_size + nmb_total * data_size

    second_relative = (
        "node_00000001/" + Path(original_fields[4]).name
    )
    second_payload = manifest.parent / second_relative
    second_payload.parent.mkdir()
    split_offset = header_size + first_blocks * data_size
    second_payload.write_bytes(original_data[:header_size] + original_data[split_offset:])
    original_payload.write_bytes(original_data[:split_offset])

    payload_records = (
        f"payload 0 {first_blocks} {header_size + first_blocks * data_size} "
        f"{original_fields[4]}",
        f"payload 1 {second_blocks} {header_size + second_blocks * data_size} "
        f"{second_relative}",
    )
    segment_records = (
        f"segment 0 0 {first_blocks} 0",
        f"segment 1 {first_blocks} {second_blocks} 0",
    )
    end_index = lines.index("end")
    lines[payload_index:end_index] = payload_records + segment_records
    payload_count_index = next(
        index for index, line in enumerate(lines)
        if line.startswith("payload_count=")
    )
    lines[payload_count_index] = "payload_count=2"
    manifest.write_text("\n".join(lines) + "\n")


def _interleave_manifest_across_two_payloads(manifest: Path):
    lines = manifest.read_text().splitlines()
    header_size = next(
        int(line.split("=", 1)[1]) for line in lines
        if line.startswith("header_size=")
    )
    data_size = next(
        int(line.split("=", 1)[1]) for line in lines
        if line.startswith("data_size=")
    )
    nmb_total = next(
        int(line.split("=", 1)[1]) for line in lines
        if line.startswith("nmb_total=")
    )
    assert nmb_total >= 2
    payload_index = next(
        index for index, line in enumerate(lines) if line.startswith("payload ")
    )
    original_fields = lines[payload_index].split()
    original_payload = manifest.parent / original_fields[4]
    original_data = original_payload.read_bytes()
    assert len(original_data) == header_size + nmb_total * data_size
    blocks = [
        original_data[
            header_size + index * data_size:header_size + (index + 1) * data_size
        ]
        for index in range(nmb_total)
    ]
    second_relative = "node_00000001/" + Path(original_fields[4]).name
    second_payload = manifest.parent / second_relative
    second_payload.parent.mkdir(exist_ok=True)
    payload_blocks = (blocks[::2], blocks[1::2])
    original_payload.write_bytes(
        original_data[:header_size] + b"".join(payload_blocks[0])
    )
    second_payload.write_bytes(original_data[:header_size] + b"".join(payload_blocks[1]))
    payload_records = (
        f"payload 0 {len(payload_blocks[0])} "
        f"{header_size + len(payload_blocks[0]) * data_size} {original_fields[4]}",
        f"payload 1 {len(payload_blocks[1])} "
        f"{header_size + len(payload_blocks[1]) * data_size} {second_relative}",
    )
    segment_records = tuple(
        f"segment {payload_id} {gid} 1 {payload_offset}"
        for gid in range(nmb_total)
        for payload_id, payload_offset in ((gid % 2, gid // 2),)
    )
    end_index = lines.index("end")
    lines[payload_index:end_index] = payload_records + segment_records
    payload_count_index = next(
        index for index, line in enumerate(lines)
        if line.startswith("payload_count=")
    )
    lines[payload_count_index] = "payload_count=2"
    manifest.write_text("\n".join(lines) + "\n")


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
        str(
            shared / "bin"
            / "io_node.density.r_2.5000000000000000e-01.00000.sph.bin"
        )
    )
    node_surface = read_sphslice(
        str(node_dir / "io_node.density.r_2.5000000000000000e-01.00000.sph.bin")
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
        timeout=90,
    )
    assert (node_dir / "io_node.full.00000.athdf").exists()


def test_node_restart_manifest_resumes_without_overwriting_terminal_checkpoint(tmp_path):
    run_dir, _ = _run(
        tmp_path,
        "resume",
        *NODE_OVERRIDES,
        "time/final_output_policy=restart_only",
    )
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
        timeout=90,
    )
    assert terminal.read_bytes() == original_manifest
    assert (run_dir / "rst" / "io_node.00002.rst").exists()
    assert not list((run_dir / "rst").rglob("*.assembled"))


def test_node_restart_generation_collision_preserves_ambient_reservation(tmp_path):
    run_dir = tmp_path / "generation_collision"
    reserve_dir = run_dir / "rst"
    reserve_dir.mkdir(parents=True)
    ambient_reservation = reserve_dir / ".io_node.00000.g42.reserve"
    ambient_reservation.write_text("ambient")
    env = os.environ.copy()
    env["ATHENAK_TEST_NODE_RESTART_GENERATION"] = "42"

    subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            *RESTART_ONLY_NODE_OVERRIDES,
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
        env=env,
    )

    manifest = reserve_dir / "io_node.00000.rst"
    assert ".g43.payload.rst" in manifest.read_text()
    assert ambient_reservation.read_text() == "ambient"
    assert sorted(reserve_dir.glob(".*.reserve")) == [ambient_reservation]


@pytest.mark.parametrize(
    ("stage", "expected"),
    (
        ("after_payload_write", "Injected node restart failure after payload write"),
        (
            "after_payload_publication",
            "Injected node restart failure after payload publication",
        ),
    ),
)
def test_injected_node_restart_failure_discards_owned_attempt_artifacts(
    tmp_path, stage, expected
):
    run_dir = tmp_path / f"injected_{stage}_failure"
    run_dir.mkdir()
    env = os.environ.copy()
    env["ATHENAK_TEST_NODE_RESTART_FAIL_STAGE"] = stage

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
            *RESTART_ONLY_NODE_OVERRIDES,
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=90,
    )

    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)
    assert not list((run_dir / "rst").rglob("*.tmp"))
    assert not list((run_dir / "rst").rglob("*.payload.rst"))
    assert not list((run_dir / "rst").glob(".*.reserve"))
    assert not (run_dir / "rst" / "io_node.00000.rst").exists()


@pytest.mark.parametrize(
    ("stage", "expected"),
    (
        (
            "after_manifest_publication",
            "Injected node restart failure after manifest publication",
        ),
        (
            "reservation_removal",
            "Injected node restart generation reservation removal failure",
        ),
    ),
)
def test_post_commit_node_restart_failure_preserves_published_checkpoint(
    tmp_path, stage, expected
):
    run_dir = tmp_path / f"injected_{stage}_failure"
    run_dir.mkdir()
    env = os.environ.copy()
    env["ATHENAK_TEST_NODE_RESTART_FAIL_STAGE"] = stage
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
            *RESTART_ONLY_NODE_OVERRIDES,
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=90,
    )

    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)
    manifest = run_dir / "rst" / "io_node.00000.rst"
    assert manifest.exists()
    assert list((run_dir / "rst").rglob("*.payload.rst"))
    assert not list((run_dir / "rst").rglob("*.tmp"))
    assert not list((run_dir / "rst").glob(".*.reserve"))
    _resume(tmp_path / f"{stage}_resume", manifest)


def test_restart_persisted_exhausted_counter_can_resume_without_publication(tmp_path):
    input_file = tmp_path / "counter_sentinel.athinput"
    input_file.write_text(
        Path(INPUT_FILE).read_text().replace(
            "<output6>\n", "<output6>\nfile_number = 2147483646\n", 1
        )
    )
    run_dir = tmp_path / "counter_sentinel"
    run_dir.mkdir()
    overrides = tuple(f"output{number}/dt=-1" for number in range(1, 6)) + (
        "time/final_output_policy=none",
    )
    subprocess.run(
        ["mpirun", "-np", "2", "./athena", "-i", str(input_file), "-d", str(run_dir),
         *overrides],
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
    )
    restart = run_dir / "rst" / "io_node.2147483646.rst"
    assert restart.exists()
    _resume(tmp_path / "counter_sentinel_resume", restart)


def test_node_restart_rejects_generated_payload_path_before_open(tmp_path):
    run_dir = tmp_path / "oversized_generated_payload_path"
    run_dir.mkdir()
    basename = "x" * 1000
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
            f"job/basename={basename}",
            "output1/dt=-1",
            "output2/dt=-1",
            "output3/dt=-1",
            "output4/dt=-1",
            "output5/dt=-1",
            "output6/single_file_per_node=true",
            "time/final_output_policy=none",
        ],
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert proc.returncode != 0
    assert "generated payload path exceeds the 1024-byte limit" in (
        proc.stdout + proc.stderr
    )
    assert not list((run_dir / "rst").rglob("*payload.rst*"))


@pytest.mark.parametrize(
    ("corruption", "expected"),
    (
        ("traversal", "payload path"),
        ("absolute", "payload path"),
        ("incomplete", "completion record"),
        ("byte_count", "payload block count or byte count"),
        ("payload_zero", "payload block count must be positive"),
        ("payload_count", "'payload_count' value must be between"),
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
        ("segment_zero", "segment block count must be positive"),
        ("segment_inventory", "too many segment records"),
        ("oversized_signature", "manifest signature exceeds"),
        ("oversized_scalar", "'complete' record exceeds"),
        ("oversized_payload", "payload inventory record exceeds"),
        ("oversized_segment", "segment inventory record exceeds"),
        ("oversized_trailing", "trailing record exceeds"),
        ("oversized_payload_path", "generated path exceeds"),
        ("sparse_total_manifest", "node restart manifest exceeds"),
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
    elif corruption == "payload_zero":
        lines = text.splitlines()
        index = next(index for index, line in enumerate(lines)
                     if line.startswith("payload "))
        fields = lines[index].split()
        fields[2] = "0"
        lines[index] = " ".join(fields)
        text = "\n".join(lines) + "\n"
    elif corruption == "payload_count":
        text = text.replace("payload_count=1", "payload_count=1000000000", 1)
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
    elif corruption == "oversized_signature":
        lines = text.splitlines()
        lines[0] = "AthenaK node restart manifest version=" + OVERSIZED_MANIFEST_LINE
        text = "\n".join(lines) + "\n"
    elif corruption == "oversized_scalar":
        text = text.replace("complete=1", "complete=" + OVERSIZED_MANIFEST_LINE, 1)
    elif corruption == "oversized_payload":
        lines = text.splitlines()
        index = next(index for index, line in enumerate(lines)
                     if line.startswith("payload "))
        lines[index] += " " + OVERSIZED_MANIFEST_LINE
        text = "\n".join(lines) + "\n"
    elif corruption == "oversized_segment":
        lines = text.splitlines()
        index = next(index for index, line in enumerate(lines)
                     if line.startswith("segment "))
        lines[index] += " " + OVERSIZED_MANIFEST_LINE
        text = "\n".join(lines) + "\n"
    elif corruption == "oversized_trailing":
        text += OVERSIZED_MANIFEST_LINE + "\n"
    elif corruption == "oversized_payload_path":
        lines = text.splitlines()
        index = next(index for index, line in enumerate(lines)
                     if line.startswith("payload "))
        fields = lines[index].split()
        directory, leaf = fields[4].split("/", 1)
        prefix = leaf.split(".g", 1)[0]
        fields[4] = (
            f"{directory}/{prefix}.g"
            + "1" * (MAX_GENERATED_PAYLOAD_PATH_BYTES + 1)
            + ".payload.rst"
        )
        lines[index] = " ".join(fields)
        text = "\n".join(lines) + "\n"
    elif corruption == "sparse_total_manifest":
        manifest.write_text(text)
        with manifest.open("r+b") as stream:
            stream.seek(MAX_NODE_MANIFEST_BYTES)
            stream.write(b"\n")
        text = None
    elif corruption == "mixed_generation":
        _split_manifest_into_two_payloads(manifest)
        lines = text.splitlines()
        lines = manifest.read_text().splitlines()
        payload_indexes = [
            index for index, line in enumerate(lines) if line.startswith("payload ")
        ]
        fields = lines[payload_indexes[1]].split()
        fields[4] = fields[4].replace(".payload.rst", "1.payload.rst", 1)
        lines[payload_indexes[1]] = " ".join(fields)
        text = "\n".join(lines) + "\n"
    elif corruption == "header_mismatch":
        _split_manifest_into_two_payloads(manifest)
        lines = manifest.read_text().splitlines()
        payload_indexes = [
            index for index, line in enumerate(lines) if line.startswith("payload ")
        ]
        fields = lines[payload_indexes[1]].split()
        header_size = next(int(line.split("=")[1]) for line in lines
                           if line.startswith("header_size="))
        duplicate = manifest.parent / fields[4]
        duplicate_bytes = duplicate.read_bytes()
        header = bytearray(duplicate_bytes[:header_size])
        header[0] ^= 1
        duplicate.write_bytes(header + duplicate_bytes[header_size:])
        text = "\n".join(lines) + "\n"
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
        elif corruption == "segment_zero":
            fields[3] = "0"
        elif corruption == "segment_inventory":
            nmb_total = next(
                int(line.split("=")[1])
                for line in lines
                if line.startswith("nmb_total=")
            )
            lines[segment_indexes[0]:segment_indexes[0]] = [
                lines[segment_indexes[0]]
            ] * nmb_total
            text = "\n".join(lines) + "\n"
            manifest.write_text(text)
            fields = None
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
    if text is not None:
        manifest.write_text(text)
    proc = subprocess.run(
        ["mpirun", "-np", "2", "./athena", "-r", str(manifest), "-d", str(run_dir)],
        capture_output=True,
        text=True,
        timeout=90,
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
    ("exact", "dot", "repeated_slash", "symlink", "temporary", "hardlink", "copied"),
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
    elif alias == "temporary":
        restart = Path(str(payload) + ".tmp")
        shutil.copyfile(payload, restart)
    elif alias == "hardlink":
        restart = tmp_path / "payload_hardlink_alias.rst"
        os.link(payload, restart)
    else:
        restart = tmp_path / "payload_copy_alias.rst"
        shutil.copyfile(payload, restart)
    proc = _resume(
        tmp_path / f"payload_entry_resume_{alias}", restart, check=False
    )
    assert proc.returncode != 0
    assert "use the public manifest path" in (proc.stdout + proc.stderr)


def test_node_restart_rejects_corrupt_payload_marker(tmp_path, node_restart_template):
    _, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "payload_marker")
    payload = _payload_path(manifest)
    data = payload.read_bytes()
    assert NODE_PAYLOAD_MARKER in data
    payload.write_bytes(
        data.replace(NODE_PAYLOAD_MARKER, b"X" * len(NODE_PAYLOAD_MARKER), 1)
    )
    proc = _resume(tmp_path / "payload_marker_resume", manifest, check=False)
    assert proc.returncode != 0
    assert "payload marker is absent or invalid" in (proc.stdout + proc.stderr)


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


def test_node_restart_one_rank_combines_interleaved_spans_from_multiple_payloads(
    tmp_path, node_restart_template
):
    _, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "interleaved")
    _interleave_manifest_across_two_payloads(manifest)
    _resume(tmp_path / "interleaved_resume", manifest, env=_env("7"), nranks=1)


def test_node_restart_supports_rank_count_change_on_local_node(
    tmp_path, node_restart_template
):
    _, manifest = _copy_node_checkpoint(tmp_path, node_restart_template, "rank_count")
    # This workstation exercises redistribution across ranks on one node. True
    # multi-node payload routing and non-owning-node qualification remain external gates.
    _resume(tmp_path / "rank_count_resume", manifest, nranks=1)


def test_existing_per_rank_restart_resumes(tmp_path):
    run_dir, _ = _run(
        tmp_path,
        "per_rank",
        "output1/single_file_per_rank=true",
        "output6/single_file_per_rank=true",
        "time/final_output_policy=restart_only",
    )
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
        timeout=90,
    )


def test_origin_main_per_rank_restart_fixture_resumes(tmp_path):
    run_dir = tmp_path / "origin_main_per_rank_resume"
    run_dir.mkdir()
    rank0 = (
        FIXTURES
        / "rst"
        / "per_rank"
        / "rank_00000000"
        / "io_legacy_per_rank.00001.rst"
    )
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
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
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
        timeout=90,
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
        timeout=90,
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
        / "io_node_example.density.r_2.5000000000000000e-01.00000.sph.bin"
    )
    assert bin_convert.read_binary(str(node_bin), assemble_shards=True)["n_mbs"] == 4
    assert (
        bin_convert.read_coarsened_binary(str(node_cbin), assemble_shards=True)["n_mbs"]
        == 4
    )
    assert read_pdf(str(node_pdf))["header"]["distribution"] == "node"
    assert read_sphslice(str(node_surface))["data"].shape == (16, 32, 1)
    assert (run_dir / "rst" / "io_node_example.00001.rst").exists()
    assert "PERFORMANCE_REGION outputs " in proc.stdout

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
        timeout=30,
    )
    assert "binary meshblocks=4" in summary.stdout

    coarsened_summary = subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"),
            "cbin",
            str(node_cbin),
            "--assemble-shards",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert "coarsened meshblocks=4" in coarsened_summary.stdout

    pdf_summary = subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"),
            "pdf",
            str(node_pdf),
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert "pdf shape=" in pdf_summary.stdout

    sphslice_summary = subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"),
            "sphslice",
            str(node_surface),
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert "sphslice shape=(16, 32, 1)" in sphslice_summary.stdout

    manifest = run_dir / "rst" / "io_node_example.00001.rst"
    _resume(tmp_path / "node_example_resume", manifest)
    assert not list(tmp_path.rglob("*.assembled"))
