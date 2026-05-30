"""Serial integration coverage for widened writer sequences and temp cleanup."""

from pathlib import Path
import subprocess

import pytest


INPUT_FILE = "inputs/io_node_sharding.athinput"
SEQUENCES = ("99999", "100000", "100001")
EXPECTED_OUTPUTS = (
    ("bin", "bin/io_node.full.{sequence}.bin"),
    ("cbin", "cbin_coarse_2/io_node.coarse.{sequence}.cbin"),
    ("pdf", "pdf_node_pdf_hydro_w_d/io_node.{sequence}.pdf"),
    ("sphslice", "bin/io_node.density.r_*.{sequence}.sph.bin"),
    ("rst", "rst/io_node.{sequence}.rst"),
)
WRITER_TARGETS = {
    "bin": (1, Path("bin/io_node.full.00000.bin")),
    "cbin": (3, Path("cbin_coarse_2/io_node.coarse.00000.cbin")),
    "pdf": (4, Path("pdf_node_pdf_hydro_w_d/io_node.00000.pdf")),
    "sphslice": (
        5,
        Path("bin/io_node.density.r_2.5000000000000000e-01.00000.sph.bin"),
    ),
    "rst": (6, Path("rst/io_node.00000.rst")),
}


def _target_only_overrides(target_number):
    return [
        *(f"output{number}/dt=-1" for number in range(1, 7) if number != target_number),
        "time/nlim=0",
        "time/tlim=0.0",
        "time/final_output_policy=none",
    ]


def _run_target(run_dir, target_number):
    return subprocess.run(
        [
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            *_target_only_overrides(target_number),
        ],
        capture_output=True,
        text=True,
        timeout=90,
    )


@pytest.fixture(scope="module", params=SEQUENCES)
def rcp02_serial_outputs(request, tmp_path_factory):
    sequence = request.param
    run_dir = tmp_path_factory.mktemp(f"rcp02_serial_outputs_{sequence}")
    input_file = run_dir / "io_node_sharding_widened.athinput"
    input_text = Path(INPUT_FILE).read_text()
    for number in range(1, 7):
        marker = f"<output{number}>\n"
        assert marker in input_text
        input_text = input_text.replace(
            marker, marker + f"file_number = {sequence}\n", 1
        )
    input_file.write_text(input_text)
    overrides = [
        "time/nlim=0",
        "time/tlim=0.0",
        "time/final_output_policy=none",
    ]
    subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir), *overrides],
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
    )
    return sequence, run_dir


@pytest.mark.parametrize(("kind", "pattern"), EXPECTED_OUTPUTS)
def test_serial_feature_writer_sequence_boundaries(
    rcp02_serial_outputs, kind, pattern
):
    sequence, run_dir = rcp02_serial_outputs
    pattern = pattern.format(sequence=sequence)
    matches = list(run_dir.glob(pattern))
    assert len(matches) == 1, f"{kind} output did not match {pattern}: {matches}"


def test_serial_feature_writer_publication_leaves_no_temporary_files(
    rcp02_serial_outputs,
):
    _, run_dir = rcp02_serial_outputs
    leftovers = [
        path
        for path in run_dir.rglob("*")
        if ".tmp" in path.name
    ]
    assert not leftovers


@pytest.mark.parametrize(("kind", "target"), WRITER_TARGETS.items())
def test_serial_feature_writer_replaces_stale_temporary_file(tmp_path, kind, target):
    target_number, relative_path = target
    run_dir = tmp_path / kind
    temporary = run_dir / f"{relative_path}.tmp"
    temporary.parent.mkdir(parents=True)
    temporary.write_text("stale temporary")

    proc = _run_target(run_dir, target_number)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (run_dir / relative_path).is_file()
    assert not temporary.exists()


@pytest.mark.parametrize(("kind", "target"), WRITER_TARGETS.items())
def test_serial_feature_writer_failed_rename_discards_temporary_file(
    tmp_path, kind, target
):
    target_number, relative_path = target
    run_dir = tmp_path / kind
    published = run_dir / relative_path
    published.mkdir(parents=True)
    temporary = run_dir / f"{relative_path}.tmp"

    proc = _run_target(run_dir, target_number)

    assert proc.returncode != 0
    assert "could not publish output file" in proc.stdout + proc.stderr
    assert not temporary.exists()


def test_serial_binary_unwritable_directory_rejects_without_temporary_file(tmp_path):
    run_dir = tmp_path / "unwritable"
    binary_dir = run_dir / "bin"
    binary_dir.mkdir(parents=True)
    binary_dir.chmod(0o555)
    probe = binary_dir / "probe"
    try:
        try:
            probe.write_text("permission probe")
        except PermissionError:
            pass
        else:
            probe.unlink()
            pytest.skip("platform permits writes to a mode-0555 directory")

        proc = _run_target(run_dir, 1)
    finally:
        binary_dir.chmod(0o755)

    assert proc.returncode != 0
    assert "could not be opened" in proc.stdout + proc.stderr
    assert not (binary_dir / "io_node.full.00000.bin.tmp").exists()
