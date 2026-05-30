"""Regression tests for final-output policy and opt-in IO timing."""

from pathlib import Path
import re
import struct
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
INPUT_FILE = "inputs/io_finalization_timing.athinput"
METADATA_INPUT_FILE = "inputs/io_restart_metadata.athinput"
INT_BYTES = struct.calcsize("=i")
REAL_BYTES = struct.calcsize("=d")
REGION_SIZE_BYTES = struct.calcsize("=9d")
REGION_INDICES_BYTES = struct.calcsize("=19i")


def _subprocess_run(*args, **kwargs):
    kwargs.setdefault("timeout", 90)
    return subprocess.run(*args, **kwargs)


def _run_case(tmp_path: Path, *overrides: str):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = _subprocess_run(
        ["./athena", "-i", INPUT_FILE, "-d", str(run_dir), *overrides],
        check=True,
        capture_output=True,
        text=True,
    )
    bin_files = sorted((run_dir / "bin").glob("*.bin"))
    rst_files = sorted((run_dir / "rst").glob("*.rst"))
    return proc.stdout, bin_files, rst_files


def test_default_policy_retains_terminal_outputs(tmp_path):
    stdout, bin_files, rst_files = _run_case(tmp_path)
    assert "[output-io]" not in stdout
    assert len(bin_files) == 2
    assert len(rst_files) == 2


def test_restart_only_suppresses_terminal_diagnostic(tmp_path):
    _, bin_files, rst_files = _run_case(
        tmp_path, "time/final_output_policy=restart_only"
    )
    assert len(bin_files) == 1
    assert len(rst_files) == 2


def test_none_suppresses_all_terminal_outputs(tmp_path):
    _, bin_files, rst_files = _run_case(tmp_path, "time/final_output_policy=none")
    assert len(bin_files) == 1
    assert len(rst_files) == 1


def test_timing_is_opt_in(tmp_path):
    stdout, _, _ = _run_case(tmp_path, "time/output_timing=true")
    assert stdout.count("[output-io] event=initial ") == 2
    assert stdout.count("[output-io] event=final ") == 2
    assert (
        "event=initial block=output1 type=bin distribution=shared "
        "elapsed_max_s=" in stdout
    )
    assert (
        "event=initial block=output2 type=rst distribution=shared "
        "elapsed_max_s=" in stdout
    )


def test_invalid_final_output_policy_is_rejected(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = _subprocess_run(
        [
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "time/final_output_policy=invalid",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "final_output_policy = 'invalid' not implemented" in proc.stdout


def test_terminal_restart_resume_advances_counter_without_overwrite(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    _subprocess_run(
        [
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    terminal_restart = run_dir / "rst" / "io_policy.00001.rst"
    saved_terminal_bytes = terminal_restart.read_bytes()
    _subprocess_run(
        [
            "./athena",
            "-r",
            str(terminal_restart),
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
    assert terminal_restart.read_bytes() == saved_terminal_bytes
    assert (run_dir / "rst" / "io_policy.00002.rst").exists()


def test_origin_main_shared_restart_fixture_resumes(tmp_path):
    run_dir = tmp_path / "origin_main_shared_resume"
    run_dir.mkdir()
    restart = FIXTURES / "rst" / "shared" / "io_legacy_shared.00001.rst"
    _subprocess_run(
        [
            "./athena",
            "-r",
            str(restart),
            "-d",
            str(run_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def test_missing_restart_path_fails_immediately(tmp_path):
    run_dir = tmp_path / "missing_restart_run"
    run_dir.mkdir()
    restart = tmp_path / "missing.rst"
    proc = _subprocess_run(
        ["./athena", "-r", str(restart), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    output = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert f"Unable to open restart file: {restart}" in output
    assert "Error opening file" not in output
    assert "could not be opened" not in output


def test_corrupt_restart_meshblock_extent_fails_before_mesh_arithmetic(tmp_path):
    restart = tmp_path / "corrupt_meshblock_extent.rst"
    shutil.copyfile(
        FIXTURES / "rst" / "shared" / "io_legacy_shared.00001.rst", restart
    )
    data = bytearray(restart.read_bytes())
    header = data.index(b"<par_end>") + len(b"<par_end>\n")
    meshblock_nx1 = header + 2 * INT_BYTES + REGION_SIZE_BYTES
    meshblock_nx1 += REGION_INDICES_BYTES + INT_BYTES
    struct.pack_into("=i", data, meshblock_nx1, 0)
    restart.write_bytes(data)
    run_dir = tmp_path / "corrupt_meshblock_extent_run"
    run_dir.mkdir()
    proc = _subprocess_run(
        ["./athena", "-r", str(restart), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    output = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "restart MeshBlock indices are inconsistent." in output
    assert "Root grid =" not in output


def _corrupt_origin_main_restart(tmp_path, name, mutate):
    restart = tmp_path / f"{name}.rst"
    shutil.copyfile(
        FIXTURES / "rst" / "shared" / "io_legacy_shared.00001.rst", restart
    )
    data = bytearray(restart.read_bytes())
    header = data.index(b"<par_end>") + len(b"<par_end>\n")
    mutate(data, header)
    restart.write_bytes(data)
    return restart


def _resume_corrupt_restart(tmp_path, restart):
    run_dir = tmp_path / f"{restart.stem}_run"
    run_dir.mkdir()
    proc = _subprocess_run(
        ["./athena", "-r", str(restart), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    return proc, proc.stdout + proc.stderr


def _restart_metadata_offsets(data):
    header = data.index(b"<par_end>") + len(b"<par_end>\n")
    nmb_total = struct.unpack_from("=i", data, header)[0]
    locations = (
        header + 3 * INT_BYTES + 2 * REAL_BYTES + REGION_SIZE_BYTES
        + 2 * REGION_INDICES_BYTES
    )
    costs = locations + nmb_total * struct.calcsize("=4i")
    return header, nmb_total, locations, costs


def _rewrite_restart_inventory(restart, locations, costs):
    data = bytearray(restart.read_bytes())
    header, old_nmb_total, old_locations, old_costs = _restart_metadata_offsets(data)
    old_end = old_costs + old_nmb_total * struct.calcsize("=f")
    encoded_locations = b"".join(struct.pack("=4i", *location) for location in locations)
    encoded_costs = b"".join(struct.pack("=f", cost) for cost in costs)
    data[old_locations:old_end] = encoded_locations + encoded_costs
    struct.pack_into("=i", data, header, len(locations))
    restart.write_bytes(data)


@pytest.fixture(scope="module")
def restart_metadata_template(tmp_path_factory):
    run_dir = tmp_path_factory.mktemp("restart_metadata_template")
    _subprocess_run(
        ["./athena", "-i", METADATA_INPUT_FILE, "-d", str(run_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    return run_dir / "rst" / "io_restart_metadata.00000.rst"


def test_corrupt_restart_mesh_spacing_fails_before_tree_reconstruction(tmp_path):
    def corrupt_dx1(data, header):
        dx1 = header + 2 * INT_BYTES + 6 * REAL_BYTES
        struct.pack_into("=d", data, dx1, 2.0 * struct.unpack_from("=d", data, dx1)[0])

    restart = _corrupt_origin_main_restart(tmp_path, "corrupt_mesh_spacing", corrupt_dx1)
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "restart mesh spacing is inconsistent with its bounds." in output
    assert "Root grid =" not in output


def test_corrupt_restart_logical_level_fails_before_tree_reconstruction(tmp_path):
    def corrupt_level(data, header):
        logical_locations = (
            header + 3 * INT_BYTES + 2 * REAL_BYTES + REGION_SIZE_BYTES
            + 2 * REGION_INDICES_BYTES
        )
        struct.pack_into("=i", data, logical_locations + 3 * INT_BYTES, -1)

    restart = _corrupt_origin_main_restart(
        tmp_path, "corrupt_logical_level", corrupt_level
    )
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "restart MeshBlock logical level is inconsistent." in output
    assert "Root grid =" not in output


def test_corrupt_restart_meshblock_cost_fails_before_tree_reconstruction(tmp_path):
    def corrupt_cost(data, header):
        nmb_total = struct.unpack_from("=i", data, header)[0]
        logical_locations = (
            header + 3 * INT_BYTES + 2 * REAL_BYTES + REGION_SIZE_BYTES
            + 2 * REGION_INDICES_BYTES
        )
        costs = logical_locations + nmb_total * struct.calcsize("=4i")
        struct.pack_into("=f", data, costs, float("nan"))

    restart = _corrupt_origin_main_restart(
        tmp_path, "corrupt_meshblock_cost", corrupt_cost
    )
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "restart MeshBlock cost is inconsistent." in output
    assert "Root grid =" not in output


def test_corrupt_restart_dimension_flags_fail_before_tree_reconstruction(tmp_path):
    def corrupt_nx2(data, header):
        mesh_indices = header + 2 * INT_BYTES + REGION_SIZE_BYTES
        struct.pack_into("=i", data, mesh_indices + 2 * INT_BYTES, 1)

    restart = _corrupt_origin_main_restart(tmp_path, "corrupt_dimension", corrupt_nx2)
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "restart mesh dimensions disagree with input parameters." in output
    assert "Root grid =" not in output


def test_corrupt_restart_inactive_coarse_axis_fails_before_tree(tmp_path):
    run_dir = tmp_path / "inactive_coarse_axis_source"
    run_dir.mkdir()
    _subprocess_run(
        ["./athena", "-i", INPUT_FILE, "-d", str(run_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    restart = run_dir / "rst" / "io_policy.00000.rst"
    data = bytearray(restart.read_bytes())
    header = data.index(b"<par_end>") + len(b"<par_end>\n")
    meshblock_indices = header + 2 * INT_BYTES + REGION_SIZE_BYTES
    meshblock_indices += REGION_INDICES_BYTES
    struct.pack_into("=i", data, meshblock_indices + 15 * INT_BYTES, 1)
    restart.write_bytes(data)
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "restart MeshBlock indices are inconsistent." in output
    assert "Root grid =" not in output


def test_corrupt_multilevel_restart_odd_meshblock_extent_fails_before_tree(tmp_path):
    def corrupt_odd_extent(data, header):
        dx1 = header + 2 * INT_BYTES + 6 * REAL_BYTES
        struct.pack_into("=d", data, dx1, 1.0 / 14.0)
        mesh_indices = header + 2 * INT_BYTES + REGION_SIZE_BYTES
        struct.pack_into("=i", data, mesh_indices + INT_BYTES, 14)
        struct.pack_into("=i", data, mesh_indices + 5 * INT_BYTES, 15)
        meshblock_indices = mesh_indices + REGION_INDICES_BYTES
        struct.pack_into("=i", data, meshblock_indices + INT_BYTES, 7)
        struct.pack_into("=i", data, meshblock_indices + 5 * INT_BYTES, 8)
        struct.pack_into("=i", data, meshblock_indices + 10 * INT_BYTES, 3)
        struct.pack_into("=i", data, meshblock_indices + 14 * INT_BYTES, 4)

    restart = _corrupt_origin_main_restart(
        tmp_path, "corrupt_odd_extent", corrupt_odd_extent
    )
    data = bytearray(restart.read_bytes())
    data = data.replace(b"refinement = none", b"refinement = static", 1)
    restart.write_bytes(data)
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "restart MeshBlock indices are inconsistent." in output
    assert "Root grid =" not in output


def test_corrupt_restart_adaptive_level_count_fails_before_tree(
    tmp_path, restart_metadata_template
):
    restart = tmp_path / "corrupt_adaptive_levels.rst"
    shutil.copyfile(restart_metadata_template, restart)
    data = restart.read_bytes()
    data, refinement_count = re.subn(
        rb"(refinement\s*=\s*)static", rb"\1adaptive", data, count=1
    )
    data, level_count = re.subn(
        rb"(num_levels\s*=\s*)1", rb"\g<1>2147483647", data, count=1
    )
    assert refinement_count == 1
    assert level_count == 1
    restart.write_bytes(data)
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert "Number of refinement levels must be between 1 and 30" in output


def test_generated_3d_restart_with_nonzero_x2_resumes(tmp_path):
    run_dir = tmp_path / "generated_3d"
    run_dir.mkdir()
    _subprocess_run(
        [
            "./athena",
            "-i",
            "inputs/io_node_sharding.athinput",
            "-d",
            str(run_dir),
            "mesh/nx2=16",
            "time/final_output_policy=none",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    restart = run_dir / "rst" / "io_node.00000.rst"
    proc, _ = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode == 0


@pytest.mark.parametrize(
    ("name", "locations", "costs", "expected"),
    (
        (
            "duplicate",
            ((0, 0, 0, 1), (0, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock inventory contains a duplicate location.",
        ),
        (
            "permuted",
            ((1, 0, 0, 1), (0, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock inventory is not in canonical order.",
        ),
        (
            "overlap",
            ((0, 0, 0, 1), (0, 0, 0, 2), (1, 0, 0, 2), (1, 0, 0, 1)),
            (1.0, 1.0, 1.0, 1.0),
            "restart MeshBlock inventory contains overlapping levels.",
        ),
        (
            "incomplete",
            ((0, 0, 0, 2), (1, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock inventory has incomplete refinement.",
        ),
        (
            "imbalanced",
            ((0, 0, 0, 2), (2, 0, 0, 3), (3, 0, 0, 3), (1, 0, 0, 1)),
            (1.0, 1.0, 1.0, 1.0),
            "Neighbor search failed",
        ),
        (
            "level31",
            ((0, 0, 0, 31), (1, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock logical level is inconsistent.",
        ),
        (
            "below_root",
            ((0, 0, 0, 0), (1, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock logical level is inconsistent.",
        ),
        (
            "active_axis_range",
            ((2, 0, 0, 1), (1, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock logical location is inconsistent.",
        ),
        (
            "inactive_axis_range",
            ((0, 1, 0, 1), (1, 0, 0, 1)),
            (1.0, 1.0),
            "restart MeshBlock logical location is inconsistent.",
        ),
        (
            "zero_cost",
            ((0, 0, 0, 1), (1, 0, 0, 1)),
            (0.0, 1.0),
            "restart MeshBlock cost is inconsistent.",
        ),
        (
            "infinite_cost",
            ((0, 0, 0, 1), (1, 0, 0, 1)),
            (float("inf"), 1.0),
            "restart MeshBlock cost is inconsistent.",
        ),
        (
            "aggregate_cost",
            ((0, 0, 0, 1), (1, 0, 0, 1)),
            (3.4e38, 3.4e38),
            "restart MeshBlock aggregate cost is inconsistent.",
        ),
    ),
)
def test_corrupt_restart_inventory_fails_before_payload_loading(
    tmp_path, restart_metadata_template, name, locations, costs, expected
):
    restart = tmp_path / f"{name}.rst"
    shutil.copyfile(restart_metadata_template, restart)
    _rewrite_restart_inventory(restart, locations, costs)
    proc, output = _resume_corrupt_restart(tmp_path, restart)
    assert proc.returncode != 0
    assert expected in output
