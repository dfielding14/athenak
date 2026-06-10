"""Focused tests for the lean direct Stage I campaign launcher."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import re
import stat
import struct
from types import SimpleNamespace
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
LAUNCHER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast.py"
EXPECTED_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
EXPECTED_SOURCE = Path("/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422")
EXPECTED_MATRIX_RELATIVE = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
EXPECTED_MATRIX_SHA256 = (
    "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
)
EXPECTED_ATHENA = (
    EXPECTED_ROOT
    / "build/frontier-hip-9e07542281e4-cpe25.09-cce20-rocm6.4.2/src/athena"
)
EXPECTED_ATHENA_SHA256 = (
    "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
)


def load_launcher():
    """Import the standalone launcher by path."""

    if not LAUNCHER.is_file():
        pytest.skip(f"forthcoming lean launcher is not present: {LAUNCHER}")
    name = "cgl_lf_stage_i_fast_for_tests"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fast():
    return load_launcher()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_history(path: Path, columns: dict[str, list[float]]) -> None:
    names = list(columns)
    lines = [
        "# " + " ".join(
            f"[{index}]={name}" for index, name in enumerate(names, start=1)
        )
    ]
    lines.extend(
        " ".join(format(value, ".17g") for value in row)
        for row in zip(*(columns[name] for name in names))
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def restart_payload(
    time: float,
    *,
    mesh_nx1: int = 1,
    local_blocks: int = 1,
) -> bytes:
    parameter_header = (
        f"<mesh>\nnx1 = {mesh_nx1}\nnx2 = 1\nnx3 = 1\n"
        "<meshblock>\nnx1 = 1\nnx2 = 1\nnx3 = 1\n"
        f"<time>\nrestart_time = {time:.17g}\n"
        "<par_end>\n"
    ).encode()
    mesh_header = bytearray(252)
    struct.pack_into("<ii", mesh_header, 0, mesh_nx1, 0)
    struct.pack_into("<d", mesh_header, 232, time)
    struct.pack_into("<d", mesh_header, 240, 0.01)
    struct.pack_into("<i", mesh_header, 248, 1)
    locations = b"".join(
        struct.pack("<4i", location, 0, 0, 0)
        for location in range(mesh_nx1)
    )
    costs = struct.pack(f"<{mesh_nx1}f", *([1.0] * mesh_nx1))
    metadata = bytearray(248)
    struct.pack_into("<i", metadata, 4, 1)
    data_size = 8
    return (
        parameter_header
        + mesh_header
        + locations
        + costs
        + metadata
        + bytes(296)
        + bytes(6 * 8)
        + struct.pack("<d", 0.0)
        + struct.pack("<18d", *([0.0] * 18))
        + struct.pack("<Q", data_size)
        + bytes(local_blocks * data_size)
    )


def binary_payload(
    time: float,
    *,
    mesh_nx1: int = 1,
    locations: tuple[int, ...] = (0,),
) -> bytes:
    parameter_header = (
        f"<mesh>\nnx1 = {mesh_nx1}\nnx2 = 1\nnx3 = 1\n"
        "<meshblock>\nnx1 = 1\nnx2 = 1\nnx3 = 1\n"
        "<par_end>\n"
    ).encode()
    header = (
        "Athena binary output version=1.1\n"
        "  size of preheader=5\n"
        f"  time={time:.17g}\n"
        "  cycle=1\n"
        "  size of location=8\n"
        "  size of variable=4\n"
        "  number of variables=1\n"
        "  variables:  dens\n"
        f"  header offset={len(parameter_header)}\n"
    ).encode() + parameter_header
    blocks = b"".join(
        struct.pack("<10i", 0, 0, 0, 0, 0, 0, location, 0, 0, 0)
        + struct.pack(
            "<6d", float(location), float(location + 1), 0.0, 1.0, 0.0, 1.0
        )
        + struct.pack("<f", 1.0)
        for location in locations
    )
    return header + blocks


def write_binary_rank_group(
    directory: Path,
    filename: str,
    rank_count: int,
    time: float,
) -> Path:
    rank_zero = directory / "rank_00000000" / filename
    for rank in range(rank_count):
        path = directory / f"rank_{rank:08d}" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            binary_payload(time, mesh_nx1=rank_count, locations=(rank,))
        )
    return rank_zero


def write_restart_rank_group(
    directory: Path,
    filename: str,
    rank_count: int,
    time: float,
) -> Path:
    rank_zero = directory / "rank_00000000" / filename
    for rank in range(rank_count):
        path = directory / f"rank_{rank:08d}" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            restart_payload(time, mesh_nx1=rank_count, local_blocks=1)
        )
    return rank_zero


def write_rank_group(
    directory: Path,
    filename: str,
    rank_count: int,
    payload: bytes,
) -> Path:
    rank_zero = directory / "rank_00000000" / filename
    for rank in range(rank_count):
        path = directory / f"rank_{rank:08d}" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    return rank_zero


def segment_fixture(
    root: Path,
    *,
    case_id: str = "R05",
    start_time: float = 0.0,
    target_time: float = 1.0,
    final_time: float = 1.0,
    strict_failure: float = 0.0,
) -> Path:
    segment = root / "fast_s000_t0_to_t1"
    output = segment / "output"
    mhd = {
        "time": [start_time, final_time],
        "mass": [1.0, 1.0],
        "lf_dfloor": [0.0, strict_failure],
        "lf_pfloor": [0.0, 0.0],
        "lf_nonfin": [0.0, 0.0],
        "lf_nonpos": [0.0, 0.0],
        "lf_hardbd": [0.0, 0.0],
    }
    user = {"time": [start_time, final_time], "mass": [1.0, 1.0]}
    write_history(output / "fixture.mhd.hst", mhd)
    write_history(output / "fixture.user.hst", user)
    write_restart_rank_group(output / "rst", "fixture.00001.rst", 2, final_time)
    write_binary_rank_group(output / "bin", "fixture.00001.bin", 2, final_time)
    write_json(
        segment / "manifest/fast_run.json",
        {
            "case_id": case_id,
            "ranks": 2,
            "start_time": start_time,
            "target_time": target_time,
        },
    )
    return segment


def submit_segment_fixture(root: Path, *, job_id: str | None = None) -> Path:
    segment = root / "fast_s000_t0_to_t10"
    write_json(
        segment / "manifest/fast_run.json",
        {"case_id": "R03", "job_id": job_id},
    )
    script = segment / "manifest/run.sbatch"
    script.write_text("#!/bin/bash\n", encoding="utf-8")
    return segment


def batch_manifest(root: Path) -> dict[str, object]:
    run_dir = root / "runs/R05/fast_s000_t0_to_t10"
    return {
        "case_id": "R05",
        "sequence": 0,
        "nodes": 4,
        "input": str(root / "input.athinput"),
        "input_sha256": "a" * 64,
        "restart": None,
        "restart_sha256": None,
        "output_dir": str(run_dir / "output"),
        "run_dir": str(run_dir),
        "slurm_log": str(root / "logs/%x.%j.log"),
        "run_basename": "E03_fast_fixture",
        "target_time": 10.0,
        "root": str(root),
        "fast_script": str(LAUNCHER),
        "fast_script_sha256": "b" * 64,
    }


def test_cli_exposes_lean_launch_status_analyze_advance(fast, monkeypatch, capsys):
    monkeypatch.setattr(fast.sys, "argv", [str(LAUNCHER), "--help"])
    with pytest.raises(SystemExit) as stopped:
        fast.main()
    assert stopped.value.code == 0
    match = re.search(r"\{([^}]+)\}", capsys.readouterr().out)
    assert match is not None
    assert set(match.group(1).split(",")) == {
        "launch",
        "status",
        "analyze",
        "advance",
    }


def test_frozen_defaults_and_resource_profiles_are_explicit(fast):
    assert fast.DEFAULT_ROOT == EXPECTED_ROOT
    assert fast.FROZEN_SOURCE == EXPECTED_SOURCE
    assert fast.MATRIX_RELATIVE == EXPECTED_MATRIX_RELATIVE
    assert fast.MATRIX_SHA256 == EXPECTED_MATRIX_SHA256
    assert fast.ATHENA == EXPECTED_ATHENA
    assert fast.ATHENA_SHA256 == EXPECTED_ATHENA_SHA256
    assert fast.TARGET_TIME == 10.0
    assert fast.CASE_NODES["R03"] == 1
    assert all(fast.CASE_NODES[f"R{number:02d}"] == 4 for number in range(4, 16))
    assert fast.CASE_NODES["R16"] == 1
    assert fast.CASE_NODES["R17"] == 8


def test_sha256_file_and_require_sha_bind_exact_bytes(fast, tmp_path):
    path = tmp_path / "frozen"
    path.write_bytes(b"frozen bytes\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()

    assert fast.sha256_file(path) == digest
    fast.require_sha(path, digest, "fixture")

    path.write_bytes(b"changed\n")
    with pytest.raises(fast.FastRunError, match="fixture checksum mismatch"):
        fast.require_sha(path, digest, "fixture")


def test_write_and_load_json_are_deterministic_and_replaceable(fast, tmp_path):
    path = tmp_path / "state/run.json"
    fast.write_json(path, {"z": 2, "a": 1})
    assert path.read_text(encoding="utf-8") == '{\n  "a": 1,\n  "z": 2\n}\n'

    fast.write_json(path, {"state": "submitted"})
    assert fast.load_json(path) == {"state": "submitted"}


def test_load_json_rejects_non_object(fast, tmp_path):
    path = tmp_path / "list.json"
    path.write_text("[]\n", encoding="utf-8")
    with pytest.raises(fast.FastRunError, match="expected JSON object"):
        fast.load_json(path)


def test_campaign_cases_authenticates_matrix_and_selects_supported_cases(
    fast, tmp_path, monkeypatch
):
    source = tmp_path / "source"
    matrix = source / EXPECTED_MATRIX_RELATIVE
    write_json(
        matrix,
        {
            "cases": [
                {"id": "R02", "name": "excluded"},
                {"id": "R03", "name": "three"},
                {"id": "R17", "name": "seventeen"},
            ]
        },
    )
    monkeypatch.setattr(fast, "FROZEN_SOURCE", source)
    monkeypatch.setattr(fast, "MATRIX_SHA256", fast.sha256_file(matrix))

    assert list(fast.campaign_cases()) == ["R03", "R17"]


def test_expand_cases_supports_ranges_commas_and_stable_deduplication(fast):
    assert fast.expand_cases(["R03,R04", "R04", "R05-R07"]) == [
        "R03",
        "R04",
        "R05",
        "R06",
        "R07",
    ]


@pytest.mark.parametrize("values", [["R02"], ["R18"], ["r03"], ["R17-R03"]])
def test_expand_cases_rejects_unknown_or_descending_requests(fast, values):
    with pytest.raises(fast.FastRunError):
        fast.expand_cases(values)


@pytest.mark.parametrize(
    ("value", "expected"),
    [(0.0, "0"), (0.5, "0p5"), (-0.25, "m0p25"), (10.0, "10")],
)
def test_time_tag_is_path_safe_and_stable(fast, value, expected):
    assert fast.time_tag(value) == expected


def test_segment_directories_filter_and_sort_by_numeric_sequence(fast, tmp_path):
    case = fast.case_root(tmp_path, "R03")
    for name in (
        "fast_s010_t1_to_t10",
        "fast_s002_t0p5_to_t10",
        "fast_s001_t0_to_t10",
        "other",
    ):
        (case / name).mkdir(parents=True)
    (case / "fast_s003_file").write_text("not a directory", encoding="utf-8")

    assert [path.name for path in fast.segment_directories(tmp_path, "R03")] == [
        "fast_s001_t0_to_t10",
        "fast_s002_t0p5_to_t10",
        "fast_s010_t1_to_t10",
    ]


def test_parse_history_returns_finite_named_columns(fast, tmp_path):
    path = tmp_path / "fixture.hst"
    write_history(path, {"time": [0.0, 1.0], "mass": [2.0, 2.0]})
    assert fast.parse_history(path) == {"time": [0.0, 1.0], "mass": [2.0, 2.0]}


@pytest.mark.parametrize("value", ["nan", "inf"])
def test_parse_history_rejects_nonfinite_values(fast, tmp_path, value):
    path = tmp_path / "fixture.hst"
    path.write_text(f"# [1]=time [2]=mass\n0 1\n1 {value}\n", encoding="utf-8")
    with pytest.raises(fast.FastRunError, match="non-finite"):
        fast.parse_history(path)


def test_restart_time_requires_one_finite_physical_time(fast, tmp_path):
    path = tmp_path / "fixture.rst"
    path.write_bytes(restart_payload(0.625))
    assert fast.restart_time(path) == 0.625

    path.write_bytes(b"<time>\ntime=0.5\nrestart_time=0.5\n<par_end>\n")
    with pytest.raises(fast.FastRunError, match="ambiguous"):
        fast.restart_time(path)

    path.write_bytes(b"<time>\nrestart_time=0.5\n<par_end>\nX\n")
    with pytest.raises(fast.FastRunError, match="mesh decomposition"):
        fast.restart_time(path)


def test_binary_time_requires_magic_and_finite_physical_time(fast, tmp_path):
    path = tmp_path / "fixture.bin"
    path.write_bytes(binary_payload(0.625))
    assert fast.binary_time(path) == 0.625

    path.write_bytes(b"not an Athena binary\n")
    with pytest.raises(fast.FastRunError, match="magic/version"):
        fast.binary_time(path)

    path.write_bytes(
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        b"  time=0.625\n"
    )
    with pytest.raises(fast.FastRunError, match="truncated cycle"):
        fast.binary_time(path)


def test_terminal_product_group_requires_a_complete_rank_group(fast, tmp_path):
    directory = tmp_path / "rst"
    terminal = write_restart_rank_group(directory, "fixture.00002.rst", 3, 2.0)

    result = fast.terminal_product_group(directory, ".rst", 3)

    assert result["rank_zero"] == str(terminal.resolve())
    assert result["rank_count"] == 3
    assert result["physical_time"] == 2.0

    (directory / "rank_00000002/fixture.00002.rst").unlink()
    with pytest.raises(fast.FastRunError, match="incomplete terminal"):
        fast.terminal_product_group(directory, ".rst", 3)


def test_terminal_restart_group_rejects_meshblock_boundary_truncation(fast, tmp_path):
    directory = tmp_path / "rst"
    write_rank_group(
        directory,
        "fixture.00001.rst",
        1,
        restart_payload(2.0, mesh_nx1=2, local_blocks=1),
    )

    with pytest.raises(fast.FastRunError, match="meshblock coverage"):
        fast.terminal_product_group(directory, ".rst", 1)


def test_terminal_restart_group_rejects_swapped_rank_block_counts(fast, tmp_path):
    directory = tmp_path / "rst"
    for rank, local_blocks in enumerate((3, 2)):
        path = directory / f"rank_{rank:08d}/fixture.00001.rst"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            restart_payload(2.0, mesh_nx1=5, local_blocks=local_blocks)
        )

    with pytest.raises(fast.FastRunError, match="meshblock coverage"):
        fast.terminal_product_group(directory, ".rst", 2)


def test_terminal_binary_group_selects_latest_physical_time(fast, tmp_path):
    directory = tmp_path / "bin"
    terminal = write_binary_rank_group(directory, "fixture.00001.bin", 3, 2.0)
    write_binary_rank_group(directory, "fixture.00002.bin", 3, 1.0)

    result = fast.terminal_product_group(directory, ".bin", 3)

    assert result["rank_zero"] == str(terminal.resolve())
    assert result["physical_time"] == 2.0

    (directory / "rank_00000002/fixture.00001.bin").write_bytes(b"corrupt")
    with pytest.raises(fast.FastRunError, match="magic/version"):
        fast.terminal_product_group(directory, ".bin", 3)


def test_terminal_binary_group_rejects_meshblock_boundary_truncation(fast, tmp_path):
    directory = tmp_path / "bin"
    write_rank_group(
        directory,
        "fixture.00001.bin",
        1,
        binary_payload(2.0, mesh_nx1=2, locations=(0,)),
    )

    with pytest.raises(fast.FastRunError, match="schema or logical coverage"):
        fast.terminal_product_group(directory, ".bin", 1)


def test_terminal_binary_group_rejects_cross_rank_schema_mismatch(fast, tmp_path):
    directory = tmp_path / "bin"
    write_binary_rank_group(directory, "fixture.00001.bin", 2, 2.0)
    rank_one = directory / "rank_00000001/fixture.00001.bin"
    rank_one.write_bytes(rank_one.read_bytes().replace(b"  cycle=1\n", b"  cycle=2\n"))

    with pytest.raises(fast.FastRunError, match="schema or logical coverage"):
        fast.terminal_product_group(directory, ".bin", 2)


def test_analyze_segment_accepts_complete_synchronized_clean_output(fast, tmp_path):
    segment = segment_fixture(tmp_path)

    result = fast.analyze_segment(segment, save=False)

    assert result["passed"] is True
    assert result["complete"] is True
    assert result["final_time"] == 1.0
    assert result["mass_relative_drift"] == 0.0
    assert result["mhd_user_mass_relative_mismatch"] == 0.0


def test_analyze_segment_reports_strict_lf_failure_without_accepting(fast, tmp_path):
    segment = segment_fixture(tmp_path, strict_failure=1.0)
    result = fast.analyze_segment(segment, save=False)
    assert result["passed"] is False
    assert result["strict_lf_failure_maxima"]["lf_dfloor"] == 1.0


def test_batch_script_pins_inputs_uses_unique_output_and_direct_auto_advance(
    fast, tmp_path
):
    manifest = batch_manifest(tmp_path)
    script = fast.batch_script_text(manifest)

    assert f'require_sha {fast.ATHENA_SHA256} "$ATHENA" executable' in script
    assert f"require_sha {manifest['input_sha256']} \"$INPUT\" input" in script
    assert f"require_sha {manifest['fast_script_sha256']} {LAUNCHER}" in script
    assert 'mkdir "$OUT_DIR"' in script
    assert f'"$ATHENA" -i {manifest["input"]} -d "$OUT_DIR"' in script
    assert 'advance --segment "$RUN_DIR" --submit' in script
    assert "cgl_lf_stage_i.py" not in script
    assert script.count("RUN_DIR=") == 1


def test_prepare_segment_creates_a_unique_no_clobber_run_directory(
    fast, tmp_path, monkeypatch
):
    monkeypatch.setattr(fast, "require_sha", lambda *_args: None)
    monkeypatch.setattr(fast, "sha256_file", lambda _path: "a" * 64)
    monkeypatch.setattr(fast, "utc_now", lambda: "2026-06-06T12:00:00+00:00")
    case = {"id": "R05", "name": "five", "input": "inputs/five.athinput"}

    segment = fast.prepare_segment(tmp_path, case, 0, 0.0, None, None)

    manifest = fast.load_json(fast.segment_manifest(segment))
    assert manifest["case_id"] == "R05"
    assert manifest["nodes"] == 4
    assert manifest["ranks"] == 32
    assert manifest["matrix_sha256"] == fast.MATRIX_SHA256
    assert manifest["executable_sha256"] == fast.ATHENA_SHA256
    assert manifest["job_id"] is None
    script = segment / "manifest/run.sbatch"
    assert script.stat().st_mode & stat.S_IXUSR
    marker = segment / "operator-note.txt"
    marker.write_text("retain\n", encoding="utf-8")

    with pytest.raises(FileExistsError):
        fast.prepare_segment(tmp_path, case, 0, 0.0, None, None)
    assert marker.read_text(encoding="utf-8") == "retain\n"


def test_submit_segment_invokes_direct_sbatch_and_records_job(
    fast, tmp_path, monkeypatch
):
    segment = submit_segment_fixture(tmp_path)
    observed: dict[str, object] = {}

    def fake_run(argv, **kwargs):
        observed["argv"] = argv
        observed["kwargs"] = kwargs
        return SimpleNamespace(stdout="12345;frontier\n")

    monkeypatch.setattr(fast.subprocess, "run", fake_run)
    monkeypatch.setattr(fast, "utc_now", lambda: "2026-06-06T12:00:00+00:00")

    assert fast.submit_segment(segment) == "12345"
    assert observed["argv"] == [
        "/usr/bin/sbatch",
        "--parsable",
        str(segment / "manifest/run.sbatch"),
    ]
    kwargs = observed["kwargs"]
    assert kwargs["check"] is True
    assert kwargs["capture_output"] is True
    assert kwargs["text"] is True
    assert kwargs.get("shell", False) is False
    manifest = fast.load_json(fast.segment_manifest(segment))
    assert manifest["job_id"] == "12345"
    assert manifest["submitted_utc"] == "2026-06-06T12:00:00+00:00"


@pytest.mark.parametrize(
    "output",
    ["", "0\n", "Submitted batch job 12345\n", "12345\n12346\n", "12;one;two\n"],
)
def test_submit_segment_rejects_ambiguous_sbatch_output(
    fast, tmp_path, monkeypatch, output
):
    segment = submit_segment_fixture(tmp_path)
    monkeypatch.setattr(
        fast.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(stdout=output),
    )
    with pytest.raises(fast.FastRunError, match="unexpected sbatch response"):
        fast.submit_segment(segment)


def test_submit_segment_refuses_resubmission(fast, tmp_path, monkeypatch):
    segment = submit_segment_fixture(tmp_path, job_id="12345")

    def forbidden(*_args, **_kwargs):
        raise AssertionError("sbatch must not run")

    monkeypatch.setattr(fast.subprocess, "run", forbidden)
    with pytest.raises(fast.FastRunError, match="already submitted"):
        fast.submit_segment(segment)


def test_job_state_prefers_live_squeue(fast, monkeypatch):
    calls: list[list[str]] = []

    def fake_run(argv, **_kwargs):
        calls.append(argv)
        return SimpleNamespace(stdout="RUNNING\n")

    monkeypatch.setattr(fast.subprocess, "run", fake_run)
    assert fast.job_state("12345") == "RUNNING"
    assert calls == [["/usr/bin/squeue", "-h", "-j", "12345", "-o", "%T"]]


def test_job_state_falls_back_to_sacct(fast, monkeypatch):
    replies = iter(["", "COMPLETED|\n"])
    calls: list[list[str]] = []

    def fake_run(argv, **_kwargs):
        calls.append(argv)
        return SimpleNamespace(stdout=next(replies))

    monkeypatch.setattr(fast.subprocess, "run", fake_run)
    assert fast.job_state("12345") == "COMPLETED"
    assert calls[1][0] == "/usr/bin/sacct"


def test_seed_for_case_uses_accepted_continuations_and_fresh_other_cases(fast):
    start, restart, digest = fast.seed_for_case("R03")
    assert start == 0.5
    assert restart == fast.CASE_SEEDS["R03"]["restart"]
    assert digest == fast.CASE_SEEDS["R03"]["sha256"]
    assert fast.seed_for_case("R12") == (0.0, None, None)


def test_next_from_segment_prepares_and_optionally_submits_clean_partial(
    fast, tmp_path, monkeypatch
):
    segment = tmp_path / "fast_s000_t0_to_t10"
    write_json(
        fast.segment_manifest(segment),
        {"root": str(tmp_path), "case_id": "R05", "sequence": 0},
    )
    restart = tmp_path / "rank_00000000/fixture.rst"
    prepared = tmp_path / "fast_s001_t3_to_t10"
    observed: dict[str, object] = {}
    monkeypatch.setattr(
        fast,
        "analyze_segment",
        lambda _segment: {
            "passed": True,
            "complete": False,
            "final_time": 3.0,
            "terminal_restart": {"rank_zero": str(restart), "sha256": "a" * 64},
        },
    )
    monkeypatch.setattr(
        fast,
        "campaign_cases",
        lambda: {"R05": {"id": "R05", "name": "five", "input": "five"}},
    )

    def fake_prepare(root, case, sequence, start_time, restart_path, restart_sha):
        observed["prepare"] = (
            root,
            case,
            sequence,
            start_time,
            restart_path,
            restart_sha,
        )
        return prepared

    monkeypatch.setattr(fast, "prepare_segment", fake_prepare)
    monkeypatch.setattr(
        fast,
        "submit_segment",
        lambda value: observed.setdefault("submitted", value) or "12345",
    )

    assert fast.next_from_segment(segment, submit=True) == prepared
    assert observed["prepare"][2:] == (1, 3.0, restart, "a" * 64)
    assert observed["submitted"] == prepared


def test_next_from_segment_stops_at_target_and_rejects_failed_science(
    fast, tmp_path, monkeypatch
):
    segment = tmp_path / "fast_s000_t0_to_t10"
    write_json(
        fast.segment_manifest(segment),
        {"root": str(tmp_path), "case_id": "R05", "sequence": 0},
    )
    monkeypatch.setattr(
        fast,
        "analyze_segment",
        lambda _segment: {"passed": True, "complete": True, "final_time": 10.0},
    )
    assert fast.next_from_segment(segment, submit=True) is None

    monkeypatch.setattr(
        fast,
        "analyze_segment",
        lambda _segment: {"passed": False, "complete": False, "final_time": 1.0},
    )
    with pytest.raises(fast.FastRunError, match="scientific sanity check failed"):
        fast.next_from_segment(segment, submit=False)
