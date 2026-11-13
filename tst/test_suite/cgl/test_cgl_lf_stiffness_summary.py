"""Stiffness summaries include every rank and retain sparse worst-cell evidence."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


SCRIPT = Path(__file__).resolve().parents[3] / "scripts/summarize_cgl_lf_stiffness.py"
spec = importlib.util.spec_from_file_location("stiffness_summary", SCRIPT)
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)


def write_map(path, values, rank=0, nranks=1, ids=None, dtype="<f8", **overrides):
    values = np.asarray(values, dtype=dtype).reshape(-1, 1, 1, 4)
    maximum = np.finfo(values.dtype).max
    meta = dict(cycle=7, phase="pre", rank=rank, nranks=nranks,
                nblocks=len(values), total_blocks=len(values)*nranks,
                nx1=4, nx2=1, nx3=1, real_bytes=values.dtype.itemsize,
                endian="big" if dtype[0] == ">" else "little",
                time=2., dt_outer=2., dt_sweep=1., sts_safety=.5,
                dt_sts=.125, nstages=7,
                raw_dt_rank=min(maximum, float(values.min())))
    meta.update(overrides)
    ids = list(range(rank*len(values), (rank+1)*len(values))) if ids is None else ids
    with path.open("wb") as stream:
        stream.write(b"# AthenaK LF stiffness v1\n")
        stream.write((" ".join(f"{key}={value}" for key, value in meta.items())+"\n")
                     .encode())
        for gid in ids:
            stream.write(f"block {gid} 0 {gid} 0 0\n".encode())
        stream.write(b"data\n")
        stream.write(values.tobytes())
    return path


@pytest.mark.parametrize("dtype", ["<f8", ">f4"])
def test_stiffness_distribution_and_block_extremes(tmp_path, dtype):
    paths = [write_map(tmp_path / f"rank{rank}", values, rank, 2, dtype=dtype)
             for rank, values in enumerate(([1., 2., 4., 8.], [.25, 16., 32., 64.]))]
    result = summary.summarize(paths[::-1])
    assert result["owned_cells"] == 8
    assert result["maximum_stiffness"] == 8
    assert result["cells_within_factor_ten"] == 3
    assert result["blocks_within_factor_ten"] == 2
    assert result["worst_block"]["gid"] == 1
    assert result["worst_block"]["cell_ijk"] == [0, 0, 0]
    assert result["next_worst_block"]["gid"] == 0
    assert result["predicted_rkl_stages"] == result["selected_rkl_stages"] == 7
    assert result["predicted_rkl_stages_without_worst_block"] == 3
    assert result["median_stiffness_interval"][0] <= .25
    assert result["median_stiffness_interval"][1] >= .25
    assert result["p99_stiffness_interval"][0] <= 8
    assert not result["warnings"]


def test_header_round_trips_at_writer_precision(tmp_path):
    minimum = np.float32(0.123456)
    path = write_map(tmp_path / "float32", [minimum, 1, 2, 3], dtype="<f4",
                     raw_dt_rank=f"{minimum:.9g}")
    assert summary.summarize([path])["worst_block"]["raw_dt"] == float(minimum)


def test_inactive_and_invalid_cells_have_distinct_meanings(tmp_path):
    path = write_map(tmp_path / "map", [0., 1., np.inf, np.inf])
    result = summary.summarize([path])
    assert result["invalid_bound_cells"] == 1
    assert result["inactive_cells"] == 2
    assert result["maximum_stiffness"] == "infinite"
    assert result["median_stiffness_interval"] == [0., 0.]
    assert result["p99_stiffness_interval"] == ["infinite", "infinite"]
    assert result["cells_within_factor_ten"] == 1
    json.dumps(result, allow_nan=False)
    path = write_map(path, [np.inf]*4)
    result = summary.summarize([path])
    assert result["maximum_stiffness"] == result["predicted_rkl_stages"] == 0
    assert result["cells_within_factor_ten"] == 0
    assert result["worst_block"]["raw_dt"] == "inactive"


def test_budget_difference_is_reported_even_when_stage_count_agrees(tmp_path):
    path = write_map(tmp_path / "map", [.25, 1., 2., 3.], dt_sts=.12)
    result = summary.summarize([path])
    assert result["selected_rkl_stages"] == result["predicted_rkl_stages"] == 7
    assert result["warnings"] == [
        "selected dt_sts differs from sts_safety*minimum raw timestep"]


@pytest.mark.parametrize("change", ["missing_rank", "duplicate_rank", "mixed_phase",
                                   "duplicate_block", "truncated", "nan", "negative",
                                   "wrong_minimum"])
def test_invalid_map_sets_fail(tmp_path, change):
    paths = [write_map(tmp_path / f"rank{rank}", [.25, 1., 2., 3.], rank, 2)
             for rank in range(2)]
    if change == "missing_rank":
        paths.pop()
    elif change == "duplicate_rank":
        paths[1] = paths[0]
    elif change == "mixed_phase":
        write_map(paths[1], [.25, 1., 2., 3.], 1, 2, phase="post")
    elif change == "duplicate_block":
        write_map(paths[1], [.25, 1., 2., 3.], 1, 2, ids=[0])
    elif change == "truncated":
        paths[0].write_bytes(paths[0].read_bytes()[:-1])
    elif change == "wrong_minimum":
        write_map(paths[0], [.25, 1., 2., 3.], 0, 2, raw_dt_rank=.5)
    else:
        value = np.nan if change == "nan" else -1.
        write_map(paths[0], [.25, value, 2., 3.], 0, 2, raw_dt_rank=.25)
    with pytest.raises(ValueError):
        summary.summarize(paths)


def test_cli_outputs_strict_json_and_stage_formula_handles_boundaries(tmp_path):
    path = write_map(tmp_path / "map", [.25, 1., 2., 3.])
    result = subprocess.run([sys.executable, str(SCRIPT), str(path)],
                            capture_output=True, text=True, check=True)
    assert json.loads(result.stdout)["maximum_stiffness"] == 8
    assert [summary.rkl_stages(value) for value in [0., 2.5, 3., 7., 7.5]] == [0, 3, 5, 5, 7]
    assert summary.rkl_stages(np.inf) == "exceeds int32"


def test_many_rank_files_use_bounded_open_file_count(tmp_path):
    resource = pytest.importorskip("resource")
    paths = [write_map(tmp_path / f"rank{rank}", [.25, 1., 2., 3.], rank, 96)
             for rank in range(96)]

    def limit_files():
        _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (32, hard))

    result = subprocess.run([sys.executable, str(SCRIPT), *map(str, paths)],
                            capture_output=True, text=True, check=True,
                            preexec_fn=limit_files)
    assert json.loads(result.stdout)["blocks"] == 96
