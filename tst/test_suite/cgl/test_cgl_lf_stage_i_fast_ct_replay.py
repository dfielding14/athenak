"""Focused tests for exact-t=9 direct-fast CT replay preparation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
REPLAY_PATH = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_ct_replay.py"


@pytest.fixture(scope="module")
def replay():
    spec = importlib.util.spec_from_file_location("cgl_ct_replay_tests", REPLAY_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")


def write_restart(path: Path, time_value: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        (
            "<time>\n"
            f"restart_time = {time_value:.17g}\n"
            "<par_end>\n"
        ).encode("ascii")
    )


def selected_case(tmp_path: Path) -> tuple[dict[str, object], str]:
    matrix_sha = "a" * 64
    segment = tmp_path / "fast_s000_t0_to_t10"
    output = segment / "output"
    manifest_path = segment / "manifest/fast_run.json"
    write_json(manifest_path, {
        "case_id": "R03",
        "case_name": "fixture",
        "matrix_sha256": matrix_sha,
        "ranks": 1,
        "output_dir": str(output.resolve()),
    })
    for index, time_value in enumerate((8.0, 8.75, 9.0001, 10.0)):
        write_restart(
            output / "rst/rank_00000000" / f"fixture.{index:05d}.rst",
            time_value,
        )
    case = {
        "case_id": "R03",
        "case_name": "fixture",
        "lineage": [{
            "kind": "fast",
            "order": 0,
            "case_id": "R03",
            "case_name": "fixture",
            "matrix_sha256": matrix_sha,
            "ranks": 1,
            "segment": segment.name,
            "output": str(output.resolve()),
            "manifest": binding(manifest_path),
        }],
    }
    return case, matrix_sha


def test_select_parent_uses_latest_complete_state_strictly_before_t9(
    replay, tmp_path
) -> None:
    case, matrix_sha = selected_case(tmp_path)

    selected = replay.select_parent_group("R03", case, matrix_sha)

    assert selected["time"] == 8.75
    assert selected["name"] == "fixture.00001.rst"
    assert selected["rank_count"] == 1
    assert len(selected["rank_files"]) == 1
    assert selected["rank_files"][0]["rank"] == 0


def test_terminal_failures_are_excluded_from_default_replay_set(replay) -> None:
    present = {f"R{number:02d}": {} for number in range(2, 18)}

    selected = replay.expand_cases([], present)

    assert "R14" not in selected
    assert "R15" not in selected
    assert selected == [
        case_id for case_id in sorted(present) if case_id not in {"R14", "R15"}
    ]
