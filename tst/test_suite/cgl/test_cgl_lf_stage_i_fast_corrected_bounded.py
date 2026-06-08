"""Tests for short native-checkpoint continuations of R14 and R15."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
LAUNCHER = (
    REPOSITORY
    / "scripts/frontier/cgl_lf_stage_i_fast_corrected_bounded.py"
)


def load_launcher():
    name = "cgl_lf_stage_i_fast_corrected_bounded_for_tests"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def bounded():
    return load_launcher()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_scope_and_walltimes_are_explicit(bounded):
    assert bounded.BOUNDED_CASES == ("R14", "R15")
    assert bounded.REQUESTED_WALLTIME == "00:40:00"
    assert bounded.ATHENA_WALLTIME == "00:30:00"
    bounded.require_bounded_case("R14")
    bounded.require_bounded_case("R15")
    with pytest.raises(bounded.BoundedContinuationError, match="only R14, R15"):
        bounded.require_bounded_case("R12")


def test_configure_installs_wrapper_identity_and_short_walltimes(
    bounded, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    monkeypatch.setattr(bounded.corrected, "CAMPAIGN_ROOT", root)

    assert bounded.configure(root) == root.resolve()
    assert bounded.corrected.SCRIPT_PATH == LAUNCHER.resolve()
    assert bounded.corrected.fast.__file__ == str(LAUNCHER.resolve())
    assert bounded.corrected.fast.REQUESTED_WALLTIME == "00:40:00"
    assert bounded.corrected.fast.ATHENA_WALLTIME == "00:30:00"
    assert bounded.corrected.fast.prepare_segment is bounded.prepare_segment


def test_prepare_marks_only_runtime_segmentation_without_changing_physics(
    bounded, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    segment = (
        root
        / bounded.corrected.RUNS_RELATIVE
        / "R14/fast_s004_t6p3_to_t10"
    )

    def fake_prepare(_root, case, sequence, start_time, restart, restart_sha):
        manifest = {
            "case_id": case["id"],
            "sequence": sequence,
            "start_time": start_time,
            "restart": str(restart),
            "restart_sha256": restart_sha,
            "fast_script": str(LAUNCHER.resolve()),
            "variant": bounded.corrected.FINITE_LIMITER_VARIANT,
            "command_line_overrides": [
                bounded.corrected.FINITE_LIMITER_OVERRIDE
            ],
        }
        write_json(bounded.corrected.fast.segment_manifest(segment), manifest)
        script = segment / "manifest/run.sbatch"
        script.write_text(
            "#SBATCH -t 00:40:00\n"
            f"{LAUNCHER.resolve()}\n"
            '"$ATHENA" -r restart -t 00:30:00 '
            "time/tlim=10.0 "
            f"{bounded.corrected.FINITE_LIMITER_OVERRIDE}\n",
            encoding="utf-8",
        )
        return segment

    monkeypatch.setattr(bounded, "_CORRECTED_PREPARE_SEGMENT", fake_prepare)
    monkeypatch.setattr(
        bounded.corrected,
        "validate_segment",
        lambda path: bounded.corrected.fast.load_json(
            bounded.corrected.fast.segment_manifest(path)
        ),
    )
    restart = root / "runs/R14/s003/output/rst/rank_00000000/r.rst"
    result = bounded.prepare_segment(
        root,
        {"id": "R14"},
        4,
        6.3,
        restart,
        "a" * 64,
    )
    manifest = bounded.validate_bounded_segment(result)

    assert manifest["runtime_segmentation_policy"] == (
        "native_wallclock_checkpoint_30m_v1"
    )
    assert manifest["runtime_segmentation_changes_physics"] is False
    assert manifest["variant"] == bounded.corrected.FINITE_LIMITER_VARIANT
    assert manifest["command_line_overrides"] == [
        bounded.corrected.FINITE_LIMITER_OVERRIDE
    ]


def test_bounded_prepare_requires_existing_continuation(bounded, tmp_path):
    with pytest.raises(
        bounded.BoundedContinuationError,
        match="requires an authenticated continuation",
    ):
        bounded.prepare_segment(
            tmp_path,
            {"id": "R14"},
            0,
            0.0,
            None,
            None,
        )
    with pytest.raises(bounded.BoundedContinuationError, match="only R14, R15"):
        bounded.prepare_segment(
            tmp_path,
            {"id": "R12"},
            1,
            1.0,
            tmp_path / "r.rst",
            "a" * 64,
        )


def test_launch_refuses_to_create_a_fresh_lineage(
    bounded, tmp_path, monkeypatch
):
    monkeypatch.setattr(
        bounded.corrected.fast,
        "segment_directories",
        lambda _root, _case_id: [],
    )
    with pytest.raises(bounded.BoundedContinuationError, match="no corrected lineage"):
        bounded.launch_case(tmp_path, "R14", submit=False)
