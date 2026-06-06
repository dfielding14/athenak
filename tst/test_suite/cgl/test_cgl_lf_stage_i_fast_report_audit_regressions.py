"""Regression tests for audited direct-fast Stage I report provenance blockers."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
REPORTER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_report.py"
CASE_ID = "R03"
CASE_NAME = "paper_standard_active_alfvenic_beta100"


def load_reporter():
    name = "cgl_lf_stage_i_fast_report_audit_regressions"
    spec = importlib.util.spec_from_file_location(name, REPORTER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def report():
    return load_reporter()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_history(path: Path, final_time: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# [1]=time [2]=mass\n"
        f"0 1\n{final_time:.17g} 1\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fast_segment(
    report,
    root: Path,
    *,
    segment_name: str,
    sequence: object,
    final_time: float,
    case_name: str = CASE_NAME,
    restart: Path | None = None,
    restart_sha256: str | None = None,
    overrides: list[str] | None = None,
) -> Path:
    segment = root / report.FAST_RUNS_RELATIVE / CASE_ID / segment_name
    output = segment / "output"
    write_history(output / "fixture.mhd.hst", final_time)
    write_history(output / "fixture.user.hst", final_time)
    write_json(
        segment / "manifest/fast_run.json",
        {
            "case_id": CASE_ID,
            "case_name": case_name,
            "sequence": sequence,
            "run_dir": str(segment.absolute()),
            "output_dir": str(output.absolute()),
            "start_time": 0.0,
            "target_time": report.TARGET_TIME,
            "restart": str(restart.absolute()) if restart is not None else None,
            "restart_sha256": restart_sha256,
            "variant": "standard",
            "command_line_overrides": overrides or [],
            "claim_scope": "standard",
            "nodes": 1,
            "ranks": 1,
            "ranks_per_node": 1,
            "input_sha256": "c" * 64,
            "matrix_sha256": "d" * 64,
            "executable_sha256": "e" * 64,
        },
    )
    (segment / "manifest/run_exit_code").write_text("0\n", encoding="utf-8")
    return segment


@pytest.mark.parametrize("invalid_chain", ["historical_seed", "fast_restart"])
def test_invalid_restart_chain_cannot_be_selected_complete(
    report, tmp_path, invalid_chain
):
    if invalid_chain == "historical_seed":
        missing = tmp_path / "missing-seed.rst"
        write_fast_segment(
            report,
            tmp_path,
            segment_name="fast_s000_t0_to_t10",
            sequence=0,
            final_time=10.0,
            restart=missing,
            restart_sha256="a" * 64,
        )
    else:
        parent = write_fast_segment(
            report,
            tmp_path,
            segment_name="fast_s000_t0_to_t10",
            sequence=0,
            final_time=1.0,
        )
        missing = parent / "output/rst/rank_00000000/missing-parent.rst"
        write_fast_segment(
            report,
            tmp_path,
            segment_name="fast_s001_t1_to_t10",
            sequence=1,
            final_time=10.0,
            restart=missing,
            restart_sha256="b" * 64,
        )

    selected, _unselected, _warnings = report.select_fast_lineage(
        tmp_path, CASE_ID, CASE_NAME
    )

    assert not selected or report.fast_candidate_state(selected[-1]) != "complete"


def test_override_order_and_repetition_are_exact_lineage_identity(report):
    def configuration(overrides: list[str]):
        return report.fast_candidate_configuration(
            {
                "manifest": {
                    "variant": "standard",
                    "command_line_overrides": overrides,
                }
            }
        )

    forward = configuration(["mhd/value=first", "mhd/value=second"])
    reversed_order = configuration(["mhd/value=second", "mhd/value=first"])
    repeated = configuration(
        ["mhd/value=first", "mhd/value=second", "mhd/value=first"]
    )

    assert forward != reversed_order
    assert forward != repeated
    assert reversed_order != repeated


def test_assembly_preserves_ordered_repeated_overrides(
    report, tmp_path, monkeypatch
):
    overrides = ["mhd/value=first", "mhd/value=second", "mhd/value=first"]
    frozen_source = tmp_path / "frozen"
    expected_input = frozen_source / "inputs/fixture.athinput"
    expected_input.parent.mkdir(parents=True)
    expected_input.write_text("<mhd>\n", encoding="utf-8")
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        report,
        "select_fast_lineage",
        lambda _root, _case_id, _case_name: (
            [
                {
                    "manifest": {
                        "restart": None,
                        "command_line_overrides": overrides,
                    }
                }
            ],
            [],
            [],
        ),
    )
    monkeypatch.setattr(
        report,
        "fast_segment_record",
        lambda _item, _order: {
            "state": "complete",
            "variant": "standard",
            "claim_scope": "standard",
            "command_line_overrides": overrides,
            "output": str(tmp_path / "unused-output"),
            "ranks": 0,
        },
    )
    monkeypatch.setattr(report, "fast_candidate_summary", lambda _item: {})

    def record_model_choices(_path, effective_overrides):
        captured["overrides"] = effective_overrides
        return {}

    monkeypatch.setattr(report, "model_choices_for_input", record_model_choices)
    monkeypatch.setattr(
        report,
        "merge_histories",
        lambda *_args, **_kwargs: {
            "available": False,
            "warnings": [],
            "errors": [],
        },
    )
    monkeypatch.setattr(
        report,
        "index_snapshots",
        lambda *_args, **_kwargs: {
            "snapshot_count": 0,
            "complete_snapshot_count": 0,
            "warnings": [],
        },
    )

    record = report.assemble_fast_case(
        tmp_path,
        frozen_source,
        tmp_path / "analysis",
        CASE_ID,
        {"name": CASE_NAME, "input": "inputs/fixture.athinput"},
    )

    assert record["lineage_command_line_overrides"] == overrides
    assert captured["overrides"] == overrides


@pytest.mark.parametrize("sequence", [1.5, 1])
def test_manifest_sequence_must_be_integral_and_match_segment_name(
    report, tmp_path, sequence
):
    write_fast_segment(
        report,
        tmp_path,
        segment_name="fast_s000_t0_to_t10",
        sequence=sequence,
        final_time=10.0,
    )
    rejections: list[dict[str, object]] = []

    candidates = report.fast_candidates(tmp_path, CASE_ID, CASE_NAME, rejections)

    assert candidates == []
    assert any("sequence" in str(item.get("reason")) for item in rejections)


def test_manifest_case_name_must_match_matrix_case(
    report, tmp_path, monkeypatch
):
    frozen_source = tmp_path / "frozen"
    expected_input = frozen_source / "inputs/fixture.athinput"
    expected_input.parent.mkdir(parents=True)
    expected_input.write_text("<mhd>\n", encoding="utf-8")
    write_fast_segment(
        report,
        tmp_path,
        segment_name="fast_s000_t0_to_t10",
        sequence=0,
        final_time=10.0,
        case_name="wrong_case_name",
    )
    monkeypatch.setattr(report, "model_choices_for_input", lambda *_args: {})

    record = report.assemble_fast_case(
        tmp_path,
        frozen_source,
        tmp_path / "analysis",
        CASE_ID,
        {"name": CASE_NAME, "input": "inputs/fixture.athinput"},
    )

    assert record["status"] != "complete" or any(
        "case_name" in error for error in record["errors"]
    )


def write_verification_fixture(report, output: Path, defect: str) -> None:
    case_dir = output / "cases" / CASE_ID
    mhd = case_dir / "history/fixture.mhd.hst"
    user = case_dir / "history/fixture.user.hst"
    write_history(mhd, 10.0)
    write_history(user, 10.0)
    lineage = {
        "case_id": CASE_ID,
        "case_name": CASE_NAME,
        "status": "complete",
        "errors": [],
        "warnings": [],
        "lineage": [{"kind": "fast", "run_exit_code": 0}],
        "histories": {
            "mhd": {
                "available": True,
                "path": str(mhd.absolute()),
                "binding": {"sha256": sha256(mhd)},
            },
            "user": {
                "available": True,
                "path": str(user.absolute()),
                "binding": {"sha256": sha256(user)},
            },
        },
    }
    diagnostics = {
        "analysis_errors": [],
        "analysis_warnings": [],
        "health": {
            "structural_errors": [],
            "structural_warnings": [],
            "numerical_warnings": [],
            "science_warnings": [],
        },
    }
    if defect == "nonzero_exit":
        lineage["lineage"][0]["run_exit_code"] = 17
    elif defect == "provenance_error":
        lineage["errors"] = ["fixture provenance identity mismatch"]
    elif defect == "diagnostic_structural_error":
        diagnostics["health"]["structural_errors"] = [
            "fixture diagnostic structural error"
        ]
    else:
        raise AssertionError(f"unknown fixture defect: {defect}")

    write_json(output / "inventory.json", {"cases": {CASE_ID: lineage}})
    write_json(case_dir / "lineage.json", lineage)
    write_json(
        case_dir / "snapshots.json",
        {
            "snapshot_count": 0,
            "complete_snapshot_count": 0,
            "snapshots": [],
        },
    )
    write_json(case_dir / "diagnostics.json", diagnostics)


@pytest.mark.parametrize(
    "defect",
    ["nonzero_exit", "provenance_error", "diagnostic_structural_error"],
)
def test_release_verification_fails_provenance_and_structural_defects(
    report, tmp_path, defect
):
    write_verification_fixture(report, tmp_path, defect)

    result = report.command_verify(
        SimpleNamespace(
            output=tmp_path,
            cases=[CASE_ID],
            require_complete=True,
        )
    )
    verification = json.loads((tmp_path / "verify.json").read_text(encoding="utf-8"))

    assert result == 1
    assert verification["result"] == "fail"
    assert verification["errors"]
