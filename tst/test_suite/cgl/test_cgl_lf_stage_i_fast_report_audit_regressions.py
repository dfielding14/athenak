"""Regression tests for audited direct-fast Stage I report provenance blockers."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
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
    case_id: str = CASE_ID,
    case_name: str = CASE_NAME,
    source_relative: Path | None = None,
    start_time: float = 0.0,
    restart: Path | None = None,
    restart_sha256: str | None = None,
    job_id: str | None = None,
    run_exit_code: int | None = 0,
    variant: str | None = "standard",
    overrides: list[str] | None = None,
    claim_scope: str = "standard",
) -> Path:
    segment = (
        root
        / (source_relative or report.FAST_RUNS_RELATIVE)
        / case_id
        / segment_name
    )
    output = segment / "output"
    write_history(output / "fixture.mhd.hst", final_time)
    write_history(output / "fixture.user.hst", final_time)
    write_json(
        segment / "manifest/fast_run.json",
        {
            "case_id": case_id,
            "case_name": case_name,
            "sequence": sequence,
            "run_dir": str(segment.absolute()),
            "output_dir": str(output.absolute()),
            "start_time": start_time,
            "target_time": report.TARGET_TIME,
            "restart": str(restart.absolute()) if restart is not None else None,
            "restart_sha256": restart_sha256,
            "variant": variant,
            "command_line_overrides": overrides or [],
            "claim_scope": claim_scope,
            "job_id": job_id,
            "nodes": 1,
            "ranks": 1,
            "ranks_per_node": 1,
            "input_sha256": "c" * 64,
            "matrix_sha256": "d" * 64,
            "executable_sha256": "e" * 64,
        },
    )
    if run_exit_code is not None:
        (segment / "manifest/run_exit_code").write_text(
            f"{run_exit_code}\n", encoding="utf-8"
        )
    return segment


def write_restart(path: Path, time: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(f"<time>\ntime={time:.17g}\n<par_end>\n".encode())


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


def test_running_race_outranks_unsubmitted_partial_continuation(
    report, tmp_path, monkeypatch
):
    parent = write_fast_segment(
        report,
        tmp_path,
        segment_name="fast_s000_t0_to_t10",
        sequence=0,
        final_time=4.0,
    )
    restart = parent / "output/rst/rank_00000000/fixture.00001.rst"
    write_restart(restart, 4.0)
    write_fast_segment(
        report,
        tmp_path,
        segment_name="fast_s001_t4_to_t10",
        sequence=1,
        start_time=4.0,
        final_time=4.7,
        restart=restart,
        restart_sha256=sha256(restart),
        run_exit_code=None,
    )
    race = write_fast_segment(
        report,
        tmp_path,
        source_relative=report.RACE_RUNS_RELATIVE,
        segment_name="fast_s000_t0_to_t10",
        sequence=0,
        final_time=1.0,
        job_id="race-job",
        run_exit_code=None,
    )
    monkeypatch.setattr(
        report,
        "slurm_job_evidence",
        lambda job_id: {
            "job_id": job_id,
            "source": "squeue",
            "state": "RUNNING",
            "exit_code": None,
        },
    )

    selected, unselected, _warnings = report.select_fast_lineage(
        tmp_path, CASE_ID, CASE_NAME
    )

    assert Path(selected[-1]["segment"]) == race
    assert report.fast_candidate_state(selected[-1]) == "in_progress"
    assert any(
        record["terminal"]["state"] == "unsubmitted_partial"
        for record in unselected
        if isinstance(record.get("terminal"), dict)
    )


def test_r15_variant_selection_retains_and_classifies_strict_failures(
    report, tmp_path, monkeypatch
):
    strict_original = write_fast_segment(
        report,
        tmp_path,
        segment_name="fast_s000_t0_to_t10",
        sequence=0,
        final_time=0.2,
        case_id="R15",
        job_id="cancelled-strict",
        run_exit_code=None,
    )
    strict_race = write_fast_segment(
        report,
        tmp_path,
        source_relative=report.RACE_RUNS_RELATIVE,
        segment_name="fast_s000_t0_to_t10",
        sequence=0,
        final_time=1.2,
        case_id="R15",
        job_id="failed-strict",
        run_exit_code=143,
    )
    variant = write_fast_segment(
        report,
        tmp_path,
        source_relative=report.RELAXED_RUNS_RELATIVE,
        segment_name="fast_s000_t0_to_t10",
        sequence=0,
        final_time=0.0,
        case_id="R15",
        job_id="pending-variant",
        run_exit_code=None,
        variant=report.NONFATAL_HARD_BOUND_VARIANT,
        overrides=["mhd/cgl_lf_strict_admissibility=false"],
        claim_scope="R15 diagnostic variant",
    )
    states = {
        "cancelled-strict": "CANCELLED",
        "failed-strict": "FAILED",
        "pending-variant": "PENDING",
    }
    monkeypatch.setattr(
        report,
        "slurm_job_evidence",
        lambda job_id: {
            "job_id": job_id,
            "source": "squeue" if states[job_id] == "PENDING" else "sacct",
            "state": states[job_id],
            "exit_code": None,
        },
    )

    selected, unselected, _warnings = report.select_fast_lineage(
        tmp_path, "R15", CASE_NAME
    )
    strict_states = {
        Path(record["terminal"]["segment"]): record["terminal"]["state"]
        for record in unselected
        if isinstance(record.get("terminal"), dict)
    }

    assert Path(selected[-1]["segment"]) == variant
    assert Path(strict_original) in strict_states
    assert Path(strict_race) in strict_states
    assert strict_states[Path(strict_original)] == "failed"
    assert strict_states[Path(strict_race)] == "failed"


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


def test_r15_nonfatal_hard_bound_warning_is_case_scoped(report):
    mhd = {
        "time": [0.0, 1.0],
        "mass": [1.0, 1.0],
        "lf_dfloor": [0.0, 0.0],
        "lf_pfloor": [0.0, 0.0],
        "lf_nonfin": [0.0, 0.0],
        "lf_nonpos": [0.0, 0.0],
        "lf_hardbd": [0.0, 2.0],
    }
    user = {
        "time": [0.0, 1.0],
        "mass": [1.0, 1.0],
        "hard_vol": [0.0, 1.0],
    }

    health = report.compute_health(
        {
            "case_id": "R15",
            "status": "complete",
            "warnings": [],
            "errors": [],
            "lineage_variants": [report.NONFATAL_HARD_BOUND_VARIANT],
        },
        mhd,
        user,
        {"limiter_hardwall": "false"},
    )

    assert health["nonfatal_hard_bound_variant"] is True
    assert health["numerical_warnings"] == []
    assert any(
        "nonfatal R15 diagnostic retained hard-bound events" in warning
        for warning in health["science_warnings"]
    )


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


def test_snapshot_worker_plan_bounds_count_cpu_and_memory(report, tmp_path, monkeypatch):
    snapshots = [tmp_path / f"snapshot-{index}.bin" for index in range(5)]
    for snapshot in snapshots:
        with snapshot.open("wb") as stream:
            stream.truncate(1024 ** 3)
    exact = {str(snapshot): [snapshot] for snapshot in snapshots}
    monkeypatch.setattr(report.os, "sched_getaffinity", lambda _pid: set(range(3)))

    cpu_limited = report.snapshot_worker_plan(8, snapshots, exact, 384.0)
    memory_limited = report.snapshot_worker_plan(8, snapshots, exact, 50.0)

    assert cpu_limited["workers"] == 3
    assert memory_limited["workers"] == 2
    assert memory_limited["estimated_peak_bytes_per_process"] == 24 * 1024 ** 3


def test_snapshot_worker_plan_rejects_budget_below_one_worker(
    report, tmp_path, monkeypatch
):
    snapshot = tmp_path / "snapshot.bin"
    with snapshot.open("wb") as stream:
        stream.truncate(1024 ** 3)
    monkeypatch.setattr(report.os, "sched_getaffinity", lambda _pid: set(range(8)))

    with pytest.raises(report.ReportError, match="cannot fit one estimated"):
        report.snapshot_worker_plan(8, [snapshot], {str(snapshot): [snapshot]}, 23.0)


def test_pinned_analyzer_loads_exact_bytes_and_rejects_source_change(report, tmp_path):
    analyzer_path = tmp_path / "analyzer.py"
    analyzer_path.write_text("VALUE = 'pinned'\n", encoding="utf-8")
    pin = report.pin_pure_analyzer(analyzer_path)

    analyzer_path.write_text("VALUE = 'mutated'\n", encoding="utf-8")
    analyzer = report.load_pure_analyzer(pin)

    assert analyzer.VALUE == "pinned"
    assert getattr(analyzer, report.ANALYZER_DIGEST_ATTRIBUTE) == pin["sha256"]
    with pytest.raises(report.ReportError, match="changed after it was pinned"):
        report.verify_pure_analyzer_pin(pin)


def test_snapshot_worker_rejects_analyzer_payload_digest_mismatch(report):
    with pytest.raises(report.ReportError, match="payload digest mismatch"):
        report.analyzer_pin_from_task({
            "analyzer_path": "/unused/analyzer.py",
            "analyzer_source": b"VALUE = 1\n",
            "analyzer_sha256": "0" * 64,
        })


class FakeSnapshotAnalyzer:
    np = np

    def __init__(self, times):
        self.times = times
        self.serial_calls = 0
        self.serial_arguments = None

    def snapshot_sibling_paths(self, path, _expected_ranks):
        return [path]

    def snapshot_time(self, path):
        return self.times[str(path)]

    def time_mask(self, times, start, end):
        return (times >= start) & (times <= end)

    def snapshot_digest_provenance(self, path, _expected_ranks, _exact_rank_set):
        return {"path": str(path), "sha256": sha256(path)}

    def read_snapshot(self, path):
        time = self.times[str(path)]
        values = np.asarray([time, time + 1.0])
        return {"values": values}, (1.0, 1.0, 1.0), time

    def pdf_fields(self, fields, _lengths):
        return {"value": fields["values"]}

    def pressure_density_fields(self, fields):
        return {"joint": (fields["values"], 2.0 * fields["values"])}

    def average_snapshot_records(self, records):
        return {"snapshot_count": len(records), "ordered_paths": list(records)}

    def analyze_snapshot_paths(
        self,
        paths,
        bins,
        alignment_shells,
        time_start,
        time_end,
        model_choices,
        eddy_samples,
        eddy_bins,
        eddy_seed,
    ):
        self.serial_calls += 1
        self.serial_arguments = (
            paths,
            bins,
            alignment_shells,
            time_start,
            time_end,
            model_choices,
            eddy_samples,
            eddy_bins,
            eddy_seed,
        )
        return {"serial": {}}, {"snapshot_count": 1, "serial": True}


def test_parallel_snapshot_orchestration_preserves_order_ranges_and_provenance(
    report, tmp_path, monkeypatch
):
    snapshots = [tmp_path / "first.bin", tmp_path / "second.bin"]
    for index, snapshot in enumerate(snapshots):
        snapshot.write_bytes(f"snapshot-{index}".encode())
    times = {str(snapshots[0]): 8.0, str(snapshots[1]): 9.0}
    analyzer = FakeSnapshotAnalyzer(times)
    monkeypatch.setattr(report.os, "sched_getaffinity", lambda _pid: set(range(8)))
    captured = {}

    def runner(tasks, workers):
        captured["tasks"] = tasks
        captured["workers"] = workers
        return [
            (task["path"], {"time": times[task["path"]], "worker": workers})
            for task in tasks
        ]

    def range_runner(tasks, workers):
        captured["range_tasks"] = tasks
        captured["range_workers"] = workers
        return [
            (
                task["path"],
                report.snapshot_range_record(
                    analyzer,
                    Path(task["path"]),
                    [Path(value) for value in task["exact_rank_set"]],
                ),
            )
            for task in tasks
        ]

    records, ensemble = report.analyze_snapshot_paths_bounded(
        analyzer,
        snapshots,
        64,
        [2, 4],
        8.0,
        10.0,
        {},
        2_000_000,
        24,
        731,
        {str(path): 1 for path in snapshots},
        2,
        384.0,
        worker_runner=runner,
        range_runner=range_runner,
    )

    assert captured["workers"] == 2
    assert captured["range_workers"] == 2
    worker_pins = {
        (task["analyzer_sha256"], task["analyzer_source"])
        for task in captured["tasks"] + captured["range_tasks"]
    }
    assert len(worker_pins) == 1
    digest, source = worker_pins.pop()
    assert digest == hashlib.sha256(source).hexdigest()
    assert captured["tasks"][0]["ranges"] == {"value": (8.0, 10.0)}
    assert captured["tasks"][0]["joint_ranges"] == {
        "joint": ((8.0, 10.0), (16.0, 20.0))
    }
    assert list(records) == [str(path) for path in snapshots]
    assert ensemble["ordered_paths"] == [str(path) for path in snapshots]
    assert records[str(snapshots[0])]["snapshot_provenance"]["sha256"] == sha256(
        snapshots[0]
    )
    assert analyzer.serial_calls == 0


def test_single_snapshot_worker_uses_exact_legacy_analyzer_path(
    report, tmp_path, monkeypatch
):
    snapshot = tmp_path / "snapshot.bin"
    snapshot.write_bytes(b"snapshot")
    analyzer = FakeSnapshotAnalyzer({str(snapshot): 8.0})
    monkeypatch.setattr(report.os, "sched_getaffinity", lambda _pid: set(range(8)))

    records, ensemble = report.analyze_snapshot_paths_bounded(
        analyzer,
        [snapshot],
        64,
        [2],
        8.0,
        10.0,
        {},
        2_000_000,
        24,
        731,
        {str(snapshot): 1},
        4,
        384.0,
    )

    assert records == {"serial": {}}
    assert ensemble == {"snapshot_count": 1, "serial": True}
    assert analyzer.serial_calls == 1
    assert len(analyzer.serial_arguments) == 9


def test_parallel_snapshot_orchestration_rejects_analyzer_source_change(
    report, tmp_path, monkeypatch
):
    snapshots = [tmp_path / "first.bin", tmp_path / "second.bin"]
    for snapshot in snapshots:
        snapshot.write_bytes(snapshot.name.encode())
    times = {str(snapshots[0]): 8.0, str(snapshots[1]): 9.0}
    analyzer = FakeSnapshotAnalyzer(times)
    analyzer_path = tmp_path / "analyzer.py"
    analyzer_path.write_text("VALUE = 'pinned'\n", encoding="utf-8")
    analyzer_pin = report.pin_pure_analyzer(analyzer_path)
    monkeypatch.setattr(report.os, "sched_getaffinity", lambda _pid: set(range(8)))

    def range_runner(tasks, _workers):
        records = [
            (
                task["path"],
                report.snapshot_range_record(
                    analyzer,
                    Path(task["path"]),
                    [Path(value) for value in task["exact_rank_set"]],
                ),
            )
            for task in tasks
        ]
        analyzer_path.write_text("VALUE = 'mutated'\n", encoding="utf-8")
        return records

    with pytest.raises(report.ReportError, match="changed after it was pinned"):
        report.analyze_snapshot_paths_bounded(
            analyzer,
            snapshots,
            64,
            [2],
            8.0,
            10.0,
            {},
            0,
            24,
            731,
            {str(path): 1 for path in snapshots},
            2,
            384.0,
            worker_runner=lambda tasks, _workers: [
                (task["path"], {"time": times[task["path"]]})
                for task in tasks
            ],
            range_runner=range_runner,
            analyzer_pin=analyzer_pin,
        )


def test_parallel_snapshot_orchestration_rejects_out_of_order_results(
    report, tmp_path, monkeypatch
):
    snapshots = [tmp_path / "first.bin", tmp_path / "second.bin"]
    for snapshot in snapshots:
        snapshot.write_bytes(snapshot.name.encode())
    times = {str(snapshots[0]): 8.0, str(snapshots[1]): 9.0}
    analyzer = FakeSnapshotAnalyzer(times)
    monkeypatch.setattr(report.os, "sched_getaffinity", lambda _pid: set(range(8)))

    with pytest.raises(ValueError, match="out of order"):
        report.analyze_snapshot_paths_bounded(
            analyzer,
            snapshots,
            64,
            [2],
            8.0,
            10.0,
            {},
            0,
            24,
            731,
            {str(path): 1 for path in snapshots},
            2,
            384.0,
            worker_runner=lambda tasks, _workers: [
                (task["path"], {"time": times[task["path"]]})
                for task in reversed(tasks)
            ],
            range_runner=lambda tasks, _workers: [
                (
                    task["path"],
                    report.snapshot_range_record(
                        analyzer,
                        Path(task["path"]),
                        [Path(value) for value in task["exact_rank_set"]],
                    ),
                )
                for task in tasks
            ],
        )
