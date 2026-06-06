"""Focused regressions for the direct-fast parallel analysis launcher."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
LAUNCHER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_analyze.py"


def load_launcher():
    name = "cgl_lf_stage_i_fast_analyze_tests"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fast_analyze():
    return load_launcher()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory_fixture(
    tmp_path: Path,
    *,
    analysis: Path | None = None,
    declared_output: Path | None = None,
) -> tuple[Path, Path, Path]:
    root = tmp_path / "campaign"
    output = analysis or tmp_path / "analysis"
    inventory = output / "inventory.json"
    write_json(
        inventory,
        {
            "schema_version": 1,
            "root": str(root),
            "output": str(declared_output or output),
            "cases": {
                "R02": {
                    "case_id": "R02",
                    "case_name": "case-r02",
                    "status": "complete",
                },
                "R03": {
                    "case_id": "R03",
                    "case_name": "case-r03",
                    "status": "in_progress",
                },
            },
        },
    )
    return root, output, inventory


def command_args(
    analysis: Path,
    jobs: Path,
    *,
    cases: list[str] | None = None,
    include_all: bool = False,
    eddy_samples: int | None = 1234,
) -> SimpleNamespace:
    return SimpleNamespace(
        analysis=analysis,
        jobs_dir=jobs,
        cases=cases or [],
        all=include_all,
        account="ast207",
        partition="batch",
        walltime="02:00:00",
        cpus_per_task=8,
        python=Path(sys.executable),
        submit=False,
        skip_snapshots=False,
        snapshot_time_start=8.0,
        snapshot_time_end=10.0,
        snapshot_workers=4,
        snapshot_memory_budget_gib=384.0,
        pdf_bins=64,
        alignment_shells="2,4,6",
        eddy_samples=eddy_samples,
        eddy_bins=24,
        eddy_seed=731,
    )


def write_complete_analysis(
    fast_analyze, analysis: Path, jobs: Path, case_id: str = "R02"
) -> tuple[Path, dict[str, object]]:
    attempt = jobs / case_id / "attempt-000"
    manifest = fast_analyze.load_json(attempt / "manifest.json")
    inventory = fast_analyze.load_json(analysis / "inventory.json")
    case = inventory["cases"][case_id]
    lineage = analysis / "cases" / case_id / "lineage.json"
    write_json(lineage, case)
    diagnostics = {
        "schema_version": 1,
        "analyzed_utc": fast_analyze.utc_now(),
        "case_id": case_id,
        "case_name": case["case_name"],
        "assembly_status": "complete",
        "analysis_status": "complete",
        "snapshot_analysis_status": "complete",
        "analysis_errors": [],
        "provenance": {
            "adapter": manifest["reporter"],
            "lineage": fast_analyze.artifact_binding(lineage),
        },
    }
    write_json(analysis / "cases" / case_id / "diagnostics.json", diagnostics)
    (attempt / "exit_code.txt").write_text("0\n", encoding="utf-8")
    return attempt, diagnostics


def test_launch_binds_inventory_and_declared_analysis_output(
    fast_analyze, tmp_path
):
    root, analysis, inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"

    assert fast_analyze.launch(command_args(inventory, jobs, cases=["R02"])) == 0

    manifest = fast_analyze.load_json(jobs / "R02/attempt-000/manifest.json")
    assert manifest["analysis_output"] == str(analysis.resolve())
    assert manifest["inventory"]["path"] == str(inventory.resolve())
    assert manifest["inventory"]["sha256"] == sha256(inventory)
    assert manifest["command"][2:5] == [
        "--output",
        str(analysis.resolve()),
        "analyze-case",
    ]
    assert manifest["command"][5] == "R02"
    assert "--snapshot-workers" in manifest["report_options"]
    assert "4" in manifest["report_options"]
    assert "--snapshot-memory-budget-gib" in manifest["report_options"]
    assert "384.0" in manifest["report_options"]
    assert manifest["nodes"] == 1
    assert manifest["job_id"] is None
    assert not (root / "runs").exists()


def test_inventory_declared_output_mismatch_rejected_without_job_writes(
    fast_analyze, tmp_path
):
    _root, analysis, _inventory = inventory_fixture(
        tmp_path,
        declared_output=tmp_path / "different-output",
    )
    jobs = tmp_path / "jobs"

    with pytest.raises(fast_analyze.AnalysisLaunchError, match="does not match"):
        fast_analyze.launch(command_args(analysis, jobs, cases=["R02"]))

    assert not jobs.exists()


@pytest.mark.parametrize("unsafe_target", ["analysis", "jobs", "jobs_symlink"])
def test_simulation_root_separation_is_enforced(
    fast_analyze, tmp_path, unsafe_target
):
    root = tmp_path / "campaign"
    analysis = tmp_path / "analysis"
    jobs = tmp_path / "jobs"
    if unsafe_target == "analysis":
        analysis = root / "runs/analysis"
    elif unsafe_target == "jobs":
        jobs = root / "runs/jobs"
    else:
        target = root / "runs/jobs"
        target.mkdir(parents=True)
        jobs.symlink_to(target, target_is_directory=True)
    inventory_fixture(tmp_path, analysis=analysis)

    with pytest.raises(fast_analyze.AnalysisLaunchError, match="outside runs"):
        fast_analyze.inventory_context(analysis, jobs)


def test_only_complete_cases_launch_unless_all_is_explicit(
    fast_analyze, tmp_path
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    complete_jobs = tmp_path / "complete-jobs"
    all_jobs = tmp_path / "all-jobs"

    fast_analyze.launch(command_args(analysis, complete_jobs))
    fast_analyze.launch(command_args(analysis, all_jobs, include_all=True))

    assert (complete_jobs / "R02/attempt-000/manifest.json").is_file()
    assert not (complete_jobs / "R03").exists()
    assert (all_jobs / "R02/attempt-000/manifest.json").is_file()
    assert (all_jobs / "R03/attempt-000/manifest.json").is_file()


def test_generated_batch_script_quotes_arguments_and_captures_exit(
    fast_analyze, tmp_path
):
    root, analysis, inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs with spaces"
    injected = tmp_path / "must-not-exist"
    hostile_option = f"2,4; touch {injected}"
    attempt = fast_analyze.prepare_attempt(
        analysis=analysis.resolve(),
        inventory_path=inventory.resolve(),
        jobs=jobs.resolve(),
        case_id="R02",
        report_options=["--alignment-shells", hostile_option],
        account="ast207",
        partition="batch",
        walltime="02:00:00",
        cpus_per_task=8,
        python=Path(sys.executable),
    )
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()
    fake_srun = fake_bin / "srun"
    fake_srun.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' \"$@\" > \"$FAKE_SRUN_ARGS\"\n"
        "exit \"$FAKE_SRUN_EXIT\"\n",
        encoding="utf-8",
    )
    fake_srun.chmod(0o755)
    captured = tmp_path / "srun-args.txt"
    environment = {
        "PATH": str(fake_bin),
        "FAKE_SRUN_ARGS": str(captured),
        "FAKE_SRUN_EXIT": "17",
    }

    completed = subprocess.run(
        ["/bin/bash", str(attempt / "run.sbatch")],
        text=True,
        capture_output=True,
        check=False,
        env=environment,
    )

    script = (attempt / "run.sbatch").read_text(encoding="utf-8")
    assert completed.returncode == 17
    assert (attempt / "exit_code.txt").read_text(encoding="utf-8") == "17\n"
    assert hostile_option in captured.read_text(encoding="utf-8").splitlines()
    assert "#SBATCH --nodes=1" in script
    assert "#SBATCH --ntasks=1" in script
    assert "cgl_lf_stage_i_fast_report.py" in script
    assert "analyze-case" in script
    assert not injected.exists()
    assert not (root / "runs").exists()


def test_retry_reuses_prepared_then_preserves_failed_attempt_options(
    fast_analyze, tmp_path
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    args = command_args(analysis, jobs, cases=["R02"], eddy_samples=9876)
    fast_analyze.launch(args)

    fast_analyze.retry(args)
    assert len(fast_analyze.attempt_directories(jobs.resolve(), "R02")) == 1

    first = jobs / "R02/attempt-000"
    (first / "exit_code.txt").write_text("9\n", encoding="utf-8")
    inherited = command_args(analysis, jobs, cases=["R02"])
    inherited.walltime = None
    inherited.cpus_per_task = None
    fast_analyze.retry(inherited)

    second = jobs / "R02/attempt-001"
    first_manifest = fast_analyze.load_json(first / "manifest.json")
    second_manifest = fast_analyze.load_json(second / "manifest.json")
    assert second_manifest["retry_of"] == str(first.resolve())
    assert second_manifest["report_options"] == first_manifest["report_options"]
    assert "--eddy-samples" in second_manifest["report_options"]
    assert "9876" in second_manifest["report_options"]
    assert second_manifest["walltime"] == first_manifest["walltime"]
    assert second_manifest["cpus_per_task"] == first_manifest["cpus_per_task"]

    (second / "exit_code.txt").write_text("9\n", encoding="utf-8")
    overridden = command_args(analysis, jobs, cases=["R02"])
    overridden.walltime = "05:30:00"
    overridden.cpus_per_task = 32
    fast_analyze.retry(overridden)
    third_manifest = fast_analyze.load_json(
        jobs / "R02/attempt-002/manifest.json"
    )
    assert third_manifest["report_options"] == first_manifest["report_options"]
    assert third_manifest["walltime"] == "05:30:00"
    assert third_manifest["cpus_per_task"] == 32


def test_retry_explicit_resources_update_unsubmitted_attempt(
    fast_analyze, tmp_path
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    args = command_args(analysis, jobs, cases=["R02"])
    fast_analyze.launch(args)
    args.walltime = "04:00:00"
    args.cpus_per_task = 24

    fast_analyze.retry(args)

    attempt = jobs / "R02/attempt-000"
    manifest = fast_analyze.load_json(attempt / "manifest.json")
    script = (attempt / "run.sbatch").read_text(encoding="utf-8")
    assert manifest["walltime"] == "04:00:00"
    assert manifest["cpus_per_task"] == 24
    assert "#SBATCH --time=04:00:00" in script
    assert "#SBATCH --cpus-per-task=24" in script


def test_zero_exit_requires_complete_bound_case_diagnostics(
    fast_analyze, tmp_path
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    fast_analyze.launch(command_args(analysis, jobs, cases=["R02"]))
    attempt = jobs / "R02/attempt-000"
    manifest = fast_analyze.load_json(attempt / "manifest.json")
    (attempt / "exit_code.txt").write_text("0\n", encoding="utf-8")

    assert fast_analyze.attempt_state(attempt, manifest) == (
        fast_analyze.FAILED_ANALYSIS_OUTPUT
    )

    attempt, _diagnostics = write_complete_analysis(
        fast_analyze, analysis, jobs
    )
    assert fast_analyze.attempt_state(attempt, manifest) == "COMPLETED"


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("analyzed_utc", "2000-01-01T00:00:00+00:00"),
        ("analysis_status", "partial"),
        ("snapshot_analysis_status", "failed"),
        ("analysis_errors", ["snapshot analysis failed"]),
    ],
)
def test_zero_exit_rejects_incomplete_analysis_fields(
    fast_analyze, tmp_path, field, value
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    fast_analyze.launch(command_args(analysis, jobs, cases=["R02"]))
    attempt, diagnostics = write_complete_analysis(fast_analyze, analysis, jobs)
    diagnostics[field] = value
    write_json(analysis / "cases/R02/diagnostics.json", diagnostics)
    manifest = fast_analyze.load_json(attempt / "manifest.json")

    assert fast_analyze.attempt_state(attempt, manifest) == (
        fast_analyze.FAILED_ANALYSIS_OUTPUT
    )


def test_zero_exit_rejects_stale_inventory_and_reporter_provenance(
    fast_analyze, tmp_path, monkeypatch
):
    reporter = tmp_path / "reporter.py"
    reporter.write_text("version = 1\n", encoding="utf-8")
    monkeypatch.setattr(fast_analyze, "REPORTER", reporter)
    _root, analysis, inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    fast_analyze.launch(command_args(analysis, jobs, cases=["R02"]))
    attempt, diagnostics = write_complete_analysis(fast_analyze, analysis, jobs)
    manifest = fast_analyze.load_json(attempt / "manifest.json")
    assert fast_analyze.attempt_state(attempt, manifest) == "COMPLETED"

    diagnostics["provenance"]["adapter"]["sha256"] = "0" * 64
    write_json(analysis / "cases/R02/diagnostics.json", diagnostics)
    assert fast_analyze.attempt_state(attempt, manifest) == (
        fast_analyze.FAILED_ANALYSIS_OUTPUT
    )

    write_complete_analysis(fast_analyze, analysis, jobs)
    original_inventory = fast_analyze.load_json(inventory)
    inventory_value = dict(original_inventory)
    inventory_value["assembled_utc"] = "changed"
    write_json(inventory, inventory_value)
    assert fast_analyze.attempt_state(attempt, manifest) == (
        fast_analyze.FAILED_ANALYSIS_OUTPUT
    )

    write_json(inventory, original_inventory)
    reporter.write_text("version = 2\n", encoding="utf-8")
    assert fast_analyze.attempt_state(attempt, manifest) == (
        fast_analyze.FAILED_ANALYSIS_OUTPUT
    )


def test_retry_parser_distinguishes_omitted_resource_overrides(
    fast_analyze, tmp_path
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)

    args = fast_analyze.parser().parse_args([str(analysis), "retry", "R02"])

    assert args.walltime is None
    assert args.cpus_per_task is None


@pytest.mark.parametrize(
    ("attribute", "value", "message"),
    [
        ("snapshot_workers", 0, "snapshot-workers"),
        ("snapshot_memory_budget_gib", 0.0, "snapshot-memory-budget-gib"),
        ("snapshot_memory_budget_gib", float("nan"), "snapshot-memory-budget-gib"),
    ],
)
def test_invalid_parallel_snapshot_limits_fail_before_job_writes(
    fast_analyze, tmp_path, attribute, value, message
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    args = command_args(analysis, jobs, cases=["R02"])
    setattr(args, attribute, value)

    with pytest.raises(fast_analyze.AnalysisLaunchError, match=message):
        fast_analyze.launch(args)

    assert not jobs.exists()


@pytest.mark.parametrize(
    ("state", "expected_attempts"),
    [
        ("RUNNING", 1),
        ("COMPLETED", 1),
        ("UNKNOWN", 1),
        ("TIMEOUT", 2),
        ("FAILED_ANALYSIS_OUTPUT", 2),
    ],
)
def test_retry_state_handling(
    fast_analyze, tmp_path, monkeypatch, state, expected_attempts
):
    _root, analysis, _inventory = inventory_fixture(tmp_path)
    jobs = tmp_path / "jobs"
    args = command_args(analysis, jobs, cases=["R02"])
    fast_analyze.launch(args)
    monkeypatch.setattr(fast_analyze, "attempt_state", lambda *_args: state)

    fast_analyze.retry(args)

    assert len(fast_analyze.attempt_directories(jobs.resolve(), "R02")) == (
        expected_attempts
    )
