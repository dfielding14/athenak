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
                "R02": {"case_id": "R02", "status": "complete"},
                "R03": {"case_id": "R03", "status": "in_progress"},
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
        pdf_bins=64,
        alignment_shells="2,4,6",
        eddy_samples=eddy_samples,
        eddy_bins=24,
        eddy_seed=731,
    )


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
    fast_analyze.retry(args)

    second = jobs / "R02/attempt-001"
    first_manifest = fast_analyze.load_json(first / "manifest.json")
    second_manifest = fast_analyze.load_json(second / "manifest.json")
    assert second_manifest["retry_of"] == str(first.resolve())
    assert second_manifest["report_options"] == first_manifest["report_options"]
    assert "--eddy-samples" in second_manifest["report_options"]
    assert "9876" in second_manifest["report_options"]


@pytest.mark.parametrize(
    ("state", "expected_attempts"),
    [("RUNNING", 1), ("COMPLETED", 1), ("UNKNOWN", 1), ("TIMEOUT", 2)],
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
