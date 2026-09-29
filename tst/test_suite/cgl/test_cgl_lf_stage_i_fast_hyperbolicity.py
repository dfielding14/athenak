"""Focused regressions for the direct-fast retained-state audit launcher."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
LAUNCHER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_hyperbolicity.py"


def load_launcher():
    name = "cgl_lf_stage_i_fast_hyperbolicity_tests"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def launcher():
    return load_launcher()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256(path),
    }


def add_snapshot(
    campaign: dict[str, object],
    case_id: str,
    time: float,
    *,
    complete: bool = True,
) -> None:
    output = Path(str(campaign["outputs"][case_id]))
    index_path = Path(str(campaign["analysis"])) / "cases" / case_id / "snapshots.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    number = len(index["snapshots"])
    records = []
    for rank in range(2):
        path = (
            output
            / "bin"
            / f"rank_{rank:08d}"
            / f"{case_id}.mhd_w_bcc.{number:05d}.bin"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"{case_id}:{time}:{rank}".encode())
        records.append({"path": str(path.resolve()), "size_bytes": path.stat().st_size})
    index["snapshots"].append(
        {
            "complete": complete,
            "expected_ranks": 2,
            "lineage_order": 0,
            "missing_rank_files": [],
            "empty_rank_files": [],
            "rank_files": records,
            "representative": records[0]["path"],
            "time": time,
        }
    )
    index["snapshot_count"] = len(index["snapshots"])
    index["complete_snapshot_count"] = sum(
        item["complete"] for item in index["snapshots"]
    )
    write_json(index_path, index)

    lineage_path = Path(str(campaign["analysis"])) / "cases" / case_id / "lineage.json"
    lineage = json.loads(lineage_path.read_text(encoding="utf-8"))
    lineage["snapshots"]["snapshot_count"] = index["snapshot_count"]
    lineage["snapshots"]["complete_snapshot_count"] = index["complete_snapshot_count"]
    write_json(lineage_path, lineage)


def refresh_inventory(campaign: dict[str, object]) -> Path:
    analysis = Path(str(campaign["analysis"]))
    inventory_path = analysis / "inventory.json"
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    inventory["cases"] = {
        case_id: json.loads(
            (analysis / "cases" / case_id / "lineage.json").read_text(encoding="utf-8")
        )
        for case_id in campaign["case_ids"]
    }
    write_json(inventory_path, inventory)
    return inventory_path


def add_case(
    campaign: dict[str, object],
    case_id: str,
    *,
    passive: bool,
    status: str = "in_progress",
    complete_snapshot: bool = True,
) -> None:
    root = Path(str(campaign["root"]))
    analysis = Path(str(campaign["analysis"]))
    matrix = Path(str(campaign["matrix"]))
    segment = root / "runs" / "selected" / case_id / "segment-000"
    output = segment / "output"
    manifest = segment / "manifest" / "fast_run.json"
    input_path = root / "inputs" / f"{case_id}.athinput"
    input_path.parent.mkdir(parents=True, exist_ok=True)
    input_path.write_text(
        f"<problem>\ncase={case_id}\npassive={str(passive).lower()}\n",
        encoding="utf-8",
    )
    write_json(
        manifest,
        {
            "schema_version": 1,
            "case_id": case_id,
            "output_dir": str(output.resolve()),
            "run_dir": str(segment.resolve()),
            "ranks": 2,
        },
    )
    index_path = analysis / "cases" / case_id / "snapshots.json"
    write_json(
        index_path,
        {
            "schema_version": 1,
            "snapshot_count": 0,
            "complete_snapshot_count": 0,
            "snapshots": [],
        },
    )
    campaign["outputs"][case_id] = output
    lineage = {
        "schema_version": 1,
        "case_id": case_id,
        "case_name": f"case_{case_id}",
        "status": status,
        "errors": [],
        "model_choices": {"passive_delta": str(passive).lower()},
        "input": binding(input_path),
        "lineage_identities": {
            "matrix_sha256": [sha256(matrix)],
            "input_sha256": [sha256(input_path)],
        },
        "lineage": [
            {
                "kind": "fast",
                "order": 0,
                "output": str(output.resolve()),
                "manifest": binding(manifest),
            }
        ],
        "snapshots": {
            "path": str(index_path.resolve()),
            "snapshot_count": 0,
            "complete_snapshot_count": 0,
        },
    }
    write_json(analysis / "cases" / case_id / "lineage.json", lineage)
    if complete_snapshot:
        add_snapshot(campaign, case_id, 0.5)
        add_snapshot(campaign, case_id, 0.75, complete=False)
    campaign["case_ids"].append(case_id)


@pytest.fixture
def campaign(tmp_path: Path) -> dict[str, object]:
    root = tmp_path / "campaign"
    analysis = tmp_path / "assembled"
    matrix = root / "matrix.json"
    adapter = root / "reporter.py"
    audit = tmp_path / "audit_cgl_hyperbolicity.py"
    matrix.parent.mkdir(parents=True)
    matrix.write_text('{"matrix": "authority"}\n', encoding="utf-8")
    adapter.write_text("# reporter authority\n", encoding="utf-8")
    audit.write_text("#!/usr/bin/env python3\n# test audit\n", encoding="utf-8")
    value: dict[str, object] = {
        "root": root,
        "analysis": analysis,
        "matrix": matrix,
        "adapter": adapter,
        "audit": audit,
        "case_ids": [],
        "outputs": {},
    }
    add_case(value, "R03", passive=False)
    add_case(value, "R06", passive=True)
    add_case(value, "R10", passive=False)
    add_case(value, "R11", passive=False, complete_snapshot=False)
    write_json(
        analysis / "inventory.json",
        {
            "schema_version": 1,
            "root": str(root.resolve()),
            "output": str(analysis.resolve()),
            "matrix": binding(matrix),
            "adapter": binding(adapter),
            "cases": {},
        },
    )
    refresh_inventory(value)
    return value


def command_args(
    campaign: dict[str, object],
    jobs: Path,
    *,
    cases: list[str] | None = None,
    exploratory: bool = False,
    include_partial: bool = True,
    snapshot_policy: str | None = None,
) -> SimpleNamespace:
    inventory = Path(str(campaign["analysis"])) / "inventory.json"
    args = SimpleNamespace(
        analysis=inventory,
        inventory_sha256=sha256(inventory),
        jobs_dir=jobs,
        audit_script=Path(str(campaign["audit"])),
        cases=cases or [],
        exploratory=exploratory,
        include_partial=include_partial,
        account="ast207",
        partition="batch",
        walltime="00:30:00",
        cpus_per_task=8,
        python=Path(sys.executable),
        submit=False,
    )
    if snapshot_policy is not None:
        args.snapshot_policy = snapshot_policy
    return args


def test_launch_authenticates_and_binds_final_retained_state(
    launcher, campaign, tmp_path
):
    jobs = tmp_path / "jobs"
    args = command_args(campaign, jobs, cases=["R03"])

    assert launcher.launch(args) == 0

    manifest = json.loads(
        (jobs / "R03/attempt-000/manifest.json").read_text(encoding="utf-8")
    )
    inventory = Path(str(campaign["analysis"])) / "inventory.json"
    assert manifest["inventory"]["sha256"] == sha256(inventory)
    assert manifest["matrix"]["sha256"] == sha256(Path(str(campaign["matrix"])))
    assert manifest["case_lineage"]["path"].endswith("/cases/R03/lineage.json")
    assert manifest["snapshot_index"]["path"].endswith("/cases/R03/snapshots.json")
    assert manifest["snapshot_policy"] == "latest"
    assert manifest["snapshot_coverage"]["selected_snapshot_count"] == 1
    assert manifest["snapshot_coverage"]["snapshot_index_complete_count"] == 1
    assert manifest["snapshot_coverage"]["all_complete_retained_snapshots_selected"] is (
        False
    )
    assert len(manifest["selected_snapshots"]) == 1
    assert manifest["selected_snapshot"]["time"] == 0.5
    assert manifest["selected_snapshot"]["expected_ranks"] == 2
    assert manifest["selected_snapshot"]["rank_inventory_sha256"] == (
        launcher.canonical_sha256(manifest["selected_snapshot"]["rank_files"])
    )
    assert manifest["audit_script"]["path"] == str(
        Path(str(campaign["audit"])).resolve()
    )
    assert manifest["command"][-3:] == ["--format", "json", "--hash-inputs"]
    assert len(manifest["command"]) == 6
    assert "rank_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]" in (
        manifest["command"][2]
    )
    assert Path(manifest["result"]["path"]).is_relative_to(jobs)
    assert not Path(manifest["result"]["path"]).is_relative_to(
        Path(str(campaign["root"])) / "runs"
    )
    batch = (jobs / "R03/attempt-000/run.sbatch").read_text(encoding="utf-8")
    assert "#SBATCH --nodes=1" in batch
    assert "#SBATCH --ntasks=1" in batch
    assert "audit_cgl_hyperbolicity.py" in batch
    assert "result.sha256" in batch


def test_literature_correct_formula_is_bound_and_result_mismatch_fails_closed(
    launcher, campaign, tmp_path
):
    jobs = tmp_path / "literature-jobs"
    args = command_args(campaign, jobs, cases=["R03"])
    args.formula = "literature-correct"

    assert launcher.launch(args) == 0

    attempt = jobs / "R03/attempt-000"
    manifest = json.loads((attempt / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["formula_id"] == "literature-correct"
    assert manifest["command"][-2:] == ["--formula", "literature-correct"]

    result = valid_result(launcher, manifest)
    result["provenance"]["formula_id"] = "literature-correct"
    write_json(attempt / "result.json", result)
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    (attempt / "exit_code.txt").write_text("0\n", encoding="utf-8")
    assert launcher.attempt_state(attempt, manifest) == "COMPLETED"

    result["provenance"]["formula_id"] = "qualified-legacy"
    write_json(attempt / "result.json", result)
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    assert launcher.attempt_state(attempt, manifest) == "INVALID_RESULT"


def test_inventory_and_lineage_provenance_fail_closed_without_job_writes(
    launcher, campaign, tmp_path
):
    inventory = Path(str(campaign["analysis"])) / "inventory.json"
    expected = sha256(inventory)
    inventory.write_text(inventory.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    jobs = tmp_path / "inventory-drift-jobs"
    args = command_args(campaign, jobs, cases=["R03"])
    args.inventory_sha256 = expected

    with pytest.raises(launcher.HyperbolicityLaunchError, match="SHA-256 differs"):
        launcher.launch(args)
    assert not jobs.exists()

    refresh_inventory(campaign)
    segment_manifest = (
        Path(str(campaign["root"]))
        / "runs/selected/R03/segment-000/manifest/fast_run.json"
    )
    segment_manifest.write_text(
        segment_manifest.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )
    jobs = tmp_path / "lineage-drift-jobs"
    args = command_args(campaign, jobs, cases=["R03"])
    with pytest.raises(launcher.HyperbolicityLaunchError, match="declared binding"):
        launcher.launch(args)
    assert not jobs.exists()


def test_scope_selects_active_partial_cases_and_gates_r10(
    launcher, campaign, tmp_path
):
    final_only = tmp_path / "final-only"
    assert (
        launcher.launch(command_args(campaign, final_only, include_partial=False))
        == 0
    )
    assert not final_only.exists()

    partial_jobs = tmp_path / "partial-jobs"
    assert launcher.launch(command_args(campaign, partial_jobs)) == 0
    assert (partial_jobs / "R03/attempt-000/manifest.json").is_file()
    assert not (partial_jobs / "R06").exists()
    assert not (partial_jobs / "R10").exists()
    assert not (partial_jobs / "R11").exists()

    marked_but_not_explicit = tmp_path / "marked-not-explicit"
    launcher.launch(command_args(campaign, marked_but_not_explicit, exploratory=True))
    assert not (marked_but_not_explicit / "R10").exists()

    explicit_not_marked = tmp_path / "explicit-not-marked"
    launcher.launch(command_args(campaign, explicit_not_marked, cases=["R10"]))
    assert not (explicit_not_marked / "R10").exists()

    exploratory_jobs = tmp_path / "exploratory-jobs"
    launcher.launch(
        command_args(campaign, exploratory_jobs, cases=["R10"], exploratory=True)
    )
    assert (exploratory_jobs / "R10/attempt-000/manifest.json").is_file()

def test_retry_reuses_prepared_then_retries_failed_and_stale_state(
    launcher, campaign, tmp_path
):
    jobs = tmp_path / "jobs"
    args = command_args(campaign, jobs, cases=["R03"])
    launcher.launch(args)

    launcher.retry(args)
    assert len(launcher.attempt_directories(jobs.resolve(), "R03")) == 1

    first = jobs / "R03/attempt-000"
    (first / "exit_code.txt").write_text("9\n", encoding="utf-8")
    launcher.retry(args)
    second = jobs / "R03/attempt-001"
    assert second.is_dir()
    first_manifest = json.loads(
        (first / "manifest.json").read_text(encoding="utf-8")
    )
    second_manifest = json.loads(
        (second / "manifest.json").read_text(encoding="utf-8")
    )
    assert second_manifest["retry_of"] == str(first.resolve())
    assert second_manifest["walltime"] == first_manifest["walltime"]
    assert second_manifest["cpus_per_task"] == first_manifest["cpus_per_task"]

    add_snapshot(campaign, "R03", 1.0)
    inventory = refresh_inventory(campaign)
    args.inventory_sha256 = sha256(inventory)
    launcher.retry(args)
    third = jobs / "R03/attempt-002"
    third_manifest = json.loads((third / "manifest.json").read_text(encoding="utf-8"))
    assert third_manifest["retry_of"] == str(second.resolve())
    assert third_manifest["selected_snapshot"]["time"] == 1.0
    assert third_manifest["selection_sha256"] != second_manifest["selection_sha256"]

    audit = Path(str(campaign["audit"]))
    audit.write_text(
        audit.read_text(encoding="utf-8") + "# revised\n", encoding="utf-8"
    )
    launcher.retry(args)
    fourth = jobs / "R03/attempt-003"
    fourth_manifest = json.loads((fourth / "manifest.json").read_text(encoding="utf-8"))
    assert fourth_manifest["retry_of"] == str(third.resolve())
    assert fourth_manifest["audit_script"]["sha256"] == sha256(audit)
    assert fourth_manifest["selection_sha256"] != third_manifest["selection_sha256"]


def valid_result(launcher, manifest: dict[str, object]) -> dict[str, object]:
    selected_snapshots = manifest.get(
        "selected_snapshots", [manifest["selected_snapshot"]]
    )
    input_patterns = []
    snapshots = []
    for snapshot_number, selected in enumerate(selected_snapshots):
        input_patterns.append(
            selected.get(
                "audit_input_pattern",
                [record["path"] for record in selected["rank_files"]],
            )
        )
        rank_files = [
            {
                **record,
                "rank_id": index,
                "sha256": f"{snapshot_number * 100 + index + 1:064x}",
            }
            for index, record in enumerate(selected["rank_files"])
        ]
        snapshots.append(
            {
                "time": selected["time"],
                "active_cgl_signal_speed": True,
                "rank_files": rank_files,
                "ranks_contiguous_from_zero": True,
                "input_inventory_sha256": launcher.canonical_sha256(rank_files),
                "aggregate": {
                    "negative": 0,
                    "nonfinite_discriminant": 0,
                },
            }
        )
    flattened_patterns = [
        item
        for value in input_patterns
        for item in (value if isinstance(value, list) else [value])
    ]
    return {
        "provenance": {
            "script_path": manifest["audit_script"]["path"],
            "script_sha256": manifest["audit_script"]["sha256"],
            "input_patterns": flattened_patterns,
            "hash_inputs": True,
        },
        "snapshots": snapshots,
    }


def test_all_snapshot_policy_binds_every_complete_snapshot_in_one_case_job(
    launcher, campaign, tmp_path
):
    add_snapshot(campaign, "R03", 1.0)
    inventory = refresh_inventory(campaign)
    jobs = tmp_path / "all-jobs"
    args = command_args(
        campaign, jobs, cases=["R03"], snapshot_policy="all"
    )
    args.inventory_sha256 = sha256(inventory)

    assert launcher.launch(args) == 0

    attempts = launcher.attempt_directories(jobs.resolve(), "R03")
    assert len(attempts) == 1
    manifest = json.loads(
        (attempts[0] / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["snapshot_policy"] == "all"
    assert manifest["snapshot_coverage"] == {
        "all_complete_retained_snapshots_selected": True,
        "selected_snapshot_count": 2,
        "selected_snapshot_positions": [0, 2],
        "selected_snapshot_times": [0.5, 1.0],
        "selected_snapshots_sha256": launcher.canonical_sha256(
            manifest["selected_snapshots"]
        ),
        "snapshot_index_complete_count": 2,
        "snapshot_policy": "all",
    }
    assert manifest["selected_snapshot"] == manifest["selected_snapshots"][-1]
    assert manifest["result"]["expected_snapshot_count"] == 2
    assert manifest["result"]["required_coverage"] == (
        "exactly_once_per_selected_snapshot"
    )
    assert manifest["command"][2:4] == [
        snapshot["audit_input_pattern"] for snapshot in manifest["selected_snapshots"]
    ]
    assert manifest["command"][-3:] == ["--format", "json", "--hash-inputs"]


def test_all_snapshot_result_requires_exactly_once_authenticated_coverage(
    launcher, campaign, tmp_path
):
    add_snapshot(campaign, "R03", 1.0)
    inventory = refresh_inventory(campaign)
    jobs = tmp_path / "all-result-jobs"
    args = command_args(
        campaign, jobs, cases=["R03"], snapshot_policy="all"
    )
    args.inventory_sha256 = sha256(inventory)
    launcher.launch(args)
    attempt = jobs / "R03/attempt-000"
    manifest = json.loads((attempt / "manifest.json").read_text(encoding="utf-8"))
    result = valid_result(launcher, manifest)
    result["snapshots"].reverse()
    write_json(attempt / "result.json", result)
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    (attempt / "exit_code.txt").write_text("0\n", encoding="utf-8")

    assert launcher.attempt_state(attempt, manifest) == "COMPLETED"
    assert launcher.audit_disposition(attempt, manifest) == "HYPERBOLIC"

    tampered = json.loads(json.dumps(manifest))
    tampered["command"][2] += ".unexpected"
    assert launcher.attempt_state(attempt, tampered) == "INVALID_RESULT"
    tampered = json.loads(json.dumps(manifest))
    tampered["selected_snapshot"]["time"] = -1.0
    assert launcher.attempt_state(attempt, tampered) == "INVALID_RESULT"

    result["snapshots"][-1]["aggregate"]["negative"] = 1
    write_json(attempt / "result.json", result)
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    assert launcher.attempt_state(attempt, manifest) == "COMPLETED"
    assert launcher.audit_disposition(attempt, manifest) == "NEGATIVE"

    result["snapshots"] = result["snapshots"][:1]
    write_json(attempt / "result.json", result)
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    assert launcher.attempt_state(attempt, manifest) == "INVALID_RESULT"


def test_retry_treats_policy_switch_and_retained_index_growth_as_stale(
    launcher, campaign, tmp_path
):
    jobs = tmp_path / "policy-retry-jobs"
    latest = command_args(campaign, jobs, cases=["R03"])
    launcher.launch(latest)

    all_snapshots = command_args(
        campaign, jobs, cases=["R03"], snapshot_policy="all"
    )
    launcher.retry(all_snapshots)
    switched = jobs / "R03/attempt-001"
    switched_manifest = json.loads(
        (switched / "manifest.json").read_text(encoding="utf-8")
    )
    assert switched_manifest["snapshot_policy"] == "all"
    assert switched_manifest["retry_of"].endswith("R03/attempt-000")

    add_snapshot(campaign, "R03", 1.0)
    inventory = refresh_inventory(campaign)
    all_snapshots.inventory_sha256 = sha256(inventory)
    launcher.retry(all_snapshots)
    grown = jobs / "R03/attempt-002"
    grown_manifest = json.loads((grown / "manifest.json").read_text(encoding="utf-8"))
    assert grown_manifest["snapshot_coverage"]["selected_snapshot_count"] == 2
    assert grown_manifest["retry_of"] == str(switched.resolve())

    add_snapshot(campaign, "R03", 1.25, complete=False)
    inventory = refresh_inventory(campaign)
    all_snapshots.inventory_sha256 = sha256(inventory)
    launcher.retry(all_snapshots)
    rebound = jobs / "R03/attempt-003"
    rebound_manifest = json.loads(
        (rebound / "manifest.json").read_text(encoding="utf-8")
    )
    assert rebound_manifest["snapshot_coverage"]["selected_snapshot_count"] == 2
    assert rebound_manifest["retry_of"] == str(grown.resolve())
    assert rebound_manifest["selection_sha256"] != grown_manifest["selection_sha256"]


def test_completed_result_is_provenance_checked_and_status_reports_it(
    launcher, campaign, tmp_path, capsys
):
    jobs = tmp_path / "jobs"
    args = command_args(campaign, jobs, cases=["R03"])
    launcher.launch(args)
    attempt = jobs / "R03/attempt-000"
    manifest = json.loads((attempt / "manifest.json").read_text(encoding="utf-8"))
    write_json(attempt / "result.json", valid_result(launcher, manifest))
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    (attempt / "exit_code.txt").write_text("0\n", encoding="utf-8")

    assert launcher.attempt_state(attempt, manifest) == "COMPLETED"
    status_args = SimpleNamespace(
        analysis=args.analysis,
        inventory_sha256=args.inventory_sha256,
        jobs_dir=args.jobs_dir,
        audit_script=args.audit_script,
        cases=["R03"],
        all_attempts=False,
    )
    launcher.status(status_args)
    output = capsys.readouterr().out
    assert "COMPLETED" in output
    assert "HYPERBOLIC" in output

    result = json.loads((attempt / "result.json").read_text(encoding="utf-8"))
    result["provenance"]["script_sha256"] = "0" * 64
    write_json(attempt / "result.json", result)
    (attempt / "result.sha256").write_text(
        f"{sha256(attempt / 'result.json')}  result.json\n",
        encoding="utf-8",
    )
    assert launcher.attempt_state(attempt, manifest) == "INVALID_RESULT"
