"""Focused tests for corrected-production downstream orchestration."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
TOOL = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_corrected_downstream.py"
IDENTITY_TOOL = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_corrected_identity.py"
REVISION = "0c406312fa35d5e1c7041d80333b0ea24b0127ae"
LEGACY_EXECUTABLE_SHA256 = (
    "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def downstream():
    return load_module("cgl_corrected_downstream_tests", TOOL)


@pytest.fixture(scope="module")
def identity_tool():
    return load_module("cgl_corrected_identity_for_downstream_tests", IDENTITY_TOOL)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_artifact(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def campaign_fixture(
    downstream,
    identity_tool,
    tmp_path: Path,
    *,
    include_passive: bool = True,
    legacy_active: str | None = None,
) -> dict[str, object]:
    authority = tmp_path / "authority"
    artifacts = {
        "audit": write_artifact(authority / "audit.py", "audit\n"),
        "eos": write_artifact(authority / "eos.hpp", "corrected ppar squared\n"),
        "executable": write_artifact(authority / "athena", "corrected executable\n"),
        "matrix": write_artifact(authority / "matrix.json", '{"cases": []}\n'),
    }
    root = tmp_path / "campaign"
    run_relatives = [
        f"runs/mks24-stage-i-fast/E03-forcing-policy/{case_id}"
        for case_id in downstream.ACTIVE_CASES
    ]
    identity_path = identity_tool.write_identity(
        "mks24-stage-i-corrected-test",
        root,
        REVISION,
        artifacts,
        run_relatives,
    )
    identity_sha = sha256(identity_path)
    corrected_executable_sha = sha256(artifacts["executable"])
    matrix_sha = sha256(artifacts["matrix"])
    cases: dict[str, object] = {}

    for case_id, relative in zip(downstream.ACTIVE_CASES, run_relatives):
        case_root = root / relative
        segment = case_root / "fast_s000_t0_to_t10"
        manifest_path = segment / "manifest/fast_run.json"
        executable_sha = (
            LEGACY_EXECUTABLE_SHA256 if case_id == legacy_active
            else corrected_executable_sha
        )
        executable_path = (
            str(tmp_path / "legacy/athena")
            if case_id == legacy_active
            else str(artifacts["executable"].resolve())
        )
        manifest = {
            "campaign_id": "mks24-stage-i-corrected-test",
            "campaign_root": str(root.resolve()),
            "root": str(root.resolve()),
            "source_revision": REVISION,
            "campaign_identity": str(identity_path),
            "campaign_identity_sha256": identity_sha,
            "matrix_sha256": matrix_sha,
            "executable": executable_path,
            "executable_sha256": executable_sha,
            "legacy_restart_permitted": False,
            "case_id": case_id,
            "sequence": 0,
            "start_time": 0.0,
            "restart": None,
            "launch_origin": "fresh_t0_corrected_production",
            "fresh_lineage_root": True,
        }
        write_json(manifest_path, manifest)
        cases[case_id] = {
            "case_id": case_id,
            "case_name": f"case-{case_id}",
            "status": "complete",
            "errors": [],
            "model_choices": {"passive_delta": "false"},
            "lineage_identities": {
                "matrix_sha256": [matrix_sha],
                "executable_sha256": [executable_sha],
            },
            "lineage": [{
                "kind": "fast",
                "order": 0,
                "segment_dir": str(segment.resolve()),
                "manifest": binding(manifest_path),
            }],
        }

    selected = list(downstream.ACTIVE_CASES)
    if include_passive:
        legacy = write_artifact(tmp_path / "legacy/athena", "legacy executable\n")
        legacy_sha = sha256(legacy)
        for case_id in downstream.PASSIVE_CASES:
            cases[case_id] = {
                "case_id": case_id,
                "case_name": f"case-{case_id}",
                "status": "complete",
                "errors": [],
                "model_choices": {"passive_delta": "true"},
                "lineage_identities": {
                    "matrix_sha256": [matrix_sha],
                    "executable_sha256": [legacy_sha],
                },
                "lineage": [{"kind": "legacy-passive"}],
            }
        selected.extend(downstream.PASSIVE_CASES)

    output = tmp_path / "analysis/corrected-composite"
    for case_id, case in cases.items():
        write_json(output / "cases" / case_id / "lineage.json", case)
    reporter = write_artifact(tmp_path / "authority/reporter.py", "reporter\n")
    inventory = {
        "schema_version": 1,
        "root": str(root.resolve()),
        "output": str(output.resolve()),
        "matrix": binding(artifacts["matrix"]),
        "adapter": binding(reporter),
        "cases": cases,
    }
    if include_passive:
        inventory.update({
            "record_type": "cgl_lf_stage_i_corrected_composite_report",
            "composite_adapter": binding(downstream.CORRECTED_REPORT_TOOL),
        })
    inventory_path = output / "inventory.json"
    write_json(inventory_path, inventory)
    return {
        "root": root,
        "identity": identity_path,
        "inventory": inventory_path,
        "inventory_sha256": sha256(inventory_path),
        "artifacts": artifacts,
        "selected": sorted(selected),
    }


def install_composite_validator(downstream, monkeypatch) -> None:
    original = downstream.load_module

    class CompositeError(RuntimeError):
        pass

    fake = SimpleNamespace(
        DEFAULT_CONFIG=object(),
        CompositeReportError=CompositeError,
        validate_composite_inventory=lambda _config, output: {
            "result": "pass",
            "output": str(output),
            "cases": list(downstream.ALL_CASES),
        },
    )
    monkeypatch.setattr(
        downstream,
        "load_module",
        lambda name, path: fake
        if Path(path) == downstream.CORRECTED_REPORT_TOOL
        else original(name, path),
    )


def test_composite_inventory_allows_legacy_passive_but_binds_every_active_corrected(
    downstream, identity_tool, tmp_path, monkeypatch
) -> None:
    install_composite_validator(downstream, monkeypatch)
    fixture = campaign_fixture(downstream, identity_tool, tmp_path)

    context = downstream.validate_inventory(
        fixture["identity"], fixture["inventory"], fixture["inventory_sha256"]
    )

    assert context["selected_cases"] == fixture["selected"]
    assert context["active_cases"] == list(downstream.ACTIVE_CASES)
    assert context["passive_cases"] == list(downstream.PASSIVE_CASES)
    assert set(context["active_manifest_bindings"]) == set(downstream.ACTIVE_CASES)
    corrected_sha = sha256(fixture["artifacts"]["executable"])
    assert all(
        corrected_sha not in hashes
        for hashes in context["passive_executable_sha256"].values()
    )


def test_isolated_corrected_active_inventory_is_supported(
    downstream, identity_tool, tmp_path
) -> None:
    fixture = campaign_fixture(
        downstream, identity_tool, tmp_path, include_passive=False
    )

    context = downstream.validate_inventory(
        fixture["identity"], fixture["inventory"], fixture["inventory_sha256"]
    )

    assert context["selected_cases"] == sorted(downstream.ACTIVE_CASES)
    assert context["passive_cases"] == []


def test_legacy_active_evidence_is_rejected(
    downstream, identity_tool, tmp_path
) -> None:
    fixture = campaign_fixture(
        downstream,
        identity_tool,
        tmp_path,
        include_passive=False,
        legacy_active="R03",
    )

    with pytest.raises(
        downstream.CorrectedDownstreamError,
        match="R03 executable lineage is not corrected",
    ):
        downstream.validate_inventory(
            fixture["identity"], fixture["inventory"], fixture["inventory_sha256"]
        )


def test_commands_are_all_snapshot_literature_correct_and_never_mix_passive_hyper(
    downstream, identity_tool, tmp_path, monkeypatch
) -> None:
    install_composite_validator(downstream, monkeypatch)
    fixture = campaign_fixture(downstream, identity_tool, tmp_path)
    context = downstream.validate_inventory(
        fixture["identity"], fixture["inventory"], fixture["inventory_sha256"]
    )
    commands = downstream.workflow_commands(
        context,
        tmp_path / "workflow",
        Path(sys.executable),
        "ast207",
        "batch",
        "02:00:00",
        56,
    )
    hyper = commands["hyperbolicity"]["command"]
    analysis = commands["analysis"]["command"]
    ct = commands["ct"]["command"]
    publication = commands["publication"]["command"]

    assert "--formula" in hyper
    assert hyper[hyper.index("--formula") + 1] == "literature-correct"
    assert hyper[hyper.index("--snapshot-policy") + 1] == "all"
    assert "--exploratory" in hyper
    assert not set(downstream.PASSIVE_CASES) & set(hyper)
    assert set(downstream.ACTIVE_CASES) <= set(hyper)
    assert set(downstream.ALL_CASES) <= set(analysis)
    assert ct[ct.index("--snapshot-policy") + 1] == "all"
    assert ",".join(sorted(downstream.ALL_CASES)) in ct
    assert "--acceptance" in publication
    assert "--submit" not in hyper
    assert "--submit" not in analysis


def test_prepare_writes_jobs_and_workflow_without_submitting(
    downstream, identity_tool, tmp_path, monkeypatch
) -> None:
    install_composite_validator(downstream, monkeypatch)
    fixture = campaign_fixture(downstream, identity_tool, tmp_path)
    context = downstream.validate_inventory(
        fixture["identity"], fixture["inventory"], fixture["inventory_sha256"]
    )
    workflow_root = tmp_path / "workflow"
    observed: list[list[str]] = []
    tools = {"orchestrator": downstream.artifact_binding(TOOL)}
    attempts = {
        "hyperbolicity": {
            case_id: str(tmp_path / f"hyper/{case_id}/attempt-000")
            for case_id in downstream.ACTIVE_CASES
        },
        "analysis": {
            case_id: str(tmp_path / f"analysis/{case_id}/attempt-000")
            for case_id in context["selected_cases"]
        },
    }
    monkeypatch.setattr(downstream, "validate_inventory", lambda *_args: context)
    monkeypatch.setattr(downstream, "tool_bindings", lambda _context: tools)
    monkeypatch.setattr(downstream, "run_checked", lambda command: observed.append(command))
    monkeypatch.setattr(
        downstream, "prepared_attempts", lambda _commands, _context: attempts
    )
    monkeypatch.setattr(downstream, "patch_hyper_attempts", lambda *_args: None)
    args = SimpleNamespace(
        identity=fixture["identity"],
        inventory=fixture["inventory"],
        inventory_sha256=fixture["inventory_sha256"],
        workflow_root=workflow_root,
        python=Path(sys.executable),
        account="ast207",
        partition="batch",
        walltime="02:00:00",
        cpus_per_task=56,
    )

    path = downstream.prepare_workflow(args)

    workflow = json.loads(path.read_text(encoding="utf-8"))
    assert len(observed) == 2
    assert all("--submit" not in command for command in observed)
    assert workflow["formula_binding"] == downstream.exact_formula_binding(context)
    assert workflow["formula_binding"]["audit_script"] == binding(
        fixture["artifacts"]["audit"]
    )
    assert workflow["formula_binding"]["corrected_eos"] == binding(
        fixture["artifacts"]["eos"]
    )
    for stage in ("ct", "publication"):
        job = workflow_root / f"jobs/{stage}/attempt-000"
        manifest = json.loads((job / "manifest.json").read_text(encoding="utf-8"))
        script = (job / "run.sbatch").read_text(encoding="utf-8")
        assert manifest["job_id"] is None
        assert manifest["formula_id"] == "literature-correct"
        assert "require_sha" in script
        assert "/usr/bin/sbatch" not in script


def test_formula_binding_drift_is_rejected(
    downstream, identity_tool, tmp_path, monkeypatch
) -> None:
    fixture = campaign_fixture(
        downstream, identity_tool, tmp_path, include_passive=False
    )
    context = downstream.validate_inventory(
        fixture["identity"], fixture["inventory"], fixture["inventory_sha256"]
    )
    workflow = {
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "corrected_executable": context["identity_artifacts"]["executable"],
        "formula_binding": {
            "audit_formula_id": "qualified-legacy",
            "executable_formula_id": "literature-correct",
            "compatibility": "incompatible",
        },
        "tools": {},
        "case_classification": {
            "corrected_active": context["active_cases"],
            "reused_passive": context["passive_cases"],
            "selected": context["selected_cases"],
        },
        "active_execution_manifests": context["active_manifest_bindings"],
    }
    monkeypatch.setattr(downstream, "validate_inventory", lambda *_args: context)

    with pytest.raises(
        downstream.CorrectedDownstreamError,
        match="formula/executable binding differs",
    ):
        downstream.workflow_context(workflow)


def test_complete_writes_immutable_corrected_only_record_and_pointer(
    downstream, tmp_path, monkeypatch
) -> None:
    root = tmp_path / "workflow"
    root.mkdir()
    inventory = write_artifact(tmp_path / "inventory.json", "{}\n")
    identity = write_artifact(tmp_path / "campaign-identity.json", "{}\n")
    executable = write_artifact(tmp_path / "athena", "corrected\n")
    workflow_file = write_artifact(root / "workflow.json", "{}\n")
    workflow = {
        "workflow_root": str(root.resolve()),
        "case_classification": {
            "corrected_active": list(downstream.ACTIVE_CASES),
            "reused_passive": list(downstream.PASSIVE_CASES),
            "selected": list(downstream.ALL_CASES),
        },
    }
    context = {
        "identity_binding": binding(identity),
        "inventory_binding": binding(inventory),
        "active_cases": list(downstream.ACTIVE_CASES),
        "identity_artifacts": {"executable": binding(executable)},
    }
    record = {
        "schema": downstream.COMPLETION_SCHEMA,
        "status": "complete",
        "corrected-only": True,
        "corrected_executable": binding(executable),
        "formula_binding": {
            "audit_formula_id": "literature-correct",
            "executable_formula_id": "literature-correct",
            "compatibility": "compatible",
            "active_cases": list(downstream.ACTIVE_CASES),
        },
    }
    monkeypatch.setattr(
        downstream,
        "load_workflow",
        lambda _root: (workflow, binding(workflow_file)),
    )
    monkeypatch.setattr(downstream, "workflow_context", lambda _workflow: context)
    monkeypatch.setattr(
        downstream,
        "completion_record",
        lambda _workflow, _binding, _context: record,
    )

    pointer_path = downstream.complete_workflow(root)
    pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
    retained = Path(pointer["completion_record"]["path"])

    assert pointer_path.name == "corrected-production-complete.json"
    assert retained.name == f"{sha256(inventory)}.json"
    assert pointer["campaign_kind"] == "corrected-production"
    assert pointer["formula_binding"]["audit_formula_id"] == "literature-correct"
    assert pointer["formula_binding"]["executable_formula_id"] == "literature-correct"
    assert json.loads(retained.read_text(encoding="utf-8")) == record
    assert downstream.complete_workflow(root) == pointer_path
