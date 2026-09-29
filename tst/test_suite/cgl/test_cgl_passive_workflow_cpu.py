"""Executable workflows skip fenced passive modes without changing historical catalogs."""

import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace


def test_passive_cases_are_recorded_but_not_executed(tmp_path, monkeypatch):
    source = Path(__file__).resolve().parents[3] / "scripts/cgl_lf_workflow.py"
    spec = importlib.util.spec_from_file_location("cgl_passive_workflow", source)
    workflow = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = workflow
    spec.loader.exec_module(workflow)
    monkeypatch.setattr(workflow, "executable_path", lambda _: (tmp_path, tmp_path / "athena"))
    monkeypatch.setattr(workflow, "ensure_executable", lambda *args: None)
    monkeypatch.setattr(workflow, "run_case", lambda case, *args: {"name": case.name})
    monkeypatch.setattr(workflow, "evaluate_manifest", lambda *args: {"passed": True})
    monkeypatch.setattr(workflow, "write_summary", lambda *args: None)
    paths = workflow.RunPaths(tmp_path)
    paths.create()
    result = workflow.execute_workflow(SimpleNamespace(workflow="paper-smoke"), paths)
    assert result == 0
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert [case["name"] for case in manifest["cases"]] == [
        "paper_smoke_active_alfvenic", "paper_smoke_active_random",
    ]
    assert manifest["disabled_cases"] == [{
        "name": "paper_smoke_passive_alfvenic",
        "reason": "Passive CGL thermal energy equation is disabled pending WO2.",
    }]
    assert len(workflow.workflow_cases("paper-mks24-stage-i")) == 16
