"""Validated smoke workflows use passive physical energy without changing historical catalogs."""

import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace


def test_passive_smoke_executes_with_unchanged_catalog(tmp_path, monkeypatch):
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
        "paper_smoke_passive_alfvenic",
    ]
    assert manifest["disabled_cases"] == []
    assert len(workflow.workflow_cases("paper-mks24-stage-i")) == 16


def test_passive_smoke_uses_physical_energy_columns(tmp_path, monkeypatch):
    source = Path(__file__).resolve().parents[3] / "scripts/cgl_lf_workflow.py"
    spec = importlib.util.spec_from_file_location("cgl_passive_energy_workflow", source)
    workflow = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = workflow
    spec.loader.exec_module(workflow)
    user = {"force_prp2": [0., 1.], "force_prl2": [0., 0.],
            "volume": [1., 1.], "b2": [1., 1.], "b4": [1., 1.],
            "beta": [10., 10.], "mirror_vol": [0., 0.],
            "fire_vol": [0., 0.], "hard_vol": [0., 0.],
            "force_work": [0., .2]}
    mhd = {"cgl-J": [100., -100.], "cgl-A": [0., 0.],
           "1-KE": [0., .1], "2-KE": [0., .1], "3-KE": [0., 0.],
           "thermal-U": [3., 3.05]}
    monkeypatch.setattr(workflow, "case_history_path",
                        lambda case, root, key=None: Path("user" if key else "mhd"))
    monkeypatch.setattr(workflow, "parse_history",
                        lambda path: user if path.name == "user" else mhd)
    result = workflow.evaluate_paper_smoke(
        {"model_choices": {"passive_delta": "true",
                           "forcing_mode": "alfvenic_z_perpendicular"}}, tmp_path)
    assert result["passed"]
    assert result["energy_measure"] == "kinetic"
    assert result["energy_delta"] == .2
    assert abs(result["passive_thermal_energy_delta"] - .05) < 1.e-14
    assert not result["energy_work_residual_required"]
