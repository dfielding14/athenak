"""Focused tests for the corrected/composite manuscript-ready data gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
TOOL = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_corrected_release.py"


def load_tool():
    name = "cgl_lf_stage_i_fast_corrected_release_tests"
    spec = importlib.util.spec_from_file_location(name, TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def gate():
    return load_tool()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    path = path.resolve()
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def write_json(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return path


def write_text(path: Path, value: str = "fixture\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8")
    return path


def snapshot_files(root: Path) -> dict[str, bytes]:
    return {
        str(path.relative_to(root)): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file()
    }


def exact_invocation(gate, analysis: Path, publication: Path, roots: list[Path]) -> list[str]:
    invocation = [
        sys.executable,
        str(gate.PUBLICATION_TOOL.resolve()),
        str(analysis.resolve()),
        "--output",
        str(publication.resolve()),
    ]
    for root in sorted(path.resolve() for path in roots):
        invocation.extend(["--acceptance", str(root)])
    return invocation


def fixture_tree(
    gate, tmp_path: Path, *, science_result: str = "pass", ct_result: str = "pass"
) -> dict[str, object]:
    identity_path = write_json(tmp_path / "identity.json", {"identity": "corrected"})
    analysis = tmp_path / "analysis/composite"
    cases: dict[str, object] = {}
    for case_id in gate.ALL_CASES:
        evidence_class = (
            gate.ACTIVE_EVIDENCE_CLASS
            if case_id in gate.ACTIVE_CASES
            else gate.PASSIVE_EVIDENCE_CLASS
        )
        cases[case_id] = {
            "case_id": case_id,
            "status": "complete",
            "evidence_class": evidence_class,
            "execution_authority": {"evidence_class": evidence_class},
        }
    inventory = {
        "record_type": "cgl_lf_stage_i_corrected_composite_report",
        "output": str(analysis.resolve()),
        "evidence_classes": {
            gate.ACTIVE_EVIDENCE_CLASS: list(gate.ACTIVE_CASES),
            gate.PASSIVE_EVIDENCE_CLASS: list(gate.PASSIVE_CASES),
        },
        "cases": cases,
    }
    inventory_path = write_json(analysis / "inventory.json", inventory)
    context = {
        "inventory_kind": "corrected-composite",
        "active_cases": list(gate.ACTIVE_CASES),
        "passive_cases": list(gate.PASSIVE_CASES),
        "selected_cases": list(gate.ALL_CASES),
        "inventory": inventory,
        "inventory_binding": binding(inventory_path),
        "identity_binding": binding(identity_path),
        "inventory_output": analysis.resolve(),
    }

    acceptance = tmp_path / "analysis/composite-acceptance"
    acceptance_provenance = write_json(
        acceptance / "provenance.json", {"record_type": "acceptance-provenance"}
    )
    acceptance_summary = write_json(
        acceptance / "summary.json",
        {
            "record_type": "cgl-lf-stage-i-direct-fast-acceptance-summary",
            "selected_cases": list(gate.ALL_CASES),
        },
    )
    acceptance_campaign = write_json(
        acceptance / "campaign_evidence.json",
        {
            "record_type": "cgl-lf-stage-i-direct-fast-campaign-evidence",
            "result": "inconclusive",
        },
    )

    science_root = tmp_path / "analysis/composite-science"
    direct_science = write_json(
        science_root / "science.json",
        {
            "record_type": gate.DIRECT_SCIENCE_RECORD_TYPE,
            "result": science_result,
            "selected_cases": list(gate.ALL_CASES),
        },
    )
    science_provenance = write_json(
        science_root / "provenance.json",
        {"record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-provenance"},
    )
    science_table = write_text(science_root / "tables/cases.csv")
    science_record = {
        "schema": gate.SCIENCE_SCHEMA,
        "schema_version": 1,
        "record_type": gate.SCIENCE_RECORD_TYPE,
        "status": "complete",
        "result": science_result,
        "campaign_kind": "corrected-composite",
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "case_classification": gate.expected_science_classification(),
        "acceptance": {
            "root": str(acceptance.resolve()),
            "provenance": binding(acceptance_provenance),
            "summary": binding(acceptance_summary),
            "campaign_evidence": binding(acceptance_campaign),
        },
        "science": {
            "root": str(science_root.resolve()),
            "provenance": binding(science_provenance),
            "record": binding(direct_science),
            "tables": {"tables/cases.csv": binding(science_table)},
        },
    }
    corrected_science = write_json(
        science_root / "corrected-composite-science.json", science_record
    )

    workflow = tmp_path / "workflow"
    ct_record = {
        "record_type": gate.CT_RECORD_TYPE,
        "result": ct_result,
        "selection": {
            "cases": list(gate.ALL_CASES),
            "snapshot_policy": "all",
        },
        "cases": {
            case_id: {
                "provenance_authenticated": True,
                "native_restart_ct": {"coverage_complete": True},
            }
            for case_id in gate.ALL_CASES
        },
    }
    ct_path = write_json(workflow / "ct/ct_audit.json", ct_record)

    hyper_root = analysis / "corrected-downstream/hyperbolicity"
    hyperbolicity: dict[str, object] = {}
    for case_id in gate.ACTIVE_CASES:
        manifest = write_json(
            hyper_root / case_id / "attempt-000/manifest.json",
            {"record_type": "hyperbolicity-manifest", "case_id": case_id},
        )
        result = write_json(
            hyper_root / case_id / "attempt-000/result.json",
            {"record_type": "hyperbolicity-result", "case_id": case_id},
        )
        hyperbolicity[case_id] = {
            "manifest": binding(manifest),
            "result": binding(result),
        }

    downstream_analysis: dict[str, object] = {}
    for case_id in gate.ALL_CASES:
        diagnostics = write_json(
            analysis / "cases" / case_id / "diagnostics.json",
            {"record_type": "case-diagnostics", "case_id": case_id},
        )
        downstream_analysis[case_id] = {"diagnostics": binding(diagnostics)}

    figure = write_text(analysis / "figures/paper/overview.pdf")
    write_json(
        analysis / "manuscript/results.json",
        {
            "case_results": {case_id: {} for case_id in gate.ALL_CASES},
            "campaign_health": {
                "partial_or_unavailable_cases": [],
                "structural_error_cases": [],
            },
        },
    )
    write_text(analysis / "manuscript/results_macros.tex")
    write_json(analysis / "manuscript/claim_evidence.json", {"claims": []})
    write_text(analysis / "manuscript/claim_evidence.md")
    write_json(
        analysis / "manuscript/figure_manifest.json",
        {"figures": [str(figure.relative_to(analysis))], "warnings": []},
    )
    write_json(
        analysis / "manuscript/provenance_manifest.json",
        {
            "inputs": [context["inventory_binding"]],
            "case_inputs": {case_id: {} for case_id in gate.ALL_CASES},
        },
    )
    write_text(analysis / "manuscript/open_questions.md")
    write_json(
        analysis / "verify.json",
        {
            "result": "pass",
            "require_complete": True,
            "cases": {case_id: {"errors": []} for case_id in gate.ALL_CASES},
            "errors": [],
            "warnings": [],
            "inventory": context["inventory_binding"],
        },
    )

    publication = workflow / "publication"
    products = [
        binding(write_text(publication / relative, f"{relative}\n"))
        for relative in gate.REQUIRED_PUBLICATION_PRODUCTS
    ]
    sources = [
        context["inventory_binding"],
        binding(ct_path),
        binding(corrected_science),
        binding(direct_science),
        binding(acceptance_provenance),
        binding(acceptance_summary),
        binding(acceptance_campaign),
        binding(science_provenance),
        *[
            record[key]
            for record in hyperbolicity.values()
            for key in ("manifest", "result")
        ],
        *[record["diagnostics"] for record in downstream_analysis.values()],
    ]
    roots = [acceptance, science_root, workflow / "ct", hyper_root]
    publication_manifest = {
        "schema_version": 2,
        "record_type": gate.PUBLICATION_RECORD_TYPE,
        "analysis_output": str(analysis.resolve()),
        "evidence_state": (
            "complete integration / partial evidence"
            if science_result == "inconclusive" or ct_result == "inconclusive"
            else "complete"
        ),
        "normalized_invocation": exact_invocation(gate, analysis, publication, roots),
        "renderer_ingestion_warnings": [],
        "reviewed_science": {
            "record_type": gate.DIRECT_SCIENCE_RECORD_TYPE,
            "result": science_result,
            "selected_cases": list(gate.ALL_CASES),
            "release_authorizing": False,
        },
        "direct_ct_audit": {
            "record_type": gate.CT_RECORD_TYPE,
            "numerical_result": ct_result,
            "selected_cases": list(gate.ALL_CASES),
            "full_stage_i_coverage": True,
            "release_authorizing": False,
        },
        "sources": sources,
        "products": products,
    }
    publication_manifest_path = write_json(
        publication / "manifest.json", publication_manifest
    )

    completion = {
        "schema": gate.DOWNSTREAM_COMPLETION_SCHEMA,
        "schema_version": 1,
        "status": "complete",
        "campaign_kind": "corrected-production",
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "case_classification": gate.expected_downstream_classification(),
        "analysis": downstream_analysis,
        "hyperbolicity": hyperbolicity,
        "ct": {"audit": binding(ct_path)},
        "publication": {"manifest": binding(publication_manifest_path)},
    }
    completion_path = write_json(
        workflow / "completion/records/inventory.json", completion
    )
    pointer_path = write_json(
        workflow / "completion/corrected-production-complete.json",
        {
            "schema": gate.DOWNSTREAM_POINTER_SCHEMA,
            "schema_version": 1,
            "status": "complete",
            "campaign_kind": "corrected-production",
            "campaign_identity": context["identity_binding"],
            "inventory": context["inventory_binding"],
            "completion_record": binding(completion_path),
        },
    )
    criteria = write_json(tmp_path / "criteria.json", {"criteria": True})
    criteria_review = write_json(tmp_path / "criteria-review.json", {"review": True})
    args = SimpleNamespace(
        identity=identity_path,
        inventory=inventory_path,
        inventory_sha256=sha256(inventory_path),
        workflow_root=workflow,
        acceptance_output=acceptance,
        science_output=science_root,
        marker=workflow / gate.RECORD_NAME,
        criteria=criteria,
        criteria_review=criteria_review,
    )
    return {
        "args": args,
        "context": context,
        "corrected_science": corrected_science,
        "publication_manifest": publication_manifest_path,
        "completion": completion_path,
        "pointer": pointer_path,
        "hyper_root": hyper_root,
        "ct": ct_path,
        "root": tmp_path,
    }


def install_dependencies(gate, monkeypatch, fixture: dict[str, object]) -> None:
    science = SimpleNamespace(
        validate_corrected_context=lambda *_args: fixture["context"],
        validate_products=lambda *_args: json.loads(
            fixture["corrected_science"].read_text(encoding="utf-8")
        ),
    )
    downstream = SimpleNamespace(validate_completion=lambda *_args: fixture["pointer"])
    monkeypatch.setattr(
        gate,
        "load_module",
        lambda _name, path: (
            science
            if Path(path) == gate.CORRECTED_SCIENCE_TOOL
            else downstream
            if Path(path) == gate.CORRECTED_DOWNSTREAM_TOOL
            else pytest.fail(f"unexpected dependency: {path}")
        ),
    )


def rebind_completion(fixture: dict[str, object]) -> None:
    completion = json.loads(fixture["completion"].read_text(encoding="utf-8"))
    completion["publication"]["manifest"] = binding(fixture["publication_manifest"])
    write_json(fixture["completion"], completion)
    pointer = json.loads(fixture["pointer"].read_text(encoding="utf-8"))
    pointer["completion_record"] = binding(fixture["completion"])
    write_json(fixture["pointer"], pointer)


def rebind_ct_evidence(fixture: dict[str, object]) -> None:
    ct_binding = binding(fixture["ct"])
    manifest = json.loads(fixture["publication_manifest"].read_text(encoding="utf-8"))
    manifest["sources"] = [
        ct_binding if value["path"] == ct_binding["path"] else value
        for value in manifest["sources"]
    ]
    write_json(fixture["publication_manifest"], manifest)
    completion = json.loads(fixture["completion"].read_text(encoding="utf-8"))
    completion["ct"]["audit"] = ct_binding
    completion["publication"]["manifest"] = binding(fixture["publication_manifest"])
    write_json(fixture["completion"], completion)
    pointer = json.loads(fixture["pointer"].read_text(encoding="utf-8"))
    pointer["completion_record"] = binding(fixture["completion"])
    write_json(fixture["pointer"], pointer)


def test_run_only_writes_immutable_manuscript_marker(gate, tmp_path, monkeypatch):
    fixture = fixture_tree(gate, tmp_path)
    install_dependencies(gate, monkeypatch, fixture)
    before = snapshot_files(tmp_path)

    path = gate.run_marker(fixture["args"])

    after = snapshot_files(tmp_path)
    assert set(after) == {*before, str(path.relative_to(tmp_path))}
    assert all(after[name] == value for name, value in before.items())
    assert "subprocess" not in gate.__dict__
    assert "run_checked" not in gate.__dict__
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["status"] == "manuscript_ready"
    assert record["gate_policy"]["adapter_mode"] == "validation-only"
    assert record["gate_policy"]["marker_is_final_initiative_release"] is False
    assert record["gate_policy"]["required_publication_product_count"] == 52
    assert record["ct_numerical_result"] == "pass"
    assert set(record["tools"]) == {
        "adapter",
        "corrected_report",
        "corrected_downstream",
        "corrected_science",
        "publication",
    }
    assert gate.validate_marker(fixture["args"]) == path


def test_existing_marker_refuses_before_any_validation(gate, tmp_path, monkeypatch):
    fixture = fixture_tree(gate, tmp_path)
    write_json(fixture["args"].marker, {"status": "already-present"})
    before = snapshot_files(tmp_path)
    monkeypatch.setattr(
        gate,
        "build_marker",
        lambda _args: pytest.fail("validation ran after existing-marker check"),
    )

    with pytest.raises(gate.ManuscriptReadyError, match="already exists"):
        gate.run_marker(fixture["args"])

    assert snapshot_files(tmp_path) == before


def test_marker_must_be_exact_workflow_path(gate, tmp_path, monkeypatch):
    fixture = fixture_tree(gate, tmp_path)
    install_dependencies(gate, monkeypatch, fixture)
    fixture["args"].marker = tmp_path / "external-manuscript-ready.json"

    with pytest.raises(gate.ManuscriptReadyError, match="exact workflow marker"):
        gate.build_marker(fixture["args"])


def test_interrupted_atomic_publish_leaves_no_canonical_marker(
    gate, tmp_path, monkeypatch
):
    fixture = fixture_tree(gate, tmp_path)
    install_dependencies(gate, monkeypatch, fixture)
    marker = fixture["args"].marker
    monkeypatch.setattr(
        gate.os,
        "link",
        lambda *_args: (_ for _ in ()).throw(OSError("injected publish failure")),
    )

    with pytest.raises(OSError, match="injected publish failure"):
        gate.run_marker(fixture["args"])

    assert not marker.exists()
    assert list(marker.parent.glob(f".{marker.name}.*.tmp")) == []


def test_completion_publication_must_equal_exact_workflow_publication(
    gate, tmp_path, monkeypatch
):
    fixture = fixture_tree(gate, tmp_path)
    install_dependencies(gate, monkeypatch, fixture)
    alternate = write_json(
        fixture["args"].workflow_root / "publication-intermediate/manifest.json",
        json.loads(fixture["publication_manifest"].read_text(encoding="utf-8")),
    )
    completion = json.loads(fixture["completion"].read_text(encoding="utf-8"))
    completion["publication"]["manifest"] = binding(alternate)
    write_json(fixture["completion"], completion)
    pointer = json.loads(fixture["pointer"].read_text(encoding="utf-8"))
    pointer["completion_record"] = binding(fixture["completion"])
    write_json(fixture["pointer"], pointer)

    with pytest.raises(gate.ManuscriptReadyError, match="exact workflow_root/publication"):
        gate.build_marker(fixture["args"])


def test_publication_requires_composite_hyperbolicity_root(
    gate, tmp_path, monkeypatch
):
    fixture = fixture_tree(gate, tmp_path)
    install_dependencies(gate, monkeypatch, fixture)
    alternate = fixture["root"] / "alternate-hyperbolicity"
    alternate.mkdir()
    manifest = json.loads(fixture["publication_manifest"].read_text(encoding="utf-8"))
    manifest["normalized_invocation"] = [
        str(alternate.resolve()) if value == str(fixture["hyper_root"].resolve()) else value
        for value in manifest["normalized_invocation"]
    ]
    write_json(fixture["publication_manifest"], manifest)
    rebind_completion(fixture)

    with pytest.raises(gate.ManuscriptReadyError, match="exact acceptance, science, CT"):
        gate.build_marker(fixture["args"])


def test_publication_requires_exact_csv_tex_product_inventory(
    gate, tmp_path, monkeypatch
):
    fixture = fixture_tree(gate, tmp_path)
    install_dependencies(gate, monkeypatch, fixture)
    manifest = json.loads(fixture["publication_manifest"].read_text(encoding="utf-8"))
    missing = "tables/direct_ct_health.tex"
    manifest["products"] = [
        value for value in manifest["products"] if not value["path"].endswith(missing)
    ]
    write_json(fixture["publication_manifest"], manifest)
    rebind_completion(fixture)

    with pytest.raises(gate.ManuscriptReadyError, match="product inventory differs"):
        gate.build_marker(fixture["args"])

    assert len(gate.REQUIRED_PUBLICATION_PRODUCTS) == 52
    assert "tables/direct_ct_health.tex" in gate.REQUIRED_PUBLICATION_PRODUCTS
    assert "tables/direct_ct_health.md" not in gate.REQUIRED_PUBLICATION_PRODUCTS


def test_inconclusive_reviewed_science_is_honestly_manuscript_ready(
    gate, tmp_path, monkeypatch
):
    fixture = fixture_tree(gate, tmp_path, science_result="inconclusive")
    install_dependencies(gate, monkeypatch, fixture)

    record = gate.build_marker(fixture["args"])

    assert record["status"] == "manuscript_ready"
    assert record["reviewed_science_result"] == "inconclusive"
    assert record["reviewed_science"]["result"] == "inconclusive"


def test_incomplete_ct_coverage_is_rejected(gate, tmp_path, monkeypatch):
    fixture = fixture_tree(gate, tmp_path, ct_result="inconclusive")
    install_dependencies(gate, monkeypatch, fixture)
    ct = json.loads(fixture["ct"].read_text(encoding="utf-8"))
    ct["cases"]["R17"]["native_restart_ct"]["coverage_complete"] = False
    write_json(fixture["ct"], ct)
    rebind_ct_evidence(fixture)

    with pytest.raises(gate.ManuscriptReadyError, match="coverage is incomplete"):
        gate.build_marker(fixture["args"])
