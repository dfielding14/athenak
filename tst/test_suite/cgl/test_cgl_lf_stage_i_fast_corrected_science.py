"""Focused tests for corrected/composite reviewed-science orchestration."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
TOOL = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_corrected_science.py"
CRITERIA = REPOSITORY / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
CRITERIA_REVIEW = REPOSITORY / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)
CURRENT_SCIENCE_SCOPE_LIMITATION = json.loads(
    CRITERIA_REVIEW.read_text(encoding="utf-8")
)["current_science_scope_limitation"]
ACTIVE_PASSIVE_INTERVENTION_SCOPE = json.loads(
    CRITERIA.read_text(encoding="utf-8")
)["family_gates"]["active_passive"]["intervention_scope"]


def load_tool():
    name = "cgl_lf_stage_i_fast_corrected_science_tests"
    spec = importlib.util.spec_from_file_location(name, TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tool():
    return load_tool()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def dependency_args(tool) -> dict[str, str]:
    return {
        argument: sha256(path)
        for argument, _label, path in tool.DEPENDENCY_PINS
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def context_fixture(tool, tmp_path: Path) -> dict[str, object]:
    report = tmp_path / "report"
    inventory_path = report / "inventory.json"
    identity_path = tmp_path / "identity.json"
    write_json(identity_path, {"identity": "corrected"})
    cases = {}
    contract_id = "legacy-passive-controls-unaffected-by-active-fastdisc-fix-v1"
    for case_id in tool.ALL_CASES:
        evidence_class = (
            tool.ACTIVE_EVIDENCE_CLASS
            if case_id in tool.ACTIVE_CASES
            else tool.PASSIVE_EVIDENCE_CLASS
        )
        cases[case_id] = {
            "case_id": case_id,
            "status": "complete",
            "evidence_class": evidence_class,
            "execution_authority": {
                "evidence_class": evidence_class,
                "compatibility_contract_id": (
                    contract_id if case_id in tool.PASSIVE_CASES else None
                ),
            },
        }
    inventory = {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_corrected_composite_report",
        "output": str(report.resolve()),
        "evidence_classes": {
            tool.ACTIVE_EVIDENCE_CLASS: list(tool.ACTIVE_CASES),
            tool.PASSIVE_EVIDENCE_CLASS: list(tool.PASSIVE_CASES),
        },
        "passive_compatibility": {
            "contract_id": contract_id,
            "result": "pass",
            "scope": list(tool.PASSIVE_CASES),
        },
        "cases": cases,
    }
    write_json(inventory_path, inventory)
    return {
        "inventory": inventory,
        "inventory_binding": binding(inventory_path),
        "identity_binding": binding(identity_path),
        "inventory_output": report,
        "inventory_kind": "corrected-composite",
        "active_cases": list(tool.ACTIVE_CASES),
        "passive_cases": list(tool.PASSIVE_CASES),
        "selected_cases": list(tool.ALL_CASES),
    }


def install_downstream(tool, monkeypatch, context: dict[str, object]) -> None:
    fake = SimpleNamespace(validate_inventory=lambda *_args: context)
    original = tool.load_module
    monkeypatch.setattr(
        tool,
        "load_module",
        lambda name, path: fake
        if Path(path) == tool.CORRECTED_DOWNSTREAM_TOOL
        else original(name, path),
    )


def test_exact_corrected_active_and_authenticated_passive_split_is_required(
    tool, tmp_path, monkeypatch
):
    context = context_fixture(tool, tmp_path)
    install_downstream(tool, monkeypatch, context)

    observed = tool.validate_corrected_context(
        tmp_path / "identity.json",
        tmp_path / "report/inventory.json",
        sha256(tmp_path / "report/inventory.json"),
    )
    assert observed is context

    context["inventory"]["cases"]["R03"]["evidence_class"] = tool.PASSIVE_EVIDENCE_CLASS
    with pytest.raises(tool.CorrectedScienceError, match="R03 evidence class"):
        tool.validate_corrected_context(
            tmp_path / "identity.json",
            tmp_path / "report/inventory.json",
            sha256(tmp_path / "report/inventory.json"),
        )


def test_uncontracted_or_extra_legacy_passive_evidence_is_rejected(
    tool, tmp_path, monkeypatch
):
    context = context_fixture(tool, tmp_path)
    install_downstream(tool, monkeypatch, context)
    context["inventory"]["cases"]["R06"]["execution_authority"][
        "compatibility_contract_id"
    ] = "different"

    with pytest.raises(tool.CorrectedScienceError, match="R06 lacks"):
        tool.validate_corrected_context(
            tmp_path / "identity.json",
            tmp_path / "report/inventory.json",
            sha256(tmp_path / "report/inventory.json"),
        )

    context = context_fixture(tool, tmp_path / "second")
    install_downstream(tool, monkeypatch, context)
    context["passive_cases"] = [*tool.PASSIVE_CASES, "R10"]
    with pytest.raises(tool.CorrectedScienceError, match="only authenticated legacy"):
        tool.validate_corrected_context(
            tmp_path / "second/identity.json",
            tmp_path / "second/report/inventory.json",
            sha256(tmp_path / "second/report/inventory.json"),
        )


def test_commands_delegate_exact_r02_r17_to_existing_tools(tool, tmp_path):
    context = context_fixture(tool, tmp_path)
    acceptance, science = tool.workflow_commands(
        context,
        tmp_path / "acceptance",
        tmp_path / "science",
        tool.DEFAULT_CRITERIA,
        tool.DEFAULT_CRITERIA_REVIEW,
        Path(sys.executable),
    )
    cases = ",".join(tool.ALL_CASES)

    assert acceptance[1] == str(tool.FAST_ACCEPTANCE_TOOL)
    assert "--inventory-only" in acceptance
    assert acceptance[acceptance.index("--cases") + 1] == cases
    assert science[1] == str(tool.FAST_SCIENCE_TOOL)
    assert science[science.index("--cases") + 1] == cases
    assert science[science.index("--acceptance") + 1] == str(tmp_path / "acceptance")


def test_dependency_pins_cover_every_direct_workflow_tool(tool):
    assert {
        path for _argument, _label, path in tool.DEPENDENCY_PINS
    } == {
        tool.CORRECTED_DOWNSTREAM_TOOL,
        tool.CORRECTED_REPORT_TOOL,
        tool.FAST_ACCEPTANCE_TOOL,
        tool.FAST_SCIENCE_TOOL,
        tool.REVIEWED_ACCEPTANCE_TOOL,
    }


def test_dependency_pin_mismatch_fails_closed(tool):
    args = SimpleNamespace(**dependency_args(tool))
    args.fast_science_sha256 = "0" * 64

    with pytest.raises(tool.CorrectedScienceError, match="fast science SHA-256 differs"):
        tool.verify_dependency_pins(args)


def product_fixture(tool, tmp_path: Path, *, complete: bool = True):
    context = context_fixture(tool, tmp_path)
    acceptance = tmp_path / "acceptance"
    science = tmp_path / "science"
    criteria = tmp_path / "criteria.json"
    review = tmp_path / "review.json"
    write_json(criteria, {"criteria": True})
    write_json(review, {"review": True})

    case_acceptance = {}
    reviewed_case_evidence = {}
    for case_id in tool.ALL_CASES:
        case_path = acceptance / "cases" / case_id / "case_acceptance.json"
        reviewed_path = acceptance / "cases" / case_id / "reviewed_case_evidence.json"
        write_json(case_path, {"case_id": case_id, "result": "pass"})
        write_json(reviewed_path, {"case_id": case_id, "result": "pass"})
        case_acceptance[case_id] = binding(case_path)
        reviewed_case_evidence[case_id] = binding(reviewed_path)
    campaign = acceptance / "campaign_evidence.json"
    summary = acceptance / "summary.json"
    write_json(
        campaign,
        {
            "result": "pass",
            "active_passive_intervention_scope": ACTIVE_PASSIVE_INTERVENTION_SCOPE,
            "current_science_scope_limitation": CURRENT_SCIENCE_SCOPE_LIMITATION,
        },
    )
    write_json(
        summary,
        {
            "record_type": "cgl-lf-stage-i-direct-fast-acceptance-summary",
            "selected_cases": list(tool.ALL_CASES),
            "case_results": {case_id: "pass" for case_id in tool.ALL_CASES},
            "active_passive_intervention_scope": ACTIVE_PASSIVE_INTERVENTION_SCOPE,
            "current_science_scope_limitation": CURRENT_SCIENCE_SCOPE_LIMITATION,
        },
    )
    acceptance_provenance = {
        "outputs": {
            "summary": binding(summary),
            "campaign_evidence": binding(campaign),
            "case_acceptance": case_acceptance,
            "reviewed_case_evidence": reviewed_case_evidence,
        }
    }
    write_json(acceptance / "provenance.json", acceptance_provenance)

    science_record = {
        "record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-comparisons",
        "selected_cases": list(tool.ALL_CASES),
        "result": "pass",
        "active_passive_intervention_scope": ACTIVE_PASSIVE_INTERVENTION_SCOPE,
        "current_science_scope_limitation": CURRENT_SCIENCE_SCOPE_LIMITATION,
        "case_dispositions": {
            case_id: {
                "inventory_status": "complete",
                "diagnostics_available": complete,
                "snapshot_products_authenticated": complete,
            }
            for case_id in tool.ALL_CASES
        },
        "provenance": {},
    }
    case_lineages = {}
    case_diagnostics = {}
    for case_id in tool.ALL_CASES:
        lineage = tmp_path / "report" / "cases" / case_id / "lineage.json"
        diagnostics = tmp_path / "report" / "cases" / case_id / "diagnostics.json"
        write_json(lineage, {"case_id": case_id})
        write_json(diagnostics, {"case_id": case_id})
        case_lineages[case_id] = binding(lineage)
        case_diagnostics[case_id] = binding(diagnostics)
    science_record["provenance"] = {
        "inventory": context["inventory_binding"],
        "acceptance_provenance": binding(acceptance / "provenance.json"),
        "acceptance_campaign_evidence": binding(campaign),
        "aggregator": binding(tool.FAST_SCIENCE_TOOL),
        "case_acceptance": case_acceptance,
        "case_lineages": case_lineages,
        "case_diagnostics": case_diagnostics,
    }
    write_json(science / "science.json", science_record)
    table_bindings = []
    for relative in tool.EXPECTED_SCIENCE_TABLES:
        path = science / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n", encoding="utf-8")
        table_bindings.append(binding(path))
    write_json(
        science / "provenance.json",
        {
            "record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-provenance",
            "inputs": science_record["provenance"],
            "outputs": {
                "science": binding(science / "science.json"),
                "tables": table_bindings,
            },
        },
    )
    return context, acceptance, science, criteria, review, acceptance_provenance


def install_science_modules(
    tool,
    monkeypatch,
    acceptance_provenance: dict[str, object],
    campaign: dict[str, object] | None = None,
):
    campaign_record = campaign or {
        "result": "pass",
        "active_passive_intervention_scope": ACTIVE_PASSIVE_INTERVENTION_SCOPE,
        "current_science_scope_limitation": CURRENT_SCIENCE_SCOPE_LIMITATION,
    }
    fast_science = SimpleNamespace(
        validate_inventory=lambda *_args: ({}, None, {}, {}),
        validate_acceptance_root=lambda *_args: (
            {case_id: {"result": "pass"} for case_id in tool.ALL_CASES},
            campaign_record,
            acceptance_provenance,
            {},
        ),
    )
    reviewed = SimpleNamespace(
        CURRENT_SCIENCE_SCOPE_LIMITATION=CURRENT_SCIENCE_SCOPE_LIMITATION,
        ACTIVE_PASSIVE_INTERVENTION_SCOPE=ACTIVE_PASSIVE_INTERVENTION_SCOPE,
        load_validated_policy=lambda *_args: {
            "current_science_scope_limitation": CURRENT_SCIENCE_SCOPE_LIMITATION,
            "criteria": {
                "family_gates": {
                    "active_passive": {
                        "intervention_scope": ACTIVE_PASSIVE_INTERVENTION_SCOPE
                    }
                }
            },
        },
        verify_evidence_digest=lambda *_args: None,
    )
    original = tool.load_module
    monkeypatch.setattr(
        tool,
        "load_module",
        lambda name, path: (
            fast_science
            if Path(path) == tool.FAST_SCIENCE_TOOL
            else reviewed
            if Path(path) == tool.REVIEWED_ACCEPTANCE_TOOL
            else original(name, path)
        ),
    )


def test_final_product_validation_binds_existing_tools_and_exact_coverage(
    tool, tmp_path, monkeypatch
):
    context, acceptance, science, criteria, review, provenance = product_fixture(
        tool, tmp_path
    )
    install_science_modules(tool, monkeypatch, provenance)

    record = tool.validate_products(context, acceptance, science, criteria, review)

    assert record["status"] == "complete"
    assert record["result"] == "pass"
    assert record["case_classification"]["corrected_active"] == list(tool.ACTIVE_CASES)
    assert record["case_classification"]["authenticated_legacy_passive"] == list(
        tool.PASSIVE_CASES
    )
    assert record["active_passive_intervention_scope"] == (
        ACTIVE_PASSIVE_INTERVENTION_SCOPE
    )
    assert record["current_science_scope_limitation"] == (
        CURRENT_SCIENCE_SCOPE_LIMITATION
    )
    assert set(record["science"]["tables"]) == set(tool.EXPECTED_SCIENCE_TABLES)


@pytest.mark.parametrize(
    ("source", "field", "message"),
    [
        (
            "campaign",
            "current_science_scope_limitation",
            "acceptance campaign current science scope limitation differs",
        ),
        (
            "campaign",
            "active_passive_intervention_scope",
            "acceptance campaign active/passive intervention scope differs",
        ),
        (
            "summary",
            "current_science_scope_limitation",
            "acceptance summary current science scope limitation differs",
        ),
        (
            "summary",
            "active_passive_intervention_scope",
            "acceptance summary active/passive intervention scope differs",
        ),
        (
            "science",
            "current_science_scope_limitation",
            "direct reviewed science current science scope limitation differs",
        ),
        (
            "science",
            "active_passive_intervention_scope",
            "direct reviewed science active/passive intervention scope differs",
        ),
    ],
)
def test_final_product_validation_requires_exact_scope_records(
    tool, tmp_path, monkeypatch, source, field, message
):
    context, acceptance, science, criteria, review, provenance = product_fixture(
        tool, tmp_path
    )
    campaign = {
        "result": "pass",
        "active_passive_intervention_scope": json.loads(
            json.dumps(ACTIVE_PASSIVE_INTERVENTION_SCOPE)
        ),
        "current_science_scope_limitation": json.loads(
            json.dumps(CURRENT_SCIENCE_SCOPE_LIMITATION)
        ),
    }
    if source == "campaign":
        value = campaign
    else:
        path = acceptance / "summary.json" if source == "summary" else science / "science.json"
        value = json.loads(path.read_text(encoding="utf-8"))
    if field == "current_science_scope_limitation":
        value[field]["full_scope_independent_review_complete"] = True
    else:
        value[field]["excluded_interpretation"] = "none"
    if source != "campaign":
        write_json(path, value)
    install_science_modules(tool, monkeypatch, provenance, campaign)

    with pytest.raises(tool.CorrectedScienceError, match=message):
        tool.validate_products(context, acceptance, science, criteria, review)


def test_incomplete_downstream_science_products_cannot_be_finalized(
    tool, tmp_path, monkeypatch
):
    context, acceptance, science, criteria, review, provenance = product_fixture(
        tool, tmp_path, complete=False
    )
    install_science_modules(tool, monkeypatch, provenance)

    with pytest.raises(tool.CorrectedScienceError, match="R02 lacks complete"):
        tool.validate_products(context, acceptance, science, criteria, review)


def test_run_orders_acceptance_then_science_and_revalidates(
    tool, tmp_path, monkeypatch
):
    context = context_fixture(tool, tmp_path)
    calls = []
    monkeypatch.setattr(tool, "validate_corrected_context", lambda *_args: context)
    monkeypatch.setattr(
        tool,
        "validate_output_layout",
        lambda *_args: (tmp_path / "acceptance", tmp_path / "science"),
    )
    monkeypatch.setattr(
        tool,
        "workflow_commands",
        lambda *_args: (["acceptance-command"], ["science-command"]),
    )
    monkeypatch.setattr(tool, "run_checked", lambda command: calls.append(command))
    monkeypatch.setattr(
        tool,
        "validate_products",
        lambda *_args: {"schema": tool.SCHEMA, "result": "pass"},
    )
    monkeypatch.setattr(
        tool,
        "load_json",
        lambda *_args: {"schema": tool.SCHEMA, "result": "pass"},
    )
    args = SimpleNamespace(
        identity=tmp_path / "identity.json",
        inventory=tmp_path / "report/inventory.json",
        inventory_sha256="0" * 64,
        acceptance_output=tmp_path / "acceptance",
        science_output=tmp_path / "science",
        criteria=tool.DEFAULT_CRITERIA,
        criteria_review=tool.DEFAULT_CRITERIA_REVIEW,
        python=Path(sys.executable),
        **dependency_args(tool),
    )

    path = tool.run_workflow(args)

    assert calls == [["acceptance-command"], ["science-command"]]
    assert path == tmp_path / "science" / tool.RECORD_NAME
