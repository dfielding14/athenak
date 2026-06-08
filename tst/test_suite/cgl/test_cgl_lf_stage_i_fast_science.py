"""Focused regressions for the direct-fast final science aggregator."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
AGGREGATOR = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_science.py"
FAST_ACCEPTANCE = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_acceptance.py"
FAST_REPORT = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_report.py"
REVIEWED = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
ANALYZER = REPOSITORY / "scripts/analyze_cgl_lf_paper.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def science():
    return load_module("cgl_lf_stage_i_fast_science_test", AGGREGATOR)


@pytest.fixture(scope="module")
def reviewed():
    return load_module("cgl_lf_stage_i_fast_science_reviewed_test", REVIEWED)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def write_history(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# [1]=time [2]=mass\n0 1\n4 1\n6 1\n8 1\n10 1\n",
        encoding="utf-8",
    )


def scalar(mean: float, *, deviation: float = 0.1, error: float = 0.01) -> dict[str, object]:
    return {
        "mean": mean,
        "standard_deviation": deviation,
        "standard_error": error,
        "block_bootstrap_interval_95": [mean - 1.96 * error, mean + 1.96 * error],
        "effective_sample_count": 4.0,
        "independent_time_block_count": 3,
        "gap_adequacy": "pass",
        "sample_count": 5,
    }


def analyzer_scalar(
    mean: float, *, deviation: float = 0.1, error: float = 0.01
) -> dict[str, object]:
    return {
        "mean": mean,
        "standard_deviation": deviation,
        "standard_error": error,
    }


def metric_record(mean: float) -> dict[str, object]:
    return {
        "history": "user",
        "column": "fixture",
        "windows": {
            "full": {"result": "available", "statistics": scalar(mean)},
            "early": {"result": "available", "statistics": scalar(mean)},
            "late": {"result": "available", "statistics": scalar(mean)},
        },
        "sampling_adequacy": "pass",
        "stationarity": {"result": "pass"},
    }


def reviewed_metrics(active_mean: float) -> dict[str, object]:
    means = {
        "abs_dp": active_mean,
        "mirror_occupancy": active_mean,
        "firehose_occupancy": 0.0,
        "unstable_occupancy": active_mean,
        "kinetic": 1.0,
        "magnetic": 1.0,
        "beta": 10.0,
        "nu_eff": 2.0,
    }
    return {
        metric: {
            window: scalar(mean)
            for window in ("full", "early", "late")
        }
        for metric, mean in means.items()
    }


def reviewed_case_evidence(
    reviewed,
    root: Path,
    case_id: str,
    case_name: str,
    mhd: Path,
    user: Path,
    policy: dict[str, object],
    active_mean: float,
) -> Path:
    path = root / "acceptance" / "cases" / case_id / "reviewed_case_evidence.json"
    write_json(
        path,
        reviewed.seal_evidence({
            "schema_version": 1,
            "record_type": "stage-i-scientific-case-evidence",
            "case_id": case_id,
            "case_name": case_name,
            "result": "pass",
            "evaluation_inputs": {
                "mhd_history": binding(mhd),
                "user_history": binding(user),
                "diagnostics": {},
            },
            "metrics": reviewed_metrics(active_mean),
            "analyzer_metrics": {"peak_alignment": analyzer_scalar(active_mean)},
            "convergence_products": {},
            "panel_products": [],
            "gates": (
                [
                    reviewed.gate(
                        "landau_fluid_activity",
                        "pass",
                        reason="fixture active LF transport",
                        observations={"lf_qface_increment": 1.0},
                    )
                ]
                if case_id in ("R12", "R02", "R13")
                else []
            ),
            "provenance": {
                "criteria": policy["criteria_binding"],
                "criteria_review": policy["review_binding"],
                "acceptance_utility": binding(REVIEWED),
                "inputs": [binding(mhd), binding(user)],
            },
        }),
    )
    return path


def build_fixture(
    science,
    reviewed,
    root: Path,
    *,
    case_results: dict[str, str],
    diagnostics_case: str | None = None,
) -> dict[str, Path]:
    policy = reviewed.load_validated_policy(
        science.DEFAULT_CRITERIA, science.DEFAULT_CRITERIA_REVIEW
    )
    manifest = policy["manifest"]
    manifest_cases = {str(item["id"]): item for item in manifest["cases"]}
    report = root / "report"
    inventory_cases: dict[str, object] = {}
    case_inputs: dict[str, tuple[Path, Path, Path]] = {}
    for case_id in case_results:
        case = manifest_cases[case_id]
        case_dir = report / "cases" / case_id
        mhd = case_dir / "history" / "fixture.mhd.hst"
        user = case_dir / "history" / "fixture.user.hst"
        write_history(mhd)
        write_history(user)
        snapshots = case_dir / "snapshots.json"
        write_json(
            snapshots,
            {"snapshot_count": 0, "complete_snapshot_count": 0, "snapshots": []},
        )
        lineage = {
            "schema_version": 1,
            "case_id": case_id,
            "case_name": case["name"],
            "status": "complete",
            "histories": {
                "mhd": {
                    "available": True,
                    "path": str(mhd.resolve()),
                    "binding": binding(mhd),
                },
                "user": {
                    "available": True,
                    "path": str(user.resolve()),
                    "binding": binding(user),
                },
            },
        }
        lineage_path = case_dir / "lineage.json"
        write_json(lineage_path, lineage)
        inventory_cases[case_id] = lineage
        case_inputs[case_id] = (lineage_path, mhd, user)
        if diagnostics_case == case_id:
            write_json(
                case_dir / "diagnostics.json",
                {
                    "schema_version": 1,
                    "case_id": case_id,
                    "case_name": case["name"],
                    "analysis_status": "complete",
                    "snapshot_analysis_status": "not_yet_available",
                    "snapshot_ensemble": {"snapshot_count": 0},
                    "snapshots": {},
                    "provenance": {
                        "lineage": binding(lineage_path),
                        "snapshot_index": binding(snapshots),
                        "merged_mhd_history": binding(mhd),
                        "merged_user_history": binding(user),
                        "analyzer": binding(ANALYZER),
                        "adapter": binding(FAST_REPORT),
                    },
                },
            )
    inventory = report / "inventory.json"
    write_json(
        inventory,
        {
            "schema_version": 1,
            "output": str(report.resolve()),
            "matrix": binding(Path(str(policy["verified_sources"]["stage_i_manifest"]["path"]))),
            "adapter": binding(FAST_REPORT),
            "cases": inventory_cases,
        },
    )

    acceptance = root / "acceptance"
    summary_bindings: dict[str, object] = {}
    reviewed_bindings: dict[str, object] = {}
    for case_id, result in case_results.items():
        lineage_path, mhd, user = case_inputs[case_id]
        active_mean = 0.0 if case_id in ("R02", "R03", "R04", "R05") else 1.0
        reviewed_path = (
            reviewed_case_evidence(
                reviewed,
                root,
                case_id,
                str(manifest_cases[case_id]["name"]),
                mhd,
                user,
                policy,
                active_mean,
            )
            if result == "pass"
            else None
        )
        summary = {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-case-acceptance",
            "case_id": case_id,
            "case_name": manifest_cases[case_id]["name"],
            "result": result,
            "reason": f"fixture {result}",
            "active_passive_intervention_scope": (
                reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
            ),
            "current_science_scope_limitation": (
                reviewed.CURRENT_SCIENCE_SCOPE_LIMITATION
            ),
            "health": {"result": "pass", "complete_to_target": True},
            "scope": {
                "classification": "standard_claim_scope",
                "campaign_interpretation_eligible": True,
            },
            "comparison_evidence": {"result": "pass"},
            "reviewed_case_evidence": (
                {"result": "pass", "binding": binding(reviewed_path)}
                if reviewed_path is not None
                else {"result": "inconclusive", "binding": None}
            ),
            "history_statistics": {
                "abs_dp": metric_record(active_mean),
                "mirror_occupancy": metric_record(active_mean),
                "firehose_occupancy": metric_record(0.0),
                "unstable_occupancy": metric_record(active_mean),
                "kinetic": metric_record(1.0),
                "magnetic": metric_record(1.0),
                "beta": metric_record(10.0),
                "nu_eff": metric_record(2.0),
            },
            "provenance": {
                "lineage": binding(lineage_path),
                "histories": {"mhd": binding(mhd), "user": binding(user)},
            },
        }
        summary_path = acceptance / "cases" / case_id / "case_acceptance.json"
        write_json(summary_path, summary)
        summary_bindings[case_id] = binding(summary_path)
        if reviewed_path is not None:
            reviewed_bindings[case_id] = binding(reviewed_path)
    campaign_path = acceptance / "campaign_evidence.json"
    write_json(
        campaign_path,
        reviewed.seal_evidence({
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-campaign-evidence",
            "result": (
                "pass" if all(value == "pass" for value in case_results.values())
                else "inconclusive"
            ),
            "case_results": case_results,
            "active_passive_intervention_scope": (
                reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
            ),
            "current_science_scope_limitation": (
                reviewed.CURRENT_SCIENCE_SCOPE_LIMITATION
            ),
            "provenance": {
                "criteria": policy["criteria_binding"],
                "criteria_review": policy["review_binding"],
                "reviewed_acceptance_utility": binding(REVIEWED),
            },
        }),
    )
    provenance = acceptance / "provenance.json"
    write_json(
        provenance,
        {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-acceptance-provenance",
            "inputs": {
                "inventory": binding(inventory),
                "criteria": policy["criteria_binding"],
                "criteria_review": policy["review_binding"],
                "reviewed_acceptance_utility": binding(REVIEWED),
                "driver": binding(FAST_ACCEPTANCE),
            },
            "outputs": {
                "campaign_evidence": binding(campaign_path),
                "case_acceptance": summary_bindings,
                "reviewed_case_evidence": reviewed_bindings,
            },
        },
    )
    return {
        "inventory": inventory,
        "acceptance": acceptance,
        "report": report,
        "criteria": science.DEFAULT_CRITERIA,
        "criteria_review": science.DEFAULT_CRITERIA_REVIEW,
    }


def aggregate(science, fixture: dict[str, Path], output: Path, cases: list[str]):
    return science.aggregate_science(
        fixture["inventory"],
        fixture["acceptance"],
        output,
        fixture["criteria"],
        fixture["criteria_review"],
        cases,
    )[0]


def test_partial_campaign_is_honest_and_cannot_pass(science, reviewed, tmp_path):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )

    result = aggregate(science, fixture, tmp_path / "science", ["R02"])

    assert result["result"] == "inconclusive"
    assert result["case_dispositions"]["R02"]["claim_eligible"] is False
    pair = result["families"]["active_passive"]["R02_R06"]
    assert pair["result"] == "inconclusive"
    assert pair["claim_eligible"] is False
    assert pair["intervention_scope"] == reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
    assert result["active_passive_intervention_scope"] == (
        reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
    )
    assert result["current_science_scope_limitation"] == (
        reviewed.CURRENT_SCIENCE_SCOPE_LIMITATION
    )
    assert (tmp_path / "science/tables/contrasts.csv").is_file()
    assert (tmp_path / "science/tables/mks24.md").is_file()


def test_passing_pair_uses_descriptive_coherent_standardized_effect(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science,
        reviewed,
        tmp_path,
        case_results={"R02": "pass", "R06": "pass"},
    )

    result = aggregate(science, fixture, tmp_path / "science", ["R02", "R06"])

    pair = result["families"]["active_passive"]["R02_R06"]
    assert pair["claim_eligible"] is True
    assert pair["result"] == "pass"
    assert pair["intervention_scope"] == reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
    available = [item for item in pair["metrics"] if item.get("available")]
    assert all(item.get("direction_coherent") is True for item in available)
    assert any(abs(float(item["standardized_effect"])) >= 0.5 for item in available)
    assert all(
        item.get("claim_scope") == "descriptive_within_realization"
        for item in available
    )
    assert all(
        item.get("intervention_scope") == reviewed.ACTIVE_PASSIVE_INTERVENTION_SCOPE
        for item in available
    )
    assert all("two_sided_p" not in item and "holm_significant" not in item for item in available)
    assert pair["population_inference"] == reviewed.POPULATION_INFERENCE_LIMITATION
    contrasts_csv = (tmp_path / "science/tables/contrasts.csv").read_text(
        encoding="utf-8"
    )
    contrasts_md = (tmp_path / "science/tables/contrasts.md").read_text(
        encoding="utf-8"
    )
    for table in (contrasts_csv, contrasts_md):
        assert "intervention_estimand" in table
        assert "total_effect_of_enabling_active_cgl" in table
        assert "anisotropic_stress_alone" in table
        assert "realized_forcing_after_trajectory_divergence" in table


def test_forged_campaign_evidence_digest_fails_closed(science, reviewed, tmp_path):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )
    campaign = fixture["acceptance"] / "campaign_evidence.json"
    value = json.loads(campaign.read_text(encoding="utf-8"))
    value["result"] = "fail"
    write_json(campaign, value)
    provenance = fixture["acceptance"] / "provenance.json"
    declared = json.loads(provenance.read_text(encoding="utf-8"))
    declared["outputs"]["campaign_evidence"] = binding(campaign)
    write_json(provenance, declared)

    with pytest.raises(science.ScienceError, match="forged"):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_campaign_current_science_scope_limitation_must_be_exact(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )
    campaign = fixture["acceptance"] / "campaign_evidence.json"
    value = json.loads(campaign.read_text(encoding="utf-8"))
    value.pop("evidence_digest")
    value["current_science_scope_limitation"][
        "full_scope_independent_review_complete"
    ] = True
    write_json(campaign, reviewed.seal_evidence(value))
    provenance = fixture["acceptance"] / "provenance.json"
    declared = json.loads(provenance.read_text(encoding="utf-8"))
    declared["outputs"]["campaign_evidence"] = binding(campaign)
    write_json(provenance, declared)

    with pytest.raises(
        science.ScienceError,
        match="campaign evidence current science scope limitation differs",
    ):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_case_active_passive_intervention_scope_must_be_exact(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )
    summary = fixture["acceptance"] / "cases/R02/case_acceptance.json"
    value = json.loads(summary.read_text(encoding="utf-8"))
    value["active_passive_intervention_scope"]["excluded_interpretation"] = "none"
    write_json(summary, value)
    provenance = fixture["acceptance"] / "provenance.json"
    declared = json.loads(provenance.read_text(encoding="utf-8"))
    declared["outputs"]["case_acceptance"]["R02"] = binding(summary)
    write_json(provenance, declared)

    with pytest.raises(
        science.ScienceError,
        match="R02 direct-fast acceptance active/passive intervention scope differs",
    ):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_stale_diagnostics_provenance_fails_closed(science, reviewed, tmp_path):
    fixture = build_fixture(
        science,
        reviewed,
        tmp_path,
        case_results={"R02": "inconclusive"},
        diagnostics_case="R02",
    )
    diagnostics = fixture["report"] / "cases/R02/diagnostics.json"
    value = json.loads(diagnostics.read_text(encoding="utf-8"))
    value["provenance"]["lineage"]["sha256"] = "0" * 64
    write_json(diagnostics, value)

    with pytest.raises(science.ScienceError, match="diagnostics lineage"):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_stale_acceptance_inventory_binding_fails_closed(science, reviewed, tmp_path):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )
    inventory = json.loads(fixture["inventory"].read_text(encoding="utf-8"))
    inventory["forged_after_acceptance"] = True
    write_json(fixture["inventory"], inventory)

    with pytest.raises(science.ScienceError, match="acceptance input inventory"):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_claim_grade_metrics_come_from_sealed_reviewed_evidence(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science,
        reviewed,
        tmp_path,
        case_results={"R02": "pass", "R06": "pass"},
    )
    summary = fixture["acceptance"] / "cases/R02/case_acceptance.json"
    value = json.loads(summary.read_text(encoding="utf-8"))
    value["history_statistics"]["abs_dp"] = metric_record(0.0)
    write_json(summary, value)
    provenance = fixture["acceptance"] / "provenance.json"
    declared = json.loads(provenance.read_text(encoding="utf-8"))
    declared["outputs"]["case_acceptance"]["R02"] = binding(summary)
    write_json(provenance, declared)

    result = aggregate(science, fixture, tmp_path / "science", ["R02", "R06"])

    pair = result["families"]["active_passive"]["R02_R06"]
    abs_dp = next(item for item in pair["metrics"] if item.get("metric") == "abs_dp")
    assert pair["result"] == "pass"
    assert abs_dp["difference"] == -1.0


def test_lf_strength_gate_excludes_r06_and_accepts_near_insensitivity(
    science, reviewed,
):
    policy = reviewed.load_validated_policy(
        science.DEFAULT_CRITERIA, science.DEFAULT_CRITERIA_REVIEW
    )

    def case(value: float) -> dict[str, object]:
        return {
            "analyzer_metrics": {"peak_alignment": analyzer_scalar(value)},
            "gates": [
                reviewed.gate(
                    "landau_fluid_activity",
                    "pass",
                    reason="fixture active LF transport",
                    observations={"lf_qface_increment": 1.0},
                )
            ],
        }

    cases = {
        "R12": case(0.5),
        "R02": case(0.5),
        "R13": case(0.5),
        "R06": case(0.1),
    }
    gate = science.lf_strength_gate(
        reviewed,
        policy,
        cases,
        {case_id: True for case_id in cases},
    )

    assert gate["result"] == "pass"
    responses = gate["observations"]["responses"]
    assert [record["pair"] for record in responses] == [
        ["R12", "R02"], ["R13", "R02"]
    ]
    assert all(record["response_classification"] == "near_insensitive" for record in responses)
    assert all("R06" not in record["pair"] for record in responses)


def test_finite_limiter_gate_is_descriptive_and_claim_scoped(
    science, reviewed,
):
    cases = {
        "R14": {"metrics": {"nu_eff": reviewed_metrics(1.0)["nu_eff"]}},
        "R15": {"metrics": {"nu_eff": reviewed_metrics(2.0)["nu_eff"]}},
    }
    gate = science.limiter_ordering_gate(
        reviewed,
        cases,
        {"R14": True, "R15": True},
        {
            "R14": {"scope": {"classification": "scoped_nonfatal_hard_bound_variant"}},
            "R15": {"scope": {"classification": "scoped_nonfatal_hard_bound_variant"}},
        },
    )

    assert gate["observations"]["claim_scope"] == "descriptive_within_trajectory"
    assert "descriptive_block_lower_95" in gate["observations"]
    assert "confidence" not in gate["reason"].lower()
    assert "signific" not in gate["reason"].lower()
    assert "population" not in str(gate).lower()


def test_stale_reviewed_evaluation_input_fails_closed(science, reviewed, tmp_path):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "pass"}
    )
    diagnostic = tmp_path / "reviewed-diagnostic.json"
    write_json(diagnostic, {"value": 1})
    evidence = fixture["acceptance"] / "cases/R02/reviewed_case_evidence.json"
    value = json.loads(evidence.read_text(encoding="utf-8"))
    value.pop("evidence_digest")
    value["evaluation_inputs"]["diagnostics"] = {"full": binding(diagnostic)}
    write_json(evidence, reviewed.seal_evidence(value))
    summary = fixture["acceptance"] / "cases/R02/case_acceptance.json"
    summary_value = json.loads(summary.read_text(encoding="utf-8"))
    summary_value["reviewed_case_evidence"]["binding"] = binding(evidence)
    write_json(summary, summary_value)
    provenance = fixture["acceptance"] / "provenance.json"
    declared = json.loads(provenance.read_text(encoding="utf-8"))
    declared["outputs"]["case_acceptance"]["R02"] = binding(summary)
    declared["outputs"]["reviewed_case_evidence"]["R02"] = binding(evidence)
    write_json(provenance, declared)
    write_json(diagnostic, {"value": 2})

    with pytest.raises(science.ScienceError, match="reviewed evaluation inputs"):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_forged_snapshot_ensemble_without_records_fails_closed(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science,
        reviewed,
        tmp_path,
        case_results={"R02": "inconclusive"},
        diagnostics_case="R02",
    )
    diagnostics = fixture["report"] / "cases/R02/diagnostics.json"
    value = json.loads(diagnostics.read_text(encoding="utf-8"))
    value["snapshot_analysis_status"] = "complete"
    value["snapshot_ensemble"] = {"snapshot_count": 1}
    value["compat"] = {"snapshot_ensemble": {"snapshot_count": 1}}
    write_json(diagnostics, value)

    with pytest.raises(science.ScienceError, match="snapshot record count"):
        aggregate(science, fixture, tmp_path / "science", ["R02"])


def test_output_cannot_overlap_report_or_acceptance(science, reviewed, tmp_path):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )

    with pytest.raises(science.ScienceError, match="separate"):
        aggregate(science, fixture, fixture["report"] / "science", ["R02"])
    with pytest.raises(science.ScienceError, match="separate"):
        aggregate(science, fixture, fixture["acceptance"] / "science", ["R02"])


def test_nested_symlink_output_escape_is_rejected_before_writing(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science, reviewed, tmp_path, case_results={"R02": "inconclusive"}
    )
    output = tmp_path / "science"
    outside = tmp_path / "outside"
    outside.mkdir()
    (output / "tables").mkdir(parents=True)
    (output / "tables/escape").symlink_to(outside, target_is_directory=True)

    with pytest.raises(science.ScienceError, match="nested symlink escape"):
        aggregate(science, fixture, output, ["R02"])

    assert not (outside / "cases.csv").exists()


def test_outputs_are_deterministic_and_inputs_remain_unchanged(
    science, reviewed, tmp_path
):
    fixture = build_fixture(
        science,
        reviewed,
        tmp_path,
        case_results={"R02": "pass", "R06": "pass"},
    )
    tracked = [
        fixture["inventory"],
        fixture["acceptance"] / "provenance.json",
        fixture["acceptance"] / "campaign_evidence.json",
    ]
    before = {path: (path.stat().st_mtime_ns, sha256(path)) for path in tracked}
    first = tmp_path / "science-a"
    second = tmp_path / "science-b"

    aggregate(science, fixture, first, ["R02", "R06"])
    aggregate(science, fixture, second, ["R02", "R06"])

    compared = [
        "science.json",
        "tables/cases.csv",
        "tables/cases.md",
        "tables/contrasts.csv",
        "tables/contrasts.md",
        "tables/gates.csv",
        "tables/gates.md",
        "tables/resolution.csv",
        "tables/resolution.md",
        "tables/mks24.csv",
        "tables/mks24.md",
    ]
    assert all((first / name).read_bytes() == (second / name).read_bytes() for name in compared)
    assert before == {path: (path.stat().st_mtime_ns, sha256(path)) for path in tracked}
