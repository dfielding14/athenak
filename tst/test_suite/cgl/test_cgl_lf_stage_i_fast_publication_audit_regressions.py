"""Regression tests for audited direct-fast Stage I publication blockers."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
RENDERER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_publication.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def publication():
    return load_module("cgl_lf_stage_i_fast_publication_regressions", RENDERER)


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
        "path": str(path.absolute()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def seal(record: dict[str, object]) -> dict[str, object]:
    body = dict(record)
    payload = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    body["evidence_digest"] = {
        "method": "sha256-canonical-json-without-evidence_digest-v1",
        "sha256": hashlib.sha256(payload).hexdigest(),
    }
    return body


def seal_ct(record: dict[str, object]) -> dict[str, object]:
    body = dict(record)
    payload = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    body["evidence_digest"] = {
        "method": "sha256-canonical-json-without-evidence-digest",
        "sha256": hashlib.sha256(payload).hexdigest(),
    }
    return body


def write_science_record(root: Path, analysis: Path) -> Path:
    inventory = analysis / "inventory.json"
    support = root / "support.txt"
    support.parent.mkdir(parents=True, exist_ok=True)
    support.write_text("reviewed science support\n", encoding="utf-8")
    provenance = {
        "inventory": binding(inventory),
        "acceptance_provenance": binding(support),
        "acceptance_campaign_evidence": binding(support),
        "criteria": binding(support),
        "criteria_review": binding(support),
        "reviewed_acceptance_utility": binding(support),
        "fast_report_utility": binding(support),
        "paper_analyzer": binding(support),
        "aggregator": binding(support),
        "case_acceptance": {},
        "case_lineages": {},
        "case_diagnostics": {},
    }
    record = seal({
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-comparisons",
        "authority": "non-authorizing-direct-fast-scientific-assessment",
        "release_authorizing": False,
        "result": "inconclusive",
        "selected_cases": ["R02", "R04"],
        "case_dispositions": {
            "R02": {"claim_eligible": True, "acceptance_result": "pass"},
            "R04": {"claim_eligible": True, "acceptance_result": "pass"},
        },
        "families": {
            "forcing": {
                "R02_R04": {
                    "left": "R02",
                    "right": "R04",
                    "result": "pass",
                    "claim_eligible": True,
                    "metrics": [{
                        "metric": "kinetic",
                        "available": True,
                        "left_mean": 1.0,
                        "right_mean": 1.4,
                        "difference": -0.4,
                        "combined_standard_error": 0.1,
                        "z_score": -4.0,
                        "two_sided_p": 0.001,
                        "standardized_effect": -1.25,
                        "holm_threshold": 0.01,
                        "holm_significant": True,
                    }],
                }
            }
        },
        "resolution": {
            "result": "inconclusive",
            "reason": "partial campaign",
            "limits": {"common_k_perp_over_pi": [4, 24]},
            "observations": [{
                "kind": "curve",
                "product": "velocity_spectrum_shape",
                "available": False,
                "reason": "R16/R17 unavailable",
            }],
        },
        "mks24": {
            "result": "inconclusive",
            "panels": {
                "fig11bottom": {
                    "result": "inconclusive",
                    "products": [{
                        "product_id": "fig11_alignment_active_alfvenic_beta10_nperp192",
                        "case_id": "R02",
                        "source": "unavailable",
                        "result": "inconclusive",
                        "reason": "snapshot analysis unavailable",
                    }],
                }
            },
        },
        "gates": [{
            "name": "forcing:R02:R04",
            "result": "pass",
            "reason": "reviewed contrast passes",
            "observations": [{"metric": "kinetic", "holm_significant": True}],
            "limits": {"alpha": 0.05},
        }, {
            "name": "R16_R02_R17_resolution_convergence",
            "result": "inconclusive",
            "reason": "partial campaign",
            "observations": [],
            "limits": {"common_k_perp_over_pi": [4, 24]},
        }],
        "provenance": provenance,
    })
    path = root / "science.json"
    write_json(path, record)
    write_json(
        root / "provenance.json",
        {
            "schema_version": 1,
            "record_type": "cgl-lf-stage-i-direct-fast-reviewed-science-provenance",
            "inputs": provenance,
            "outputs": {"science": binding(path), "tables": []},
        },
    )
    return path


def write_ct_record(root: Path, analysis: Path) -> Path:
    support = root / "ct_support.txt"
    support.parent.mkdir(parents=True, exist_ok=True)
    support.write_text("direct CT support\n", encoding="utf-8")
    record = seal_ct({
        "schema_version": 2,
        "record_type": "stage-i-direct-fast-ct-audit",
        "result": "pass",
        "inventory": binding(analysis / "inventory.json"),
        "source_bindings": {"audit_utility": binding(support)},
        "claim_boundary": {
            "campaign_authority_eligible": False,
            "release_authorizing": False,
        },
        "selection": {"cases": ["R02"]},
        "summary": {"requested_case_count": 1, "ct_pass_case_count": 1},
        "cases": {
            "R02": {
                "ct_result": "pass",
                "ct_evidence_available": True,
                "ct_claim_supported": True,
                "provenance_authenticated": True,
                "reason": "sampled native restart CT-divB is below threshold",
                "native_restart_ct": {
                    "coverage_complete": True,
                    "ct_evidence_available": True,
                    "ct_claim_supported": True,
                    "maximum_normalized_ct_divb": 2.4e-13,
                    "normalized_ct_divb_lt": 1.0e-10,
                    "campaign_authority_eligible": False,
                    "release_authorizing": False,
                },
            }
        },
    })
    path = root / "ct_audit.json"
    write_json(path, record)
    return path


def write_history(path: Path, kinetic: float = 1.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# [1]=time [2]=volume [3]=kinetic\n"
        f"0 1 {kinetic}\n4 1 {kinetic}\n10 1 {kinetic}\n",
        encoding="utf-8",
    )


def case_fixture(
    analysis: Path, case_id: str, *, kinetic: float = 1.0
) -> tuple[Path, Path, Path]:
    case_dir = analysis / "cases" / case_id
    mhd = case_dir / "history" / "fixture.mhd.hst"
    user = case_dir / "history" / "fixture.user.hst"
    write_history(mhd, kinetic)
    write_history(user, kinetic)
    lineage = case_dir / "lineage.json"
    write_json(
        lineage,
        {
            "case_id": case_id,
            "status": "complete",
            "final_time": 10.0,
            "lineage_variants": [],
            "model_choices": {"cgl_lf_strict_admissibility": "true"},
            "histories": {
                "mhd": {"path": str(mhd.absolute())},
                "user": {"path": str(user.absolute())},
            },
        },
    )
    write_json(
        case_dir / "diagnostics.json",
        {
            "health": {"result": "clean", "final_time": 10.0},
            "windows": {
                "steady": {
                    "history": {
                        "analysis_window": {
                            "kinetic_mean": kinetic,
                            "abs_dp_mean": kinetic,
                        }
                    }
                }
            },
        },
    )
    return lineage, mhd, user


def write_direct_case(
    acceptance: Path,
    case_id: str,
    lineage: Path,
    mhd: Path,
    user: Path,
    *,
    result: str,
    health: str,
) -> Path:
    path = acceptance / "cases" / case_id / "case_acceptance.json"
    write_json(
        path,
        {
            "record_type": "cgl-lf-stage-i-direct-fast-case-acceptance",
            "case_id": case_id,
            "result": result,
            "health": {"result": health},
            "provenance": {
                "lineage": binding(lineage),
                "histories": {"mhd": binding(mhd), "user": binding(user)},
            },
        },
    )
    write_json(
        acceptance / "provenance.json",
        {
            "record_type": "cgl-lf-stage-i-direct-fast-acceptance-provenance",
            "outputs": {"case_acceptance": {case_id: binding(path)}},
        },
    )
    return path


def empty_data(publication, analysis: Path):
    return publication.PublicationData(
        analysis=analysis,
        cases={
            case_id: publication.CaseRecord(case_id=case_id)
            for case_id in publication.CASE_IDS
        },
        aggregate=None,
        comparisons=None,
        campaign_acceptance=None,
        science_record=None,
        ct_audit_record=None,
        acceptance_records=[],
        audit_records=[],
        source_paths=set(),
        ingestion_warnings=[],
    )


def test_stale_native_pass_cannot_override_current_direct_inconclusive(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    lineage, mhd, user = case_fixture(analysis, "R03")
    acceptance = tmp_path / "acceptance"
    write_direct_case(
        acceptance, "R03", lineage, mhd, user,
        result="inconclusive", health="inconclusive",
    )

    stale_mhd = tmp_path / "stale" / "old.mhd.hst"
    stale_user = tmp_path / "stale" / "old.user.hst"
    write_history(stale_mhd, 999.0)
    write_history(stale_user, 999.0)
    write_json(
        acceptance / "cases" / "R03" / "stale_native.json",
        seal({
            "record_type": "stage-i-scientific-case-evidence",
            "case_id": "R03",
            "result": "pass",
            "evaluation_inputs": {
                "mhd_history": binding(stale_mhd),
                "user_history": binding(stale_user),
            },
            "provenance": {"inputs": [binding(stale_mhd), binding(stale_user)]},
        }),
    )

    data = publication.discover_data(analysis, [acceptance])

    assert data.cases["R03"].acceptance is None
    assert data.cases["R03"].direct_acceptance is not None
    assert publication.acceptance_status(data, data.cases["R03"]) == "inconclusive"
    assert any("stale_native.json" in warning and "is stale" in warning
               for warning in data.ingestion_warnings)


def test_forged_evidence_digest_is_rejected(publication, tmp_path):
    analysis = tmp_path / "analysis"
    _, mhd, user = case_fixture(analysis, "R04")
    acceptance = tmp_path / "acceptance"
    forged = seal({
        "record_type": "stage-i-scientific-case-evidence",
        "case_id": "R04",
        "result": "pass",
        "evaluation_inputs": {
            "mhd_history": binding(mhd),
            "user_history": binding(user),
        },
        "provenance": {"inputs": [binding(mhd), binding(user)]},
    })
    forged["result"] = "fail"
    write_json(acceptance / "forged.json", forged)

    data = publication.discover_data(analysis, [acceptance])

    assert data.cases["R04"].acceptance is None
    assert any("evidence self-digest differs" in warning
               for warning in data.ingestion_warnings)


def test_self_sealed_record_without_required_provenance_is_rejected(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    _, mhd, user = case_fixture(analysis, "R05")
    acceptance = tmp_path / "acceptance"
    write_json(
        acceptance / "self_sealed_forgery.json",
        seal({
            "record_type": "stage-i-scientific-case-evidence",
            "case_id": "R05",
            "result": "pass",
            "evaluation_inputs": {
                "mhd_history": binding(mhd),
                "user_history": binding(user),
            },
            "provenance": {"inputs": [binding(mhd), binding(user)]},
        }),
    )

    data = publication.discover_data(analysis, [acceptance])

    assert data.cases["R05"].acceptance is None
    assert any("lacks required provenance binding criteria" in warning
               for warning in data.ingestion_warnings)


def test_real_hyperbolicity_snapshot_aggregate_schema_is_consumed(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    script = tmp_path / "audit.py"
    helper = tmp_path / "bin_convert.py"
    rank = analysis / "cases" / "R10" / "snapshots" / "rank_00000000" / "state.bin"
    for path, payload in ((script, "audit"), (helper, "helper"), (rank, "state")):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(payload, encoding="utf-8")
    rank_record = {
        **binding(rank),
        "mtime_ns": rank.stat().st_mtime_ns,
        "rank_id": 0,
    }
    audit = {
        "provenance": {
            "script_path": str(script),
            "script_sha256": sha256(script),
            "bin_convert_path": str(helper),
            "bin_convert_sha256": sha256(helper),
            "input_patterns": [str(analysis / "cases" / "R10")],
        },
        "snapshots": [{
            "rank_files": [rank_record],
            "input_inventory_sha256": hashlib.sha256(
                json.dumps(
                    [rank_record], sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
            ).hexdigest(),
            "aggregate": {
                "evaluated": 1000,
                "negative": 90,
                "negative_fraction": 0.09,
                "nonfinite_discriminant": 0,
                "minimum": -9.5,
            },
        }],
    }
    write_json(analysis / "audits" / "R10" / "hyperbolicity.json", audit)

    data = publication.discover_data(analysis, [])
    result = publication.hyperbolicity_diagnostics(data, "R10")

    assert len(data.audit_records) == 1
    assert result["result"] == "fail"
    assert result["negative_discriminant_fraction"] == pytest.approx(0.09)
    assert result["negative_discriminant_count"] == 90
    assert result["cell_direction_evaluations"] == 1000
    assert result["minimum_discriminant"] == pytest.approx(-9.5)


def test_inconclusive_or_failed_health_cannot_populate_response_products(
    publication, tmp_path
):
    data = empty_data(publication, tmp_path)
    for case_id in ("R02", "R14", "R16"):
        data.cases[case_id].diagnostics = {
            "health": {"result": "clean", "final_time": 10.0},
            "windows": {
                "steady": {
                    "history": {"analysis_window": {"kinetic_mean": 7.0, "abs_dp_mean": 7.0}},
                    "lf_history": {"applied_heat_flux_work": {"total": 8.0}},
                }
            },
        }
        data.cases[case_id].direct_acceptance = {
            "result": "inconclusive",
            "health": {"result": "inconclusive" if case_id != "R14" else "fail"},
        }
    data.cases["R04"].diagnostics = {
        "health": {"result": "clean", "final_time": 10.0},
        "windows": {"steady": {"history": {"analysis_window": {"kinetic_mean": 2.0}}}},
    }
    data.comparisons = {
        "resolution_detail": {
            "products": {"fixture": {"log_rms_R16_R02": 1.0, "log_rms_R02_R17": 2.0}}
        }
    }

    assert publication.metric_value(data.cases["R02"], "kinetic") == 7.0
    assert publication.response_metric_value(data, data.cases["R02"], "kinetic") is None
    robustness = publication.robustness_rows(data)
    row = next(
        value for value in robustness
        if value["contrast"] == "forcing A, beta=10" and value["metric"] == "kinetic"
    )
    assert row["reference_value"] is None
    limiter = next(
        value for value in publication.limiter_heat_flux_rows(data)
        if value["case_id"] == "R14" and value["family"] == "limiter"
    )
    assert limiter["abs_dp"] is None
    assert limiter["applied_heat_flux_work_abs"] is None
    resolution = publication.resolution_rows(data)
    assert next(
        value for value in resolution
        if value["record_type"] == "case_metric"
        and value["case_or_product"] == "R16"
        and value["metric"] == "kinetic"
    )["value"] is None
    assert not any(value["record_type"] == "fast_report_distance" for value in resolution)
    assert publication.publication_evidence_state(data) == "partial/transient"
    assert "Evidence state: **partial/transient**" in publication.report_markdown(
        data, [], tmp_path
    )


def test_r14_scope_prioritizes_hard_bound_and_never_invents_variant(
    publication, tmp_path
):
    data = empty_data(publication, tmp_path)
    data.cases["R10"].lineage = {"lineage_variants": []}
    data.cases["R14"].lineage = {}
    data.cases["R14"].diagnostics = {
        "health": {
            "result": "warnings",
            "final_time": 10.0,
            "hard_bound_diagnostic_maximum": 2.77e11,
        }
    }
    data.audit_records = [{
        "_publication_case_ids": ["R14"],
        "snapshots": [{
            "aggregate": {
                "evaluated": 300,
                "negative": 0,
                "negative_fraction": 0.0,
                "nonfinite_discriminant": 0,
                "minimum": 1.0,
            }
        }],
    }]

    observed, status = publication.scope_observed_diagnostic(data, "R14")
    rows = {row["case_id"]: row for row in publication.scope_rows(data)}

    assert publication.hyperbolicity_status(data, "R14") == "pass"
    assert observed == "hard-bound=2.77e+11"
    assert status == "warning"
    assert rows["R14"]["observed_diagnostic"] == observed
    assert rows["R14"]["variant"] is None
    assert rows["R10"]["variant"] is None
    assert "cannot be verified" in rows["R14"]["claim_scope"]


def test_authenticated_science_and_ct_are_integrated_but_non_authorizing(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    science = write_science_record(tmp_path / "science", analysis)
    ct = write_ct_record(tmp_path / "ct", analysis)

    data = publication.discover_data(analysis, [science, ct])
    contrast = publication.science_contrast_rows(data)[0]
    resolution = publication.science_resolution_rows(data)[0]
    mks24 = publication.science_mks24_rows(data)[0]
    ct_row = next(
        row for row in publication.ct_health_rows(data) if row["case_id"] == "R02"
    )

    assert data.science_record is not None
    assert data.ct_audit_record is not None
    assert contrast["standardized_effect"] == pytest.approx(-1.25)
    assert contrast["holm_significant"] is True
    assert contrast["release_authorizing"] is False
    assert resolution["available"] is False
    assert resolution["passed"] is None
    assert mks24["result"] == "inconclusive"
    assert ct_row["ct_result"] == "pass"
    assert ct_row["campaign_authority_eligible"] is False
    assert ct_row["release_authorizing"] is False
    assert publication.health_rows(data)[0]["direct_ct_numerical"] == "pass"
    report = publication.report_markdown(data, [], tmp_path)
    assert "non-authorizing-direct-fast-scientific-assessment" in report
    assert "campaign_authority_eligible=false" in report


def test_rendered_products_expose_authenticated_science_and_ct(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    science = write_science_record(tmp_path / "science", analysis)
    ct = write_ct_record(tmp_path / "ct", analysis)
    data = publication.discover_data(analysis, [science, ct])
    output = tmp_path / "publication"

    products = publication.render_products(data, output)

    assert output / "figures/fig07_reviewed_science_ct_summary.pdf" in products
    assert (output / "figures/fig07_reviewed_science_ct_summary.pdf").is_file()
    contrasts = (
        output / "tables/reviewed_science_contrasts.csv"
    ).read_text(encoding="utf-8")
    ct_health = (output / "tables/direct_ct_health.csv").read_text(encoding="utf-8")
    assert "standardized_effect" in contrasts
    assert "-1.25" in contrasts
    assert "R02,pass,true,true,true" in ct_health
    assert "release_authorizing" in ct_health


def test_forged_or_stale_science_and_ct_records_are_rejected(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    science = write_science_record(tmp_path / "science", analysis)
    science_record = json.loads(science.read_text(encoding="utf-8"))
    science_record["result"] = "pass"
    write_json(science, science_record)

    ct = write_ct_record(tmp_path / "ct", analysis)
    ct_record = json.loads(ct.read_text(encoding="utf-8"))
    stale_inventory = tmp_path / "stale_inventory.json"
    write_json(stale_inventory, {"record_type": "stale-inventory"})
    ct_record["inventory"] = binding(stale_inventory)
    write_json(ct, seal_ct({
        key: value for key, value in ct_record.items() if key != "evidence_digest"
    }))

    data = publication.discover_data(analysis, [science, ct])

    assert data.science_record is None
    assert data.ct_audit_record is None
    assert any(
        "science.json" in warning and "evidence self-digest differs" in warning
        for warning in data.ingestion_warnings
    )
    assert any(
        "ct_audit.json" in warning and "inventory is stale" in warning
        for warning in data.ingestion_warnings
    )


def test_r15_selected_nonfatal_variant_preserves_strict_failure_scope(
    publication, tmp_path
):
    data = empty_data(publication, tmp_path)
    data.cases["R15"].lineage = {
        "lineage_variants": ["finite_limiter_hard_bound_diagnostic_nonfatal"],
        "strict_failure": {
            "result": "fail",
            "time": 1.275643,
            "hard_bound": 2,
            "job_id": 4771183,
        },
    }
    lineage_path = tmp_path / "cases/R15/lineage.json"
    write_json(lineage_path, data.cases["R15"].lineage)
    data.cases["R15"].lineage_path = lineage_path
    data.cases["R15"].model = {"cgl_lf_strict_admissibility": "false"}
    data.cases["R15"].diagnostics = {
        "health": {
            "result": "warnings",
            "hard_bound_diagnostic_maximum": 2.0,
        }
    }
    data.science_record = {
        "_publication_evidence_validated": True,
        "result": "inconclusive",
        "case_dispositions": {
            "R15": {"claim_eligible": True, "acceptance_result": "pass"}
        },
        "provenance": {"case_lineages": {"R15": binding(lineage_path)}},
        "gates": [{
            "name": "finite_limiter_ordering:R15_gt_R14",
            "result": "pass",
            "reason": "diagnostic ordering",
        }],
    }

    scope = {row["case_id"]: row for row in publication.scope_rows(data)}["R15"]
    gate = publication.science_gate_rows(data)[0]

    assert publication.scope_status("R15") == "restricted"
    assert publication.science_case_status(data, "R15") == "restricted"
    assert scope["strict_run_disposition"] == (
        "fail; t=1.275643; hard_bound=2; job=4771183"
    )
    assert "never a strict-admissibility success" in scope["claim_scope"]
    assert "restricted nonfatal-hard-bound diagnostic" in gate["claim_scope"]
    assert publication.scope_observed_diagnostic(data, "R15") == (
        "hard-bound=2", "warning"
    )


def test_r15_strict_failure_details_are_never_hardcoded(publication, tmp_path):
    data = empty_data(publication, tmp_path)
    data.cases["R15"].lineage = {
        "lineage_variants": ["finite_limiter_hard_bound_diagnostic_nonfatal"],
        "strict_failure": {
            "result": "fail",
            "time": 1.275643,
            "hard_bound": 2,
            "job_id": 4771183,
        },
    }

    row = {value["case_id"]: value for value in publication.scope_rows(data)}["R15"]

    assert publication.r15_strict_failure_disposition(data) == "unavailable/inconclusive"
    assert row["strict_run_disposition"] == "unavailable/inconclusive"
    assert "1.275643" not in publication.report_markdown(data, [], tmp_path)
