"""Regression tests for audited direct-fast Stage I acceptance blockers."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
ADAPTER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_acceptance.py"
PUBLICATION = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_publication.py"
WORKFLOW = REPOSITORY / "scripts/cgl_lf_workflow.py"
MATRIX = REPOSITORY / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
HISTORY_LABEL = re.compile(r"\[(\d+)\]=(\S+)")


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fast_acceptance():
    return load_module("cgl_lf_stage_i_fast_acceptance_audit_regressions", ADAPTER)


@pytest.fixture(scope="module")
def publication():
    return load_module("cgl_lf_stage_i_fast_publication_audit_regressions", PUBLICATION)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_cases() -> dict[str, dict[str, object]]:
    manifest = json.loads(MATRIX.read_text(encoding="utf-8"))
    return {str(case["id"]): case for case in manifest["cases"]}


def expected_model_choices(input_path: Path, overrides: list[str]) -> dict[str, str]:
    workflow = load_module("cgl_lf_stage_i_fast_acceptance_test_workflow", WORKFLOW)
    source = input_path.read_text(encoding="utf-8")
    choices = {
        str(key): str(value)
        for key, value in workflow.model_choices(source, overrides).items()
    }
    strict = workflow.input_block_value(
        source, "mhd", "cgl_lf_strict_admissibility"
    ) or "false"
    for override in overrides:
        prefix = "mhd/cgl_lf_strict_admissibility="
        if override.startswith(prefix):
            strict = override[len(prefix):]
    choices["cgl_lf_strict_admissibility"] = strict
    return choices


def write_history(path: Path, columns: dict[str, list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = list(columns)
    rows = zip(*(columns[label] for label in labels))
    text = "# " + " ".join(
        f"[{index}]={label}" for index, label in enumerate(labels, start=1)
    ) + "\n"
    text += "\n".join(
        " ".join(format(value, ".17g") for value in row) for row in rows
    )
    path.write_text(text + "\n", encoding="utf-8")


def history_fixture(
    root: Path,
    *,
    counters: bool = True,
    final_time: float = 10.0,
) -> tuple[Path, Path]:
    times = [0.0, 4.0, 6.0, 8.0, final_time]
    mhd: dict[str, list[float]] = {"time": times}
    if counters:
        for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd"):
            mhd[name] = [0.0] * len(times)
    user = {
        "time": times,
        "kinetic": [1.0, 1.1, 1.2, 1.1, 1.0],
    }
    mhd_path = root / "fixture.mhd.hst"
    user_path = root / "fixture.user.hst"
    write_history(mhd_path, mhd)
    write_history(user_path, user)
    return mhd_path, user_path


def history_record(path: Path, *, declared_sha256: str | None = None) -> dict[str, object]:
    return {
        "available": True,
        "path": str(path.absolute()),
        "binding": {
            "path": str(path.absolute()),
            "size_bytes": path.stat().st_size,
            "sha256": declared_sha256 or sha256(path),
        },
    }


def lineage_fixture(
    root: Path,
    *,
    case_id: str = "R02",
    status: str = "complete",
    counters: bool = True,
    final_time: float = 10.0,
    strict: bool = True,
    variants: list[str] | None = None,
    overrides: list[str] | None = None,
) -> dict[str, object]:
    mhd, user = history_fixture(root, counters=counters, final_time=final_time)
    selected_overrides = overrides or []
    matrix_case = matrix_cases()[case_id]
    input_path = REPOSITORY / str(matrix_case["input"])
    model_choices = expected_model_choices(input_path, selected_overrides)
    model_choices["cgl_lf_strict_admissibility"] = str(strict).lower()
    return {
        "case_id": case_id,
        "case_name": matrix_case["name"],
        "matrix_case": matrix_case,
        "status": status,
        "target_time": 10.0,
        "errors": [],
        "lineage_variants": variants or [],
        "lineage_command_line_overrides": selected_overrides,
        "lineage_claim_scopes": ["standard"],
        "model_choices": model_choices,
        "input": {
            "path": str(input_path.resolve()),
            "size_bytes": input_path.stat().st_size,
            "sha256": sha256(input_path),
        },
        "lineage_identities": {
            "input_sha256": [sha256(input_path)],
            "matrix_sha256": [sha256(MATRIX)],
        },
        "histories": {
            "mhd": history_record(mhd),
            "user": history_record(user),
        },
    }


def policy_fixture(required: list[str] | None = None) -> dict[str, object]:
    return {
        "criteria_binding": {"path": "criteria.json", "sha256": "c" * 64},
        "review_binding": {"path": "review.json", "sha256": "d" * 64},
        "manifest": json.loads(MATRIX.read_text(encoding="utf-8")),
        "verified_sources": {
            "stage_i_manifest": {
                "path": str(MATRIX.resolve()),
                "size_bytes": MATRIX.stat().st_size,
                "sha256": sha256(MATRIX),
            },
        },
        "criteria": {
            "required_cases": required or ["R02"],
            "analysis_windows": {
                "full": [4.0, 10.0],
                "early": [4.0, 8.0],
                "late": [6.0, 10.0],
            },
            "statistics_policy": {
                "bootstrap_replicates": 8,
                "gap_policy": {
                    "expected_history_cadence": 0.25,
                    "maximum_gap_expected_cadence_multiplier": 8.0,
                    "maximum_gap_forcing_tcorr_fraction": 0.25,
                },
            },
            "case_metrics": {
                "kinetic": {
                    "history": "user",
                    "column": "kinetic",
                    "stationarity_kind": "energy",
                },
            },
            "family_gates": {
                "active_passive": {"pairs": []},
            },
        },
    }


class FakeAcceptanceError(RuntimeError):
    """Fixture equivalent of the reviewed utility's input error."""


class FakeAcceptance:
    """Small reviewed-acceptance API stub for adapter contract tests."""

    AcceptanceError = FakeAcceptanceError

    def __init__(self, policy: dict[str, object]):
        self.policy = policy

    def load_validated_policy(self, _criteria: Path, _review: Path) -> dict[str, object]:
        return self.policy

    @staticmethod
    def case_name(policy: dict[str, object], case_id: str) -> str:
        for case in policy["manifest"]["cases"]:
            if case["id"] == case_id:
                return str(case["name"])
        raise FakeAcceptanceError(f"missing fixture case: {case_id}")

    @staticmethod
    def canonical_bundle_path(case_id: str) -> Path:
        return Path(f"/nonexistent/canonical/{case_id}/manifest.json")

    def load_history(
        self, path: Path, _label: str
    ) -> tuple[dict[str, list[float]], dict[str, object]]:
        payload = path.read_text(encoding="utf-8")
        labels: list[str] | None = None
        rows: list[list[float]] = []
        for line in payload.splitlines():
            if line.startswith("#"):
                found = HISTORY_LABEL.findall(line)
                if found:
                    labels = [
                        name for _, name in sorted(found, key=lambda item: int(item[0]))
                    ]
            elif line.strip():
                rows.append([float(value) for value in line.split()])
        assert labels is not None
        return (
            {
                label: [row[index] for row in rows]
                for index, label in enumerate(labels)
            },
            {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            },
        )

    @staticmethod
    def _forcing_tcorr_from_bundle(value: object) -> float:
        if isinstance(value, dict):
            candidate = value.get("forcing_tcorr")
            return float(candidate) if candidate is not None else 0.0
        if isinstance(value, (str, Path)) and Path(value).is_file():
            record = json.loads(Path(value).read_text(encoding="utf-8"))
            for key in ("forcing_tcorr", "model_choices"):
                candidate = record.get(key)
                if isinstance(candidate, dict):
                    candidate = candidate.get("forcing_tcorr")
                if candidate is not None:
                    return float(candidate)
        return 0.0

    def evaluate_case(self, *args, **kwargs) -> dict[str, object]:
        case_id = str(args[1])
        bundle = kwargs.get("bundle_path", args[6] if len(args) > 6 else None)
        tcorr = self._forcing_tcorr_from_bundle(bundle)
        statistics = self._statistics(tcorr)
        return {
            "record_type": "stage-i-scientific-case-evidence",
            "case_id": case_id,
            "result": "pass",
            "metrics": {
                "kinetic": {
                    "full": statistics,
                    "early": statistics,
                    "late": statistics,
                    "sampling_adequacy": "pass",
                    "stationarity": {"result": "pass"},
                },
            },
        }

    @staticmethod
    def _statistics(minimum_block_duration: float) -> dict[str, object]:
        return {
            "mean": 1.0,
            "standard_deviation": 0.1,
            "standard_error": 0.05,
            "confidence_interval_95": [0.9, 1.1],
            "effective_sample_count": 4.0,
            "independent_time_block_count": 2,
            "gap_adequacy": "pass",
            "sample_count": 5,
            "method": {
                "required_block_duration": minimum_block_duration,
                "minimum_block_duration": minimum_block_duration,
            },
        }

    def window_statistics(
        self,
        _times: list[float],
        _values: list[float],
        _start: float,
        _end: float,
        *,
        minimum_block_duration: float,
        **_kwargs,
    ) -> dict[str, object]:
        return self._statistics(minimum_block_duration)

    def metric_statistics(
        self,
        _history: dict[str, list[float]],
        _values: list[float],
        _metric: str,
        _policy: dict[str, object],
        *,
        kind: str,
        minimum_block_duration: float,
    ) -> dict[str, object]:
        statistics = self._statistics(minimum_block_duration)
        return {
            "full": statistics,
            "early": statistics,
            "late": statistics,
            "sampling_adequacy": "pass",
            "stationarity": {"result": "pass", "kind": kind},
        }

    @staticmethod
    def gate(
        name: str,
        result: str,
        *,
        reason: str,
        observations: object = None,
    ) -> dict[str, object]:
        return {
            "name": name,
            "result": result,
            "reason": reason,
            "observations": observations,
        }

    @staticmethod
    def aggregate_gate_result(gates: list[dict[str, object]]) -> str:
        results = [str(gate["result"]) for gate in gates]
        if "fail" in results:
            return "fail"
        if "inconclusive" in results:
            return "inconclusive"
        return "pass"

    @staticmethod
    def seal_evidence(value: dict[str, object]) -> dict[str, object]:
        return value

    @staticmethod
    def pair_contrast(
        _active: dict[str, object],
        _passive: dict[str, object],
        _policy: dict[str, object],
    ) -> dict[str, object]:
        return {"result": "pass", "reason": "fixture pair passed", "metrics": {}}

    @staticmethod
    def scalar_from_case(
        _case: dict[str, object], _metric: str, _window: str
    ) -> dict[str, float]:
        return {"mean": 2.0, "standard_error": 0.01, "standard_deviation": 0.1}


def summarize(
    fast_acceptance,
    tmp_path: Path,
    case_id: str,
    lineage: dict[str, object],
    *,
    policy: dict[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, object] | None]:
    selected_policy = policy or policy_fixture([case_id])
    return fast_acceptance.summarize_case(
        FakeAcceptance(selected_policy),
        selected_policy,
        case_id,
        lineage,
        {"path": "fixture-lineage.json", "sha256": "a" * 64},
        "fixture",
        tmp_path / "acceptance" / "cases" / case_id,
    )


def test_wrong_case_identity_is_rejected(fast_acceptance, tmp_path):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R02")

    summary, reviewed = summarize(fast_acceptance, tmp_path, "R17", lineage)

    assert summary["result"] != "pass"
    assert summary["health"]["result"] == "fail"
    assert any(
        "case_id differs" in error
        for error in summary["health"]["structural_errors"]
    )
    assert reviewed is None


@pytest.mark.parametrize(
    ("case_id", "strict", "variants", "overrides"),
    [
        ("R03", False, ["standard"], ["mhd/cgl_lf_strict_admissibility=false"]),
        (
            "R14",
            True,
            ["finite_limiter_hard_bound_diagnostic_nonfatal"],
            [
                "mhd/cgl_lf_strict_admissibility=false",
                "mhd/cgl_lf_strict_admissibility=true",
            ],
        ),
    ],
)
def test_strict_policy_and_configuration_mismatches_are_ineligible(
    fast_acceptance,
    tmp_path,
    case_id,
    strict,
    variants,
    overrides,
):
    lineage = lineage_fixture(
        tmp_path / case_id,
        case_id=case_id,
        strict=strict,
        variants=variants,
        overrides=overrides,
    )

    errors = fast_acceptance.lineage_identity_errors(
        FakeAcceptance(policy_fixture([case_id])),
        policy_fixture([case_id]),
        case_id,
        lineage,
        {"path": "fixture-lineage.json", "sha256": "a" * 64},
    )

    assert errors
    assert any(
        "strict-admissibility" in error or "override" in error
        for error in errors
    )
    if case_id == "R14":
        scope = fast_acceptance.scientific_scope(case_id, lineage)
        assert scope["campaign_interpretation_eligible"] is False
        assert scope["uniform_strict_diagnostics_eligible"] is False


def test_missing_strict_model_choice_legacy_exception_is_r02_only(
    fast_acceptance, tmp_path
):
    policy = policy_fixture(["R02", "R03"])
    acceptance = FakeAcceptance(policy)
    r02 = lineage_fixture(tmp_path / "R02", case_id="R02")
    r03 = lineage_fixture(tmp_path / "R03", case_id="R03")
    del r02["model_choices"]["cgl_lf_strict_admissibility"]
    del r03["model_choices"]["cgl_lf_strict_admissibility"]

    r02_errors = fast_acceptance.lineage_identity_errors(
        acceptance, policy, "R02", r02, {"path": "R02-lineage.json"}
    )
    r03_errors = fast_acceptance.lineage_identity_errors(
        acceptance, policy, "R03", r03, {"path": "R03-lineage.json"}
    )

    assert not r02_errors
    assert any("model choices" in error for error in r03_errors)


def test_spoofed_model_choices_and_forcing_tcorr_are_rejected(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R03")
    lineage["model_choices"]["forcing_tcorr"] = "0.000001"
    lineage["model_choices"]["beta0"] = "999"

    summary, comparison = summarize(fast_acceptance, tmp_path, "R03", lineage)

    assert summary["health"]["result"] == "fail"
    assert comparison is None
    assert any(
        "model choices differ" in error
        for error in summary["health"]["structural_errors"]
    )


def test_bound_input_must_match_matrix_authoritative_input(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", case_id="R03")
    original = Path(str(lineage["input"]["path"]))
    attacker_input = tmp_path / "attacker.athinput"
    attacker_input.write_text(
        original.read_text(encoding="utf-8").replace(
            "tcorr = 2.0", "tcorr = 0.000001"
        ),
        encoding="utf-8",
    )
    lineage["input"] = {
        "path": str(attacker_input.resolve()),
        "size_bytes": attacker_input.stat().st_size,
        "sha256": sha256(attacker_input),
    }
    lineage["lineage_identities"]["input_sha256"] = [sha256(attacker_input)]
    lineage["model_choices"] = expected_model_choices(attacker_input, [])

    summary, comparison = summarize(fast_acceptance, tmp_path, "R03", lineage)

    assert summary["health"]["result"] == "fail"
    assert comparison is None
    assert any(
        "matrix-authoritative input" in error
        for error in summary["health"]["structural_errors"]
    )


def test_declared_history_binding_mismatch_is_rejected(fast_acceptance, tmp_path):
    lineage = lineage_fixture(tmp_path / "lineage")
    lineage["histories"]["mhd"]["binding"]["sha256"] = "0" * 64
    reviewed = FakeAcceptance(policy_fixture())

    histories, _bindings, errors = fast_acceptance.load_histories(reviewed, lineage)

    assert "mhd" not in histories
    assert any("sha256" in error or "binding" in error for error in errors)


def test_failed_complete_cases_are_excluded_from_reviewed_comparisons(
    fast_acceptance, tmp_path
):
    policy = policy_fixture(["R14", "R15"])
    acceptance = FakeAcceptance(policy)
    lineage = lineage_fixture(
        tmp_path / "R14",
        case_id="R14",
        strict=False,
        variants=["finite_limiter_hard_bound_diagnostic_nonfatal"],
        overrides=["mhd/cgl_lf_strict_admissibility=false"],
    )
    histories, _bindings, errors = fast_acceptance.load_histories(acceptance, lineage)
    scope = fast_acceptance.scientific_scope("R14", lineage)
    failed_health = fast_acceptance.numerical_health(
        "R14",
        lineage,
        scope,
        histories,
        errors,
        ["fixture provenance failure"],
    )
    reviewed, _reviewed_reason = fast_acceptance.reviewed_complete_case_evidence(
        acceptance, policy, "R14", lineage, failed_health
    )
    comparison, _comparison_reason = fast_acceptance.comparison_case_evidence(
        acceptance,
        policy,
        "R14",
        lineage,
        scope,
        failed_health,
        histories,
    )
    failed_summary = {
        "result": "fail",
        "health": failed_health,
        "scope": scope,
    }
    passed_summary = {
        "result": "pass",
        "health": {"result": "pass", "complete_to_target": True},
        "scope": {
            "classification": "standard_claim_scope",
            "campaign_interpretation_eligible": True,
        },
    }
    assert reviewed is None
    assert comparison is None

    campaign = fast_acceptance.build_campaign_evidence(
        acceptance,
        policy,
        {"R14": failed_summary, "R15": passed_summary},
        {"R15": {"result": "pass"}},
    )
    gate_by_name = {gate["name"]: gate for gate in campaign["gates"]}

    reviewed_gate = gate_by_name["scoped_exact_window_comparison_evidence"]
    limiter_gate = gate_by_name["finite_limiter_ordering:R15_gt_R14"]
    assert "R14" not in reviewed_gate["observations"]["available_cases"]
    assert limiter_gate["result"] == "inconclusive"


def test_comparison_gates_require_pass_case_comparison_evidence(fast_acceptance):
    policy = policy_fixture(["R02", "R06", "R14", "R15"])
    policy["criteria"]["family_gates"]["active_passive"]["pairs"] = [["R02", "R06"]]
    acceptance = FakeAcceptance(policy)
    standard_scope = {
        "classification": "standard_claim_scope",
        "campaign_interpretation_eligible": True,
    }
    r14_scope = {
        "classification": "scoped_nonfatal_hard_bound_variant",
        "campaign_interpretation_eligible": True,
    }
    summaries = {
        case_id: {
            "result": "inconclusive",
            "health": {"result": "pass", "complete_to_target": True},
            "scope": r14_scope if case_id == "R14" else standard_scope,
        }
        for case_id in ("R02", "R06", "R14", "R15")
    }
    comparison_cases = {
        case_id: {"case_id": case_id, "result": "inconclusive"}
        for case_id in summaries
    }

    campaign = fast_acceptance.build_campaign_evidence(
        acceptance, policy, summaries, comparison_cases
    )
    gates = {gate["name"]: gate for gate in campaign["gates"]}

    assert gates["active_passive_pair:R02:R06"]["result"] == "inconclusive"
    assert gates["finite_limiter_ordering:R15_gt_R14"]["result"] == "inconclusive"


def test_complete_case_statistics_use_authenticated_forcing_tcorr(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage")

    summary, comparison = summarize(fast_acceptance, tmp_path, "R02", lineage)

    assert comparison is not None
    assert comparison["minimum_block_duration"] == pytest.approx(2.0)
    statistics = summary["history_statistics"]["kinetic"]["windows"]["full"][
        "statistics"
    ]
    assert statistics["method"]["required_block_duration"] == pytest.approx(2.0)


def test_rerun_removes_stale_reviewed_case_evidence(fast_acceptance, tmp_path):
    output = tmp_path / "acceptance" / "cases" / "R02"
    stale = output / "reviewed_case_evidence.json"
    write_json(stale, {"record_type": "stage-i-scientific-case-evidence"})
    assert stale.is_file()
    policy = policy_fixture()
    reviewed = FakeAcceptance(policy)
    partial = lineage_fixture(
        tmp_path / "partial",
        status="partial",
        final_time=2.0,
    )

    fast_acceptance.summarize_case(
        reviewed,
        policy,
        "R02",
        partial,
        {"path": "fixture-lineage.json", "sha256": "a" * 64},
        "fixture",
        output,
    )

    assert not stale.exists()


def test_missing_histories_are_inconclusive_for_in_progress_case(
    fast_acceptance,
):
    scope = {
        "hard_bound_is_fatal": True,
        "campaign_interpretation_eligible": True,
    }
    lineage = {"status": "partial", "target_time": 10.0, "errors": []}

    health = fast_acceptance.numerical_health(
        "R02",
        lineage,
        scope,
        {},
        [
            "evidence gap: assembled mhd history is unavailable",
            "evidence gap: assembled user history is unavailable",
        ],
        [],
    )

    assert health["result"] == "inconclusive"


def test_missing_fatal_counters_are_inconclusive_for_complete_case(
    fast_acceptance, tmp_path
):
    lineage = lineage_fixture(tmp_path / "lineage", counters=False)
    reviewed = FakeAcceptance(policy_fixture())
    histories, _bindings, errors = fast_acceptance.load_histories(reviewed, lineage)
    scope = fast_acceptance.scientific_scope("R02", lineage)

    health = fast_acceptance.numerical_health(
        "R02", lineage, scope, histories, errors, []
    )

    assert health["result"] == "inconclusive"
    assert all(
        health["fatal_counter_maxima"][name] is None
        for name in fast_acceptance.FATAL_COUNTERS
    )


def test_symlink_output_cannot_write_inside_fast_report_source(
    fast_acceptance, tmp_path, monkeypatch
):
    source = tmp_path / "fast-report"
    source.mkdir()
    inventory = source / "inventory.json"
    lineage = lineage_fixture(source / "histories")
    write_json(
        inventory,
        {
            "output": str(source),
            "cases": {"R02": lineage},
        },
    )
    output_link = tmp_path / "acceptance-output"
    output_link.symlink_to(source / "nested-acceptance", target_is_directory=True)
    policy = policy_fixture()
    monkeypatch.setattr(
        fast_acceptance,
        "load_acceptance_module",
        lambda: FakeAcceptance(policy),
    )

    with pytest.raises(SystemExit):
        fast_acceptance.main(
            [
                "--inventory",
                str(inventory),
                "--output",
                str(output_link),
                "--cases",
                "R02",
            ]
        )

    assert not (source / "nested-acceptance").exists()


def test_nested_symlink_output_escape_is_rejected_before_writing(
    fast_acceptance, tmp_path, monkeypatch
):
    source = tmp_path / "fast-report"
    source.mkdir()
    inventory = source / "inventory.json"
    lineage = lineage_fixture(source / "histories")
    write_json(inventory, {"output": str(source), "cases": {"R02": lineage}})
    escaped = source / "escaped-acceptance"
    escaped.mkdir()
    output = tmp_path / "acceptance-output"
    (output / "cases").mkdir(parents=True)
    (output / "cases" / "R02").symlink_to(escaped, target_is_directory=True)
    policy = policy_fixture()
    monkeypatch.setattr(
        fast_acceptance,
        "load_acceptance_module",
        lambda: FakeAcceptance(policy),
    )

    with pytest.raises(SystemExit):
        fast_acceptance.main(
            [
                "--inventory",
                str(inventory),
                "--output",
                str(output),
                "--cases",
                "R02",
            ]
        )

    assert not (escaped / "case_acceptance.json").exists()


def test_publication_discovers_scoped_direct_fast_campaign_record(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    acceptance = tmp_path / "acceptance"
    write_json(
        acceptance / "campaign_evidence.json",
        {
            "record_type": "cgl-lf-stage-i-direct-fast-campaign-evidence",
            "result": "inconclusive",
            "case_results": {"R10": "inconclusive", "R14": "pass"},
            "claim_scope": {
                "R10": {"classification": "exploratory_only"},
                "R14": {
                    "classification": "scoped_nonfatal_hard_bound_variant"
                },
            },
            "gates": [
                {
                    "name": "explicit_claim_scope",
                    "result": "pass",
                    "reason": "fixture",
                }
            ],
        },
    )

    data = publication.discover_data(analysis, [acceptance])

    assert data.campaign_acceptance is not None
    assert (
        data.campaign_acceptance["record_type"]
        == "cgl-lf-stage-i-direct-fast-campaign-evidence"
    )
    assert any(
        record.get("record_type")
        == "cgl-lf-stage-i-direct-fast-campaign-evidence"
        for record in data.acceptance_records
    )
