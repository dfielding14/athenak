"""Regression tests for audited direct-fast Stage I publication blockers."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import threading

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


def test_publication_facing_language_is_descriptive():
    source = RENDERER.read_text(encoding="utf-8")
    for prohibited in (
        "Authenticated Holm-corrected comparison",
        "Developed-window robustness summary",
        "Matched active/passive causal-mechanism evidence",
        "Common-scale convergence evidence",
        "Reviewed convergence result:",
        "Populated paired developed-window robustness cells",
    ):
        assert prohibited not in source
    for required in (
        "Authenticated descriptive comparison",
        "Developed-window measured-sensitivity summary",
        "Matched active/passive mechanism diagnostics",
        "Common-scale consistency evidence",
        "Reviewed common-scale consistency result:",
        "Populated paired developed-window sensitivity cells",
    ):
        assert required in source


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


def write_science_record(root: Path, analysis: Path, publication_module) -> Path:
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
        "current_science_scope_limitation": (
            publication_module.CURRENT_SCIENCE_SCOPE_LIMITATION
        ),
        "active_passive_intervention_scope": (
            publication_module.ACTIVE_PASSIVE_INTERVENTION_SCOPE
        ),
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
                        "standardized_effect": -1.25,
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
            "observations": [{"metric": "kinetic", "descriptive_direction": "left_lt_right"}],
            "limits": {"inference_scope": "descriptive_within_realization"},
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


def rewrite_science_record(path: Path, record: dict[str, object]) -> None:
    body = {key: value for key, value in record.items() if key != "evidence_digest"}
    write_json(path, seal(body))
    external_path = path.parent / "provenance.json"
    external = json.loads(external_path.read_text(encoding="utf-8"))
    external["inputs"] = body["provenance"]
    external["outputs"]["science"] = binding(path)
    write_json(external_path, external)


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
                    "result": "pass",
                    "numerical_result": "pass",
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


def rewrite_ct_record(path: Path, record: dict[str, object]) -> None:
    body = {key: value for key, value in record.items() if key != "evidence_digest"}
    write_json(path, seal_ct(body))


def strict_failure_record(
    root: Path, case_id: str = "R15", *, failure_time: float = 1.275643,
    hard_bound: int = 2, job_id: str = "4771183",
) -> dict[str, object]:
    manifest = root / "fast_run.json"
    exit_code = root / "run_exit_code"
    slurm_log = root / "strict.log"
    write_json(manifest, {"case_id": case_id, "job_id": job_id})
    exit_code.write_text("1\n", encoding="utf-8")
    slurm_log.write_text("strict hard-bound failure\n", encoding="utf-8")
    return {
        "schema_version": 1,
        "record_type": "cgl-lf-stage-i-retained-strict-failure-evidence",
        "case_id": case_id,
        "job_id": job_id,
        "result": "fail",
        "failure_time": failure_time,
        "failure_counters": {
            "lf_dfloor": 0,
            "lf_pfloor": 0,
            "lf_nonfin": 0,
            "lf_nonpos": 0,
            "lf_hardbd": hard_bound,
        },
        "strict_admissibility_evidence": True,
        "provenance": {
            "manifest": binding(manifest),
            "run_exit_code": binding(exit_code),
            "slurm_log": binding(slurm_log),
        },
    }


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


def tree_bytes(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def staged_publication(
    publication,
    staging: Path,
    output: Path,
    publisher: str,
    values: dict[str, str],
) -> tuple[list[Path], dict[str, object]]:
    staging.mkdir()
    products = []
    for relative, value in values.items():
        path = staging / relative
        publication.write_text(path, value)
        products.append(path)
    return products, {
        "record_type": "cgl_lf_stage_i_fast_publication_products",
        "publisher": publisher,
        "analysis_output": str((output.parent / "analysis").absolute()),
        "normalized_invocation": ["renderer", "--output", str(output.absolute())],
        "sources": [],
        "products": publication.canonical_product_bindings(
            products, staging, output
        ),
    }


def test_staging_failure_preserves_previous_canonical_authority(
    publication, tmp_path, monkeypatch
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    output = tmp_path / "publication"
    publication.write_text(output / "tables/old.csv", "old\n")
    publication.write_json(output / "manifest.json", {"authority": "old"})
    previous = tree_bytes(output)

    monkeypatch.setattr(
        publication,
        "discover_data",
        lambda analysis_path, _acceptance: empty_data(publication, analysis_path),
    )

    def interrupted_render(_data, staging):
        assert staging.parent == output.parent
        assert staging != output
        assert staging.name.startswith(".publication.staging-")
        publication.write_text(staging / "tables/new.csv", "new\n")
        raise RuntimeError("interrupted render")

    monkeypatch.setattr(publication, "render_products", interrupted_render)

    with pytest.raises(RuntimeError, match="interrupted render"):
        publication.main([str(analysis), "--output", str(output)])

    assert tree_bytes(output) == previous
    assert list(tmp_path.glob(".publication.staging-*")) == []


def test_interrupted_promotion_withdraws_manifest_and_rerun_recovers(
    publication, tmp_path, monkeypatch
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    output = tmp_path / "publication"

    monkeypatch.setattr(
        publication,
        "discover_data",
        lambda analysis_path, _acceptance: empty_data(publication, analysis_path),
    )

    def fixed_render(_data, staging):
        first = staging / "a.txt"
        second = staging / "b.txt"
        publication.write_text(first, "new a\n")
        publication.write_text(second, "new b\n")
        return [first, second]

    monkeypatch.setattr(publication, "render_products", fixed_render)
    real_promote_file = publication.promote_file
    product_promotions = 0

    def interrupt_after_first_product(staged, relative, output_descriptor):
        nonlocal product_promotions
        real_promote_file(staged, relative, output_descriptor)
        if relative.name != "manifest.json":
            product_promotions += 1
            if product_promotions == 1:
                assert not (output / "manifest.json").exists()
                raise RuntimeError("interrupted promotion")

    monkeypatch.setattr(publication, "promote_file", interrupt_after_first_product)

    with pytest.raises(RuntimeError, match="interrupted promotion"):
        publication.main([str(analysis), "--output", str(output)])

    ownership = publication.publication_ownership_path(output)
    assert publication.valid_publication_ownership(output)
    assert ownership.is_file()
    assert not (output / "manifest.json").exists()
    assert (output / "a.txt").read_text(encoding="utf-8") == "new a\n"
    assert not (output / "b.txt").exists()
    assert list(tmp_path.glob(".publication.staging-*")) == []

    monkeypatch.setattr(publication, "promote_file", real_promote_file)
    publication.main([str(analysis), "--output", str(output)])

    recovered_manifest = json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    )
    assert (
        recovered_manifest["record_type"]
        == "cgl_lf_stage_i_fast_publication_products"
    )
    assert (output / "a.txt").read_text(encoding="utf-8") == "new a\n"
    assert (output / "b.txt").read_text(encoding="utf-8") == "new b\n"
    for product in recovered_manifest["products"]:
        path = Path(product["path"])
        assert path.is_relative_to(output.absolute())
        assert publication.source_binding(path) == product


def test_post_manifest_failure_withdraws_authenticated_authority(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    staging = tmp_path / ".publication.staging-post-manifest-failure"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "post-manifest-failure",
        {"a.txt": "new a\n", "b.txt": "new b\n"},
    )
    real_promote_file = publication.promote_file

    def corrupt_product_after_manifest(staged, relative, output_descriptor):
        real_promote_file(staged, relative, output_descriptor)
        if relative == Path("manifest.json"):
            publication.write_text(output / "a.txt", "tampered after manifest\n")

    monkeypatch.setattr(
        publication, "promote_file", corrupt_product_after_manifest
    )

    with pytest.raises(
        publication.PublicationError, match="differs from staged binding"
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert publication.valid_publication_ownership(output)
    assert not (output / "manifest.json").exists()
    assert (output / "a.txt").read_text(encoding="utf-8") == (
        "tampered after manifest\n"
    )
    assert (output / "b.txt").read_text(encoding="utf-8") == "new b\n"


def test_post_manifest_failure_preserves_unrecognized_replacement(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    staging = tmp_path / ".publication.staging-unrecognized-manifest"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "unrecognized-manifest",
        {"a.txt": "new a\n"},
    )
    real_promote_file = publication.promote_file
    replacement = b'{"authority":"unrecognized replacement"}\n'

    def replace_manifest_then_fail(staged, relative, output_descriptor):
        real_promote_file(staged, relative, output_descriptor)
        if relative == Path("manifest.json"):
            (output / "manifest.json").write_bytes(replacement)
            raise RuntimeError("failure after unrecognized manifest replacement")

    monkeypatch.setattr(publication, "promote_file", replace_manifest_then_fail)

    with pytest.raises(
        RuntimeError, match="failure after unrecognized manifest replacement"
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert (output / "manifest.json").read_bytes() == replacement
    assert publication.valid_publication_ownership(output)


def test_authenticated_withdrawal_preserves_replacement_after_quarantine_rename(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    output.mkdir()
    expected = b'{"authority":"authenticated"}\n'
    replacement = b'{"authority":"replacement-after-quarantine"}\n'
    (output / "manifest.json").write_bytes(expected)
    output_descriptor = publication.open_output_directory(output)
    real_read_regular_file_at = publication.read_regular_file_at
    replacement_installed = False

    def replace_after_quarantine(output_descriptor_value, relative):
        nonlocal replacement_installed
        payload = real_read_regular_file_at(output_descriptor_value, relative)
        if (
            output_descriptor_value != output_descriptor
            and relative == Path("manifest.json")
            and not replacement_installed
        ):
            staged_replacement = output / ".replacement-manifest"
            staged_replacement.write_bytes(replacement)
            os.replace(staged_replacement, output / "manifest.json")
            replacement_installed = True
        return payload

    monkeypatch.setattr(
        publication, "read_regular_file_at", replace_after_quarantine
    )
    try:
        assert publication.withdraw_authenticated_canonical_manifest(
            output_descriptor, expected
        )
    finally:
        os.close(output_descriptor)

    assert replacement_installed
    assert (output / "manifest.json").read_bytes() == replacement
    quarantines = list(
        tmp_path.glob(
            f".publication{publication.PUBLICATION_MANIFEST_QUARANTINE_PREFIX}*"
        )
    )
    assert len(quarantines) == 1
    assert (quarantines[0] / "manifest.json").read_bytes() == expected


def test_quarantine_name_collision_preserves_unknown_directory(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    output.mkdir()
    expected = b'{"authority":"authenticated"}\n'
    (output / "manifest.json").write_bytes(expected)
    existing = (
        tmp_path
        / f".publication{publication.PUBLICATION_MANIFEST_QUARANTINE_PREFIX}collision"
    )
    existing.mkdir()
    (existing / "unknown.txt").write_text("preserve\n", encoding="utf-8")
    tokens = iter(("collision", "unique"))
    monkeypatch.setattr(publication.secrets, "token_hex", lambda _size: next(tokens))
    output_descriptor = publication.open_output_directory(output)

    try:
        assert publication.withdraw_authenticated_canonical_manifest(
            output_descriptor, expected
        )
    finally:
        os.close(output_descriptor)

    assert (existing / "unknown.txt").read_text(encoding="utf-8") == "preserve\n"
    assert not (output / "manifest.json").exists()
    unique = (
        tmp_path
        / f".publication{publication.PUBLICATION_MANIFEST_QUARANTINE_PREFIX}unique"
    )
    assert (unique / "manifest.json").read_bytes() == expected


def test_authenticated_withdrawal_never_deletes_quarantine_replacement(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    output.mkdir()
    expected = b'{"authority":"authenticated"}\n'
    replacement = b'{"authority":"unknown-quarantine-replacement"}\n'
    (output / "manifest.json").write_bytes(expected)
    output_descriptor = publication.open_output_directory(output)
    real_read_regular_file_at = publication.read_regular_file_at
    replacement_installed = False

    def replace_quarantined_manifest(descriptor, relative):
        nonlocal replacement_installed
        payload = real_read_regular_file_at(descriptor, relative)
        if (
            descriptor != output_descriptor
            and relative == Path("manifest.json")
            and not replacement_installed
        ):
            staged = ".unknown-replacement"
            flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
            replacement_descriptor = os.open(staged, flags, 0o600, dir_fd=descriptor)
            try:
                os.write(replacement_descriptor, replacement)
            finally:
                os.close(replacement_descriptor)
            os.replace(
                staged,
                "manifest.json",
                src_dir_fd=descriptor,
                dst_dir_fd=descriptor,
            )
            replacement_installed = True
        return payload

    monkeypatch.setattr(
        publication, "read_regular_file_at", replace_quarantined_manifest
    )
    try:
        assert publication.withdraw_authenticated_canonical_manifest(
            output_descriptor, expected
        )
    finally:
        os.close(output_descriptor)

    quarantines = list(
        tmp_path.glob(
            f".publication{publication.PUBLICATION_MANIFEST_QUARANTINE_PREFIX}*"
        )
    )
    assert replacement_installed
    assert len(quarantines) == 1
    assert (quarantines[0] / "manifest.json").read_bytes() == replacement
    assert not (output / "manifest.json").exists()


def test_all_staged_products_are_fsynced_before_first_promotion(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    staging = tmp_path / ".publication.staging-fsync-products"
    staging.mkdir()
    text = staging / "tables/result.txt"
    publication.write_text(text, "result\n")
    pdf = staging / "figures/result.pdf"
    pdf.parent.mkdir(parents=True)
    pdf.write_bytes(b"%PDF-1.4\nfixture\n")
    products = [text, pdf]
    manifest = {
        "record_type": "cgl_lf_stage_i_fast_publication_products",
        "publisher": "fsync-products",
        "normalized_invocation": ["renderer", "--output", str(output.absolute())],
        "products": publication.canonical_product_bindings(
            products, staging, output
        ),
    }
    real_fsync_regular_file = publication.fsync_regular_file
    real_promote_file = publication.promote_file
    synced: list[Path] = []

    def observed_fsync(path):
        real_fsync_regular_file(path)
        synced.append(path.absolute())

    def observed_promote(staged, relative, output_descriptor):
        expected = {
            text.absolute(),
            pdf.absolute(),
            (staging / "manifest.json").absolute(),
        }
        assert expected <= set(synced)
        real_promote_file(staged, relative, output_descriptor)

    monkeypatch.setattr(publication, "fsync_regular_file", observed_fsync)
    monkeypatch.setattr(publication, "promote_file", observed_promote)

    publication.promote_staged_publication(
        staging,
        output,
        products,
        manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )

    assert pdf.absolute() in synced
    assert json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    ) == manifest


def test_semantically_equivalent_unfsynced_manifest_replacement_is_rejected(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    staging = tmp_path / ".publication.staging-manifest-exact-bytes"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "manifest-exact-bytes",
        {"a.txt": "new a\n"},
    )
    real_fsync_regular_file = publication.fsync_regular_file
    replaced = False

    def replace_after_fsync(path):
        nonlocal replaced
        real_fsync_regular_file(path)
        if path == staging / "manifest.json" and not replaced:
            value = json.loads(path.read_text(encoding="utf-8"))
            path.write_text(json.dumps(value, separators=(",", ":")), encoding="utf-8")
            replaced = True

    monkeypatch.setattr(publication, "fsync_regular_file", replace_after_fsync)

    with pytest.raises(
        publication.PublicationError,
        match="staged publication manifest differs from exact expected bytes",
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert replaced
    assert not (output / "manifest.json").exists()


def test_staged_pdf_fsync_failure_preserves_previous_canonical_authority(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    initial_staging = tmp_path / ".publication.staging-fsync-initial"
    initial_products, initial_manifest = staged_publication(
        publication,
        initial_staging,
        output,
        "fsync-initial",
        {"figures/result.pdf": "%PDF-1.4\nold\n"},
    )
    publication.promote_staged_publication(
        initial_staging,
        output,
        initial_products,
        initial_manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )
    previous = tree_bytes(output)
    staging = tmp_path / ".publication.staging-fsync-failure"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "fsync-failure",
        {"figures/result.pdf": "%PDF-1.4\nnew\n"},
    )
    real_fsync_regular_file = publication.fsync_regular_file

    def fail_pdf_fsync(path):
        if path.suffix == ".pdf":
            raise publication.PublicationError("simulated staged PDF fsync failure")
        real_fsync_regular_file(path)

    monkeypatch.setattr(publication, "fsync_regular_file", fail_pdf_fsync)

    with pytest.raises(
        publication.PublicationError, match="simulated staged PDF fsync failure"
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert tree_bytes(output) == previous
    assert json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    ) == initial_manifest


def test_main_commit_is_idempotent_and_uses_canonical_paths(
    publication, tmp_path, monkeypatch
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    acceptance = tmp_path / "acceptance"
    output = tmp_path / "publication"

    monkeypatch.setattr(
        publication,
        "discover_data",
        lambda analysis_path, _acceptance: empty_data(publication, analysis_path),
    )
    invocation = [
        str(analysis),
        "--output",
        str(output),
        "--acceptance",
        str(acceptance),
    ]

    publication.main(invocation)
    first = tree_bytes(output)
    manifest = json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    )

    assert manifest["normalized_invocation"] == publication.normalized_invocation(
        analysis.absolute(), output.absolute(), [acceptance]
    )
    assert manifest["publication_commit_contract"]["stable_parent_boundary"] == (
        "the output parent pathname must continue to name the held locked directory "
        "through commit return"
    )
    assert ".staging-" not in json.dumps(manifest)
    for product in manifest["products"]:
        path = Path(product["path"])
        assert path.is_relative_to(output.absolute())
        assert publication.source_binding(path) == product

    publication.main(invocation)

    assert tree_bytes(output) == first
    assert list(tmp_path.glob(".publication.staging-*")) == []


def test_analysis_as_output_is_rejected_without_deleting_evidence(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "evidence"})
    publication.write_text(analysis / "keep.txt", "keep evidence\n")
    previous = tree_bytes(analysis)

    with pytest.raises(
        publication.PublicationError, match="equal or contain the analysis root"
    ):
        publication.main([str(analysis), "--output", str(analysis)])

    assert tree_bytes(analysis) == previous
    assert not publication.publication_ownership_path(analysis).exists()
    assert not publication.publication_lock_path(analysis).exists()


def test_output_containing_discovered_evidence_is_rejected_without_deletion(
    publication, tmp_path, monkeypatch
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    output = tmp_path / "publication"
    evidence = output / "acceptance/evidence.json"
    write_json(evidence, {"record_type": "discovered-evidence"})
    publication.write_text(output / "keep.txt", "keep output evidence\n")
    previous = tree_bytes(output)
    data = empty_data(publication, analysis)
    data.source_paths.add(evidence.absolute())
    monkeypatch.setattr(
        publication, "discover_data", lambda _analysis, _acceptance: data
    )

    with pytest.raises(
        publication.PublicationError,
        match="contains discovered source or acceptance evidence",
    ):
        publication.main([str(analysis), "--output", str(output)])

    assert tree_bytes(output) == previous
    assert not publication.publication_ownership_path(output).exists()
    assert not publication.publication_lock_path(output).exists()


@pytest.mark.parametrize("nested_output", [False, True], ids=["equal", "contains"])
def test_explicit_empty_acceptance_root_equal_or_contains_output_is_rejected(
    publication, tmp_path, nested_output
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    acceptance = tmp_path / "acceptance"
    acceptance.mkdir()
    output = acceptance / "publication" if nested_output else acceptance

    with pytest.raises(
        publication.PublicationError,
        match="overlaps an acceptance root",
    ):
        publication.main([
            str(analysis),
            "--output",
            str(output),
            "--acceptance",
            str(acceptance),
        ])

    assert tree_bytes(acceptance) == {}
    assert not publication.publication_ownership_path(output).exists()
    assert not publication.publication_lock_path(output).exists()
    if nested_output:
        assert not output.exists()


@pytest.mark.parametrize(
    "root_name",
    ["scientific-acceptance", "acceptance", "audits", "scientific-audits"],
)
@pytest.mark.parametrize("nested_output", [False, True], ids=["equal", "contains"])
def test_builtin_acceptance_root_equal_or_containing_output_is_rejected(
    publication, tmp_path, root_name, nested_output
):
    analysis = tmp_path / "analysis"
    root = analysis / root_name
    root.mkdir(parents=True)
    output = root / "publication" if nested_output else root

    with pytest.raises(
        publication.PublicationError, match="overlaps an acceptance root"
    ):
        publication.main([str(analysis), "--output", str(output)])

    assert tree_bytes(root) == {}
    assert output.exists() is not nested_output
    assert not publication.publication_ownership_path(output).exists()
    assert not publication.publication_lock_path(output).exists()


def test_arbitrary_nonempty_unowned_output_is_rejected(
    publication, tmp_path, monkeypatch
):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    output = tmp_path / "publication"
    publication.write_text(output / "new.txt", "do not delete\n")
    previous = tree_bytes(output)

    monkeypatch.setattr(
        publication,
        "discover_data",
        lambda analysis_path, _acceptance: empty_data(publication, analysis_path),
    )

    def fixed_render(_data, staging):
        product = staging / "new.txt"
        publication.write_text(product, "new publication\n")
        return [product]

    monkeypatch.setattr(publication, "render_products", fixed_render)

    with pytest.raises(publication.PublicationError, match="nonempty unowned output"):
        publication.main([str(analysis), "--output", str(output)])

    assert tree_bytes(output) == previous
    assert not publication.publication_ownership_path(output).exists()
    assert list(tmp_path.glob(".publication.staging-*")) == []


def test_default_output_descendant_of_analysis_works_without_evidence_overlap(
    publication, tmp_path, monkeypatch
):
    analysis = tmp_path / "analysis"
    inventory = analysis / "inventory.json"
    write_json(inventory, {"record_type": "fixture-inventory"})
    inventory_bytes = inventory.read_bytes()
    output = analysis / "publication-products"

    def fixed_render(_data, staging):
        product = staging / "tables/summary.csv"
        publication.write_text(product, "result\npass\n")
        return [product]

    monkeypatch.setattr(publication, "render_products", fixed_render)

    publication.main([str(analysis)])

    assert inventory.read_bytes() == inventory_bytes
    assert publication.valid_publication_ownership(output)
    assert (output / "tables/summary.csv").read_text(encoding="utf-8") == (
        "result\npass\n"
    )
    manifest = json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["analysis_output"] == str(analysis.absolute())
    assert manifest["normalized_invocation"] == publication.normalized_invocation(
        analysis.absolute(), output.absolute(), []
    )
    assert all(
        Path(product["path"]).is_relative_to(output.absolute())
        for product in manifest["products"]
    )
    first = tree_bytes(output)

    publication.main([str(analysis)])

    assert tree_bytes(output) == first


def test_valid_existing_manifest_migrates_ownership_before_interruption(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    old_products = [output / "a.txt", output / "b.txt"]
    for product in old_products:
        publication.write_text(product, f"old {product.stem}\n")
    publication.write_json(
        output / "manifest.json",
        {
            "record_type": "cgl_lf_stage_i_fast_publication_products",
            "normalized_invocation": [
                "renderer",
                "--output",
                str(output.absolute()),
            ],
            "products": [
                publication.source_binding(product) for product in old_products
            ],
        },
    )
    assert publication.valid_existing_publication_manifest(output)
    ownership = publication.publication_ownership_path(output)
    write_json(ownership, {
        "schema_version": 2,
        "record_type": publication.PUBLICATION_OWNERSHIP_RECORD_TYPE,
        "output": str(output.absolute()),
        "directory_identity": {
            "device": output.stat().st_dev,
            "inode": output.stat().st_ino,
        },
    })
    assert not publication.valid_publication_ownership(output)
    staging = tmp_path / ".publication.staging-existing-manifest"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "existing-manifest",
        {"a.txt": "new a\n", "b.txt": "new b\n"},
    )
    real_promote_file = publication.promote_file
    product_promotions = 0

    def interrupt_after_migrated_ownership(
        staged_path, relative, output_descriptor
    ):
        nonlocal product_promotions
        real_promote_file(staged_path, relative, output_descriptor)
        if relative.name != "manifest.json":
            product_promotions += 1
            if product_promotions == 1:
                assert publication.valid_publication_ownership(output)
                token = (
                    output / publication.PUBLICATION_AUTHORITY_TOKEN_NAME
                ).read_bytes()
                owner = json.loads(ownership.read_text(encoding="utf-8"))
                assert owner["schema_version"] == 3
                assert owner["token_sha256"] == hashlib.sha256(token).hexdigest()
                assert "directory_identity" not in owner
                assert not (output / "manifest.json").exists()
                raise RuntimeError("interrupted migrated publication")

    monkeypatch.setattr(
        publication, "promote_file", interrupt_after_migrated_ownership
    )

    with pytest.raises(RuntimeError, match="interrupted migrated publication"):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert publication.valid_publication_ownership(output)
    ownership_bytes = ownership.read_bytes()
    assert not (output / "manifest.json").exists()
    assert (output / "a.txt").read_text(encoding="utf-8") == "new a\n"
    assert (output / "b.txt").read_text(encoding="utf-8") == "old b\n"

    monkeypatch.setattr(publication, "promote_file", real_promote_file)
    recovery_staging = tmp_path / ".publication.staging-existing-manifest-recovery"
    recovery_products, recovery_manifest = staged_publication(
        publication,
        recovery_staging,
        output,
        "existing-manifest",
        {"a.txt": "new a\n", "b.txt": "new b\n"},
    )
    publication.promote_staged_publication(
        recovery_staging,
        output,
        recovery_products,
        recovery_manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )

    assert publication.valid_publication_ownership(output)
    assert ownership.read_bytes() == ownership_bytes
    assert (output / "a.txt").read_text(encoding="utf-8") == "new a\n"
    assert (output / "b.txt").read_text(encoding="utf-8") == "new b\n"
    assert json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    ) == recovery_manifest


def test_deleted_recreated_output_cannot_reuse_persistent_authority(
    publication, tmp_path
):
    output = tmp_path / "publication"
    staging = tmp_path / ".publication.staging-owned-instance"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "owned-instance",
        {"owned.txt": "owned publication\n"},
    )
    publication.promote_staged_publication(
        staging,
        output,
        products,
        manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )
    assert publication.valid_publication_ownership(output)
    ownership = json.loads(
        publication.publication_ownership_path(output).read_text(encoding="utf-8")
    )
    assert ownership["schema_version"] == 3
    assert "token_sha256" in ownership
    assert "directory_identity" not in ownership

    displaced = tmp_path / "displaced-publication"
    output.rename(displaced)
    output.mkdir()
    replacement = tree_bytes(output)

    assert not publication.valid_publication_ownership(output)
    replacement_staging = tmp_path / ".publication.staging-replacement"
    replacement_products, replacement_manifest = staged_publication(
        publication,
        replacement_staging,
        output,
        "replacement",
        {"new.txt": "new publication\n"},
    )

    with pytest.raises(
        publication.PublicationError,
        match="does not bind the current output directory instance",
    ):
        publication.promote_staged_publication(
            replacement_staging,
            output,
            replacement_products,
            replacement_manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert tree_bytes(output) == replacement
    assert (displaced / "owned.txt").read_text(encoding="utf-8") == (
        "owned publication\n"
    )


def test_rerender_preserves_unrecognized_preexisting_regular_manifest(
    publication, tmp_path
):
    output = tmp_path / "publication"
    initial_staging = tmp_path / ".publication.staging-recognized-manifest"
    initial_products, initial_manifest = staged_publication(
        publication,
        initial_staging,
        output,
        "recognized-manifest",
        {"a.txt": "old a\n"},
    )
    publication.promote_staged_publication(
        initial_staging,
        output,
        initial_products,
        initial_manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )
    unknown_manifest = b'{"authority":"unrecognized-preexisting"}\n'
    replacement = output / ".replacement-manifest"
    replacement.write_bytes(unknown_manifest)
    os.replace(replacement, output / "manifest.json")
    previous = tree_bytes(output)
    staging = tmp_path / ".publication.staging-rerender-unrecognized-manifest"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "rerender-unrecognized-manifest",
        {"a.txt": "new a\n"},
    )

    with pytest.raises(
        publication.PublicationError,
        match="manifest is unrecognized; refusing to replace it",
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert tree_bytes(output) == previous
    assert (output / "manifest.json").read_bytes() == unknown_manifest


def test_rerender_preserves_manifest_replaced_after_recognition_before_withdrawal(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    initial_staging = tmp_path / ".publication.staging-recognition-race-initial"
    initial_products, initial_manifest = staged_publication(
        publication,
        initial_staging,
        output,
        "recognition-race-initial",
        {"a.txt": "old a\n"},
    )
    publication.promote_staged_publication(
        initial_staging,
        output,
        initial_products,
        initial_manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )
    staging = tmp_path / ".publication.staging-recognition-race-rerender"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "recognition-race-rerender",
        {"a.txt": "new a\n"},
    )
    unknown_manifest = b'{"authority":"replacement-after-recognition"}\n'
    real_fsync_regular_file = publication.fsync_regular_file
    replacement_installed = False

    def replace_after_recognition(path):
        nonlocal replacement_installed
        real_fsync_regular_file(path)
        if path == staging / "manifest.json" and not replacement_installed:
            replacement = output / ".replacement-manifest"
            replacement.write_bytes(unknown_manifest)
            os.replace(replacement, output / "manifest.json")
            replacement_installed = True

    monkeypatch.setattr(
        publication, "fsync_regular_file", replace_after_recognition
    )

    with pytest.raises(
        publication.PublicationError,
        match="changed before authority withdrawal",
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert replacement_installed
    assert (output / "manifest.json").read_bytes() == unknown_manifest
    assert (output / "a.txt").read_text(encoding="utf-8") == "old a\n"
    assert not any(
        path.name.startswith(publication.PUBLICATION_MANIFEST_QUARANTINE_PREFIX)
        for path in output.iterdir()
    )


@pytest.mark.parametrize(
    ("replace_after", "expected_error"),
    [
        ("a.txt", "changed during product promotion"),
        ("manifest.json", "changed before commit return"),
    ],
)
def test_canonical_replacement_during_promotion_never_blesses_replacement(
    publication, tmp_path, monkeypatch, replace_after, expected_error,
):
    output = tmp_path / "publication"
    initial_staging = tmp_path / ".publication.staging-race-initial"
    initial_products, initial_manifest = staged_publication(
        publication,
        initial_staging,
        output,
        "race-initial",
        {"a.txt": "old a\n", "b.txt": "old b\n"},
    )
    publication.promote_staged_publication(
        initial_staging,
        output,
        initial_products,
        initial_manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )
    staging = tmp_path / ".publication.staging-race-replacement"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "race-replacement",
        {"a.txt": "new a\n", "b.txt": "new b\n"},
    )
    real_promote_file = publication.promote_file
    displaced = tmp_path / "displaced-during-promotion"
    replaced = False

    def replace_canonical_after_first_product(
        staged, relative, output_descriptor
    ):
        nonlocal replaced
        real_promote_file(staged, relative, output_descriptor)
        if relative.name == replace_after and not replaced:
            output.rename(displaced)
            output.mkdir()
            publication.write_text(output / "unrelated.txt", "replacement\n")
            replaced = True

    monkeypatch.setattr(
        publication, "promote_file", replace_canonical_after_first_product
    )

    with pytest.raises(
        publication.PublicationError,
        match=expected_error,
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert tree_bytes(output) == {"unrelated.txt": b"replacement\n"}
    assert not (output / "manifest.json").exists()
    assert not (displaced / "manifest.json").exists()
    assert (displaced / "a.txt").read_text(encoding="utf-8") == "new a\n"
    assert (displaced / "b.txt").read_text(encoding="utf-8") == "new b\n"


def test_concurrent_promotions_survive_sibling_lock_path_replacement(
    publication, tmp_path, monkeypatch
):
    output = tmp_path / "publication"
    first_staging = tmp_path / ".publication.staging-first"
    second_staging = tmp_path / ".publication.staging-second"
    first_products, first_manifest = staged_publication(
        publication,
        first_staging,
        output,
        "first",
        {"a.txt": "first a\n", "b.txt": "first b\n"},
    )
    second_products, second_manifest = staged_publication(
        publication,
        second_staging,
        output,
        "second",
        {"a.txt": "second a\n", "b.txt": "second b\n"},
    )

    real_promote_file = publication.promote_file
    first_paused = threading.Event()
    release_first = threading.Event()
    second_started = threading.Event()
    second_finished = threading.Event()
    order = []
    errors = []

    def observed_promote(staged, relative, output_descriptor):
        publisher = (
            "first" if staged.is_relative_to(first_staging) else "second"
        )
        order.append((publisher, relative.name))
        real_promote_file(staged, relative, output_descriptor)
        if publisher == "first" and relative.name == "a.txt":
            first_paused.set()
            assert release_first.wait(10)

    def run_promotion(products, manifest, staging, started=None, finished=None):
        try:
            if started is not None:
                started.set()
            publication.promote_staged_publication(
                staging,
                output,
                products,
                manifest,
                analysis=output.parent / "analysis",
                evidence_paths=[],
            )
        except BaseException as error:
            errors.append(error)
        finally:
            if finished is not None:
                finished.set()

    monkeypatch.setattr(publication, "promote_file", observed_promote)
    first = threading.Thread(
        target=run_promotion,
        args=(first_products, first_manifest, first_staging),
    )
    second = threading.Thread(
        target=run_promotion,
        args=(
            second_products,
            second_manifest,
            second_staging,
            second_started,
            second_finished,
        ),
    )

    first.start()
    assert first_paused.wait(10)
    lock_path = publication.publication_lock_path(output)
    displaced_lock = tmp_path / "displaced-publication.lock"
    lock_path.rename(displaced_lock)
    publication.write_text(lock_path, "replacement lock pathname\n")
    second.start()
    assert second_started.wait(10)
    try:
        assert not second_finished.wait(0.5)
        assert all(publisher == "first" for publisher, _ in order)
    finally:
        release_first.set()
    first.join(10)
    second.join(10)

    assert not first.is_alive()
    assert not second.is_alive()
    assert errors == []
    assert publication.publication_lock_path(output).is_file()
    assert not publication.publication_lock_path(output).is_symlink()
    assert displaced_lock.is_file()
    assert lock_path.read_text(encoding="utf-8") == "replacement lock pathname\n"
    assert [publisher for publisher, _ in order] == [
        "first",
        "first",
        "first",
        "second",
        "second",
        "second",
    ]
    assert json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    ) == second_manifest
    assert (output / "a.txt").read_text(encoding="utf-8") == "second a\n"
    assert (output / "b.txt").read_text(encoding="utf-8") == "second b\n"
    for product in second_manifest["products"]:
        assert publication.source_binding(Path(product["path"])) == product


def test_parent_directory_replacement_invalidates_trusted_lock_boundary(
    publication, tmp_path
):
    parent = tmp_path / "parent"
    parent.mkdir()
    output = parent / "publication"

    with publication.exclusive_publication_lock(output) as parent_descriptor:
        assert publication.canonical_parent_resolves_to_held_parent(
            output, parent_descriptor
        )
        displaced = tmp_path / "displaced-parent"
        parent.rename(displaced)
        parent.mkdir()

        assert not publication.canonical_parent_resolves_to_held_parent(
            output, parent_descriptor
        )
        with pytest.raises(
            publication.PublicationError,
            match="held trusted stable parent boundary",
        ):
            publication.require_stable_publication_parent(
                output, parent_descriptor, "during adversarial test"
            )


def test_parent_replacement_cannot_redirect_ownership_write(
    publication, tmp_path, monkeypatch
):
    parent = tmp_path / "parent"
    parent.mkdir()
    output = parent / "publication"
    output.mkdir()
    output_descriptor = publication.open_output_directory(output)
    token = publication.create_publication_authority_token(output_descriptor)
    real_write = publication.write_durable_regular_file_at
    replacement_sentinel = b'{"authority":"unrelated-sentinel"}\n'
    displaced = tmp_path / "displaced-parent"
    replaced = False

    def replace_parent_then_write(parent_descriptor, name, payload):
        nonlocal replaced
        parent.rename(displaced)
        parent.mkdir()
        (parent / name).write_bytes(replacement_sentinel)
        replaced = True
        real_write(parent_descriptor, name, payload)

    monkeypatch.setattr(
        publication, "write_durable_regular_file_at", replace_parent_then_write
    )
    parent_descriptor = os.open(parent, publication.directory_open_flags())
    try:
        with pytest.raises(
            publication.PublicationError,
            match="held trusted stable parent boundary",
        ):
            publication.write_durable_publication_ownership(
                output, parent_descriptor, output_descriptor, token
            )
    finally:
        os.close(parent_descriptor)
        os.close(output_descriptor)

    assert replaced
    ownership_name = publication.publication_ownership_path(output).name
    assert (parent / ownership_name).read_bytes() == replacement_sentinel
    assert json.loads((displaced / ownership_name).read_text(encoding="utf-8")) == (
        publication.publication_ownership_record(output, token)
    )


def test_output_root_symlink_escape_is_rejected(publication, tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    output = tmp_path / "publication"
    output.symlink_to(outside, target_is_directory=True)
    staging = tmp_path / ".publication.staging-root-symlink"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "symlink-root",
        {"tables/new.csv": "new\n"},
    )

    with pytest.raises(publication.PublicationError, match="output root.*symlink"):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert tree_bytes(outside) == {}


def test_output_descendant_symlink_escape_is_rejected(publication, tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    output = tmp_path / "publication"
    output.mkdir()
    publication.write_json(output / "manifest.json", {"authority": "old"})
    (output / "tables").symlink_to(outside, target_is_directory=True)
    staging = tmp_path / ".publication.staging-tables-symlink"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "symlink-tables",
        {"tables/new.csv": "new\n"},
    )

    with pytest.raises(
        publication.PublicationError, match="descendant.*symlink"
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    ) == {"authority": "old"}
    assert tree_bytes(outside) == {}


def test_output_descendant_file_symlink_is_rejected(publication, tmp_path):
    outside = tmp_path / "outside.txt"
    outside.write_text("outside\n", encoding="utf-8")
    output = tmp_path / "publication"
    output.mkdir()
    publication.write_json(output / "manifest.json", {"authority": "old"})
    (output / "tables").mkdir()
    (output / "tables/old.csv").symlink_to(outside)
    staging = tmp_path / ".publication.staging-file-symlink"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "symlink-file",
        {"tables/new.csv": "new\n"},
    )

    with pytest.raises(
        publication.PublicationError, match="descendant.*symlink"
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert json.loads(
        (output / "manifest.json").read_text(encoding="utf-8")
    ) == {"authority": "old"}
    assert outside.read_text(encoding="utf-8") == "outside\n"


def test_rerender_fails_closed_and_preserves_stale_canonical_paths(
    publication, tmp_path
):
    output = tmp_path / "publication"
    initial_staging = tmp_path / ".publication.staging-initial-tree"
    initial_products, initial_manifest = staged_publication(
        publication,
        initial_staging,
        output,
        "initial-tree",
        {
            "tables/current.csv": "old current\n",
            "tables/obsolete.csv": "obsolete\n",
            "obsolete/old.txt": "obsolete\n",
        },
    )
    publication.promote_staged_publication(
        initial_staging,
        output,
        initial_products,
        initial_manifest,
        analysis=tmp_path / "analysis",
        evidence_paths=[],
    )
    staging = tmp_path / ".publication.staging-exact-tree"
    products, manifest = staged_publication(
        publication,
        staging,
        output,
        "exact-tree",
        {"tables/current.csv": "new current\n"},
    )
    previous = tree_bytes(output)

    with pytest.raises(
        publication.PublicationError, match="stale or unrecognized paths"
    ):
        publication.promote_staged_publication(
            staging,
            output,
            products,
            manifest,
            analysis=tmp_path / "analysis",
            evidence_paths=[],
        )

    assert tree_bytes(output) == previous
    assert (output / "tables/obsolete.csv").read_text(encoding="utf-8") == (
        "obsolete\n"
    )
    assert (output / "obsolete/old.txt").read_text(encoding="utf-8") == "obsolete\n"


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
    assert result["result"] == "negative"
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
    assert limiter["applied_heat_flux_work_total"] is None
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


def test_headline_sensitivity_rows_use_reviewed_science_not_fast_report(
    publication, tmp_path
):
    data = empty_data(publication, tmp_path)
    data.cases["R02"].diagnostics = {
        "windows": {
            "steady": {
                "history": {"analysis_window": {"kinetic_mean": 10.0}}
            }
        }
    }
    data.cases["R04"].diagnostics = {
        "windows": {
            "steady": {
                "history": {"analysis_window": {"kinetic_mean": 20.0}}
            }
        }
    }
    data.science_record = {
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
                        "pooled_within_realization_standard_deviation": 0.32,
                        "standardized_effect": -1.25,
                        "standardized_effect_scope": (
                            "descriptive_within_realization"
                        ),
                        "claim_scope": "descriptive_within_realization",
                    }],
                }
            }
        }
    }

    row = next(
        value for value in publication.robustness_rows(data)
        if value["contrast"] == "forcing A, beta=10"
        and value["metric"] == "kinetic"
    )

    assert row["reference_value"] == pytest.approx(1.0)
    assert row["variant_value"] == pytest.approx(1.4)
    assert row["variant_minus_reference"] == pytest.approx(0.4)
    assert row["reviewed_left_minus_right"] == pytest.approx(-0.4)
    assert row["standardized_effect"] == pytest.approx(-1.25)
    assert row["authority"] == publication.SCIENCE_AUTHORITY
    assert row["release_authorizing"] is False


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

    assert publication.hyperbolicity_status(data, "R14") == (
        "nonnegative_discriminant"
    )
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
    science = write_science_record(tmp_path / "science", analysis, publication)
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
    assert contrast["inference_scope"] == "descriptive_within_realization"
    assert contrast["full_scope_independent_review_complete"] is False
    assert contrast["release_authorizing"] is False
    assert resolution["available"] is False
    assert resolution["passed"] is None
    assert mks24["result"] == "inconclusive"
    assert ct_row["ct_result"] == "pass"
    assert ct_row["campaign_authority_eligible"] is False
    assert ct_row["release_authorizing"] is False
    assert publication.health_rows(data)[0]["direct_ct_numerical"] == "pass"
    report = publication.report_markdown(data, [], tmp_path)
    assert len(publication.integrated_audit_records(data)) == 1
    assert publication.integrated_audit_records(data)[0] is data.ct_audit_record
    assert "- Audit records discovered: `1`" in report
    assert "non-authorizing-direct-fast-scientific-assessment" in report
    assert "campaign_authority_eligible=false" in report


@pytest.mark.parametrize(
    "mutation",
    (
        "missing_current_scope",
        "missing_intervention_scope",
        "deprecated_inference",
    ),
)
def test_publication_rejects_unscoped_or_inferential_science(
    publication, tmp_path, mutation
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    science = write_science_record(tmp_path / "science", analysis, publication)
    record = json.loads(science.read_text(encoding="utf-8"))
    if mutation == "missing_current_scope":
        record.pop("current_science_scope_limitation")
    elif mutation == "missing_intervention_scope":
        record.pop("active_passive_intervention_scope")
    else:
        record["families"]["forcing"]["R02_R04"]["metrics"][0]["two_sided_p"] = 0.01
    rewrite_science_record(science, record)

    data = publication.discover_data(analysis, [science])

    assert data.science_record is None
    assert any(
        "reviewed science aggregate" in warning for warning in data.ingestion_warnings
    )


def test_science_pass_allows_explicit_ineligible_r10_and_preserves_available_contrast(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    science = write_science_record(tmp_path / "science", analysis, publication)
    record = json.loads(science.read_text(encoding="utf-8"))
    record["result"] = "pass"
    record["selected_cases"] = list(publication.CASE_IDS)
    record["case_dispositions"] = {
        case_id: {
            "claim_eligible": case_id != "R10",
            "acceptance_result": "inconclusive" if case_id == "R10" else "pass",
        }
        for case_id in publication.CASE_IDS
    }
    record["families"]["descriptive"] = {
        "R02_R03": {
            "left": "R02",
            "right": "R03",
            "result": "available",
            "claim_eligible": True,
            "metrics": [{"metric": "kinetic", "available": True}],
        }
    }
    record["gates"] = [{
        "name": "reviewed_campaign_science",
        "result": "pass",
        "reason": "all required science gates passed",
        "observations": [],
    }]
    rewrite_science_record(science, record)

    data = publication.discover_data(analysis, [science])
    rows = {
        row["contrast"]: row for row in publication.science_contrast_rows(data)
    }

    assert data.science_record is not None
    assert publication.aggregate_science_status(data) == "pass"
    assert publication.science_case_status(data, "R10") == "inconclusive"
    assert rows["R02_R03"]["result"] == "available"


def test_rendered_products_expose_authenticated_science_and_ct(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    science = write_science_record(tmp_path / "science", analysis, publication)
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
    science = write_science_record(tmp_path / "science", analysis, publication)
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


def test_ct_audit_rejects_unselected_case_records(publication, tmp_path):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    ct = write_ct_record(tmp_path / "ct", analysis)
    record = json.loads(ct.read_text(encoding="utf-8"))
    record["cases"]["R03"] = dict(record["cases"]["R02"])
    rewrite_ct_record(ct, record)

    data = publication.discover_data(analysis, [ct])

    assert data.ct_audit_record is None
    assert any(
        "cases differ from selected-case inventory" in warning
        for warning in data.ingestion_warnings
    )


def test_ct_audit_rejects_pass_inconsistent_with_numerical_evidence(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_json(analysis / "inventory.json", {"record_type": "fixture-inventory"})
    ct = write_ct_record(tmp_path / "ct", analysis)
    record = json.loads(ct.read_text(encoding="utf-8"))
    record["cases"]["R02"]["native_restart_ct"]["maximum_normalized_ct_divb"] = 2.0e-10
    rewrite_ct_record(ct, record)

    data = publication.discover_data(analysis, [ct])

    assert data.ct_audit_record is None
    assert any(
        "CT result inconsistent with its numerical evidence" in warning
        for warning in data.ingestion_warnings
    )


def test_ct_campaign_coverage_requires_complete_native_coverage(publication, tmp_path):
    data = empty_data(publication, tmp_path)
    data.ct_audit_record = {
        "selection": {"cases": list(publication.CASE_IDS)},
        "cases": {
            case_id: {"native_restart_ct": {"coverage_complete": True}}
            for case_id in publication.CASE_IDS
        },
    }

    assert publication.ct_campaign_coverage_complete(data) is True
    data.ct_audit_record["cases"]["R17"]["native_restart_ct"][
        "coverage_complete"
    ] = False
    assert publication.ct_campaign_coverage_complete(data) is False


def test_r15_selected_nonfatal_variant_preserves_strict_failure_scope(
    publication, tmp_path
):
    data = empty_data(publication, tmp_path)
    data.cases["R15"].lineage = {
        "lineage_variants": ["finite_limiter_hard_bound_diagnostic_nonfatal"],
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
    science_bound_acceptance = tmp_path / "science-bound/R15/case_acceptance.json"
    write_json(science_bound_acceptance, {
        "scope": {
            "retained_strict_failure_evidence": [
                strict_failure_record(tmp_path / "strict-failure")
            ]
        },
    })
    data.science_record = {
        "_publication_evidence_validated": True,
        "result": "inconclusive",
        "case_dispositions": {
            "R15": {"claim_eligible": True, "acceptance_result": "pass"}
        },
        "provenance": {
            "case_acceptance": {"R15": binding(science_bound_acceptance)},
            "case_lineages": {"R15": binding(lineage_path)},
        },
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


def test_r15_strict_failure_rejects_stale_bound_provenance(publication, tmp_path):
    data = empty_data(publication, tmp_path)
    failure = strict_failure_record(tmp_path / "strict-failure")
    data.cases["R15"].direct_acceptance = {
        "_publication_evidence_validated": True,
        "scope": {"retained_strict_failure_evidence": [failure]},
    }
    Path(failure["provenance"]["slurm_log"]["path"]).write_text(
        "changed after binding\n", encoding="utf-8"
    )

    assert publication.r15_strict_failure_disposition(data) == "unavailable/inconclusive"


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


def material_diagnostics(case_id: str) -> dict[str, object]:
    offset = int(case_id[1:]) / 100.0
    return {
        "health": {
            "result": "clean",
            "final_time": 10.0,
            "mass_relative_drift": {"mhd": 2.0e-16, "user": 3.0e-16},
            "mhd_user_mass_relative_mismatch": 0.0,
        },
        "snapshot_analysis_status": "complete",
        "windows": {
            "steady": {
                "lf_history": {
                    "available": True,
                    "applied_heat_flux_work": {
                        "signed": True,
                        "parallel": 2.0 + offset,
                        "perpendicular": -0.5 - offset,
                        "total": 1.5,
                    },
                    "applied_pressure_work": {
                        "signed": True,
                        "total": -1.0 - offset,
                        "anisotropic": -0.8 - offset,
                    },
                    "heat_flux_cap_fractions": {
                        "parallel_over_1": 1.0e-4 + offset * 1.0e-4,
                        "parallel_over_10": 0.0,
                        "perpendicular_over_1": 2.0e-4 + offset * 1.0e-4,
                        "perpendicular_over_10": 0.0,
                    },
                }
            }
        },
        "snapshot_ensemble": {
            "snapshot_count": 3,
            "pressure_work_decomposition": {
                "available": True,
                "snapshot_count": 3,
                "applied_to_flow": True,
                "isotropic_perpendicular_pressure_power_mean": 0.5 + offset,
                "anisotropic_stress_power_mean": -0.2 - offset,
                "total_cgl_pressure_power_mean": 0.3,
                "parallel_strain_rms_mean": 0.7 + offset,
                "time_integral_estimate": {
                    "available": True,
                    "anisotropic_stress_power_integral": -1.2 - offset,
                },
            },
            "pressure_balance": {
                "available": True,
                "snapshot_count": 3,
                "correlation_mean": -0.9 + offset,
                "normalized_residual_variance_mean": 0.1 + offset,
            },
            "spectral_scalar_diagnostics": {
                "compressive_velocity_power_fraction": {
                    "available": True,
                    "snapshot_count": 3,
                    "fraction_mean": 0.2 - offset,
                },
            },
            "heat_flux_transport_proxy": {
                "available": True,
                "snapshot_count": 3,
                "regularized_total_power_mean": 0.4 + offset,
                "unlimited_total_power_mean": 0.6 + offset,
                "parallel_cap_active_volume_fraction_mean": 0.01 + offset,
                "perpendicular_cap_active_volume_fraction_mean": 0.02 + offset,
                "time_integral_estimate": {
                    "available": True,
                    "regularized_total_power_integral": 2.4 + offset,
                    "unlimited_total_power_integral": 3.6 + offset,
                },
            },
            "pdf": {
                "bb_grad_velocity": {
                    "edges": [-1.0, 0.0, 1.0],
                    "density": [0.4 + offset, 0.6 - offset],
                }
            },
            "spectra": {
                "velocity": {
                    "dk": 3.141592653589793,
                    "k": [
                        12.566370614359172,
                        25.132741228718345,
                        50.26548245743669,
                        75.39822368615503,
                    ],
                    "power_per_dk": [4.0 + offset, 3.0, 2.0, 1.0],
                },
                "magnetic_fluctuation": {
                    "k": [
                        12.566370614359172,
                        25.132741228718345,
                        50.26548245743669,
                        75.39822368615503,
                    ],
                    "power_per_dk": [3.0 + offset, 2.5, 1.5, 0.8],
                },
            },
            "alignment": {
                shell: {"edges": [0.0, 0.5, 1.0], "density": [1.0, 2.0]}
                for shell in ("4", "8", "16", "24")
            },
        },
    }


def install_material_case(publication, data, analysis: Path, case_id: str) -> dict:
    diagnostics = material_diagnostics(case_id)
    diagnostics_path = analysis / "cases" / case_id / "diagnostics.json"
    write_json(diagnostics_path, diagnostics)
    case = data.cases[case_id]
    case.diagnostics = diagnostics
    case.lineage = {
        "status": "complete",
        "final_time": 10.0,
        "selected_fast_lineage": {
            "reason": "highest_ranked_restart_linked_lineage",
            "terminal": {
                "state": "complete",
                "source_family": "original",
                "job_id": f"selected-{case_id}",
                "observed_final_time": 10.0,
                "restart_link_valid": True,
                "segment": f"/selected/{case_id}",
            },
            "segments": [{"segment": f"/selected/{case_id}"}],
        },
        "unselected_lineages": [{
            "reason": "lower_ranked_restart_linked_lineage",
            "terminal": {
                "state": "failed",
                "source_family": "race",
                "job_id": f"failed-{case_id}",
                "observed_final_time": 1.0,
                "run_exit_code": 143,
                "restart_link_valid": True,
                "segment": f"/failed/{case_id}",
            },
            "segments": [{"segment": f"/failed/{case_id}"}],
        }],
    }
    lineage_path = analysis / "cases" / case_id / "lineage.json"
    write_json(lineage_path, case.lineage)
    case.lineage_path = lineage_path
    case.user_history = {
        "time": [0.0, 4.0, 10.0],
        "volume": [1.0, 1.0, 1.0],
        "b2": [2.0, 2.0, 2.0],
        "b4": [4.08, 4.16, 4.2],
    }
    user_history_path = analysis / "cases" / case_id / "history/material.user.hst"
    user_history_path.parent.mkdir(parents=True, exist_ok=True)
    user_history_path.write_text(
        "# [1]=time [2]=volume [3]=b2 [4]=b4\n"
        "0 1 2 4.08\n4 1 2 4.16\n10 1 2 4.2\n",
        encoding="utf-8",
    )
    case.history_paths["user"] = user_history_path
    case.direct_acceptance = {
        "_publication_evidence_validated": True,
        "result": "pass",
        "health": {
            "result": "pass",
            "complete_to_target": True,
            "observed_final_time": 10.0,
            "fatal_counter_maxima": {
                "lf_dfloor": 0,
                "lf_pfloor": 0,
                "lf_nonfin": 0,
                "lf_nonpos": 0,
            },
        },
        "scope": {"classification": "standard_claim_scope"},
        "provenance": {
            "lineage": binding(lineage_path),
            "histories": {"user": binding(user_history_path)},
        },
        "history_statistics": {
            metric: {
                "history": "user",
                "column": metric,
                "sampling_adequacy": "pass",
                "stationarity": {"result": "pass"},
                "windows": {
                    "full": {
                        "statistics": {
                            "mean": 1.0 + int(case_id[1:]) / 100.0,
                            "standard_error": 0.05,
                            "confidence_interval_95": [0.9, 1.1],
                            "effective_sample_count": 3.0,
                        }
                    }
                },
            }
            for metric in publication.PRIMARY_SCALAR_METRICS
        },
    }
    case.acceptance = {
        "_publication_evidence_validated": True,
        "result": "pass",
        "gates": [{
            "name": "active_energy_closure",
            "result": "pass",
            "observations": {
                "windows": {
                    "whole_lineage": {
                        "increment_normalized_residual": 3.0e-12,
                        "state_normalized_mismatch": 2.0e-13,
                    },
                    "developed": {
                        "increment_normalized_residual": 4.0e-12,
                        "state_normalized_mismatch": 3.0e-13,
                    },
                }
            },
        }],
    }
    return binding(diagnostics_path)


def install_material_science(publication, data, bindings: dict[str, dict]) -> None:
    selected = sorted(bindings)
    data.science_record = {
        "_publication_evidence_validated": True,
        "result": "pass",
        "selected_cases": selected,
        "current_science_scope_limitation": publication.CURRENT_SCIENCE_SCOPE_LIMITATION,
        "active_passive_intervention_scope": publication.ACTIVE_PASSIVE_INTERVENTION_SCOPE,
        "case_dispositions": {
            case_id: {"claim_eligible": True, "acceptance_result": "pass"}
            for case_id in selected
        },
        "families": {
            "active_passive": {
                "R02_R06": {
                    "active": "R02",
                    "passive": "R06",
                    "result": "pass",
                    "claim_eligible": True,
                    "intervention_scope": publication.ACTIVE_PASSIVE_INTERVENTION_SCOPE,
                    "metrics": [{
                        "metric": "abs_dp",
                        "available": True,
                        "active_mean": 1.0,
                        "passive_mean": 2.0,
                        "difference": -1.0,
                        "pooled_within_realization_standard_deviation": 2.0 / 3.0,
                        "standardized_effect": -1.5,
                        "standardized_effect_scope": (
                            "descriptive_within_realization"
                        ),
                        "expected_direction": "active_lower",
                        "direction_coherent": True,
                        "large_direction_coherent_effect": True,
                        "claim_scope": "descriptive_within_realization",
                    }],
                }
            }
        },
        "resolution": {
            "result": "pass",
            "limits": {"common_k_perp_over_pi": [4.0, 24.0]},
            "observations": [{
                "kind": "curve",
                "product": "velocity_spectrum_shape",
                "available": True,
                "R02_R17_distance": 0.1,
                "R02_R17_limit": 0.2,
                "improved": True,
                "passed": True,
            }],
        },
        "gates": [],
        "provenance": {"case_diagnostics": bindings},
    }


def write_all_snapshot_hyperbolicity_evidence(
    publication, analysis: Path, root: Path, case_id: str = "R02",
    aggregates: tuple[dict[str, object], ...] | None = None,
    formula_id: str | None = "qualified-legacy",
    formula_provenance: object | None = None,
    executable_formula_id: str | None = None,
    formula_executable_compatibility: object | None = None,
) -> tuple[Path, Path]:
    if aggregates is None:
        aggregates = (
            {
                "evaluated": 1000,
                "negative": 0,
                "nonfinite_discriminant": 0,
                "minimum": 0.25,
            },
            {
                "evaluated": 1000,
                "negative": 0,
                "nonfinite_discriminant": 0,
                "minimum": 0.25,
            },
        )
    assert len(aggregates) == 2
    case_dir = analysis / "cases" / case_id
    rank_records = []
    selected = []
    result_snapshots = []
    for index, time in enumerate((4.0, 10.0)):
        rank = case_dir / "snapshots" / f"snapshot-{index}" / "rank_00000000.bin"
        rank.parent.mkdir(parents=True, exist_ok=True)
        rank.write_text(f"rank {index}\n", encoding="utf-8")
        rank_record = {
            **binding(rank),
            "mtime_ns": rank.stat().st_mtime_ns,
            "rank_id": 0,
        }
        rank_records.append(rank_record)
        selected.append({
            "index_position": index,
            "time": time,
            "audit_input_pattern": str(rank),
            "rank_files": [rank_record],
        })
        result_snapshots.append({
            "time": time,
            "active_cgl_signal_speed": True,
            "rank_files": [rank_record],
            "ranks_contiguous_from_zero": True,
            "input_inventory_sha256": hashlib.sha256(
                json.dumps(
                    [rank_record], sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
            ).hexdigest(),
            "aggregate": aggregates[index],
        })
    snapshot_index = case_dir / "snapshots.json"
    write_json(snapshot_index, {
        "complete_snapshot_count": 2,
        "snapshots": [{"time": 4.0}, {"time": 10.0}],
    })
    lineage = case_dir / "lineage.json"
    write_json(lineage, {
        "case_id": case_id,
        "snapshots": {
            "path": str(snapshot_index.absolute()),
            "complete_snapshot_count": 2,
        },
    })
    audit_script = root / "audit.py"
    bin_convert = root / "bin_convert.py"
    launcher = root / "launcher.py"
    for path in (audit_script, bin_convert, launcher):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(path.name + "\n", encoding="utf-8")
    provenance = {
        "script_path": str(audit_script.absolute()),
        "script_sha256": sha256(audit_script),
        "bin_convert_path": str(bin_convert.absolute()),
        "bin_convert_sha256": sha256(bin_convert),
        "input_patterns": [str(record["audit_input_pattern"]) for record in selected],
        "hash_inputs": True,
        "formula_reference": str(audit_script.absolute()),
    }
    if formula_id is not None:
        provenance["formula_id"] = formula_id
    if formula_provenance is None and formula_id is not None:
        formula_provenance = {"fixture_formula_id": formula_id}
    if formula_provenance is not None:
        provenance["formula_provenance"] = formula_provenance
    if executable_formula_id is not None:
        provenance["executable_formula_id"] = executable_formula_id
    if formula_executable_compatibility is not None:
        provenance["formula_executable_compatibility"] = (
            formula_executable_compatibility
        )
    result = {
        "provenance": provenance,
        "snapshots": result_snapshots,
    }
    attempt = root / case_id / "attempt-000"
    result_path = attempt / "result.json"
    write_json(result_path, result)
    result_sha = attempt / "result.sha256"
    result_sha.write_text(f"{sha256(result_path)}  result.json\n", encoding="utf-8")
    selected_digest = hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    manifest = {
        "schema_version": 1,
        "record_type": publication.HYPERBOLICITY_JOB_RECORD_TYPE,
        "case_id": case_id,
        "attempt": 0,
        "case_lineage": binding(lineage),
        "snapshot_index": binding(snapshot_index),
        "audit_script": binding(audit_script),
        "launcher": binding(launcher),
        "snapshot_policy": "all",
        "snapshot_coverage": {
            "snapshot_policy": "all",
            "snapshot_index_complete_count": 2,
            "selected_snapshot_count": 2,
            "selected_snapshot_positions": [0, 1],
            "selected_snapshot_times": [4.0, 10.0],
            "selected_snapshots_sha256": selected_digest,
            "all_complete_retained_snapshots_selected": True,
        },
        "selected_snapshots": selected,
        "selected_snapshot": selected[-1],
        "result": {
            "path": str(result_path.absolute()),
            "sha256_path": str(result_sha.absolute()),
            "required_coverage": "exactly_once_per_selected_snapshot",
            "expected_snapshot_count": 2,
            "expected_selected_snapshots_sha256": selected_digest,
        },
    }
    if formula_id is not None:
        manifest["formula_id"] = formula_id
    manifest_path = attempt / "manifest.json"
    write_json(manifest_path, manifest)
    return manifest_path, Path(str(rank_records[-1]["path"]))


def test_material_tables_authenticate_health_energy_and_r14_r15_failures(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    data = empty_data(publication, analysis)
    bindings = {
        case_id: install_material_case(publication, data, analysis, case_id)
        for case_id in ("R02", "R06", "R14", "R15")
    }
    install_material_science(publication, data, bindings)
    for case_id, failure_time, hard_bound, job_id in (
        ("R14", 0.1326939, 167, "4770854"),
        ("R15", 1.275643, 2, "4771183"),
    ):
        data.cases[case_id].direct_acceptance["scope"] = {
            "classification": "scoped_nonfatal_hard_bound_variant",
            "retained_strict_failure_evidence": [
                strict_failure_record(
                    tmp_path / f"strict-{case_id}", case_id,
                    failure_time=failure_time, hard_bound=hard_bound, job_id=job_id,
                )
            ],
        }

    health = {
        row["case_id"]: row
        for row in publication.numerical_health_provenance_rows(data)
    }
    scalars = publication.primary_full_window_scalar_rows(data)

    assert health["R02"]["completion"] == "pass"
    assert health["R02"]["mass_relative_drift_maximum"] == pytest.approx(3.0e-16)
    assert health["R02"]["active_energy_closure"] == "pass"
    assert health["R02"]["energy_increment_residual_maximum"] == pytest.approx(4.0e-12)
    assert health["R02"]["reporter_diagnostics_provenance"] == "authenticated"
    assert health["R06"]["active_energy_closure"] == "not_applicable"
    assert "hard_bound=167" in health["R14"]["strict_failure"]
    assert "hard_bound=2" in health["R15"]["strict_failure"]
    kinetic = next(
        row for row in scalars if row["case_id"] == "R02" and row["metric"] == "kinetic"
    )
    assert kinetic["availability"] == "available"
    assert kinetic["history"] == "user"
    assert kinetic["column"] == "kinetic"
    assert kinetic["standard_error"] == pytest.approx(0.05)
    assert kinetic["effective_sample_count"] == pytest.approx(3.0)
    assert kinetic["stationarity"] == "pass"

    diagnostics_path = analysis / "cases/R02/diagnostics.json"
    diagnostics_path.write_text("{}\n", encoding="utf-8")
    unauthenticated = {
        row["case_id"]: row
        for row in publication.numerical_health_provenance_rows(data)
    }
    assert unauthenticated["R02"]["mass_relative_drift_maximum"] is None
    assert unauthenticated["R02"]["reporter_diagnostics_provenance"] == "inconclusive"


def test_material_figures_and_tables_are_integrated(publication, tmp_path):
    analysis = tmp_path / "analysis"
    data = empty_data(publication, analysis)
    bindings = {
        case_id: install_material_case(publication, data, analysis, case_id)
        for case_id in ("R02", "R06", "R16", "R17")
    }
    install_material_science(publication, data, bindings)
    output = tmp_path / "publication"

    products = publication.render_products(data, output)

    for relative in (
        "figures/fig08_causal_mechanism.pdf",
        "figures/fig09_resolution_curves.pdf",
        "figures/fig10_hyperbolicity_coverage.pdf",
        "tables/numerical_health_provenance.csv",
        "tables/numerical_health_provenance.tex",
        "tables/primary_full_window_scalars.csv",
        "tables/primary_full_window_scalars.tex",
        "tables/hyperbolicity_all_snapshot_coverage.csv",
        "tables/signed_lf_cap_work_ledger.csv",
        "tables/mks24_panel_dispositions.csv",
        "tables/lineage_dispositions.csv",
        "tables/coherent_direction_mechanism.csv",
    ):
        assert output / relative in products
        assert (output / relative).is_file()
    assert publication.reviewed_pair_effect_rows(data, "R02", "R06")[0][
        "standardized_effect"
    ] == pytest.approx(-1.5)
    assert publication.normalized_history_series(data.cases["R02"], "c_b2")[1][
        -1
    ] == pytest.approx(0.05)
    assert publication.normalized_resolution_spectrum(data, "R17", "velocity")
    assert publication.resolution_alignment_curve(data, "R16")
    health = (output / "tables/numerical_health_provenance.csv").read_text(
        encoding="utf-8"
    )
    scalar = (output / "tables/primary_full_window_scalars.csv").read_text(
        encoding="utf-8"
    )
    hyperbolicity_header = (
        output / "tables/hyperbolicity_all_snapshot_coverage.csv"
    ).read_text(encoding="utf-8").splitlines()[0].split(",")
    assert "floor_margin" not in health
    assert "mass_relative_drift_maximum" in health
    assert "strict_failure" in health
    assert "effective_sample_count" in scalar
    assert "experiment_scope" in hyperbolicity_header
    assert "formula_id" in hyperbolicity_header
    assert "formula_disposition_family" in hyperbolicity_header
    assert "formula_provenance" in hyperbolicity_header
    assert "formula_executable_compatibility" in hyperbolicity_header
    assert "legacy_implementation_disposition" in hyperbolicity_header
    assert "literature_correct_disposition" in hyperbolicity_header
    assert "strict_hyperbolic_claim_status" in hyperbolicity_header
    assert "strict_hyperbolic_claim_eligible" in hyperbolicity_header
    assert "strict_hyperbolic_claim_reason" in hyperbolicity_header
    assert "audit_state_scope" in hyperbolicity_header
    assert "audit_direction_scope" in hyperbolicity_header
    assert "claim_scope" not in hyperbolicity_header


def test_all_snapshot_hyperbolicity_requires_authenticated_exact_coverage(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    _, rank = write_all_snapshot_hyperbolicity_evidence(
        publication, analysis, tmp_path / "hyper",
        formula_id=None,
        aggregates=(
            {
                "evaluated": 870580224,
                "negative": 0,
                "nonfinite_discriminant": 0,
                "minimum": 0.25,
            },
            {
                "evaluated": 870580224,
                "negative": 739552,
                "nonfinite_discriminant": 0,
                "minimum": -16.118885826958632,
            },
        ),
    )

    data = publication.discover_data(analysis, [tmp_path / "hyper"])
    row = {
        value["case_id"]: value
        for value in publication.hyperbolicity_coverage_rows(data)
    }["R02"]

    assert row["coverage_result"] == "pass"
    assert row["numerical_result"] == "negative"
    assert row["snapshot_policy"] == "all"
    assert row["complete_retained_snapshot_count"] == 2
    assert row["selected_snapshot_count"] == 2
    assert row["audited_snapshot_count"] == 2
    assert row["negative_discriminant_count"] == 739552
    assert row["nonfinite_discriminant_count"] == 0
    assert row["cell_direction_evaluations"] == 1741160448
    assert row["minimum_discriminant"] == -16.118885826958632
    assert row["experiment_scope"] == "standard"
    assert row["formula_id"] == "qualified-legacy"
    assert row["formula_disposition_family"] == "legacy_implementation"
    assert row["formula_id_provenance"] == "authenticated_legacy_default_contract"
    assert row["formula_provenance_status"] == "authenticated"
    assert row["executable_formula_id"] is None
    assert row["formula_executable_compatibility"] == "inconclusive"
    assert row["legacy_implementation_disposition"] == "negative"
    assert row["literature_correct_disposition"] == "inconclusive"
    assert (
        row["strict_hyperbolic_claim_status"]
        == "inconclusive_legacy_implementation"
    )
    assert row["strict_hyperbolic_claim_eligible"] is None
    assert "neither supports nor refutes" in row["strict_hyperbolic_claim_reason"]
    assert row["selection_provenance"] == "authenticated"
    assert row["result_provenance"] == "authenticated"
    assert {
        value["case_id"]: value
        for value in publication.hyperbolicity_coverage_rows(data)
    }["R06"]["coverage_result"] == "not_applicable"
    csv_path, _ = publication.write_table(
        tmp_path / "tables",
        "exact_hyperbolicity",
        ["cell_direction_evaluations", "minimum_discriminant"],
        [row],
    )
    assert csv_path.read_text(encoding="utf-8").splitlines()[1] == (
        "1741160448,-16.118885826958632"
    )

    rank.write_text("changed after authenticated audit\n", encoding="utf-8")
    stale = {
        value["case_id"]: value
        for value in publication.hyperbolicity_coverage_rows(data)
    }["R02"]
    assert stale["coverage_result"] == "inconclusive"
    assert stale["numerical_result"] == "inconclusive"
    assert stale["formula_id"] is None
    assert stale["legacy_implementation_disposition"] == "inconclusive"
    assert stale["literature_correct_disposition"] == "inconclusive"
    assert stale["strict_hyperbolic_claim_status"] == "inconclusive"
    assert stale["strict_hyperbolic_claim_eligible"] is None
    assert stale["selection_provenance"] == "inconclusive"
    assert stale["result_provenance"] == "inconclusive"


def test_directional_discriminant_dispositions_never_overclaim_strict_hyperbolicity(
    publication,
):
    hyperbolic = publication.snapshot_hyperbolicity_summary([{
        "aggregate": {
            "evaluated": 10,
            "negative": 0,
            "nonfinite_discriminant": 0,
            "minimum": 0.125,
        }
    }])
    nonfinite = publication.snapshot_hyperbolicity_summary([{
        "aggregate": {
            "evaluated": 10,
            "negative": 2,
            "nonfinite_discriminant": 1,
            "minimum": -4.0,
        }
    }])

    assert hyperbolic["result"] == "nonnegative_discriminant"
    assert nonfinite["result"] == "nonfinite"
    assert publication.strict_hyperbolic_claim_summary(
        "standard",
        "pass",
        "nonnegative_discriminant",
        "literature-correct",
        "authenticated",
        "compatible",
    ) == {
        "status": "inconclusive_retained_coordinate_discriminant_scope",
        "eligible": None,
        "reason": (
            "compatible authenticated literature-correct discriminants are "
            "nonnegative only for retained cell-centered snapshots in the three "
            "coordinate-normal directions; this does not test reconstructed faces, "
            "intermediate states, oblique directions, full or strict hyperbolicity, "
            "or absence of sqrt(|D|) fallback"
        ),
    }
    restricted = publication.strict_hyperbolic_claim_summary(
        "restricted",
        "pass",
        "nonnegative_discriminant",
        "literature-correct",
        "authenticated",
        "compatible",
    )
    assert restricted["status"] == "excluded_experiment_scope"
    assert restricted["eligible"] is False
    legacy = publication.strict_hyperbolic_claim_summary(
        "standard",
        "pass",
        "negative",
        "qualified-legacy",
        "authenticated",
        "compatible",
    )
    assert legacy["status"] == "inconclusive_legacy_implementation"
    assert legacy["eligible"] is None
    incompatible = publication.strict_hyperbolic_claim_summary(
        "standard",
        "pass",
        "nonnegative_discriminant",
        "literature-correct",
        "authenticated",
        "incompatible",
    )
    assert incompatible["status"] == "excluded_formula_executable_incompatible"
    assert incompatible["eligible"] is False
    no_compatibility = publication.strict_hyperbolic_claim_summary(
        "standard",
        "pass",
        "nonnegative_discriminant",
        "literature-correct",
        "authenticated",
        "inconclusive",
    )
    assert (
        no_compatibility["status"]
        == "inconclusive_formula_executable_compatibility"
    )
    assert no_compatibility["eligible"] is None


def test_formula_identity_routes_dispositions_and_gates_strict_claim(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    write_all_snapshot_hyperbolicity_evidence(
        publication,
        analysis,
        tmp_path / "literature",
        formula_id="literature-correct",
        formula_provenance={"reference": "independently-qualified expression"},
        executable_formula_id="literature-correct",
    )

    data = publication.discover_data(analysis, [tmp_path / "literature"])
    rows = publication.hyperbolicity_coverage_rows(data)
    row = {
        value["case_id"]: value
        for value in rows
    }["R02"]

    assert row["numerical_result"] == "nonnegative_discriminant"
    assert row["legacy_implementation_disposition"] == "inconclusive"
    assert row["literature_correct_disposition"] == "nonnegative_discriminant"
    assert row["formula_executable_compatibility"] == "compatible"
    assert (
        row["formula_executable_compatibility_provenance"]
        == "derived_from_authenticated_formula_ids"
    )
    assert (
        row["strict_hyperbolic_claim_status"]
        == "inconclusive_retained_coordinate_discriminant_scope"
    )
    assert row["strict_hyperbolic_claim_eligible"] is None
    assert row["audit_state_scope"] == "retained_cell_centered_snapshots"
    assert row["audit_direction_scope"] == "three_coordinate_normal_directions"
    assert "reconstructed faces" in row["strict_hyperbolic_claim_reason"]
    assert "sqrt(|D|) fallback" in row["strict_hyperbolic_claim_reason"]
    assert all(
        value["strict_hyperbolic_claim_eligible"] is not True for value in rows
    )

    incompatible_analysis = tmp_path / "incompatible-analysis"
    write_all_snapshot_hyperbolicity_evidence(
        publication,
        incompatible_analysis,
        tmp_path / "incompatible",
        formula_id="literature-correct",
        formula_provenance={"reference": "independently-qualified expression"},
        executable_formula_id="qualified-legacy",
    )
    incompatible_data = publication.discover_data(
        incompatible_analysis, [tmp_path / "incompatible"]
    )
    incompatible_row = {
        value["case_id"]: value
        for value in publication.hyperbolicity_coverage_rows(incompatible_data)
    }["R02"]
    assert incompatible_row["numerical_result"] == "nonnegative_discriminant"
    assert (
        incompatible_row["literature_correct_disposition"]
        == "nonnegative_discriminant"
    )
    assert incompatible_row["formula_executable_compatibility"] == "incompatible"
    assert (
        incompatible_row["strict_hyperbolic_claim_status"]
        == "excluded_formula_executable_incompatible"
    )
    assert incompatible_row["strict_hyperbolic_claim_eligible"] is False

    known_legacy = publication.authenticated_formula_summary({
        "provenance": {
            "script_sha256": next(
                iter(publication.KNOWN_LEGACY_HYPERBOLICITY_AUDIT_SCRIPTS)
            ),
            "script_version": "1.0.0",
            "formula_reference": "qualified/source/src/eos/eos.hpp:95-99",
        },
        "snapshots": [{
            "aggregate": {
                "evaluated": 10,
                "negative": 1,
                "nonfinite_discriminant": 0,
                "minimum": -1.0,
            },
        }],
    })
    assert known_legacy["formula_id"] == "qualified-legacy"
    assert (
        known_legacy["formula_id_provenance"]
        == "authenticated_known_legacy_audit_script"
    )
    assert known_legacy["legacy_implementation_disposition"] == "negative"
    assert known_legacy["literature_correct_disposition"] == "inconclusive"


def test_formula_id_manifest_result_mismatch_fails_closed(publication, tmp_path):
    analysis = tmp_path / "analysis"
    manifest_path, _ = write_all_snapshot_hyperbolicity_evidence(
        publication, analysis, tmp_path / "hyper"
    )
    result_path = manifest_path.parent / "result.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["provenance"]["formula_id"] = "literature-correct"
    write_json(result_path, result)
    (manifest_path.parent / "result.sha256").write_text(
        f"{sha256(result_path)}  result.json\n", encoding="utf-8"
    )

    data = publication.discover_data(analysis, [tmp_path / "hyper"])
    row = {
        value["case_id"]: value
        for value in publication.hyperbolicity_coverage_rows(data)
    }["R02"]

    assert row["coverage_result"] == "inconclusive"
    assert row["formula_id"] is None
    assert row["strict_hyperbolic_claim_eligible"] is None
    assert any(
        "result provenance differs from selected coverage" in warning
        for warning in data.ingestion_warnings
    )


def test_final_evidence_tables_preserve_semantics_and_fail_closed(
    publication, tmp_path
):
    analysis = tmp_path / "analysis"
    data = empty_data(publication, analysis)
    bindings = {
        case_id: install_material_case(publication, data, analysis, case_id)
        for case_id in ("R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09")
    }
    install_material_science(publication, data, bindings)
    data.science_record["mks24"] = {
        "result": "pass",
        "panels": {
            "fig2b": {
                "result": "pass",
                "reason": "all admitted products passed",
                "products": [{
                    "product_id": "fig2b_active",
                    "case_id": "R02",
                    "source": "authenticated_direct_fast_recomputation",
                    "result": "pass",
                }],
            }
        },
    }
    data.acceptance_records = [{
        "_publication_evidence_validated": True,
        "record_type": "stage-i-scientific-campaign-evidence",
        "gates": [{
            "name": "panel:fig3external",
            "result": "blocked_out_of_scope",
            "reason": "panel is explicitly blocked or external",
        }],
    }]

    ledgers = {
        row["case_id"]: row for row in publication.signed_lf_cap_work_ledger_rows(data)
    }
    assert ledgers["R02"]["applied_heat_flux_parallel"] == pytest.approx(2.02)
    assert ledgers["R02"]["applied_heat_flux_perpendicular"] == pytest.approx(-0.52)
    assert ledgers["R02"]["applied_pressure_work_total"] == pytest.approx(-1.02)
    assert ledgers["R02"]["reconstructed_anisotropic_stress_power_mean"] == pytest.approx(
        -0.22
    )
    assert "not applied accounting" in ledgers["R02"]["semantics"]

    panels = publication.mks24_panel_disposition_rows(data)
    assert next(row for row in panels if row["panel"] == "fig2b")[
        "disposition"
    ] == "admitted"
    assert next(row for row in panels if row["panel"] == "fig3external")[
        "disposition"
    ] == "blocked_or_external"

    lineage = [
        row for row in publication.lineage_disposition_rows(data)
        if row["case_id"] == "R02"
    ]
    assert [row["disposition"] for row in lineage] == ["selected", "failed"]
    assert lineage[1]["run_exit_code"] == 143

    directions = {
        row["metric"]: row
        for row in publication.coherent_direction_mechanism_rows(data)
    }
    assert directions["applied_pressure_work_total"]["descriptive_direction"] == (
        "active_gt_passive"
    )
    assert directions["parallel_strain_rms_mean"]["descriptive_direction"] == (
        "active_lt_passive"
    )
    assert directions["pressure_balance_correlation"]["descriptive_direction"] == (
        "active_lt_passive"
    )
    assert directions[
        "pressure_balance_normalized_residual_variance"
    ]["descriptive_direction"] == "active_lt_passive"
    assert directions[
        "compressive_velocity_power_fraction"
    ]["descriptive_direction"] == "active_gt_passive"
    assert directions["c_b2_full_window_mean"]["descriptive_direction"] == "equal"
    assert directions[
        "reviewed_abs_dp_signed_standardized_active_minus_passive_effect"
    ][
        "descriptive_direction"
    ] == "inconclusive"
    assert directions["applied_pressure_work_total"]["inference_scope"] == (
        "descriptive_only_no_preregistered_pass_gate"
    )

    data.cases["R02"].history_paths["user"].write_text(
        "changed after binding\n", encoding="utf-8"
    )
    stale_direction = {
        row["metric"]: row
        for row in publication.coherent_direction_mechanism_rows(data)
    }["c_b2_full_window_mean"]
    assert stale_direction["descriptive_direction"] == "inconclusive"
    assert stale_direction["available_pair_count"] == 3

    diagnostics_path = analysis / "cases/R02/diagnostics.json"
    diagnostics_path.write_text("{}\n", encoding="utf-8")
    stale_ledger = {
        row["case_id"]: row for row in publication.signed_lf_cap_work_ledger_rows(data)
    }["R02"]
    assert stale_ledger["availability"] == "inconclusive"
    assert stale_ledger["applied_pressure_work_total"] is None
    lineage_path = analysis / "cases/R02/lineage.json"
    lineage_path.write_text("{}\n", encoding="utf-8")
    stale_lineage = [
        row for row in publication.lineage_disposition_rows(data)
        if row["case_id"] == "R02"
    ]
    assert len(stale_lineage) == 1
    assert stale_lineage[0]["disposition"] == "inconclusive"
