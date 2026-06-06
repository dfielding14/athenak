"""Focused tests for corrected/legacy composite Stage I report provenance."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from dataclasses import replace

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
ADAPTER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_corrected_report.py"


def load_adapter():
    name = "cgl_lf_stage_i_fast_corrected_report_tests"
    spec = importlib.util.spec_from_file_location(name, ADAPTER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def adapter():
    return load_adapter()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def commit_source(source: Path) -> str:
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    subprocess.run(["git", "-C", str(source), "config", "user.name", "Fixture"], check=True)
    subprocess.run(
        ["git", "-C", str(source), "config", "user.email", "fixture@example.invalid"],
        check=True,
    )
    subprocess.run(["git", "-C", str(source), "add", "."], check=True)
    subprocess.run(["git", "-C", str(source), "commit", "-qm", "fixture"], check=True)
    return subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()


@pytest.fixture
def campaign(adapter, tmp_path: Path) -> dict[str, object]:
    corrected_root = tmp_path / "corrected-campaign"
    legacy_root = tmp_path / "legacy-campaign"
    corrected_source = tmp_path / "corrected-source"
    legacy_source = tmp_path / "legacy-source"
    authority_dir = tmp_path / "authority"
    audit = authority_dir / "audit.py"
    audit.parent.mkdir()
    audit.write_text("fixture audit\n", encoding="utf-8")

    cases = []
    input_hashes = {}
    for case_id in adapter.ALL_CASES:
        passive = case_id in adapter.PASSIVE_CONTROLS
        relative = f"inputs/{case_id}.athinput"
        payload = (
            "<mhd>\n"
            f"passive = {'true' if passive else 'false'}\n"
            "cgl_lf_strict_admissibility = true\n"
        )
        target_sources = (corrected_source, legacy_source) if passive else (corrected_source,)
        for source in target_sources:
            path = source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
        cases.append({"id": case_id, "name": f"fixture_{case_id}", "input": relative})
        if passive:
            input_hashes[case_id] = sha256(corrected_source / relative)

    matrix = {"cases": cases}
    for source in (corrected_source, legacy_source):
        write_json(source / adapter.MATRIX_RELATIVE, matrix)
        rsolver = source / adapter.PASSIVE_RSOLVER_RELATIVE
        rsolver.parent.mkdir(parents=True, exist_ok=True)
        rsolver.write_bytes(b"prefix\n" + adapter.PASSIVE_BRANCH + b"\nsuffix\n")

    legacy_eos = legacy_source / adapter.EOS_RELATIVE
    corrected_eos = corrected_source / adapter.EOS_RELATIVE
    legacy_eos.parent.mkdir(parents=True, exist_ok=True)
    corrected_eos.parent.mkdir(parents=True, exist_ok=True)
    legacy_eos.write_bytes(b"prefix\n" + adapter.LEGACY_EOS_TERM + b"\nsuffix\n")
    corrected_eos.write_bytes(b"prefix\n" + adapter.CORRECTED_EOS_TERM + b"\nsuffix\n")

    corrected_revision = commit_source(corrected_source)
    legacy_revision = commit_source(legacy_source)
    corrected_executable = authority_dir / "corrected-athena"
    legacy_executable = authority_dir / "legacy-athena"
    corrected_executable.write_bytes(b"corrected executable\n")
    legacy_executable.write_bytes(b"legacy executable\n")
    for root in (corrected_root, legacy_root):
        for case_id in adapter.ALL_CASES:
            if (
                root == corrected_root
                and case_id in adapter.ACTIVE_CASES
                or root == legacy_root
                and case_id in adapter.PASSIVE_CONTROLS
            ):
                (root / adapter.RUNS_RELATIVE / case_id).mkdir(parents=True)

    corrected = adapter.ExecutionAuthority(
        campaign_id="fixture-corrected",
        evidence_class="corrected_active_production",
        root=corrected_root,
        source=corrected_source,
        source_revision=corrected_revision,
        executable=corrected_executable,
        executable_sha256=sha256(corrected_executable),
        matrix=corrected_source / adapter.MATRIX_RELATIVE,
        matrix_sha256=sha256(corrected_source / adapter.MATRIX_RELATIVE),
    )
    legacy = adapter.ExecutionAuthority(
        campaign_id="fixture-legacy",
        evidence_class="authenticated_legacy_passive_control",
        root=legacy_root,
        source=legacy_source,
        source_revision=legacy_revision,
        executable=legacy_executable,
        executable_sha256=sha256(legacy_executable),
        matrix=legacy_source / adapter.MATRIX_RELATIVE,
        matrix_sha256=sha256(legacy_source / adapter.MATRIX_RELATIVE),
    )
    identity_path = adapter.identity_tool.write_identity(
        "fixture-corrected",
        corrected_root,
        corrected_revision,
        {
            "audit": audit,
            "eos": corrected_eos,
            "executable": corrected_executable,
            "matrix": corrected.matrix,
        },
        [corrected.run_root / case_id for case_id in adapter.ACTIVE_CASES],
    )
    config = adapter.CompositeConfig(
        corrected=corrected,
        legacy=legacy,
        corrected_identity=identity_path,
        corrected_identity_sha256=sha256(identity_path),
        compatibility=adapter.PassiveCompatibility(
            rsolver_sha256=sha256(corrected_source / adapter.PASSIVE_RSOLVER_RELATIVE),
            legacy_eos_sha256=sha256(legacy_eos),
            corrected_eos_sha256=sha256(corrected_eos),
            input_sha256=tuple(sorted(input_hashes.items())),
        ),
    )
    return {
        "config": config,
        "cases": {case["id"]: case for case in cases},
        "output": tmp_path / "composite-output",
    }


def fake_assemble(adapter, campaign: dict[str, object], monkeypatch) -> None:
    config = campaign["config"]
    cases = campaign["cases"]

    def assemble(root, source, output, case_id, case):
        authority = (
            config.corrected if case_id in adapter.ACTIVE_CASES else config.legacy
        )
        assert Path(root) == authority.root
        assert Path(source) == authority.source
        segment = authority.run_root / case_id / "fast_s000_t0_to_t10"
        manifest_path = segment / "manifest/fast_run.json"
        manifest = {
            "case_id": case_id,
            "case_name": case["name"],
            "sequence": 0,
            "root": str(authority.root.resolve()),
            "run_dir": str(segment.resolve()),
            "output_dir": str((segment / "output").resolve()),
            "input": str((authority.source / case["input"]).resolve()),
            "input_sha256": sha256(authority.source / case["input"]),
            "matrix_sha256": authority.matrix_sha256,
            "executable": str(authority.executable.resolve()),
            "executable_sha256": authority.executable_sha256,
            "start_time": 0.0,
            "target_time": 10.0,
            "restart": None,
            "restart_sha256": None,
        }
        if case_id in adapter.ACTIVE_CASES:
            manifest.update(
                {
                    "campaign_id": config.corrected.campaign_id,
                    "campaign_root": str(config.corrected.root.resolve()),
                    "campaign_identity": str(config.corrected_identity.resolve()),
                    "campaign_identity_sha256": config.corrected_identity_sha256,
                    "source": str(config.corrected.source.resolve()),
                    "source_revision": config.corrected.source_revision,
                    "legacy_restart_permitted": False,
                }
            )
        write_json(manifest_path, manifest)
        (segment / "output").mkdir()
        return {
            "schema_version": 1,
            "case_id": case_id,
            "case_name": case["name"],
            "status": "complete",
            "final_time": 10.0,
            "errors": [],
            "warnings": [],
            "lineage_identities": {
                "input_sha256": [manifest["input_sha256"]],
                "matrix_sha256": [authority.matrix_sha256],
                "executable_sha256": [authority.executable_sha256],
            },
            "lineage": [
                {
                    "kind": "fast",
                    "segment_dir": str(segment.resolve()),
                    "source_root": str(authority.run_root.resolve()),
                    "input": manifest["input"],
                    "input_sha256": manifest["input_sha256"],
                    "matrix_sha256": authority.matrix_sha256,
                    "executable": manifest["executable"],
                    "executable_sha256": authority.executable_sha256,
                    "manifest": {
                        "path": str(manifest_path.resolve()),
                        "sha256": sha256(manifest_path),
                    },
                }
            ],
            "histories": {},
            "snapshots": {"path": str((output / "cases" / case_id / "snapshots.json"))},
        }

    monkeypatch.setattr(adapter.report, "assemble_fast_case", assemble)


def test_preflight_authenticates_explicit_passive_compatibility(adapter, campaign):
    validation = adapter.validate_authorities(campaign["config"])
    compatibility = validation["passive_compatibility"]

    assert compatibility["result"] == "pass"
    assert compatibility["scope"] == list(adapter.PASSIVE_CONTROLS)
    assert compatibility["disposition_basis"]["commit"] == adapter.COMPATIBILITY_COMMIT
    assert compatibility["disposition_basis"]["test"]["sha256"] == (
        adapter.COMPATIBILITY_TEST_SHA256
    )
    assert compatibility["disposition_basis"]["reviewed_result"] == {
        "result": "pass",
        "summary": "1 passed",
        "elapsed_seconds": 233.45,
    }
    assert set(compatibility["supporting_evidence"]["inputs"]) == set(
        adapter.PASSIVE_CONTROLS
    )


def test_composite_assembly_preserves_actual_source_executable_and_evidence_class(
    adapter, campaign, monkeypatch
):
    fake_assemble(adapter, campaign, monkeypatch)
    inventory_path = adapter.assemble_composite(campaign["config"], campaign["output"])
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    active = inventory["cases"]["R02"]
    passive = inventory["cases"]["R06"]

    assert active["evidence_class"] == "corrected_active_production"
    assert active["execution_source"]["revision"] == (
        campaign["config"].corrected.source_revision
    )
    assert active["execution_executable"]["sha256"] == (
        campaign["config"].corrected.executable_sha256
    )
    assert passive["evidence_class"] == "authenticated_legacy_passive_control"
    assert passive["execution_source"]["revision"] == campaign["config"].legacy.source_revision
    assert passive["execution_executable"]["sha256"] == (
        campaign["config"].legacy.executable_sha256
    )
    assert passive["execution_authority"]["compatibility_contract_id"] == (
        inventory["passive_compatibility"]["contract_id"]
    )
    assert not any(path.is_symlink() for path in campaign["output"].rglob("*"))
    assert adapter.validate_composite_inventory(campaign["config"], campaign["output"])[
        "result"
    ] == "pass"


def test_executable_hash_drift_fails_before_assembly(adapter, campaign, monkeypatch):
    calls = []
    monkeypatch.setattr(adapter.report, "assemble_fast_case", lambda *_args: calls.append(True))
    campaign["config"].legacy.executable.write_bytes(b"drifted executable\n")

    with pytest.raises(adapter.CompositeReportError, match="executable SHA-256 mismatch"):
        adapter.assemble_composite(campaign["config"], campaign["output"])
    assert calls == []
    assert not campaign["output"].exists()


def test_passive_control_cannot_be_relabelled_with_corrected_root(
    adapter, campaign, monkeypatch
):
    fake_assemble(adapter, campaign, monkeypatch)
    original = adapter.report.assemble_fast_case

    def relabelled(root, source, output, case_id, case):
        record = original(root, source, output, case_id, case)
        if case_id == "R06":
            manifest_path = Path(record["lineage"][0]["manifest"]["path"])
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["root"] = str(campaign["config"].corrected.root.resolve())
            write_json(manifest_path, manifest)
            record["lineage"][0]["manifest"]["sha256"] = sha256(manifest_path)
        return record

    monkeypatch.setattr(adapter.report, "assemble_fast_case", relabelled)
    with pytest.raises(adapter.CompositeReportError, match="R06 root mismatch"):
        adapter.assemble_composite(campaign["config"], campaign["output"])
    assert not campaign["output"].exists()


def test_identity_root_mismatch_fails_closed(adapter, campaign):
    identity = json.loads(
        campaign["config"].corrected_identity.read_text(encoding="utf-8")
    )
    identity["campaign_root"] = str(campaign["config"].legacy.root.resolve())
    write_json(campaign["config"].corrected_identity, identity)
    changed = adapter.CompositeConfig(
        corrected=campaign["config"].corrected,
        legacy=campaign["config"].legacy,
        corrected_identity=campaign["config"].corrected_identity,
        corrected_identity_sha256=sha256(campaign["config"].corrected_identity),
        compatibility=campaign["config"].compatibility,
    )

    with pytest.raises(adapter.CompositeReportError, match="identity"):
        adapter.validate_authorities(changed)


def test_existing_or_symlinked_output_is_never_replaced(adapter, campaign, tmp_path):
    output = campaign["output"]
    output.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (output / "legacy-control").symlink_to(outside, target_is_directory=True)

    with pytest.raises(adapter.CompositeReportError, match="forbidden symlink"):
        adapter.assemble_composite(campaign["config"], output)


def test_only_explicit_legacy_race_root_is_accepted(adapter, campaign):
    legacy = campaign["config"].legacy
    race_relative = Path("runs/legacy-passive-races")
    authorized = legacy.root / race_relative / "R06/fast_s000_t0_to_t10"
    unauthorized = legacy.root / "runs/legacy-passive-relaxed/R06/fast_s000_t0_to_t10"
    authorized.mkdir(parents=True)
    unauthorized.mkdir(parents=True)
    authority = replace(
        legacy, run_relatives=(*legacy.run_relatives, race_relative)
    )

    resolved, root = adapter.selected_authority_run_root(
        authority, authorized, "authorized race"
    )
    assert resolved == authorized.resolve()
    assert root == (legacy.root / race_relative).resolve()
    with pytest.raises(adapter.CompositeReportError, match="matched 0 authorized"):
        adapter.selected_authority_run_root(authority, unauthorized, "relaxed run")


def test_final_composite_rejects_partial_case(adapter, campaign, monkeypatch):
    fake_assemble(adapter, campaign, monkeypatch)
    original = adapter.report.assemble_fast_case

    def partial(root, source, output, case_id, case):
        record = original(root, source, output, case_id, case)
        if case_id == "R17":
            record["status"] = "partial"
            record["final_time"] = 9.5
        return record

    monkeypatch.setattr(adapter.report, "assemble_fast_case", partial)
    with pytest.raises(adapter.CompositeReportError, match="R17 is not complete"):
        adapter.assemble_composite(campaign["config"], campaign["output"])
    assert not campaign["output"].exists()
