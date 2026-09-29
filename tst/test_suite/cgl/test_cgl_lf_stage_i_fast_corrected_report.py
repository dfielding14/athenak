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
            "backup_limiters = false\n"
            "mirror_limiter = true\n"
            "firehose_limiter = true\n"
            f"limiter_nu_coll = {'200.0' if case_id == 'R15' else '20.0'}\n"
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


def history_payload(times: list[float]) -> str:
    rows = "".join(f"{time:.16e} 1.0000000000000000e+00\n" for time in times)
    return f"# [1]=time [2]=value\n{rows}"


def install_r14_failed_replay_fixture(
    adapter,
    campaign: dict[str, object],
    monkeypatch,
    mutate=None,
    case_id: str = "R14",
) -> None:
    """Replace one terminal-case fixture with two exact failed replays."""

    fake_assemble(adapter, campaign, monkeypatch)
    original = adapter.report.assemble_fast_case
    config = campaign["config"]

    def assemble(root, source, output, assembled_case_id, case):
        record = original(root, source, output, assembled_case_id, case)
        if assembled_case_id != case_id:
            return record

        parent = Path(record["lineage"][0]["segment_dir"])
        parent_manifest_path = parent / "manifest/fast_run.json"
        parent_manifest = json.loads(parent_manifest_path.read_text(encoding="utf-8"))
        parent_time = 4.0
        failure_time = 4.25
        overrides = [adapter.TERMINAL_STRICT_OVERRIDE]
        parent_manifest.update(
            {
                "variant": adapter.TERMINAL_APPROVED_VARIANT,
                "command_line_overrides": overrides,
                "strict_admissibility": False,
                "continuation_policy": "finite_progress_complete_terminal_products",
                "runtime_segmentation_changes_physics": False,
                "job_id": "9000",
                "slurm_log": str(
                    config.corrected.root / "logs/slurm-fast/%x.%j.log"
                ),
            }
        )
        write_json(parent_manifest_path, parent_manifest)
        parent_output = parent / "output"
        parent_mhd = parent_output / "fixture_parent.mhd.hst"
        parent_user = parent_output / "fixture_parent.user.hst"
        parent_mhd.write_text(history_payload([0.0, parent_time]), encoding="utf-8")
        parent_user.write_text(history_payload([0.0, parent_time]), encoding="utf-8")
        (parent / "manifest/run_exit_code").write_text("0\n", encoding="utf-8")
        restart = parent_output / "rst/rank_00000000/fixture.00001.rst"
        restart.parent.mkdir(parents=True)
        restart.write_text(
            f"<time>\ntime = {parent_time:.16e}\n<par_end>\n",
            encoding="utf-8",
        )
        restart_sha = sha256(restart)

        attempts = []
        logs = []
        for sequence, job_id in ((1, "9001"), (2, "9002")):
            segment = (
                config.corrected.run_root
                / case_id
                / f"fast_s{sequence:03d}_t4_to_t10"
            )
            manifest_path = segment / "manifest/fast_run.json"
            output_dir = segment / "output"
            output_dir.mkdir(parents=True)
            manifest = dict(parent_manifest)
            manifest.update(
                {
                    "sequence": sequence,
                    "run_dir": str(segment.resolve()),
                    "output_dir": str(output_dir.resolve()),
                    "run_basename": f"fixture_{case_id}_s{sequence:03d}",
                    "start_time": parent_time,
                    "restart": str(restart.resolve()),
                    "restart_sha256": restart_sha,
                    "job_id": job_id,
                }
            )
            write_json(manifest_path, manifest)
            (segment / "manifest/run_exit_code").write_text("1\n", encoding="utf-8")
            mhd = output_dir / f"fixture_{case_id}_s{sequence:03d}.mhd.hst"
            user = output_dir / f"fixture_{case_id}_s{sequence:03d}.user.hst"
            mhd.write_text(
                history_payload([parent_time, failure_time]), encoding="utf-8"
            )
            user.write_text(
                history_payload([parent_time, failure_time]), encoding="utf-8"
            )
            log = (
                config.corrected.root
                / "logs/slurm-fast"
                / f"cglc_{case_id}_s{sequence:03d}.{job_id}.log"
            )
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(
                f"time={failure_time:.16e}\n{adapter.TERMINAL_FATAL_SIGNATURE}\n",
                encoding="utf-8",
            )
            attempts.append(
                {
                    "sequence": sequence,
                    "segment": segment,
                    "manifest": manifest_path,
                    "mhd": mhd,
                    "user": user,
                    "log": log,
                }
            )
            logs.append(log)

        selected = attempts[-1]
        parent_lineage = record["lineage"][0]
        parent_lineage.update(
            {
                "manifest": {
                    "path": str(parent_manifest_path.resolve()),
                    "sha256": sha256(parent_manifest_path),
                },
                "observed_final_time": parent_time,
                "state": "exited_success_partial",
                "run_exit_code": 0,
                "variant": parent_manifest["variant"],
                "command_line_overrides": overrides,
                "mhd_history": str(parent_mhd.resolve()),
                "user_history": str(parent_user.resolve()),
            }
        )
        selected_manifest = json.loads(
            selected["manifest"].read_text(encoding="utf-8")
        )
        selected_lineage = {
            "kind": "fast",
            "segment_dir": str(selected["segment"].resolve()),
            "source_root": str(config.corrected.run_root.resolve()),
            "input": selected_manifest["input"],
            "input_sha256": selected_manifest["input_sha256"],
            "matrix_sha256": selected_manifest["matrix_sha256"],
            "executable": selected_manifest["executable"],
            "executable_sha256": selected_manifest["executable_sha256"],
            "manifest": {
                "path": str(selected["manifest"].resolve()),
                "sha256": sha256(selected["manifest"]),
            },
            "observed_final_time": failure_time,
            "state": "failed",
            "run_exit_code": 1,
            "restart": str(restart.resolve()),
            "restart_sha256": restart_sha,
            "variant": selected_manifest["variant"],
            "command_line_overrides": overrides,
            "mhd_history": str(selected["mhd"].resolve()),
            "user_history": str(selected["user"].resolve()),
        }
        record["lineage"] = [parent_lineage, selected_lineage]
        record["status"] = "failed_partial"
        record["final_time"] = failure_time
        record["model_choices"] = adapter.report.model_choices_for_input(
            config.corrected.source / case["input"], overrides
        )
        history_dir = output / f"cases/{case_id}/history"
        history_dir.mkdir(parents=True, exist_ok=True)
        merged_mhd = history_dir / f"{case['name']}.mhd.hst"
        merged_user = history_dir / f"{case['name']}.user.hst"
        for path in (merged_mhd, merged_user):
            path.write_text(
                history_payload([0.0, parent_time, failure_time]), encoding="utf-8"
            )
        record["histories"] = {
            "mhd": {
                "available": True,
                "path": str(merged_mhd.resolve()),
                "binding": adapter.report.artifact_binding(merged_mhd),
                "time_final": failure_time,
                "errors": [],
            },
            "user": {
                "available": True,
                "path": str(merged_user.resolve()),
                "binding": adapter.report.artifact_binding(merged_user),
                "time_final": failure_time,
                "errors": [],
            },
        }
        context = {
            "record": record,
            "parent_restart": restart,
            "attempts": attempts,
            "logs": logs,
            "failure_time": failure_time,
        }
        if mutate is not None:
            mutate(context)
        return record

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


def test_selected_attempt_sequences_allow_preserved_failed_attempt_gaps(adapter):
    assert adapter.validate_selected_attempt_sequences(
        [{"sequence": 0}, {"sequence": 2}, {"sequence": 5}], "R14"
    ) == [0, 2, 5]
    with pytest.raises(
        adapter.CompositeReportError, match="not strictly increasing"
    ):
        adapter.validate_selected_attempt_sequences(
            [{"sequence": 0}, {"sequence": 0}], "R14"
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


@pytest.mark.parametrize("status", ["partial", "failed_partial"])
def test_final_composite_rejects_non_r14_partial_case(
    adapter, campaign, monkeypatch, status
):
    fake_assemble(adapter, campaign, monkeypatch)
    original = adapter.report.assemble_fast_case

    def partial(root, source, output, case_id, case):
        record = original(root, source, output, case_id, case)
        if case_id == "R17":
            record["status"] = status
            record["final_time"] = 9.5
        return record

    monkeypatch.setattr(adapter.report, "assemble_fast_case", partial)
    with pytest.raises(adapter.CompositeReportError, match="R17 is not complete"):
        adapter.assemble_composite(campaign["config"], campaign["output"])
    assert not campaign["output"].exists()


def diverge_first_replay_restart(context):
    original = context["parent_restart"]
    replacement = original.with_name("different-checkpoint.rst")
    replacement.write_text(
        "<time>\ntime = 4.0000000000000000e+00\n<par_end>\n",
        encoding="utf-8",
    )
    manifest_path = context["attempts"][0]["manifest"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["restart"] = str(replacement.resolve())
    manifest["restart_sha256"] = sha256(replacement)
    write_json(manifest_path, manifest)


def test_r14_reproducible_failure_is_stamped_and_revalidated(
    adapter, campaign, monkeypatch
):
    install_r14_failed_replay_fixture(adapter, campaign, monkeypatch)

    inventory_path = adapter.assemble_composite(
        campaign["config"], campaign["output"]
    )
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    disposition = inventory["cases"]["R14"]["terminal_disposition"]

    assert disposition["disposition"] == (
        "reproducible_finite_time_model_runtime_failure"
    )
    assert disposition["attempt_count"] == 2
    assert [attempt["sequence"] for attempt in disposition["attempts"]] == [1, 2]
    assert disposition["physics"]["cgl_lf_strict_admissibility"] is False
    assert disposition["physics"]["variant"] == adapter.R14_APPROVED_VARIANT
    assert disposition["physics"]["command_line_overrides"] == (
        adapter.R14_APPROVED_OVERRIDES
    )
    assert disposition["physics"]["backup_limiters"] is False
    assert disposition["failure_time"] == pytest.approx(4.25, abs=1.0e-12)
    assert disposition["replay_history_consensus"]["mhd_sha256"]
    assert disposition["replay_history_consensus"]["user_sha256"]
    retained = json.loads(
        (
            campaign["output"] / "cases/R14/lineage.json"
        ).read_text(encoding="utf-8")
    )
    assert retained["terminal_disposition"] == disposition
    assert adapter.validate_composite_inventory(
        campaign["config"], campaign["output"]
    )["result"] == "pass"


def test_r15_reproducible_failure_uses_its_exact_limiter_contract(
    adapter, campaign, monkeypatch
):
    install_r14_failed_replay_fixture(
        adapter, campaign, monkeypatch, case_id="R15"
    )

    inventory_path = adapter.assemble_composite(
        campaign["config"], campaign["output"]
    )
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    disposition = inventory["cases"]["R15"]["terminal_disposition"]

    assert disposition["case_id"] == "R15"
    assert disposition["attempt_count"] == 2
    assert disposition["physics"]["limiter_nu_coll"] == 200.0
    assert disposition["physics"]["command_line_overrides"] == (
        adapter.TERMINAL_APPROVED_OVERRIDES
    )
    assert adapter.validate_composite_inventory(
        campaign["config"], campaign["output"]
    )["result"] == "pass"


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda context: context["logs"][0].write_text(
                "different fatal event\n", encoding="utf-8"
            ),
            "at least two distinct authenticated failed replay attempts",
        ),
        (
            lambda context: context["attempts"][0]["mhd"].write_text(
                history_payload([4.0, 4.250001]), encoding="utf-8"
            ),
            "at least two distinct authenticated failed replay attempts",
        ),
        (
            lambda context: context["attempts"][0]["mhd"].write_text(
                history_payload([4.0, 4.1, 4.25]), encoding="utf-8"
            ),
            "histories are not byte-identical",
        ),
        (
            lambda context: context["attempts"][0]["user"].write_text(
                history_payload([4.0, 4.1, 4.25]), encoding="utf-8"
            ),
            "histories are not byte-identical",
        ),
        (
            diverge_first_replay_restart,
            "at least two distinct authenticated failed replay attempts",
        ),
        (
            lambda context: context["record"]["model_choices"].update(
                {"cgl_lf_strict_admissibility": "true"}
            ),
            "approved physics contract",
        ),
        (
            lambda context: context["record"]["model_choices"].update(
                {"backup_limiters": "true"}
            ),
            "approved physics contract",
        ),
        (
            lambda context: context["record"]["histories"]["mhd"].update(
                {"available": False}
            ),
            "retained mhd history is unavailable",
        ),
    ],
)
def test_r14_failure_disposition_fails_closed_on_incomplete_evidence(
    adapter, campaign, monkeypatch, mutation, message
):
    install_r14_failed_replay_fixture(
        adapter, campaign, monkeypatch, mutate=mutation
    )

    with pytest.raises(adapter.CompositeReportError, match=message):
        adapter.assemble_composite(campaign["config"], campaign["output"])
    assert not campaign["output"].exists()


def test_r14_retained_disposition_is_rederived_from_source_attempts(
    adapter, campaign, monkeypatch
):
    install_r14_failed_replay_fixture(adapter, campaign, monkeypatch)
    inventory_path = adapter.assemble_composite(
        campaign["config"], campaign["output"]
    )
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    inventory["cases"]["R14"]["terminal_disposition"]["failure_time"] = 9.0
    write_json(inventory_path, inventory)
    lineage_path = campaign["output"] / "cases/R14/lineage.json"
    lineage = json.loads(lineage_path.read_text(encoding="utf-8"))
    lineage["terminal_disposition"]["failure_time"] = 9.0
    write_json(lineage_path, lineage)

    with pytest.raises(
        adapter.CompositeReportError,
        match="terminal disposition differs from source attempts",
    ):
        adapter.validate_composite_inventory(
            campaign["config"], campaign["output"]
        )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda manifest: manifest.update({"variant": "unreviewed_variant"}),
        lambda manifest: manifest.update(
            {
                "command_line_overrides": [
                    "mhd/cgl_lf_strict_admissibility=false",
                    "mhd/backup_limiters=true",
                ]
            }
        ),
        lambda manifest: manifest.update(
            {"runtime_segmentation_changes_physics": True}
        ),
    ],
)
def test_r14_failure_rejects_unapproved_manifest_physics(
    adapter, campaign, monkeypatch, mutation
):
    def mutate(context):
        for attempt in context["attempts"]:
            manifest = json.loads(
                attempt["manifest"].read_text(encoding="utf-8")
            )
            mutation(manifest)
            write_json(attempt["manifest"], manifest)
        selected = context["record"]["lineage"][-1]
        selected_manifest = json.loads(
            context["attempts"][-1]["manifest"].read_text(encoding="utf-8")
        )
        selected["manifest"]["sha256"] = sha256(context["attempts"][-1]["manifest"])
        selected["variant"] = selected_manifest["variant"]
        selected["command_line_overrides"] = selected_manifest[
            "command_line_overrides"
        ]

    install_r14_failed_replay_fixture(
        adapter, campaign, monkeypatch, mutate=mutate
    )

    with pytest.raises(adapter.CompositeReportError):
        adapter.assemble_composite(campaign["config"], campaign["output"])


def corrected_verification_fixture(
    adapter, tmp_path: Path, *, extra_error: str | None = None
) -> tuple[Path, dict[str, object]]:
    output = tmp_path / "report"
    segment = tmp_path / "runs/R14/replay"
    terminal_error = "selected lineage segment 0 has nonzero run exit code: 1"
    incomplete_error = "case is not complete: failed_partial"
    errors = [terminal_error, incomplete_error]
    if extra_error is not None:
        errors.append(extra_error)
    disposition = {
        "record_type": "cgl_lf_stage_i_terminal_disposition",
        "case_id": "R14",
        "status": "failed_partial",
        "disposition": "reproducible_finite_time_model_runtime_failure",
        "selected_terminal_segment": str(segment.resolve()),
        "attempts": [{
            "segment": str(segment.resolve()),
            "run_exit_code": 1,
        }],
    }
    write_json(
        output / "inventory.json",
        {
            "cases": {
                "R14": {
                    "lineage": [{
                        "state": "failed",
                        "segment_dir": str(segment.resolve()),
                        "run_exit_code": 1,
                    }]
                }
            }
        },
    )
    write_json(
        output / "verify.base.json",
        {
            "result": "fail",
            "require_complete": True,
            "adapter": {"fixture": True},
            "cases": {"R14": {"errors": errors}},
            "errors": [f"R14: {error}" for error in errors],
            "warnings": [],
        },
    )
    return output, {"terminal_dispositions": {"R14": disposition}}


def test_corrected_verification_removes_only_authenticated_r14_errors(
    adapter, tmp_path
):
    output, validation = corrected_verification_fixture(adapter, tmp_path)

    record = adapter.corrected_verification(output, validation)

    assert record["result"] == "pass"
    assert record["errors"] == []
    assert record["cases"]["R14"]["errors"] == []
    assert record["record_type"] == (
        "cgl_lf_stage_i_corrected_composite_verification"
    )
    assert len(record["accepted_base_errors"]) == 2


def test_corrected_verification_preserves_unrelated_errors(adapter, tmp_path):
    output, validation = corrected_verification_fixture(
        adapter, tmp_path, extra_error="history authentication differs"
    )

    record = adapter.corrected_verification(output, validation)

    assert record["result"] == "fail"
    assert record["errors"] == ["R14: history authentication differs"]
    assert record["cases"]["R14"]["errors"] == [
        "history authentication differs"
    ]


def test_corrected_verification_accepts_complete_campaign_without_exception(
    adapter, tmp_path
):
    output = tmp_path / "report"
    write_json(output / "inventory.json", {"cases": {"R14": {}}})
    write_json(
        output / "verify.base.json",
        {
            "result": "pass",
            "require_complete": True,
            "adapter": {"fixture": True},
            "cases": {"R14": {"errors": []}},
            "errors": [],
            "warnings": [],
        },
    )

    record = adapter.corrected_verification(
        output, {"terminal_dispositions": {}}
    )

    assert record["result"] == "pass"
    assert record["terminal_dispositions"] == {}
    assert record["accepted_base_errors"] == []
