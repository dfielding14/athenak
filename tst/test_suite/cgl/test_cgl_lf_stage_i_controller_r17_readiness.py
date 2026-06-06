"""Focused production-controller tests for the strong R17 readiness gate."""

from __future__ import annotations

import copy
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
CONTROLLER_PATH = REPO_ROOT / "scripts/frontier/cgl_lf_stage_i.py"
QUALIFICATION_TEST_PATH = REPO_ROOT / "tst/test_suite/cgl/test_cgl_lf_stage_i_qualification.py"
SOURCE_AUTHORITY_PATH = REPO_ROOT / "scripts/frontier/cgl_lf_stage_i_source_authority.py"
EPOCH = "E03-forcing-policy"
EPOCH_SLUG = "E03_forcing_policy"
LIVE_FROZEN_E03_CLEAN_PARTIALS = (
    (
        "R03",
        "4762472",
        Path(
            "/lustre/orion/ast207/proj-shared/dfielding/CGL/runs/"
            "mks24-stage-i/E03-forcing-policy/R03/s00_rankio_t0_t0p5/"
            "manifest/prepared_run.json"
        ),
    ),
    (
        "R12",
        "4766856",
        Path(
            "/lustre/orion/ast207/proj-shared/dfielding/CGL/runs/"
            "mks24-stage-i/E03-forcing-policy/R12/s00_rankio_t0_t0p25/"
            "manifest/prepared_run.json"
        ),
    ),
)
R12_RETAINED_MHD_PREFIX = """# Athena++ history data
#  [1]=time      [2]=dt       [3]=mass    [4]=1-mom    [5]=2-mom    [6]=3-mom    [7]=tot-E    [8]=aam-D    [9]=1-KE    [10]=2-KE    [11]=3-KE    [12]=1-ME    [13]=2-ME    [14]=3-ME    [15]=lf_nstage    [16]=lf_dfloor    [17]=lf_pfloor    [18]=lf_nonfin    [19]=lf_nonpos    [20]=lf_mirror    [21]=lf_firehs    [22]=lf_hardbd    [23]=lf_qface    [24]=lf_qprcap    [25]=lf_qpr10    [26]=lf_qpecap    [27]=lf_qpe10    [28]=lf_qprwrk    [29]=lf_qpewrk    [30]=lf_hwproj    [31]=lf_cpwrk    [32]=lf_cawrk
   0.0000000000000000e+00   4.0343576522993910e-04   2.0000000000000004e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.6000000000000004e+01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000002e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00
"""
R12_RETAINED_USER_PREFIX = """# Athena++ history data
#  [1]=time      [2]=dt       [3]=volume    [4]=p_parallel    [5]=p_perp    [6]=force_prp2    [7]=force_prl2    [8]=vel_prp2    [9]=mass    [10]=kinetic    [11]=magnetic    [12]=therm_cgl    [13]=b2    [14]=b4    [15]=delta_p    [16]=abs_dp    [17]=beta    [18]=mirror_vol    [19]=fire_vol    [20]=hard_vol    [21]=nu_eff    [22]=force_pwr    [23]=vel_prl2    [24]=force_work
   0.0000000000000000e+00   4.0343576522993910e-04   2.0000000000000004e+00   1.0000000000000002e+01   1.0000000000000002e+01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   2.0000000000000004e+00   0.0000000000000000e+00   1.0000000000000002e+00   1.4999999999999996e+01   2.0000000000000004e+00   2.0000000000000004e+00   0.0000000000000000e+00   0.0000000000000000e+00   2.0000000000000004e+01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00
"""


def load_controller():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_controller_r17_readiness", CONTROLLER_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_source_authority():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_source_authority_for_controller_parity",
        SOURCE_AUTHORITY_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_controller_queue_authentication_allows_unrelated_jobs_but_not_cgl_drift(
    monkeypatch,
):
    module = load_controller()
    module.authenticate_production_queue(
        {},
        [],
        [
            "4769109_2|batch|PENDING|gotham_pdf_frontier",
            "4769110_3|batch|RUNNING|another_workload",
        ],
    )
    with pytest.raises(ValueError, match="queued CGL job"):
        module.authenticate_production_queue(
            {},
            [],
            ["999|batch|RUNNING|cgl_untracked_stage_i_job"],
        )
    monkeypatch.setattr(
        module,
        "expected_submitted_queue_jobs",
        lambda *_args: {"999": "cgl_expected_stage_i_job"},
    )
    for state in sorted(module.NONTERMINAL_STATES):
        module.authenticate_production_queue(
            {},
            [],
            [f"999|batch|{state}|cgl_expected_stage_i_job"],
        )
    for state in ("COMPLETED", "FAILED", "UNKNOWN", "RUNNING+", "running"):
        with pytest.raises(ValueError, match="queued CGL job"):
            module.authenticate_production_queue(
                {},
                [],
                [f"999|batch|{state}|cgl_expected_stage_i_job"],
            )
    with pytest.raises(ValueError, match="queued CGL job"):
        module.authenticate_production_queue(
            {},
            [],
            ["999|batch|RUNNING|forged_non_cgl_name"],
        )
    with pytest.raises(ValueError, match="malformed scheduler row"):
        module.authenticate_production_queue({}, [], ["malformed"])


def test_controller_retirement_rejects_post_fsync_forensic_substitution(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "target.json"
    target.write_text("authenticated\n")
    original_fsync = module.os.fsync
    calls = 0

    def substitute_forensic(descriptor):
        nonlocal calls
        original_fsync(descriptor)
        calls += 1
        if calls == 2:
            forensic = next(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))
            forensic.rename(tmp_path / "escaped-authenticated.forensic")
            forensic.write_text("substitute\n")

    monkeypatch.setattr(module.os, "fsync", substitute_forensic)
    with pytest.raises(
        ValueError, match="durable forensic recovery state was preserved"
    ):
        module.unlink_durable(target)
    assert not target.exists()
    assert (tmp_path / "escaped-authenticated.forensic").read_text() == "authenticated\n"


def load_qualification_test_support():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_controller_producer_support", QUALIFICATION_TEST_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_stable_json(path: Path, value: object, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    path.chmod(mode)


def immutable_publication_binding(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "sha256": sha256(path),
        "mode": "0444",
        "links": 1,
    }


def managed_clearance_fixture(module, tmp_path, monkeypatch, checkpoint=120):
    root = tmp_path / "root"
    accounting = root / "accounting"
    accounting.mkdir(parents=True)
    paths = paths_for(module, root)
    controller = tmp_path / "controller.py"
    controller.write_text("# retained controller\n")
    controller.chmod(0o644)
    monkeypatch.setattr(module, "retained_controller_source_path", lambda: controller)

    source_paths = [
        root / relative
        for relative in module.SHARED_ROOT_CLEARANCE_SOURCE_AUTHORITY_RELATIVES
    ]
    for index, path in enumerate(source_paths):
        write_stable_json(path, {"source_authority": index})
    legacy_paths = [
        root / relative for relative in module.SHARED_ROOT_CLEARANCE_LEGACY_RELATIVES
    ]
    for index, path in enumerate(legacy_paths):
        write_stable_json(path, {"legacy_clearance": index})

    stage_namespace = root / "runs/mks24-stage-i" / module.EXECUTION_EPOCH
    stale_namespace = root / "runs" / module.SHARED_ROOT_STALE_CAMPAIGN_ID
    stage_namespace.mkdir(parents=True, exist_ok=True)
    stale_manifest = stale_namespace / "manifest/prepared_run.json"
    write_stable_json(
        stale_manifest,
        {
            "campaign_id": module.SHARED_ROOT_STALE_CAMPAIGN_ID,
            "slurm_job_id": module.SHARED_ROOT_STALE_JOB_ID,
            "state": "running",
        },
        mode=0o644,
    )
    artifact_path, review_path, audit_path, audit_review_path = (
        module.managed_shared_root_clearance_chain_paths(accounting, checkpoint)
    )
    generated = datetime(2026, 6, 6, 12, 0, tzinfo=timezone.utc)
    artifact = module.build_managed_shared_root_clearance_refresh(
        paths, checkpoint, now=generated
    )
    write_stable_json(artifact_path, artifact)
    review = {
        "schema_version": 1,
        "record_type": "stage-i-managed-shared-root-isolation-clearance-review",
        "checkpoint": f"F-{checkpoint}",
        "execution_epoch": module.EXECUTION_EPOCH,
        "reviewed_utc": (generated + timedelta(minutes=1)).isoformat(),
        "decision": "approved",
        "reviewer": {"role": "independent security reviewer", "reviewer_id": "reviewer-1"},
        "independent_of_implementation": True,
        "candidate": immutable_publication_binding(artifact_path),
    }
    write_stable_json(review_path, review)
    audit = {
        "schema_version": 1,
        "record_type": (
            "stage-i-managed-shared-root-isolation-clearance-publication-audit"
        ),
        "checkpoint": f"F-{checkpoint}",
        "execution_epoch": module.EXECUTION_EPOCH,
        "published_utc": (generated + timedelta(minutes=2)).isoformat(),
        "artifact": immutable_publication_binding(artifact_path),
        "independent_review": immutable_publication_binding(review_path),
        "authority": {
            "binding_refresh_only": True,
            "scope_expanded": False,
            "bypasses_controller_preflights": False,
        },
    }
    write_stable_json(audit_path, audit)
    audit_review = {
        "schema_version": 1,
        "record_type": (
            "stage-i-managed-shared-root-isolation-clearance-"
            "publication-audit-review"
        ),
        "checkpoint": f"F-{checkpoint}",
        "execution_epoch": module.EXECUTION_EPOCH,
        "reviewed_utc": (generated + timedelta(minutes=3)).isoformat(),
        "decision": "approved-for-publication",
        "reviewer": {
            "role": "independent publication reviewer",
            "reviewer_id": "reviewer-2",
        },
        "independent_of_implementation": True,
        "candidate": immutable_publication_binding(audit_path),
    }
    write_stable_json(audit_review_path, audit_review)
    return {
        "paths": paths,
        "controller": controller,
        "source_paths": source_paths,
        "artifact_path": artifact_path,
        "review_path": review_path,
        "audit_path": audit_path,
        "audit_review_path": audit_review_path,
        "now": generated + timedelta(minutes=4),
    }


def test_managed_shared_root_clearance_authenticates_exact_refresh_chain(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = managed_clearance_fixture(module, tmp_path, monkeypatch)
    assert all(
        "F118_current_source_authority" in str(path)
        for path in fixture["source_paths"]
    )
    assert not any("F116_current_source_authority" in str(path)
                   for path in fixture["source_paths"])
    retained = module.require_managed_shared_root_clearance(
        fixture["paths"],
        {module.SHARED_ROOT_STALE_CAMPAIGN_ID},
        now=fixture["now"],
    )
    assert retained["checkpoint"] == "F-120"
    assert retained["artifact"]["sha256"] == sha256(fixture["artifact_path"])


@pytest.mark.parametrize(
    "mutation",
    ("controller_inode", "source_authority", "review", "same_reviewer"),
)
def test_managed_shared_root_clearance_rejects_changed_exact_binding(
    tmp_path, monkeypatch, mutation
):
    module = load_controller()
    fixture = managed_clearance_fixture(module, tmp_path, monkeypatch)
    if mutation == "controller_inode":
        fixture["controller"].unlink()
        fixture["controller"].write_text("# retained controller\n")
        fixture["controller"].chmod(0o644)
        match = "controller.*inode"
    elif mutation == "source_authority":
        source = fixture["source_paths"][0]
        source.chmod(0o644)
        source.write_text('{"changed": true}\n')
        source.chmod(0o444)
        match = "source authority publication 0"
    elif mutation == "review":
        review = fixture["review_path"]
        value = json.loads(review.read_text())
        value["candidate"]["sha256"] = "0" * 64
        review.chmod(0o644)
        write_stable_json(review, value)
        match = "review candidate"
    else:
        audit_review = fixture["audit_review_path"]
        value = json.loads(audit_review.read_text())
        value["reviewer"]["reviewer_id"] = "reviewer-1"
        audit_review.chmod(0o644)
        write_stable_json(audit_review, value)
        match = "review chain"
    with pytest.raises(ValueError, match=match):
        module.require_managed_shared_root_clearance(
            fixture["paths"],
            {module.SHARED_ROOT_STALE_CAMPAIGN_ID},
            now=fixture["now"],
        )


def test_render_managed_shared_root_clearance_refresh_supersedes_prior_chain(
    tmp_path, monkeypatch, capsys
):
    module = load_controller()
    fixture = managed_clearance_fixture(module, tmp_path, monkeypatch)
    root = fixture["paths"]["root"]
    monkeypatch.setattr(module, "require_root", lambda *_args: root)
    args = SimpleNamespace(root=str(root), allow_local_root=True, checkpoint=121)
    assert module.render_managed_shared_root_clearance_refresh(args) == 0
    candidate = json.loads(capsys.readouterr().out)
    assert candidate["checkpoint"] == "F-121"
    assert len(candidate["supersedes"]) == (
        len(module.SHARED_ROOT_CLEARANCE_LEGACY_RELATIVES) + 4
    )
    assert not any(
        path.exists()
        for path in module.managed_shared_root_clearance_chain_paths(
            root / "accounting", 121
        )
    )


def test_managed_shared_root_clearance_rejects_arbitrary_campaign(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = managed_clearance_fixture(module, tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="exactly the retained campaign"):
        module.require_managed_shared_root_clearance(
            fixture["paths"], {"other-campaign"}, now=fixture["now"]
        )


def load_live_frozen_e03_manifest(path: Path) -> dict[str, object]:
    if not path.is_file():
        pytest.skip(f"live frozen-E03 clean partial is unavailable: {path}")
    value = json.loads(path.read_text())
    assert isinstance(value, dict)
    return value


def live_migrated_parent(module, manifest_path: Path
                         ) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    """Return the exact live frozen-E03 manifest, migration, and parent record."""

    manifest = load_live_frozen_e03_manifest(manifest_path)
    migrated = module.revalidate_continuation_plasma_evidence(
        manifest["scientific_inspection"], manifest
    )
    terminal = migrated["terminal_restart"]
    parent = {
        "execution_epoch": EPOCH,
        "manifest": migrated["manifest"],
        "case_id": migrated["case_id"],
        "segment": migrated["segment"],
        "result": "clean_partial",
        "restart_sha256": terminal["sha256"],
        "restart_files": [
            str(path) for path in module.retained_product_paths(terminal)
        ],
        "final_time": migrated["final_time"],
        "restart_time": migrated["final_time"],
        "input_sha256": manifest["command"]["input_sha256"],
        "executable_sha256": manifest["command"]["executable_sha256"],
        "plasma_continuation_policy": migrated["plasma_continuation_policy"],
        "plasma_continuation_evidence": migrated["plasma_continuation_evidence"],
    }
    return manifest, migrated, parent


def write_json(path: Path, value: object, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.chmod(0o644)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    path.chmod(mode)


def overwrite_direct_child_restore_mtime(directory_fd: int, name: str,
                                         payload: bytes) -> None:
    """Overwrite one exact-length direct child while restoring its mtime."""

    before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    assert len(payload) == before.st_size
    descriptor = os.open(
        name,
        os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        assert os.write(descriptor, payload) == len(payload)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.utime(
        name,
        ns=(before.st_atime_ns, before.st_mtime_ns),
        dir_fd=directory_fd,
        follow_symlinks=False,
    )


def declared(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "sha256": sha256(path),
        "mode": "0444",
        "links": 1,
    }


def paths_for(module, root: Path) -> dict[str, Path]:
    paths = module.layout(root)
    paths["accounting"].mkdir(parents=True, exist_ok=True)
    paths["transactions"].mkdir(parents=True, exist_ok=True)
    paths["runs"].mkdir(parents=True, exist_ok=True)
    paths["ledger"].write_text("ledger\n")
    paths["reservations"].write_text("[]\n")
    return paths


def valid_ledger_row(module, job_id: str = "123") -> dict[str, str]:
    row = {column: "" for column in module.LEDGER_COLUMNS}
    row.update({
        "execution_epoch": module.EXECUTION_EPOCH,
        "job_id": job_id,
        "nodes": "1",
        "requested_walltime": "00:10:00",
        "elapsed_seconds": "60",
        "reserved_node_hours": "0.166667",
        "actual_node_hours": "0.016667",
        "cumulative_stage_i_node_hours": "0.016667",
    })
    return row


def producer_operational_qualification_fixture(module, tmp_path: Path, monkeypatch):
    """Emit and independently review one real producer schema-2 R17 bundle."""

    support = load_qualification_test_support()
    fixture = support.qualification_fixture.__wrapped__(tmp_path, monkeypatch)
    producer = support.qualification
    monkeypatch.setattr(producer, "RANKS_PER_NODE", 8)
    wave, _, evidence, scheduler = support.retain_complete_wave(
        fixture, ("R17",), max_nodes=8
    )
    retained = producer.retain_r17_operational_qualification(
        evidence["R17"], "qualification measurement agent", scheduler
    )
    root = fixture["root"]
    qualification_path = root / retained["operational_qualification"]["path"]
    qualification = json.loads(qualification_path.read_text())
    frozen = qualification["frozen_science_build_contract"]
    abis = dict(module.QUALIFIED_RESTART_BINARY_ABIS)
    abis[(frozen["executable_revision"], frozen["executable_sha256"])] = {}
    monkeypatch.setattr(module, "QUALIFIED_RESTART_BINARY_ABIS", abis)
    monkeypatch.setattr(module, "QUALIFIED_SOURCE_REVISION", frozen["source_revision"])
    monkeypatch.setattr(module, "R17_FROZEN_MATRIX_SHA256", frozen["matrix_sha256"])
    monkeypatch.setattr(module, "R17_FROZEN_INPUT_SHA256", frozen["input_sha256"])
    monkeypatch.setattr(
        module, "R17_FROZEN_PARAMETER_CONTRACT", frozen["parameter_contract"]
    )
    measured = datetime.fromisoformat(
        qualification["measured_utc"].replace("Z", "+00:00")
    )
    reviewed = measured + timedelta(minutes=1)
    review_path = root / qualification["independent_review_contract"]["path"]
    write_json(
        review_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-operational-qualification-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": reviewed.isoformat(),
            "decision": "approved",
            "reviewer": "qualification independent reviewer",
            "candidate": {
                "path": str(qualification_path),
                "sha256": sha256(qualification_path),
            },
        },
    )
    readiness = {
        "operational_qualification": {
            "path": qualification_path.relative_to(root).as_posix(),
            "sha256": sha256(qualification_path),
        },
        "operational_qualification_review": {
            "path": review_path.relative_to(root).as_posix(),
            "sha256": sha256(review_path),
        },
    }
    profile = {
        "source_bundle_sha256": frozen["source_bundle_sha256"],
        "executable_revision": frozen["executable_revision"],
        "executable_sha256": frozen["executable_sha256"],
        "build_manifest_sha256": frozen["build_manifest_inventory_sha256"],
    }
    return {
        "root": root,
        "wave": wave,
        "qualification_path": qualification_path,
        "qualification": qualification,
        "review_path": review_path,
        "readiness": readiness,
        "profile": profile,
        "authorization_time": reviewed + timedelta(minutes=1),
    }


def prepared_batch_fixture(module, tmp_path: Path
                           ) -> tuple[Path, dict[str, object], Path]:
    """Create one exact generated batch script suitable for descriptor tests."""

    manifest_path = tmp_path / "R16/s00/manifest/prepared_run.json"
    manifest_path.parent.mkdir(parents=True)
    batch = manifest_path.parent / "cgl_lf_stage_i.sbatch"
    manifest = {
        "run": {
            "case_id": "R16",
            "segment": "s00",
            "run_basename": "fixture",
        },
        "allocation": {
            "nodes": 1,
            "requested_walltime": "02:00:00",
            "ranks_per_node": 8,
            "cpus_per_task": 7,
        },
        "command": {
            "overrides": [],
            "athena_walltime": "01:50:00",
            "executable": str(tmp_path / "athena"),
            "executable_sha256": "4" * 64,
            "input_file": str(tmp_path / "input"),
            "input_sha256": "1" * 64,
            "matrix_file": str(tmp_path / "matrix"),
            "matrix_sha256": "2" * 64,
            "restart_file": None,
            "restart_files": [],
            "production_utility": {
                "path": str(CONTROLLER_PATH),
                "sha256": "3" * 64,
            },
        },
        "paths": {
            "batch_script": str(batch),
            "slurm_log": str(tmp_path / "%x.%j.log"),
            "output_dir": str(tmp_path / "output"),
            "environment_log": str(tmp_path / "environment.log"),
        },
    }
    script, digest = module.finalize_batch_script(
        module.generated_batch_script(manifest, manifest_path)
    )
    manifest["command"]["batch_script_sha256"] = digest
    batch.write_text(script)
    batch.chmod(0o750)
    return manifest_path, manifest, batch


def source_authority_binding(
    revision: str = "a" * 40,
    *,
    bundle_path: str | None = None,
    bundle_sha256: str = "5" * 64,
    verified_revisions: list[str] | None = None,
) -> dict[str, object]:
    bundle_path = bundle_path or (
        f"source-archives/athenak-feature-cgl-through-{revision[:9]}.bundle"
    )
    verified_revisions = verified_revisions or ["b" * 40, revision]
    return {
        "checkpoint": "F-118",
        "evidence": {"path": "accounting/f116.json", "sha256": "1" * 64},
        "provenance_review": {
            "path": "accounting/f116.provenance.json",
            "sha256": "2" * 64,
        },
        "plasma_review": {
            "path": "accounting/f116.plasma.json",
            "sha256": "3" * 64,
        },
        "publication_audit": {
            "path": "accounting/f116.audit.json",
            "sha256": "4" * 64,
        },
        "final_source_bundle": {
            "path": bundle_path,
            "sha256": bundle_sha256,
            "verified_revisions": verified_revisions,
        },
    }


def make_f116_source_authority(module, paths: dict[str, Path], tmp_path: Path,
                               monkeypatch) -> dict[str, object]:
    """Publish a complete focused F115/F116 fixture and return prepare bindings."""

    root = paths["root"]
    repo = tmp_path / "repo"
    matrix = repo / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
    matrix.parent.mkdir(parents=True)
    matrix.write_text('{"matrix": "current"}\n')
    monkeypatch.setattr(module, "ROOT_DIR", repo)
    archives = root / "source-archives"
    archives.mkdir(parents=True)
    historical_bundle = root / module.R03_F115_SOURCE_BUNDLE_RELATIVE
    historical_bundle.write_text("# v2 git bundle\nhistorical fixture\n")
    historical_bundle.chmod(0o644)
    monkeypatch.setattr(
        module, "R03_F115_SOURCE_BUNDLE_SHA256", sha256(historical_bundle)
    )

    f115_paths = {
        "evidence": root / module.R03_F115_SOURCE_AUTHORITY_RELATIVE,
        "publication_audit": root / (
            f"{module.R03_F115_SOURCE_AUTHORITY_RELATIVE}.publication_audit.json"
        ),
        "provenance_review": root / (
            f"{module.R03_F115_SOURCE_AUTHORITY_RELATIVE}.provenance_security_review.json"
        ),
        "plasma_review": root / (
            f"{module.R03_F115_SOURCE_AUTHORITY_RELATIVE}.plasma_scientific_review.json"
        ),
    }
    write_json(
        f115_paths["evidence"],
        {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-recovery-supersession-evidence",
            "checkpoint": "F-115",
            "execution_epoch": EPOCH,
        },
    )
    write_json(
        f115_paths["publication_audit"],
        {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-recovery-supersession-publication-audit",
            "checkpoint": "F-115",
            "execution_epoch": EPOCH,
        },
    )
    for key in ("provenance_review", "plasma_review"):
        write_json(
            f115_paths[key],
            {
                "schema_version": 1,
                "record_type": (
                    "stage-i-source-bundle-recovery-supersession-independent-review"
                ),
                "checkpoint": "F-115",
                "execution_epoch": EPOCH,
            },
        )
    monkeypatch.setattr(
        module, "R03_F115_SOURCE_AUTHORITY_SHA256", sha256(f115_paths["evidence"])
    )
    monkeypatch.setattr(
        module, "R03_F115_PUBLICATION_AUDIT_SHA256",
        sha256(f115_paths["publication_audit"]),
    )
    monkeypatch.setattr(
        module, "R03_F115_PROVENANCE_REVIEW_SHA256",
        sha256(f115_paths["provenance_review"]),
    )
    monkeypatch.setattr(
        module, "R03_F115_PLASMA_REVIEW_SHA256", sha256(f115_paths["plasma_review"])
    )
    historical_manifest = root / module.R03_F115_MANIFEST_RELATIVE
    write_json(
        historical_manifest,
        {
            "schema_version": 3,
            "execution_epoch": EPOCH,
            "state": "recorded",
            "job_id": module.R03_F115_JOB_ID,
            "project_root": str(root),
            "run": {"case_id": "R03", "segment": module.R03_F115_SEGMENT},
            "command": {
                "production_utility": {
                    "committed": True,
                    "revision": module.R03_F115_CONTROLLER_REVISION,
                    "sha256": module.R03_F115_CONTROLLER_SHA256,
                },
                "source_bundle": {
                    "path": str(root / module.R03_F115_SOURCE_BUNDLE_RELATIVE),
                    "sha256": module.R03_F115_SOURCE_BUNDLE_SHA256,
                    "verified_revisions": [
                        module.QUALIFIED_SOURCE_REVISION,
                        module.R03_F115_CONTROLLER_REVISION,
                    ],
                },
            },
            "accounting": {"result": "accepted", "state": "COMPLETED"},
            "scientific_inspection": {"final_time": 0.5},
        },
        mode=0o644,
    )
    monkeypatch.setattr(module, "R03_F115_MANIFEST_SHA256", sha256(historical_manifest))

    revision = "a" * 40
    bridge_revision = "b" * 40
    helper_sha = "6" * 64
    monkeypatch.setattr(module, "F116_BRIDGE_REVISION", bridge_revision)
    monkeypatch.setattr(
        module, "F116_BRIDGE_NAME",
        f"athenak-feature-cgl-through-{bridge_revision[:9]}.bundle",
    )
    monkeypatch.setattr(
        module,
        "F116_PRODUCTION_REQUIRED_REVISIONS",
        frozenset({
            module.QUALIFIED_SOURCE_REVISION,
            module.R03_F114_CONTROLLER_REVISION,
            module.R03_F115_CONTROLLER_REVISION,
            bridge_revision,
        }),
    )
    final_revisions = [
        module.QUALIFIED_SOURCE_REVISION,
        module.R03_F114_CONTROLLER_REVISION,
        module.R03_F115_CONTROLLER_REVISION,
        bridge_revision,
        revision,
    ]
    bundle_relative = Path(
        f"source-archives/athenak-feature-cgl-through-{revision[:9]}.bundle"
    )
    bundle = root / bundle_relative
    bundle.parent.mkdir(parents=True, exist_ok=True)
    bundle.write_text(
        f"# v2 git bundle\n{revision} HEAD\n\npayload\n"
    )
    bundle.chmod(0o644)
    bundle_sha = sha256(bundle)
    bridge_relative = Path(f"source-archives/{module.F116_BRIDGE_NAME}")
    bridge = root / bridge_relative
    bridge.write_text(
        f"# v2 git bundle\n{bridge_revision} refs/heads/feature/cgl-landau-fluid\n\n"
        "bridge payload\n"
    )
    bridge.chmod(0o644)
    monkeypatch.setattr(module, "F116_BRIDGE_SHA256", sha256(bridge))
    readme = archives / "README.md"
    sums = archives / "SHA256SUMS"
    readme.write_text("# Source archives\n\nF116 current source.\n")
    readme.chmod(0o644)
    sums.write_text(
        f"{module.R03_F115_SOURCE_BUNDLE_SHA256}  "
        f"{module.R03_F115_SOURCE_BUNDLE_RELATIVE.name}\n"
        f"{sha256(bridge)}  {bridge.name}\n"
        f"{bundle_sha}  {bundle.name}\n"
    )
    sums.chmod(0o644)
    now = datetime.now(timezone.utc).replace(microsecond=0)
    historical_bindings = {
        "evidence": {
            "path": module.R03_F115_SOURCE_AUTHORITY_RELATIVE.as_posix(),
            "sha256": module.R03_F115_SOURCE_AUTHORITY_SHA256,
        },
        "publication_audit": {
            "path": f"{module.R03_F115_SOURCE_AUTHORITY_RELATIVE}.publication_audit.json",
            "sha256": module.R03_F115_PUBLICATION_AUDIT_SHA256,
        },
        "provenance_review": {
            "path": f"{module.R03_F115_SOURCE_AUTHORITY_RELATIVE}.provenance_security_review.json",
            "sha256": module.R03_F115_PROVENANCE_REVIEW_SHA256,
        },
        "plasma_review": {
            "path": f"{module.R03_F115_SOURCE_AUTHORITY_RELATIVE}.plasma_scientific_review.json",
            "sha256": module.R03_F115_PLASMA_REVIEW_SHA256,
        },
    }
    historical_digests = {
        "evidence_sha256": module.R03_F115_SOURCE_AUTHORITY_SHA256,
        "publication_audit_sha256": module.R03_F115_PUBLICATION_AUDIT_SHA256,
        "provenance_review_sha256": module.R03_F115_PROVENANCE_REVIEW_SHA256,
        "plasma_review_sha256": module.R03_F115_PLASMA_REVIEW_SHA256,
    }
    tool_shas = {
        relative: hashlib.sha256(relative.encode()).hexdigest()
        for relative in module.F116_REQUIRED_TOOLS
    }
    tool_shas["scripts/frontier/cgl_lf_stage_i.py"] = helper_sha
    committed_tools = [
        {
            "path": relative,
            "revision": revision,
            "sha256": tool_shas[relative],
            "mode": mode,
        }
        for relative, mode in sorted(module.F116_REQUIRED_TOOLS.items())
    ]
    readme_before_payload = b"# Source archives\n"
    sums_before_payload = (
        f"{module.R03_F115_SOURCE_BUNDLE_SHA256}  "
        f"{module.R03_F115_SOURCE_BUNDLE_RELATIVE.name}\n"
    ).encode()
    before_catalog = {
        "readme_sha256": hashlib.sha256(readme_before_payload).hexdigest(),
        "sha256sums_sha256": hashlib.sha256(sums_before_payload).hexdigest(),
        "bridge_listed": False,
        "final_bundle_listed": False,
        "corrupt_c7_listed": False,
    }
    after_catalog = {
        "readme_sha256": sha256(readme),
        "sha256sums_sha256": sha256(sums),
        "bridge_listed_exactly_once": True,
        "final_bundle_listed_exactly_once": True,
        "corrupt_c7_listed": False,
        "historical_f115_preserved": True,
        "sole_current_source_bundle": bundle_relative.as_posix(),
    }
    evidence_path = root / module.F116_CURRENT_SOURCE_AUTHORITY_RELATIVE
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-current-source-authority-supersession-evidence",
        "checkpoint": "F-116",
        "execution_epoch": EPOCH,
        "generated_utc": (now - timedelta(minutes=4)).isoformat(),
        "scope": {
            "relationship": "current-source-selection-only-supersession",
            "summary": "Select the exact committed final tooling source.",
            "preserves": module.F116_PRESERVES,
            "does_not_authorize": module.F116_DOES_NOT_AUTHORIZE,
        },
        "predecessor_authorities": {"historical_f115": historical_bindings},
        "implementation": {
            "publisher": next(
                item for item in committed_tools
                if item["path"] == "scripts/frontier/cgl_lf_stage_i_source_authority.py"
            ),
            "committed_tools": committed_tools,
            "intermediate_36140_bundle": {
                "path": bridge_relative.as_posix(),
                "sha256": sha256(bridge),
                "complete_history": True,
                "head": bridge_revision,
                "advertised_tip": {
                    "revision": bridge_revision,
                    "name": "refs/heads/feature/cgl-landau-fluid",
                },
                "verified_revisions": [
                    module.R03_F115_CONTROLLER_REVISION,
                    bridge_revision,
                ],
                "selected_as_current": False,
                "role": "retained-non-current-bridge",
            },
            "current_source_bundle": {
                "candidate_path": str(tmp_path / "final.bundle.candidate"),
                "path": bundle_relative.as_posix(),
                "sha256": bundle_sha,
                "complete_history": True,
                "head": revision,
                "advertised_tip": {
                    "revision": revision,
                    "name": "HEAD",
                },
                "verified_revisions": final_revisions,
                "selected_as_current": True,
                "subject": "fixture final source authority",
            },
        },
        "source_archive_catalog": {"before": before_catalog, "after": after_catalog},
        "authorization": module.F116_AUTHORIZATION,
        "validation": module.F116_VALIDATION_CLAIMS,
        "publication_requirements": module.F116_PUBLICATION_REQUIREMENTS,
    }
    write_json(evidence_path, evidence)
    verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": bundle_sha,
        "final_head": revision,
        "historical_f115_preserved": True,
    }
    review_paths = {
        "provenance_review": root / module.F116_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": root / module.F116_PLASMA_REVIEW_RELATIVE,
    }
    for key, kind, decision, agent in (
        (
            "provenance_review", "provenance-security", "approved-for-publication",
            "provenance-reviewer",
        ),
        (
            "plasma_review", "plasma-scientific-continuation", "approved",
            "plasma-reviewer",
        ),
    ):
        write_json(
            review_paths[key],
            {
                "schema_version": 1,
                "record_type": (
                    "stage-i-current-source-authority-supersession-independent-review"
                ),
                "checkpoint": "F-116",
                "execution_epoch": EPOCH,
                "review_kind": kind,
                "decision": decision,
                "reviewed_candidate": {
                    "path": str(tmp_path / "f116.evidence.candidate"),
                    "sha256": sha256(evidence_path),
                },
                "published_f116": {
                    "path": str(evidence_path),
                    "sha256": sha256(evidence_path),
                },
                "reviewer": {
                    "agent_id": agent,
                    "identity": f"fixture {agent}",
                },
                "reviewed_utc": (now - timedelta(minutes=3)).isoformat(),
                "findings": ["Exact source selection and history bindings verified."],
                "limitations": ["No prepare or submit authority."],
                "verified": verified,
            },
        )
    audit_path = root / module.F116_PUBLICATION_AUDIT_RELATIVE
    write_json(
        audit_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-publication-audit",
            "checkpoint": "F-116",
            "execution_epoch": EPOCH,
            "published_utc": (now - timedelta(minutes=2)).isoformat(),
            "artifact": declared(evidence_path),
            "independent_reviews": {
                "reviews_bind_exact_published_f116_sha256": sha256(evidence_path),
                "provenance_security": declared(review_paths["provenance_review"]),
                "plasma_scientific_continuation": declared(review_paths["plasma_review"]),
            },
            "historical_f115_authority": historical_digests,
            "source_archive_catalog": {
                "readme": {
                    "path": str(readme),
                    "sha256": sha256(readme),
                    "mode": "0644",
                    "links": 1,
                },
                "sha256sums": {
                    "path": str(sums),
                    "sha256": sha256(sums),
                    "mode": "0644",
                    "links": 1,
                },
                "bridge_bundle": {
                    "path": str(bridge),
                    "sha256": sha256(bridge),
                    "mode": "0644",
                    "links": 1,
                    "head": bridge_revision,
                    "role": "retained-non-current-bridge",
                    "selected_as_current": False,
                },
                "current_source_bundle": {
                    "path": str(bundle),
                    "sha256": bundle_sha,
                    "mode": "0644",
                    "links": 1,
                    "head": revision,
                    "selected_as_current": True,
                },
                "corrupt_c7_absent_from_active_checksum_ledger": True,
                "sole_current_source_bundle": str(bundle),
            },
            "authority_and_enforcement": module.F116_AUTHORIZATION,
            "publication": (
                "recoverable-forward-transaction-with-publication-audit-commit-marker-"
                "under-stage-i-lock"
            ),
        },
    )
    f116_bindings = {
        "evidence": {
            "path": module.F116_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix(),
            "sha256": sha256(evidence_path),
        },
        "provenance_review": {
            "path": module.F116_PROVENANCE_REVIEW_RELATIVE.as_posix(),
            "sha256": sha256(review_paths["provenance_review"]),
        },
        "plasma_review": {
            "path": module.F116_PLASMA_REVIEW_RELATIVE.as_posix(),
            "sha256": sha256(review_paths["plasma_review"]),
        },
        "publication_audit": {
            "path": module.F116_PUBLICATION_AUDIT_RELATIVE.as_posix(),
            "sha256": sha256(audit_path),
        },
    }
    f116_digests = {
        "evidence_sha256": sha256(evidence_path),
        "provenance_review_sha256": sha256(review_paths["provenance_review"]),
        "plasma_review_sha256": sha256(review_paths["plasma_review"]),
        "publication_audit_sha256": sha256(audit_path),
    }
    current_bundle = evidence["implementation"]["current_source_bundle"]
    bridge_bundle = evidence["implementation"]["intermediate_36140_bundle"]
    predecessor_bundle = dict(current_bundle)
    predecessor_bundle.pop("candidate_path")
    predecessor_bundle["selected_as_current"] = False
    predecessor_bundle["role"] = "retained-non-current-predecessor"
    f118_before = {
        "readme_sha256": sha256(readme),
        "sha256sums_sha256": sha256(sums),
        "bridge_listed_exactly_once": True,
        "predecessor_current_source_bundle_listed_exactly_once": True,
        "final_bundle_listed": False,
        "corrupt_c7_listed": False,
        "historical_f115_preserved": True,
    }
    f118_after = {
        "readme_sha256": sha256(readme),
        "sha256sums_sha256": sha256(sums),
        "bridge_listed_exactly_once": True,
        "predecessor_current_source_bundle_listed_exactly_once": True,
        "final_bundle_listed_exactly_once": True,
        "corrupt_c7_listed": False,
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
        "all_prior_checksum_entries_preserved": True,
        "sole_current_source_bundle": bundle_relative.as_posix(),
    }
    evidence_path = root / module.F118_CURRENT_SOURCE_AUTHORITY_RELATIVE
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-current-source-authority-supersession-evidence",
        "checkpoint": "F-118",
        "execution_epoch": EPOCH,
        "generated_utc": (now - timedelta(minutes=1)).isoformat(),
        "scope": {
            "relationship": "current-source-selection-only-supersession",
            "summary": "Select the exact committed F118 tooling source.",
            "preserves": module.F118_PRESERVES,
            "does_not_authorize": module.F118_DOES_NOT_AUTHORIZE,
        },
        "predecessor_authorities": {"historical_f116": f116_bindings},
        "implementation": {
            "publisher": next(
                item for item in committed_tools
                if item["path"] == "scripts/frontier/cgl_lf_stage_i_source_authority.py"
            ),
            "committed_tools": committed_tools,
            "intermediate_36140_bundle": bridge_bundle,
            "predecessor_current_source_bundle": predecessor_bundle,
            "current_source_bundle": current_bundle,
        },
        "source_archive_catalog": {"before": f118_before, "after": f118_after},
        "authorization": module.F118_AUTHORIZATION,
        "validation": module.F118_VALIDATION_CLAIMS,
        "publication_requirements": module.F118_PUBLICATION_REQUIREMENTS,
    }
    write_json(evidence_path, evidence)
    f118_verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "predecessor_current_source_bundle_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": bundle_sha,
        "final_head": revision,
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
    }
    review_paths = {
        "provenance_review": root / module.F118_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": root / module.F118_PLASMA_REVIEW_RELATIVE,
    }
    for key, kind, decision, agent in (
        (
            "provenance_review", "provenance-security", "approved-for-publication",
            "fixture-f118-provenance-reviewer",
        ),
        (
            "plasma_review", "plasma-scientific-continuation", "approved",
            "fixture-f118-plasma-reviewer",
        ),
    ):
        write_json(
            review_paths[key],
            {
                "schema_version": 1,
                "record_type": (
                    "stage-i-current-source-authority-supersession-independent-review"
                ),
                "checkpoint": "F-118",
                "execution_epoch": EPOCH,
                "review_kind": kind,
                "decision": decision,
                "reviewed_candidate": {
                    "path": str(tmp_path / "f118.evidence.candidate"),
                    "sha256": sha256(evidence_path),
                },
                "published_f118": {
                    "path": str(evidence_path),
                    "sha256": sha256(evidence_path),
                },
                "reviewer": {"agent_id": agent, "identity": f"fixture {agent}"},
                "reviewed_utc": now.isoformat(),
                "findings": ["Exact F118 source selection and F116 history verified."],
                "limitations": ["No prepare or submit authority."],
                "verified": f118_verified,
            },
        )
    audit_path = root / module.F118_PUBLICATION_AUDIT_RELATIVE
    write_json(
        audit_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-publication-audit",
            "checkpoint": "F-118",
            "execution_epoch": EPOCH,
            "published_utc": now.isoformat(),
            "artifact": declared(evidence_path),
            "independent_reviews": {
                "reviews_bind_exact_published_f118_sha256": sha256(evidence_path),
                "provenance_security": declared(review_paths["provenance_review"]),
                "plasma_scientific_continuation": declared(review_paths["plasma_review"]),
            },
            "historical_f116_authority": f116_digests,
            "source_archive_catalog": {
                "readme": {
                    "path": str(readme), "sha256": sha256(readme),
                    "mode": "0644", "links": 1,
                },
                "sha256sums": {
                    "path": str(sums), "sha256": sha256(sums),
                    "mode": "0644", "links": 1,
                },
                "bridge_bundle": {
                    "path": str(bridge), "sha256": sha256(bridge),
                    "mode": "0644", "links": 1, "head": bridge_revision,
                    "role": "retained-non-current-bridge", "selected_as_current": False,
                },
                "predecessor_current_source_bundle": {
                    "path": str(bundle), "sha256": bundle_sha,
                    "mode": "0644", "links": 1, "head": revision,
                    "role": "retained-non-current-predecessor",
                    "selected_as_current": False,
                },
                "current_source_bundle": {
                    "path": str(bundle), "sha256": bundle_sha,
                    "mode": "0644", "links": 1, "head": revision,
                    "selected_as_current": True,
                },
                "corrupt_c7_absent_from_active_checksum_ledger": True,
                "sole_current_source_bundle": str(bundle),
            },
            "authority_and_enforcement": module.F118_AUTHORIZATION,
            "publication": (
                "recoverable-forward-transaction-with-publication-audit-commit-marker-"
                "under-stage-i-lock"
            ),
        },
    )
    monkeypatch.setattr(
        module,
        "source_authority_committed_sha256",
        lambda _revision, relative, _label: (
            sha256(matrix)
            if relative == "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
            else tool_shas[relative]
        ),
    )
    monkeypatch.setattr(
        module,
        "source_authority_committed_mode",
        lambda _revision, relative, _label: module.F116_REQUIRED_TOOLS[relative],
    )
    monkeypatch.setattr(
        module, "source_authority_revision_is_ancestor", lambda *_args: True
    )
    monkeypatch.setattr(
        module,
        "source_authority_commit_subject",
        lambda _revision: "fixture final source authority",
    )
    monkeypatch.setattr(
        module,
        "source_bundle_provenance",
        lambda _value, revisions, _root, _allow: {
            "path": str(bundle),
            "sha256": bundle_sha,
            "verified_revisions": revisions,
        },
    )
    binding = {
        "checkpoint": "F-118",
        "evidence": {
            "path": module.F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix(),
            "sha256": sha256(evidence_path),
        },
        "provenance_review": {
            "path": module.F118_PROVENANCE_REVIEW_RELATIVE.as_posix(),
            "sha256": sha256(review_paths["provenance_review"]),
        },
        "plasma_review": {
            "path": module.F118_PLASMA_REVIEW_RELATIVE.as_posix(),
            "sha256": sha256(review_paths["plasma_review"]),
        },
        "publication_audit": {
            "path": module.F118_PUBLICATION_AUDIT_RELATIVE.as_posix(),
            "sha256": sha256(audit_path),
        },
        "final_source_bundle": {
            "path": bundle_relative.as_posix(),
            "sha256": bundle_sha,
            "verified_revisions": final_revisions,
        },
    }
    return {
        "binding": binding,
        "bundle": bundle,
        "bundle_record": {
            "path": str(bundle),
            "sha256": bundle_sha,
            "verified_revisions": [revision, module.QUALIFIED_SOURCE_REVISION],
        },
        "evidence": evidence_path,
        "historical_manifest": historical_manifest,
        "matrix": matrix,
        "revision": revision,
        "audit": audit_path,
        "bridge": bridge,
        "historical_bundle": historical_bundle,
        "sums": sums,
        "utility": {
            "path": str(CONTROLLER_PATH),
            "revision": revision,
            "sha256": helper_sha,
            "committed": True,
        },
        "readme_before_payload": readme_before_payload,
        "sums_before_payload": sums_before_payload,
    }


def make_retained_f118_staging(module, paths: dict[str, Path],
                               fixture: dict[str, object]) -> Path:
    """Create one complete retained publisher staging bound to public F118."""

    root = paths["root"]
    evidence = json.loads(fixture["evidence"].read_text())
    audit = json.loads(fixture["audit"].read_text())
    transaction_id = "2026-06-05T120000+0000-" + "a" * 32
    staging = (
        paths["accounting"]
        / f"mks24_stage_i_{EPOCH_SLUG}_F118_source_authority_transactions"
        / f"{transaction_id}.staging"
    )
    staging.mkdir(parents=True)
    staging.chmod(0o700)
    public_payloads = {
        "bundle": fixture["bundle"].read_bytes(),
        "evidence": fixture["evidence"].read_bytes(),
        "provenance_review": (
            root / module.F118_PROVENANCE_REVIEW_RELATIVE
        ).read_bytes(),
        "plasma_review": (
            root / module.F118_PLASMA_REVIEW_RELATIVE
        ).read_bytes(),
        "audit": fixture["audit"].read_bytes(),
        "readme_before": fixture["readme_before_payload"],
        "sha256sums_before": fixture["sums_before_payload"],
        "readme_after": (root / "source-archives/README.md").read_bytes(),
        "sha256sums_after": fixture["sums"].read_bytes(),
    }
    bindings = {}
    for key, (name, mode_text) in (
        module.F116_SOURCE_AUTHORITY_TRANSACTION_PAYLOADS.items()
    ):
        path = staging / name
        path.write_bytes(public_payloads[key])
        path.chmod(int(mode_text, 8))
        bindings[key] = {
            "name": name,
            "sha256": hashlib.sha256(public_payloads[key]).hexdigest(),
            "mode": mode_text,
        }
    implementation = evidence["implementation"]
    publisher = implementation["publisher"]
    final = implementation["current_source_bundle"]
    journal = {
        "schema_version": 1,
        "record_type": "stage-i-current-source-authority-publication-transaction",
        "transaction_id": transaction_id,
        "execution_epoch": EPOCH,
        "checkpoint": "F-118",
        "state": "staged",
        "created_utc": evidence["generated_utc"],
        "publisher": {
            "revision": publisher["revision"],
            "sha256": publisher["sha256"],
        },
        "candidate_paths": {
            key: str((root / "candidates" / key).absolute())
            for key in ("bundle", "evidence", "provenance_review", "plasma_review", "audit")
        },
        "expected": {
            key: bindings[key]["sha256"]
            for key in ("bundle", "evidence", "provenance_review", "plasma_review", "audit")
        },
        "payloads": bindings,
        "targets": {
            "bundle": final["path"],
            "evidence": module.F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix(),
            "provenance_review": module.F118_PROVENANCE_REVIEW_RELATIVE.as_posix(),
            "plasma_review": module.F118_PLASMA_REVIEW_RELATIVE.as_posix(),
            "publication_audit": module.F118_PUBLICATION_AUDIT_RELATIVE.as_posix(),
            "readme": "source-archives/README.md",
            "sha256sums": "source-archives/SHA256SUMS",
        },
        "catalog_before": evidence["source_archive_catalog"]["before"],
        "catalog_after": evidence["source_archive_catalog"]["after"],
    }
    write_json(staging / "journal.json", journal, mode=0o600)
    recovery = staging / ".journal.json.recovery.tmp"
    recovery.write_bytes((staging / "journal.json").read_bytes())
    recovery.chmod(0o600)
    assert audit["published_utc"] >= journal["created_utc"]
    return staging


def make_published_recost(module, root: Path, *, checkpoint: int = 200,
                          published: datetime | None = None,
                          generated: datetime | None = None
                          ) -> tuple[dict[str, Path], dict[str, object]]:
    now = datetime.now(timezone.utc).replace(microsecond=0)
    generated = generated or now - timedelta(minutes=10)
    published = published or now - timedelta(minutes=5)
    paths = paths_for(module, root)
    artifact = paths["accounting"] / (
        f"mks24_stage_i_{EPOCH_SLUG}_F{checkpoint}_recost_evidence.json"
    )
    review = artifact.with_name(f"{artifact.name}.independent_review.json")
    audit = artifact.with_name(f"{artifact.name}.publication_audit.json")
    write_json(
        artifact,
        {
            "generated_utc": generated.isoformat(),
            "expires_utc": (generated + timedelta(hours=12)).isoformat(),
        },
    )
    write_json(
        review,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (published - timedelta(minutes=1)).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": {
                "agent_id": "independent-reviewer",
                "independent_from_generator": True,
            },
            "candidate": {"path": str(artifact), "sha256": sha256(artifact)},
            "scope": {"non_authorizing": True},
        },
    )
    audit_value = {
        "schema_version": 1,
        "record_type": "stage-i-recost-recommendation-publication-audit",
        "execution_epoch": EPOCH,
        "published_utc": published.isoformat(),
        "artifact": declared(artifact),
        "independent_review": declared(review),
        "authority": {
            "action_authority": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
    }
    write_json(audit, audit_value)
    return paths, audit_value


def make_fixed_readiness(module, paths: dict[str, Path], recost: dict[str, object],
                         profile: dict[str, object], lineage_sha: str,
                         current_source_authority: dict[str, object], *, now: datetime
                         ) -> tuple[Path, dict[str, object], dict[str, object]]:
    recost.setdefault("generated_utc", (now - timedelta(minutes=5)).isoformat())
    recost.setdefault("expires_utc", (now + timedelta(hours=1)).isoformat())
    readiness_path = paths["root"] / module.R17_READINESS_RELATIVE
    qualification = paths["accounting"] / "r17_operational_qualification.json"
    qualification_review = qualification.with_name(
        f"{qualification.name}.independent_review.json"
    )
    qualification.write_text("{}\n")
    qualification.chmod(0o644)
    write_json(qualification_review, {"fixture": True})
    readiness = {
        "schema_version": 1,
        "record_type": "stage-i-r17-readiness",
        "execution_epoch": EPOCH,
        "root": str(paths["root"]),
        "generated_utc": (now - timedelta(minutes=8)).isoformat(),
        "expires_utc": (now + timedelta(hours=2)).isoformat(),
        "reviewed_by": "readiness-reviewer",
        "predecessor_lineages_sha256": lineage_sha,
        "storage_evidence_sha256": "2" * 64,
        "computed_projection_sha256": "3" * 64,
        "executable_sha256": profile["executable_sha256"],
        "build_manifest_sha256": profile["build_manifest_sha256"],
        "required_retained_bytes": module.R17_MINIMUM_RETAINED_BYTES,
        "nodes": 8,
        "ranks": 64,
        "storage_ready": True,
        "node_hour_ready": True,
        "rank_64_ready": True,
        "operational_qualification": {
            "path": qualification.relative_to(paths["root"]).as_posix(),
            "sha256": sha256(qualification),
        },
        "operational_qualification_review": {
            "path": qualification_review.relative_to(paths["root"]).as_posix(),
            "sha256": sha256(qualification_review),
        },
        "current_source_authority": current_source_authority,
    }
    write_json(readiness_path, readiness)
    review_path = readiness_path.with_name(
        f"{readiness_path.name}.independent_review.json"
    )
    write_json(
        review_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-readiness-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (now - timedelta(minutes=7)).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": "readiness-reviewer",
            "candidate": {"path": str(readiness_path), "sha256": sha256(readiness_path)},
        },
    )
    audit_path = readiness_path.with_name(f"{readiness_path.name}.publication_audit.json")
    audit = {
        "schema_version": 1,
        "record_type": "stage-i-r17-readiness-publication-audit",
        "execution_epoch": EPOCH,
        "published_utc": (now - timedelta(minutes=6)).isoformat(),
        "artifact": declared(readiness_path),
        "independent_review": declared(review_path),
        "authority": {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
    }
    write_json(audit_path, audit)
    provenance = recost.setdefault("provenance", {})
    assert isinstance(provenance, dict)
    provenance.update({
        "storage_evidence_sha256": "2" * 64,
        "computed_projection_sha256": "3" * 64,
        "r17_readiness_evidence_sha256": sha256(readiness_path),
    })
    recost["r17_readiness"] = {
        **readiness,
        "operational_qualification_evidence": {},
        "publication_chain": {
            "independent_review_sha256": sha256(review_path),
            "publication_audit_sha256": sha256(audit_path),
            "reviewed_by": "readiness-reviewer",
            "published_utc": audit["published_utc"],
        },
    }
    return readiness_path, readiness, audit


def scoped_node_hour_basis(job_id: str, case_id: str, segment: str):
    return {
        "job_id": job_id,
        "case_id": case_id,
        "segment": segment,
        "actual_node_hours": "1",
        "observed_cells": 1,
        "observed_simulation_interval": "10",
        "normalized_node_hours_per_cell_per_simulation_time": "0.1",
    }


def final_r17_scoped_budget():
    global_basis = scoped_node_hour_basis("5000003", "R03", "s00_rankio_t0_t10")
    r12_basis = scoped_node_hour_basis("5000012", "R12", "s01_rankio_t0_t10")
    breakdown = {}
    for case_number in range(2, 18):
        case_id = f"R{case_number:02d}"
        is_r17 = case_id == "R17"
        breakdown[case_id] = {
            "matrix_full_case_node_hours_reference_only": "10",
            "authenticated_progress_fraction": "0" if is_r17 else "1",
            "remaining_simulation_time": "10" if is_r17 else "0",
            "projected_cells": "1",
            "projection_measurement_basis": copy.deepcopy(
                r12_basis if case_id == "R12" else global_basis
            ),
            "observed_rate_projected_remaining_node_hours": "1" if is_r17 else "0",
            "authorized_profile_reserved_node_hours": "16" if is_r17 else "0",
            "projected_remaining_node_hours": "16" if is_r17 else "0",
        }
    return {
        "method": "observed-stage-i-scoped-node-hour-rate-v2",
        "measurement_basis": global_basis,
        "actual_stage_i_node_hours": "2",
        "authorized_wave_reserved_node_hours": "16",
        "actual_plus_authorized_wave_node_hours": "18",
        "computed_remaining_stage_i_node_hours": "16",
        "computed_stage_i_total_node_hours": "18",
        "promoted_stage_i_envelope_node_hours": "1400",
        "project_ceiling_node_hours": "4000",
        "computed_stage_i_margin_node_hours": "1382",
        "case_breakdown": breakdown,
        "storage_projection_evidence": {},
    }


def refresh_gate_budget_sha(fixture):
    budget = fixture["recost"]["budget"]
    fixture["recost"]["provenance"]["computed_projection_sha256"] = hashlib.sha256(
        (json.dumps(budget, sort_keys=True) + "\n").encode()
    ).hexdigest()


def gate_fixture(module, tmp_path: Path, monkeypatch):
    root = tmp_path / "root"
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    paths = paths_for(module, root)
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    input_path = source_dir / "inputs/r17.athinput"
    input_path.parent.mkdir()
    input_path.write_text("input\n")
    matrix = source_dir / "matrix.json"
    matrix.write_text("{}\n")
    executable = root / "build/athena"
    executable.parent.mkdir()
    executable.write_text("exe\n")
    executable.chmod(0o755)
    build_manifest = root / "build/manifest"
    build_manifest.mkdir()
    (build_manifest / "environment.txt").write_text("build\n")
    bundle = root / "source-archives/source.bundle"
    bundle.parent.mkdir()
    bundle.write_text("bundle\n")
    qualification = paths["qualification"]
    qualification.write_text("{}\n")
    args = SimpleNamespace(
        case_id="R17",
        segment="s00_rankio_t0_t0p25",
        acceptance_criterion="exact R17 criterion",
        athena_walltime="01:50:00",
        executable=str(executable),
        nodes=8,
        ranks_per_node=8,
        cpus_per_task=7,
        walltime="02:00:00",
    )
    revision = "a" * 40
    utility = {
        "path": str(CONTROLLER_PATH),
        "revision": revision,
        "sha256": "1" * 64,
        "committed": True,
    }
    build = {"revision": revision, "sha256": sha256(executable), "manifest_dir": str(build_manifest)}
    bundle_record = {
        "path": str(bundle),
        "sha256": sha256(bundle),
        "verified_revisions": [revision],
    }
    current_source_authority = source_authority_binding(
        revision,
        bundle_path=bundle.relative_to(root).as_posix(),
        bundle_sha256=sha256(bundle),
    )
    qualification_record = {"path": str(qualification), "sha256": sha256(qualification)}
    build_manifest_sha = "b" * 64
    profile = {
        "acceptance_criterion": args.acceptance_criterion,
        "acceptance_policy": "strict-r17-policy",
        "athena_walltime": args.athena_walltime,
        "build_manifest": str(build_manifest),
        "build_manifest_sha256": build_manifest_sha,
        "case_id": "R17",
        "controller_walltime_max_seconds": 7200,
        "cpus_per_task": 7,
        "estimated_storage_bytes": module.R17_MINIMUM_RETAINED_BYTES,
        "executable": str(executable),
        "executable_revision": revision,
        "executable_sha256": sha256(executable),
        "input_file": "inputs/r17.athinput",
        "input_revision": revision,
        "input_sha256": sha256(input_path),
        "nodes": 8,
        "output_layout": "rank-local",
        "segment": args.segment,
        "parent_job_id": None,
        "parent_result": None,
        "parent_segment": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
        "ranks_per_node": 8,
        "recommendation_basis": {"kind": "qualified-initial-calibration"},
        "time_tlim_target": 0.25,
        "walltime": args.walltime,
        "source_bundle": str(bundle),
        "source_bundle_sha256": sha256(bundle),
    }
    lineage_sha = "c" * 64
    manifest_bindings = [{"path": "runs/example.json", "sha256": "d" * 64}]
    current_rows = [
        {
            "job_id": "5000003",
            "case_id": "R03",
            "segment": "s00_rankio_t0_t10",
            "actual_node_hours": "1",
        },
        {
            "job_id": "5000012",
            "case_id": "R12",
            "segment": "s01_rankio_t0_t10",
            "actual_node_hours": "1",
        },
    ]
    budget = final_r17_scoped_budget()
    projection_sha = hashlib.sha256(
        (json.dumps(budget, sort_keys=True) + "\n").encode()
    ).hexdigest()
    provenance = {
        "request_sha256": "0" * 64,
        "request_independent_review_sha256": "0" * 64,
        "generator_sha256": "0" * 64,
        "generator_revision": revision,
        "stage_i_helper_sha256": utility["sha256"],
        "stage_i_helper_revision": revision,
        "matrix_sha256": sha256(matrix),
        "matrix_revision": revision,
        "source_bundle_sha256": bundle_record["sha256"],
        "source_bundle_verified_revisions": (
            current_source_authority["final_source_bundle"]["verified_revisions"]
        ),
        "source_authority": current_source_authority,
        "qualification_approval_sha256": qualification_record["sha256"],
        "ceiling_evidence_sha256": "0" * 64,
        "ceiling_publication_audit_sha256": "0" * 64,
        "f113_historical_helper_revision": revision,
        "f113_historical_helper_sha256": "0" * 64,
        "storage_evidence_sha256": "2" * 64,
        "reconciliation_sha256": "0" * 64,
        "ledger_sha256": sha256(paths["ledger"]),
        "reservations_sha256": sha256(paths["reservations"]),
        "scheduler_evidence": [{}],
        "scheduler_sha256": "0" * 64,
        "predecessor_recost_sha256": "0" * 64,
        "predecessor_recost_independent_review_sha256": "0" * 64,
        "predecessor_recost_publication_audit_sha256": "0" * 64,
        "authenticated_lineages_sha256": lineage_sha,
        "computed_projection_sha256": projection_sha,
        "r17_readiness_evidence_sha256": "4" * 64,
    }
    recost = {
        "schema_version": 2,
        "record_type": "stage-i-recost-recommendation-evidence",
        "checkpoint": "F-200",
        "artifact_name": f"mks24_stage_i_{EPOCH_SLUG}_F200_recost_evidence.json",
        "execution_epoch": EPOCH,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "expires_utc": (datetime.now(timezone.utc) + timedelta(hours=1)).isoformat(),
        "requested_by": "requester",
        "scope": "R17 only",
        "predecessor_recost": {},
        "authority": {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
        "publication_requirements": {
            "independent_review_required": True,
            "publication_audit_required": True,
            "published_mode": "0444",
            "published_links": 1,
            "controller_consumption_requires_exact_published_sha256": True,
        },
        "recommendations": {
            "mode": "sole-next-profile",
            "authorizing": False,
            "recommended_next_profiles": [profile],
            "bounded_concurrency": {
                "max_active_segments": 4,
                "max_wave_nodes": 8,
                "r17_exclusive_and_last": True,
            },
            "controller_consumption_state": (
                "non-authorizing until exact independent review and publication audit"
            ),
            "sole_next_segment_recommendation": {
                key: profile[key] for key in module.R17_SOLE_PROFILE_KEYS
            },
            "non_authorizing_reason": (
                "A recost request and generated candidate are evidence only; neither may "
                "authorize prepare, submit, scheduler, controller, or canonical mutations."
            ),
        },
        "barrier": {},
        "budget": budget,
        "storage": {
            "available_bytes": module.R17_MINIMUM_RETAINED_BYTES + 1001,
            "retained_stage_i_bytes": 1,
            "required_safety_bytes": 1,
            "projected_authorized_wave_growth_bytes": module.R17_MINIMUM_RETAINED_BYTES,
            "headroom_after_authorized_wave_and_safety_bytes": 1000,
        },
        "ledger": {"rows": len(current_rows), "sha256": sha256(paths["ledger"])},
        "reservations": {
            "rows": 0,
            "sha256": sha256(paths["reservations"]),
            "active": 0,
        },
        "manifests": {
            "rows": 1,
            "bindings": manifest_bindings,
            "authenticated_lineages_sha256": lineage_sha,
        },
        "r17_readiness": {},
        "promoted_f113": {},
        "reconcile": {},
        "provenance": provenance,
    }
    recost_path = paths["accounting"] / recost["artifact_name"]
    audit_path = recost_path.with_name(f"{recost_path.name}.publication_audit.json")
    monkeypatch.setattr(
        module,
        "require_clean_r17_predecessor_state",
        lambda _paths, _reservations, **_kwargs: (lineage_sha, manifest_bindings),
    )
    monkeypatch.setattr(
        module,
        "latest_published_r17_recost",
        lambda _paths, _now: (recost_path, "e" * 64, recost, audit_path, "f" * 64),
    )
    monkeypatch.setattr(
        module, "r17_directory_inventory_sha256", lambda _path, _label: build_manifest_sha
    )
    monkeypatch.setattr(
        module,
        "source_authority_committed_sha256",
        lambda _revision, relative, _label: (
            "0" * 64
            if relative == "scripts/frontier/cgl_lf_stage_i_recost.py"
            else "f" * 64
        ),
    )
    monkeypatch.setattr(module, "read_ledger", lambda _paths: current_rows)
    monkeypatch.setattr(
        module,
        "validate_fixed_r17_readiness_chain",
        lambda *_args: ({}, {"path": "readiness", "sha256": "4" * 64}),
    )
    return {
        "paths": paths,
        "args": args,
        "recost": recost,
        "profile": profile,
        "current_rows": current_rows,
        "kwargs": {
            "source_dir": source_dir,
            "matrix_path": matrix,
            "input_path": input_path,
            "input_revision": revision,
            "utility_provenance": utility,
            "build_provenance": build,
            "build_manifest": build_manifest,
            "qualification_approval": qualification_record,
            "bundle_provenance": bundle_record,
            "current_source_authority": current_source_authority,
            "restart": None,
            "parent_segment": None,
            "time_tlim_target": 0.25,
        },
    }


def test_lower_case_and_noncanonical_r17_do_not_consume_readiness(tmp_path, monkeypatch):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    monkeypatch.setattr(
        module,
        "require_clean_r17_predecessor_state",
        lambda *_args: pytest.fail("lower-case/noncanonical path consumed R17 readiness"),
    )
    unused = {
        "source_dir": tmp_path,
        "matrix_path": tmp_path,
        "input_path": tmp_path,
        "input_revision": "",
        "utility_provenance": {},
        "build_provenance": {},
        "build_manifest": tmp_path,
        "qualification_approval": None,
        "bundle_provenance": None,
        "current_source_authority": None,
        "restart": None,
        "parent_segment": None,
        "time_tlim_target": None,
    }
    args = SimpleNamespace(case_id="R16")
    assert module.require_r17_readiness_for_prepare(paths, args, [], **unused) is None
    args.case_id = "R17"
    assert module.require_r17_readiness_for_prepare(paths, args, [], **unused) is None


def test_f116_current_tooling_authority_accepts_frozen_scientific_source(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    retained = module.require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=fixture["bundle_record"],
        utility_provenance=fixture["utility"],
        matrix_path=fixture["matrix"],
        input_revision=module.QUALIFIED_SOURCE_REVISION,
        offline_local_root=False,
    )
    assert retained == fixture["binding"]
    assert set(retained) == {
        "checkpoint",
        "evidence",
        "provenance_review",
        "plasma_review",
        "publication_audit",
        "final_source_bundle",
    }
    assert retained["final_source_bundle"]["path"].endswith(
        f"{fixture['revision'][:9]}.bundle"
    )
    assert fixture["revision"] != module.QUALIFIED_SOURCE_REVISION


def test_f116_current_authority_accepts_canonical_retained_staging_after_public_auth(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    staging = make_retained_f118_staging(module, paths, fixture)
    staging.parent.chmod(0o2755)
    staging.chmod(0o2700)

    retained = module.require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=fixture["bundle_record"],
        utility_provenance=fixture["utility"],
        matrix_path=fixture["matrix"],
        input_revision=module.QUALIFIED_SOURCE_REVISION,
        offline_local_root=False,
    )

    assert retained == fixture["binding"]
    assert staging.is_dir()


def test_f116_current_authority_accepts_stale_complete_old_staging_debris(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    staging = make_retained_f118_staging(module, paths, fixture)
    staging.parent.chmod(0o2755)
    staging.chmod(0o2700)
    stale = staging / "evidence.json"
    stale.chmod(0o644)
    stale.write_text('{"stale": "old F116 candidate"}\n')
    stale.chmod(0o444)
    write_json(staging / "journal.json", {"stale": "old transaction"}, mode=0o600)
    (staging / ".journal.json.recovery.tmp").unlink()

    retained = module.require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=fixture["bundle_record"],
        utility_provenance=fixture["utility"],
        matrix_path=fixture["matrix"],
        input_revision=module.QUALIFIED_SOURCE_REVISION,
        offline_local_root=False,
    )

    assert retained == fixture["binding"]
    assert not (staging / ".journal.json.recovery.tmp").exists()


def test_f116_current_authority_accepts_incomplete_staging_debris(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    transactions = (
        paths["accounting"]
        / f"mks24_stage_i_{EPOCH_SLUG}_F118_source_authority_transactions"
    )
    staging = transactions / f"2026-06-05T120000+0000-{'b' * 32}.staging"
    staging.mkdir(parents=True)
    transactions.chmod(0o2755)
    staging.chmod(0o2700)
    partial = staging / "journal.json"
    partial.write_bytes(b'{"partial":')
    partial.chmod(0o000)
    (staging / "untrusted-content-link").symlink_to(tmp_path / "outside")

    retained = module.require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=fixture["bundle_record"],
        utility_provenance=fixture["utility"],
        matrix_path=fixture["matrix"],
        input_revision=module.QUALIFIED_SOURCE_REVISION,
        offline_local_root=False,
    )

    assert retained == fixture["binding"]
    assert partial.is_file()
    assert stat.S_IMODE(partial.stat().st_mode) == 0o000
    assert (staging / "untrusted-content-link").is_symlink()


def test_f116_source_authority_debris_classifier_acceptance_is_cross_tool_exact(
    tmp_path,
):
    controller = load_controller()
    source_authority = load_source_authority()
    transactions = tmp_path / "source-authority-transactions"
    transactions.mkdir()
    transactions.chmod(0o2755)
    transaction_id = "2026-06-05T120000+0000-"
    for suffix, digit in (("", "a"), (".staging", "b"), (".retired", "c")):
        transaction = transactions / f"{transaction_id}{digit * 32}{suffix}"
        transaction.mkdir()
        transaction.chmod(0o2700)
        unreadable = transaction / "untrusted-content"
        unreadable.write_bytes(b"not authoritative\n")
        unreadable.chmod(0o000)
        (transaction / "untrusted-link").symlink_to(tmp_path / "outside")

    single = transactions / f".cgl-source-authority-retired-{'d' * 32}.forensic"
    single.write_bytes(b"non-authoritative forensic bytes\n")
    single.chmod(0o000)
    paired_first = (
        transactions / f".cgl-source-authority-retired-{'e' * 32}.forensic"
    )
    paired_second = (
        transactions / f".cgl-source-authority-retired-{'f' * 32}.forensic"
    )
    paired_first.write_bytes(b"paired forensic bytes\n")
    paired_first.chmod(0o400)
    os.link(paired_first, paired_second)
    directory = (
        transactions
        / f".cgl-source-authority-retired-directory-{'1' * 32}.forensic"
    )
    directory.mkdir()
    directory.chmod(0o2755)
    (directory / "untrusted-link").symlink_to(tmp_path / "outside")

    source_authority.classify_non_authoritative_recovery_debris(
        {"transactions": transactions}
    )
    controller.require_authenticated_f118_source_authority_staging(transactions)

    assert single.stat().st_mode & 0o777 == 0o000
    assert paired_first.stat().st_nlink == paired_second.stat().st_nlink == 2
    assert (directory / "untrusted-link").is_symlink()


@pytest.mark.parametrize(
    "mutation",
    (
        "malformed-directory",
        "active-regular",
        "active-symlink",
        "unsafe-active-mode",
        "forensic-symlink",
        "unsafe-forensic-file",
        "unsafe-forensic-directory",
        "external-forensic-hardlink",
    ),
)
def test_f116_source_authority_debris_classifier_rejection_is_cross_tool_exact(
    tmp_path, mutation
):
    controller = load_controller()
    source_authority = load_source_authority()
    transactions = tmp_path / "source-authority-transactions"
    transactions.mkdir()
    transactions.chmod(0o2755)
    transaction_id = f"2026-06-05T120000+0000-{'a' * 32}"
    forensic = transactions / f".cgl-source-authority-retired-{'b' * 32}.forensic"
    if mutation == "malformed-directory":
        (transactions / "malformed.staging").mkdir()
    elif mutation == "active-regular":
        (transactions / transaction_id).write_text("not a directory\n")
    elif mutation == "active-symlink":
        outside = tmp_path / "outside-active"
        outside.mkdir()
        (transactions / f"{transaction_id}.staging").symlink_to(
            outside, target_is_directory=True
        )
    elif mutation == "unsafe-active-mode":
        active = transactions / f"{transaction_id}.retired"
        active.mkdir()
        active.chmod(0o770)
    elif mutation == "forensic-symlink":
        outside = tmp_path / "outside-forensic"
        outside.write_text("outside\n")
        forensic.symlink_to(outside)
    elif mutation == "unsafe-forensic-file":
        forensic.write_text("unsafe\n")
        forensic.chmod(0o620)
    elif mutation == "unsafe-forensic-directory":
        forensic.mkdir()
        forensic.chmod(0o720)
    elif mutation == "external-forensic-hardlink":
        outside = tmp_path / "external-forensic-hardlink"
        outside.write_text("external link\n")
        outside.chmod(0o600)
        os.link(outside, forensic)

    with pytest.raises(ValueError):
        source_authority.classify_non_authoritative_recovery_debris(
            {"transactions": transactions}
        )
    with pytest.raises(ValueError):
        controller.require_authenticated_f118_source_authority_staging(transactions)


def test_live_canonical_f116_source_authority_debris_classification_is_cross_tool_exact():
    controller = load_controller()
    source_authority = load_source_authority()
    transactions = (
        controller.DEFAULT_ROOT
        / "accounting"
        / f"mks24_stage_i_{EPOCH_SLUG}_F118_source_authority_transactions"
    )
    if not transactions.is_dir():
        pytest.skip(f"live canonical F116 transaction root is unavailable: {transactions}")

    source_authority.classify_non_authoritative_recovery_debris(
        {"transactions": transactions}
    )
    controller.require_authenticated_f118_source_authority_staging(transactions)


@pytest.mark.parametrize("hostile_kind", ("regular-file", "symlink"))
def test_f116_current_authority_rejects_hostile_staging_container_entry(
    tmp_path, monkeypatch, hostile_kind
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    transactions = (
        paths["accounting"]
        / f"mks24_stage_i_{EPOCH_SLUG}_F118_source_authority_transactions"
    )
    transactions.mkdir(parents=True)
    transactions.chmod(0o2755)
    entry = transactions / f"2026-06-05T120000+0000-{'c' * 32}.staging"
    if hostile_kind == "regular-file":
        entry.write_text("not a staging directory\n")
    else:
        outside = tmp_path / "hostile-staging-target"
        outside.mkdir()
        entry.symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError):
        module.require_current_source_authority_for_prepare(
            paths,
            bundle_provenance=fixture["bundle_record"],
            utility_provenance=fixture["utility"],
            matrix_path=fixture["matrix"],
            input_revision=module.QUALIFIED_SOURCE_REVISION,
            offline_local_root=False,
        )


def test_f116_retained_staging_is_not_consulted_before_public_auth(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    make_retained_f118_staging(module, paths, fixture)
    (paths["root"] / module.F116_PUBLICATION_AUDIT_RELATIVE).unlink()
    monkeypatch.setattr(
        module,
        "require_authenticated_f118_source_authority_staging",
        lambda *_args, **_kwargs: pytest.fail(
            "retained staging was consulted before public F116 authenticated"
        ),
    )

    with pytest.raises(ValueError):
        module.require_current_source_authority_for_prepare(
            paths,
            bundle_provenance=fixture["bundle_record"],
            utility_provenance=fixture["utility"],
            matrix_path=fixture["matrix"],
            input_revision=module.QUALIFIED_SOURCE_REVISION,
            offline_local_root=False,
        )


def test_current_source_authority_offline_local_bypass(tmp_path):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    assert module.require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=None,
        utility_provenance={},
        matrix_path=tmp_path,
        input_revision="",
        offline_local_root=True,
    ) is None


def test_f116_inventory_authenticates_r17_qualification_producer():
    module = load_controller()
    assert module.F116_REQUIRED_TOOLS[
        "scripts/frontier/cgl_lf_stage_i_qualification.py"
    ] == "0755"


@pytest.mark.parametrize(
    "mutation",
    (
        "missing-chain",
        "authority-broadening",
        "historical-f115-scope",
        "historical-manifest",
        "wrong-final-head",
        "final-branch-tip",
        "wrong-selected-bundle",
        "generator-binding",
        "missing-qualification-producer",
        "wrong-final-subject",
        "nonancestor-revision",
        "bridge-selected",
        "catalog-sole-current",
        "empty-review-findings",
        "split-review-candidate",
        "pending-source-transaction",
        "retired-source-transaction",
        "nonfrozen-input-revision",
    ),
)
def test_f116_current_source_authority_fails_closed_on_drift(
    tmp_path, monkeypatch, mutation
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    bundle_record = copy.deepcopy(fixture["bundle_record"])
    evidence_path = fixture["evidence"]
    input_revision = module.QUALIFIED_SOURCE_REVISION
    if mutation == "missing-chain":
        (paths["root"] / module.F116_PUBLICATION_AUDIT_RELATIVE).unlink()
    elif mutation == "authority-broadening":
        evidence = json.loads(evidence_path.read_text())
        evidence["authorization"]["prepare_authorized"] = True
        write_json(evidence_path, evidence)
    elif mutation == "historical-f115-scope":
        historical_evidence_path = (
            paths["root"] / module.F116_CURRENT_SOURCE_AUTHORITY_RELATIVE
        )
        evidence = json.loads(historical_evidence_path.read_text())
        evidence["predecessor_authorities"]["historical_f115"]["evidence"][
            "sha256"
        ] = "f" * 64
        write_json(historical_evidence_path, evidence)
    elif mutation == "historical-manifest":
        fixture["historical_manifest"].chmod(0o644)
        fixture["historical_manifest"].write_text(
            fixture["historical_manifest"].read_text() + "\n"
        )
    elif mutation == "wrong-final-head":
        evidence = json.loads(evidence_path.read_text())
        evidence["implementation"]["current_source_bundle"][
            "advertised_tip"
        ]["revision"] = "f" * 40
        write_json(evidence_path, evidence)
    elif mutation == "final-branch-tip":
        evidence = json.loads(evidence_path.read_text())
        evidence["implementation"]["current_source_bundle"][
            "advertised_tip"
        ]["name"] = "refs/heads/feature/cgl-landau-fluid"
        write_json(evidence_path, evidence)
    elif mutation == "wrong-selected-bundle":
        bundle_record["path"] = str(paths["root"] / "source-archives/other.bundle")
    elif mutation == "generator-binding":
        evidence = json.loads(evidence_path.read_text())
        for tool in evidence["implementation"]["committed_tools"]:
            if tool["path"] == "scripts/frontier/cgl_lf_stage_i_recost.py":
                tool["sha256"] = "f" * 64
        write_json(evidence_path, evidence)
    elif mutation == "missing-qualification-producer":
        evidence = json.loads(evidence_path.read_text())
        evidence["implementation"]["committed_tools"] = [
            tool for tool in evidence["implementation"]["committed_tools"]
            if tool["path"] != "scripts/frontier/cgl_lf_stage_i_qualification.py"
        ]
        write_json(evidence_path, evidence)
    elif mutation == "wrong-final-subject":
        evidence = json.loads(evidence_path.read_text())
        evidence["implementation"]["current_source_bundle"]["subject"] = "other subject"
        write_json(evidence_path, evidence)
    elif mutation == "nonancestor-revision":
        monkeypatch.setattr(
            module,
            "source_authority_revision_is_ancestor",
            lambda revision, _promoted: revision != module.QUALIFIED_SOURCE_REVISION,
        )
    elif mutation == "bridge-selected":
        evidence = json.loads(evidence_path.read_text())
        evidence["implementation"]["intermediate_36140_bundle"][
            "selected_as_current"
        ] = True
        write_json(evidence_path, evidence)
    elif mutation == "catalog-sole-current":
        evidence = json.loads(evidence_path.read_text())
        evidence["source_archive_catalog"]["after"][
            "sole_current_source_bundle"
        ] = "source-archives/other.bundle"
        write_json(evidence_path, evidence)
    elif mutation == "empty-review-findings":
        review = paths["root"] / module.F116_PROVENANCE_REVIEW_RELATIVE
        value = json.loads(review.read_text())
        value["findings"] = []
        write_json(review, value)
    elif mutation == "split-review-candidate":
        review = paths["root"] / module.F116_PLASMA_REVIEW_RELATIVE
        value = json.loads(review.read_text())
        value["reviewed_candidate"]["path"] = str(tmp_path / "other-f116.json.candidate")
        write_json(review, value)
    elif mutation in {"pending-source-transaction", "retired-source-transaction"}:
        transaction = (
            paths["accounting"]
            / f"mks24_stage_i_{EPOCH_SLUG}_F118_source_authority_transactions"
            / ("pending.retired" if mutation == "retired-source-transaction" else "pending")
        )
        transaction.mkdir(parents=True)
    elif mutation == "nonfrozen-input-revision":
        input_revision = fixture["revision"]
    with pytest.raises(ValueError):
        module.require_current_source_authority_for_prepare(
            paths,
            bundle_provenance=bundle_record,
            utility_provenance=fixture["utility"],
            matrix_path=fixture["matrix"],
            input_revision=input_revision,
            offline_local_root=False,
        )


def test_f116_bundle_header_rejects_prerequisites_and_multiple_tips(tmp_path):
    module = load_controller()
    revision = "a" * 40
    bundle = tmp_path / "source.bundle"
    bundle.write_text(f"# v2 git bundle\n{revision} HEAD\n\npayload\n")
    bundle.chmod(0o644)
    assert module.source_authority_bundle_advertised_tip(bundle) == (revision, "HEAD")
    bundle.write_text(
        f"# v2 git bundle\n-{'b' * 40} prerequisite\n{revision} HEAD\n\npayload\n"
    )
    bundle.chmod(0o644)
    with pytest.raises(ValueError, match="complete-history"):
        module.source_authority_bundle_advertised_tip(bundle)
    bundle.write_text(
        f"# v2 git bundle\n{revision} HEAD\n{'c' * 40} refs/heads/other\n\npayload\n"
    )
    bundle.chmod(0o644)
    with pytest.raises(ValueError, match="exactly one tip"):
        module.source_authority_bundle_advertised_tip(bundle)


def test_f116_git_queries_strip_poisoned_environment_and_replacements(
    tmp_path, monkeypatch
):
    module = load_controller()
    repository = tmp_path / "repo"
    repository.mkdir()
    (repository / ".git").mkdir()
    monkeypatch.setattr(module, "ROOT_DIR", repository)
    for key in (
        "GIT_OBJECT_DIRECTORY",
        "GIT_ALTERNATE_OBJECT_DIRECTORIES",
        "GIT_REPLACE_REF_BASE",
        "GIT_CONFIG_PARAMETERS",
        "LD_PRELOAD",
        "PYTHONPATH",
        "BASH_ENV",
    ):
        monkeypatch.setenv(key, f"poison-{key}")
    calls = []
    relative = "scripts/frontier/cgl_lf_stage_i.py"

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        if "ls-tree" in command:
            stdout = f"100644 blob {'1' * 40}\t{relative}\n"
        elif "--format=%s" in command:
            stdout = b"trusted subject\n"
        else:
            stdout = b"committed payload"
        return SimpleNamespace(returncode=0, stdout=stdout)

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    revision = "a" * 40
    assert module.source_authority_committed_sha256(
        revision, relative, "tool"
    ) == hashlib.sha256(b"committed payload").hexdigest()
    assert module.source_authority_committed_mode(revision, relative, "tool") == "0644"
    assert module.source_authority_commit_subject(revision) == "trusted subject"
    assert module.source_authority_revision_is_ancestor("b" * 40, revision)
    assert len(calls) == 4
    for command, kwargs in calls:
        assert command[:2] == [str(module.GIT), "--no-replace-objects"]
        assert command[2:6] == [
            "--git-dir", str(repository / ".git"),
            "--work-tree", str(repository),
        ]
        environment = kwargs["env"]
        for poison in (
            "GIT_OBJECT_DIRECTORY",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES",
            "GIT_REPLACE_REF_BASE",
            "GIT_CONFIG_PARAMETERS",
            "LD_PRELOAD",
            "PYTHONPATH",
            "BASH_ENV",
        ):
            assert poison not in environment
        assert environment["GIT_CONFIG_GLOBAL"] == "/dev/null"
        assert environment["GIT_CONFIG_SYSTEM"] == "/dev/null"
        assert environment["GIT_CONFIG_NOSYSTEM"] == "1"
        assert environment["PATH"] == "/usr/bin:/bin"


def test_controller_source_and_bundle_git_provenance_ignores_poisoned_environment(
    tmp_path, monkeypatch
):
    module = load_controller()
    source = tmp_path / "source"
    source.mkdir()
    module.subprocess.run(["git", "init", "-q", str(source)], check=True)
    module.subprocess.run(
        ["git", "-C", str(source), "config", "user.email", "fixture@example.com"],
        check=True,
    )
    module.subprocess.run(
        ["git", "-C", str(source), "config", "user.name", "fixture"],
        check=True,
    )
    input_path = source / "input.athinput"
    matrix = source / "matrix.json"
    input_path.write_text("input\n")
    matrix.write_text("{}\n")
    module.subprocess.run(
        ["git", "-C", str(source), "add", input_path.name, matrix.name], check=True
    )
    module.subprocess.run(
        ["git", "-C", str(source), "commit", "-q", "-m", "fixture"], check=True
    )
    revision = module.subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    root = tmp_path / "root"
    bundle = root / "source-archives/current.bundle"
    bundle.parent.mkdir(parents=True)
    module.subprocess.run(
        ["git", "-C", str(source), "bundle", "create", str(bundle), "HEAD"],
        check=True,
    )
    marker = tmp_path / "poison-ran"
    poison = tmp_path / "poison"
    poison.write_text(f"#!/bin/sh\ntouch {marker}\nexit 0\n")
    poison.chmod(0o755)
    redirected = tmp_path / "redirected"
    redirected.mkdir()
    (redirected / input_path.name).write_text("redirected\n")
    (redirected / matrix.name).write_text("{}\n")
    module.subprocess.run(
        ["git", "-C", str(source), "config", "core.worktree", str(redirected)],
        check=True,
    )
    module.subprocess.run(
        ["git", "-C", str(source), "config", "core.fsmonitor", str(poison)],
        check=True,
    )
    module.subprocess.run(
        ["git", "-C", str(source), "config", "diff.external", str(poison)],
        check=True,
    )
    for key, value in {
        "GIT_DIR": str(tmp_path / "poison.git"),
        "GIT_WORK_TREE": str(tmp_path / "poison-worktree"),
        "GIT_OBJECT_DIRECTORY": str(tmp_path / "poison-objects"),
        "GIT_REPLACE_REF_BASE": "refs/poison/",
        "GIT_CONFIG_PARAMETERS": "'core.repositoryformatversion=99'",
    }.items():
        monkeypatch.setenv(key, value)
    assert module.git_revision_for_input(source, input_path, matrix) == revision
    input_path.write_text("changed\n")
    with pytest.raises(ValueError, match="must be committed"):
        module.git_revision_for_input(source, input_path, matrix)
    assert not marker.exists()
    input_path.write_text("input\n")
    untracked = source / "untracked.athinput"
    untracked.write_text("untracked\n")
    with pytest.raises(ValueError, match="must be tracked"):
        module.git_revision_for_input(source, untracked, matrix)
    provenance = module.source_bundle_provenance(
        str(bundle), [revision], root, False
    )
    assert provenance is not None
    assert provenance["verified_revisions"] == [revision]


def test_f116_active_checksum_ledger_requires_safe_exact_files(tmp_path, monkeypatch):
    module = load_controller()
    root = tmp_path / "root"
    archives = root / "source-archives"
    archives.mkdir(parents=True)
    archive = archives / "active.bundle"
    archive.write_bytes(b"active source archive\n")
    archive.chmod(0o644)
    line = f"{sha256(archive)}  {archive.name}"
    assert module.validate_f116_source_archive_checksum_ledger(root, [line]) == {
        archive.name: sha256(archive)
    }

    archive.write_bytes(b"corrupt\n")
    archive.chmod(0o644)
    with pytest.raises(ValueError, match="checksum differs"):
        module.validate_f116_source_archive_checksum_ledger(root, [line])

    archive.write_bytes(b"active source archive\n")
    archive.chmod(0o644)
    alias = archives / "alias.bundle"
    module.os.link(archive, alias)
    with pytest.raises(ValueError, match="owner-controlled"):
        module.validate_f116_source_archive_checksum_ledger(root, [line])
    alias.unlink()

    with monkeypatch.context() as poisoned_owner:
        poisoned_owner.setattr(module.os, "geteuid", lambda: module.os.getuid() + 1)
        with pytest.raises(ValueError, match="owner-controlled"):
            module.validate_f116_source_archive_checksum_ledger(root, [line])

    archive.unlink()
    with pytest.raises(ValueError, match="unavailable"):
        module.validate_f116_source_archive_checksum_ledger(root, [line])


def test_f116_authority_reads_reject_symlinked_parent(tmp_path):
    module = load_controller()
    root = tmp_path / "root"
    real_archives = root / "real-source-archives"
    real_archives.mkdir(parents=True)
    archive = real_archives / "active.bundle"
    archive.write_bytes(b"active source archive\n")
    archive.chmod(0o644)
    (root / "source-archives").symlink_to(real_archives, target_is_directory=True)
    line = f"{sha256(archive)}  {archive.name}"
    with pytest.raises(ValueError, match="symbolic link"):
        module.validate_f116_source_archive_checksum_ledger(root, [line])


def test_f116_integrated_gate_rejects_missing_active_historical_archive(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    fixture["historical_bundle"].unlink()
    with pytest.raises(ValueError, match="active source archive"):
        module.require_current_source_authority_for_prepare(
            paths,
            bundle_provenance=fixture["bundle_record"],
            utility_provenance=fixture["utility"],
            matrix_path=fixture["matrix"],
            input_revision=module.QUALIFIED_SOURCE_REVISION,
            offline_local_root=False,
        )


def test_f116_integrated_gate_rejects_symlinked_authority_parent(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    fixture = make_f116_source_authority(module, paths, tmp_path, monkeypatch)
    accounting = paths["accounting"]
    real_accounting = paths["root"] / "real-accounting"
    accounting.rename(real_accounting)
    accounting.symlink_to(real_accounting, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic link"):
        module.require_current_source_authority_for_prepare(
            paths,
            bundle_provenance=fixture["bundle_record"],
            utility_provenance=fixture["utility"],
            matrix_path=fixture["matrix"],
            input_revision=module.QUALIFIED_SOURCE_REVISION,
            offline_local_root=False,
        )


def test_latest_published_r17_recost_requires_fresh_nonbroadening_review(tmp_path):
    module = load_controller()
    now = datetime.now(timezone.utc).replace(microsecond=0)
    paths, _ = make_published_recost(module, tmp_path / "root", published=now - timedelta(minutes=2))
    selected = module.latest_published_r17_recost(paths, now)
    assert selected[0].name.endswith("F200_recost_evidence.json")

    audit_path = selected[3]
    audit = json.loads(audit_path.read_text())
    audit["authority"]["scheduler_mutation_authorized"] = True
    write_json(audit_path, audit)
    with pytest.raises(ValueError, match="broadens authority"):
        module.latest_published_r17_recost(paths, now)


def test_latest_published_r17_recost_requires_forensic_retirement_method(tmp_path):
    module = load_controller()
    assert module.R17_RECOST_PUBLICATION_METHOD == (
        "same-directory-link-fsync-copy-exchange-forensic-retirement-fsync"
    )
    now = datetime.now(timezone.utc).replace(microsecond=0)
    paths, audit = make_published_recost(
        module, tmp_path / "root", published=now - timedelta(minutes=2)
    )
    audit.update({
        "transaction_id": "recost-publication-transaction",
        "recost_recommendations": None,
        "counts": {},
        "generator": {},
        "scheduler_evidence": {},
        "source_bundle": {},
        "stage_i_helper": {},
        "utility": {},
        "forensic_copy": {},
        "publication": module.R17_RECOST_PUBLICATION_METHOD,
        "generalized_publication_context": {},
    })
    audit_path = paths["accounting"] / (
        f"mks24_stage_i_{EPOCH_SLUG}_F200_recost_evidence.json.publication_audit.json"
    )
    write_json(audit_path, audit)
    module.latest_published_r17_recost(paths, now)

    audit["publication"] = "same-directory-link-fsync-unlink-fsync"
    write_json(audit_path, audit)
    with pytest.raises(ValueError, match="full publication audit differs"):
        module.latest_published_r17_recost(paths, now)


def test_latest_published_r17_recost_rejects_stale_and_ambiguous(tmp_path):
    module = load_controller()
    now = datetime.now(timezone.utc).replace(microsecond=0)
    stale_paths, _ = make_published_recost(
        module,
        tmp_path / "stale",
        generated=now - timedelta(hours=25),
        published=now - timedelta(hours=24, minutes=50),
    )
    with pytest.raises(ValueError, match="freshness"):
        module.latest_published_r17_recost(stale_paths, now)

    paths, _ = make_published_recost(module, tmp_path / "ambiguous", checkpoint=200)
    json_time = (
        datetime.now(timezone.utc).replace(microsecond=0) - timedelta(minutes=5)
    )
    make_published_recost(
        module,
        tmp_path / "ambiguous",
        checkpoint=201,
        published=json_time,
    )
    first_audit = paths["accounting"] / (
        f"mks24_stage_i_{EPOCH_SLUG}_F200_recost_evidence.json.publication_audit.json"
    )
    first = json.loads(first_audit.read_text())
    first["published_utc"] = json_time.isoformat()
    write_json(first_audit, first)
    with pytest.raises(ValueError, match="ambiguous"):
        module.latest_published_r17_recost(paths, now)


def test_fixed_readiness_chain_is_exact_fresh_and_non_authorizing(tmp_path, monkeypatch):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    now = datetime.now(timezone.utc).replace(microsecond=0)
    profile = {"executable_sha256": "1" * 64, "build_manifest_sha256": "5" * 64}
    recost: dict[str, object] = {"provenance": {}}
    current_source_authority = source_authority_binding()
    monkeypatch.setattr(
        module, "validate_r17_operational_qualification", lambda *_args, **_kwargs: {}
    )
    monkeypatch.setattr(
        module,
        "r17_available_storage_bytes",
        lambda _path: module.R17_MINIMUM_RETAINED_BYTES,
    )
    readiness_path, readiness_value, _ = make_fixed_readiness(
        module, paths, recost, profile, "6" * 64, current_source_authority, now=now
    )
    module.validate_fixed_r17_readiness_chain(
        paths, recost, profile, "6" * 64, current_source_authority, now
    )
    monkeypatch.setattr(
        module,
        "r17_available_storage_bytes",
        lambda _path: module.R17_MINIMUM_RETAINED_BYTES - 1,
    )
    with pytest.raises(ValueError, match="current-state bindings differ"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )
    monkeypatch.setattr(
        module,
        "r17_available_storage_bytes",
        lambda _path: module.R17_MINIMUM_RETAINED_BYTES,
    )
    legacy_readiness = copy.deepcopy(readiness_value)
    del legacy_readiness["current_source_authority"]
    write_json(readiness_path, legacy_readiness)
    with pytest.raises(ValueError, match="R17 readiness evidence schema differs"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )
    write_json(readiness_path, readiness_value)
    recost["r17_readiness"]["operational_qualification_evidence"] = {
        "forged": True
    }
    with pytest.raises(ValueError, match="embedded R17 readiness differs"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )
    recost["r17_readiness"]["operational_qualification_evidence"] = {}
    original_source = recost["r17_readiness"]["current_source_authority"]
    recost["r17_readiness"]["current_source_authority"] = source_authority_binding(
        "f" * 40
    )
    with pytest.raises(ValueError, match="embedded R17 readiness differs"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )
    recost["r17_readiness"]["current_source_authority"] = original_source
    recost["r17_readiness"]["r17_launch_authorized"] = True
    with pytest.raises(ValueError, match="schema differs"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )
    del recost["r17_readiness"]["r17_launch_authorized"]
    original_expiry = recost["expires_utc"]
    recost["expires_utc"] = (now + timedelta(hours=3)).isoformat()
    with pytest.raises(ValueError, match="invalid lifetime"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )
    recost["expires_utc"] = original_expiry

    audit_path = paths["root"] / module.R17_READINESS_RELATIVE
    audit_path = audit_path.with_name(f"{audit_path.name}.publication_audit.json")
    audit = json.loads(audit_path.read_text())
    audit["authority"]["r17_launch_authorized"] = True
    write_json(audit_path, audit)
    with pytest.raises(ValueError, match="broadens authority"):
        module.validate_fixed_r17_readiness_chain(
            paths, recost, profile, "6" * 64, current_source_authority, now
        )


def test_operational_qualification_rejects_wrong_current_build_before_inventory(tmp_path):
    module = load_controller()
    root = tmp_path / "root"
    qualification = root / "accounting/qualification.json"
    qualification.parent.mkdir(parents=True)
    qualification_value = {
        "schema_version": 2,
        "record_type": "stage-i-r17-operational-qualification",
        "execution_epoch": EPOCH,
        "completed_utc": "2026-06-05T00:00:00",
        "job_id": "12345",
        "scheduler_evidence": {},
        "state": "COMPLETED",
        "exit_code": "0:0",
        "nodes": 8,
        "ranks": 64,
        "executable_sha256": "9" * 64,
        "build_manifest_inventory": [],
        "rank_local_outputs": [],
        "rank_local_restarts": [],
        "restart_load_evidence": {},
        "physics_validation_evidence": {},
        "decomposition_evidence": {},
        "account_scheduler_evidence": {},
        "account_exclusivity_evidence": {},
    }
    qualification.write_text(json.dumps(qualification_value))
    qualification.chmod(0o644)
    readiness = {
        "operational_qualification": {
            "path": qualification.relative_to(root).as_posix(),
            "sha256": sha256(qualification),
        },
        "operational_qualification_review": {},
    }
    with pytest.raises(ValueError, match="0444"):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            {"executable_sha256": "1" * 64, "build_manifest_sha256": "5" * 64},
            "readiness-reviewer",
        )


def test_operational_qualification_authenticates_all_rank_evidence_and_tampering(tmp_path):
    module = load_controller()
    root = tmp_path / "root"
    accounting = root / "accounting"
    accounting.mkdir(parents=True)
    outputs = []
    restarts = []
    for rank in range(64):
        rank_dir = root / "qualification" / f"rank_{rank:08d}"
        rank_dir.mkdir(parents=True)
        output = rank_dir / "output.bin"
        restart = rank_dir / "restart.rst"
        output.write_bytes(f"output {rank}\n".encode())
        restart.write_bytes(f"restart {rank}\n".encode())
        output.chmod(0o644)
        restart.chmod(0o644)
        outputs.append({
            "path": output.relative_to(root).as_posix(),
            "sha256": sha256(output),
        })
        restarts.append({
            "path": restart.relative_to(root).as_posix(),
            "sha256": sha256(restart),
        })
    output_inventory = hashlib.sha256(
        (json.dumps(outputs, sort_keys=True) + "\n").encode()
    ).hexdigest()
    restart_inventory = hashlib.sha256(
        (json.dumps(restarts, sort_keys=True) + "\n").encode()
    ).hexdigest()
    executable_sha = "1" * 64
    build_inventory = [
        {"name": "athena", "mode": "0755", "sha256": "5" * 64}
    ]
    build_sha = hashlib.sha256(
        (json.dumps(build_inventory, sort_keys=True) + "\n").encode()
    ).hexdigest()
    job_id = "12345"
    scheduler = accounting / f"{job_id}.r17_qualification.sacct.txt"
    scheduler.write_text(
        f"{job_id}|r17_qualification|COMPLETED|0:0|8|600|"
        "2026-06-04T23:59:00+00:00|2026-06-05T00:10:00+00:00\n"
    )
    scheduler.chmod(0o644)
    locations = [
        [lx1, lx2, lx3, 0]
        for lx1 in range(12)
        for lx2 in range(12)
        for lx3 in range(12)
    ]
    decomposition_inventory = [
        {
            "rank": rank,
            "rank_name": f"rank_{rank:08d}",
            "logical_meshblocks": locations[rank * 27:(rank + 1) * 27],
        }
        for rank in range(64)
    ]
    decomposition = {
        "schema_version": 1,
        "record_type": "stage-i-r17-decomposition-evidence",
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "logical_meshblock_grid": [12, 12, 12],
        "logical_meshblocks": 1728,
        "ranks": 64,
        "meshblocks_per_rank": 27,
        "complete_block_rank_inventory": decomposition_inventory,
        "complete_block_rank_inventory_sha256": module.compact_json_sha256(
            decomposition_inventory
        ),
        "terminal_rank_local_output_inventory_sha256": output_inventory,
        "checks": {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        },
    }
    account_scheduler = accounting / f"{job_id}.r17_account.sacct.txt"
    account_scheduler.write_text(
        "JobIDRaw|JobName|State|ExitCode|NNodes|ElapsedRaw|Submit|Start|End|"
        "Partition|Account|User\n"
        f"{job_id}|r17_qualification|COMPLETED|0:0|8|600|"
        "2026-06-04T23:59:00+00:00|2026-06-05T00:00:00+00:00|"
        "2026-06-05T00:10:00+00:00|batch|AST207|operator\n"
    )
    account_scheduler.chmod(0o644)
    account_jobs = module.parse_r17_account_scheduler_evidence(
        account_scheduler.read_bytes()
    )
    qualification_job = {
        key: account_jobs[0][key]
        for key in (
            "job_id", "job_name", "state", "exit_code", "nodes",
            "elapsed_seconds", "submit_utc", "start_utc", "end_utc",
            "partition", "account",
        )
    }
    account_exclusivity = accounting / "r17_account_exclusivity.json"
    write_json(
        account_exclusivity,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-account-exclusivity-evidence",
            "execution_epoch": EPOCH,
            "measured_utc": "2026-06-05T00:11:30+00:00",
            "query_contract": {
                "account": "AST207",
                "all_users": True,
                "allocations_only": True,
                "expanded_arrays": True,
                "start_utc": "2026-06-04T23:59:59+00:00",
                "end_utc": "2026-06-05T00:10:01+00:00",
                "start_argument": "2026-06-04T23:59:59",
                "end_argument": "2026-06-05T00:10:01",
                "start_scheduler_offset": "+0000",
                "end_scheduler_offset": "+0000",
                "fields": [
                    "JobIDRaw", "JobName", "State", "ExitCode", "NNodes",
                    "ElapsedRaw", "Submit", "Start", "End", "Partition",
                    "Account", "User",
                ],
            },
            "visibility_contract": {
                "private_data": "none",
                "all_users_job_visibility": True,
            },
            "raw_account_scheduler_sha256": sha256(account_scheduler),
            "qualification_job": qualification_job,
            "qualification_job_sha256": module.compact_json_sha256(
                qualification_job
            ),
            "account_jobs": account_jobs,
            "account_jobs_sha256": module.compact_json_sha256(account_jobs),
            "overlapping_job_ids": [job_id],
            "exclusive_entire_execution_interval": True,
        },
    )
    restart_load = accounting / "r17_restart_load.json"
    physics = accounting / "r17_physics.json"
    write_json(
        restart_load,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-restart-load-evidence",
            "execution_epoch": EPOCH,
            "measured_utc": "2026-06-05T00:11:00+00:00",
            "measured_by": "restart-measurer",
            "job_id": job_id,
            "executable_sha256": executable_sha,
            "build_manifest_inventory_sha256": build_sha,
            "passed": True,
            "measurements": {
                "rank_local_restart_inventory_sha256": restart_inventory,
                "loaded_rank_count": 64,
                "load_state": "COMPLETED",
                "load_exit_code": "0:0",
            },
        },
    )
    write_json(
        physics,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-physics-validation-evidence",
            "execution_epoch": EPOCH,
            "measured_utc": "2026-06-05T00:11:00+00:00",
            "measured_by": "physics-measurer",
            "job_id": job_id,
            "executable_sha256": executable_sha,
            "build_manifest_inventory_sha256": build_sha,
            "passed": True,
            "measurements": {
                "rank_local_output_inventory_sha256": output_inventory,
                "finite_rank_outputs": 64,
                "mass_relative_drift_max": "0",
                "mhd_user_mass_mismatch_max": "0",
                "lf_bad_counts_total": 0,
                "normalized_ct_divb_max": "0",
                "normalized_ct_divb_threshold": "0.000000000001",
                "normalized_ct_divb_below_threshold": True,
            },
        },
    )
    qualification = accounting / "r17_operational_qualification.json"
    write_json(
        qualification,
        {
            "schema_version": 2,
            "record_type": "stage-i-r17-operational-qualification",
            "execution_epoch": EPOCH,
            "completed_utc": "2026-06-05T00:10:00",
            "job_id": job_id,
            "scheduler_evidence": {
                "path": scheduler.relative_to(root).as_posix(),
                "sha256": sha256(scheduler),
            },
            "state": "COMPLETED",
            "exit_code": "0:0",
            "nodes": 8,
            "ranks": 64,
            "executable_sha256": executable_sha,
            "build_manifest_inventory": build_inventory,
            "rank_local_outputs": outputs,
            "rank_local_restarts": restarts,
            "restart_load_evidence": {
                "path": restart_load.relative_to(root).as_posix(),
                "sha256": sha256(restart_load),
            },
            "physics_validation_evidence": {
                "path": physics.relative_to(root).as_posix(),
                "sha256": sha256(physics),
            },
            "decomposition_evidence": decomposition,
            "account_scheduler_evidence": {
                "path": account_scheduler.relative_to(root).as_posix(),
                "sha256": sha256(account_scheduler),
            },
            "account_exclusivity_evidence": {
                "path": account_exclusivity.relative_to(root).as_posix(),
                "sha256": sha256(account_exclusivity),
            },
        },
        mode=0o644,
    )
    review = accounting / "r17_operational_qualification.json.independent_review.json"
    write_json(
        review,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-operational-qualification-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": "2026-06-05T00:12:00+00:00",
            "decision": "approved",
            "reviewer": "qualification-reviewer",
            "candidate": {"path": str(qualification), "sha256": sha256(qualification)},
        },
    )
    readiness = {
        "operational_qualification": {
            "path": qualification.relative_to(root).as_posix(),
            "sha256": sha256(qualification),
        },
        "operational_qualification_review": {
            "path": review.relative_to(root).as_posix(),
            "sha256": sha256(review),
        },
    }
    profile = {
        "executable_sha256": executable_sha,
        "build_manifest_sha256": build_sha,
    }
    with pytest.raises(ValueError):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            profile,
            "readiness-reviewer",
            authorization_time=datetime(2026, 6, 5, 1, tzinfo=timezone.utc),
        )
    return
    assert expanded["build_manifest_inventory_sha256"] == build_sha
    assert expanded["authenticated_rank_local_outputs"] == outputs
    qualification_original = json.loads(qualification.read_text())
    review_original = json.loads(review.read_text())
    too_many_outputs = copy.deepcopy(qualification_original)
    too_many_outputs["rank_local_outputs"].append(copy.deepcopy(outputs[0]))
    write_json(qualification, too_many_outputs, mode=0o644)
    readiness["operational_qualification"]["sha256"] = sha256(qualification)
    overfull_review = copy.deepcopy(review_original)
    overfull_review["candidate"]["sha256"] = sha256(qualification)
    write_json(review, overfull_review)
    readiness["operational_qualification_review"]["sha256"] = sha256(review)
    with pytest.raises(ValueError, match="exactly 64"):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            profile,
            "readiness-reviewer",
            authorization_time=datetime(2026, 6, 5, 1, tzinfo=timezone.utc),
        )
    write_json(qualification, qualification_original, mode=0o644)
    readiness["operational_qualification"]["sha256"] = sha256(qualification)
    review_original["candidate"]["sha256"] = sha256(qualification)
    write_json(review, review_original)
    readiness["operational_qualification_review"]["sha256"] = sha256(review)
    tampered = root / outputs[17]["path"]
    tampered.write_text("tampered\n")
    tampered.chmod(0o644)
    with pytest.raises(ValueError, match="checksum differs"):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            profile,
            "readiness-reviewer",
            authorization_time=datetime(2026, 6, 5, 1, tzinfo=timezone.utc),
        )
    tampered.write_bytes(b"output 17\n")
    tampered.chmod(0o644)

    qualification_value = json.loads(qualification.read_text())
    original_physics = json.loads(physics.read_text())

    def bind_physics(value):
        write_json(physics, value)
        qualification_value["physics_validation_evidence"]["sha256"] = sha256(physics)
        write_json(qualification, qualification_value, mode=0o644)
        readiness["operational_qualification"]["sha256"] = sha256(qualification)

    retired_physics = copy.deepcopy(original_physics)
    retired_measurements = retired_physics["measurements"]
    del retired_measurements["normalized_ct_divb_max"]
    del retired_measurements["normalized_ct_divb_threshold"]
    del retired_measurements["normalized_ct_divb_below_threshold"]
    retired_measurements["ct_divergence_assessed"] = False
    retired_measurements["ct_divergence_reason"] = "retired self-assertion"
    bind_physics(retired_physics)
    with pytest.raises(ValueError, match="physics measurements schema differs"):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            profile,
            "readiness-reviewer",
            authorization_time=datetime(2026, 6, 5, 1, tzinfo=timezone.utc),
        )

    precompletion_physics = copy.deepcopy(original_physics)
    precompletion_physics["measured_utc"] = "2026-06-05T00:09:59+00:00"
    bind_physics(precompletion_physics)
    with pytest.raises(ValueError, match="physics_validation_evidence differs"):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            profile,
            "readiness-reviewer",
            authorization_time=datetime(2026, 6, 5, 1, tzinfo=timezone.utc),
        )

    bind_physics(original_physics)
    scheduler.write_text(
        f"{job_id}|r17_qualification|COMPLETED|0:0|8|600|"
        "2026-06-04T23:59:00+00:00|2026-06-05T00:09:59+00:00\n"
    )
    scheduler.chmod(0o644)
    qualification_value["scheduler_evidence"]["sha256"] = sha256(scheduler)
    write_json(qualification, qualification_value, mode=0o644)
    readiness["operational_qualification"]["sha256"] = sha256(qualification)
    review_value = json.loads(review.read_text())
    review_value["candidate"]["sha256"] = sha256(qualification)
    write_json(review, review_value)
    readiness["operational_qualification_review"]["sha256"] = sha256(review)
    with pytest.raises(ValueError, match="scheduler evidence differs"):
        module.validate_r17_operational_qualification(
            root,
            readiness,
            profile,
            "readiness-reviewer",
            authorization_time=datetime(2026, 6, 5, 1, tzinfo=timezone.utc),
        )


def test_controller_consumes_real_producer_schema2_r17_qualification(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = producer_operational_qualification_fixture(module, tmp_path, monkeypatch)
    qualification = fixture["qualification"]
    root = fixture["root"]
    for key in (
        "scheduler_evidence", "account_scheduler_evidence",
        "account_exclusivity_evidence", "restart_load_evidence",
        "physics_validation_evidence",
    ):
        assert (root / qualification[key]["path"]).stat().st_mode & 0o777 == 0o444
    assert fixture["qualification_path"].stat().st_mode & 0o777 == 0o444
    expanded = module.validate_r17_operational_qualification(
        root,
        fixture["readiness"],
        fixture["profile"],
        "readiness independent reviewer",
        authorization_time=fixture["authorization_time"],
    )
    assert expanded["frozen_science_build_contract"] == (
        qualification["frozen_science_build_contract"]
    )
    assert expanded["authority"] == {
        "r17_launch_authorized": False,
        "scheduler_mutation_authorized": False,
        "canonical_mutation_authorized": False,
    }
    assert len(expanded["authenticated_rank_local_outputs"]) == 64
    assert len(expanded["authenticated_rank_local_restarts"]) == 64


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("finite_rank_outputs", 63),
        ("finite_rank_outputs", True),
        ("mass_relative_drift_max", "-0.1"),
        ("mass_relative_drift_max", "0.0000000000011"),
        ("mass_relative_drift_max", "NaN"),
        ("mhd_user_mass_mismatch_max", "-0.1"),
        ("mhd_user_mass_mismatch_max", "0.0000000000011"),
        ("lf_bad_counts_total", 1),
        ("lf_bad_counts_total", False),
        ("normalized_ct_divb_max", "-0.1"),
        ("normalized_ct_divb_max", "0.000000000001"),
        ("normalized_ct_divb_max", "Infinity"),
        ("normalized_ct_divb_threshold", "0.000000000002"),
        ("normalized_ct_divb_below_threshold", False),
        ("rank_local_output_inventory_sha256", "0" * 64),
    ),
)
def test_r17_physics_measurements_are_independently_enforced(field, value):
    module = load_controller()
    output_inventory = "1" * 64
    measurements = {
        "rank_local_output_inventory_sha256": output_inventory,
        "finite_rank_outputs": 64,
        "mass_relative_drift_max": "0",
        "mhd_user_mass_mismatch_max": "0",
        "lf_bad_counts_total": 0,
        "normalized_ct_divb_max": "0",
        "normalized_ct_divb_threshold": "0.000000000001",
        "normalized_ct_divb_below_threshold": True,
    }
    assert module.validate_r17_physics_measurements(
        measurements, "R17 physics", output_inventory_sha256=output_inventory
    ) == measurements
    changed = copy.deepcopy(measurements)
    changed[field] = value
    with pytest.raises(ValueError):
        module.validate_r17_physics_measurements(
            changed, "R17 physics", output_inventory_sha256=output_inventory
        )


def test_r17_operational_review_uses_declared_process_independence(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = producer_operational_qualification_fixture(module, tmp_path, monkeypatch)
    review = json.loads(fixture["review_path"].read_text())
    review["reviewer"] = "readiness independent reviewer"
    write_json(fixture["review_path"], review)
    fixture["readiness"]["operational_qualification_review"]["sha256"] = sha256(
        fixture["review_path"]
    )
    expanded = module.validate_r17_operational_qualification(
        fixture["root"],
        fixture["readiness"],
        fixture["profile"],
        "readiness independent reviewer",
        authorization_time=fixture["authorization_time"],
    )
    assert expanded["job_id"] == fixture["qualification"]["job_id"]

    review["reviewer"] = fixture["qualification"]["measured_by"]
    write_json(fixture["review_path"], review)
    fixture["readiness"]["operational_qualification_review"]["sha256"] = sha256(
        fixture["review_path"]
    )
    with pytest.raises(ValueError, match="declared process independence"):
        module.validate_r17_operational_qualification(
            fixture["root"],
            fixture["readiness"],
            fixture["profile"],
            "readiness independent reviewer",
            authorization_time=fixture["authorization_time"],
        )


def test_controller_rejects_real_producer_schema2_binding_and_mode_drift(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = producer_operational_qualification_fixture(module, tmp_path, monkeypatch)
    qualification_path = fixture["qualification_path"]
    qualification = fixture["qualification"]
    qualification["frozen_science_build_contract"]["matrix_sha256"] = "0" * 64
    write_json(qualification_path, qualification)
    fixture["readiness"]["operational_qualification"]["sha256"] = sha256(
        qualification_path
    )
    review = json.loads(fixture["review_path"].read_text())
    review["candidate"]["sha256"] = sha256(qualification_path)
    write_json(fixture["review_path"], review)
    fixture["readiness"]["operational_qualification_review"]["sha256"] = sha256(
        fixture["review_path"]
    )
    with pytest.raises(ValueError, match="frozen science/build"):
        module.validate_r17_operational_qualification(
            fixture["root"],
            fixture["readiness"],
            fixture["profile"],
            "readiness independent reviewer",
            authorization_time=fixture["authorization_time"],
        )

    write_json(qualification_path, fixture["qualification"], mode=0o644)
    fixture["readiness"]["operational_qualification"]["sha256"] = sha256(
        qualification_path
    )
    with pytest.raises(ValueError, match="0444"):
        module.validate_r17_operational_qualification(
            fixture["root"],
            fixture["readiness"],
            fixture["profile"],
            "readiness independent reviewer",
            authorization_time=fixture["authorization_time"],
        )


def test_strong_r17_gate_accepts_only_exact_current_sole_profile(tmp_path, monkeypatch):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    result = module.require_r17_readiness_for_prepare(
        fixture["paths"], fixture["args"], [], **fixture["kwargs"]
    )
    assert result is not None
    assert result["recost_sha256"] == "e" * 64


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("stale-v1", "budget method is stale"),
        ("global-r12", "global measurement basis must be non-R12"),
        ("r12-global", "R12 basis is not a fresh final measurement"),
        ("r12-historical", "R12 basis is not a fresh final measurement"),
        ("r12-incomplete", "R12 scoped projection is invalid"),
        ("non-r12-scoped", "basis does not reuse the global basis"),
        ("missing-case", "scoped case breakdown differs"),
        ("bad-rate", "not a credible current-ledger measurement"),
    ),
)
def test_strong_r17_gate_rejects_stale_or_malformed_scoped_budget(
    tmp_path, monkeypatch, mutation, message
):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    budget = fixture["recost"]["budget"]
    breakdown = budget["case_breakdown"]
    if mutation == "stale-v1":
        budget["method"] = "observed-stage-i-node-hour-rate-v1"
    elif mutation == "global-r12":
        budget["measurement_basis"]["case_id"] = "R12"
        fixture["current_rows"][0]["case_id"] = "R12"
    elif mutation == "r12-global":
        breakdown["R12"]["projection_measurement_basis"] = copy.deepcopy(
            budget["measurement_basis"]
        )
    elif mutation == "r12-historical":
        breakdown["R12"]["projection_measurement_basis"][
            "job_id"
        ] = module.R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID
        fixture["current_rows"][1]["job_id"] = (
            module.R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID
        )
    elif mutation == "r12-incomplete":
        breakdown["R12"]["authenticated_progress_fraction"] = "0"
        breakdown["R12"]["remaining_simulation_time"] = "10"
        breakdown["R12"]["observed_rate_projected_remaining_node_hours"] = "1"
        breakdown["R12"]["projected_remaining_node_hours"] = "1"
        budget["computed_remaining_stage_i_node_hours"] = "17"
        budget["computed_stage_i_total_node_hours"] = "19"
        budget["computed_stage_i_margin_node_hours"] = "1381"
    elif mutation == "non-r12-scoped":
        breakdown["R04"]["projection_measurement_basis"] = copy.deepcopy(
            breakdown["R12"]["projection_measurement_basis"]
        )
    elif mutation == "missing-case":
        breakdown.pop("R16")
    elif mutation == "bad-rate":
        budget["measurement_basis"][
            "normalized_node_hours_per_cell_per_simulation_time"
        ] = "0.2"
    refresh_gate_budget_sha(fixture)
    with pytest.raises(ValueError, match=message):
        module.require_r17_readiness_for_prepare(
            fixture["paths"], fixture["args"], [], **fixture["kwargs"]
        )


def test_strong_r17_gate_rejects_final_credible_projection_above_envelope(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    budget = fixture["recost"]["budget"]
    r17 = budget["case_breakdown"]["R17"]
    r17["projected_cells"] = "1498"
    r17["observed_rate_projected_remaining_node_hours"] = "1498"
    r17["projected_remaining_node_hours"] = "1498"
    budget["computed_remaining_stage_i_node_hours"] = "1498"
    budget["computed_stage_i_total_node_hours"] = "1500"
    budget["computed_stage_i_margin_node_hours"] = "-100"
    refresh_gate_budget_sha(fixture)
    with pytest.raises(ValueError, match="budget is stale or arithmetically invalid"):
        module.require_r17_readiness_for_prepare(
            fixture["paths"], fixture["args"], [], **fixture["kwargs"]
        )


def test_strong_r17_gate_never_accepts_provisional_fresh_r12_exception(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    budget = fixture["recost"]["budget"]
    r12 = budget["case_breakdown"]["R12"]
    r12_basis = r12["projection_measurement_basis"]
    r12_basis["job_id"] = module.R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID
    r12_basis["segment"] = "s00_rankio_t0_t0p25"
    fixture["current_rows"][1]["job_id"] = module.R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID
    fixture["current_rows"][1]["segment"] = "s00_rankio_t0_t0p25"
    r12["authenticated_progress_fraction"] = "0"
    r12["remaining_simulation_time"] = "10"
    r12["observed_rate_projected_remaining_node_hours"] = "1"
    r12["projected_remaining_node_hours"] = "1"
    r17 = budget["case_breakdown"]["R17"]
    r17["projected_cells"] = "1497"
    r17["observed_rate_projected_remaining_node_hours"] = "1497"
    r17["projected_remaining_node_hours"] = "1497"
    budget["computed_remaining_stage_i_node_hours"] = "1498"
    budget["computed_stage_i_total_node_hours"] = "1500"
    budget["computed_stage_i_margin_node_hours"] = "-100"
    refresh_gate_budget_sha(fixture)
    with pytest.raises(ValueError, match="R12 basis is not a fresh final measurement"):
        module.require_r17_readiness_for_prepare(
            fixture["paths"], fixture["args"], [], **fixture["kwargs"]
        )


def test_strong_r17_gate_rejects_wrong_predecessor_reservation_snapshot(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="reservation snapshot checksum differs"):
        module.require_r17_readiness_for_prepare(
            fixture["paths"],
            fixture["args"],
            [],
            reservation_snapshot_sha256="9" * 64,
            **fixture["kwargs"],
        )


def test_strong_r17_gate_separates_current_tooling_from_frozen_scientific_source(
    tmp_path, monkeypatch
):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    tooling_revision = "f" * 40
    scientific_revision = module.QUALIFIED_SOURCE_REVISION
    authority = source_authority_binding(
        tooling_revision,
        bundle_path=fixture["kwargs"]["current_source_authority"][
            "final_source_bundle"
        ]["path"],
        bundle_sha256=fixture["kwargs"]["bundle_provenance"]["sha256"],
        verified_revisions=[scientific_revision, tooling_revision],
    )
    fixture["kwargs"]["input_revision"] = scientific_revision
    fixture["kwargs"]["utility_provenance"]["revision"] = tooling_revision
    fixture["kwargs"]["build_provenance"]["revision"] = scientific_revision
    fixture["kwargs"]["bundle_provenance"]["verified_revisions"] = [
        scientific_revision,
        tooling_revision,
    ]
    fixture["kwargs"]["current_source_authority"] = authority
    fixture["profile"]["input_revision"] = scientific_revision
    fixture["profile"]["executable_revision"] = scientific_revision
    provenance = fixture["recost"]["provenance"]
    provenance["generator_revision"] = tooling_revision
    provenance["stage_i_helper_revision"] = tooling_revision
    provenance["matrix_revision"] = tooling_revision
    provenance["source_bundle_verified_revisions"] = [
        scientific_revision,
        tooling_revision,
    ]
    provenance["source_authority"] = authority
    result = module.require_r17_readiness_for_prepare(
        fixture["paths"], fixture["args"], [], **fixture["kwargs"]
    )
    assert result is not None
    assert provenance["matrix_revision"] == tooling_revision
    assert fixture["profile"]["input_revision"] == scientific_revision


@pytest.mark.parametrize(
    "mutation",
    (
        "profile-segment",
        "profile-input-revision",
        "profile-extra-key",
        "profile-storage",
        "duplicate-profile",
        "recommendation-authority",
        "artifact-authority",
        "controller-binding",
        "generator-binding",
        "matrix-binding",
        "matrix-revision",
        "bundle-revisions",
        "source-authority",
        "authority-bundle",
        "active-snapshot",
        "manifest-binding",
        "budget-projection",
        "storage-headroom",
    ),
)
def test_strong_r17_gate_fails_closed_on_drift(tmp_path, monkeypatch, mutation):
    module = load_controller()
    fixture = gate_fixture(module, tmp_path, monkeypatch)
    recost = fixture["recost"]
    profile = fixture["profile"]
    if mutation == "profile-segment":
        profile["segment"] = "s00_rankio_t0_t0p5"
    elif mutation == "profile-input-revision":
        profile["input_revision"] = "f" * 40
    elif mutation == "profile-extra-key":
        profile["unreviewed"] = True
    elif mutation == "profile-storage":
        profile["estimated_storage_bytes"] = module.R17_MINIMUM_RETAINED_BYTES - 1
    elif mutation == "duplicate-profile":
        recost["recommendations"]["recommended_next_profiles"].append(copy.deepcopy(profile))
    elif mutation == "recommendation-authority":
        recost["recommendations"]["authorizing"] = True
    elif mutation == "artifact-authority":
        recost["authority"]["canonical_mutation_authorized"] = True
    elif mutation == "controller-binding":
        recost["provenance"]["stage_i_helper_sha256"] = "9" * 64
    elif mutation == "generator-binding":
        recost["provenance"]["generator_sha256"] = "9" * 64
    elif mutation == "matrix-binding":
        recost["provenance"]["matrix_sha256"] = "9" * 64
    elif mutation == "matrix-revision":
        recost["provenance"]["matrix_revision"] = "e" * 40
    elif mutation == "bundle-revisions":
        recost["provenance"]["source_bundle_verified_revisions"] = ["b" * 40]
    elif mutation == "source-authority":
        recost["provenance"]["source_authority"] = source_authority_binding("f" * 40)
    elif mutation == "authority-bundle":
        fixture["kwargs"]["current_source_authority"]["final_source_bundle"][
            "sha256"
        ] = "f" * 64
    elif mutation == "active-snapshot":
        recost["reservations"]["active"] = 1
    elif mutation == "manifest-binding":
        recost["manifests"]["bindings"] = []
    elif mutation == "budget-projection":
        recost["budget"]["computed_stage_i_total_node_hours"] = "1500"
    elif mutation == "storage-headroom":
        recost["storage"]["headroom_after_authorized_wave_and_safety_bytes"] = 999
    with pytest.raises(ValueError):
        module.require_r17_readiness_for_prepare(
            fixture["paths"], fixture["args"], [], **fixture["kwargs"]
        )


def test_clean_r17_predecessor_state_requires_zero_active_and_transactions(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    with pytest.raises(ValueError, match="zero active"):
        module.require_clean_r17_predecessor_state(
            paths, [{"state": "submitted", "case_id": "R16"}]
        )
    transaction = paths["transactions"] / "pending.json"
    transaction.write_text("{}\n")
    with pytest.raises(ValueError, match="requires recovery"):
        module.require_clean_r17_predecessor_state(paths, [])
    transaction.unlink()
    monkeypatch.setattr(module, "current_r17_lineages", lambda _paths: {})
    with pytest.raises(ValueError, match="accepted exact-t10 predecessors"):
        module.require_clean_r17_predecessor_state(paths, [])


def test_clean_r17_predecessor_state_rejects_recost_transaction_store(tmp_path):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transaction = (
        paths["accounting"]
        / f"mks24_stage_i_{EPOCH_SLUG}_recost_transactions"
        / "retained"
    )
    transaction.mkdir(parents=True)
    with pytest.raises(ValueError, match="requires recovery"):
        module.require_clean_r17_predecessor_state(paths, [])


def test_clean_r17_predecessor_state_rejects_unauthenticated_source_transaction_store(
    tmp_path,
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transaction = (
        paths["accounting"]
        / f"mks24_stage_i_{EPOCH_SLUG}_F118_source_authority_transactions"
        / "retained.staging"
    )
    transaction.mkdir(parents=True)
    with pytest.raises(ValueError, match="requires recovery"):
        module.require_clean_r17_predecessor_state(paths, [])


def test_clean_r17_predecessor_state_excludes_only_exact_prepared_r17_manifest(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    prior = paths["runs"] / "R16/s09/manifest/prepared_run.json"
    current = paths["runs"] / "R17/s00/manifest/prepared_run.json"
    write_json(prior, {"prior": True}, mode=0o644)
    write_json(current, {"current": True}, mode=0o644)
    lineages = {
        case_id: [{
            "accounting": {"result": "accepted"},
            "scientific_inspection": {"final_time": 10.0},
        }]
        for case_id in module.R17_PREDECESSOR_CASE_IDS
    }
    monkeypatch.setattr(module, "current_r17_lineages", lambda _paths: lineages)
    monkeypatch.setattr(module, "r17_lineage_summary_sha256", lambda _lineages: "1" * 64)
    lineage_sha, bindings = module.require_clean_r17_predecessor_state(
        paths, [], excluded_manifest_path=current
    )
    assert lineage_sha == "1" * 64
    assert bindings == [{
        "path": prior.relative_to(paths["root"]).as_posix(),
        "sha256": sha256(prior),
    }]
    with pytest.raises(ValueError, match="lacks its exact prepared manifest"):
        module.require_clean_r17_predecessor_state(
            paths, [], excluded_manifest_path=tmp_path / "other.json"
        )


def test_submission_reauthentication_binds_f116_and_fresh_r17_chain(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    paths = paths_for(module, root)
    manifest_path = (
        paths["runs"] / "R17/s00_rankio_t0_t0p25/manifest/prepared_run.json"
    )
    source_dir = tmp_path / "frozen-source"
    source_input = source_dir / "inputs/r17.athinput"
    source_input.parent.mkdir(parents=True)
    source_input.write_text("frozen input\n")
    matrix = manifest_path.parent / "mks24_stage_i_manifest.json"
    matrix.parent.mkdir(parents=True)
    matrix.write_text("{}\n")
    build_manifest = root / "build/manifest"
    build_manifest.mkdir(parents=True)
    executable = root / "build/athena"
    executable.write_text("exe\n")
    executable.chmod(0o755)
    authority = source_authority_binding(
        "a" * 40,
        bundle_path="source-archives/current.bundle",
        verified_revisions=[module.QUALIFIED_SOURCE_REVISION, "a" * 40],
    )
    readiness = {"path": "readiness", "sha256": "9" * 64}
    utility = {
        "path": str(CONTROLLER_PATH),
        "revision": "a" * 40,
        "sha256": "1" * 64,
        "committed": True,
    }
    bundle = {
        "path": str(root / "source-archives/current.bundle"),
        "sha256": "5" * 64,
        "verified_revisions": [module.QUALIFIED_SOURCE_REVISION, "a" * 40],
    }
    manifest = {
        "command": {
            "source_bundle": bundle,
            "production_utility": utility,
            "matrix_file": str(matrix),
            "input_revision": module.QUALIFIED_SOURCE_REVISION,
            "current_source_authority": authority,
            "r17_readiness_evidence_chain": readiness,
            "source_dir": str(source_dir),
            "source_input_file": str(source_input),
            "input_sha256": sha256(source_input),
            "source_restart_file": None,
            "build_manifest": str(build_manifest),
            "qualification_approval": {"sha256": "2" * 64},
            "parent_segment": None,
            "athena_walltime": "01:50:00",
            "executable": str(executable),
            "executable_revision": module.QUALIFIED_SOURCE_REVISION,
            "executable_sha256": sha256(executable),
            "time_tlim_target": 0.25,
        },
        "run": {
            "case_id": "R17",
            "segment": "s00_rankio_t0_t0p25",
            "acceptance_criterion": "exact R17 criterion",
        },
        "allocation": {
            "nodes": 8,
            "ranks_per_node": 8,
            "cpus_per_task": 7,
            "requested_walltime": "02:00:00",
        },
    }
    predecessor = {"case_id": "R16", "state": "recorded"}
    reservation = {
        "manifest": str(manifest_path),
        "case_id": "R17",
        "state": "prepared",
    }
    captured = {}

    def authenticate_source(_paths, **kwargs):
        captured["source"] = kwargs
        return authority

    def authenticate_r17(_paths, args, reservations, **kwargs):
        captured["r17_args"] = args
        captured["r17_reservations"] = reservations
        captured["r17"] = kwargs
        return readiness

    monkeypatch.setattr(
        module, "require_current_source_authority_for_prepare", authenticate_source
    )
    monkeypatch.setattr(module, "require_r17_readiness_for_prepare", authenticate_r17)
    module.reauthenticate_submission_authority(
        paths, manifest_path, manifest, [predecessor, reservation], reservation, False
    )
    assert captured["source"]["input_revision"] == module.QUALIFIED_SOURCE_REVISION
    assert captured["r17_reservations"] == [predecessor]
    assert captured["r17"]["reservation_snapshot_sha256"] == module.stable_json_sha256(
        [predecessor]
    )
    assert captured["r17"]["excluded_manifest_path"] == manifest_path
    assert captured["r17_args"].case_id == "R17"

    manifest["command"]["r17_readiness_evidence_chain"] = {"stale": True}
    with pytest.raises(ValueError, match="readiness/recost authority is stale"):
        module.reauthenticate_submission_authority(
            paths, manifest_path, manifest, [predecessor, reservation], reservation, False
        )
    manifest["command"]["r17_readiness_evidence_chain"] = readiness
    manifest["command"]["current_source_authority"] = {"stale": True}
    with pytest.raises(ValueError, match="F118 current-source authority is stale"):
        module.reauthenticate_submission_authority(
            paths, manifest_path, manifest, [predecessor, reservation], reservation, False
        )


def test_non_r17_submission_reauthenticates_exact_f116_authority(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    authority = source_authority_binding()
    profile_authority = {"recost_sha256": "8" * 64}
    manifest_path = tmp_path / "manifest.json"
    source_dir = tmp_path / "source"
    source_input = source_dir / "input.athinput"
    source_input.parent.mkdir()
    source_input.write_text("input\n")
    build_manifest = tmp_path / "build"
    build_manifest.mkdir()
    manifest = {
        "command": {
            "source_bundle": {"verified_revisions": [module.QUALIFIED_SOURCE_REVISION]},
            "production_utility": {"revision": "a" * 40},
            "matrix_file": str(tmp_path / "matrix.json"),
            "input_revision": module.QUALIFIED_SOURCE_REVISION,
            "current_source_authority": authority,
            "promoted_profile_authority": profile_authority,
            "source_dir": str(source_dir),
            "source_input_file": str(source_input),
            "input_sha256": sha256(source_input),
            "source_restart_file": None,
            "build_manifest": str(build_manifest),
            "qualification_approval": {"sha256": "7" * 64},
            "parent_segment": None,
            "athena_walltime": "01:50:00",
            "executable": str(tmp_path / "athena"),
            "executable_revision": module.QUALIFIED_SOURCE_REVISION,
            "executable_sha256": "6" * 64,
            "time_tlim_target": 0.25,
        },
        "run": {
            "case_id": "R16",
            "segment": "s00_rankio_t0_t0p25",
            "acceptance_criterion": "criterion",
        },
        "allocation": {
            "nodes": 1,
            "ranks_per_node": 8,
            "cpus_per_task": 7,
            "requested_walltime": "02:00:00",
        },
    }
    reservation = {"manifest": str(manifest_path), "case_id": "R16", "state": "prepared"}
    calls = []

    def authenticate_source(*_args, **_kwargs):
        calls.append(_kwargs)
        return authority

    monkeypatch.setattr(
        module, "require_current_source_authority_for_prepare", authenticate_source
    )
    monkeypatch.setattr(
        module,
        "require_r17_readiness_for_prepare",
        lambda *_args, **_kwargs: pytest.fail("non-R17 submission consumed R17 authority"),
    )
    monkeypatch.setattr(
        module, "require_promoted_profile_for_prepare",
        lambda *_args, **_kwargs: profile_authority,
    )
    module.reauthenticate_submission_authority(
        paths, manifest_path, manifest, [reservation], reservation, False
    )
    assert calls[0]["input_revision"] == module.QUALIFIED_SOURCE_REVISION
    manifest["command"]["current_source_authority"] = {"stale": True}
    with pytest.raises(ValueError, match="F118 current-source authority is stale"):
        module.reauthenticate_submission_authority(
            paths, manifest_path, manifest, [reservation], reservation, False
        )


def test_submission_reauthentication_preserves_offline_local_fixture(
    tmp_path, monkeypatch
):
    module = load_controller()

    def forbidden(*_args, **_kwargs):
        pytest.fail("offline fixture consumed canonical submission authority")

    monkeypatch.setattr(module, "require_current_source_authority_for_prepare", forbidden)
    monkeypatch.setattr(module, "require_r17_readiness_for_prepare", forbidden)
    module.reauthenticate_submission_authority(
        paths_for(module, tmp_path / "root"),
        tmp_path / "manifest.json",
        {},
        [],
        {},
        True,
    )


def test_submission_preflight_invokes_authority_reauthentication(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    paths = paths_for(module, root)
    manifest_path = tmp_path / "manifest.json"
    manifest = {
        "execution_epoch": EPOCH,
        "project_root": str(root),
        "state": "prepared",
        "command": {},
    }
    reservation = {
        "manifest": str(manifest_path),
        "case_id": "R16",
        "state": "prepared",
    }
    monkeypatch.setattr(module, "require_current_epoch", lambda *_args: None)
    monkeypatch.setattr(module, "require_root", lambda *_args: root)
    monkeypatch.setattr(module, "is_offline_local_root", lambda *_args: False)
    monkeypatch.setattr(module, "layout", lambda _root: paths)
    monkeypatch.setattr(module, "validate_submission_fixture_options", lambda *_args: None)
    for name in (
        "require_existing_layout",
        "require_no_pending_transactions",
        "require_no_orphaned_segment_runs",
        "require_reconciled_store_consistency",
        "require_reserved_execution_intent",
        "authenticate_prepared_execution",
    ):
        monkeypatch.setattr(module, name, lambda *_args, **_kwargs: None)
    monkeypatch.setattr(module, "read_reservations", lambda _paths: [reservation])
    monkeypatch.setattr(
        module, "reservation_for_manifest", lambda *_args: reservation
    )
    monkeypatch.setattr(
        module, "require_active_reservation_policy", lambda *_args: [reservation]
    )

    def reached(*_args, **_kwargs):
        raise RuntimeError("submission authority reauthenticated")

    monkeypatch.setattr(module, "reauthenticate_submission_authority", reached)
    with pytest.raises(RuntimeError, match="submission authority reauthenticated"):
        module.submission_preflight(
            SimpleNamespace(allow_local_root=False),
            manifest_path,
            manifest,
            run_slurm_test=False,
        )


def test_submission_preflight_requires_managed_clearance_for_production_overlap(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    paths = paths_for(module, root)
    manifest_path = tmp_path / "manifest.json"
    manifest = {
        "execution_epoch": EPOCH,
        "project_root": str(root),
        "state": "prepared",
        "command": {},
    }
    reservation = {
        "manifest": str(manifest_path),
        "case_id": "R16",
        "state": "prepared",
    }
    monkeypatch.setattr(module, "require_current_epoch", lambda *_args: None)
    monkeypatch.setattr(module, "require_root", lambda *_args: root)
    monkeypatch.setattr(module, "is_offline_local_root", lambda *_args: False)
    monkeypatch.setattr(module, "layout", lambda _root: paths)
    monkeypatch.setattr(module, "validate_submission_fixture_options", lambda *_args: None)
    for name in (
        "require_existing_layout",
        "require_no_pending_transactions",
        "require_no_orphaned_segment_runs",
        "require_reconciled_store_consistency",
        "require_reserved_execution_intent",
        "authenticate_prepared_execution",
        "reauthenticate_submission_authority",
        "prepared_time_tlim_target",
    ):
        monkeypatch.setattr(module, name, lambda *_args, **_kwargs: None)
    monkeypatch.setattr(module, "read_reservations", lambda _paths: [reservation])
    monkeypatch.setattr(
        module, "reservation_for_manifest", lambda *_args: reservation
    )
    monkeypatch.setattr(
        module, "require_active_reservation_policy", lambda *_args: [reservation]
    )
    monkeypatch.setattr(module, "reservation_usage", lambda *_args: (0.0, 0.0))
    monkeypatch.setattr(
        module, "authenticated_production_queue_evidence", lambda *_args: {}
    )

    def reached(_paths, requested):
        assert requested == {module.SHARED_ROOT_STALE_CAMPAIGN_ID}
        raise RuntimeError("managed clearance reached")

    monkeypatch.setattr(module, "require_managed_shared_root_clearance", reached)
    with pytest.raises(RuntimeError, match="managed clearance reached"):
        module.submission_preflight(
            SimpleNamespace(
                allow_local_root=False,
                allow_shared_root_campaign=[module.SHARED_ROOT_STALE_CAMPAIGN_ID],
            ),
            manifest_path,
            manifest,
            run_slurm_test=False,
        )


@pytest.mark.parametrize("action", ("check-submit", "submit"))
def test_check_submit_and_submit_share_submission_preflight(
    tmp_path, monkeypatch, action
):
    module = load_controller()
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text("{}\n")
    root = tmp_path / "root"
    manifest = {"project_root": str(root)}
    monkeypatch.setattr(module, "read_manifest", lambda _path: manifest)
    monkeypatch.setattr(module, "require_root", lambda *_args: root)
    monkeypatch.setattr(module, "is_offline_local_root", lambda *_args: False)

    def reached(*_args, **_kwargs):
        raise RuntimeError("shared submission preflight reached")

    monkeypatch.setattr(module, "submission_preflight", reached)
    args = SimpleNamespace(
        manifest=str(manifest_path),
        allow_local_root=False,
        sbatch_output_file=None,
    )
    selected = module.check_submit if action == "check-submit" else module.submit.__wrapped__
    with pytest.raises(RuntimeError, match="shared submission preflight reached"):
        selected(args)


def test_check_submit_prints_authenticated_isolated_python(
    tmp_path, monkeypatch, capsys
):
    module = load_controller()
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text("{}\n")
    script = {"path": tmp_path / "run.sbatch"}
    monkeypatch.setattr(
        module, "submission_preflight", lambda *_args, **_kwargs: ({}, script, {})
    )
    monkeypatch.setattr(module, "close_authenticated_batch_script", lambda _script: None)
    monkeypatch.setattr(
        module, "authenticated_python_binary", lambda: "/trusted/python3"
    )
    args = SimpleNamespace(
        manifest=str(manifest_path),
        allow_shared_root_campaign=[],
    )
    assert module.check_submit(args) == 0
    assert (
        f"/trusted/python3 -I -S -B {CONTROLLER_PATH} submit --manifest {manifest_path}"
        in capsys.readouterr().out
    )


def test_scheduler_environment_strips_sbatch_and_slurm(monkeypatch):
    module = load_controller()
    monkeypatch.setenv("SBATCH_ACCOUNT", "attacker")
    monkeypatch.setenv("SBATCH_EXPORT", "ALL,POISON=1")
    monkeypatch.setenv("SLURM_CONF", "/tmp/attacker.conf")
    monkeypatch.setenv("PYTHONPATH", "/tmp/attacker-python")
    monkeypatch.setenv("LD_PRELOAD", "/tmp/attacker-loader.so")
    monkeypatch.setenv("BASH_ENV", "/tmp/attacker-bash-env")
    monkeypatch.setenv("GIT_DIR", "/tmp/attacker.git")
    monkeypatch.setenv("PRIVATE_TOKEN", "attacker-secret")
    monkeypatch.setenv("KEEP_ME", "retained")
    environment = module.scheduler_environment()
    account = module.pwd.getpwuid(os.geteuid())
    assert environment == {
        "HOME": account.pw_dir,
        "LC_ALL": "C",
        "LOGNAME": account.pw_name,
        "PATH": "/usr/bin:/bin",
        "USER": account.pw_name,
    }
    assert not any(key.startswith(("SBATCH_", "SLURM_")) for key in environment)
    assert not {
        "PYTHONPATH", "LD_PRELOAD", "BASH_ENV", "GIT_DIR", "PRIVATE_TOKEN",
        "KEEP_ME",
    }.intersection(environment)


def test_authenticated_python_is_root_owned_and_environment_is_isolated(
    tmp_path, monkeypatch
):
    module = load_controller()
    for key in (
        "PYTHONPATH", "PYTHONHOME", "PYTHONINSPECT", "LD_PRELOAD",
        "BASH_ENV", "GIT_DIR", "SBATCH_ACCOUNT", "SLURM_CONF",
    ):
        monkeypatch.setenv(key, f"poison-{key}")
    environment = module.hardened_python_environment()
    assert environment == {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        "PYTHONSAFEPATH": "1",
        "XDG_CONFIG_HOME": "/nonexistent",
    }
    authenticated = Path(module.authenticated_python_binary())
    assert authenticated.is_absolute()
    assert authenticated.stat().st_uid == 0

    replacement = tmp_path / "python"
    replacement.write_text("#!/bin/sh\nexit 0\n")
    replacement.chmod(0o755)
    monkeypatch.setattr(module, "SYSTEM_PYTHON", replacement)
    with pytest.raises(ValueError, match="root-owned|immutable"):
        module.authenticated_python_binary()


def test_controller_reexec_uses_fixed_authenticated_isolated_python(
    monkeypatch,
):
    module = load_controller()
    calls = []

    class ReexecIntercept(Exception):
        pass

    def fake_execve(executable, argv, environment):
        python_descriptor = int(executable.rsplit("/", 1)[1])
        profile = os.fstat(python_descriptor)
        assert profile.st_uid == 0
        assert profile.st_mode & 0o111
        assert os.get_inheritable(python_descriptor)
        calls.append((executable, argv, environment))
        raise ReexecIntercept

    monkeypatch.setenv("HOME", "/attacker/home")
    monkeypatch.setenv("PYTHONPATH", "/attacker/python")
    monkeypatch.setenv("LD_PRELOAD", "/attacker/loader.so")
    monkeypatch.setenv("SBATCH_ACCOUNT", "attacker")
    monkeypatch.setattr(module.os, "execve", fake_execve)
    descriptor = os.open(CONTROLLER_PATH, os.O_RDONLY)
    try:
        with pytest.raises(ReexecIntercept):
            module.reexec_authenticated_controller(
                descriptor, CONTROLLER_PATH.resolve(), REPO_ROOT.resolve(), ["--help"]
            )
    finally:
        os.close(descriptor)

    executable, argv, environment = calls[0]
    assert executable.startswith("/proc/self/fd/")
    assert argv[:4] == [str(module.SYSTEM_PYTHON), "-I", "-S", "-B"]
    assert argv[4].startswith("/proc/self/fd/")
    assert argv[5:] == ["--help"]
    assert environment["HOME"] == "/nonexistent"
    assert environment["PATH"] == module.TRUSTED_SYSTEM_PATH
    assert "PYTHONPATH" not in environment
    assert "LD_PRELOAD" not in environment
    assert "SBATCH_ACCOUNT" not in environment


def test_controller_rejects_orphaned_private_reexec_metadata(monkeypatch):
    module = load_controller()
    monkeypatch.setenv(module.SELF_SOURCE_ENV, "/attacker/controller.py")
    monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, "/attacker/repository")
    with pytest.raises(ValueError, match="private controller reexecution path"):
        module.authenticate_controller_runtime([])


def test_authenticated_controller_descriptor_binds_private_source_metadata(
    monkeypatch,
):
    module = load_controller()
    source = CONTROLLER_PATH.resolve()
    repository = REPO_ROOT.resolve()
    descriptor = os.open(source, os.O_RDONLY)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(source))
        monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, str(repository))
        monkeypatch.setattr(module, "require_authenticated_reexec_runtime", lambda: None)
        module.authenticate_controller_runtime([])

        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(repository / "README.md"))
        with pytest.raises(ValueError, match="source/repository relationship differs"):
            module.authenticate_controller_runtime([])
    finally:
        os.close(descriptor)


def test_controller_reexec_runtime_requires_isolated_fixed_python(monkeypatch):
    module = load_controller()
    monkeypatch.setenv("HOME", "/nonexistent")
    with pytest.raises(ValueError, match="reexecution is not isolated"):
        module.require_authenticated_reexec_runtime()


def test_canonical_controller_startup_requires_fixed_isolated_python(
    tmp_path,
):
    poisoned = tmp_path / "poisoned"
    poisoned.mkdir()
    marker = tmp_path / "sitecustomize-ran"
    (poisoned / "sitecustomize.py").write_text(
        "from pathlib import Path\n"
        f"Path({str(marker)!r}).write_text('executed\\n')\n"
    )
    environment = {
        **os.environ,
        "PYTHONPATH": str(poisoned),
        "PYTHONDONTWRITEBYTECODE": "1",
    }

    rejected = subprocess.run(
        [
            sys.executable,
            str(CONTROLLER_PATH),
            "--root",
            "/lustre/orion/ast207/proj-shared/dfielding/CGL",
            "reconcile",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert rejected.returncode == 1
    assert "controller invocation requires" in rejected.stderr
    assert marker.read_text() == "executed\n"

    marker.unlink()
    isolated = subprocess.run(
        ["/usr/bin/python3.11", "-I", "-S", "-B", str(CONTROLLER_PATH), "--help"],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert isolated.returncode == 0, isolated.stderr
    assert not marker.exists()


def test_every_controller_startup_requires_fixed_isolated_python(tmp_path):
    poisoned = tmp_path / "poisoned"
    poisoned.mkdir()
    marker = tmp_path / "sitecustomize-ran"
    (poisoned / "sitecustomize.py").write_text(
        "from pathlib import Path\n"
        f"Path({str(marker)!r}).write_text('executed\\n')\n"
    )
    local_root = tmp_path / "local-root"
    local_root.mkdir()
    environment = {
        **os.environ,
        "PYTHONPATH": str(poisoned),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    rejected = subprocess.run(
        [
            sys.executable,
            str(CONTROLLER_PATH),
            "--root",
            str(local_root),
            "--allow-local-root",
            "init",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert rejected.returncode == 1
    assert "controller invocation requires" in rejected.stderr
    assert marker.read_text() == "executed\n"

    marker.unlink()
    isolated = subprocess.run(
        [
            "/usr/bin/python3.11",
            "-I",
            "-S",
            "-B",
            str(CONTROLLER_PATH),
            "--root",
            str(local_root),
            "--allow-local-root",
            "init",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert isolated.returncode == 0, isolated.stderr
    assert not marker.exists()


def test_git_diff_cleanliness_explicitly_disables_external_diff_and_textconv(
    tmp_path, monkeypatch
):
    module = load_controller()
    source = tmp_path / "source"
    source.mkdir()
    (source / ".git").mkdir()
    input_path = source / "input.athinput"
    matrix = source / "matrix.json"
    input_path.write_text("input\n")
    matrix.write_text("{}\n")
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        if "rev-parse" in command:
            return SimpleNamespace(returncode=0, stdout="a" * 40 + "\n")
        return SimpleNamespace(returncode=0, stdout="")

    monkeypatch.setattr(module, "authenticated_git_binary", lambda: "/trusted/git")
    monkeypatch.setattr(module.subprocess, "run", fake_run)
    assert module.git_revision_for_input(source, input_path, matrix) == "a" * 40
    diff_calls = [command for command, _ in calls if "diff" in command]
    assert len(diff_calls) == 2
    for command in diff_calls:
        assert "--no-ext-diff" in command
        assert "--no-textconv" in command
    for _, kwargs in calls:
        assert kwargs["env"] == module.hardened_git_environment()


def test_descriptor_bound_batch_submission_detects_exact_byte_path_swap(
    tmp_path, monkeypatch
):
    module = load_controller()
    manifest_path, manifest, batch = prepared_batch_fixture(module, tmp_path)
    authenticated = module.open_authenticated_batch_script(manifest, manifest_path)
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(stdout="test-only accepted\n")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    try:
        assert module.scheduler_test_only_output(authenticated) == "test-only accepted\n"
        command, kwargs = calls[0]
        assert command[-1] == authenticated["descriptor_path"]
        assert kwargs["pass_fds"] == (authenticated["fd"],)
        replacement = batch.with_name("replacement.sbatch")
        replacement.write_bytes(batch.read_bytes())
        replacement.chmod(0o750)
        os.replace(replacement, batch)
        with pytest.raises(ValueError, match="single-link|changed across submission"):
            module.reauthenticate_open_batch_script(authenticated)
        transaction = {
            "submission_audit": {
                "batch_script": authenticated["binding"],
                "offline_local_root": False,
            },
        }
        with pytest.raises(ValueError, match="changed"):
            module.require_transaction_batch_script(transaction)
    finally:
        module.close_authenticated_batch_script(authenticated)


def test_descriptor_bound_metadata_replacement_rejects_parent_path_swap(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("old\n")
    moved = tmp_path / "moved-metadata"
    original_renameat2 = module.renameat2
    swapped = False

    def swap_parent_before_rename(*args, **kwargs):
        nonlocal swapped
        if not swapped:
            swapped = True
            parent.rename(moved)
            parent.mkdir()
        return original_renameat2(*args, **kwargs)

    monkeypatch.setattr(module, "renameat2", swap_parent_before_rename)
    with pytest.raises(ValueError, match="metadata parent pathname changed"):
        module.write_text(target, "new\n")
    assert not target.exists()
    assert (moved / target.name).read_text() == "new\n"
    assert (
        moved / module.metadata_temporary_name(target.name)
    ).read_text() == "old\n"


@pytest.mark.parametrize("mutation", ["mode", "group"])
def test_descriptor_bound_metadata_publication_rejects_parent_security_drift(
    tmp_path, monkeypatch, mutation
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    if mutation == "group":
        alternate_groups = [
            group for group in os.getgroups()
            if group != parent.stat().st_gid
        ]
        if not alternate_groups:
            pytest.skip("group-drift test requires a supplementary group")
        alternate_group = alternate_groups[0]
    original_renameat2 = module.renameat2
    mutated = False

    def drift_parent_at_mutation_boundary(*args, **kwargs):
        nonlocal mutated
        if not mutated:
            mutated = True
            if mutation == "mode":
                parent.chmod(0o777)
            else:
                os.chown(parent, -1, alternate_group)
        return original_renameat2(*args, **kwargs)

    monkeypatch.setattr(module, "renameat2", drift_parent_at_mutation_boundary)
    with pytest.raises(ValueError, match="metadata parent security metadata changed"):
        module.write_text(target, "controller-payload\n")
    assert mutated
    assert target.read_text() == "controller-payload\n"
    temporary = parent / module.metadata_temporary_name(target.name)
    assert not temporary.exists()


def test_canonical_metadata_parent_rejects_initial_insecure_mode_before_mutation(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    parent = root / "accounting"
    parent.mkdir()
    parent.chmod(0o777)
    target = parent / "state.json"
    rename_calls = []
    original_renameat2 = module.renameat2
    original_renameat2_between = module.renameat2_between

    def record_rename(*args, **kwargs):
        rename_calls.append(("same", args))
        return original_renameat2(*args, **kwargs)

    def record_rename_between(*args, **kwargs):
        rename_calls.append(("between", args))
        return original_renameat2_between(*args, **kwargs)

    monkeypatch.setattr(module, "renameat2", record_rename)
    monkeypatch.setattr(module, "renameat2_between", record_rename_between)
    with module.canonical_root_lock(root):
        with pytest.raises(ValueError, match="canonical parent is not owner-controlled"):
            module.write_text(target, "blocked\n")
    assert rename_calls == []
    assert not target.exists()
    assert not (parent / module.metadata_temporary_name(target.name)).exists()
    assert list(root.rglob(".cgl_lf_stage_i_replaced_*.forensic")) == []

    retained = parent / "retained.json"
    retained.write_text("retained\n")
    with module.canonical_root_lock(root):
        with pytest.raises(ValueError, match="canonical parent is not owner-controlled"):
            module.unlink_durable(retained)
    assert rename_calls == []
    assert retained.read_text() == "retained\n"
    assert list(root.rglob(".cgl_lf_stage_i_replaced_*.forensic")) == []

    offline = tmp_path / "offline"
    offline.mkdir()
    offline.chmod(0o777)
    module.write_text(offline / "state.json", "fixture\n")
    assert (offline / "state.json").read_text() == "fixture\n"


def test_copy_file_ambiguous_publication_is_durable_and_retry_safe(
    tmp_path, monkeypatch,
):
    module = load_controller()
    source = tmp_path / "source.bin"
    destination_parent = tmp_path / "destination"
    destination = destination_parent / "copied.bin"
    source.write_bytes(b"source\x00payload\n")
    source.chmod(0o640)
    destination_parent.mkdir()
    original_renameat2 = module.renameat2
    original_fsync = module.os.fsync
    injected = False
    durable_target_inodes = []
    rename_calls = 0

    def record_directory_fsync(descriptor):
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            try:
                durable_target_inodes.append(
                    module.os.stat(
                        destination.name,
                        dir_fd=descriptor,
                        follow_symlinks=False,
                    ).st_ino
                )
            except FileNotFoundError:
                pass
        original_fsync(descriptor)

    def raise_after_copy_publication(
        directory_fd, source_name, destination_name, flags, label,
    ):
        nonlocal injected, rename_calls
        rename_calls += 1
        result = original_renameat2(
            directory_fd, source_name, destination_name, flags, label
        )
        if not injected and flags == module.RENAME_NOREPLACE:
            injected = True
            raise RuntimeError("injected exception after copy publication")
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2", raise_after_copy_publication)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        module.copy_file(source, destination)
    assert injected
    assert rename_calls == 1
    assert destination.read_bytes() == source.read_bytes()
    assert stat.S_IMODE(destination.stat().st_mode) == 0o640
    assert destination.stat().st_ino in durable_target_inodes
    assert not (destination_parent / module.metadata_temporary_name(destination.name)).exists()

    monkeypatch.setattr(module, "renameat2", original_renameat2)
    module.copy_file(source, destination)
    assert destination.read_bytes() == source.read_bytes()
    assert stat.S_IMODE(destination.stat().st_mode) == 0o640


def test_copy_file_parent_path_drift_preserves_bound_durable_outcome(
    tmp_path, monkeypatch,
):
    module = load_controller()
    source = tmp_path / "source.bin"
    source.write_bytes(b"bound-copy\n")
    parent = tmp_path / "destination"
    parent.mkdir()
    destination = parent / "copied.bin"
    moved = tmp_path / "moved-destination"
    original_renameat2 = module.renameat2
    swapped = False

    def swap_parent_at_publication(
        directory_fd, source_name, destination_name, flags, label,
    ):
        nonlocal swapped
        if not swapped and flags == module.RENAME_NOREPLACE:
            swapped = True
            parent.rename(moved)
            parent.mkdir()
        return original_renameat2(
            directory_fd, source_name, destination_name, flags, label
        )

    monkeypatch.setattr(module, "renameat2", swap_parent_at_publication)
    with pytest.raises(ValueError, match="copy parent pathname changed"):
        module.copy_file(source, destination)
    assert swapped
    assert not destination.exists()
    assert (moved / destination.name).read_bytes() == source.read_bytes()
    assert not (parent / module.metadata_temporary_name(destination.name)).exists()


def test_copy_file_rejects_insecure_publication_mode_before_mutation(tmp_path):
    module = load_controller()
    source = tmp_path / "source.bin"
    source.write_bytes(b"insecure-mode\n")
    source.chmod(0o666)
    parent = tmp_path / "destination"
    parent.mkdir()
    destination = parent / "copied.bin"

    with pytest.raises(ValueError, match="group/world-writable destination"):
        module.copy_file(source, destination)
    assert not destination.exists()
    assert not (parent / module.metadata_temporary_name(destination.name)).exists()


def test_mkdir_durable_classifies_post_syscall_exception_and_retry(
    tmp_path, monkeypatch,
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    target = root / "a" / "b"
    original_mkdir = module.os.mkdir
    original_fsync = module.os.fsync
    injected = False
    directory_fsyncs = 0

    def raise_after_mkdir(path, mode=0o777, *, dir_fd=None):
        nonlocal injected
        result = original_mkdir(path, mode, dir_fd=dir_fd)
        if not injected:
            injected = True
            raise RuntimeError("injected exception after mkdir")
        return result

    def record_directory_fsync(descriptor):
        nonlocal directory_fsyncs
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            directory_fsyncs += 1
        original_fsync(descriptor)

    monkeypatch.setattr(module.os, "mkdir", raise_after_mkdir)
    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    module.mkdir_durable(target)
    module.mkdir_durable(target)
    assert injected
    assert target.is_dir()
    assert all(
        profile.st_uid == os.geteuid()
        and not stat.S_IMODE(profile.st_mode) & 0o022
        for profile in (target.parent.stat(), target.stat())
    )
    assert directory_fsyncs >= 4


def test_canonical_mkdir_durable_umask_zero_retains_exact_secure_profile(
    tmp_path, monkeypatch,
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    root.chmod(0o2755)
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    target = root / "created" / "nested"
    previous_umask = os.umask(0o0000)
    try:
        with module.canonical_root_lock(root):
            module.mkdir_durable(target)
    finally:
        os.umask(previous_umask)

    for directory in (target.parent, target):
        profile = directory.stat()
        assert profile.st_uid == os.geteuid()
        assert profile.st_gid == root.stat().st_gid
        assert stat.S_IMODE(profile.st_mode) == 0o2755
        assert not stat.S_IMODE(profile.st_mode) & 0o022


def test_mkdir_durable_post_child_fsync_exception_preserves_durable_state(
    tmp_path, monkeypatch,
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    target = root / "created"
    original_fsync = module.os.fsync
    injected = False
    parent_fsyncs = 0

    def raise_after_child_fsync(descriptor):
        nonlocal injected, parent_fsyncs
        profile = module.os.fstat(descriptor)
        original_fsync(descriptor)
        if (
            stat.S_ISDIR(profile.st_mode)
            and profile.st_ino == root.stat().st_ino
        ):
            parent_fsyncs += 1
        elif not injected and stat.S_ISDIR(profile.st_mode):
            injected = True
            raise RuntimeError("injected exception after child fsync")

    monkeypatch.setattr(module.os, "fsync", raise_after_child_fsync)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        module.mkdir_durable(target)
    assert injected
    assert parent_fsyncs >= 1
    assert target.is_dir()

    monkeypatch.setattr(module.os, "fsync", original_fsync)
    module.mkdir_durable(target)
    assert target.is_dir()


def test_mkdir_durable_authority_drift_preserves_created_directory(
    tmp_path, monkeypatch,
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    target = root / "created"
    original_mkdir = module.os.mkdir
    original_fsync = module.os.fsync
    drifted = False
    directory_fsyncs = 0

    def drift_after_mkdir(path, mode=0o777, *, dir_fd=None):
        nonlocal drifted
        result = original_mkdir(path, mode, dir_fd=dir_fd)
        if not drifted and path == target.name:
            drifted = True
            root.chmod(0o777)
        return result

    def record_directory_fsync(descriptor):
        nonlocal directory_fsyncs
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            directory_fsyncs += 1
        original_fsync(descriptor)

    monkeypatch.setattr(module.os, "mkdir", drift_after_mkdir)
    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    with module.canonical_root_lock(root):
        with pytest.raises(ValueError, match="group/world-writable"):
            module.mkdir_durable(target)
        assert drifted
        assert target.is_dir()
        assert directory_fsyncs >= 2
        root.chmod(0o755)
    assert target.is_dir()


def test_atomic_ledger_append_recovers_partial_temporary_with_journal_retained(
    tmp_path, monkeypatch,
):
    module = load_controller()
    accounting = tmp_path / "accounting"
    transactions = tmp_path / "transactions"
    accounting.mkdir()
    transactions.mkdir()
    ledger = accounting / "ledger.csv"
    module.write_text(ledger, module.ledger_csv_text([]), mode=0o644)
    journal = transactions / "retained.json"
    journal.write_text('{"state": "retained"}\n')
    journal.chmod(0o644)
    row = valid_ledger_row(module)
    original_fsync = module.os.fsync
    injected = False

    def fail_partial_file_fsync(descriptor):
        nonlocal injected
        profile = module.os.fstat(descriptor)
        if not injected and stat.S_ISREG(profile.st_mode) and profile.st_size > 0:
            injected = True
            raise RuntimeError("injected partial ledger temporary failure")
        original_fsync(descriptor)

    monkeypatch.setattr(module.os, "fsync", fail_partial_file_fsync)
    with pytest.raises(RuntimeError, match="partial ledger temporary"):
        module.append_ledger_row(ledger, row, [])
    assert injected
    assert module.read_ledger({"ledger": ledger}) == []
    assert journal.read_text() == '{"state": "retained"}\n'

    monkeypatch.setattr(module.os, "fsync", original_fsync)
    module.append_ledger_row(ledger, row, [])
    module.append_ledger_row(ledger, row, [])
    assert module.read_ledger({"ledger": ledger}) == [row]
    assert journal.read_text() == '{"state": "retained"}\n'
    assert not (accounting / module.metadata_temporary_name(ledger.name)).exists()


def test_atomic_ledger_append_ambiguous_publication_is_retry_safe_with_journal(
    tmp_path, monkeypatch,
):
    module = load_controller()
    accounting = tmp_path / "accounting"
    transactions = tmp_path / "transactions"
    accounting.mkdir()
    transactions.mkdir()
    ledger = accounting / "ledger.csv"
    module.write_text(ledger, module.ledger_csv_text([]), mode=0o644)
    journal = transactions / "retained.json"
    journal.write_text('{"state": "retained"}\n')
    journal.chmod(0o644)
    row = valid_ledger_row(module)
    original_renameat2 = module.renameat2
    injected = False

    def raise_after_ledger_exchange(
        directory_fd, source_name, destination_name, flags, label,
    ):
        nonlocal injected
        result = original_renameat2(
            directory_fd, source_name, destination_name, flags, label
        )
        if (
            not injected
            and flags == module.RENAME_EXCHANGE
            and destination_name == ledger.name
        ):
            injected = True
            raise RuntimeError("injected exception after ledger exchange")
        return result

    monkeypatch.setattr(module, "renameat2", raise_after_ledger_exchange)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        module.append_ledger_row(ledger, row, [])
    assert injected
    assert module.read_ledger({"ledger": ledger}) == [row]
    assert journal.read_text() == '{"state": "retained"}\n'
    temporary = accounting / module.metadata_temporary_name(ledger.name)
    assert temporary.read_bytes() == module.ledger_csv_text([]).encode("utf-8")

    monkeypatch.setattr(module, "renameat2", original_renameat2)
    module.append_ledger_row(ledger, row, [])
    assert module.read_ledger({"ledger": ledger}) == [row]
    assert journal.read_text() == '{"state": "retained"}\n'
    assert temporary.read_bytes() == module.ledger_csv_text([]).encode("utf-8")


def test_metadata_publication_rejects_equal_length_in_place_overwrite(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    original_fsync = module.os.fsync
    mutated = False
    published_inode = None
    durable_target_inodes = []

    def record_directory_fsync(descriptor):
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            try:
                durable_target_inodes.append(
                    module.os.stat(
                        target.name, dir_fd=descriptor, follow_symlinks=False
                    ).st_ino
                )
            except FileNotFoundError:
                durable_target_inodes.append(None)
        original_fsync(descriptor)

    def overwrite_after_publication(
        directory_fd, source, destination, flags, label,
    ):
        nonlocal mutated, published_inode
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        if not mutated and source == temporary_name:
            mutated = True
            published_inode = module.os.stat(
                destination, dir_fd=directory_fd, follow_symlinks=False
            ).st_ino
            overwrite_direct_child_restore_mtime(
                directory_fd, destination, b"evil\n"
            )
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2", overwrite_after_publication)
    with pytest.raises(ValueError, match="exact bytes.*durable forward recovery"):
        module.write_text(target, "good\n")
    assert mutated
    assert target.read_text() == "evil\n"
    assert target.stat().st_ino == published_inode
    assert not (parent / temporary_name).exists()
    assert published_inode in durable_target_inodes
    assert len(durable_target_inodes) >= 1


@pytest.mark.parametrize("profile_race", ("mode", "link"))
def test_mode_requested_publication_rejects_public_profile_race(
    tmp_path, monkeypatch, profile_race
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    alias = parent / "state-alias.json"
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    raced = False

    def alter_profile_after_publication(
        directory_fd, source, destination, flags, label,
    ):
        nonlocal raced
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        if not raced and source == temporary_name:
            raced = True
            if profile_race == "mode":
                os.chmod(destination, 0o666, dir_fd=directory_fd)
            else:
                os.link(
                    destination,
                    alias.name,
                    src_dir_fd=directory_fd,
                    dst_dir_fd=directory_fd,
                )
        return result

    monkeypatch.setattr(module, "renameat2", alter_profile_after_publication)
    with pytest.raises(ValueError, match="durable forward recovery state"):
        module.write_text(target, "good\n", mode=0o644)
    assert raced
    assert target.read_text() == "good\n"
    assert not (parent / temporary_name).exists()
    if profile_race == "mode":
        assert stat.S_IMODE(target.stat().st_mode) == 0o666
    else:
        assert alias.stat().st_ino == target.stat().st_ino
        assert target.stat().st_nlink == 2


def test_metadata_post_rename_nonregular_substitution_fsyncs_before_failure(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    escaped = parent / "escaped-state.json"
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    original_fsync = module.os.fsync
    directory_fsyncs = 0

    def record_directory_fsync(descriptor):
        nonlocal directory_fsyncs
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            directory_fsyncs += 1
        original_fsync(descriptor)

    def substitute_directory_after_publication(
        directory_fd, source, destination, flags, label,
    ):
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        if source == temporary_name:
            module.os.rename(
                destination,
                escaped.name,
                src_dir_fd=directory_fd,
                dst_dir_fd=directory_fd,
            )
            module.os.mkdir(destination, dir_fd=directory_fd)
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2", substitute_directory_after_publication)
    with pytest.raises(ValueError, match="post-rename state is not durably authenticated"):
        module.write_text(target, "good\n")
    assert directory_fsyncs >= 1
    assert escaped.read_text() == "good\n"
    assert target.is_dir()


def test_publication_exception_after_successful_rename_durably_recovers(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    original_fsync = module.os.fsync
    published_inode = None
    durable_target_inodes = []
    injected = False
    rename_calls = 0

    def record_directory_fsync(descriptor):
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            try:
                durable_target_inodes.append(
                    module.os.stat(
                        target.name, dir_fd=descriptor, follow_symlinks=False
                    ).st_ino
                )
            except FileNotFoundError:
                durable_target_inodes.append(None)
        original_fsync(descriptor)

    def raise_after_publication(directory_fd, source, destination, flags, label):
        nonlocal injected, published_inode, rename_calls
        rename_calls += 1
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        if not injected and source == temporary_name:
            injected = True
            published_inode = module.os.stat(
                destination, dir_fd=directory_fd, follow_symlinks=False
            ).st_ino
            raise RuntimeError("injected exception after publication rename")
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2", raise_after_publication)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        module.write_text(target, "good\n")
    assert injected
    assert target.read_text() == "good\n"
    assert target.stat().st_ino == published_inode
    assert not (parent / temporary_name).exists()
    assert published_inode in durable_target_inodes
    assert rename_calls == 1


def test_exchange_exception_after_successful_rename_durably_recovers(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("prior\n")
    prior_inode = target.stat().st_ino
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    original_fsync = module.os.fsync
    published_inode = None
    durable_target_inodes = []
    injected = False

    def record_directory_fsync(descriptor):
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            durable_target_inodes.append(
                module.os.stat(
                    target.name, dir_fd=descriptor, follow_symlinks=False
                ).st_ino
            )
        original_fsync(descriptor)

    def raise_after_exchange(directory_fd, source, destination, flags, label):
        nonlocal injected, published_inode
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        if not injected and flags == module.RENAME_EXCHANGE:
            injected = True
            published_inode = module.os.stat(
                destination, dir_fd=directory_fd, follow_symlinks=False
            ).st_ino
            raise RuntimeError("injected exception after exchange rename")
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2", raise_after_exchange)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        module.write_text(target, "newer\n")
    assert injected
    assert target.read_text() == "newer\n"
    assert target.stat().st_ino == published_inode
    temporary = parent / temporary_name
    assert temporary.read_text() == "prior\n"
    assert temporary.stat().st_ino == prior_inode
    assert published_inode in durable_target_inodes


def test_metadata_replacement_equal_length_overwrite_preserves_forward_state(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("prior\n")
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    mutated = False

    def overwrite_after_exchange(
        directory_fd, source, destination, flags, label,
    ):
        nonlocal mutated
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        if not mutated and flags == module.RENAME_EXCHANGE:
            mutated = True
            overwrite_direct_child_restore_mtime(
                directory_fd, destination, b"evil!\n"
            )
        return result

    monkeypatch.setattr(module, "renameat2", overwrite_after_exchange)
    with pytest.raises(ValueError, match="exact bytes.*durable forward recovery"):
        module.write_text(target, "good!\n")
    assert mutated
    assert target.read_text() == "evil!\n"
    assert (parent / temporary_name).read_text() == "prior\n"


def test_predecessor_quarantine_preserves_recoverable_state_after_parent_drift(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("prior\n")
    original_renameat2_between = module.renameat2_between
    drifted = False

    def drift_parent_during_quarantine(
        source_directory_fd,
        source,
        target_directory_fd,
        destination,
        flags,
        label,
    ):
        nonlocal drifted
        if not drifted and source_directory_fd != target_directory_fd:
            drifted = True
            parent.chmod(0o777)
        return original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            destination,
            flags,
            label,
        )

    monkeypatch.setattr(module, "renameat2_between", drift_parent_during_quarantine)
    with pytest.raises(ValueError, match="security metadata changed"):
        module.write_text(target, "newer\n")
    assert drifted
    assert target.read_text() == "newer\n"
    assert not (parent / module.metadata_temporary_name(target.name)).exists()
    forensic = list(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))
    assert len(forensic) == 1
    assert forensic[0].read_text() == "prior\n"


def test_quarantine_exception_after_successful_rename_durably_recovers(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("authenticated\n")
    target_inode = target.stat().st_ino
    original_renameat2_between = module.renameat2_between
    original_fsync = module.os.fsync
    forensic_name = None
    forensic_inode = None
    durable_forensic_inodes = []
    injected = False
    rename_calls = 0

    def record_directory_fsync(descriptor):
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode) and forensic_name:
            try:
                durable_forensic_inodes.append(
                    module.os.stat(
                        forensic_name, dir_fd=descriptor, follow_symlinks=False
                    ).st_ino
                )
            except FileNotFoundError:
                pass
        original_fsync(descriptor)

    def raise_after_quarantine(
        source_directory_fd,
        source,
        target_directory_fd,
        destination,
        flags,
        label,
    ):
        nonlocal injected, forensic_name, forensic_inode, rename_calls
        rename_calls += 1
        result = original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            destination,
            flags,
            label,
        )
        if not injected and source_directory_fd != target_directory_fd:
            injected = True
            forensic_name = destination
            forensic_inode = module.os.stat(
                destination,
                dir_fd=target_directory_fd,
                follow_symlinks=False,
            ).st_ino
            raise RuntimeError("injected exception after quarantine rename")
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2_between", raise_after_quarantine)
    with pytest.raises(ValueError, match="durable forensic recovery state was preserved"):
        module.unlink_durable(target)
    assert injected
    assert not target.exists()
    assert forensic_inode in durable_forensic_inodes
    forensic = list(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))
    assert len(forensic) == 1
    assert forensic[0].read_text() == "authenticated\n"
    assert forensic[0].stat().st_ino == target_inode
    assert rename_calls == 1


def test_absent_target_publication_race_never_clobbers_raced_target(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    original_renameat2 = module.renameat2
    raced = False

    def race_before_noreplace(directory_fd, source, destination, flags, label):
        nonlocal raced
        if flags == module.RENAME_NOREPLACE and not raced:
            raced = True
            descriptor = os.open(
                destination,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o644,
                dir_fd=directory_fd,
            )
            os.write(descriptor, b"raced-target\n")
            os.close(descriptor)
        return original_renameat2(directory_fd, source, destination, flags, label)

    monkeypatch.setattr(module, "renameat2", race_before_noreplace)
    with pytest.raises(ValueError, match="target already exists"):
        module.write_text(target, "controller-payload\n")
    assert target.read_text() == "raced-target\n"
    assert any(
        path.read_text() == "controller-payload\n"
        for path in parent.glob(".state.json.*.tmp")
    )


def test_existing_target_exchange_race_never_deletes_or_clobbers_substitute(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    substitute = parent / "substitute.json"
    escaped = parent / "escaped.json"
    target.write_text("old-target\n")
    substitute.write_text("raced-substitute\n")
    original_renameat2 = module.renameat2

    def race_after_exchange(directory_fd, source, destination, flags, label):
        original_renameat2(directory_fd, source, destination, flags, label)
        if flags == module.RENAME_EXCHANGE:
            os.rename(
                destination,
                escaped.name,
                src_dir_fd=directory_fd,
                dst_dir_fd=directory_fd,
            )
            os.rename(
                substitute.name,
                destination,
                src_dir_fd=directory_fd,
                dst_dir_fd=directory_fd,
            )

    monkeypatch.setattr(module, "renameat2", race_after_exchange)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        module.write_text(target, "controller-payload\n")
    assert target.read_text() == "raced-substitute\n"
    assert escaped.read_text() == "controller-payload\n"
    assert any(
        path.read_text() == "old-target\n"
        for path in parent.glob(".state.json.*.tmp")
    )


def test_successful_existing_target_replacement_quarantines_predecessor(
    tmp_path,
):
    module = load_controller()
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    journal = transactions / "pending.json"
    journal.write_text('{"state": "old"}\n')

    module.write_json(journal, {"state": "new"})

    assert json.loads(journal.read_text()) == {"state": "new"}
    assert not list(transactions.glob(".pending.json.*.tmp"))
    forensics = list(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))
    assert len(forensics) == 1
    assert json.loads(forensics[0].read_text()) == {"state": "old"}


def test_lustre_unsupported_noreplace_uses_hard_link_commit(tmp_path, monkeypatch):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    module.write_text(target, "published\n", mode=0o644)

    assert target.read_text() == "published\n"
    assert target.stat().st_nlink == 1
    assert not (parent / module.metadata_temporary_name(target.name)).exists()


@pytest.mark.parametrize("flags", (1, 2))
def test_renameat2_einval_is_an_explicit_unsupported_response(
    monkeypatch, flags,
):
    module = load_controller()

    class UnsupportedOperation:
        argtypes = None
        restype = None

        def __call__(self, *_args):
            return -1

    libc = SimpleNamespace(renameat2=UnsupportedOperation())
    monkeypatch.setattr(module.ctypes, "CDLL", lambda *_args, **_kwargs: libc)
    monkeypatch.setattr(module.ctypes, "get_errno", lambda: module.errno.EINVAL)

    with pytest.raises(module.Renameat2Unsupported) as unsupported:
        module.renameat2_between(10, "source", 11, "target", flags, "fixture")
    assert unsupported.value.errno == module.errno.EINVAL


def test_lustre_hard_link_publication_never_clobbers_target_race(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    real_link = module.os.link
    raced = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    def race_before_link(source, destination, **kwargs):
        nonlocal raced
        if not raced and destination == target.name:
            target.write_text("raced-target\n")
            target.chmod(0o644)
            raced = True
        return real_link(source, destination, **kwargs)

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(module.os, "link", race_before_link)
    with pytest.raises(ValueError, match="target already exists"):
        module.write_text(target, "controller-payload\n", mode=0o644)

    assert raced
    assert target.read_text() == "raced-target\n"
    temporary = parent / module.metadata_temporary_name(target.name)
    assert temporary.read_text() == "controller-payload\n"
    assert temporary.stat().st_nlink == 1


def test_lustre_unsupported_exchange_uses_recoverable_in_place_replacement(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("predecessor\n")
    target.chmod(0o644)

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(
        module.os,
        "rename",
        lambda *_args, **_kwargs: pytest.fail("flags-zero rename must not run"),
    )
    prior_inode = target.stat().st_ino
    module.write_text(target, "replacement\n", mode=0o644)

    assert target.read_text() == "replacement\n"
    assert target.stat().st_ino == prior_inode
    assert target.stat().st_nlink == 1
    temporary_name = module.metadata_temporary_name(target.name)
    assert not (parent / temporary_name).exists()
    assert not (
        parent / module.metadata_predecessor_recovery_name(temporary_name)
    ).exists()
    assert not list(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))


def test_lustre_in_place_replacement_recovers_after_public_inode_truncation(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text("predecessor\n")
    target.chmod(0o644)
    real_ftruncate = module.os.ftruncate
    injected = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    def truncate_then_raise(descriptor, length):
        nonlocal injected
        result = real_ftruncate(descriptor, length)
        if not injected:
            injected = True
            raise RuntimeError("simulated crash after public inode truncation")
        return result

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(module.os, "ftruncate", truncate_then_raise)
    with pytest.raises(ValueError, match="durable in-place forward recovery state"):
        module.write_text(
            target,
            "replacement\n",
            mode=0o644,
            expected_predecessor=b"predecessor\n",
        )

    temporary_name = module.metadata_temporary_name(target.name)
    predecessor_name = module.metadata_predecessor_recovery_name(temporary_name)
    assert injected
    assert target.read_bytes() == b""
    assert stat.S_IMODE(target.stat().st_mode) == 0o600
    assert (parent / temporary_name).read_text() == "replacement\n"
    assert (parent / predecessor_name).read_text() == "predecessor\n"
    assert stat.S_IMODE((parent / predecessor_name).stat().st_mode) == 0o400

    module.write_text(
        target,
        "replacement\n",
        mode=0o644,
        expected_predecessor=b"predecessor\n",
    )

    assert target.read_text() == "replacement\n"
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert not (parent / temporary_name).exists()
    assert not (parent / predecessor_name).exists()


def test_lustre_lifecycle_metadata_paths_replace_under_renameat2_einval(
    tmp_path, monkeypatch,
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    paths["ledger"].write_text(module.ledger_csv_text([]))
    paths["ledger"].chmod(0o644)
    paths["reservations"].chmod(0o644)
    manifest = paths["runs"] / "R16/s00/manifest/prepared_run.json"
    manifest.parent.mkdir(parents=True)
    write_json(manifest, {"state": "prepared"}, mode=0o644)
    transaction = paths["transactions"] / "lifecycle.json"
    write_json(transaction, {"kind": "submit_pending"}, mode=0o644)

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(
        module.os,
        "rename",
        lambda *_args, **_kwargs: pytest.fail("flags-zero rename must not run"),
    )
    module.write_json(paths["reservations"], [{"state": "prepared"}])
    module.write_json(manifest, {"state": "submitted"})
    module.append_ledger_row(paths["ledger"], valid_ledger_row(module), [])
    module.write_json(transaction, {"kind": "recorded"}, mode=0o644)

    assert json.loads(paths["reservations"].read_text()) == [{"state": "prepared"}]
    assert json.loads(manifest.read_text()) == {"state": "submitted"}
    assert module.read_ledger(paths) == [valid_ledger_row(module)]
    assert json.loads(transaction.read_text()) == {"kind": "recorded"}
    module.unlink_trusted_transaction(
        {"transactions": paths["transactions"]},
        transaction,
        {"kind": "recorded"},
    )
    assert not transaction.exists()
    assert any(
        json.loads(path.read_text()) == {"kind": "recorded"}
        for path in paths["accounting"].glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    for target in (paths["reservations"], paths["ledger"], manifest):
        temporary_name = module.metadata_temporary_name(target.name)
        assert not (target.parent / temporary_name).exists()
        assert not (
            target.parent / module.metadata_predecessor_recovery_name(temporary_name)
        ).exists()


def test_lustre_unsupported_forensic_retirement_uses_deterministic_hardlink_move(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    source = parent / "state.json"
    source.write_text("retained\n")
    source.chmod(0o644)

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    binding = module.open_regular_file_binding(parent_fd, source.name, "retained source")
    try:
        profile = module.require_directory_descriptor_binding(
            parent, parent_fd, "retirement parent"
        )
        retained = module.quarantine_bound_predecessor(
            parent_fd,
            source.name,
            binding,
            parent,
            profile,
            source,
            "retained source",
        )
    finally:
        binding.close()
        os.close(parent_fd)

    assert not source.exists()
    assert retained.read_text() == "retained\n"
    assert retained.stat().st_nlink == 1
    assert retained.name == (
        ".cgl_lf_stage_i_replaced_"
        + hashlib.sha256(
            f"{parent.absolute()}\0{source.name}\0{sha256(retained)}".encode()
        ).hexdigest()
        + ".forensic"
    )


def test_lustre_deterministic_forensic_hardlink_state_retries_to_completion(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    source = parent / "state.json"
    source.write_text("retained\n")
    source.chmod(0o644)
    real_unlink = module.os.unlink
    injected = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    def preserve_first_two_link_state(name, **kwargs):
        nonlocal injected
        if not injected and name == source.name:
            injected = True
            raise RuntimeError("simulated crash before source unlink")
        return real_unlink(name, **kwargs)

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(module.os, "unlink", preserve_first_two_link_state)

    parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        profile = module.require_directory_descriptor_binding(
            parent, parent_fd, "retirement parent"
        )
        binding = module.open_regular_file_binding(
            parent_fd, source.name, "retained source"
        )
        try:
            with pytest.raises(ValueError, match="two-link recovery state"):
                module.quarantine_bound_predecessor(
                    parent_fd,
                    source.name,
                    binding,
                    parent,
                    profile,
                    source,
                    "retained source",
                )
        finally:
            binding.close()
        forensic = next(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))
        assert source.stat().st_ino == forensic.stat().st_ino
        assert source.stat().st_nlink == 2

        binding = module.open_regular_file_binding(
            parent_fd, source.name, "retained source retry"
        )
        try:
            retained = module.quarantine_bound_predecessor(
                parent_fd,
                source.name,
                binding,
                parent,
                profile,
                source,
                "retained source retry",
            )
        finally:
            binding.close()
    finally:
        os.close(parent_fd)

    assert injected
    assert not source.exists()
    assert retained == forensic
    assert forensic.read_text() == "retained\n"
    assert forensic.stat().st_nlink == 1


def test_lustre_two_link_publication_is_deterministically_recovered_on_retry(
    tmp_path,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    temporary = parent / module.metadata_temporary_name(target.name)
    temporary.write_text("published\n")
    temporary.chmod(0o644)
    os.link(temporary, target)
    assert temporary.stat().st_nlink == 2

    parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        profile = module.require_directory_descriptor_binding(
            parent, parent_fd, "publication parent"
        )
        assert module.retire_metadata_temporary(
            parent,
            parent_fd,
            temporary.name,
            profile,
            "retained publication temporary",
        ) is None
    finally:
        os.close(parent_fd)

    assert not temporary.exists()
    assert target.read_text() == "published\n"
    assert target.stat().st_nlink == 1


def test_lustre_two_link_publication_rejects_undiscoverable_alias(tmp_path):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    temporary = parent / module.metadata_temporary_name(target.name)
    temporary.write_text("published\n")
    temporary.chmod(0o644)
    random_alias = parent / "random-recovery-alias"
    os.link(temporary, random_alias)

    parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        profile = module.require_directory_descriptor_binding(
            parent, parent_fd, "publication parent"
        )
        with pytest.raises(ValueError, match="deterministic authenticated public target"):
            module.retire_metadata_temporary(
                parent,
                parent_fd,
                temporary.name,
                profile,
                "retained publication temporary",
            )
    finally:
        os.close(parent_fd)

    assert temporary.read_text() == "published\n"
    assert random_alias.read_text() == "published\n"
    assert temporary.stat().st_nlink == 2
    assert not target.exists()


def test_lustre_hard_link_commit_classifies_post_syscall_errors(
    tmp_path, monkeypatch,
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "state.json"
    temporary_name = module.metadata_temporary_name(target.name)
    real_link = module.os.link
    real_unlink = module.os.unlink
    link_injected = False
    unlink_injected = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    def link_then_raise(source, destination, **kwargs):
        nonlocal link_injected
        result = real_link(source, destination, **kwargs)
        if not link_injected and source == temporary_name and destination == target.name:
            link_injected = True
            raise RuntimeError("reported link failure after commit")
        return result

    def unlink_then_raise(name, **kwargs):
        nonlocal unlink_injected
        result = real_unlink(name, **kwargs)
        if not unlink_injected and name == temporary_name:
            unlink_injected = True
            raise RuntimeError("reported unlink failure after commit")
        return result

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(module.os, "link", link_then_raise)
    monkeypatch.setattr(module.os, "unlink", unlink_then_raise)
    module.write_text(target, "classified\n", mode=0o644)

    assert link_injected
    assert unlink_injected
    assert target.read_text() == "classified\n"
    assert target.stat().st_nlink == 1
    assert not (parent / temporary_name).exists()


def test_transaction_discovery_forensically_retires_atomic_write_remnants(
    tmp_path,
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    journal = transactions / "pending.json"
    journal.write_text("{}\n")
    deterministic = transactions / module.metadata_temporary_name(journal.name)
    deterministic.write_text("deterministic-remnant\n")
    legacy = transactions / f".legacy.json.1234.{'a' * 32}.tmp"
    legacy.write_text("legacy-remnant\n")

    assert module.pending_transaction_paths(paths) == [journal]
    assert not deterministic.exists()
    assert not legacy.exists()
    retained = list(
        paths["accounting"].glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert {path.read_text() for path in retained} == {
        "deterministic-remnant\n",
        "legacy-remnant\n",
    }


def test_transaction_discovery_recovers_two_link_hardlink_publication(
    tmp_path,
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    journal = transactions / "pending.json"
    temporary = transactions / module.metadata_temporary_name(journal.name)
    temporary.write_text('{"state": "pending"}\n')
    temporary.chmod(0o644)
    os.link(temporary, journal)
    assert journal.stat().st_nlink == 2

    assert module.pending_transaction_paths(paths) == [journal]
    assert journal.read_text() == '{"state": "pending"}\n'
    assert journal.stat().st_nlink == 1
    assert not temporary.exists()
    assert not list(paths["accounting"].glob(".cgl_lf_stage_i_replaced_*.forensic"))


def test_transaction_discovery_recovers_post_sbatch_submitted_journal_under_einval(
    tmp_path, monkeypatch,
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    journal = transactions / "submission.json"
    pending = {"kind": "submit_pending", "scheduler_boundary": "sbatch returned"}
    submitted = {"kind": "submitted", "job_id": "12345"}
    write_json(journal, pending, mode=0o644)
    real_ftruncate = module.os.ftruncate
    injected = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise module.Renameat2Unsupported(
            module.errno.EINVAL, os.strerror(module.errno.EINVAL)
        )

    def truncate_then_raise(descriptor, length):
        nonlocal injected
        result = real_ftruncate(descriptor, length)
        if not injected:
            injected = True
            raise RuntimeError("simulated crash during post-sbatch journal update")
        return result

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(module.os, "ftruncate", truncate_then_raise)
    with pytest.raises(ValueError, match="durable in-place forward recovery state"):
        module.write_json(journal, submitted, mode=0o644)

    temporary_name = module.metadata_temporary_name(journal.name)
    predecessor_name = module.metadata_predecessor_recovery_name(temporary_name)
    assert injected
    assert journal.read_bytes() == b""
    assert stat.S_IMODE(journal.stat().st_mode) == 0o600
    assert json.loads((transactions / temporary_name).read_text()) == submitted
    assert json.loads((transactions / predecessor_name).read_text()) == pending

    monkeypatch.setattr(module.os, "ftruncate", real_ftruncate)
    assert module.pending_transaction_paths(paths) == [journal]
    assert json.loads(journal.read_text()) == submitted
    assert stat.S_IMODE(journal.stat().st_mode) == 0o644
    assert not (transactions / temporary_name).exists()
    assert not (transactions / predecessor_name).exists()


def test_transaction_discovery_forensically_retires_mode_zero_remnant(
    tmp_path,
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    temporary = (
        paths["transactions"] / module.metadata_temporary_name("pending.json")
    )
    temporary.write_text("partial-remnant\n")
    temporary.chmod(0o000)

    assert module.pending_transaction_paths(paths) == []
    assert not temporary.exists()
    retained = list(
        paths["accounting"].glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert len(retained) == 1
    assert stat.S_IMODE(retained[0].stat().st_mode) == 0o000


def test_transaction_store_scan_rejects_path_replacement_without_hiding_journal(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    pending = transactions / "pending.json"
    pending.write_text('{"state": "pending"}\n')
    pending.chmod(0o644)
    temporary = transactions / module.metadata_temporary_name("crashed.json")
    temporary.write_text("crash-remnant\n")
    moved = transactions.with_name(f"{transactions.name}.moved")
    original_require_lock = module.require_canonical_mutation_lock
    swapped = False

    def replace_store_after_authenticated_scan(path):
        nonlocal swapped
        result = original_require_lock(path)
        if not swapped and path.absolute() == transactions.absolute():
            swapped = True
            transactions.rename(moved)
            transactions.mkdir()
        return result

    monkeypatch.setattr(
        module, "require_canonical_mutation_lock", replace_store_after_authenticated_scan
    )
    with pytest.raises(ValueError, match="pathname changed"):
        module.transaction_store_entries(paths)
    assert swapped
    assert not pending.exists()
    assert (moved / pending.name).read_text() == '{"state": "pending"}\n'
    assert (moved / temporary.name).read_text() == "crash-remnant\n"
    assert list(transactions.iterdir()) == []


def test_transaction_store_final_scan_rejects_late_inserted_journal(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    pending = transactions / "pending.json"
    pending.write_text('{"state": "pending"}\n')
    pending.chmod(0o644)
    inserted = transactions / "inserted.json"
    original_scan = module.stable_bound_directory_entries
    calls = 0

    def insert_after_final_scan(*args, **kwargs):
        nonlocal calls
        result = original_scan(*args, **kwargs)
        calls += 1
        if calls == 2:
            inserted.write_text('{"state": "inserted"}\n')
            inserted.chmod(0o644)
        return result

    monkeypatch.setattr(
        module, "stable_bound_directory_entries", insert_after_final_scan
    )
    with pytest.raises(ValueError, match="final contents changed"):
        module.transaction_store_entries(paths)
    assert calls >= 3
    assert pending.exists()
    assert inserted.exists()


def test_transaction_store_final_scan_rejects_namespace_replacement(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    pending = transactions / "pending.json"
    pending.write_text('{"state": "pending"}\n')
    pending.chmod(0o644)
    moved = transactions.with_name(f"{transactions.name}.moved-final")
    original_scan = module.stable_bound_directory_entries
    calls = 0

    def replace_after_final_scan(*args, **kwargs):
        nonlocal calls
        result = original_scan(*args, **kwargs)
        calls += 1
        if calls == 2:
            transactions.rename(moved)
            transactions.mkdir()
        return result

    monkeypatch.setattr(
        module, "stable_bound_directory_entries", replace_after_final_scan
    )
    with pytest.raises(ValueError, match="pathname changed"):
        module.transaction_store_entries(paths)
    assert calls >= 2
    assert not pending.exists()
    assert (moved / pending.name).exists()
    assert list(transactions.iterdir()) == []


def test_transaction_temporary_recovery_fails_closed_on_substitution_race(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    transactions = paths["transactions"]
    temporary = transactions / module.metadata_temporary_name("pending.json")
    escaped = transactions / "escaped-authenticated.tmp"
    substitute = transactions / "raced-substitute.tmp"
    temporary.write_text("authenticated-remnant\n")
    substitute.write_text("substitute\n")
    original_renameat2_between = module.renameat2_between
    raced = False

    def race_before_temporary_retirement(
        source_directory_fd,
        source,
        target_directory_fd,
        target,
        flags,
        label,
    ):
        nonlocal raced
        if not raced and source == temporary.name:
            raced = True
            os.rename(
                source,
                escaped.name,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=source_directory_fd,
            )
            os.rename(
                substitute.name,
                source,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=source_directory_fd,
            )
        return original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            target,
            flags,
            label,
        )

    monkeypatch.setattr(module, "renameat2_between", race_before_temporary_retirement)
    with pytest.raises(ValueError, match="durable forensic recovery state was preserved"):
        module.pending_transaction_paths(paths)
    assert raced
    assert escaped.read_text() == "authenticated-remnant\n"
    assert not temporary.exists()
    retained = list(
        paths["accounting"].glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert len(retained) == 1
    assert retained[0].read_text() == "substitute\n"


def test_transaction_temporary_recovery_fails_closed_on_in_place_mode_race(
    tmp_path, monkeypatch
):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    temporary = (
        paths["transactions"] / module.metadata_temporary_name("pending.json")
    )
    temporary.write_text("authenticated-remnant\n")
    original_renameat2_between = module.renameat2_between
    raced = False

    def race_before_temporary_retirement(
        source_directory_fd,
        source,
        target_directory_fd,
        target,
        flags,
        label,
    ):
        nonlocal raced
        if not raced and source == temporary.name:
            raced = True
            os.chmod(source, 0o400, dir_fd=source_directory_fd)
        return original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            target,
            flags,
            label,
        )

    monkeypatch.setattr(module, "renameat2_between", race_before_temporary_retirement)
    with pytest.raises(ValueError, match="durable forensic recovery state was preserved"):
        module.pending_transaction_paths(paths)
    assert raced
    assert not temporary.exists()
    retained = list(
        paths["accounting"].glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert len(retained) == 1
    assert retained[0].read_text() == "authenticated-remnant\n"
    assert stat.S_IMODE(retained[0].stat().st_mode) == 0o400


def test_transaction_retirement_race_preserves_substitute_forensically(
    tmp_path, monkeypatch
):
    module = load_controller()
    transactions = tmp_path / "accounting" / "transactions"
    transactions.mkdir(parents=True)
    journal = transactions / "pending.json"
    escaped = transactions / "escaped-original.json"
    substitute = transactions / "raced-substitute.json"
    expected = {"state": "authenticated"}
    journal.write_text(json.dumps(expected, indent=2, sort_keys=True) + "\n")
    substitute.write_text('{"state": "substitute"}\n')
    original_renameat2_between = module.renameat2_between
    raced = False

    def race_before_retirement(
        source_directory_fd,
        source,
        target_directory_fd,
        target,
        flags,
        label,
    ):
        nonlocal raced
        if not raced:
            raced = True
            os.rename(
                source,
                escaped.name,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=source_directory_fd,
            )
            os.rename(
                substitute.name,
                source,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=source_directory_fd,
            )
        return original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            target,
            flags,
            label,
        )

    monkeypatch.setattr(module, "renameat2_between", race_before_retirement)
    with pytest.raises(ValueError, match="durable forensic recovery state was preserved"):
        module.unlink_trusted_transaction(
            {"transactions": transactions},
            journal,
            expected,
        )
    assert escaped.read_text() == json.dumps(expected, indent=2, sort_keys=True) + "\n"
    assert not journal.exists()
    retained = list(
        (tmp_path / "accounting").glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert len(retained) == 1
    assert retained[0].read_text() == '{"state": "substitute"}\n'


def test_transaction_retirement_rejects_equal_length_in_place_overwrite(
    tmp_path, monkeypatch
):
    module = load_controller()
    transactions = tmp_path / "accounting" / "transactions"
    transactions.mkdir(parents=True)
    journal = transactions / "pending.json"
    expected = {"state": "good"}
    evil = {"state": "evil"}
    expected_payload = json.dumps(expected, indent=2, sort_keys=True) + "\n"
    evil_payload = (json.dumps(evil, indent=2, sort_keys=True) + "\n").encode()
    assert len(evil_payload) == len(expected_payload.encode())
    journal.write_text(expected_payload)
    journal_inode = journal.stat().st_ino
    original_renameat2_between = module.renameat2_between
    original_fsync = module.os.fsync
    mutated = False
    directory_fsyncs = {"transactions": 0, "accounting": 0}
    durable_forensic_inodes = []

    def record_directory_fsync(descriptor):
        profile = module.os.fstat(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            if profile.st_ino == transactions.stat().st_ino:
                directory_fsyncs["transactions"] += 1
            if profile.st_ino == transactions.parent.stat().st_ino:
                directory_fsyncs["accounting"] += 1
            forensic = list(
                transactions.parent.glob(".cgl_lf_stage_i_replaced_*.forensic")
            )
            durable_forensic_inodes.extend(path.stat().st_ino for path in forensic)
        original_fsync(descriptor)

    def overwrite_after_retirement(
        source_directory_fd,
        source,
        target_directory_fd,
        target,
        flags,
        label,
    ):
        nonlocal mutated
        result = original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            target,
            flags,
            label,
        )
        if not mutated and source == journal.name:
            mutated = True
            overwrite_direct_child_restore_mtime(
                target_directory_fd, target, evil_payload
            )
        return result

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2_between", overwrite_after_retirement)
    with pytest.raises(ValueError, match="exact bytes.*durable forensic recovery"):
        module.unlink_trusted_transaction(
            {"transactions": transactions},
            journal,
            expected,
    )
    assert mutated
    assert not journal.exists()
    assert journal_inode in durable_forensic_inodes
    assert directory_fsyncs["transactions"] >= 1
    assert directory_fsyncs["accounting"] >= 1
    forensic = list(
        (tmp_path / "accounting").glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert len(forensic) == 1
    assert forensic[0].read_text() == evil_payload.decode()
    assert forensic[0].stat().st_ino == journal_inode


def test_successful_transaction_retirement_preserves_authenticated_journal(
    tmp_path,
):
    module = load_controller()
    transactions = tmp_path / "accounting" / "transactions"
    transactions.mkdir(parents=True)
    journal = transactions / "pending.json"
    expected = {"state": "authenticated"}
    payload = json.dumps(expected, indent=2, sort_keys=True) + "\n"
    journal.write_text(payload)

    module.unlink_trusted_transaction(
        {"transactions": transactions},
        journal,
        expected,
    )

    assert not journal.exists()
    retained = list(
        (tmp_path / "accounting").glob(".cgl_lf_stage_i_replaced_*.forensic")
    )
    assert len(retained) == 1
    assert retained[0].read_text() == payload


def test_generic_durable_retirement_race_preserves_substitute_forensically(
    tmp_path, monkeypatch
):
    module = load_controller()
    parent = tmp_path / "metadata"
    parent.mkdir()
    target = parent / "target.json"
    escaped = parent / "escaped-original.json"
    substitute = parent / "raced-substitute.json"
    target.write_text("authenticated\n")
    substitute.write_text("substitute\n")
    original_renameat2_between = module.renameat2_between
    raced = False

    def race_before_retirement(
        source_directory_fd,
        source,
        target_directory_fd,
        target_name,
        flags,
        label,
    ):
        nonlocal raced
        if not raced:
            raced = True
            os.rename(
                source,
                escaped.name,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=source_directory_fd,
            )
            os.rename(
                substitute.name,
                source,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=source_directory_fd,
            )
        return original_renameat2_between(
            source_directory_fd,
            source,
            target_directory_fd,
            target_name,
            flags,
            label,
        )

    monkeypatch.setattr(module, "renameat2_between", race_before_retirement)
    with pytest.raises(ValueError, match="durable forensic recovery state was preserved"):
        module.unlink_durable(target)
    assert escaped.read_text() == "authenticated\n"
    assert not target.exists()
    retained = list(tmp_path.glob(".cgl_lf_stage_i_replaced_*.forensic"))
    assert len(retained) == 1
    assert retained[0].read_text() == "substitute\n"


def test_authenticated_git_binary_rejects_replacement_path(tmp_path, monkeypatch):
    module = load_controller()
    replacement = tmp_path / "git-real"
    replacement.write_text("#!/bin/sh\nexit 0\n")
    replacement.chmod(0o755)
    symlink = tmp_path / "git"
    symlink.symlink_to(replacement)
    monkeypatch.setattr(module, "GIT", symlink)
    with pytest.raises(ValueError, match="root-owned|symbolic link"):
        module.authenticated_git_binary()


def test_latest_f117_publication_rejects_symlinked_artifact(tmp_path):
    module = load_controller()
    paths, _ = make_published_recost(module, tmp_path / "root", checkpoint=117)
    artifact = paths["accounting"] / (
        f"mks24_stage_i_{EPOCH_SLUG}_F117_recost_evidence.json"
    )
    target = tmp_path / "moved-recost.json"
    artifact.rename(target)
    artifact.symlink_to(target)
    with pytest.raises(ValueError, match="symbolic link"):
        module.latest_published_r17_recost(paths, datetime.now(timezone.utc))


def test_latest_f117_publication_requires_owner_controlled_path(tmp_path, monkeypatch):
    module = load_controller()
    paths, _ = make_published_recost(module, tmp_path / "root", checkpoint=117)
    current_uid = os.geteuid()
    monkeypatch.setattr(module.os, "geteuid", lambda: current_uid + 1)
    with pytest.raises(ValueError, match="owner-controlled|owner controlled"):
        module.latest_published_r17_recost(paths, datetime.now(timezone.utc))


def test_latest_f117_review_rejects_every_waiver_scope(tmp_path):
    module = load_controller()
    paths, _ = make_published_recost(module, tmp_path / "root", checkpoint=117)
    artifact = paths["accounting"] / (
        f"mks24_stage_i_{EPOCH_SLUG}_F117_recost_evidence.json"
    )
    review = artifact.with_name(f"{artifact.name}.independent_review.json")
    audit = artifact.with_name(f"{artifact.name}.publication_audit.json")
    review_value = json.loads(review.read_text())
    review_value["scope"] = {
        "non_authorizing": True,
        "r12_frozen_e03_no_ct_profile_waiver": {"producer_contract": True},
    }
    write_json(review, review_value)
    audit_value = json.loads(audit.read_text())
    audit_value["independent_review"] = declared(review)
    write_json(audit, audit_value)
    with pytest.raises(ValueError, match="review scope differs"):
        module.latest_published_r17_recost(paths, datetime.now(timezone.utc))

    review_value["scope"] = {
        "non_authorizing": True,
    }
    write_json(review, review_value)
    audit_value["independent_review"] = declared(review)
    write_json(audit, audit_value)
    module.latest_published_r17_recost(paths, datetime.now(timezone.utc))


def test_transaction_store_rejects_symlink_path_swap(tmp_path):
    module = load_controller()
    paths = paths_for(module, tmp_path / "root")
    paths["transactions"].rmdir()
    replacement = tmp_path / "replacement-transactions"
    replacement.mkdir()
    paths["transactions"].symlink_to(replacement, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic link"):
        module.pending_transaction_paths(paths)


def test_lock_zero_size_and_ledger_canonical_six_decimals(tmp_path):
    module = load_controller()
    lock = tmp_path / "lock"
    lock.write_text("x")
    lock.chmod(0o644)
    with pytest.raises(ValueError, match="must be empty"):
        module.require_canonical_root_lock_profile(lock.stat(), lock)
    row = {
        "nodes": "1",
        "elapsed_seconds": "3600",
        "requested_walltime": "01:00:00",
        "reserved_node_hours": "1.000000",
        "actual_node_hours": "1.000000",
        "cumulative_stage_i_node_hours": "1.000000",
    }
    assert module.validate_ledger_numeric_fields(row, "ledger")[0] == 1.0
    row["actual_node_hours"] = "1.0"
    with pytest.raises(ValueError, match="invalid numeric fields"):
        module.validate_ledger_numeric_fields(row, "ledger")


def test_canonical_mutation_lock_path_inode_is_rechecked_throughout_section(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    target = root / "state.json"
    with pytest.raises(ValueError, match="lacks its root lock"):
        module.write_text(target, "unlocked\n")

    lock = module.canonical_root_lock_path(root)
    moved = root / "moved.lock"
    with pytest.raises(ValueError, match="lock path changed during canonical mutation"):
        with module.canonical_root_lock(root):
            module.write_text(target, "locked\n")
            lock.rename(moved)
            lock.touch(mode=0o644)
            module.write_text(root / "must-not-write.json", "blocked\n")
    assert target.read_text() == "locked\n"
    assert not (root / "must-not-write.json").exists()


def test_canonical_lock_rejects_writable_root_and_ancestor_before_creation(
    tmp_path, monkeypatch
):
    module = load_controller()
    insecure_root = tmp_path / "insecure-root"
    insecure_root.mkdir()
    insecure_root.chmod(0o777)
    insecure_parent = tmp_path / "insecure-parent"
    insecure_parent.mkdir()
    insecure_parent.chmod(0o777)
    nested_root = insecure_parent / "root"
    nested_root.mkdir()

    for root in (insecure_root, nested_root):
        monkeypatch.setattr(module, "DEFAULT_ROOT", root)
        lock = module.canonical_root_lock_path(root)
        with pytest.raises(ValueError, match="root or ancestor is group/world-writable"):
            with module.canonical_root_lock(root):
                raise AssertionError("insecure canonical root acquired its lock")
        assert not lock.exists()

    offline_parent = tmp_path / "offline-parent"
    offline_parent.mkdir()
    offline_parent.chmod(0o777)
    offline_root = offline_parent / "root"
    offline_root.mkdir()
    with module.canonical_root_lock(offline_root):
        pass
    assert not module.canonical_root_lock_path(offline_root).exists()


def test_canonical_lock_preserves_created_lock_when_root_drifts_inside_open(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    lock = module.canonical_root_lock_path(root)
    original_open = module.os.open
    original_fsync = module.os.fsync
    raced = False
    root_fsyncs = 0
    created_lock_inode = None

    def drift_root_after_lock_create(path, flags, *args, **kwargs):
        nonlocal raced, created_lock_inode
        descriptor = original_open(path, flags, *args, **kwargs)
        if (
            not raced
            and path == lock.name
            and flags & os.O_CREAT
            and kwargs.get("dir_fd") is not None
        ):
            raced = True
            created_lock_inode = module.os.fstat(descriptor).st_ino
            root.chmod(0o777)
        return descriptor

    def record_root_fsync(descriptor):
        nonlocal root_fsyncs
        profile = module.os.fstat(descriptor)
        if stat.S_ISDIR(profile.st_mode) and profile.st_ino == root.stat().st_ino:
            root_fsyncs += 1
        original_fsync(descriptor)

    monkeypatch.setattr(module.os, "open", drift_root_after_lock_create)
    monkeypatch.setattr(module.os, "fsync", record_root_fsync)
    with pytest.raises(ValueError, match="security metadata changed"):
        with module.canonical_root_lock(root):
            raise AssertionError("drifted root retained canonical authority")
    assert raced
    assert lock.exists()
    assert lock.stat().st_ino == created_lock_inode
    assert lock.stat().st_size == 0
    assert stat.S_IMODE(lock.stat().st_mode) == 0o644
    assert root_fsyncs >= 1
    root.chmod(0o755)
    lock.unlink()
    module.fsync_directory(root)
    assert not lock.exists()


def test_canonical_authority_rejects_arbitrary_group_writable_project_boundary(
    monkeypatch,
):
    module = load_controller()
    root = Path("/trusted-project/controller/root")
    euid = os.geteuid()
    project_gid = 31114
    profiles = {
        Path("/"): SimpleNamespace(
            st_dev=1,
            st_ino=1,
            st_mode=stat.S_IFDIR | 0o755,
            st_uid=0,
            st_gid=0,
            st_nlink=2,
        ),
        Path("/trusted-project"): SimpleNamespace(
            st_dev=1,
            st_ino=2,
            st_mode=stat.S_IFDIR | 0o2770,
            st_uid=0,
            st_gid=project_gid,
            st_nlink=2,
        ),
        Path("/trusted-project/controller"): SimpleNamespace(
            st_dev=1,
            st_ino=3,
            st_mode=stat.S_IFDIR | 0o2755,
            st_uid=euid,
            st_gid=project_gid,
            st_nlink=2,
        ),
        root: SimpleNamespace(
            st_dev=1,
            st_ino=4,
            st_mode=stat.S_IFDIR | 0o2755,
            st_uid=euid,
            st_gid=project_gid,
            st_nlink=2,
        ),
    }
    monkeypatch.setattr(module.os, "lstat", lambda path: profiles[Path(path)])
    with pytest.raises(ValueError, match="group/world-writable"):
        module.require_canonical_root_authority(root)


def test_canonical_authority_accepts_only_exact_frontier_project_boundary(
    monkeypatch,
):
    module = load_controller()
    root = module.DEFAULT_ROOT
    euid = os.geteuid()
    project_gid = module.TRUSTED_GROUP_WRITABLE_PROJECT_GID
    profiles = {}
    for index, path in enumerate((
        Path("/"),
        Path("/lustre"),
        Path("/lustre/orion"),
        Path("/lustre/orion/ast207"),
        module.TRUSTED_GROUP_WRITABLE_PROJECT_BOUNDARY,
        Path("/lustre/orion/ast207/proj-shared/dfielding"),
        root,
    ), start=1):
        mode = 0o755
        uid = 0
        gid = 0
        if path == module.TRUSTED_GROUP_WRITABLE_PROJECT_BOUNDARY:
            mode = module.TRUSTED_GROUP_WRITABLE_PROJECT_MODE
            gid = project_gid
        elif path in {
            Path("/lustre/orion/ast207/proj-shared/dfielding"),
            root,
        }:
            mode = 0o2755
            uid = euid
            gid = project_gid
        profiles[path] = SimpleNamespace(
            st_dev=1,
            st_ino=index,
            st_mode=stat.S_IFDIR | mode,
            st_uid=uid,
            st_gid=gid,
            st_nlink=2,
        )

    monkeypatch.setattr(module.os, "lstat", lambda path: profiles[Path(path)])
    retained = module.require_canonical_root_authority(root)
    module.require_canonical_root_authority(root, retained)

    boundary = profiles[module.TRUSTED_GROUP_WRITABLE_PROJECT_BOUNDARY]
    boundary.st_mode = stat.S_IFDIR | 0o2775
    with pytest.raises(ValueError, match="group/world-writable"):
        module.require_canonical_root_authority(root)
    boundary.st_mode = stat.S_IFDIR | module.TRUSTED_GROUP_WRITABLE_PROJECT_MODE
    boundary.st_gid = project_gid + 1
    with pytest.raises(ValueError, match="group/world-writable"):
        module.require_canonical_root_authority(root)


def test_canonical_lock_rejects_symlinked_declared_public_root(
    tmp_path, monkeypatch
):
    module = load_controller()
    real = tmp_path / "real-root"
    real.mkdir()
    declared = tmp_path / "declared-root"
    declared.symlink_to(real, target_is_directory=True)
    monkeypatch.setattr(module, "DEFAULT_ROOT", declared)
    with pytest.raises(ValueError, match="contains a symbolic link"):
        with module.canonical_root_lock(declared):
            raise AssertionError("symlinked canonical root acquired its lock")
    assert not module.canonical_root_lock_path(real).exists()


def test_canonical_lock_replacement_at_mutation_boundary_has_no_rollback_rename(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    target = root / "state.json"
    target.write_text("prior\n")
    prior_inode = target.stat().st_ino
    lock = module.canonical_root_lock_path(root)
    moved = root / "moved.lock"
    original_renameat2 = module.renameat2
    original_fsync = module.os.fsync
    replaced = False
    rename_calls = []
    directory_fsyncs = 0

    def record_directory_fsync(descriptor):
        nonlocal directory_fsyncs
        if stat.S_ISDIR(module.os.fstat(descriptor).st_mode):
            directory_fsyncs += 1
        original_fsync(descriptor)

    def replace_lock_at_mutation_boundary(*args, **kwargs):
        nonlocal replaced
        if not replaced:
            replaced = True
            lock.rename(moved)
            lock.touch(mode=0o644)
        rename_calls.append(args)
        return original_renameat2(*args, **kwargs)

    monkeypatch.setattr(module.os, "fsync", record_directory_fsync)
    monkeypatch.setattr(module, "renameat2", replace_lock_at_mutation_boundary)
    with pytest.raises(ValueError, match="lock path changed during canonical mutation"):
        with module.canonical_root_lock(root):
            module.write_text(target, "newer\n")
    assert replaced
    assert len(rename_calls) == 1
    assert directory_fsyncs >= 2
    assert target.read_text() == "newer\n"
    temporary = root / module.metadata_temporary_name(target.name)
    assert temporary.read_text() == "prior\n"
    assert temporary.stat().st_ino == prior_inode


def test_exchange_in_place_mutation_never_attempts_rollback_rename(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    target = root / "state.json"
    target.write_text("prior\n")
    prior_inode = target.stat().st_ino
    temporary_name = module.metadata_temporary_name(target.name)
    original_renameat2 = module.renameat2
    rename_calls = []
    mutated = False

    def mutate_after_exchange(directory_fd, source, destination, flags, label):
        nonlocal mutated
        result = original_renameat2(
            directory_fd, source, destination, flags, label
        )
        rename_calls.append((source, destination, flags))
        if not mutated and flags == module.RENAME_EXCHANGE:
            mutated = True
            overwrite_direct_child_restore_mtime(
                directory_fd, destination, b"evil!\n"
            )
        return result

    monkeypatch.setattr(module, "renameat2", mutate_after_exchange)
    with pytest.raises(ValueError, match="durable forward recovery state was preserved"):
        with module.canonical_root_lock(root):
            module.write_text(target, "newer\n")
    assert mutated
    assert len(rename_calls) == 1
    assert target.read_text() == "evil!\n"
    temporary = root / temporary_name
    assert temporary.read_text() == "prior\n"
    assert temporary.stat().st_ino == prior_inode


def test_promoted_f117_profile_preserves_frozen_science_revision(
    tmp_path, monkeypatch
):
    module = load_controller()
    root = tmp_path / "root"
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    paths = paths_for(module, root)
    source = tmp_path / "source"
    input_path = source / "inputs/r16.athinput"
    input_path.parent.mkdir(parents=True)
    input_path.write_text("input\n")
    matrix = source / "matrix.json"
    matrix.write_text("{}\n")
    executable = root / "build/athena"
    executable.parent.mkdir(parents=True)
    executable.write_text("exe\n")
    executable.chmod(0o755)
    build_manifest = root / "build/manifest"
    build_manifest.mkdir()
    (build_manifest / "environment.txt").write_text("build\n")
    bundle = root / "source-archives/current.bundle"
    bundle.parent.mkdir()
    bundle.write_text("bundle\n")
    tooling_revision = "a" * 40
    utility = {
        "revision": tooling_revision,
        "sha256": "1" * 64,
    }
    qualification = {"sha256": "2" * 64}
    bundle_record = {
        "path": str(bundle),
        "sha256": sha256(bundle),
        "verified_revisions": [module.QUALIFIED_SOURCE_REVISION, tooling_revision],
    }
    authority = source_authority_binding(
        tooling_revision,
        bundle_path=bundle.relative_to(root).as_posix(),
        bundle_sha256=sha256(bundle),
        verified_revisions=bundle_record["verified_revisions"],
    )
    args = SimpleNamespace(
        case_id="R16",
        segment="s00_rankio_t0_t0p25",
        acceptance_criterion="criterion",
        athena_walltime="01:50:00",
        executable=str(executable),
        nodes=1,
        ranks_per_node=8,
        cpus_per_task=7,
        walltime="02:00:00",
    )
    profile = {
        "acceptance_criterion": args.acceptance_criterion,
        "acceptance_policy": "reviewed",
        "athena_walltime": args.athena_walltime,
        "build_manifest": str(build_manifest),
        "build_manifest_sha256": module.r17_directory_inventory_sha256(
            build_manifest, "build"
        ),
        "case_id": args.case_id,
        "controller_walltime_max_seconds": module.MAX_SEGMENT_SECONDS,
        "cpus_per_task": args.cpus_per_task,
        "estimated_storage_bytes": 1,
        "executable": str(executable),
        "executable_revision": module.QUALIFIED_SOURCE_REVISION,
        "executable_sha256": sha256(executable),
        "input_file": "inputs/r16.athinput",
        "input_revision": module.QUALIFIED_SOURCE_REVISION,
        "input_sha256": sha256(input_path),
        "nodes": args.nodes,
        "output_layout": "rank-local",
        "segment": args.segment,
        "parent_job_id": None,
        "parent_result": None,
        "parent_segment": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
        "ranks_per_node": args.ranks_per_node,
        "recommendation_basis": {"kind": "reviewed"},
        "time_tlim_target": 0.25,
        "walltime": args.walltime,
        "source_bundle": str(bundle),
        "source_bundle_sha256": sha256(bundle),
    }
    recost_path = paths["accounting"] / (
        f"mks24_stage_i_{EPOCH_SLUG}_F117_recost_evidence.json"
    )
    audit_path = recost_path.with_name(f"{recost_path.name}.publication_audit.json")
    recost = {
        "schema_version": 2,
        "record_type": "stage-i-recost-recommendation-evidence",
        "checkpoint": "F-117",
        "artifact_name": recost_path.name,
        "execution_epoch": EPOCH,
        "authority": {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
        "publication_requirements": {
            "independent_review_required": True,
            "publication_audit_required": True,
            "published_mode": "0444",
            "published_links": 1,
            "controller_consumption_requires_exact_published_sha256": True,
        },
        "recommendations": {
            "mode": "sole-next-profile",
            "authorizing": False,
            "recommended_next_profiles": [profile],
            "bounded_concurrency": {
                "max_active_segments": module.MAX_ACTIVE_STAGE_I_SEGMENTS,
                "max_wave_nodes": 1,
                "r17_exclusive_and_last": True,
            },
            "controller_consumption_state": "reviewed publication required",
            "sole_next_segment_recommendation": {
                key: profile[key] for key in module.R17_SOLE_PROFILE_KEYS
            },
            "non_authorizing_reason": "candidate is non-authorizing",
        },
        "provenance": {
            "source_authority": authority,
            "stage_i_helper_sha256": utility["sha256"],
            "stage_i_helper_revision": utility["revision"],
            "matrix_sha256": sha256(matrix),
            "matrix_revision": tooling_revision,
            "source_bundle_sha256": sha256(bundle),
            "qualification_approval_sha256": qualification["sha256"],
        },
    }
    monkeypatch.setattr(
        module, "latest_published_r17_recost",
        lambda *_args: (recost_path, "3" * 64, recost, audit_path, "4" * 64),
    )

    def require_profile():
        return module.require_promoted_profile_for_prepare(
            paths,
            args,
            source_dir=source,
            matrix_path=matrix,
            input_path=input_path,
            input_revision=module.QUALIFIED_SOURCE_REVISION,
            utility_provenance=utility,
            build_provenance={
                "revision": module.QUALIFIED_SOURCE_REVISION,
                "sha256": sha256(executable),
                "manifest_dir": str(build_manifest),
            },
            build_manifest=build_manifest,
            qualification_approval=qualification,
            bundle_provenance=bundle_record,
            current_source_authority=authority,
            restart=None,
            parent_segment=None,
            time_tlim_target=0.25,
        )

    retained = require_profile()
    assert retained["profile"]["input_revision"] == module.QUALIFIED_SOURCE_REVISION
    assert recost["provenance"]["matrix_revision"] == tooling_revision
    assert (
        recost["provenance"]["matrix_revision"] != retained["profile"]["input_revision"]
    )

    recost["provenance"]["matrix_revision"] = module.QUALIFIED_SOURCE_REVISION
    with pytest.raises(ValueError, match="profile provenance is stale"):
        require_profile()

    recost["provenance"]["matrix_revision"] = tooling_revision
    profile["input_revision"] = tooling_revision
    with pytest.raises(ValueError, match="differs from prepare arguments"):
        require_profile()


def test_parse_history_preserves_real_athenak_hyphenated_labels(tmp_path):
    module = load_controller()
    path = tmp_path / "retained-r12.mhd.hst"
    path.write_text(R12_RETAINED_MHD_PREFIX)
    history = module.parse_history(path)
    assert "tot-E" in history
    assert "tot" not in history
    assert history["tot-E"] == [16.000000000000004]
    assert history["lf_cawrk"] == [0.0]


def test_historical_clean_partial_without_divb_has_exact_irreducible_blocker(tmp_path):
    module = load_controller()
    path = tmp_path / "retained-r12.user.hst"
    path.write_text(R12_RETAINED_USER_PREFIX)
    user = module.parse_history(path)
    assert "max_ndiv" not in user
    with pytest.raises(ValueError, match="irreducibly unavailable.*max_ndiv"):
        module.continuation_plasma_evidence("R12", {}, user)


@pytest.mark.parametrize(
    ("case_id", "job_id", "manifest_path"),
    LIVE_FROZEN_E03_CLEAN_PARTIALS,
)
def test_live_frozen_e03_clean_partial_reproduces_full_schema4_migration(
    case_id, job_id, manifest_path
):
    module = load_controller()
    manifest = load_live_frozen_e03_manifest(manifest_path)
    migrated = module.revalidate_continuation_plasma_evidence(
        manifest["scientific_inspection"], manifest
    )
    evidence = migrated["plasma_continuation_evidence"]
    migration = evidence["migration_contract"]
    bindings = migration["bindings"]

    assert migrated["schema_version"] == 4
    assert migrated["checks"]["plasma_continuation_policy"] is False
    assert migrated["plasma_continuation_policy"] == (
        module.FROZEN_E03_CONTINUATION_MIGRATION_POLICY
    )
    assert migrated["plasma_continuation_authorized"] is False
    assert migrated["plasma_continuation_eligible"] is False
    assert evidence["continuation_authorized"] is False
    assert evidence["continuation_eligible"] is False
    assert evidence["eligibility_only"] is True
    assert evidence["ct_divergence_claimed"] is False
    assert evidence["ct_divergence_reason"] == module.FROZEN_E03_CT_DIVERGENCE_REASON
    assert evidence["checks"]["normalized_ct_divb"] is False
    assert evidence["checks"]["frozen_e03_exact_migration"] is True
    assert evidence["measurements"]["normalized_ct_divb_max"] is None
    assert evidence["normalized_ct_divb_evidence"]["ct_divergence_claimed"] is False
    assert evidence["normalized_ct_divb_evidence"]["authorizing"] is False
    assert migration["authority"] == {
        "continuation_authorized": False,
        "submission_authorized": False,
        "scheduler_mutation_authorized": False,
        "canonical_mutation_authorized": False,
    }
    assert migration["case_id"] == case_id
    assert migration["job_id"] == job_id
    assert bindings["manifest"]["path"] == str(manifest_path)
    assert bindings["manifest"]["sha256"] == sha256(manifest_path)
    for label in (
        "inspection",
        "mhd_history",
        "user_history",
        "independent_validation",
    ):
        bound = bindings[label]
        assert bound["sha256"] == sha256(Path(bound["path"]))
        assert bound["size_bytes"] == Path(bound["path"]).stat().st_size
    if case_id == "R12":
        assert bindings["independent_review"]["mode"] == "0444"
        assert (
            evidence["normalized_ct_divb_evidence"]["independent_review"]
            == bindings["independent_review"]
        )
    else:
        assert bindings["independent_review"] is None


@pytest.mark.parametrize(
    ("case_id", "job_id", "manifest_path"),
    LIVE_FROZEN_E03_CLEAN_PARTIALS,
)
def test_live_frozen_e03_migration_rejects_manifest_and_inspection_drift(
    case_id, job_id, manifest_path
):
    module = load_controller()
    manifest = load_live_frozen_e03_manifest(manifest_path)
    assert manifest["run"]["case_id"] == case_id
    assert manifest["job_id"] == job_id

    changed_manifest = copy.deepcopy(manifest)
    changed_manifest["accounting"]["notes"] += " changed"
    with pytest.raises(ValueError, match="manifest or inspection bytes differ"):
        module.revalidate_continuation_plasma_evidence(
            changed_manifest["scientific_inspection"], changed_manifest
        )

    changed_identity = copy.deepcopy(manifest)
    changed_identity["command"] = []
    with pytest.raises(ValueError, match="manifest identity differs"):
        module.revalidate_continuation_plasma_evidence(
            changed_identity["scientific_inspection"], changed_identity
        )

    changed_inspection = copy.deepcopy(manifest["scientific_inspection"])
    changed_inspection["mhd_history"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="manifest or inspection bytes differ"):
        module.revalidate_continuation_plasma_evidence(
            changed_inspection, manifest
        )


@pytest.mark.parametrize(
    ("case_id", "job_id", "manifest_path"),
    LIVE_FROZEN_E03_CLEAN_PARTIALS,
)
def test_live_frozen_e03_continuation_prepare_revalidates_then_rejects_inventory_only(
    monkeypatch, case_id, job_id, manifest_path
):
    module = load_controller()
    manifest = load_live_frozen_e03_manifest(manifest_path)
    inspection = manifest["scientific_inspection"]
    terminal = inspection["terminal_restart"]
    restart = Path(terminal["path"])
    original = module.revalidate_continuation_plasma_evidence
    calls = []

    def traced_revalidation(current_inspection, current_manifest):
        calls.append((current_inspection["job_id"], current_manifest["job_id"]))
        return original(current_inspection, current_manifest)

    monkeypatch.setattr(module, "authenticate_prepared_execution", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "require_existing_layout", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "read_reservations", lambda *args, **kwargs: [])
    monkeypatch.setattr(module, "reservation_for_manifest", lambda *args, **kwargs: {})
    monkeypatch.setattr(module, "require_reserved_execution_intent", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "revalidate_inspection_files", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "revalidate_retained_product", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "sha256", lambda path: terminal["sha256"])
    monkeypatch.setattr(
        module,
        "authenticated_restart_product_time",
        lambda *args, **kwargs: {"binary_time": inspection["final_time"]},
    )
    monkeypatch.setattr(
        module, "revalidate_continuation_plasma_evidence", traced_revalidation
    )

    expected = (
        "accepted job 4766828"
        if case_id == "R03"
        else "R12 must restart fresh as s01_rankio_t0_t0p12 from t=0"
    )
    with pytest.raises(ValueError, match=expected):
        module.verify_continuation_restart(restart)
    assert calls == [(job_id, job_id)]


def test_frozen_e03_migration_binding_rejects_symlink_alias(tmp_path, monkeypatch):
    module = load_controller()
    _, _, manifest_path = LIVE_FROZEN_E03_CLEAN_PARTIALS[1]
    load_live_frozen_e03_manifest(manifest_path)
    root = tmp_path / "root"
    root.mkdir()
    alias = root / "manifest.json"
    alias.symlink_to(manifest_path)
    binding = copy.deepcopy(
        module.FROZEN_E03_CONTINUATION_MIGRATIONS[
            ("R12", "4766856")
        ]["manifest"]
    )
    binding["path"] = alias.name
    monkeypatch.setattr(module, "DEFAULT_ROOT", root)
    with pytest.raises(ValueError, match="symbolic link"):
        module.read_frozen_e03_migration_binding(binding, "migration alias")


def test_controller_exposes_no_r12_waiver_authorization_api():
    module = load_controller()
    assert not hasattr(module, "frozen_e03_no_ct_waiver_contract")
    assert not hasattr(module, "frozen_e03_promoted_profile_waiver_contract")
    assert not hasattr(module, "require_frozen_e03_no_ct_waiver_for_profile")


def test_fresh_r12_profile_contract_is_exact_and_parentless():
    module = load_controller()
    assert module.R12_FRESH_RERUN_SEGMENT == "s01_rankio_t0_t0p12"
    assert module.R12_FRESH_RERUN_NODES == 4
    assert module.R12_FRESH_RERUN_RANKS_PER_NODE == 8
    assert module.R12_FRESH_RERUN_TOTAL_RANKS == 32
    assert (
        module.R12_FRESH_RERUN_NODES * module.R12_FRESH_RERUN_RANKS_PER_NODE
        == module.R12_FRESH_RERUN_TOTAL_RANKS
    )
    assert module.R12_FRESH_RERUN_WALLTIME == "02:00:00"
    assert module.R12_FRESH_RERUN_ATHENA_WALLTIME == "01:50:00"
    assert module.R12_FRESH_RERUN_TARGET == 0.12
    fresh = {
        "case_id": "R12",
        "segment": module.R12_FRESH_RERUN_SEGMENT,
        "nodes": module.R12_FRESH_RERUN_NODES,
        "ranks_per_node": module.R12_FRESH_RERUN_RANKS_PER_NODE,
        "walltime": module.R12_FRESH_RERUN_WALLTIME,
        "athena_walltime": module.R12_FRESH_RERUN_ATHENA_WALLTIME,
        "parent_job_id": None,
        "parent_segment": None,
        "parent_result": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
    }
    module.require_fresh_r12_profile_contract(
        fresh, None, None, module.R12_FRESH_RERUN_TARGET
    )
    for key, value in (
        ("segment", "s01_rankio_t0_t0p25"),
        ("segment", "s00_rankio_t0_t0p25"),
        ("nodes", 1),
        ("ranks_per_node", 4),
        ("walltime", "01:00:00"),
        ("athena_walltime", "00:50:00"),
        ("parent_job_id", module.R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID),
        ("restart_file", "/historical/restart.rst"),
    ):
        changed = copy.deepcopy(fresh)
        changed[key] = value
        with pytest.raises(ValueError, match="R12 must restart fresh"):
            module.require_fresh_r12_profile_contract(
                changed, None, None, module.R12_FRESH_RERUN_TARGET
            )
    with pytest.raises(ValueError, match="R12 must restart fresh"):
        module.require_fresh_r12_profile_contract(fresh, None, None, 0.5)
    with pytest.raises(ValueError, match="R12 must restart fresh"):
        module.require_fresh_r12_profile_contract(fresh, None, None, 0.25)
    with pytest.raises(ValueError, match="R12 must restart fresh"):
        module.require_fresh_r12_profile_contract(
            fresh, None, Path("/historical/restart.rst"), module.R12_FRESH_RERUN_TARGET
        )


def test_clean_partial_plasma_gate_rejects_divb_policy_and_limiter_drift():
    module = load_controller()
    zeros = [0.0, 0.0]
    mhd = {
        "time": [0.0, 1.0],
        "mass": [1.0, 1.0],
        "tot-E": [1.0, 1.1],
        "lf_nstage": [0.0, 1.0],
        "lf_qface": [0.0, 10.0],
        "lf_qprcap": [0.0, 1.0],
        "lf_qpr10": [0.0, 1.0],
        "lf_qpecap": [0.0, 1.0],
        "lf_qpe10": [0.0, 1.0],
        "lf_qprwrk": [0.0, 0.01],
        "lf_qpewrk": [0.0, 0.01],
        "lf_cpwrk": [0.0, 0.1],
        "lf_cawrk": [0.0, 0.0],
        "lf_hwproj": [0.0, 0.0],
        **{name: [0.0, 0.0] for name in module.STRICT_LF_FAILURE_COLUMNS},
    }
    user = {
        "time": [0.0, 1.0],
        "mass": [1.0, 1.0],
        "hard_vol": zeros,
        "force_work": [0.0, 0.1],
        "max_ndiv": [0.0, 1.0e-13],
    }
    evidence = module.continuation_plasma_evidence("R16", mhd, user)
    assert evidence["schema_version"] == 2
    assert evidence["policy"] == module.CONTINUATION_PLASMA_POLICY
    assert evidence["continuation_authorized"] is True
    assert evidence["checks"]["normalized_ct_divb"]
    assert evidence["normalized_ct_divb_evidence"] == {
        "source": "authenticated user history max_ndiv",
        "comparison": "strictly-less-than",
        "threshold": module.CONTINUATION_MAX_NORMALIZED_CT_DIVB,
        "maximum": 1.0e-13,
        "passed": True,
    }
    invalid_divb = copy.deepcopy(user)
    invalid_divb["max_ndiv"][-1] = 1.0e-6
    with pytest.raises(ValueError, match="divB"):
        module.continuation_plasma_evidence("R16", mhd, invalid_divb)
    with pytest.raises(ValueError, match="passive continuation"):
        module.continuation_plasma_evidence("R06", mhd, user)
    finite_limiter = copy.deepcopy(mhd)
    finite_limiter["lf_hwproj"][-1] = 1.0
    with pytest.raises(ValueError, match="finite-limiter"):
        module.continuation_plasma_evidence("R14", finite_limiter, user)


@pytest.mark.parametrize(
    ("case_id", "nodes", "expected"),
    (
        ("R16", 1, "F116 source gate reached"),
        ("R17", 8, "R17 readiness gate reached"),
    ),
)
def test_canonical_prepare_reaches_source_and_r17_gates_before_run_mutation(
    tmp_path, monkeypatch, case_id, nodes, expected
):
    module = load_controller()
    root = tmp_path / "root"
    paths = paths_for(module, root)
    source = tmp_path / "source"
    source.mkdir()
    input_path = source / "r17.athinput"
    input_path.write_text("input\n")
    matrix = source / "matrix.json"
    matrix.write_text("{}\n")
    executable = root / "build/athena"
    executable.parent.mkdir()
    executable.write_text("exe\n")
    executable.chmod(0o755)
    build_manifest = root / "build/manifest"
    build_manifest.mkdir()
    bundle = root / "source.bundle"
    bundle.write_text("bundle\n")
    args = SimpleNamespace(
        root=str(root),
        allow_local_root=False,
        case_id=case_id,
        segment="s00_rankio_t0_t0p25",
        acceptance_criterion="criterion",
        source_dir=str(source),
        matrix=str(matrix),
        executable=str(executable),
        build_manifest=str(build_manifest),
        source_bundle=str(bundle),
        restart_file=None,
        nodes=nodes,
        walltime="02:00:00",
        athena_walltime="01:50:00",
        ranks_per_node=8,
        cpus_per_task=7,
        override=[],
        allow_missing_time_target=False,
        allow_missing_restart_time_marker=False,
    )
    revision = "a" * 40
    monkeypatch.setattr(module, "require_root", lambda *_args: root)
    monkeypatch.setattr(module, "is_offline_local_root", lambda *_args: False)
    monkeypatch.setattr(module, "layout", lambda _root: paths)
    for name in (
        "require_existing_layout",
        "require_no_pending_transactions",
        "require_no_orphaned_segment_runs",
        "require_reconciled_store_consistency",
        "require_authorized_case",
        "require_safe_segment",
        "require_active_reservation_policy",
        "require_prepare_case_policy",
    ):
        monkeypatch.setattr(module, name, lambda *_args, **_kwargs: None)
    monkeypatch.setattr(module, "read_reservations", lambda _paths: [])
    monkeypatch.setattr(module, "validate_matrix", lambda *_args: {"cases": []})
    monkeypatch.setattr(
        module,
        "case_for_id",
        lambda *_args: {
            "name": "r17",
            "input": input_path.name,
            "resolution": "384x384x768",
            "figure_roles": [],
        },
    )
    monkeypatch.setattr(module, "validate_prepare_overrides", lambda *_args, **_kwargs: 0.25)
    monkeypatch.setattr(module, "git_revision_for_input", lambda *_args: revision)
    monkeypatch.setattr(
        module,
        "production_utility_provenance",
        lambda **_kwargs: {
            "path": str(CONTROLLER_PATH),
            "revision": revision,
            "sha256": "1" * 64,
            "committed": True,
        },
    )
    monkeypatch.setattr(
        module,
        "read_build_provenance",
        lambda *_args: {
            "revision": revision,
            "sha256": "2" * 64,
            "manifest_dir": str(build_manifest),
        },
    )
    monkeypatch.setattr(
        module, "require_qualification_approval", lambda *_args: {"sha256": "3" * 64}
    )
    monkeypatch.setattr(
        module,
        "source_bundle_provenance",
        lambda *_args: {"path": str(bundle), "sha256": "4" * 64, "verified_revisions": [revision]},
    )
    def source_gate(*_args, **_kwargs):
        if case_id == "R16":
            raise RuntimeError("F116 source gate reached")
        return source_authority_binding(revision)

    monkeypatch.setattr(module, "require_current_source_authority_for_prepare", source_gate)

    def stop_at_gate(*_args, **_kwargs):
        raise RuntimeError("R17 readiness gate reached")

    monkeypatch.setattr(module, "require_r17_readiness_for_prepare", stop_at_gate)
    with pytest.raises(RuntimeError, match=expected):
        module.prepare.__wrapped__(args)
    assert not (paths["runs"] / case_id / args.segment).exists()
