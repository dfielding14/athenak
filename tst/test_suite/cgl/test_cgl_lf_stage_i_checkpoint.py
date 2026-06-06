"""Local-fixture regressions for retained Stage I recost publication."""

from __future__ import annotations

import ast
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import errno
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
STAGE_I = REPOSITORY / "scripts/frontier/cgl_lf_stage_i.py"
CHECKPOINT = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_checkpoint.py"
RECOST = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_recost.py"
RECOST_TEST = REPOSITORY / "tst/test_suite/cgl/test_cgl_lf_stage_i_recost.py"
PRODUCTION_ROOT = "/lustre/orion/ast207/proj-shared/dfielding/CGL"
ARTIFACT = "mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json"
V2_ARTIFACT = "mks24_stage_i_E03_forcing_policy_F114_recost_evidence.json"
EPOCH = "E03-forcing-policy"
EPOCH_SLUG = "E03_forcing_policy"


def sha256(path: Path) -> str:
    """Return one retained fixture digest."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    """Write stable fixture JSON."""

    if path.exists():
        path.chmod(0o644)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def load_checkpoint_module():
    """Load the retained companion without executing its CLI."""

    spec = importlib.util.spec_from_file_location("cgl_lf_stage_i_checkpoint", CHECKPOINT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def force_checkpoint_renameat2_einval(path: Path) -> None:
    """Make one fixture checkpoint model Orion renameat2 flag rejection."""

    retained = path.read_text()
    marker = (
        '    """Perform one descriptor-relative Linux renameat2 operation."""\n\n'
        "    source = require_entry_name(source, label)\n"
    )
    replacement = (
        '    """Perform one descriptor-relative Linux renameat2 operation."""\n\n'
        "    raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))\n"
    )
    assert retained.count(marker) == 1
    path.write_text(retained.replace(marker, replacement))
    path.chmod(0o755)


def leave_lustre_json_post_link_state(module, target: Path, value: object,
                                      monkeypatch, *, mode: int = 0o644) -> Path:
    """Interrupt one forced-EINVAL JSON publication after link and before unlink."""

    real_unlink = module.os.unlink
    interrupted = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def interrupt_temporary_unlink(name, *args, **kwargs):
        nonlocal interrupted
        directory_descriptor = kwargs.get("dir_fd")
        if (
            not interrupted
            and isinstance(name, str)
            and directory_descriptor is not None
            and module.json_temporary_target_name(name) == target.name
        ):
            temporary_profile = os.stat(
                name, dir_fd=directory_descriptor, follow_symlinks=False
            )
            public_profile = os.stat(
                target.name, dir_fd=directory_descriptor, follow_symlinks=False
            )
            if (
                module.profile_identity(temporary_profile)
                == module.profile_identity(public_profile)
                and temporary_profile.st_nlink == public_profile.st_nlink == 2
            ):
                interrupted = True
                raise OSError(errno.EIO, os.strerror(errno.EIO))
        return real_unlink(name, *args, **kwargs)

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    monkeypatch.setattr(module.os, "unlink", interrupt_temporary_unlink)
    with pytest.raises(OSError, match=os.strerror(errno.EIO)):
        module.write_json(target, value, mode=mode)
    monkeypatch.setattr(module.os, "unlink", real_unlink)
    assert interrupted
    temporaries = module.json_temporary_entries(target)
    assert len(temporaries) == 1
    return temporaries[0]


def load_current_recost_test_module():
    """Load the current generator's fixture builders for an end-to-end probe."""

    name = f"_cgl_lf_stage_i_current_recost_test_{os.getpid()}_{id(object())}"
    spec = importlib.util.spec_from_file_location(name, RECOST_TEST)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def recost_fixture(tmp_path):
    """Create one staged-only local recost boundary."""

    repository = tmp_path / "repository"
    frontier = repository / "scripts/frontier"
    frontier.mkdir(parents=True)
    checkpoint = frontier / CHECKPOINT.name
    checkpoint.write_bytes(CHECKPOINT.read_bytes())
    checkpoint.chmod(0o755)
    stage_i = frontier / STAGE_I.name
    stage_i.write_text(
        STAGE_I.read_text().replace(
            PRODUCTION_ROOT,
            str(tmp_path / "local-canonical-sentinel"),
        )
    )
    stage_i.chmod(0o644)
    subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
    subprocess.run(["git", "add", "."], cwd=repository, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=CGL fixture",
            "-c",
            "user.email=cgl-fixture@example.invalid",
            "commit",
            "-q",
            "-m",
            "Create local checkpoint fixture",
        ],
        cwd=repository,
        check=True,
    )
    root = tmp_path / "root"
    subprocess.run(
        [
            "/usr/bin/python3.11",
            "-I",
            "-S",
            "-B",
            str(stage_i),
            "--root",
            str(root),
            "--allow-local-root",
            "init",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    accounting = root / "accounting"
    generator = accounting / "utilities/generate_test_recost.py"
    generator.parent.mkdir()
    generator.write_text("#!/usr/bin/env python3\n# retained test generator\n")
    generator.chmod(0o755)
    scheduler = accounting / "12345.stage_i.sacct.txt"
    scheduler.write_text("12345|COMPLETED|0:0|1|2340\n")
    scheduler.chmod(0o644)
    queue = tmp_path / "squeue.txt"
    queue.write_text("")
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    source_bundle = root / "source-archives/athenak-legacy-fixture.bundle"
    source_bundle.parent.mkdir()
    subprocess.run(
        ["git", "bundle", "create", str(source_bundle), "--all"],
        cwd=repository,
        check=True,
    )
    source_bundle.chmod(0o644)
    profile = {
        "athena_timeout": "00:55:00",
        "nodes": 1,
        "segment": "R02/s19_rankio_t7p25_t7p5",
        "slurm_walltime": "01:05:00",
        "source_bundle": str(source_bundle),
        "source_bundle_sha256": sha256(source_bundle),
        "threshold_seconds": 2700,
    }
    counts = {
        "transactions": 0,
        "reservations": 0,
        "active_reservations": 0,
        "ledger_rows": 0,
        "manifests": 0,
    }
    payload = {
        "execution_epoch": EPOCH,
        "authorization": {"sole_next_segment_profile": profile},
        "reconcile": {"counts": counts},
        "provenance": {
            "stage_i_helper_sha256": sha256(stage_i),
            "generator_sha256": sha256(generator),
            "scheduler_sha256": sha256(scheduler),
        },
    }
    staged = accounting / f"{ARTIFACT}.staged"
    write_json(staged, payload)
    staged.chmod(0o644)
    return {
        "root": root,
        "repository": repository,
        "checkpoint": checkpoint,
        "stage_i": stage_i,
        "accounting": accounting,
        "generator": generator,
        "scheduler": scheduler,
        "queue": queue,
        "source_bundle": source_bundle,
        "source_bundle_revisions": [revision],
        "profile": profile,
        "counts": counts,
        "staged": staged,
        "artifact_sha256": sha256(staged),
        "canonical": accounting / ARTIFACT,
        "audit": accounting / f"{ARTIFACT}.publication_audit.json",
        "stage_i_transactions": (
            accounting / f"mks24_stage_i_{EPOCH_SLUG}_transactions"
        ),
        "recost_transactions": (
            accounting / f"mks24_stage_i_{EPOCH_SLUG}_recost_transactions"
        ),
        "recost_forensics": (
            accounting / f"mks24_stage_i_{EPOCH_SLUG}_recost_forensics"
            / ARTIFACT
        ),
        "recost_forensics_root": (
            accounting / f"mks24_stage_i_{EPOCH_SLUG}_recost_forensics"
        ),
        "lock": root / f".mks24_stage_i_{EPOCH_SLUG}.lock",
    }


def fake_v2_generator_source() -> str:
    """Return a small generator-interface double that reauthenticates all inputs."""

    return textwrap.dedent(
        """\
        #!/usr/bin/env python3
        import hashlib
        import json
        from pathlib import Path


        def digest(path):
            return hashlib.sha256(Path(path).read_bytes()).hexdigest()


        class Tracker:
            def __init__(self, bindings):
                self.bindings = bindings

            def reauthenticate_all(self):
                for path, expected in self.bindings:
                    if digest(path) != expected:
                        raise ValueError(f"tracked V2 input checksum changed: {path}")


        class Build:
            def __init__(self, payload, tracker, reconcile, helper_path, helper_sha256,
                         helper_revision, matrix_path, matrix_sha256, matrix_revision,
                         storage):
                self.payload = payload
                self.tracker = tracker
                self.reconcile = reconcile
                self.helper_path = helper_path
                self.helper_sha256 = helper_sha256
                self.helper_revision = helper_revision
                self.matrix_path = matrix_path
                self.matrix_sha256 = matrix_sha256
                self.matrix_revision = matrix_revision
                self.storage_available_bytes = storage["available_bytes"]
                self.storage_retained_stage_i_bytes = storage["retained_stage_i_bytes"]
                self.storage_required_safety_bytes = storage["required_safety_bytes"]
                self.projected_storage_bytes = storage[
                    "projected_authorized_wave_growth_bytes"
                ]
                self.storage_measurements = ("fixture-directory-measurement",)


        def build_payload(args, root, source_path, repository, generator_sha256, *,
                          stage_i_lock_held=False):
            request = json.loads(args.request.read_text())
            inputs = request["inputs"]
            bindings = [
                (args.request, args.expected_request_sha256),
                (source_path, generator_sha256),
            ]
            helper_path = repository / "scripts/frontier/cgl_lf_stage_i.py"
            helper = inputs["stage_i_helper"]
            bindings.append((helper_path, helper["sha256"]))
            matrix = inputs["matrix"]
            matrix_path = repository / matrix["path"]
            bindings.append((matrix_path, matrix["sha256"]))
            for key in (
                "reconciliation", "ledger", "reservations", "storage_evidence",
                "ceiling_evidence", "ceiling_publication_audit",
                "predecessor_recost", "predecessor_recost_publication_audit",
            ):
                item = inputs[key]
                bindings.append((root / item["path"], item["sha256"]))
            bundle = inputs["source_bundle"]
            bindings.append((root / bundle["path"], bundle["sha256"]))
            for key in ("manifests", "scheduler_evidence"):
                for item in inputs[key]:
                    bindings.append((root / item["path"], item["sha256"]))
            readiness = inputs["r17_readiness_evidence"]
            if readiness is not None:
                bindings.append((root / readiness["path"], readiness["sha256"]))
            artifact = root / "accounting" / f"{request['artifact_name']}.staged"
            if not artifact.exists():
                artifact = root / "accounting" / request["artifact_name"]
            reconcile = json.loads((root / inputs["reconciliation"]["path"]).read_text())
            return Build(
                artifact.read_bytes(),
                Tracker(bindings),
                reconcile,
                helper_path,
                helper["sha256"],
                helper["revision"],
                matrix_path,
                matrix["sha256"],
                matrix["revision"],
                json.loads(artifact.read_text())["storage"],
            )


        def run_authenticated_reconcile(helper_path, helper_sha256, root, *,
                                         stage_i_lock_held=False):
            return json.loads((root / "accounting/v2_reconciliation.json").read_text())


        def require_empty_transaction_stores(root):
            return None


        def require_live_storage_boundary(root, available_bytes, retained_stage_i_bytes,
                                          required_safety_bytes, projected_growth_bytes):
            return None


        def require_directory_measurement_boundaries(measurements):
            if measurements != ("fixture-directory-measurement",):
                raise ValueError("fixture directory measurement boundary differs")
        """
    )


def bounded_profile(repository: Path, source_bundle: Path, revision: str,
                    case_id: str, segment: str, nodes: int) -> dict[str, object]:
    """Return one exact fresh V2 bounded-wave profile."""

    executable = repository / "build/bin/athena"
    build_manifest = repository / "build/build_manifest.json"
    input_file = repository / "inputs/cgl_lf_paper/v2_fixture.in"
    return {
        "acceptance_criterion": f"{case_id} fixture acceptance",
        "acceptance_policy": f"{case_id} fixture policy",
        "athena_walltime": "01:50:00",
        "build_manifest": str(build_manifest),
        "build_manifest_sha256": sha256(build_manifest),
        "case_id": case_id,
        "controller_walltime_max_seconds": 7200,
        "cpus_per_task": 7,
        "estimated_storage_bytes": 1024,
        "executable": str(executable),
        "executable_revision": revision,
        "executable_sha256": sha256(executable),
        "input_file": str(input_file),
        "input_revision": revision,
        "input_sha256": sha256(input_file),
        "nodes": nodes,
        "output_layout": "rank-local",
        "segment": segment,
        "parent_job_id": None,
        "parent_result": None,
        "parent_segment": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
        "ranks_per_node": 8,
        "time_tlim_target": (
            0.12
            if case_id == "R12" and segment == "s01_rankio_t0_t0p12"
            else 0.25
        ),
        "walltime": "02:00:00",
        "source_bundle": str(source_bundle),
        "source_bundle_sha256": sha256(source_bundle),
    }


@pytest.fixture
def bounded_recost_fixture(recost_fixture):
    """Replace the sole-profile fixture with one exact generalized V2 packet."""

    fixture = recost_fixture
    root = fixture["root"]
    repository = fixture["repository"]
    matrix = repository / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
    matrix.parent.mkdir(parents=True)
    write_json(matrix, {"fixture": "V2 matrix"})
    executable = repository / "build/bin/athena"
    executable.parent.mkdir(parents=True)
    executable.write_text("fixture executable\n")
    build_manifest = repository / "build/build_manifest.json"
    write_json(build_manifest, {"fixture": "build"})
    input_file = repository / "inputs/cgl_lf_paper/v2_fixture.in"
    input_file.write_text("<problem>\nfixture = true\n")
    live_recost = fixture["repository"] / "scripts/frontier/cgl_lf_stage_i_recost.py"
    live_recost.write_text(fake_v2_generator_source())
    live_recost.chmod(0o755)
    subprocess.run(["git", "add", "."], cwd=fixture["repository"], check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=CGL fixture",
            "-c",
            "user.email=cgl-fixture@example.invalid",
            "commit",
            "-q",
            "-m",
            "Add bounded recost generator",
        ],
        cwd=fixture["repository"],
        check=True,
    )
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=fixture["repository"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    fixture["generator"].write_bytes(live_recost.read_bytes())
    fixture["generator"].chmod(0o755)
    source_bundle = root / "source-archives/athenak-fixture.bundle"
    source_bundle.parent.mkdir(exist_ok=True)
    subprocess.run(
        ["git", "bundle", "create", str(source_bundle), "--all"],
        cwd=fixture["repository"],
        check=True,
    )
    source_bundle.chmod(0o644)
    fixture["scheduler"].write_text(
        "12345|cgl_mks24_E03_forcing_policy_R02_s28_rankio_t9p75_t10|"
        "COMPLETED|0:0|1|2340|2026-06-04T01:00:00+00:00|"
        "2026-06-04T01:39:00+00:00\n"
    )
    second_scheduler = fixture["accounting"] / "23456.stage_i.sacct.txt"
    second_scheduler.write_text(
        "23456|cgl_mks24_E03_forcing_policy_R03_s00_rankio_t0_t0p5|"
        "COMPLETED|0:0|1|1800|2026-06-04T02:00:00+00:00|"
        "2026-06-04T02:30:00+00:00\n"
    )
    second_scheduler.chmod(0o644)
    r12_scheduler = fixture["accounting"] / "4766856.stage_i.sacct.txt"
    r12_scheduler.write_text(
        "4766856|cgl_mks24_E03_forcing_policy_R12_s00_rankio_t0_t0p25|"
        "COMPLETED|0:0|4|7200|2026-06-04T03:00:00+00:00|"
        "2026-06-04T05:00:00+00:00\n"
    )
    r12_scheduler.chmod(0o644)
    scheduler_evidence = [
        {
            "path": str(fixture["scheduler"].relative_to(root)),
            "sha256": sha256(fixture["scheduler"]),
            "job_id": "12345",
            "job_name": "cgl_mks24_E03_forcing_policy_R02_s28_rankio_t9p75_t10",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "nodes": 1,
            "elapsed_seconds": 2340,
            "submitted_utc": "2026-06-04T01:00:00+00:00",
            "completed_utc": "2026-06-04T01:39:00+00:00",
        },
        {
            "path": str(second_scheduler.relative_to(root)),
            "sha256": sha256(second_scheduler),
            "job_id": "23456",
            "job_name": "cgl_mks24_E03_forcing_policy_R03_s00_rankio_t0_t0p5",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "nodes": 1,
            "elapsed_seconds": 1800,
            "submitted_utc": "2026-06-04T02:00:00+00:00",
            "completed_utc": "2026-06-04T02:30:00+00:00",
        },
        {
            "path": str(r12_scheduler.relative_to(root)),
            "sha256": sha256(r12_scheduler),
            "job_id": "4766856",
            "job_name": "cgl_mks24_E03_forcing_policy_R12_s00_rankio_t0_t0p25",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "nodes": 4,
            "elapsed_seconds": 7200,
            "submitted_utc": "2026-06-04T03:00:00+00:00",
            "completed_utc": "2026-06-04T05:00:00+00:00",
        },
    ]
    profiles = [
        bounded_profile(
            repository, source_bundle, revision, "R04", "s00_rankio_t0_t0p25", 4
        ),
        bounded_profile(
            repository, source_bundle, revision, "R12", "s01_rankio_t0_t0p12", 4
        ),
    ]
    authorization = {
        "mode": "bounded-wave",
        "authorizing": False,
        "authorized_next_profiles": profiles,
        "bounded_concurrency": {
            "max_active_segments": 4,
            "max_wave_nodes": 8,
            "r17_exclusive_and_last": True,
        },
        "controller_consumption_state": (
            "non-authorizing advisory bounded wave pending promoted controller "
            "consumption and generalized checkpoint support"
        ),
        "non_authorizing_reason": (
            "Controller consumption remains a separate promoted transition."
        ),
    }
    reconcile_completed = subprocess.run(
        [
            "/usr/bin/python3.11",
            "-I",
            "-S",
            "-B",
            str(fixture["stage_i"]),
            "--root",
            str(root),
            "--allow-local-root",
            "reconcile",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    reconcile = json.loads(reconcile_completed.stdout)
    reconciliation = fixture["accounting"] / "v2_reconciliation.json"
    write_json(reconciliation, reconcile)
    ledger = fixture["accounting"] / f"mks24_stage_i_{EPOCH_SLUG}_node_hours.csv"
    reservations = fixture["accounting"] / f"mks24_stage_i_{EPOCH_SLUG}_reservations.json"
    reservations.chmod(0o644)
    storage_evidence = fixture["accounting"] / "v2_storage_evidence.json"
    write_json(storage_evidence, {"fixture": "storage"})
    f113 = fixture["accounting"] / (
        "mks24_stage_i_E03_forcing_policy_F113_controller_transition_evidence.json"
    )
    write_json(f113, {"fixture": "F113 ceiling"})
    f113_audit = f113.with_name(f"{f113.name}.publication_audit.json")
    f113_audit_value = {"fixture": "F113 publication audit"}
    write_json(f113_audit, f113_audit_value)
    predecessor = fixture["accounting"] / (
        "mks24_stage_i_E03_forcing_policy_F113_recost_evidence.json"
    )
    generated = datetime.now(timezone.utc).replace(microsecond=0)
    predecessor_generated = generated - timedelta(hours=1)
    predecessor_published = generated - timedelta(minutes=30)
    write_json(
        predecessor,
        {
            "checkpoint": "F-113",
            "execution_epoch": EPOCH,
            "generated_utc": predecessor_generated.isoformat(),
            "record_type": "stage-i-barrier-recost-checkpoint",
        },
    )
    predecessor_audit = predecessor.with_name(
        f"{predecessor.name}.publication_audit.json"
    )
    write_json(
        predecessor_audit,
        {
            "artifact": {
                "links": 1,
                "mode": "0644",
                "path": str(predecessor),
                "sha256": sha256(predecessor),
            },
            "execution_epoch": EPOCH,
            "published_utc": predecessor_published.isoformat(),
            "record_type": "observed-publication",
        },
    )
    inputs = {
        "reconciliation": {
            "path": str(reconciliation.relative_to(root)),
            "sha256": sha256(reconciliation),
        },
        "ledger": {
            "path": str(ledger.relative_to(root)),
            "sha256": sha256(ledger),
        },
        "reservations": {
            "path": str(reservations.relative_to(root)),
            "sha256": sha256(reservations),
        },
        "manifests": [],
        "scheduler_evidence": [
            {"path": item["path"], "sha256": item["sha256"]}
            for item in scheduler_evidence
        ],
        "storage_evidence": {
            "path": str(storage_evidence.relative_to(root)),
            "sha256": sha256(storage_evidence),
        },
        "source_bundle": {
            "path": str(source_bundle.relative_to(root)),
            "sha256": sha256(source_bundle),
            "verified_revisions": [revision],
        },
        "matrix": {
            "path": str(matrix.relative_to(repository)),
            "revision": revision,
            "sha256": sha256(matrix),
        },
        "stage_i_helper": {
            "revision": revision,
            "sha256": sha256(fixture["stage_i"]),
        },
        "ceiling_evidence": {
            "path": str(f113.relative_to(root)),
            "sha256": sha256(f113),
        },
        "ceiling_publication_audit": {
            "path": str(f113_audit.relative_to(root)),
            "sha256": sha256(f113_audit),
        },
        "predecessor_recost": {
            "path": str(predecessor.relative_to(root)),
            "sha256": sha256(predecessor),
        },
        "predecessor_recost_publication_audit": {
            "path": str(predecessor_audit.relative_to(root)),
            "sha256": sha256(predecessor_audit),
        },
        "r17_readiness_evidence": None,
    }
    request = fixture["accounting"] / "v2_recost_request.json"
    request_value = {
        "schema_version": 1,
        "record_type": "stage-i-barrier-recost-request",
        "checkpoint": "F-114",
        "artifact_name": V2_ARTIFACT,
        "execution_epoch": EPOCH,
        "generated_utc": generated.isoformat(),
        "expires_utc": (generated + timedelta(hours=1)).isoformat(),
        "scope": "fixture V2 bounded-wave recost",
        "barrier": {
            "recorded_segments": [
                {
                    "case_id": "R02",
                    "segment": "s28_rankio_t9p75_t10",
                    "job_id": "12345",
                    "result": "accepted",
                },
                {
                    "case_id": "R03",
                    "segment": "s00_rankio_t0_t0p5",
                    "job_id": "23456",
                    "result": "clean_partial",
                },
                {
                    "case_id": "R12",
                    "segment": "s00_rankio_t0_t0p25",
                    "job_id": "4766856",
                    "result": "clean_partial",
                },
            ],
        },
        "inputs": inputs,
        "authorization": {
            "mode": "bounded-wave",
            "max_wave_nodes": 8,
            "profiles": profiles,
        },
    }
    write_json(request, request_value)
    predecessor_binding = {
        "artifact_name": predecessor.name,
        "path": str(predecessor),
        "sha256": sha256(predecessor),
        "publication_audit_path": str(predecessor_audit),
        "publication_audit_sha256": sha256(predecessor_audit),
        "checkpoint": "F-113",
        "generated_utc": predecessor_generated.isoformat(),
        "published_utc": predecessor_published.isoformat(),
    }
    budget = {
        "actual_plus_authorized_wave_node_hours": "12",
        "actual_stage_i_node_hours": "0",
        "authorized_wave_reserved_node_hours": "12",
        "case_breakdown": {},
        "computed_remaining_stage_i_node_hours": "12",
        "computed_stage_i_margin_node_hours": "88",
        "computed_stage_i_total_node_hours": "12",
        "method": "authenticated fixture projection",
        "project_ceiling_node_hours": "100",
        "promoted_stage_i_envelope_node_hours": "100",
    }
    projection_sha256 = hashlib.sha256(
        (json.dumps(budget, sort_keys=True) + "\n").encode()
    ).hexdigest()
    lineage_sha256 = hashlib.sha256(b"authenticated fixture lineages").hexdigest()
    projected_storage = sum(profile["estimated_storage_bytes"] for profile in profiles)
    provenance = {
        "request_sha256": sha256(request),
        "generator_sha256": sha256(fixture["generator"]),
        "generator_revision": revision,
        "stage_i_helper_sha256": sha256(fixture["stage_i"]),
        "stage_i_helper_revision": revision,
        "matrix_sha256": sha256(matrix),
        "matrix_revision": revision,
        "source_bundle_sha256": sha256(source_bundle),
        "source_bundle_verified_revisions": [revision],
        "ceiling_evidence_sha256": sha256(f113),
        "ceiling_publication_audit_sha256": sha256(f113_audit),
        "storage_evidence_sha256": sha256(storage_evidence),
        "reconciliation_sha256": sha256(reconciliation),
        "ledger_sha256": sha256(ledger),
        "reservations_sha256": sha256(reservations),
        "scheduler_evidence": scheduler_evidence,
        "predecessor_recost_sha256": sha256(predecessor),
        "predecessor_recost_publication_audit_sha256": sha256(predecessor_audit),
        "authenticated_lineages_sha256": lineage_sha256,
        "computed_projection_sha256": projection_sha256,
        "r17_readiness_evidence_sha256": None,
    }
    payload = {
        "schema_version": 1,
        "record_type": "stage-i-barrier-recost-checkpoint",
        "checkpoint": "F-114",
        "artifact_name": V2_ARTIFACT,
        "execution_epoch": EPOCH,
        "generated_utc": request_value["generated_utc"],
        "expires_utc": request_value["expires_utc"],
        "scope": request_value["scope"],
        "predecessor_recost": predecessor_binding,
        "authorization": authorization,
        "barrier": {
            "job_ids": ["12345", "23456", "4766856"],
            "recorded_segments": request_value["barrier"]["recorded_segments"],
            "scheduler_evidence": scheduler_evidence,
        },
        "budget": budget,
        "storage": {
            "available_bytes": 10000,
            "retained_stage_i_bytes": 0,
            "required_safety_bytes": 1000,
            "projected_authorized_wave_growth_bytes": projected_storage,
            "headroom_after_authorized_wave_and_safety_bytes": (
                10000 - 1000 - projected_storage
            ),
        },
        "ledger": {
            "rows": 0,
            "sha256": sha256(ledger),
            "cumulative_stage_i_node_hours": "0",
        },
        "reservations": {
            "rows": 0,
            "sha256": sha256(reservations),
            "active": 0,
        },
        "manifests": {
            "rows": 0,
            "bindings": [],
            "authenticated_lineages_sha256": lineage_sha256,
        },
        "r17_readiness": None,
        "promoted_f113": {
            "path": str(f113),
            "sha256": sha256(f113),
            "publication_audit_path": str(f113_audit),
            "publication_audit_sha256": sha256(f113_audit),
            "publication_audit": f113_audit_value,
        },
        "reconcile": reconcile,
        "provenance": provenance,
    }
    fixture["staged"].unlink()
    fixture["staged"] = fixture["accounting"] / f"{V2_ARTIFACT}.staged"
    write_json(fixture["staged"], payload)
    fixture["canonical"] = fixture["accounting"] / V2_ARTIFACT
    fixture["audit"] = fixture["accounting"] / f"{V2_ARTIFACT}.publication_audit.json"
    fixture["recost_forensics"] = fixture["recost_forensics_root"] / V2_ARTIFACT
    fixture.update(
        {
            "artifact_name": V2_ARTIFACT,
            "authorization": authorization,
            "scheduler_evidence": scheduler_evidence,
            "scheduler_files": [fixture["scheduler"], second_scheduler, r12_scheduler],
            "source_bundle": source_bundle,
            "source_bundle_sha256": sha256(source_bundle),
            "stage_i_revision": revision,
            "generator_revision": revision,
            "source_bundle_revisions": [revision],
            "request": request,
            "request_sha256": sha256(request),
            "inputs": inputs,
            "evidence_files": {
                "reconciliation": reconciliation,
                "ledger": ledger,
                "reservations": reservations,
                "storage": storage_evidence,
                "f113": f113,
                "f113_audit": f113_audit,
                "predecessor": predecessor,
                "predecessor_audit": predecessor_audit,
                "matrix": matrix,
            },
            "artifact_sha256": sha256(fixture["staged"]),
        }
    )
    return fixture


def replace_bounded_artifact(fixture, payload: dict[str, object]) -> None:
    """Replace one staged bounded artifact and refresh its external digest."""

    write_json(fixture["staged"], payload)
    fixture["artifact_sha256"] = sha256(fixture["staged"])


def mutate_bounded_artifact(fixture, mutate) -> None:
    """Mutate one V2 artifact while retaining an exact external digest."""

    payload = json.loads(fixture["staged"].read_text())
    mutate(payload)
    replace_bounded_artifact(fixture, payload)


def mutate_bounded_authorization(fixture, mutate) -> None:
    """Mutate the expected and retained bounded authorization together."""

    authorization = json.loads(json.dumps(fixture["authorization"]))
    mutate(authorization)
    payload = json.loads(fixture["staged"].read_text())
    payload["authorization"] = authorization
    fixture["authorization"] = authorization
    replace_bounded_artifact(fixture, payload)


def mutate_bounded_scheduler_evidence(fixture, mutate) -> None:
    """Mutate both retained scheduler-evidence lists and the expected packet."""

    evidence = json.loads(json.dumps(fixture["scheduler_evidence"]))
    mutate(evidence)
    payload = json.loads(fixture["staged"].read_text())
    payload["provenance"]["scheduler_evidence"] = evidence
    payload["barrier"]["scheduler_evidence"] = evidence
    fixture["scheduler_evidence"] = evidence
    replace_bounded_artifact(fixture, payload)


def checkpoint_command(fixture, action: str, *extra: str,
                       utility_sha256: str | None = None,
                       artifact_sha256: str | None = None,
                       root: Path | None = None,
                       queue_file: Path | None = None,
                       include_queue_fixture: bool = True,
                       generator_relative_path: str | None = None) -> list[str]:
    """Build one fully bound local companion invocation."""

    counts = fixture["counts"]
    if "recommendations" in fixture:
        packet = [
            "--recost-recommendations-json",
            json.dumps(fixture["recommendations"], sort_keys=True),
            "--v2-scheduler-evidence-json",
            json.dumps(fixture["scheduler_evidence"], sort_keys=True),
            "--recost-request-relative-path",
            str(fixture["request"].relative_to(fixture["root"])),
            "--expected-request-sha256",
            fixture["request_sha256"],
            "--independent-review-relative-path",
            str(fixture["independent_review"].relative_to(fixture["root"])),
            "--expected-independent-review-sha256",
            fixture["independent_review_sha256"],
            "--source-bundle-relative-path",
            str(fixture["source_bundle"].relative_to(fixture["root"])),
            "--expected-source-bundle-sha256",
            fixture["source_bundle_sha256"],
            "--expected-stage-i-revision",
            fixture["stage_i_revision"],
            "--expected-generator-revision",
            fixture["generator_revision"],
            "--expected-source-bundle-verified-revisions-json",
            json.dumps(fixture["source_bundle_revisions"]),
            "--expected-artifact-mode",
            "0444",
            "--artifact-authorization-pointer",
            "/recommendations",
        ]
    elif "authorization" in fixture:
        packet = [
            "--authorized-v2-json",
            json.dumps(fixture["authorization"], sort_keys=True),
            "--v2-scheduler-evidence-json",
            json.dumps(fixture["scheduler_evidence"], sort_keys=True),
            "--recost-request-relative-path",
            str(fixture["request"].relative_to(fixture["root"])),
            "--expected-request-sha256",
            fixture["request_sha256"],
            "--source-bundle-relative-path",
            str(fixture["source_bundle"].relative_to(fixture["root"])),
            "--expected-source-bundle-sha256",
            fixture["source_bundle_sha256"],
            "--expected-stage-i-revision",
            fixture["stage_i_revision"],
            "--expected-generator-revision",
            fixture["generator_revision"],
            "--expected-source-bundle-verified-revisions-json",
            json.dumps(fixture["source_bundle_revisions"]),
        ]
    else:
        packet = [
            "--scheduler-relative-path",
            str(fixture["scheduler"].relative_to(fixture["root"])),
            "--expected-scheduler-sha256",
            sha256(fixture["scheduler"]),
            "--authorized-next-segment-profile-json",
            json.dumps(fixture["profile"], sort_keys=True),
            "--expected-source-bundle-verified-revisions-json",
            json.dumps(fixture["source_bundle_revisions"]),
        ]
    command = [
        sys.executable,
        str(fixture["checkpoint"]),
        "--root",
        str(root or fixture["root"]),
        "--allow-local-root",
        "--expected-utility-sha256",
        utility_sha256 or sha256(fixture["checkpoint"]),
        "--expected-stage-i-sha256",
        sha256(fixture["stage_i"]),
        action,
        "--artifact-name",
        fixture.get("artifact_name", ARTIFACT),
        "--expected-artifact-sha256",
        artifact_sha256 or fixture["artifact_sha256"],
        "--generator-relative-path",
        generator_relative_path or str(fixture["generator"].relative_to(fixture["root"])),
        "--expected-generator-sha256",
        sha256(fixture["generator"]),
        *packet,
        "--expected-transactions",
        str(counts["transactions"]),
        "--expected-reservations",
        str(counts["reservations"]),
        "--expected-active-reservations",
        str(counts["active_reservations"]),
        "--expected-ledger-rows",
        str(counts["ledger_rows"]),
        "--expected-manifests",
        str(counts["manifests"]),
        *extra,
    ]
    if include_queue_fixture:
        command[9:9] = ["--squeue-file", str(queue_file or fixture["queue"])]
    return command


def run_checkpoint(fixture, action: str, *extra: str, **kwargs):
    """Run one local companion command without raising on rejection."""

    return subprocess.run(
        checkpoint_command(fixture, action, *extra, **kwargs),
        check=False,
        capture_output=True,
        text=True,
    )


def assert_rejected(completed: subprocess.CompletedProcess, pattern: str) -> None:
    """Require one fail-closed companion rejection."""

    assert completed.returncode == 1
    assert pattern in completed.stderr


def schema2_review_probe(tmp_path: Path, reviewed_utc: str):
    """Create one minimal authenticated schema-2 independent-review probe."""

    root = tmp_path / "review-root"
    accounting = root / "accounting"
    accounting.mkdir(parents=True)
    artifact_name = "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json"
    artifact_sha256 = "a" * 64
    review = accounting / f"{artifact_name}.independent_review.json"
    write_json(
        review,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": reviewed_utc,
            "decision": "approved-for-publication",
            "reviewer": {
                "agent_id": "fixture-independent-schema2-reviewer",
                "independent_from_generator": True,
            },
            "candidate": {
                "path": str(accounting / artifact_name),
                "sha256": artifact_sha256,
            },
            "scope": {"non_authorizing": True},
        },
    )
    review.chmod(0o444)
    args = SimpleNamespace(
        independent_review_relative_path=str(review.relative_to(root)),
        expected_independent_review_sha256=sha256(review),
        expected_artifact_sha256=artifact_sha256,
    )
    return load_checkpoint_module(), root, artifact_name, args


def f116_committed_tools_probe(tmp_path: Path, monkeypatch):
    """Create one exact seven-tool Git commitment for F116 consumer tests."""

    module = load_checkpoint_module()
    repository = tmp_path / "repository"
    repository.mkdir()
    for relative, mode_text in module.F116_REQUIRED_TOOLS.items():
        path = repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"fixture committed tool {relative}\n".encode())
        path.chmod(int(mode_text, 8))
    subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
    subprocess.run(["git", "add", "."], cwd=repository, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=F116 checkpoint fixture",
            "-c",
            "user.email=f116-checkpoint@example.invalid",
            "commit",
            "-q",
            "-m",
            "Commit exact F116 tools",
        ],
        cwd=repository,
        check=True,
    )
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    tools = [
        {
            "path": relative,
            "revision": head,
            "sha256": sha256(repository / relative),
            "mode": module.F116_REQUIRED_TOOLS[relative],
        }
        for relative in sorted(module.F116_REQUIRED_TOOLS)
    ]
    current = {
        "path": "source-archives/final.bundle",
        "sha256": "f" * 64,
        "complete_history": True,
        "head": head,
        "advertised_tip": {"revision": head, "name": "feature/cgl-landau-fluid"},
        "verified_revisions": [head],
        "selected_as_current": True,
        "candidate_path": str(tmp_path / "final.bundle"),
        "subject": "Commit exact F116 tools",
    }
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-current-source-authority-supersession-evidence",
        "checkpoint": "F-116",
        "execution_epoch": EPOCH,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {},
        "predecessor_authorities": {},
        "implementation": {
            "publisher": next(
                item for item in tools if item["path"] == module.F116_PUBLISHER_RELATIVE
            ),
            "committed_tools": tools,
            "intermediate_36140_bundle": {},
            "current_source_bundle": current,
        },
        "source_archive_catalog": {},
        "authorization": {},
        "validation": {},
        "publication_requirements": {},
    }
    final_binding = {
        "path": current["path"],
        "sha256": current["sha256"],
        "verified_revisions": current["verified_revisions"],
    }
    monkeypatch.setattr(
        module,
        "initial_source_path",
        lambda: repository / "scripts/frontier/cgl_lf_stage_i_checkpoint.py",
    )
    return module, evidence, final_binding


def test_checkpoint_rejects_schema2_review_predating_artifact_before_publication(
    tmp_path,
):
    generated = datetime.now(timezone.utc).replace(microsecond=0) - timedelta(hours=1)
    module, root, artifact_name, args = schema2_review_probe(
        tmp_path, (generated - timedelta(seconds=1)).isoformat()
    )
    with pytest.raises(ValueError, match="predates artifact generation"):
        module.schema2_independent_review(
            root,
            artifact_name,
            args,
            artifact_generated_utc=generated.isoformat(),
        )


def test_checkpoint_rejects_future_schema2_review_before_publication(tmp_path):
    now = datetime.now(timezone.utc).replace(microsecond=0)
    generated = now - timedelta(hours=1)
    module, root, artifact_name, args = schema2_review_probe(
        tmp_path, (now + timedelta(hours=1)).isoformat()
    )
    with pytest.raises(ValueError, match="independent review is in the future"):
        module.schema2_independent_review(
            root,
            artifact_name,
            args,
            artifact_generated_utc=generated.isoformat(),
        )


def test_checkpoint_rejects_publication_predating_schema2_review(tmp_path):
    generated = datetime.now(timezone.utc).replace(microsecond=0) - timedelta(hours=1)
    reviewed = generated + timedelta(minutes=30)
    module, root, artifact_name, args = schema2_review_probe(
        tmp_path, reviewed.isoformat()
    )
    with pytest.raises(ValueError, match="publication predates independent review"):
        module.schema2_independent_review(
            root,
            artifact_name,
            args,
            artifact_generated_utc=generated.isoformat(),
            publication_published_utc=(reviewed - timedelta(seconds=1)).isoformat(),
        )


def test_checkpoint_accepts_exact_schema2_review_chronology_boundaries(tmp_path):
    retained = datetime.now(timezone.utc).replace(microsecond=0) - timedelta(minutes=1)
    module, root, artifact_name, args = schema2_review_probe(
        tmp_path, retained.isoformat()
    )
    binding = module.schema2_independent_review(
        root,
        artifact_name,
        args,
        artifact_generated_utc=retained.isoformat(),
        publication_published_utc=retained.isoformat(),
    )
    assert binding["sha256"] == args.expected_independent_review_sha256


def test_checkpoint_declared_independence_disclaims_cryptographic_identity(tmp_path):
    retained = datetime.now(timezone.utc).replace(microsecond=0) - timedelta(minutes=1)
    module, root, artifact_name, args = schema2_review_probe(
        tmp_path, retained.isoformat()
    )
    assurance = module.validate_schema2_independent_review(
        (root / args.independent_review_relative_path).read_bytes(),
        root,
        artifact_name,
        args,
        artifact_generated_utc=retained.isoformat(),
        publication_published_utc=retained.isoformat(),
    )
    assert assurance["cryptographic_identity_verified"] is False
    assert (
        assurance["non_cryptographic_limitation"]
        == module.INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION
    )


def test_checkpoint_accepts_exact_f116_committed_tools(tmp_path, monkeypatch):
    module, evidence, final_binding = f116_committed_tools_probe(tmp_path, monkeypatch)
    module.validate_f116_committed_tools(
        (json.dumps(evidence, indent=2, sort_keys=True) + "\n").encode(), final_binding
    )


def test_checkpoint_f116_required_tools_match_source_authority_contract():
    checkpoint = load_checkpoint_module()
    source = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_source_authority.py"
    spec = importlib.util.spec_from_file_location(
        "_cgl_lf_stage_i_source_authority_contract", source
    )
    assert spec is not None and spec.loader is not None
    authority = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(authority)
    assert checkpoint.F116_REQUIRED_TOOLS == authority.REQUIRED_TOOLS


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "wrong-mode", "wrong-sha"])
def test_checkpoint_rejects_malformed_f116_committed_tools(
    tmp_path, monkeypatch, mutation
):
    module, evidence, final_binding = f116_committed_tools_probe(tmp_path, monkeypatch)
    tools = evidence["implementation"]["committed_tools"]
    if mutation == "missing":
        tools.pop()
    elif mutation == "duplicate":
        tools[1] = dict(tools[0])
    elif mutation == "wrong-mode":
        tools[0]["mode"] = "0755" if tools[0]["mode"] == "0644" else "0644"
    else:
        tools[0]["sha256"] = "0" * 64

    with pytest.raises(ValueError, match="schema-2 F116 committed"):
        module.validate_f116_committed_tools(
            (json.dumps(evidence, indent=2, sort_keys=True) + "\n").encode(),
            final_binding,
        )


def test_checkpoint_rejects_artifact_review_candidate_path_race(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    candidate = tmp_path / "independent-review.candidate.json"
    candidate.write_bytes(b'{"review": "exact candidate"}\n')
    candidate.chmod(0o444)
    expected = sha256(candidate)
    with module.authenticated_artifact_review_candidate(
        candidate, root, expected
    ) as (_, descriptor, retained):
        retired = candidate.with_suffix(".retired")
        candidate.rename(retired)
        candidate.write_bytes(retained)
        candidate.chmod(0o444)
        with pytest.raises(ValueError, match="pathname changed during installation"):
            module.require_artifact_review_candidate_stable(
                candidate, descriptor, expected
            )


def test_checkpoint_rejects_artifact_review_candidate_checksum_drift(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    candidate = tmp_path / "independent-review.candidate.json"
    candidate.write_bytes(b'{"review": "unexpected candidate"}\n')
    candidate.chmod(0o444)
    with pytest.raises(ValueError, match="candidate checksum has changed"):
        with module.authenticated_artifact_review_candidate(
            candidate, root, "0" * 64
        ):
            pass


def test_checkpoint_rejects_relative_artifact_review_candidate(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    with pytest.raises(ValueError, match="absolute normalized path"):
        with module.authenticated_artifact_review_candidate(
            Path("independent-review.candidate.json"), root, "0" * 64
        ):
            pass


def test_checkpoint_rejects_in_root_artifact_review_candidate(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    candidate = root / "independent-review.candidate.json"
    candidate.write_bytes(b'{"review": "in-root candidate"}\n')
    candidate.chmod(0o444)
    with pytest.raises(ValueError, match="must be external to the Stage I root"):
        with module.authenticated_artifact_review_candidate(
            candidate, root, sha256(candidate)
        ):
            pass


def test_checkpoint_artifact_review_scope_is_exactly_non_authorizing(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    artifact_name = "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json"
    expected = "1" * 64
    reviewed = datetime.now(timezone.utc) - timedelta(seconds=1)
    review = {
        "schema_version": 1,
        "record_type": "stage-i-recost-recommendation-independent-review",
        "execution_epoch": EPOCH,
        "reviewed_utc": reviewed.isoformat(),
        "decision": "approved-for-publication",
        "reviewer": {
            "agent_id": "independent-checkpoint-reviewer",
            "independent_from_generator": True,
        },
        "candidate": {
            "path": str(root / "accounting" / artifact_name),
            "sha256": expected,
        },
        "scope": {"non_authorizing": True},
    }
    args = SimpleNamespace(
        expected_artifact_sha256=expected,
        expected_generator_revision="2" * 40,
    )
    retained = (json.dumps(review, sort_keys=True) + "\n").encode()
    assurance = module.validate_schema2_independent_review(
        retained,
        root,
        artifact_name,
        args,
        artifact_generated_utc=(reviewed - timedelta(seconds=1)).isoformat(),
    )
    assert assurance["strict_distinct_role_and_agent_declarations"] is True

    review["scope"]["artifact_review_waiver"] = "forbidden"
    with pytest.raises(ValueError, match="independent review differs"):
        module.validate_schema2_independent_review(
            (json.dumps(review, sort_keys=True) + "\n").encode(),
            root,
            artifact_name,
            args,
            artifact_generated_utc=(reviewed - timedelta(seconds=1)).isoformat(),
        )


def test_checkpoint_rejects_artifact_review_target_create_race(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    target = accounting / "artifact.independent_review.json"
    retained = b'{"review": "exact candidate"}\n'
    expected = hashlib.sha256(retained).hexdigest()
    descriptor = os.open(accounting, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    real_renameat2 = module.renameat2
    raced = False

    def racing_renameat2(parent, source, selected_target, flags, label):
        nonlocal raced
        if (
            not raced
            and parent == descriptor
            and selected_target == target.name
            and flags == module.RENAME_NOREPLACE
        ):
            raced = True
            injected = os.open(
                target.name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o444,
                dir_fd=descriptor,
            )
            os.write(injected, b"raced target\n")
            os.close(injected)
        return real_renameat2(parent, source, selected_target, flags, label)

    monkeypatch.setattr(module, "renameat2", racing_renameat2)
    try:
        with pytest.raises(ValueError, match="already exists or changed"):
            module.create_artifact_review(descriptor, target, retained, expected)
    finally:
        os.close(descriptor)
    assert raced
    assert target.read_bytes() == b"raced target\n"


def test_checkpoint_cleans_only_exact_failed_artifact_review_create(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    target = accounting / "artifact.independent_review.json"
    retained = b'{"review": "exact candidate"}\n'
    expected = hashlib.sha256(retained).hexdigest()
    descriptor = os.open(accounting, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    observed = []

    def fail_write(_descriptor, _retained):
        entries = list(accounting.iterdir())
        observed.append((target.exists(), [entry.lstat().st_mode & 0o777 for entry in entries]))
        raise OSError("simulated artifact-review write failure")

    monkeypatch.setattr(module, "write_descriptor_bytes", fail_write)
    try:
        with pytest.raises(OSError, match="simulated artifact-review write failure"):
            module.create_artifact_review(descriptor, target, retained, expected)
    finally:
        os.close(descriptor)
    assert observed == [(False, [0o600])]
    assert not target.exists()
    assert list(accounting.iterdir()) == []


def test_checkpoint_atomic_retirement_retains_substituted_inode(tmp_path, monkeypatch):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    escaped = tmp_path / "escaped"
    substitute = tmp_path / "substitute"
    victim.write_bytes(b"authenticated victim\n")
    substitute.write_bytes(b"raced substitute\n")
    real_renameat2_between = module.renameat2_between
    raced = False
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    with module.bound_parent_descriptor(victim, "retirement race") as parent:
        expected = os.stat(victim.name, dir_fd=parent, follow_symlinks=False)

        def racing_renameat2_between(
            selected_parent, source, target_parent, target, flags, label
        ):
            nonlocal raced
            if (
                not raced
                and selected_parent == parent
                and source == victim.name
                and flags == module.RENAME_NOREPLACE
                and "retirement" in label
            ):
                raced = True
                os.rename(
                    victim.name,
                    escaped.name,
                    src_dir_fd=parent,
                    dst_dir_fd=parent,
                )
                os.rename(
                    substitute.name,
                    victim.name,
                    src_dir_fd=parent,
                    dst_dir_fd=parent,
                )
            return real_renameat2_between(
                selected_parent, source, target_parent, target, flags, label
            )

        monkeypatch.setattr(module, "renameat2_between", racing_renameat2_between)
        with pytest.raises(ValueError, match="retained as .cgl-checkpoint-retired"):
            module.unlink_bound_entry(parent, victim.name, expected, "raced victim")

    assert raced
    assert not victim.exists()
    assert escaped.read_bytes() == b"authenticated victim\n"
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    assert next(iter(retired)).read_bytes() == b"raced substitute\n"


def test_checkpoint_atomic_retirement_has_no_python_unlink_race_hook(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    victim.write_bytes(b"authenticated victim\n")

    def forbidden_unlink(*_args, **_kwargs):
        raise AssertionError("retirement must not expose a Python unlink race hook")

    monkeypatch.setattr(module.os, "unlink", forbidden_unlink)
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))
    with module.bound_parent_descriptor(victim, "retirement") as parent:
        expected = os.stat(victim.name, dir_fd=parent, follow_symlinks=False)
        module.unlink_bound_entry(parent, victim.name, expected, "authenticated victim")

    assert not victim.exists()
    assert list(tmp_path.iterdir()) == []
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    assert next(iter(retired)).read_bytes() == b"authenticated victim\n"


@pytest.mark.parametrize("drift", ["mode", "link"])
def test_checkpoint_atomic_retirement_rejects_same_inode_security_drift(
    tmp_path, monkeypatch, drift,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    linked = tmp_path / "linked"
    victim.write_bytes(b"authenticated victim\n")
    real_renameat2_between = module.renameat2_between
    raced = False
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    with module.bound_parent_descriptor(victim, "retirement security drift") as parent:
        expected = os.stat(victim.name, dir_fd=parent, follow_symlinks=False)

        def drift_before_retirement(
            source_parent, source, target_parent, target, flags, label
        ):
            nonlocal raced
            if (
                not raced
                and source == victim.name
                and flags == module.RENAME_NOREPLACE
                and "retirement" in label
            ):
                raced = True
                if drift == "mode":
                    os.chmod(victim.name, 0o666, dir_fd=source_parent)
                else:
                    os.link(
                        victim.name,
                        linked.name,
                        src_dir_fd=source_parent,
                        dst_dir_fd=source_parent,
                    )
            return real_renameat2_between(
                source_parent, source, target_parent, target, flags, label
            )

        monkeypatch.setattr(module, "renameat2_between", drift_before_retirement)
        with pytest.raises(ValueError, match="changed during atomic retirement"):
            module.unlink_bound_entry(parent, victim.name, expected, "drifted victim")

    assert raced
    assert not victim.exists()
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    retained = next(iter(retired))
    assert retained.read_bytes() == b"authenticated victim\n"
    if drift == "mode":
        assert retained.stat().st_mode & 0o777 == 0o666
    else:
        assert retained.stat().st_nlink == 2
        assert linked.stat().st_ino == retained.stat().st_ino


def test_checkpoint_atomic_retirement_durably_records_public_name_reappearance(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    victim.write_bytes(b"authenticated victim\n")
    real_renameat2_between = module.renameat2_between
    real_fsync = module.os.fsync
    fsynced_directories = []
    reappeared = False
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced_directories.append(module.profile_identity(profile))
        return result

    with module.bound_parent_descriptor(victim, "retirement reappearance") as parent:
        expected = os.stat(victim.name, dir_fd=parent, follow_symlinks=False)
        parent_identity = module.profile_identity(os.fstat(parent))
        forensic_parent_identity = module.profile_identity(tmp_path.parent.stat())

        def reappear_after_retirement(
            source_parent, source, target_parent, target, flags, label
        ):
            nonlocal reappeared
            result = real_renameat2_between(
                source_parent, source, target_parent, target, flags, label
            )
            if (
                not reappeared
                and source == victim.name
                and flags == module.RENAME_NOREPLACE
                and "retirement" in label
            ):
                reappeared = True
                descriptor = os.open(
                    source,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o644,
                    dir_fd=source_parent,
                )
                os.write(descriptor, b"reappeared public name\n")
                os.close(descriptor)
            return result

        monkeypatch.setattr(module, "renameat2_between", reappear_after_retirement)
        monkeypatch.setattr(module.os, "fsync", observe_fsync)
        with pytest.raises(ValueError, match="public name reappeared during retirement"):
            module.unlink_bound_entry(parent, victim.name, expected, "reappeared victim")

    assert reappeared
    assert victim.read_bytes() == b"reappeared public name\n"
    assert parent_identity in fsynced_directories
    assert forensic_parent_identity in fsynced_directories
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    retained = next(iter(retired))
    assert module.profile_identity(retained.stat()) == module.profile_identity(expected)
    assert retained.read_bytes() == b"authenticated victim\n"


def test_checkpoint_atomic_retirement_durably_records_post_rename_exception(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    victim.write_bytes(b"authenticated victim\n")
    real_renameat2_between = module.renameat2_between
    real_fsync = module.os.fsync
    fsynced_directories = []
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced_directories.append(module.profile_identity(profile))
        return result

    with module.bound_parent_descriptor(victim, "retirement post-rename failure") as parent:
        expected = os.stat(victim.name, dir_fd=parent, follow_symlinks=False)
        parent_identity = module.profile_identity(os.fstat(parent))
        forensic_parent_identity = module.profile_identity(tmp_path.parent.stat())

        def fail_after_retirement(
            source_parent, source, target_parent, target, flags, label
        ):
            real_renameat2_between(
                source_parent, source, target_parent, target, flags, label
            )
            raise OSError("simulated post-rename retirement failure")

        monkeypatch.setattr(module, "renameat2_between", fail_after_retirement)
        monkeypatch.setattr(module.os, "fsync", observe_fsync)
        with pytest.raises(OSError, match="simulated post-rename retirement failure"):
            module.unlink_bound_entry(parent, victim.name, expected, "failed victim")

    assert not victim.exists()
    assert parent_identity in fsynced_directories
    assert forensic_parent_identity in fsynced_directories
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    retained = next(iter(retired))
    assert module.profile_identity(retained.stat()) == module.profile_identity(expected)
    assert retained.read_bytes() == b"authenticated victim\n"


def test_checkpoint_absent_json_publication_does_not_clobber_raced_target(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "target.json"
    real_renameat2 = module.renameat2
    raced = False

    def racing_renameat2(parent, source, selected_target, flags, label):
        nonlocal raced
        if (
            not raced
            and selected_target == target.name
            and flags == module.RENAME_NOREPLACE
            and label == "JSON atomic write"
        ):
            raced = True
            descriptor = os.open(
                target.name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o644,
                dir_fd=parent,
            )
            os.write(descriptor, b"raced target\n")
            os.close(descriptor)
        return real_renameat2(parent, source, selected_target, flags, label)

    monkeypatch.setattr(module, "renameat2", racing_renameat2)
    with pytest.raises(ValueError, match="JSON atomic write target already exists"):
        module.write_json(target, {"must_not": "clobber"})

    assert raced
    assert target.read_bytes() == b"raced target\n"
    temporaries = list(tmp_path.glob(f".{target.name}.*.tmp"))
    assert len(temporaries) == 1
    assert json.loads(temporaries[0].read_text()) == {"must_not": "clobber"}


def test_checkpoint_json_write_rejects_initially_unsafe_parent_before_mutation(
    tmp_path,
):
    module = load_checkpoint_module()
    parent = tmp_path / "unsafe-parent"
    parent.mkdir()
    parent.chmod(0o777)
    target = parent / "state.json"

    with pytest.raises(ValueError, match="exceeds trusted profile 0755"):
        module.write_json(target, {"must_not": "be-created"})

    assert list(parent.iterdir()) == []


def test_checkpoint_json_publication_rejects_forged_same_inode_bytes(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "target.json"
    real_renameat2 = module.renameat2
    forged = None

    def forge_before_rename(parent, source, selected_target, flags, label):
        nonlocal forged
        if forged is None and label == "JSON atomic write":
            descriptor = os.open(source, os.O_WRONLY | os.O_NOFOLLOW, dir_fd=parent)
            try:
                forged = b"x" * os.fstat(descriptor).st_size
                os.pwrite(descriptor, forged, 0)
            finally:
                os.close(descriptor)
        return real_renameat2(parent, source, selected_target, flags, label)

    monkeypatch.setattr(module, "renameat2", forge_before_rename)
    with pytest.raises(ValueError, match="content digest changed during mutation"):
        module.write_json(target, {"must": "remain authenticated"})

    assert forged is not None
    assert target.read_bytes() == forged


def test_checkpoint_json_temporary_cleanup_retires_mode_zero_crash_remnant(tmp_path):
    module = load_checkpoint_module()
    target = tmp_path / "state.json"
    temporary = tmp_path / f".{target.name}.123.{'0' * 32}.tmp"
    temporary.write_bytes(b"interrupted JSON payload\n")
    temporary.chmod(0o000)
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    module.remove_json_temporaries(
        module.json_temporary_entries(target),
        "JSON crash temporary",
    )

    assert not temporary.exists()
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    assert next(iter(retired)).stat().st_mode & 0o777 == 0o000


def test_checkpoint_json_temporary_cleanup_fails_closed_on_substitution(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "state.json"
    temporary = tmp_path / f".{target.name}.123.{'0' * 32}.tmp"
    escaped = tmp_path / "escaped"
    substitute = tmp_path / "substitute"
    temporary.write_bytes(b"authenticated JSON temporary\n")
    temporary.chmod(0o000)
    substitute.write_bytes(b"raced substitute\n")
    substitute.chmod(0o000)
    real_renameat2_between = module.renameat2_between
    raced = False
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    def racing_renameat2_between(
        source_parent, source, target_parent, selected_target, flags, label
    ):
        nonlocal raced
        if (
            not raced
            and source == temporary.name
            and flags == module.RENAME_NOREPLACE
            and "retirement" in label
        ):
            raced = True
            os.rename(
                temporary.name,
                escaped.name,
                src_dir_fd=source_parent,
                dst_dir_fd=source_parent,
            )
            os.rename(
                substitute.name,
                temporary.name,
                src_dir_fd=source_parent,
                dst_dir_fd=source_parent,
            )
        return real_renameat2_between(
            source_parent, source, target_parent, selected_target, flags, label
        )

    monkeypatch.setattr(module, "renameat2_between", racing_renameat2_between)
    with pytest.raises(ValueError, match="changed during atomic retirement"):
        module.remove_json_temporaries(
            module.json_temporary_entries(target),
            "JSON crash temporary",
        )

    assert raced
    assert not temporary.exists()
    assert escaped.exists()
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    assert next(iter(retired)).stat().st_mode & 0o777 == 0o000


def test_checkpoint_json_forward_replacement_does_not_clobber_reappeared_target(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "target.json"
    substitute = tmp_path / "substitute"
    target.write_bytes(b"authenticated predecessor\n")
    substitute.write_bytes(b"raced substitute\n")
    predecessor = target.read_bytes()
    real_rename_bound_noreplace = module.rename_bound_noreplace
    raced = False

    def race_before_publication(parent, source, selected_target, expected, label, **kwargs):
        nonlocal raced
        if not raced and selected_target == target.name and label == "JSON atomic write":
            raced = True
            os.rename(
                substitute.name,
                target.name,
                src_dir_fd=parent,
                dst_dir_fd=parent,
            )
        return real_rename_bound_noreplace(
            parent, source, selected_target, expected, label, **kwargs
        )

    monkeypatch.setattr(module, "rename_bound_noreplace", race_before_publication)
    with pytest.raises(ValueError, match="target already exists"):
        module.write_json(target, {"new": "payload"})

    assert raced
    assert target.read_bytes() == b"raced substitute\n"
    temporaries = list(tmp_path.glob(f".{target.name}.*.tmp"))
    assert len(temporaries) == 1
    assert json.loads(temporaries[0].read_text()) == {"new": "payload"}
    retired = list(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))
    assert any(path.read_bytes() == predecessor for path in retired)


def test_checkpoint_json_forward_replacement_operates_when_renameat2_is_unsupported(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "target.json"
    target.write_bytes(b"authenticated predecessor\n")
    predecessor = target.read_bytes()

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    module.write_json(target, {"new": "payload"})

    assert json.loads(target.read_text()) == {"new": "payload"}
    assert not list(tmp_path.glob(f".{target.name}.*.tmp"))
    retired = list(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))
    assert any(path.read_bytes() == predecessor for path in retired)


def test_checkpoint_forensic_publication_does_not_clobber_raced_target(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "forensic"
    source.write_bytes(b"authenticated forensic bytes\n")
    source.chmod(0o444)
    expected = sha256(source)
    real_renameat2 = module.renameat2
    raced = False

    def racing_renameat2(parent, temporary, selected_target, flags, label):
        nonlocal raced
        if (
            not raced
            and selected_target == target.name
            and flags == module.RENAME_NOREPLACE
            and label == "recost forensic publication"
        ):
            raced = True
            descriptor = os.open(
                target.name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o444,
                dir_fd=parent,
            )
            os.write(descriptor, b"raced forensic target\n")
            os.close(descriptor)
        return real_renameat2(parent, temporary, selected_target, flags, label)

    monkeypatch.setattr(module, "renameat2", racing_renameat2)
    with pytest.raises(ValueError, match="recost forensic publication target already exists"):
        module.copy_forensic(source, target, expected, expected_mode=0o444)

    assert raced
    assert target.read_bytes() == b"raced forensic target\n"
    retained = list(tmp_path.glob(".cgl-checkpoint-forensic-*"))
    assert len(retained) == 1
    assert retained[0].read_bytes() == source.read_bytes()


def test_checkpoint_forensic_publication_rejects_forged_same_inode_bytes(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "forensic"
    source.write_bytes(b"authenticated forensic bytes\n")
    source.chmod(0o444)
    expected = sha256(source)
    real_renameat2 = module.renameat2
    forged = None

    def forge_before_rename(parent, temporary, selected_target, flags, label):
        nonlocal forged
        if forged is None and label == "recost forensic publication":
            os.chmod(temporary, 0o644, dir_fd=parent)
            descriptor = os.open(temporary, os.O_WRONLY | os.O_NOFOLLOW, dir_fd=parent)
            try:
                forged = b"x" * os.fstat(descriptor).st_size
                os.pwrite(descriptor, forged, 0)
            finally:
                os.close(descriptor)
                os.chmod(temporary, 0o444, dir_fd=parent)
        return real_renameat2(parent, temporary, selected_target, flags, label)

    monkeypatch.setattr(module, "renameat2", forge_before_rename)
    with pytest.raises(ValueError, match="content digest changed during mutation"):
        module.copy_forensic(source, target, expected, expected_mode=0o444)

    assert forged is not None
    assert target.read_bytes() == forged


def test_checkpoint_forensic_copy_recovers_deterministic_owner_only_temporary(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "forensic"
    source.write_bytes(b"authenticated forensic bytes\n")
    source.chmod(0o444)
    expected = sha256(source)
    real_fchmod = module.os.fchmod

    def interrupt_before_forensic_mode_publication(descriptor, mode):
        if mode == 0o444:
            raise OSError("simulated forensic interruption at mode 0600")
        return real_fchmod(descriptor, mode)

    with monkeypatch.context() as patch:
        patch.setattr(module.os, "fchmod", interrupt_before_forensic_mode_publication)
        with pytest.raises(OSError, match="simulated forensic interruption"):
            module.copy_forensic(source, target, expected, expected_mode=0o444)

    temporary = tmp_path / module.forensic_temporary_name(target)
    assert temporary.exists()
    assert temporary.stat().st_mode & 0o777 == 0o600

    module.copy_forensic(source, target, expected, expected_mode=0o444)

    assert target.read_bytes() == source.read_bytes()
    assert target.stat().st_mode & 0o777 == 0o444
    assert not temporary.exists()


def test_checkpoint_forensic_recovery_rejects_same_inode_forgery_during_fsync(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    destination = tmp_path / "forensic"
    source.write_bytes(b"authenticated forensic bytes\n")
    source.chmod(0o444)
    destination.write_bytes(source.read_bytes())
    destination.chmod(0o444)
    expected = sha256(source)
    identity = module.profile_identity(destination.stat())
    real_fsync = module.os.fsync
    forged = False

    def forge_during_directory_fsync(descriptor):
        nonlocal forged
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode) and not forged:
            forged = True
            destination.chmod(0o644)
            writable = os.open(destination, os.O_WRONLY | os.O_NOFOLLOW)
            try:
                os.pwrite(writable, b"x" * destination.stat().st_size, 0)
            finally:
                os.close(writable)
                destination.chmod(0o444)
        return result

    monkeypatch.setattr(module.os, "fsync", forge_during_directory_fsync)
    with pytest.raises(ValueError, match="content digest changed during mutation"):
        module.recover_or_copy_forensic(
            source,
            destination,
            expected,
            expected_mode=0o444,
            simulate_interruption_before_directory_fsync=False,
        )

    assert forged
    assert module.profile_identity(destination.stat()) == identity
    assert sha256(destination) != expected


def test_checkpoint_forensic_recovery_rejects_moved_namespace_during_fsync(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    sources = tmp_path / "sources"
    sources.mkdir()
    source = sources / "source"
    source.write_bytes(b"authenticated forensic bytes\n")
    source.chmod(0o444)
    live = tmp_path / "live"
    live.mkdir()
    destination = live / "forensic"
    destination.write_bytes(source.read_bytes())
    destination.chmod(0o444)
    detached = tmp_path / "detached"
    real_fsync = module.os.fsync
    moved = False

    def move_during_directory_fsync(descriptor):
        nonlocal moved
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode) and not moved:
            moved = True
            live.rename(detached)
            live.mkdir()
        return result

    monkeypatch.setattr(module.os, "fsync", move_during_directory_fsync)
    with pytest.raises(ValueError, match="parent path changed during mutation"):
        module.recover_or_copy_forensic(
            source,
            destination,
            sha256(source),
            expected_mode=0o444,
            simulate_interruption_before_directory_fsync=False,
        )

    assert moved
    assert not destination.exists()
    assert (detached / destination.name).read_bytes() == source.read_bytes()
    assert list(live.iterdir()) == []


def test_checkpoint_forensic_temporary_cleanup_fails_closed_on_substitution(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    destination = tmp_path / "forensic"
    temporary = tmp_path / module.forensic_temporary_name(destination)
    escaped = tmp_path / "escaped"
    substitute = tmp_path / "substitute"
    temporary.write_bytes(b"authenticated forensic temporary\n")
    temporary.chmod(0o000)
    substitute.write_bytes(b"raced substitute\n")
    substitute.chmod(0o000)
    real_renameat2_between = module.renameat2_between
    raced = False
    before = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))

    def racing_renameat2_between(
        source_parent, source, target_parent, selected_target, flags, label
    ):
        nonlocal raced
        if (
            not raced
            and source == temporary.name
            and flags == module.RENAME_NOREPLACE
            and "retirement" in label
        ):
            raced = True
            os.rename(
                temporary.name,
                escaped.name,
                src_dir_fd=source_parent,
                dst_dir_fd=source_parent,
            )
            os.rename(
                substitute.name,
                temporary.name,
                src_dir_fd=source_parent,
                dst_dir_fd=source_parent,
            )
        return real_renameat2_between(
            source_parent, source, target_parent, selected_target, flags, label
        )

    monkeypatch.setattr(module, "renameat2_between", racing_renameat2_between)
    with pytest.raises(ValueError, match="changed during atomic retirement"):
        module.remove_forensic_temporaries(tmp_path, "recost forensic temporary")

    assert raced
    assert not temporary.exists()
    assert escaped.exists()
    retired = set(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    assert next(iter(retired)).stat().st_mode & 0o777 == 0o000


def expose_legacy_canonical(fixture) -> None:
    """Replace one staged fixture with an unattested canonical-only artifact."""

    fixture["staged"].rename(fixture["canonical"])


def bind_legacy_source_bundle(fixture) -> Path:
    """Embed one exact F114-compatible source bundle in a sole-profile fixture."""

    return fixture["source_bundle"]


def test_checkpoint_promotes_and_audits_local_recost(recost_fixture):
    fixture = recost_fixture
    staged = run_checkpoint(fixture, "verify-staged-recost")
    assert staged.returncode == 0, staged.stderr
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    assert fixture["canonical"].is_file()
    assert not fixture["staged"].exists()
    assert fixture["audit"].is_file()
    audit = json.loads(fixture["audit"].read_text())
    assert audit["utility"]["committed"] is False
    assert list(fixture["recost_transactions"].iterdir()) == []
    forensic = list(fixture["recost_forensics"].iterdir())
    assert len(forensic) == 1
    assert sha256(forensic[0]) == sha256(fixture["canonical"])
    verified = run_checkpoint(fixture, "verify-promoted-recost")
    assert verified.returncode == 0, verified.stderr
    audited = run_checkpoint(
        fixture, "audit-recost", "--publication-state", "promoted"
    )
    assert audited.returncode == 0, audited.stderr


def test_checkpoint_authenticates_legacy_f114_source_bundle(recost_fixture):
    fixture = recost_fixture
    bind_legacy_source_bundle(fixture)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "missing",
    [
        ("source_bundle",),
        ("source_bundle_sha256",),
        ("source_bundle", "source_bundle_sha256"),
    ],
)
def test_checkpoint_requires_complete_legacy_f114_source_bundle(
    recost_fixture,
    missing,
):
    fixture = recost_fixture
    for key in missing:
        fixture["profile"].pop(key)
    payload = json.loads(fixture["staged"].read_text())
    payload["authorization"]["sole_next_segment_profile"] = fixture["profile"]
    write_json(fixture["staged"], payload)
    fixture["artifact_sha256"] = sha256(fixture["staged"])
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(
        completed,
        "sole-profile mode requires a complete legacy source-bundle binding",
    )


def test_checkpoint_requires_legacy_f114_source_bundle_revisions(recost_fixture):
    fixture = recost_fixture
    command = checkpoint_command(fixture, "verify-staged-recost")
    option = command.index("--expected-source-bundle-verified-revisions-json")
    del command[option:option + 2]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    assert_rejected(
        completed,
        "sole-profile mode requires source-bundle verified revisions",
    )


def test_checkpoint_rejects_legacy_source_revisions_without_helper_binding(
    recost_fixture,
):
    fixture = recost_fixture
    fixture["source_bundle_revisions"] = ["f" * 40]
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(
        completed,
        "legacy source-bundle revisions do not bind the Stage I helper",
    )


def test_checkpoint_rejects_legacy_f114_source_bundle_digest_drift(recost_fixture):
    fixture = recost_fixture
    bundle = bind_legacy_source_bundle(fixture)
    bundle.write_bytes(bundle.read_bytes() + b"\nforged drift\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "legacy sole-profile source bundle checksum has changed")


def test_checkpoint_rejects_invalid_legacy_f114_git_bundle(recost_fixture):
    fixture = recost_fixture
    bundle = bind_legacy_source_bundle(fixture)
    bundle.write_text("not a git bundle\n")
    fixture["profile"]["source_bundle_sha256"] = sha256(bundle)
    payload = json.loads(fixture["staged"].read_text())
    payload["authorization"]["sole_next_segment_profile"] = fixture["profile"]
    write_json(fixture["staged"], payload)
    fixture["artifact_sha256"] = sha256(fixture["staged"])
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "source bundle header is invalid or unexpectedly large")


def test_checkpoint_rejects_legacy_bundle_without_isolated_revision_coverage(
    recost_fixture,
    tmp_path,
):
    fixture = recost_fixture
    unrelated = tmp_path / "unrelated-repository"
    unrelated.mkdir()
    (unrelated / "unrelated.txt").write_text("unrelated history\n")
    subprocess.run(["git", "init", "-q"], cwd=unrelated, check=True)
    subprocess.run(["git", "add", "."], cwd=unrelated, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=CGL fixture",
            "-c",
            "user.email=cgl-fixture@example.invalid",
            "commit",
            "-q",
            "-m",
            "Create unrelated bundle history",
        ],
        cwd=unrelated,
        check=True,
    )
    bundle = fixture["root"] / "source-archives/unrelated.bundle"
    subprocess.run(
        ["git", "bundle", "create", str(bundle), "--all"],
        cwd=unrelated,
        check=True,
    )
    bundle.chmod(0o644)
    fixture["profile"]["source_bundle"] = str(bundle)
    fixture["profile"]["source_bundle_sha256"] = sha256(bundle)
    payload = json.loads(fixture["staged"].read_text())
    payload["authorization"]["sole_next_segment_profile"] = fixture["profile"]
    write_json(fixture["staged"], payload)
    fixture["artifact_sha256"] = sha256(fixture["staged"])
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "source bundle does not contain requested revision")


def test_checkpoint_rejects_legacy_bundle_path_replacement_during_verification(
    recost_fixture,
    monkeypatch,
):
    fixture = recost_fixture
    replacement = fixture["root"] / "source-archives/replacement.bundle"
    subprocess.run(
        ["git", "bundle", "create", str(replacement), "--all"],
        cwd=fixture["repository"],
        check=True,
    )
    replacement.chmod(0o644)
    module = load_checkpoint_module()
    args = module.parser().parse_args(
        checkpoint_command(fixture, "verify-staged-recost")[2:]
    )
    original = module.git_descriptor_run
    swapped = False

    def replace_path(repository, descriptor, arguments, *, capture_output=False):
        nonlocal swapped
        if not swapped:
            os.replace(replacement, fixture["source_bundle"])
            swapped = True
        return original(
            repository,
            descriptor,
            arguments,
            capture_output=capture_output,
        )

    monkeypatch.setattr(module, "git_descriptor_run", replace_path)
    monkeypatch.setattr(module, "initial_source_path", lambda: fixture["checkpoint"])
    monkeypatch.setattr(module, "repository_root", lambda _source: fixture["repository"])
    with pytest.raises(ValueError, match="source bundle pathname changed"):
        module.authenticate_legacy_profile_source_bundle(fixture["root"], args)


def test_checkpoint_rejects_legacy_bundle_mode_change_during_verification(
    recost_fixture,
    monkeypatch,
):
    fixture = recost_fixture
    module = load_checkpoint_module()
    args = module.parser().parse_args(
        checkpoint_command(fixture, "verify-staged-recost")[2:]
    )
    original = module.git_descriptor_run
    changed = False

    def change_mode(repository, descriptor, arguments, *, capture_output=False):
        nonlocal changed
        if not changed:
            fixture["source_bundle"].chmod(0o600)
            changed = True
        return original(
            repository,
            descriptor,
            arguments,
            capture_output=capture_output,
        )

    monkeypatch.setattr(module, "git_descriptor_run", change_mode)
    monkeypatch.setattr(module, "initial_source_path", lambda: fixture["checkpoint"])
    monkeypatch.setattr(module, "repository_root", lambda _source: fixture["repository"])
    with pytest.raises(ValueError, match="source bundle mode is 0600, expected 0644"):
        module.authenticate_legacy_profile_source_bundle(fixture["root"], args)


def test_checkpoint_adopts_and_audits_legacy_canonical_recost(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    adopted = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert adopted.returncode == 0, adopted.stderr
    audit = json.loads(fixture["audit"].read_text())
    assert audit["record_type"] == "legacy-canonical-adoption"
    assert audit["original_publication_transition_observed"] is False
    assert audit["original_publication_method"] == "unknown"
    assert audit["original_publisher"] == "unknown"
    assert audit["present_authentication"] == {
        "canonical_artifact": "authenticated-under-stage-i-lock",
        "no_queued_cgl_jobs": "required-under-stage-i-lock",
        "reconciliation": "clean-required-under-stage-i-lock",
    }
    assert list(fixture["recost_transactions"].iterdir()) == []
    forensic = list(fixture["recost_forensics"].iterdir())
    assert len(forensic) == 1
    assert sha256(forensic[0]) == sha256(fixture["canonical"])
    verified = run_checkpoint(fixture, "verify-promoted-recost")
    assert verified.returncode == 0, verified.stderr


@pytest.mark.parametrize(
    ("hook", "pattern"),
    [
        (
            "--simulate-adoption-interruption-after-journal",
            "simulated interruption after legacy adoption journal",
        ),
        (
            "--simulate-adoption-interruption-before-forensic-directory-fsync",
            "simulated interruption before forensic directory fsync",
        ),
        (
            "--simulate-adoption-interruption-after-forensic",
            "simulated interruption after legacy adoption forensic copy",
        ),
        (
            "--simulate-adoption-interruption-before-audit-directory-fsync",
            "simulated interruption before JSON directory fsync",
        ),
        (
            "--simulate-adoption-interruption-after-audit",
            "simulated interruption after legacy adoption audit",
        ),
    ],
)
def test_checkpoint_resumes_interrupted_legacy_canonical_adoption(recost_fixture,
                                                                  hook, pattern):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    interrupted = run_checkpoint(fixture, "adopt-legacy-canonical", hook)
    assert_rejected(interrupted, pattern)
    assert len(list(fixture["recost_transactions"].iterdir())) == 1
    resumed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert resumed.returncode == 0, resumed.stderr
    assert list(fixture["recost_transactions"].iterdir()) == []
    assert json.loads(fixture["audit"].read_text())["record_type"] == (
        "legacy-canonical-adoption"
    )


@pytest.mark.parametrize(
    ("hook", "pattern", "linked_record"),
    [
        (
            "--simulate-adoption-interruption-after-journal",
            "simulated interruption after legacy adoption journal",
            "journal",
        ),
        (
            "--simulate-adoption-interruption-after-forensic",
            "simulated interruption after legacy adoption forensic copy",
            "journal",
        ),
        (
            "--simulate-adoption-interruption-before-audit-directory-fsync",
            "simulated interruption before JSON directory fsync",
            "audit",
        ),
    ],
)
def test_checkpoint_lustre_recovers_exact_linked_json_lifecycle_state(
    recost_fixture, hook, pattern, linked_record,
):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    force_checkpoint_renameat2_einval(fixture["checkpoint"])
    interrupted = run_checkpoint(fixture, "adopt-legacy-canonical", hook)
    assert_rejected(interrupted, pattern)

    if linked_record == "audit":
        public = fixture["audit"]
    else:
        public = next(
            path
            for path in fixture["recost_transactions"].iterdir()
            if path.name.endswith(".json")
        )
    temporary = public.parent / f".{public.name}.123.{'0' * 32}.tmp"
    os.link(public, temporary)
    assert temporary.stat().st_ino == public.stat().st_ino
    assert temporary.stat().st_nlink == public.stat().st_nlink == 2

    resumed = run_checkpoint(fixture, "adopt-legacy-canonical")

    assert resumed.returncode == 0, resumed.stderr
    assert not temporary.exists()
    assert fixture["audit"].stat().st_nlink == 1
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_recovers_mode_zero_prejournal_legacy_adoption_temporary(
    recost_fixture,
):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    fixture["recost_transactions"].mkdir()
    fixture["recost_forensics"].mkdir(parents=True)
    temporary = (
        fixture["recost_transactions"]
        / f".interrupted.json.123.{'0' * 32}.tmp"
    )
    temporary.write_bytes(b"partial adoption journal")
    temporary.chmod(0o000)

    resumed = run_checkpoint(fixture, "adopt-legacy-canonical")

    assert resumed.returncode == 0, resumed.stderr
    assert list(fixture["recost_transactions"].iterdir()) == []
    assert json.loads(fixture["audit"].read_text())["record_type"] == (
        "legacy-canonical-adoption"
    )


def test_checkpoint_repairs_interrupted_partial_legacy_adoption_forensic_copy(
    recost_fixture,
):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    interrupted = run_checkpoint(
        fixture,
        "adopt-legacy-canonical",
        "--simulate-adoption-interruption-after-journal",
    )
    assert_rejected(interrupted, "simulated interruption after legacy adoption journal")
    journal = next(fixture["recost_transactions"].iterdir())
    transaction_id = json.loads(journal.read_text())["transaction_id"]
    forensic = (
        fixture["recost_forensics"]
        / f"{transaction_id}.{fixture['canonical'].name}.forensic"
    )
    forensic.write_bytes(b"partial forensic copy")
    forensic.chmod(0o444)
    resumed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert resumed.returncode == 0, resumed.stderr
    assert sha256(forensic) == sha256(fixture["canonical"])


def test_checkpoint_rejects_untrusted_partial_legacy_adoption_forensic_profile(
    recost_fixture,
):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    interrupted = run_checkpoint(
        fixture,
        "adopt-legacy-canonical",
        "--simulate-adoption-interruption-after-journal",
    )
    assert_rejected(interrupted, "simulated interruption after legacy adoption journal")
    journal = next(fixture["recost_transactions"].iterdir())
    transaction_id = json.loads(journal.read_text())["transaction_id"]
    forensic = (
        fixture["recost_forensics"]
        / f"{transaction_id}.{fixture['canonical'].name}.forensic"
    )
    forensic.write_bytes(b"partial forensic copy")
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "recost forensic copy mode is 0644, expected 0444")


def test_checkpoint_legacy_adoption_rejects_staged_twin(recost_fixture):
    fixture = recost_fixture
    os.link(fixture["staged"], fixture["canonical"])
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "recost artifact has 2 links, expected 1")
    assert not fixture["audit"].exists()


def test_checkpoint_legacy_adoption_rejects_nonempty_queue(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    fixture["queue"].write_text("67890|cgl_other_root_writer|RUNNING\n")
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "another CGL job is queued")
    assert not fixture["audit"].exists()


def test_checkpoint_permits_unrelated_account_queue_job(recost_fixture):
    fixture = recost_fixture
    fixture["queue"].write_text("67890|pic_unrelated|RUNNING\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "row",
    (
        "67890||RUNNING\n",
        "|pic_unrelated|RUNNING\n",
        "67890|pic_unrelated|\n",
        "67890| pic_unrelated|RUNNING\n",
        " 67890|pic_unrelated|RUNNING\n",
        "67890|pic_unrelated|RUNNING \n",
    ),
)
def test_checkpoint_rejects_malformed_queue_row(recost_fixture, row):
    fixture = recost_fixture
    fixture["queue"].write_text(row)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "squeue output row has")


def test_checkpoint_legacy_adoption_rejects_audit_relabeling(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    adopted = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert adopted.returncode == 0, adopted.stderr
    audit = json.loads(fixture["audit"].read_text())
    audit["record_type"] = "observed-publication"
    write_json(fixture["audit"], audit)
    completed = run_checkpoint(fixture, "verify-promoted-recost")
    assert_rejected(completed, "recost publication audit schema differs")


def test_checkpoint_legacy_adoption_rejects_preexisting_audit(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    fixture["audit"].write_text("{}\n")
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "recost namespace is not canonical-pending-audit-only")


def test_checkpoint_legacy_adoption_rejects_preexisting_forensic_root(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    fixture["recost_forensics_root"].mkdir()
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "requires an absent recost forensic root")
    assert fixture["recost_transactions"].is_dir()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_legacy_adoption_rejects_forged_journal(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    interrupted = run_checkpoint(
        fixture,
        "adopt-legacy-canonical",
        "--simulate-adoption-interruption-after-journal",
    )
    assert_rejected(interrupted, "simulated interruption after legacy adoption journal")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["operation"] = "observed-publication"
    write_json(journal, record)
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "legacy adoption journal operation binding differs")


def test_checkpoint_legacy_adoption_rejects_altered_forensic_copy(recost_fixture):
    fixture = recost_fixture
    expose_legacy_canonical(fixture)
    interrupted = run_checkpoint(
        fixture,
        "adopt-legacy-canonical",
        "--simulate-adoption-interruption-after-forensic",
    )
    assert_rejected(interrupted, "simulated interruption after legacy adoption forensic copy")
    forensic = next(fixture["recost_forensics"].iterdir())
    forensic.chmod(0o600)
    completed = run_checkpoint(fixture, "adopt-legacy-canonical")
    assert_rejected(completed, "recost forensic copy mode is 0600, expected 0444")


def test_checkpoint_rejects_self_authentication_failure(recost_fixture):
    completed = run_checkpoint(
        recost_fixture, "verify-staged-recost", utility_sha256="0" * 64
    )
    assert_rejected(completed, "utility checksum has changed")


def test_checkpoint_scrubs_git_repository_override_environment(recost_fixture):
    environment = dict(os.environ)
    environment["GIT_DIR"] = "/tmp/checkpoint-fixture-nonexistent-git-dir"
    completed = subprocess.run(
        checkpoint_command(recost_fixture, "verify-staged-recost"),
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert completed.returncode == 0, completed.stderr


def test_checkpoint_scrubs_slurm_routing_environment(monkeypatch):
    monkeypatch.setenv("SLURM_CLUSTERS", "alternate")
    monkeypatch.setenv("SLURM_CONF", "/tmp/untrusted-slurm.conf")
    monkeypatch.setenv("CGL_CHECKPOINT_SENTINEL", "retained")
    monkeypatch.setenv("LD_PRELOAD", "/tmp/untrusted-loader.so")
    monkeypatch.setenv("PYTHONPATH", "/tmp/untrusted-python")
    monkeypatch.setenv("GIT_DIR", "/tmp/untrusted-git")
    environment = load_checkpoint_module().scheduler_environment()
    assert not any(key.startswith("SLURM_") for key in environment)
    assert environment == {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def test_checkpoint_scheduler_child_is_descriptor_bound_and_sanitized(
    monkeypatch, tmp_path
):
    module = load_checkpoint_module()
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    module.require_empty_queue(tmp_path, False, None)

    command, kwargs = calls[0]
    descriptor = int(kwargs["executable"].removeprefix("/proc/self/fd/"))
    assert descriptor in kwargs["pass_fds"]
    assert command[0] == str(module.SQUEUE)
    assert kwargs["env"] == module.scheduler_environment()
    assert kwargs["stdin"] == subprocess.DEVNULL
    assert kwargs["timeout"] == 120


def test_checkpoint_rechecks_queue_directly_before_link(recost_fixture, tmp_path):
    late_queue = tmp_path / "late-squeue.txt"
    late_queue.write_text("67890|cgl_late_root_writer|RUNNING\n")
    completed = run_checkpoint(
        recost_fixture,
        "promote-recost",
        "--pre-link-squeue-file",
        str(late_queue),
    )
    assert_rejected(completed, "another CGL job is queued")
    assert recost_fixture["staged"].is_file()
    assert not recost_fixture["canonical"].exists()
    assert list(recost_fixture["recost_transactions"].iterdir()) == []
    assert list(recost_fixture["recost_forensics"].iterdir()) == []


def test_checkpoint_rejects_shared_queue_fixture_without_opening_it(recost_fixture):
    completed = run_checkpoint(
        recost_fixture,
        "verify-staged-recost",
        queue_file=Path("/lustre/checkpoint-fixture-must-not-open"),
    )
    assert_rejected(completed, "queue fixture must be under local fixture storage /tmp")


def test_checkpoint_rejects_nonlocal_offline_root_without_opening_it(recost_fixture):
    completed = run_checkpoint(
        recost_fixture,
        "verify-staged-recost",
        root=Path("/autofs/checkpoint-fixture-must-not-open"),
    )
    assert_rejected(completed, "offline fixture root must be under local fixture storage /tmp")


def test_checkpoint_requires_offline_queue_fixture(recost_fixture):
    completed = run_checkpoint(
        recost_fixture,
        "verify-staged-recost",
        include_queue_fixture=False,
    )
    assert_rejected(completed, "offline fixture root requires --squeue-file")


def test_checkpoint_rejects_queue_fixture_symlink(recost_fixture, tmp_path):
    queue = tmp_path / "shared-queue-link"
    queue.symlink_to("/lustre/checkpoint-fixture-must-not-resolve")
    completed = run_checkpoint(
        recost_fixture,
        "verify-staged-recost",
        queue_file=queue,
    )
    assert completed.returncode == 1
    assert not recost_fixture["canonical"].exists()


def test_checkpoint_rejects_offline_root_symlink_without_resolving_it(recost_fixture,
                                                                      tmp_path):
    root = tmp_path / "shared-root-link"
    root.symlink_to("/lustre/checkpoint-fixture-must-not-resolve", target_is_directory=True)
    completed = run_checkpoint(recost_fixture, "verify-staged-recost", root=root)
    assert completed.returncode == 1
    assert not recost_fixture["canonical"].exists()


def test_checkpoint_rejects_lock_contention(recost_fixture):
    recost_fixture["lock"].write_text("")
    recost_fixture["lock"].chmod(0o644)
    with recost_fixture["lock"].open("a+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        completed = run_checkpoint(recost_fixture, "promote-recost")
        assert_rejected(completed, "another Stage I mutation holds")
    assert recost_fixture["staged"].is_file()
    assert not recost_fixture["canonical"].exists()


def test_checkpoint_rejects_nonempty_stage_i_lock(recost_fixture):
    fixture = recost_fixture
    fixture["lock"].write_text("forged retained lock content\n")
    fixture["lock"].chmod(0o644)
    completed = run_checkpoint(fixture, "promote-recost")
    assert_rejected(completed, "Stage I lock must be empty")
    assert fixture["staged"].is_file()
    assert not fixture["canonical"].exists()


def test_checkpoint_rejects_stage_i_lock_path_replacement_after_flock(
    monkeypatch,
    tmp_path,
):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock.write_text("")
    lock.chmod(0o644)
    real_flock = module.fcntl.flock
    replaced = False

    def replace_named_lock_after_acquisition(descriptor, operation):
        nonlocal replaced
        real_flock(descriptor, operation)
        if operation == module.fcntl.LOCK_EX | module.fcntl.LOCK_NB and not replaced:
            replaced = True
            lock.unlink()
            lock.write_text("")
            lock.chmod(0o644)

    monkeypatch.setattr(module.fcntl, "flock", replace_named_lock_after_acquisition)
    with pytest.raises(ValueError, match="Stage I lock path changed while locking"):
        with module.promotion_lock({"root": root, "lock": lock}):
            raise AssertionError("replaced lock must not enter the mutation boundary")


def test_checkpoint_lock_replacement_blocks_next_write(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock.write_text("")
    lock.chmod(0o644)
    target = root / "must-not-be-created.json"

    with pytest.raises(ValueError, match="Stage I lock path changed while mutation is active"):
        with module.promotion_lock({"root": root, "lock": lock}):
            lock.rename(root / "retired-lock")
            lock.write_text("")
            lock.chmod(0o644)
            module.write_json(target, {"forbidden": True})

    assert not target.exists()


def test_checkpoint_lustre_noreplace_fallback_uses_hard_link(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"published through hard link\n")
    source_identity = module.profile_identity(source.stat())
    real_linkat = module.linkat_descriptor_noreplace
    linked = []

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def observe_linkat(descriptor, target_parent, target_name, label):
        linked.append((target_name, label))
        return real_linkat(descriptor, target_parent, target_name, label)

    def forbidden_posix_rename(*_args, **_kwargs):
        raise AssertionError("file no-replace fallback must not use POSIX rename")

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "linkat_descriptor_noreplace", observe_linkat)
    monkeypatch.setattr(module.os, "rename", forbidden_posix_rename)
    with module.bound_parent_descriptor(source, "Lustre no-replace") as parent:
        module.rename_bound_noreplace(
            parent,
            source.name,
            target.name,
            source.stat(),
            "Lustre no-replace",
        )

    assert linked == [(target.name, "Lustre no-replace")]
    assert not source.exists()
    assert module.profile_identity(target.stat()) == source_identity
    assert target.stat().st_nlink == 1
    assert target.read_bytes() == b"published through hard link\n"


def test_checkpoint_lustre_noreplace_fallback_rejects_target_race(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"authenticated source\n")
    real_linkat = module.linkat_descriptor_noreplace
    raced = False

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def race_target(descriptor, target_parent, target_name, label):
        nonlocal raced
        if not raced:
            raced = True
            injected = os.open(
                target_name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o644,
                dir_fd=target_parent,
            )
            os.write(injected, b"raced target\n")
            os.close(injected)
        return real_linkat(descriptor, target_parent, target_name, label)

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "linkat_descriptor_noreplace", race_target)
    with module.bound_parent_descriptor(source, "Lustre no-replace race") as parent:
        with pytest.raises(ValueError, match="target already exists"):
            module.rename_bound_noreplace(
                parent,
                source.name,
                target.name,
                source.stat(),
                "Lustre no-replace race",
            )

    assert raced
    assert source.read_bytes() == b"authenticated source\n"
    assert target.read_bytes() == b"raced target\n"


def test_checkpoint_lustre_noreplace_reconciles_post_unlink_exception(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"durably published\n")
    source_identity = module.profile_identity(source.stat())
    real_unlink = module.os.unlink

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def unlink_then_raise(name, *args, **kwargs):
        real_unlink(name, *args, **kwargs)
        if name == source.name:
            raise OSError(errno.EIO, os.strerror(errno.EIO))

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module.os, "unlink", unlink_then_raise)
    with module.bound_parent_descriptor(source, "Lustre unlink reconcile") as parent:
        module.rename_bound_noreplace(
            parent,
            source.name,
            target.name,
            source.stat(),
            "Lustre unlink reconcile",
        )

    assert not source.exists()
    assert module.profile_identity(target.stat()) == source_identity
    assert target.stat().st_nlink == 1


def test_checkpoint_lustre_noreplace_retries_exact_two_link_state(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"retryable publication\n")
    real_unlink = module.os.unlink
    reject_unlink = True

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def interrupt_source_unlink(name, *args, **kwargs):
        if name == source.name and reject_unlink:
            raise OSError(errno.EIO, os.strerror(errno.EIO))
        return real_unlink(name, *args, **kwargs)

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module.os, "unlink", interrupt_source_unlink)
    with module.bound_parent_descriptor(source, "retryable no-replace") as parent:
        with pytest.raises(OSError, match=os.strerror(errno.EIO)):
            module.rename_bound_noreplace(
                parent,
                source.name,
                target.name,
                source.stat(),
                "retryable no-replace",
            )
        assert source.stat().st_nlink == target.stat().st_nlink == 2
        assert source.stat().st_ino == target.stat().st_ino
        reject_unlink = False
        module.rename_bound_noreplace(
            parent,
            source.name,
            target.name,
            source.stat(),
            "retryable no-replace",
        )

    assert not source.exists()
    assert target.stat().st_nlink == 1
    assert target.read_bytes() == b"retryable publication\n"


def test_checkpoint_lustre_noreplace_rejects_same_inode_byte_drift(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"authenticated bytes\n")
    real_linkat = module.linkat_descriptor_noreplace

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def mutate_after_link(descriptor, target_parent, target_name, label):
        result = real_linkat(descriptor, target_parent, target_name, label)
        opened = os.open(source.name, os.O_WRONLY | os.O_NOFOLLOW, dir_fd=target_parent)
        try:
            os.pwrite(opened, b"x" * source.stat().st_size, 0)
        finally:
            os.close(opened)
        return result

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "linkat_descriptor_noreplace", mutate_after_link)
    with module.bound_parent_descriptor(source, "drifted no-replace") as parent:
        with pytest.raises(ValueError, match="content digest changed"):
            module.rename_bound_noreplace(
                parent,
                source.name,
                target.name,
                source.stat(),
                "drifted no-replace",
            )

    assert source.stat().st_ino == target.stat().st_ino
    assert source.stat().st_nlink == target.stat().st_nlink == 2


@pytest.mark.parametrize(
    ("name", "mode"),
    [
        ("state.json", 0o644),
        ("artifact.publication_audit.json", 0o444),
    ],
)
def test_checkpoint_write_json_retries_exact_lustre_post_link_state(
    tmp_path, monkeypatch, name, mode,
):
    module = load_checkpoint_module()
    target = tmp_path / name
    value = {"state": "durable-publication"}

    temporary = leave_lustre_json_post_link_state(
        module, target, value, monkeypatch, mode=mode
    )

    assert temporary.stat().st_ino == target.stat().st_ino
    assert temporary.stat().st_nlink == target.stat().st_nlink == 2
    module.write_json(target, value, mode=mode)

    assert not temporary.exists()
    assert target.stat().st_nlink == 1
    assert target.stat().st_mode & 0o777 == mode
    assert json.loads(target.read_text()) == value


@pytest.mark.parametrize(
    ("state", "pattern"),
    [
        (
            "wrong-name",
            "requires exactly one correctly named temporary",
        ),
        (
            "different-inode",
            "do not select the same exact two-link inode",
        ),
        (
            "public-absent",
            "linked temporary without its exact public name",
        ),
        (
            "extra-link",
            "public name has 3 links",
        ),
    ],
)
def test_checkpoint_write_json_recovery_rejects_nonexact_two_link_states(
    tmp_path, state, pattern,
):
    module = load_checkpoint_module()
    target = tmp_path / "state.json"
    value = {"state": "authenticated"}
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    temporary = tmp_path / f".{target.name}.123.{'0' * 32}.tmp"
    extra = tmp_path / "extra-link"

    if state == "wrong-name":
        target.write_bytes(payload)
        os.link(target, tmp_path / "wrong-temporary")
    elif state == "different-inode":
        target.write_bytes(payload)
        os.link(target, extra)
        temporary.write_bytes(payload)
    elif state == "public-absent":
        temporary.write_bytes(payload)
        os.link(temporary, extra)
    else:
        target.write_bytes(payload)
        os.link(target, temporary)
        os.link(target, extra)

    with pytest.raises(ValueError, match=pattern):
        module.write_json(target, value)

    if target.exists():
        assert target.read_bytes() == payload
    if temporary.exists():
        assert temporary.read_bytes() == payload


def test_checkpoint_lustre_retirement_is_deterministic_and_retryable(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    victim.write_bytes(b"retryable retirement\n")
    real_unlink = module.os.unlink
    reject_unlink = True

    def unsupported_renameat2_between(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def interrupt_source_unlink(name, *args, **kwargs):
        if name == victim.name and reject_unlink:
            raise OSError(errno.EIO, os.strerror(errno.EIO))
        return real_unlink(name, *args, **kwargs)

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2_between)
    monkeypatch.setattr(module.os, "unlink", interrupt_source_unlink)
    with module.bound_parent_descriptor(victim, "retryable retirement") as parent:
        expected = victim.stat()
        retired_name = module.deterministic_retirement_name(
            parent, victim.name, expected
        )
        retired = tmp_path.parent / retired_name
        with pytest.raises(OSError, match=os.strerror(errno.EIO)):
            module.unlink_bound_entry(
                parent, victim.name, expected, "retryable retirement"
            )
        assert victim.stat().st_ino == retired.stat().st_ino
        assert victim.stat().st_nlink == retired.stat().st_nlink == 2
        reject_unlink = False
        module.unlink_bound_entry(
            parent, victim.name, victim.stat(), "retryable retirement"
        )

    assert not victim.exists()
    assert retired.stat().st_nlink == 1
    assert retired.read_bytes() == b"retryable retirement\n"


def test_checkpoint_lustre_retirement_rejects_deterministic_target_collision(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    victim.write_bytes(b"authenticated retirement\n")

    def unsupported_renameat2_between(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2_between)
    with module.bound_parent_descriptor(victim, "retirement collision") as parent:
        expected = victim.stat()
        retired = tmp_path.parent / module.deterministic_retirement_name(
            parent, victim.name, expected
        )
        retired.write_bytes(b"collision\n")
        with pytest.raises(ValueError, match="deterministic forensic name"):
            module.unlink_bound_entry(
                parent, victim.name, expected, "retirement collision"
            )

    assert victim.read_bytes() == b"authenticated retirement\n"
    assert retired.read_bytes() == b"collision\n"


def test_checkpoint_lustre_retirement_mode_zero_fails_closed(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    victim = tmp_path / "victim"
    victim.write_bytes(b"mode-zero forensic bytes\n")
    victim.chmod(0o000)

    def unsupported_renameat2_between(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2_between)
    with module.bound_parent_descriptor(victim, "mode-zero retirement") as parent:
        expected = victim.stat()
        retired = tmp_path.parent / module.deterministic_retirement_name(
            parent, victim.name, expected
        )
        with pytest.raises(ValueError, match="exact readable-byte authentication"):
            module.unlink_bound_entry(
                parent, victim.name, expected, "mode-zero retirement"
            )

    assert victim.exists()
    assert victim.stat().st_mode & 0o777 == 0o000
    assert not retired.exists()


def test_checkpoint_lustre_exchange_fails_closed_without_fallback_mutation(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"new\n")
    target.write_bytes(b"old\n")
    source_identity = module.profile_identity(source.stat())
    target_identity = module.profile_identity(target.stat())

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    def forbidden_mutation(*_args, **_kwargs):
        raise AssertionError("unsupported exchange must not issue a fallback mutation")

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "linkat_descriptor_noreplace", forbidden_mutation)
    monkeypatch.setattr(module.os, "rename", forbidden_mutation)
    with module.bound_parent_descriptor(source, "Lustre exchange") as parent:
        with pytest.raises(ValueError, match="requires RENAME_EXCHANGE support"):
            module.exchange_bound_entries(
                parent,
                source.name,
                target.name,
                source.stat(),
                target.stat(),
                "Lustre exchange",
            )

    assert module.profile_identity(source.stat()) == source_identity
    assert module.profile_identity(target.stat()) == target_identity
    assert source.read_bytes() == b"new\n"
    assert target.read_bytes() == b"old\n"
    assert not list(tmp_path.glob(".cgl-checkpoint-replaced-*"))


def test_checkpoint_exchange_race_does_not_clobber_substitute(tmp_path, monkeypatch):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped"
    source.write_bytes(b"new\n")
    target.write_bytes(b"old\n")
    substitute.write_bytes(b"substitute\n")
    real_renameat2 = module.renameat2

    with module.bound_parent_descriptor(source, "exchange") as parent:
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        target_profile = os.stat(target.name, dir_fd=parent, follow_symlinks=False)

        def race_after_exchange(
            selected_parent, source_name, target_name, flags, label
        ):
            real_renameat2(selected_parent, source_name, target_name, flags, label)
            if flags == module.RENAME_EXCHANGE and "rollback" not in label:
                os.rename(
                    target_name,
                    escaped.name,
                    src_dir_fd=selected_parent,
                    dst_dir_fd=selected_parent,
                )
                os.rename(
                    substitute.name,
                    target_name,
                    src_dir_fd=selected_parent,
                    dst_dir_fd=selected_parent,
                )

        monkeypatch.setattr(module, "renameat2", race_after_exchange)
        with pytest.raises(ValueError, match="names changed after atomic exchange"):
            module.exchange_bound_entries(
                parent,
                source.name,
                target.name,
                source_profile,
                target_profile,
                "raced exchange",
            )

    assert target.read_bytes() == b"substitute\n"
    assert escaped.read_bytes() == b"new\n"
    assert source.read_bytes() == b"old\n"


@pytest.mark.parametrize("operation", ["noreplace", "exchange"])
def test_checkpoint_post_rename_exception_is_durably_recorded(
    tmp_path, monkeypatch, operation,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"new\n")
    if operation == "exchange":
        target.write_bytes(b"old\n")
    real_renameat2 = module.renameat2
    real_fsync = module.os.fsync
    fsynced_directories = []

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced_directories.append(module.profile_identity(profile))
        return result

    def fail_after_successful_rename(parent, source_name, target_name, flags, label):
        real_renameat2(parent, source_name, target_name, flags, label)
        raise OSError("simulated exception after successful namespace mutation")

    monkeypatch.setattr(module.os, "fsync", observe_fsync)
    monkeypatch.setattr(module, "renameat2", fail_after_successful_rename)
    with module.bound_parent_descriptor(source, "post-rename durability") as parent:
        parent_identity = module.profile_identity(os.fstat(parent))
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        if operation == "noreplace":
            with pytest.raises(OSError, match="after successful namespace mutation"):
                module.rename_bound_noreplace(
                    parent, source.name, target.name, source_profile, "durable rename"
                )
        else:
            target_profile = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
            with pytest.raises(OSError, match="after successful namespace mutation"):
                module.exchange_bound_entries(
                    parent,
                    source.name,
                    target.name,
                    source_profile,
                    target_profile,
                    "durable exchange",
                )

    assert parent_identity in fsynced_directories
    assert target.read_bytes() == b"new\n"
    if operation == "noreplace":
        assert not source.exists()
    else:
        assert source.read_bytes() == b"old\n"


def test_checkpoint_exchange_never_rolls_back_after_lock_loss(tmp_path, monkeypatch):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock.write_text("")
    lock.chmod(0o644)
    source = root / "source"
    target = root / "target"
    source.write_bytes(b"new\n")
    target.write_bytes(b"old\n")
    real_renameat2_between = module.renameat2_between
    mutations = []

    def lose_lock_after_exchange(
        source_parent, source_name, target_parent, target_name, flags, label
    ):
        mutations.append(label)
        result = real_renameat2_between(
            source_parent, source_name, target_parent, target_name, flags, label
        )
        if len(mutations) == 1:
            lock.rename(root / "retired-lock")
            lock.write_text("")
            lock.chmod(0o644)
        return result

    monkeypatch.setattr(module, "renameat2_between", lose_lock_after_exchange)
    with pytest.raises(ValueError, match="Stage I lock path changed while mutation is active"):
        with module.promotion_lock({"root": root, "lock": lock}):
            with module.bound_parent_descriptor(source, "lock-loss exchange") as parent:
                module.exchange_bound_entries(
                    parent,
                    source.name,
                    target.name,
                    source.stat(),
                    target.stat(),
                    "lock-loss exchange",
                )

    assert mutations == ["lock-loss exchange"]
    assert source.read_bytes() == b"old\n"
    assert target.read_bytes() == b"new\n"


def test_checkpoint_exchange_never_rolls_back_after_parent_loss(tmp_path, monkeypatch):
    module = load_checkpoint_module()
    live = tmp_path / "live"
    live.mkdir()
    source = live / "source"
    target = live / "target"
    source.write_bytes(b"new\n")
    target.write_bytes(b"old\n")
    detached = tmp_path / "detached"
    real_renameat2_between = module.renameat2_between
    mutations = []

    def detach_parent_after_exchange(
        source_parent, source_name, target_parent, target_name, flags, label
    ):
        mutations.append(label)
        result = real_renameat2_between(
            source_parent, source_name, target_parent, target_name, flags, label
        )
        if len(mutations) == 1:
            live.rename(detached)
            live.mkdir()
        return result

    monkeypatch.setattr(module, "renameat2_between", detach_parent_after_exchange)
    with pytest.raises(ValueError, match="parent path changed during mutation"):
        with module.bound_parent_descriptor(source, "parent-loss exchange") as parent:
            module.exchange_bound_entries(
                parent,
                source.name,
                target.name,
                source.stat(),
                target.stat(),
                "parent-loss exchange",
            )

    assert mutations == ["parent-loss exchange"]
    assert (detached / source.name).read_bytes() == b"old\n"
    assert (detached / target.name).read_bytes() == b"new\n"
    assert list(live.iterdir()) == []


def test_checkpoint_publication_link_uses_authenticated_source_descriptor(tmp_path):
    module = load_checkpoint_module()
    staged = tmp_path / "staged"
    retired = tmp_path / "retired"
    canonical = tmp_path / "canonical"
    staged.write_bytes(b"authenticated staged bytes\n")

    with module.bound_parent_descriptor(staged, "publication") as parent:
        descriptor = os.open(staged.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
        try:
            identity = module.profile_identity(os.fstat(descriptor))
            staged.rename(retired)
            staged.write_bytes(b"raced staged pathname\n")
            module.link_descriptor_noreplace(
                descriptor,
                parent,
                canonical.name,
                identity,
                "canonical publication",
            )
        finally:
            os.close(descriptor)

    assert canonical.read_bytes() == b"authenticated staged bytes\n"
    assert staged.read_bytes() == b"raced staged pathname\n"


def test_checkpoint_rejects_stage_i_lock_through_symlinked_parent(tmp_path):
    module = load_checkpoint_module()
    real_parent = tmp_path / "real-parent"
    real_root = real_parent / "root"
    real_root.mkdir(parents=True)
    alias = tmp_path / "alias"
    alias.symlink_to(real_parent, target_is_directory=True)
    root = alias / "root"
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    with pytest.raises(OSError):
        with module.promotion_lock({"root": root, "lock": lock}):
            raise AssertionError("symlinked lock path must not enter the mutation boundary")
    assert not (real_root / lock.name).exists()


def test_checkpoint_locked_reconcile_uses_in_process_report(monkeypatch, tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock.write_text("")
    lock.chmod(0o644)
    helper = tmp_path / "stage_i_helper.py"
    helper.write_text(
        """
import fcntl
import json
from pathlib import Path
import sys

def reconcile_report(root):
    return {
        "execution_epoch": "E03-forcing-policy",
        "consistent": True,
        "counts": {
            "transactions": 0,
            "reservations": 0,
            "active_reservations": 0,
            "ledger_rows": 0,
            "manifests": 0,
        },
        "issues": [],
    }

if __name__ == "__main__":
    root = Path(sys.argv[sys.argv.index("--root") + 1])
    with (root / ".mks24_stage_i_E03_forcing_policy.lock").open("r+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    print(json.dumps(reconcile_report(root)))
""".lstrip()
    )
    helper.chmod(0o644)

    @contextmanager
    def descriptor(*_args):
        value = os.open(helper, os.O_RDONLY)
        try:
            yield value
        finally:
            os.close(value)

    monkeypatch.setattr(module, "stage_i_descriptor", descriptor)
    args = type(
        "Args",
        (),
        {
            "expected_stage_i_sha256": sha256(helper),
            "expected_transactions": 0,
            "expected_reservations": 0,
            "expected_active_reservations": 0,
            "expected_ledger_rows": 0,
            "expected_manifests": 0,
        },
    )()
    with lock.open("r+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        report = module.run_reconcile(
            tmp_path, root, False, args, stage_i_lock_held=True
        )
    assert report["consistent"] is True


def test_checkpoint_reconcile_child_is_isolated_and_sanitized(monkeypatch, tmp_path):
    module = load_checkpoint_module()
    helper = tmp_path / "stage_i_helper.py"
    helper.write_text("# authenticated fixture helper\n")
    helper.chmod(0o644)
    calls = []

    @contextmanager
    def descriptor(*_args):
        value = os.open(helper, os.O_RDONLY)
        try:
            yield value
        finally:
            os.close(value)

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=json.dumps(
                {
                    "execution_epoch": EPOCH,
                    "consistent": True,
                    "counts": {
                        "transactions": 0,
                        "reservations": 0,
                        "active_reservations": 0,
                        "ledger_rows": 0,
                        "manifests": 0,
                    },
                }
            ),
            stderr="",
        )

    monkeypatch.setattr(module, "stage_i_descriptor", descriptor)
    monkeypatch.setattr(module.subprocess, "run", fake_run)
    args = SimpleNamespace(
        expected_stage_i_sha256=sha256(helper),
        expected_transactions=0,
        expected_reservations=0,
        expected_active_reservations=0,
        expected_ledger_rows=0,
        expected_manifests=0,
    )
    module.run_reconcile(tmp_path, tmp_path, True, args)

    command, kwargs = calls[0]
    assert command[1:4] == ["-I", "-S", "-B"]
    python_descriptor = int(kwargs["executable"].removeprefix("/proc/self/fd/"))
    assert python_descriptor in kwargs["pass_fds"]
    helper_descriptor = next(
        descriptor for descriptor in kwargs["pass_fds"] if descriptor != python_descriptor
    )
    assert kwargs["env"] == module.stage_i_reconcile_environment(
        tmp_path, helper_descriptor, python_descriptor
    )
    assert kwargs["env"][module.STAGE_I_SOURCE_ENV] == str(tmp_path / module.STAGE_I_RELATIVE)
    assert kwargs["env"][module.STAGE_I_REPOSITORY_ROOT_ENV] == str(tmp_path)
    assert kwargs["stdin"] == subprocess.DEVNULL
    assert kwargs["timeout"] == 120


def test_checkpoint_normalizes_restricted_lock_after_interruption(recost_fixture):
    fixture = recost_fixture
    fixture["lock"].write_text("")
    fixture["lock"].chmod(0o600)
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    assert fixture["lock"].stat().st_mode & 0o777 == 0o644


def test_checkpoint_creates_exact_profiles_under_restrictive_umask(recost_fixture):
    fixture = recost_fixture
    promoted = subprocess.run(
        checkpoint_command(fixture, "promote-recost"),
        check=False,
        capture_output=True,
        text=True,
        preexec_fn=lambda: os.umask(0o777),
    )
    assert promoted.returncode == 0, promoted.stderr
    assert fixture["lock"].stat().st_mode & 0o777 == 0o644
    assert fixture["audit"].stat().st_mode & 0o777 == 0o644
    forensic = next(fixture["recost_forensics"].iterdir())
    assert forensic.stat().st_mode & 0o777 == 0o444


@pytest.mark.parametrize(
    "fixture_key",
    ["recost_transactions", "recost_forensics_root", "recost_forensics"],
)
def test_checkpoint_rejects_writable_managed_directory(recost_fixture, fixture_key):
    fixture = recost_fixture
    directory = fixture[fixture_key]
    directory.mkdir(parents=True, exist_ok=True)
    directory.chmod(0o777)
    completed = run_checkpoint(fixture, "promote-recost")
    assert_rejected(completed, "exceeds trusted profile 0755")
    assert fixture["staged"].is_file()
    assert not fixture["canonical"].exists()


@pytest.mark.parametrize("fixture_key", ["root", "accounting", "stage_i_transactions"])
def test_checkpoint_rejects_writable_trust_boundary_directory(recost_fixture,
                                                              fixture_key):
    fixture = recost_fixture
    fixture[fixture_key].chmod(0o777)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "exceeds trusted profile 0755")
    assert fixture["staged"].is_file()
    assert not fixture["canonical"].exists()


def test_checkpoint_rejects_namespace_ambiguity_and_hash_mismatch(recost_fixture):
    ambiguous = recost_fixture["accounting"] / f".{ARTIFACT}.unexpected"
    ambiguous.write_text("ambiguous\n")
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert_rejected(completed, "recost namespace is not staged-only")
    ambiguous.unlink()
    recost_fixture["staged"].write_text("{}\n")
    completed = run_checkpoint(
        recost_fixture,
        "verify-staged-recost",
        artifact_sha256="f" * 64,
    )
    assert_rejected(completed, "recost artifact checksum has changed")


def test_checkpoint_rejects_exact_hidden_twin(recost_fixture):
    hidden = recost_fixture["accounting"] / f".{ARTIFACT}"
    hidden.write_text("ambiguous\n")
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert_rejected(completed, "recost namespace is not staged-only")


def test_checkpoint_rejects_external_evidence_parent_symlink(recost_fixture):
    fixture = recost_fixture
    real = fixture["accounting"] / "real-utilities"
    real.mkdir()
    escaped = real / "generate_test_recost.py"
    escaped.write_bytes(fixture["generator"].read_bytes())
    escaped.chmod(0o755)
    fixture["generator"].unlink()
    fixture["generator"].parent.rmdir()
    fixture["generator"].parent.symlink_to(real, target_is_directory=True)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert completed.returncode == 1
    assert not fixture["canonical"].exists()


def test_checkpoint_rejects_stage_i_transaction_and_bad_artifact_mode(recost_fixture):
    transaction = recost_fixture["stage_i_transactions"] / "unexpected.json"
    transaction.write_text("{}\n")
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert_rejected(completed, "Stage I transaction directory is not empty")
    transaction.unlink()
    recost_fixture["staged"].chmod(0o600)
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert_rejected(completed, "recost artifact mode is 0600, expected 0644")


@pytest.mark.parametrize(
    ("option", "mode", "pattern"),
    [
        ("--expected-artifact-mode", "0666", "artifact mode must not exceed 0644"),
        ("--expected-generator-mode", "0777", "generator mode must not exceed 0755"),
        (
            "--expected-scheduler-mode",
            "0666",
            "scheduler evidence mode must not exceed 0644",
        ),
    ],
)
def test_checkpoint_rejects_unsafe_retained_mode_argument(recost_fixture,
                                                         option,
                                                         mode,
                                                         pattern):
    completed = run_checkpoint(recost_fixture, "verify-staged-recost", option, mode)
    assert completed.returncode == 2
    assert pattern in completed.stderr


def test_checkpoint_rejects_nonzero_expected_active_reservations(recost_fixture):
    recost_fixture["counts"]["active_reservations"] = 1
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert_rejected(completed, "expected active reservation count must be zero")


def test_checkpoint_rejects_lexical_evidence_traversal(recost_fixture):
    completed = run_checkpoint(
        recost_fixture,
        "verify-staged-recost",
        generator_relative_path="../outside-generator.py",
    )
    assert_rejected(completed, "must be a relative path without '..'")


def test_checkpoint_rejects_nonempty_recost_transaction_directory(recost_fixture):
    recost_fixture["recost_transactions"].mkdir()
    transaction = recost_fixture["recost_transactions"] / "unexpected.json"
    transaction.write_text("{}\n")
    completed = run_checkpoint(recost_fixture, "promote-recost")
    assert_rejected(completed, "recost transaction directory is not empty")
    assert recost_fixture["staged"].is_file()
    assert not recost_fixture["canonical"].exists()


def test_checkpoint_rejects_staged_artifact_extra_hardlink(recost_fixture, tmp_path):
    os.link(recost_fixture["staged"], tmp_path / "unexpected-artifact-link.json")
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert_rejected(completed, "recost artifact has 2 links, expected 1")
    assert recost_fixture["staged"].is_file()
    assert not recost_fixture["canonical"].exists()


def test_checkpoint_rejects_scheduler_leaf_symlink(recost_fixture):
    scheduler = recost_fixture["scheduler"]
    retained = scheduler.with_name(f"{scheduler.name}.retained")
    scheduler.rename(retained)
    scheduler.symlink_to(retained.name)
    completed = run_checkpoint(recost_fixture, "verify-staged-recost")
    assert completed.returncode == 1
    assert "symbolic links" in completed.stderr
    assert recost_fixture["staged"].is_file()
    assert not recost_fixture["canonical"].exists()


def test_checkpoint_retains_forensics_and_finalizes_linked_pair(recost_fixture):
    failed = run_checkpoint(
        recost_fixture,
        "promote-recost",
        "--simulate-post-link-failure",
    )
    assert_rejected(failed, "simulated post-link publication failure")
    staged = recost_fixture["staged"].stat()
    canonical = recost_fixture["canonical"].stat()
    assert (staged.st_dev, staged.st_ino) == (canonical.st_dev, canonical.st_ino)
    assert staged.st_nlink == 2
    journals = list(recost_fixture["recost_transactions"].iterdir())
    forensics = list(recost_fixture["recost_forensics"].iterdir())
    assert len(journals) == 1
    assert len(forensics) == 1
    assert sha256(forensics[0]) == sha256(recost_fixture["staged"])
    journal = json.loads(journals[0].read_text())
    assert journal["state"] == "ambiguous-after-link-attempt"
    assert journal["staged_inode_identity"] == {
        "device": staged.st_dev,
        "inode": staged.st_ino,
    }

    finalized = run_checkpoint(recost_fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr
    assert recost_fixture["canonical"].is_file()
    assert recost_fixture["canonical"].stat().st_nlink == 1
    assert not recost_fixture["staged"].exists()
    assert recost_fixture["audit"].is_file()
    assert list(recost_fixture["recost_transactions"].iterdir()) == []
    verified = run_checkpoint(recost_fixture, "verify-promoted-recost")
    assert verified.returncode == 0, verified.stderr


def test_checkpoint_retains_link_return_interruption_for_recovery(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-link-return-failure")
    assert_rejected(failed, "simulated failure immediately after recost publication link")
    staged = fixture["staged"].stat()
    canonical = fixture["canonical"].stat()
    assert (staged.st_dev, staged.st_ino) == (canonical.st_dev, canonical.st_ino)
    assert staged.st_nlink == 2
    journals = list(fixture["recost_transactions"].iterdir())
    assert len(journals) == 1
    journal = json.loads(journals[0].read_text())
    assert journal["state"] == "ambiguous-after-link-attempt"
    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr


def test_checkpoint_recovers_interruption_after_single_link_exchange(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(
        fixture,
        "promote-recost",
        "--simulate-single-link-post-exchange-failure",
    )
    assert_rejected(failed, "simulated interruption after canonical copy publication")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    assert record["state"] == "single-link-copy-prepared"
    assert fixture["canonical"].stat().st_nlink == 1
    replacement = fixture["accounting"] / record["single_link_replacement_name"]
    assert fixture["staged"].stat().st_nlink in {2, 3}
    assert not replacement.exists()

    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr
    assert fixture["canonical"].stat().st_nlink == 1
    assert not fixture["staged"].exists()
    assert not replacement.exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_lustre_recovers_interrupted_canonical_publication(recost_fixture):
    fixture = recost_fixture
    force_checkpoint_renameat2_einval(fixture["checkpoint"])
    failed = run_checkpoint(
        fixture,
        "promote-recost",
        "--simulate-single-link-post-exchange-failure",
    )
    assert_rejected(failed, "simulated interruption after canonical copy publication")

    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr
    verified = run_checkpoint(fixture, "verify-promoted-recost")
    assert verified.returncode == 0, verified.stderr
    assert fixture["canonical"].stat().st_nlink == 1
    assert not fixture["staged"].exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


@pytest.mark.parametrize("partial", [False, True])
def test_checkpoint_recovers_prejournal_single_link_copy(recost_fixture, partial):
    fixture = recost_fixture
    failed = run_checkpoint(
        fixture,
        "promote-recost",
        "--simulate-post-link-failure",
    )
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    replacement = fixture["accounting"] / record["single_link_replacement_name"]
    if partial:
        replacement.write_bytes(b"partial canonical copy")
        replacement.chmod(0o000)
    else:
        module = load_checkpoint_module()
        paths = module.layout(fixture["root"], ARTIFACT)
        args = module.parser().parse_args(
            checkpoint_command(fixture, "finalize-linked-pair")[2:]
        )
        with module.bound_parent_descriptor(
            paths["canonical"], "fixture prejournal canonical copy"
        ) as directory_descriptor:
            module.prepare_single_link_copy(
                directory_descriptor,
                paths,
                args,
                replacement.name,
                expected_staged_identity=module.recovery_staged_identity(record),
            )
    assert "canonical_inode_identity" not in record

    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr
    assert fixture["canonical"].stat().st_nlink == 1
    assert not fixture["staged"].exists()
    assert not replacement.exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_json_write_rejects_detached_parent(tmp_path, monkeypatch):
    module = load_checkpoint_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    target = parent / "state.json"
    target.write_text('{"old": true}\n')
    detached = tmp_path / "accounting-detached"
    original_replacement = module.replace_bound_entry_forward
    raced = False

    def detach_parent(*args, **kwargs):
        nonlocal raced
        result = original_replacement(*args, **kwargs)
        if not raced:
            raced = True
            parent.rename(detached)
            parent.mkdir()
        return result

    monkeypatch.setattr(module, "replace_bound_entry_forward", detach_parent)
    with pytest.raises(ValueError, match="parent path changed during mutation"):
        module.write_json(target, {"new": True})
    assert not target.exists()
    assert json.loads((detached / target.name).read_text()) == {"new": True}


def test_checkpoint_recovery_rejects_changed_authorization_context(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    fixture["staged"].chmod(0o666)
    completed = run_checkpoint(
        fixture,
        "finalize-linked-pair",
        "--expected-artifact-mode",
        "0666",
    )
    assert completed.returncode == 2
    assert "artifact mode must not exceed 0644" in completed.stderr
    assert fixture["staged"].exists()
    assert fixture["canonical"].exists()
    assert len(list(fixture["recost_transactions"].iterdir())) == 1


def test_checkpoint_recovery_rejects_journal_hardlink(recost_fixture, tmp_path):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    os.link(journal, tmp_path / "unexpected-journal-link.json")
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost transaction journal has 2 links, expected 1")
    assert fixture["staged"].exists()
    assert fixture["canonical"].exists()


def test_checkpoint_recovery_rejects_forensic_mode_change(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    forensic = next(fixture["recost_forensics"].iterdir())
    forensic.chmod(0o600)
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost forensic copy mode is 0600, expected 0444")
    assert fixture["staged"].exists()
    assert fixture["canonical"].exists()


def test_checkpoint_recovery_rejects_unsafe_accounting_before_orphan_cleanup(
    recost_fixture,
):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    orphan = (
        fixture["recost_transactions"]
        / f".orphan.json.123.{'0' * 32}.tmp"
    )
    orphan.write_text("partial")
    orphan.chmod(0o000)
    fixture["accounting"].chmod(0o777)
    before = set(fixture["accounting"].iterdir())

    completed = run_checkpoint(fixture, "finalize-linked-pair")

    assert_rejected(completed, "exceeds trusted profile 0755")
    assert orphan.exists()
    assert set(fixture["accounting"].iterdir()) == before


@pytest.mark.parametrize("journal_state", ["preparing", "link-pending"])
def test_checkpoint_retires_interrupted_prepublication_transaction(recost_fixture,
                                                                    journal_state):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["state"] = journal_state
    record.pop("ambiguity_recorded_utc")
    write_json(journal, record)
    fixture["canonical"].unlink()
    if journal_state == "preparing":
        next(fixture["recost_forensics"].iterdir()).chmod(0o400)
    retired = run_checkpoint(fixture, "retire-preparing")
    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert fixture["staged"].stat().st_nlink == 1
    assert not fixture["canonical"].exists()
    assert list(fixture["recost_transactions"].iterdir()) == []
    assert list(fixture["recost_forensics"].iterdir()) == []


def test_checkpoint_retires_prejournal_partial_copy_from_staged_state(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["state"] = "link-pending"
    record.pop("ambiguity_recorded_utc")
    write_json(journal, record)
    fixture["canonical"].unlink()
    replacement = fixture["accounting"] / record["single_link_replacement_name"]
    replacement.write_bytes(b"partial canonical copy")
    replacement.chmod(0o000)

    retired = run_checkpoint(fixture, "retire-preparing")
    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert fixture["staged"].stat().st_nlink == 1
    assert not replacement.exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


@pytest.mark.parametrize("mode", [0o000, 0o600])
def test_checkpoint_retires_orphan_journal_temporary(recost_fixture, mode):
    fixture = recost_fixture
    fixture["recost_transactions"].mkdir()
    fixture["recost_forensics"].mkdir(parents=True)
    temporary = (
        fixture["recost_transactions"]
        / f".orphan.json.123.{'0' * 32}.tmp"
    )
    temporary.write_text("partial")
    temporary.chmod(mode)
    retired = run_checkpoint(fixture, "retire-preparing")
    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert list(fixture["recost_transactions"].iterdir()) == []


@pytest.mark.parametrize("legacy_name", [False, True])
def test_checkpoint_retires_mode_zero_forensic_copy_temporary(
    recost_fixture, legacy_name,
):
    fixture = recost_fixture
    fixture["recost_transactions"].mkdir()
    fixture["recost_forensics"].mkdir(parents=True)
    if legacy_name:
        name = f".cgl-checkpoint-forensic-{'0' * 32}"
    else:
        destination = fixture["recost_forensics"] / "interrupted.forensic"
        name = (
            ".cgl-checkpoint-forensic-"
            f"{hashlib.sha256(destination.name.encode()).hexdigest()}.tmp"
        )
    temporary = fixture["recost_forensics"] / name
    temporary.write_bytes(b"partial forensic copy")
    temporary.chmod(0o000)

    retired = run_checkpoint(fixture, "retire-preparing")

    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert list(fixture["recost_transactions"].iterdir()) == []
    assert list(fixture["recost_forensics"].iterdir()) == []


def test_checkpoint_retires_link_pending_after_forensic_cleanup_interruption(
    recost_fixture,
):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["state"] = "link-pending"
    record.pop("ambiguity_recorded_utc")
    write_json(journal, record)
    fixture["canonical"].unlink()
    next(fixture["recost_forensics"].iterdir()).unlink()
    retired = run_checkpoint(fixture, "retire-preparing")
    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_rejects_forensic_directory_symlink_before_finalize(
    recost_fixture,
):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    retained = fixture["recost_forensics"].with_name(
        f"{fixture['recost_forensics'].name}.retained"
    )
    fixture["recost_forensics"].rename(retained)
    fixture["recost_forensics"].symlink_to(retained.name, target_is_directory=True)
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost forensic directory is not a directory")
    staged = fixture["staged"].stat()
    canonical = fixture["canonical"].stat()
    assert (staged.st_dev, staged.st_ino) == (canonical.st_dev, canonical.st_ino)
    assert staged.st_nlink == 2
    assert len(list(fixture["recost_transactions"].iterdir())) == 1


def test_checkpoint_rejects_extra_forensic_entry(recost_fixture):
    fixture = recost_fixture
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    extra = fixture["recost_forensics"] / "unexpected"
    extra.write_text("unexpected\n")
    completed = run_checkpoint(fixture, "verify-promoted-recost")
    assert_rejected(completed, "recost forensic directory entries differ")


def test_checkpoint_rechecks_queue_directly_before_recovery_unlink(recost_fixture,
                                                                   tmp_path):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    late_queue = tmp_path / "late-unlink-squeue.txt"
    late_queue.write_text("67890|cgl_late_root_writer|RUNNING\n")
    completed = run_checkpoint(
        fixture,
        "finalize-linked-pair",
        "--pre-unlink-squeue-file",
        str(late_queue),
    )
    assert_rejected(completed, "another CGL job is queued")
    assert fixture["staged"].is_file()
    assert fixture["canonical"].is_file()
    assert len(list(fixture["recost_transactions"].iterdir())) == 1


def test_checkpoint_recovery_revalidates_staged_name_before_unlink(recost_fixture):
    fixture = recost_fixture
    module = load_checkpoint_module()
    paths = module.layout(fixture["root"], ARTIFACT)
    args = module.parser().parse_args(
        checkpoint_command(fixture, "finalize-linked-pair")[2:]
    )
    os.link(fixture["staged"], fixture["canonical"])
    staged_identity = (
        fixture["staged"].stat().st_dev,
        fixture["staged"].stat().st_ino,
    )
    replacement = module.single_link_replacement_name("fixture-staged-race")

    def replace_staged_name():
        retained = fixture["staged"].read_bytes()
        fixture["staged"].unlink()
        fixture["staged"].write_bytes(retained)
        fixture["staged"].chmod(0o644)

    with module.bound_parent_descriptor(
        paths["canonical"], "fixture staged-name race"
    ) as directory_descriptor:
        canonical_identity = module.prepare_single_link_copy(
            directory_descriptor,
            paths,
            args,
            replacement,
            expected_staged_identity=staged_identity,
        )
        replace_staged_name()
        with pytest.raises(
            ValueError,
            match="staged recost publication link does not select its journaled inode",
        ):
            module.complete_single_link_transition(
                directory_descriptor,
                paths,
                args,
                replacement,
                staged_identity,
                canonical_identity,
            )
    assert fixture["staged"].is_file()
    assert fixture["canonical"].is_file()
    assert fixture["staged"].stat().st_ino != fixture["canonical"].stat().st_ino


def test_checkpoint_recovery_rechecks_accounting_profile_before_unlink(recost_fixture):
    fixture = recost_fixture
    module = load_checkpoint_module()
    paths = module.layout(fixture["root"], ARTIFACT)
    args = module.parser().parse_args(
        checkpoint_command(fixture, "finalize-linked-pair")[2:]
    )
    os.link(fixture["staged"], fixture["canonical"])
    staged_identity = (
        fixture["staged"].stat().st_dev,
        fixture["staged"].stat().st_ino,
    )
    replacement = module.single_link_replacement_name("fixture-accounting-race")
    with module.bound_parent_descriptor(
        paths["canonical"], "fixture accounting-profile race"
    ) as directory_descriptor:
        canonical_identity = module.prepare_single_link_copy(
            directory_descriptor,
            paths,
            args,
            replacement,
            expected_staged_identity=staged_identity,
        )
    record = {
        "staged_inode_identity": module.identity_binding(staged_identity),
        "canonical_inode_identity": module.identity_binding(canonical_identity),
        "single_link_replacement_name": replacement,
    }
    fixture["accounting"].chmod(0o777)
    with pytest.raises(ValueError, match="exceeds trusted profile 0755"):
        module.complete_journaled_single_link_transition(
            paths,
            fixture["root"],
            True,
            str(fixture["queue"]),
            args,
            record,
        )
    assert fixture["staged"].is_file()
    assert fixture["canonical"].is_file()


def test_checkpoint_single_link_transition_rejects_accounting_mode_drift(
    recost_fixture, monkeypatch,
):
    fixture = recost_fixture
    module = load_checkpoint_module()
    paths = module.layout(fixture["root"], ARTIFACT)
    args = module.parser().parse_args(
        checkpoint_command(fixture, "finalize-linked-pair")[2:]
    )
    os.link(fixture["staged"], fixture["canonical"])
    staged_identity = (
        fixture["staged"].stat().st_dev,
        fixture["staged"].stat().st_ino,
    )
    replacement = module.single_link_replacement_name("fixture-mode-drift")
    original_publish = module.rename_bound_noreplace

    def weaken_after_publication(*publish_args, **publish_kwargs):
        result = original_publish(*publish_args, **publish_kwargs)
        fixture["accounting"].chmod(0o777)
        return result

    monkeypatch.setattr(module, "rename_bound_noreplace", weaken_after_publication)
    with pytest.raises(ValueError, match="parent path changed during mutation"):
        with module.bound_parent_descriptor(
            paths["canonical"], "fixture accounting mode drift"
        ) as directory_descriptor:
            canonical_identity = module.prepare_single_link_copy(
                directory_descriptor,
                paths,
                args,
                replacement,
                expected_staged_identity=staged_identity,
            )
            module.complete_single_link_transition(
                directory_descriptor,
                paths,
                args,
                replacement,
                staged_identity,
                canonical_identity,
            )


def test_checkpoint_recovery_rejects_non_utc_journal_timestamp(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["created_utc"] = "2026-06-02T01:00:00+05:00"
    write_json(journal, record)
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost journal creation timestamp must use UTC")


def test_checkpoint_finalizes_after_orphan_audit_temporary(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    temporary = (
        fixture["audit"].parent
        / f".{fixture['audit'].name}.123.{'0' * 32}.tmp"
    )
    temporary.write_text("partial")
    temporary.chmod(0o600)
    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr
    assert not temporary.exists()
    assert fixture["audit"].is_file()


def test_checkpoint_finalizes_after_mode_zero_orphan_audit_temporary(recost_fixture):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    temporary = (
        fixture["audit"].parent
        / f".{fixture['audit'].name}.123.{'0' * 32}.tmp"
    )
    temporary.write_text("partial")
    temporary.chmod(0o000)

    finalized = run_checkpoint(fixture, "finalize-linked-pair")

    assert finalized.returncode == 0, finalized.stderr
    assert not temporary.exists()
    assert fixture["audit"].is_file()


@pytest.mark.parametrize(
    ("fixture_key", "mode_option"),
    [
        ("generator", "--expected-generator-mode"),
        ("scheduler", "--expected-scheduler-mode"),
    ],
)
def test_checkpoint_promoted_audit_rejects_evidence_mode_drift(recost_fixture,
                                                               fixture_key,
                                                               mode_option):
    fixture = recost_fixture
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    fixture[fixture_key].chmod(0o600)
    completed = run_checkpoint(
        fixture,
        "verify-promoted-recost",
        mode_option,
        "0600",
    )
    assert completed.returncode == 1
    assert "binding differs" in completed.stderr


@pytest.mark.parametrize("journal_state", ["ambiguous-after-link-attempt", "link-pending"])
def test_checkpoint_finalizes_canonical_only_recovery(recost_fixture, journal_state):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["state"] = journal_state
    if journal_state == "link-pending":
        record.pop("ambiguity_recorded_utc")
    write_json(journal, record)
    fixture["staged"].unlink()
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert completed.returncode == 0, completed.stderr
    assert fixture["canonical"].is_file()
    assert not fixture["staged"].exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_rejects_replacement_canonical_inode_during_recovery(
    recost_fixture,
):
    fixture = recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    retained = fixture["canonical"].read_bytes()
    fixture["staged"].unlink()
    fixture["canonical"].unlink()
    fixture["canonical"].write_bytes(retained)
    fixture["canonical"].chmod(0o644)
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost artifact inode identity has changed")
    assert not fixture["audit"].exists()


def test_checkpoint_rejects_forged_descriptor_source_profile(recost_fixture, tmp_path):
    fixture = recost_fixture
    root = tmp_path / "forged-source"
    source = root / "scripts/frontier/cgl_lf_stage_i_checkpoint.py"
    source.parent.mkdir(parents=True)
    source.write_bytes(CHECKPOINT.read_bytes())
    source.chmod(0o644)
    descriptor = os.memfd_create("forged-cgl-checkpoint")
    python_descriptor = os.open(Path("/proc/self/exe").resolve(), os.O_RDONLY)
    try:
        os.write(descriptor, CHECKPOINT.read_bytes())
        os.lseek(descriptor, 0, os.SEEK_SET)
        environment = dict(os.environ)
        environment["_CGL_LF_RECOST_UTILITY_DESCRIPTOR"] = str(descriptor)
        environment["_CGL_LF_RECOST_UTILITY_PYTHON_DESCRIPTOR"] = str(python_descriptor)
        environment["_CGL_LF_RECOST_UTILITY_SOURCE"] = str(source)
        environment["_CGL_LF_RECOST_REPOSITORY_ROOT"] = str(root)
        completed = subprocess.run(
            [
                sys.executable,
                "-I",
                f"/proc/self/fd/{descriptor}",
                *checkpoint_command(fixture, "verify-staged-recost")[2:],
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
            pass_fds=(descriptor, python_descriptor),
        )
    finally:
        os.close(descriptor)
        os.close(python_descriptor)
    assert_rejected(completed, "retained utility mode is 0644, expected 0755")


def test_checkpoint_promotes_and_audits_bounded_wave(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    verified = run_checkpoint(fixture, "verify-staged-recost")
    assert verified.returncode == 0, verified.stderr
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    audit = json.loads(fixture["audit"].read_text())
    assert audit["authorized_bounded_wave"] == fixture["authorization"]
    assert "authorized_sole_next_segment_profile" not in audit
    assert len(audit["scheduler_evidence"]) == 3
    assert audit["source_bundle"]["sha256"] == sha256(fixture["source_bundle"])
    assert audit["source_bundle"]["verified_revisions"] == fixture["source_bundle_revisions"]
    assert audit["generator"]["revision"] == fixture["generator_revision"]
    assert audit["stage_i_helper"]["revision"] == fixture["stage_i_revision"]
    context = audit["generalized_publication_context"]
    assert context["artifact"]["basename"] == V2_ARTIFACT
    assert context["checkpoint"] == "F-114"
    assert context["request"]["sha256"] == fixture["request_sha256"]
    assert context["ledger_tail_job_ids"] == ["12345", "23456", "4766856"]
    assert context["controller_enforcement"]["bounded_wave_authorizing"] is False
    assert audit["forensic_copy"]["generalized_vectors"] == "exact-artifact-payload"
    forensic = Path(audit["forensic_copy"]["path"])
    assert forensic.read_bytes() == fixture["canonical"].read_bytes()
    assert list(fixture["recost_transactions"].iterdir()) == []
    completed = run_checkpoint(fixture, "verify-promoted-recost")
    assert completed.returncode == 0, completed.stderr


def test_checkpoint_finalizes_bounded_wave_linked_pair(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    assert record["authorized_bounded_wave"] == fixture["authorization"]
    assert record["bounded_wave_scheduler_evidence"] == fixture["scheduler_evidence"]
    assert record["generalized_publication_context"]["checkpoint"] == "F-114"
    assert (
        record["generalized_publication_context"]["inputs"]
        == json.loads(fixture["request"].read_text())["inputs"]
    )
    assert record["forensic_generalized_vectors"] == "exact-artifact-payload"
    assert "authorized_sole_next_segment_profile" not in record
    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr
    assert fixture["canonical"].is_file()
    assert not fixture["staged"].exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_v2_post_link_recovery_survives_request_expiry(
    bounded_recost_fixture,
    monkeypatch,
):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")

    module = load_checkpoint_module()
    paths = module.layout(fixture["root"], fixture["artifact_name"])
    args = module.parser().parse_args(
        checkpoint_command(fixture, "finalize-linked-pair")[2:]
    )

    class FutureDateTime(datetime):
        @classmethod
        def now(cls, tz=None):
            retained = datetime.now(timezone.utc) + timedelta(days=2)
            return retained if tz is not None else retained.replace(tzinfo=None)

    monkeypatch.setattr(module, "datetime", FutureDateTime)
    monkeypatch.setattr(module, "initial_source_path", lambda: fixture["checkpoint"])
    monkeypatch.setattr(module, "repository_root", lambda _source: fixture["repository"])
    module.finalize_linked_pair(
        paths,
        fixture["repository"],
        fixture["root"],
        True,
        args,
    )
    assert fixture["canonical"].is_file()
    assert not fixture["staged"].exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_v2_post_link_recovery_survives_unrelated_head_movement(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    unrelated = fixture["repository"] / "unrelated.txt"
    unrelated.write_text("unrelated later repository change\n")
    subprocess.run(["git", "add", unrelated.name], cwd=fixture["repository"], check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=CGL fixture",
            "-c",
            "user.email=cgl-fixture@example.invalid",
            "commit",
            "-q",
            "-m",
            "Advance unrelated fixture HEAD",
        ],
        cwd=fixture["repository"],
        check=True,
    )
    finalized = run_checkpoint(fixture, "finalize-linked-pair")
    assert finalized.returncode == 0, finalized.stderr


def make_v2_transaction_provably_prelink(fixture) -> None:
    """Convert one simulated post-link interruption to exact staged-only state."""

    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["state"] = "link-pending"
    record.pop("ambiguity_recorded_utc")
    write_json(journal, record)
    fixture["canonical"].unlink()


def test_checkpoint_v2_prelink_retirement_survives_request_expiry(
    bounded_recost_fixture,
    monkeypatch,
):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    make_v2_transaction_provably_prelink(fixture)

    module = load_checkpoint_module()
    paths = module.layout(fixture["root"], fixture["artifact_name"])
    args = module.parser().parse_args(
        checkpoint_command(fixture, "retire-preparing")[2:]
    )

    class FutureDateTime(datetime):
        @classmethod
        def now(cls, tz=None):
            retained = datetime.now(timezone.utc) + timedelta(days=2)
            return retained if tz is not None else retained.replace(tzinfo=None)

    monkeypatch.setattr(module, "datetime", FutureDateTime)
    monkeypatch.setattr(module, "initial_source_path", lambda: fixture["checkpoint"])
    monkeypatch.setattr(module, "repository_root", lambda _source: fixture["repository"])
    module.retire_preparing(
        paths,
        fixture["repository"],
        fixture["root"],
        True,
        args,
    )
    assert fixture["staged"].is_file()
    assert not fixture["canonical"].exists()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_v2_prelink_retirement_survives_unrelated_head_movement(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    make_v2_transaction_provably_prelink(fixture)
    unrelated = fixture["repository"] / "unrelated-retirement.txt"
    unrelated.write_text("unrelated later repository change\n")
    subprocess.run(["git", "add", unrelated.name], cwd=fixture["repository"], check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=CGL fixture",
            "-c",
            "user.email=cgl-fixture@example.invalid",
            "commit",
            "-q",
            "-m",
            "Advance unrelated retirement fixture HEAD",
        ],
        cwd=fixture["repository"],
        check=True,
    )
    retired = run_checkpoint(fixture, "retire-preparing")
    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert list(fixture["recost_transactions"].iterdir()) == []


def test_checkpoint_v2_pre_link_still_rejects_unrelated_head_movement(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    unrelated = fixture["repository"] / "unrelated.txt"
    unrelated.write_text("unrelated prepublication repository change\n")
    subprocess.run(["git", "add", unrelated.name], cwd=fixture["repository"], check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=CGL fixture",
            "-c",
            "user.email=cgl-fixture@example.invalid",
            "commit",
            "-q",
            "-m",
            "Advance prepublication fixture HEAD",
        ],
        cwd=fixture["repository"],
        check=True,
    )
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(
        completed,
        "Stage I recost generator revision differs from repository HEAD",
    )


def test_checkpoint_v2_pre_link_still_rejects_expired_request(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    now = datetime.now(timezone.utc).replace(microsecond=0)
    request = json.loads(fixture["request"].read_text())
    request["generated_utc"] = (now - timedelta(hours=2)).isoformat()
    request["expires_utc"] = (now - timedelta(hours=1)).isoformat()
    write_json(fixture["request"], request)
    fixture["request_sha256"] = sha256(fixture["request"])

    def bind_expired_request(payload):
        payload["generated_utc"] = request["generated_utc"]
        payload["expires_utc"] = request["expires_utc"]
        payload["provenance"]["request_sha256"] = fixture["request_sha256"]

    mutate_bounded_artifact(fixture, bind_expired_request)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 recost request has expired")


def test_checkpoint_bounded_recovery_rejects_changed_packet_context(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    fixture["authorization"]["authorized_next_profiles"][0][
        "estimated_storage_bytes"
    ] += 1
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost artifact V2 authorization differs")
    assert fixture["staged"].is_file()
    assert fixture["canonical"].is_file()


def test_checkpoint_rejects_bounded_wave_duplicate_lane(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def duplicate_lane(authorization):
        authorization["authorized_next_profiles"][1]["case_id"] = "R04"

    mutate_bounded_authorization(fixture, duplicate_lane)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 profiles duplicate lane R04")


def test_checkpoint_rejects_bounded_wave_above_four_lanes(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def add_lanes(authorization):
        profiles = authorization["authorized_next_profiles"]
        for case_id in ("R05", "R06", "R07"):
            retained = json.loads(json.dumps(profiles[0]))
            retained["case_id"] = case_id
            retained["nodes"] = 1
            profiles.append(retained)
        profiles.sort(key=lambda item: (item["case_id"], item["segment"]))
        authorization["bounded_concurrency"]["max_wave_nodes"] = 9

    mutate_bounded_authorization(fixture, add_lanes)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded-wave authorization requires 2-4 profiles")


def test_checkpoint_rejects_bounded_wave_above_ten_nodes(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def exceed_nodes(authorization):
        authorization["authorized_next_profiles"][0]["nodes"] = 6
        authorization["authorized_next_profiles"][1]["nodes"] = 5
        authorization["bounded_concurrency"]["max_wave_nodes"] = 11

    mutate_bounded_authorization(fixture, exceed_nodes)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 authorization exceeds the 10-node ceiling")


def test_checkpoint_rejects_r17_in_bounded_wave(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def insert_r17(authorization):
        authorization["authorized_next_profiles"][0]["nodes"] = 1
        authorization["authorized_next_profiles"][1]["case_id"] = "R17"
        authorization["authorized_next_profiles"][1]["nodes"] = 8
        authorization["bounded_concurrency"]["max_wave_nodes"] = 9

    mutate_bounded_authorization(fixture, insert_r17)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "R17 is forbidden in bounded-wave mode")


def test_checkpoint_rejects_bounded_wave_without_r17_last_binding(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def remove_last_policy(authorization):
        authorization["bounded_concurrency"]["r17_exclusive_and_last"] = False

    mutate_bounded_authorization(fixture, remove_last_policy)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "does not preserve R17 exclusive/last")


def test_checkpoint_rejects_bounded_scheduler_row_mismatch(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def change_elapsed(evidence):
        evidence[0]["elapsed_seconds"] += 1

    mutate_bounded_scheduler_evidence(fixture, change_elapsed)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded scheduler evidence differs for job 12345")


def test_checkpoint_scheduler_schema_preserves_noncompleted_terminal_identity():
    module = load_checkpoint_module()
    evidence = {
        "path": "accounting/12345.stage_i.sacct.txt",
        "sha256": "1" * 64,
        "job_id": "12345",
        "job_name": "cgl_mks24_E03_forcing_policy_R03_s00_rankio_t0_t0p5",
        "state": "FAILED",
        "exit_code": "1:0",
        "nodes": 1,
        "elapsed_seconds": 60,
        "submitted_utc": "2026-06-04T01:00:00+00:00",
        "completed_utc": "2026-06-04T01:01:00+00:00",
    }
    assert module.validate_bounded_scheduler_evidence([evidence]) == [evidence]
    evidence["state"] = "RUNNING"
    with pytest.raises(ValueError, match="state is not terminal"):
        module.validate_bounded_scheduler_evidence([evidence])


def test_checkpoint_scheduler_schema_preserves_authenticated_barrier_order():
    module = load_checkpoint_module()
    first = {
        "path": "accounting/23456.stage_i.sacct.txt",
        "sha256": "1" * 64,
        "job_id": "23456",
        "job_name": "cgl_mks24_E03_forcing_policy_R12_s00_rankio_t0_t0p25",
        "state": "COMPLETED",
        "exit_code": "0:0",
        "nodes": 4,
        "elapsed_seconds": 60,
        "submitted_utc": "2026-06-04T01:00:00+00:00",
        "completed_utc": "2026-06-04T01:01:00+00:00",
    }
    second = {
        "path": "accounting/12345.stage_i.sacct.txt",
        "sha256": "2" * 64,
        "job_id": "12345",
        "job_name": "cgl_mks24_E03_forcing_policy_R03_s02_rankio_t0p312823_t0p5",
        "state": "COMPLETED",
        "exit_code": "0:0",
        "nodes": 1,
        "elapsed_seconds": 60,
        "submitted_utc": "2026-06-04T01:00:00+00:00",
        "completed_utc": "2026-06-04T01:01:00+00:00",
    }
    scheduler = module.validate_bounded_scheduler_evidence([first, second])
    barrier = [
        {
            "case_id": "R12",
            "segment": "s00_rankio_t0_t0p25",
            "job_id": "23456",
            "result": "clean_partial",
        },
        {
            "case_id": "R03",
            "segment": "s02_rankio_t0p312823_t0p5",
            "job_id": "12345",
            "result": "accepted",
        },
    ]
    module.validate_bounded_barrier(barrier, scheduler)
    with pytest.raises(ValueError, match="job order differs"):
        module.validate_bounded_barrier(list(reversed(barrier)), scheduler)


def test_checkpoint_allows_exact_fresh_r12_after_historical_clean_partial(
    bounded_recost_fixture,
):
    module = load_checkpoint_module()
    fixture = bounded_recost_fixture
    authorization = json.loads(json.dumps(fixture["authorization"]))
    profile = authorization["authorized_next_profiles"][1]
    profile.update(
        {
            "segment": "s01_rankio_t0_t0p12",
            "nodes": 4,
            "parent_job_id": None,
            "parent_result": None,
            "parent_segment": None,
            "restart_file": None,
            "restart_file_sha256": None,
            "restart_time": None,
            "time_tlim_target": 0.12,
        }
    )
    authorization["bounded_concurrency"]["max_wave_nodes"] = 8
    args = SimpleNamespace(
        source_bundle_relative_path=str(
            fixture["source_bundle"].relative_to(fixture["root"])
        ),
        expected_source_bundle_sha256=fixture["source_bundle_sha256"],
    )

    retained = module.validate_bounded_wave_authorization(
        authorization, fixture["root"], args
    )
    module.validate_fresh_r12_rerun_evidence(
        retained["authorized_next_profiles"],
        [module.HISTORICAL_R12_CLEAN_PARTIAL],
    )


@pytest.mark.parametrize(
    ("mutation", "error"),
    [
        ({"nodes": 2}, "fresh R12 rerun profile differs"),
        ({"time_tlim_target": 0.5}, "fresh R12 rerun profile differs"),
        ({"parent_job_id": "4766856"}, "fresh R12 rerun profile differs"),
    ],
)
def test_checkpoint_rejects_drifted_fresh_r12_decision(
    bounded_recost_fixture, mutation, error,
):
    module = load_checkpoint_module()
    profile = json.loads(
        json.dumps(bounded_recost_fixture["authorization"]["authorized_next_profiles"][1])
    )
    profile.update(
        {
            "segment": "s01_rankio_t0_t0p12",
            "nodes": 4,
            "parent_job_id": None,
            "parent_result": None,
            "parent_segment": None,
            "restart_file": None,
            "restart_file_sha256": None,
            "restart_time": None,
            "time_tlim_target": 0.12,
            **mutation,
        }
    )
    with pytest.raises(ValueError, match=error):
        module.validate_fresh_r12_rerun_evidence(
            [profile],
            [module.HISTORICAL_R12_CLEAN_PARTIAL],
        )


def test_checkpoint_rejects_fresh_r12_without_exact_historical_inventory(
    bounded_recost_fixture,
):
    module = load_checkpoint_module()
    profile = json.loads(
        json.dumps(bounded_recost_fixture["authorization"]["authorized_next_profiles"][1])
    )
    profile.update(
        {
            "segment": "s01_rankio_t0_t0p12",
            "nodes": 4,
            "parent_job_id": None,
            "parent_result": None,
            "parent_segment": None,
            "restart_file": None,
            "restart_file_sha256": None,
            "restart_time": None,
            "time_tlim_target": 0.12,
        }
    )
    forged = dict(module.HISTORICAL_R12_CLEAN_PARTIAL)
    forged["result"] = "accepted"
    with pytest.raises(ValueError, match="exact historical s00/4766856 clean_partial"):
        module.validate_fresh_r12_rerun_evidence([profile], [forged])


def test_checkpoint_rejects_historical_r12_clean_partial_as_parent(
    bounded_recost_fixture,
):
    module = load_checkpoint_module()
    fixture = bounded_recost_fixture
    authorization = json.loads(json.dumps(fixture["authorization"]))
    profile = authorization["authorized_next_profiles"][1]
    profile.update(
        {
            "segment": "s01_rankio_t0p137193_t0p25",
            "nodes": 4,
            "parent_job_id": "4766856",
            "parent_result": "clean_partial",
            "parent_segment": "s00_rankio_t0_t0p25",
            "restart_file": str(fixture["root"] / "runs/R12/restart.rst"),
            "restart_file_sha256": "1" * 64,
            "restart_time": 0.1371931229426507,
            "time_tlim_target": 0.25,
        }
    )
    authorization["bounded_concurrency"]["max_wave_nodes"] = 8
    args = SimpleNamespace(
        source_bundle_relative_path=str(
            fixture["source_bundle"].relative_to(fixture["root"])
        ),
        expected_source_bundle_sha256=fixture["source_bundle_sha256"],
    )
    with pytest.raises(ValueError, match="historical R12 s00/4766856.*non-authorizing"):
        module.validate_bounded_wave_authorization(authorization, fixture["root"], args)


def test_checkpoint_rejects_historical_r12_s00_as_new_profile(
    bounded_recost_fixture,
):
    module = load_checkpoint_module()
    fixture = bounded_recost_fixture
    authorization = json.loads(json.dumps(fixture["authorization"]))
    authorization["authorized_next_profiles"][1]["segment"] = "s00_rankio_t0_t0p25"
    args = SimpleNamespace(
        source_bundle_relative_path=str(
            fixture["source_bundle"].relative_to(fixture["root"])
        ),
        expected_source_bundle_sha256=fixture["source_bundle_sha256"],
    )
    with pytest.raises(ValueError, match="inventory-only"):
        module.validate_bounded_wave_authorization(authorization, fixture["root"], args)


def test_checkpoint_rejects_unattached_private_reexecution_environment(
    monkeypatch,
    tmp_path,
):
    module = load_checkpoint_module()
    forged_source = tmp_path / "forged-checkpoint.py"
    forged_source.write_text("raise RuntimeError('forged checkpoint executed')\n")
    forged_root = tmp_path / "forged-repository"
    monkeypatch.delenv(module.SELF_DESCRIPTOR_ENV, raising=False)
    monkeypatch.setenv(module.SELF_SOURCE_ENV, str(forged_source))
    monkeypatch.setenv(module.ROOT_DIR_ENV, str(forged_root))

    with pytest.raises(ValueError, match="forbidden without an authenticated descriptor"):
        module.authenticate_self(sha256(CHECKPOINT))


def test_checkpoint_git_environment_isolated_from_caller_configuration(monkeypatch):
    module = load_checkpoint_module()
    monkeypatch.setenv("PATH", "/tmp/forged-path")
    monkeypatch.setenv("LD_PRELOAD", "/tmp/forged-loader.so")
    monkeypatch.setenv("PYTHONPATH", "/tmp/forged-python")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/tmp/forged-global-config")
    monkeypatch.setenv("GIT_CONFIG_SYSTEM", "/tmp/forged-system-config")
    monkeypatch.setenv("GIT_OBJECT_DIRECTORY", "/tmp/forged-objects")
    environment = module.hardened_git_environment()

    assert environment["GIT_CONFIG_GLOBAL"] == "/dev/null"
    assert environment["GIT_CONFIG_NOSYSTEM"] == "1"
    assert environment["GIT_CONFIG_SYSTEM"] == "/dev/null"
    assert environment["GIT_EXEC_PATH"] == str(module.GIT_EXEC_PATH)
    assert environment["GIT_OPTIONAL_LOCKS"] == "0"
    assert environment["LC_ALL"] == "C"
    assert environment["PATH"] == module.TRUSTED_SYSTEM_PATH
    assert "LD_PRELOAD" not in environment
    assert "PYTHONPATH" not in environment
    assert "GIT_OBJECT_DIRECTORY" not in environment


def test_checkpoint_reexec_environment_strips_private_loader_and_interpreter_state(
    monkeypatch,
):
    module = load_checkpoint_module()
    for name, value in {
        module.SELF_DESCRIPTOR_ENV: "99",
        module.SELF_SOURCE_ENV: "/tmp/forged-source",
        module.ROOT_DIR_ENV: "/tmp/forged-repository",
        "PATH": "/tmp/forged-path",
        "LD_PRELOAD": "/tmp/forged-loader.so",
        "PYTHONPATH": "/tmp/forged-python",
        "GIT_OBJECT_DIRECTORY": "/tmp/forged-objects",
    }.items():
        monkeypatch.setenv(name, value)

    environment = module.reexec_environment(7, 8, CHECKPOINT, REPOSITORY)

    assert environment[module.SELF_DESCRIPTOR_ENV] == "7"
    assert environment[module.PYTHON_DESCRIPTOR_ENV] == "8"
    assert environment[module.SELF_SOURCE_ENV] == str(CHECKPOINT)
    assert environment[module.ROOT_DIR_ENV] == str(REPOSITORY)
    assert environment["PATH"] == module.TRUSTED_SYSTEM_PATH
    assert environment["PYTHONDONTWRITEBYTECODE"] == "1"
    assert environment["HOME"] == "/nonexistent"
    assert environment["XDG_CONFIG_HOME"] == "/nonexistent"
    assert set(environment) == {
        module.SELF_DESCRIPTOR_ENV,
        module.PYTHON_DESCRIPTOR_ENV,
        module.SELF_SOURCE_ENV,
        module.ROOT_DIR_ENV,
        "HOME",
        "LC_ALL",
        "PATH",
        "PYTHONDONTWRITEBYTECODE",
        "XDG_CONFIG_HOME",
    }
    assert "LD_PRELOAD" not in environment
    assert "PYTHONPATH" not in environment
    assert "GIT_OBJECT_DIRECTORY" not in environment


def test_authenticated_checkpoint_descriptor_controls_private_metadata(monkeypatch):
    module = load_checkpoint_module()
    descriptor = os.open(CHECKPOINT, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(), os.O_RDONLY)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(CHECKPOINT))
        monkeypatch.setenv(module.ROOT_DIR_ENV, str(REPOSITORY))
        monkeypatch.setattr(module.sys, "flags", type("Flags", (), {"isolated": 1})())

        source, repository = module.authenticate_self(sha256(CHECKPOINT))
    finally:
        os.close(descriptor)
        os.close(python_descriptor)

    assert source == CHECKPOINT
    assert repository == REPOSITORY


def test_authenticated_checkpoint_rejects_nonisolated_python(monkeypatch):
    module = load_checkpoint_module()
    descriptor = os.open(CHECKPOINT, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(), os.O_RDONLY)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(CHECKPOINT))
        monkeypatch.setenv(module.ROOT_DIR_ENV, str(REPOSITORY))
        monkeypatch.setattr(module.sys, "flags", type("Flags", (), {"isolated": 0})())

        with pytest.raises(ValueError, match="interpreter is not isolated"):
            module.authenticate_self(sha256(CHECKPOINT))
    finally:
        os.close(descriptor)
        os.close(python_descriptor)


def test_checkpoint_rejects_different_root_executable_as_python_descriptor(monkeypatch):
    module = load_checkpoint_module()
    descriptor = os.open(module.SQUEUE, os.O_RDONLY)
    try:
        monkeypatch.setenv(module.PYTHON_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setattr(module.sys, "flags", type("Flags", (), {"isolated": 1})())
        with pytest.raises(ValueError, match="is not this interpreter"):
            module.require_authenticated_python_descriptor()
    finally:
        os.close(descriptor)


def test_checkpoint_git_execution_is_descriptor_bound(monkeypatch):
    module = load_checkpoint_module()
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    module.git_run(REPOSITORY, ["diff", "--quiet", "--", CHECKPOINT.name])

    assert len(calls) == 1
    command, kwargs = calls[0]
    descriptor = int(str(kwargs["executable"]).removeprefix("/proc/self/fd/"))
    assert descriptor in kwargs["pass_fds"]
    assert command[:4] == [
        str(module.GIT), "--no-replace-objects", "-C", str(REPOSITORY)
    ]
    assert "--no-ext-diff" in command
    assert "--no-textconv" in command
    assert f"core.hooksPath={os.devnull}" in command
    assert kwargs["stdin"] == subprocess.DEVNULL
    assert kwargs["env"] == module.hardened_git_environment()


def test_checkpoint_canonical_actions_require_canonical_repository(tmp_path):
    module = load_checkpoint_module()
    with pytest.raises(ValueError, match="canonical use requires"):
        module.require_canonical_repository(tmp_path, False)
    module.require_canonical_repository(module.CANONICAL_REPOSITORY_ROOT, False)
    module.require_canonical_repository(tmp_path, True)


def test_checkpoint_rejects_bounded_controller_walltime_maximum_drift(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def change_maximum(authorization):
        authorization["authorized_next_profiles"][0][
            "controller_walltime_max_seconds"
        ] = 7199

    mutate_bounded_authorization(fixture, change_maximum)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 profile 0 controller walltime maximum differs")


def test_checkpoint_rejects_bounded_scheduler_timestamp_drift(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def change_timestamp(evidence):
        evidence[0]["completed_utc"] = "2026-06-04T01:40:00+00:00"

    mutate_bounded_scheduler_evidence(fixture, change_timestamp)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded scheduler evidence timestamp differs for job 12345")


def test_checkpoint_rejects_bounded_scheduler_invalid_chronology(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def reverse_chronology(evidence):
        evidence[0]["completed_utc"] = evidence[0]["submitted_utc"]

    mutate_bounded_scheduler_evidence(fixture, reverse_chronology)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded scheduler evidence 0 chronology is invalid")


def test_checkpoint_rejects_bounded_scheduler_barrier_mismatch(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    payload = json.loads(fixture["staged"].read_text())
    payload["barrier"]["recorded_segments"].pop()
    replace_bounded_artifact(fixture, payload)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded scheduler evidence job order differs from the barrier")


def test_checkpoint_rejects_bounded_scheduler_checksum_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    fixture["scheduler_files"][1].write_text("forged scheduler evidence\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded scheduler evidence 1 checksum has changed")


def test_checkpoint_rejects_bounded_source_bundle_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    fixture["source_bundle"].write_text("forged source bundle\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "source bundle checksum has changed")


def test_checkpoint_rejects_bounded_uncovered_source_revision(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    fixture["source_bundle_revisions"].append("f" * 40)
    payload = json.loads(fixture["staged"].read_text())
    payload["provenance"]["source_bundle_verified_revisions"] = (
        fixture["source_bundle_revisions"]
    )
    replace_bounded_artifact(fixture, payload)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "source bundle does not contain requested revision")


def test_checkpoint_rejects_bounded_live_generator_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    live = fixture["repository"] / "scripts/frontier/cgl_lf_stage_i_recost.py"
    live.write_text(live.read_text() + "\n# forged generator drift\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "Stage I recost generator checksum has changed")


def test_checkpoint_rejects_bounded_restart_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    restart = fixture["root"] / "runs/restart.rst"
    restart.parent.mkdir(parents=True, exist_ok=True)
    restart.write_text("authenticated restart\n")

    def add_parent(authorization):
        profile = authorization["authorized_next_profiles"][0]
        profile.update(
            {
                "segment": "s01_rankio_t0p25_t0p5",
                "parent_job_id": "34567",
                "parent_result": "clean_partial",
                "parent_segment": "s00_rankio_t0_t0p25",
                "restart_file": str(restart),
                "restart_file_sha256": sha256(restart),
                "restart_time": 0.25,
                "time_tlim_target": 0.5,
            }
        )

    mutate_bounded_authorization(fixture, add_parent)
    restart.write_text("forged restart\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded profile 0 restart checksum has changed")


def test_checkpoint_rejects_bounded_artifact_count_detail_drift(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    payload = json.loads(fixture["staged"].read_text())
    payload["reservations"]["rows"] = 1
    replace_bounded_artifact(fixture, payload)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "recost artifact reservation rows differs")


def test_checkpoint_rejects_mixed_bounded_and_sole_scheduler_bindings(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    completed = run_checkpoint(
        fixture,
        "verify-staged-recost",
        "--scheduler-relative-path",
        str(fixture["scheduler"].relative_to(fixture["root"])),
        "--expected-scheduler-sha256",
        sha256(fixture["scheduler"]),
    )
    assert_rejected(completed, "V2 mode rejects sole-profile scheduler bindings")


def test_checkpoint_rejects_bounded_noncanonical_json_pointer(
    bounded_recost_fixture,
):
    completed = run_checkpoint(
        bounded_recost_fixture,
        "verify-staged-recost",
        "--artifact-authorization-pointer",
        "/alternate_authorization",
    )
    assert_rejected(completed, "V2 mode requires canonical JSON pointer /authorization")


def test_checkpoint_v2_matches_current_recost_generator_interface():
    tree = ast.parse(RECOST.read_text())
    functions = {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert [argument.arg for argument in functions["build_payload"].args.args] == [
        "args",
        "root",
        "source_path",
        "repository",
        "generator_sha256",
    ]
    assert "run_authenticated_reconcile" in functions
    assert "require_empty_transaction_stores" in functions
    assert "require_live_storage_boundary" in functions
    assert "require_directory_measurement_boundaries" in functions
    build_result = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "BuildResult"
    )
    fields = {
        node.target.id
        for node in build_result.body
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
    }
    assert fields == {
        "payload",
        "tracker",
        "reconcile",
        "helper_path",
        "helper_sha256",
        "helper_revision",
        "matrix_path",
        "matrix_sha256",
        "matrix_revision",
        "storage_available_bytes",
        "storage_retained_stage_i_bytes",
        "storage_required_safety_bytes",
        "projected_storage_bytes",
        "storage_measurements",
    }
    checkpoint_tree = ast.parse(CHECKPOINT.read_text())
    reexecute = next(
        node
        for node in checkpoint_tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "reexecute_v2_recost"
    )
    assert any(
        isinstance(call.func, ast.Attribute)
        and call.func.attr == "require_directory_measurement_boundaries"
        and len(call.args) == 1
        and isinstance(call.args[0], ast.Attribute)
        and isinstance(call.args[0].value, ast.Name)
        and call.args[0].value.id == "build"
        and call.args[0].attr == "storage_measurements"
        for call in ast.walk(reexecute)
        if isinstance(call, ast.Call)
    )


def test_checkpoint_current_recost_generator_end_to_end(tmp_path, monkeypatch):
    recost_test = load_current_recost_test_module()
    original_write_immutable_json = recost_test.write_immutable_json

    def write_exact_f116_publisher(path, value):
        if path.name == recost_test.F116_NAME:
            value = json.loads(json.dumps(value))
            tools = value["implementation"]["committed_tools"]
            value["implementation"]["publisher"] = next(
                item
                for item in tools
                if item["path"] == "scripts/frontier/cgl_lf_stage_i_source_authority.py"
            )
        original_write_immutable_json(path, value)

    monkeypatch.setattr(
        recost_test, "write_immutable_json", write_exact_f116_publisher
    )
    generated = recost_test.recost_fixture.__wrapped__(tmp_path)
    repository = generated["repository"]
    root = generated["root"]
    accounting = generated["accounting"]
    assert isinstance(repository, Path)
    assert isinstance(root, Path)
    assert isinstance(accounting, Path)

    checkpoint = repository / "scripts/frontier/cgl_lf_stage_i_checkpoint.py"
    checkpoint.write_bytes(CHECKPOINT.read_bytes())
    checkpoint.chmod(0o755)
    force_checkpoint_renameat2_einval(checkpoint)
    retained_generator = accounting / "utilities/cgl_lf_stage_i_recost.py"
    retained_generator.parent.mkdir()
    retained_generator.write_bytes(RECOST.read_bytes())
    retained_generator.chmod(0o755)
    completed = recost_test.run_generator(generated)
    assert completed.returncode == 0, completed.stderr
    staged = generated["output"]
    request = generated["request"]
    source_bundle = generated["source_bundle"]
    helper = generated["helper"]
    revision = generated["revision"]
    queue = generated["queue"]
    assert all(
        isinstance(item, Path)
        for item in (staged, request, source_bundle, helper, queue)
    )
    assert isinstance(revision, str)
    artifact = json.loads(staged.read_text())
    assert set(artifact["provenance"]["source_authority"]) == {
        "checkpoint",
        "evidence",
        "provenance_review",
        "plasma_review",
        "publication_audit",
        "final_source_bundle",
    }
    artifact_name = staged.name.removesuffix(".staged")
    independent_review = accounting / f"{artifact_name}.independent_review.json"
    review_candidate = tmp_path / f"{artifact_name}.independent_review.candidate.json"
    write_json(
        review_candidate,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": datetime.now(timezone.utc).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": {
                "agent_id": "fixture-independent-schema2-reviewer",
                "independent_from_generator": True,
            },
            "candidate": {
                "path": str(accounting / artifact_name),
                "sha256": sha256(staged),
            },
            "scope": {"non_authorizing": True},
        },
    )
    review_candidate.chmod(0o444)
    fixture = {
        "root": root,
        "repository": repository,
        "checkpoint": checkpoint,
        "stage_i": helper,
        "accounting": accounting,
        "generator": retained_generator,
        "queue": queue,
        "counts": artifact["reconcile"]["counts"],
        "staged": staged,
        "artifact_name": artifact_name,
        "artifact_sha256": sha256(staged),
        "canonical": accounting / artifact_name,
        "audit": accounting / f"{artifact_name}.publication_audit.json",
        "stage_i_transactions": generated["transaction_store"],
        "recost_transactions": generated["recost_transaction_store"],
        "recost_forensics_root": (
            accounting / f"mks24_stage_i_{EPOCH_SLUG}_recost_forensics"
        ),
        "recost_forensics": (
            accounting / f"mks24_stage_i_{EPOCH_SLUG}_recost_forensics" / artifact_name
        ),
        "lock": root / f".mks24_stage_i_{EPOCH_SLUG}.lock",
        "recommendations": artifact["recommendations"],
        "scheduler_evidence": artifact["provenance"]["scheduler_evidence"],
        "source_bundle": source_bundle,
        "source_bundle_sha256": sha256(source_bundle),
        "stage_i_revision": artifact["provenance"]["stage_i_helper_revision"],
        "generator_revision": artifact["provenance"]["generator_revision"],
        "source_bundle_revisions": artifact["provenance"][
            "source_bundle_verified_revisions"
        ],
        "request": request,
        "request_sha256": sha256(request),
        "independent_review": independent_review,
        "independent_review_sha256": sha256(review_candidate),
    }
    candidate_symlink = tmp_path / "artifact-review-candidate-symlink.json"
    candidate_symlink.symlink_to(review_candidate)
    rejected = run_checkpoint(
        fixture,
        "install-artifact-review",
        "--artifact-review-candidate",
        str(candidate_symlink),
    )
    assert rejected.returncode == 1
    assert not independent_review.exists()

    independent_review.symlink_to(review_candidate)
    rejected = run_checkpoint(
        fixture,
        "install-artifact-review",
        "--artifact-review-candidate",
        str(review_candidate),
    )
    assert_rejected(rejected, "artifact-review target already exists or changed")
    independent_review.unlink()

    checkpoint_module = load_checkpoint_module()
    with checkpoint_module.promotion_lock(fixture):
        rejected = run_checkpoint(
            fixture,
            "install-artifact-review",
            "--artifact-review-candidate",
            str(review_candidate),
        )
    assert_rejected(rejected, "another Stage I mutation holds")
    installed = run_checkpoint(
        fixture,
        "install-artifact-review",
        "--artifact-review-candidate",
        str(review_candidate),
    )
    assert installed.returncode == 0, installed.stderr
    assert independent_review.read_bytes() == review_candidate.read_bytes()
    assert independent_review.stat().st_mode & 0o777 == 0o444
    assert independent_review.stat().st_nlink == 1
    assert review_candidate.stat().st_nlink == 1
    repeated = run_checkpoint(
        fixture,
        "install-artifact-review",
        "--artifact-review-candidate",
        str(review_candidate),
    )
    assert repeated.returncode == 0, repeated.stderr

    retained = staged.read_bytes()
    forged = json.loads(retained)
    authority = forged["provenance"]["source_authority"]
    forged["provenance"]["source_authority"] = {
        "checkpoint": "F-116",
        "evidence_sha256": authority["evidence"]["sha256"],
        "publication_audit_sha256": authority["publication_audit"]["sha256"],
        "current_source_bundle": authority["final_source_bundle"],
    }
    write_json(staged, forged)
    staged.chmod(0o444)
    fixture["artifact_sha256"] = sha256(staged)
    rejected = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(rejected, "schema-2 artifact F116 source authority binding differs")
    staged.chmod(0o644)
    staged.write_bytes(retained)
    staged.chmod(0o444)
    fixture["artifact_sha256"] = sha256(staged)

    forged = json.loads(retained)
    forged["authority"]["authorizing"] = True
    write_json(staged, forged)
    staged.chmod(0o444)
    fixture["artifact_sha256"] = sha256(staged)
    rejected = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(rejected, "improperly grants authority")
    staged.chmod(0o644)
    staged.write_bytes(retained)
    staged.chmod(0o444)
    fixture["artifact_sha256"] = sha256(staged)
    verified = run_checkpoint(fixture, "verify-staged-recost")
    assert verified.returncode == 0, verified.stderr
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    verified = run_checkpoint(fixture, "verify-promoted-recost")
    assert verified.returncode == 0, verified.stderr
    publication_audit = json.loads(fixture["audit"].read_text())
    assert publication_audit["record_type"] == (
        "stage-i-recost-recommendation-publication-audit"
    )
    assert publication_audit["authority"]["action_authority"] is False
    assert publication_audit["generalized_publication_context"][
        "controller_enforcement"
    ]["launch_authority"] is False
    assert publication_audit["generalized_publication_context"]["provenance"][
        "source_authority"
    ] == artifact["provenance"]["source_authority"]
    assert fixture["audit"].stat().st_mode & 0o777 == 0o444
    reviewed = datetime.fromisoformat(
        json.loads(independent_review.read_text())["reviewed_utc"].replace("Z", "+00:00")
    )
    publication_audit["published_utc"] = (
        reviewed - timedelta(microseconds=1)
    ).isoformat()
    write_json(fixture["audit"], publication_audit)
    fixture["audit"].chmod(0o444)
    rejected = run_checkpoint(fixture, "verify-promoted-recost")
    assert_rejected(rejected, "publication predates independent review")


def test_checkpoint_v2_preserves_sole_profile_compatibility(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    profile = fixture["authorization"]["authorized_next_profiles"][0]
    compatibility_keys = (
        "athena_walltime",
        "case_id",
        "nodes",
        "parent_job_id",
        "parent_result",
        "parent_segment",
        "restart_file",
        "restart_time",
        "segment",
        "source_bundle",
        "source_bundle_sha256",
        "time_tlim_target",
        "walltime",
    )
    authorization = {
        "mode": "sole-next-profile",
        "authorizing": True,
        "authorized_next_profiles": [profile],
        "bounded_concurrency": {
            "max_active_segments": 4,
            "max_wave_nodes": profile["nodes"],
            "r17_exclusive_and_last": True,
        },
        "controller_consumption_state": (
            "sole-profile-compatible-with-existing-checkpoint; "
            "controller consumption remains pending"
        ),
        "sole_next_segment_profile": {
            key: profile[key] for key in compatibility_keys
        },
    }
    scheduler = fixture["scheduler_evidence"][0]
    request = json.loads(fixture["request"].read_text())
    request["authorization"] = {
        "mode": "sole-next-profile",
        "max_wave_nodes": profile["nodes"],
        "profiles": [profile],
    }
    request["barrier"]["recorded_segments"] = request["barrier"]["recorded_segments"][:1]
    request["inputs"]["scheduler_evidence"] = [
        {"path": scheduler["path"], "sha256": scheduler["sha256"]}
    ]
    write_json(fixture["request"], request)
    fixture["request_sha256"] = sha256(fixture["request"])
    fixture["authorization"] = authorization
    fixture["scheduler_evidence"] = [scheduler]

    def mutate(payload):
        payload["authorization"] = authorization
        payload["barrier"]["job_ids"] = [scheduler["job_id"]]
        payload["barrier"]["recorded_segments"] = request["barrier"]["recorded_segments"]
        payload["barrier"]["scheduler_evidence"] = [scheduler]
        payload["provenance"]["request_sha256"] = fixture["request_sha256"]
        payload["provenance"]["scheduler_evidence"] = [scheduler]
        payload["provenance"]["scheduler_sha256"] = scheduler["sha256"]
        payload["storage"]["projected_authorized_wave_growth_bytes"] = (
            profile["estimated_storage_bytes"]
        )
        payload["storage"]["headroom_after_authorized_wave_and_safety_bytes"] = (
            payload["storage"]["available_bytes"]
            - payload["storage"]["required_safety_bytes"]
            - profile["estimated_storage_bytes"]
        )

    mutate_bounded_artifact(fixture, mutate)
    verified = run_checkpoint(fixture, "verify-staged-recost")
    assert verified.returncode == 0, verified.stderr
    promoted = run_checkpoint(fixture, "promote-recost")
    assert promoted.returncode == 0, promoted.stderr
    audit = json.loads(fixture["audit"].read_text())
    assert (
        audit["generalized_publication_context"]["authorization"][
            "sole_next_segment_profile"
        ]
        == authorization["sole_next_segment_profile"]
    )


def test_checkpoint_v2_rejects_request_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    request = json.loads(fixture["request"].read_text())
    request["scope"] = "forged scope"
    write_json(fixture["request"], request)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 recost request checksum has changed")


@pytest.mark.parametrize(
    "evidence",
    [
        "reconciliation",
        "ledger",
        "reservations",
        "storage",
        "f113",
        "f113_audit",
        "predecessor",
        "predecessor_audit",
        "matrix",
    ],
)
def test_checkpoint_v2_rejects_bound_evidence_drift(
    bounded_recost_fixture,
    evidence,
):
    fixture = bounded_recost_fixture
    path = fixture["evidence_files"][evidence]
    path.write_bytes(path.read_bytes() + b"\nforged drift\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "checksum has changed")


def test_checkpoint_v2_rejects_helper_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    fixture["stage_i"].write_text(fixture["stage_i"].read_text() + "\n# forged drift\n")
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "Stage I helper must be committed before recost publication")


def test_checkpoint_v2_rejects_profile_schema_extension(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def add_field(authorization):
        authorization["authorized_next_profiles"][0]["unreviewed"] = True

    mutate_bounded_authorization(fixture, add_field)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 profile 0 schema differs")


@pytest.mark.parametrize(
    ("field", "replacement", "error"),
    [
        (
            "artifact_name",
            "mks24_stage_i_E03_forcing_policy_F115_recost_evidence.json",
            "V2 artifact basename binding differs",
        ),
        ("checkpoint", "F-115", "V2 artifact checkpoint binding differs"),
        ("generated_utc", "2026-06-04T00:00:00+00:00", "generation timestamp differs"),
        ("expires_utc", "2026-06-04T00:00:01+00:00", "expiry timestamp differs"),
    ],
)
def test_checkpoint_v2_rejects_artifact_request_identity_drift(
    bounded_recost_fixture,
    field,
    replacement,
    error,
):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload[field] = replacement

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, error)


def test_checkpoint_v2_rejects_barrier_job_vector_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload["barrier"]["job_ids"].pop()

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 artifact barrier ledger-tail bindings differ")


def test_checkpoint_v2_rejects_predecessor_digest_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload["predecessor_recost"]["sha256"] = "0" * 64

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 predecessor recost binding differs")


def test_checkpoint_v2_rejects_predecessor_artifact_name_drift(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload["predecessor_recost"]["artifact_name"] = (
            "mks24_stage_i_E03_forcing_policy_F112_recost_evidence.json"
        )

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 predecessor artifact-name binding differs")


def test_checkpoint_v2_rejects_authenticated_lineage_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload["manifests"]["authenticated_lineages_sha256"] = "0" * 64

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 authenticated-lineage digest differs")


def test_checkpoint_v2_rejects_budget_projection_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload["budget"]["computed_stage_i_total_node_hours"] = "99"

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 artifact budget projection digest differs")


def test_checkpoint_v2_rejects_storage_arithmetic_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture

    def mutate(payload):
        payload["storage"]["headroom_after_authorized_wave_and_safety_bytes"] += 1

    mutate_bounded_artifact(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "V2 artifact storage arithmetic differs")


def test_checkpoint_v2_rejects_authorizing_bounded_artifact(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def mutate(authorization):
        authorization["authorizing"] = True

    mutate_bounded_authorization(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "bounded-wave V2 artifact must remain non-authorizing")


def test_checkpoint_v2_rejects_missing_controller_transition_disclosure(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture

    def mutate(authorization):
        authorization["controller_consumption_state"] = "ready"

    mutate_bounded_authorization(fixture, mutate)
    completed = run_checkpoint(fixture, "verify-staged-recost")
    assert_rejected(completed, "does not disclose pending controller enforcement")


def test_checkpoint_v2_recovery_rejects_generalized_context_tamper(
    bounded_recost_fixture,
):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    journal = next(fixture["recost_transactions"].iterdir())
    record = json.loads(journal.read_text())
    record["generalized_publication_context"]["checkpoint"] = "F-999"
    write_json(journal, record)
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(
        completed,
        "recost recovery journal generalized_publication_context binding differs",
    )


def test_checkpoint_v2_recovery_rejects_forensic_drift(bounded_recost_fixture):
    fixture = bounded_recost_fixture
    failed = run_checkpoint(fixture, "promote-recost", "--simulate-post-link-failure")
    assert_rejected(failed, "simulated post-link publication failure")
    forensic = next(fixture["recost_forensics"].iterdir())
    forensic.chmod(0o644)
    forensic.write_text("forged forensic artifact\n")
    forensic.chmod(0o444)
    completed = run_checkpoint(fixture, "finalize-linked-pair")
    assert_rejected(completed, "recost forensic copy checksum has changed")


@pytest.mark.parametrize(
    ("mode", "message"),
    [(0o770, "group-writable and untrusted"), (0o777, "world-writable and untrusted")],
)
def test_checkpoint_lock_rejects_unsafe_writable_ancestor_before_mutation(
    tmp_path, mode, message,
):
    module = load_checkpoint_module()
    unsafe = tmp_path / "unsafe"
    unsafe.mkdir()
    unsafe.chmod(mode)
    root = unsafe / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"

    with pytest.raises(ValueError, match=message):
        with module.promotion_lock({"root": root, "lock": lock}):
            raise AssertionError("unsafe ancestor must not enter the mutation boundary")

    assert not lock.exists()


def test_checkpoint_lock_binds_public_root_identity_before_public_write(tmp_path):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    detached = tmp_path / "detached-root"
    target = root / "forbidden.json"

    with pytest.raises(ValueError, match="Stage I public namespace changed"):
        with module.promotion_lock({"root": root, "lock": lock}):
            root.rename(detached)
            root.mkdir()
            module.write_json(target, {"forbidden": True})

    assert not target.exists()
    assert not (detached / target.name).exists()


def test_checkpoint_noreplace_preserves_durable_state_after_parent_authority_loss(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    live = tmp_path / "live"
    live.mkdir()
    source = live / "source"
    target = live / "target"
    source.write_bytes(b"published\n")
    detached = tmp_path / "detached"
    real_renameat2_between = module.renameat2_between
    mutations = []

    def detach_after_rename(
        source_parent, source_name, target_parent, target_name, flags, label
    ):
        mutations.append(label)
        result = real_renameat2_between(
            source_parent, source_name, target_parent, target_name, flags, label
        )
        live.rename(detached)
        live.mkdir()
        return result

    monkeypatch.setattr(module, "renameat2_between", detach_after_rename)
    with pytest.raises(ValueError, match="parent path changed during mutation"):
        with module.bound_parent_descriptor(source, "detached no-replace") as parent:
            module.rename_bound_noreplace(
                parent,
                source.name,
                target.name,
                source.stat(),
                "detached no-replace",
            )

    assert mutations == ["detached no-replace"]
    assert (detached / target.name).read_bytes() == b"published\n"
    assert list(live.iterdir()) == []


def test_checkpoint_retirement_preserves_durable_state_after_parent_authority_loss(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    live = tmp_path / "live"
    live.mkdir()
    victim = live / "victim"
    victim.write_bytes(b"retire me\n")
    detached = tmp_path / "detached"
    real_renameat2_between = module.renameat2_between
    mutations = []
    before = set(tmp_path.glob(".cgl-checkpoint-retired-*.forensic"))

    def detach_after_retirement(
        source_parent, source_name, target_parent, target_name, flags, label
    ):
        mutations.append(label)
        result = real_renameat2_between(
            source_parent, source_name, target_parent, target_name, flags, label
        )
        live.rename(detached)
        live.mkdir()
        return result

    monkeypatch.setattr(module, "renameat2_between", detach_after_retirement)
    with pytest.raises(ValueError, match="parent path changed during mutation"):
        with module.bound_parent_descriptor(victim, "detached retirement") as parent:
            module.unlink_bound_entry(
                parent, victim.name, victim.stat(), "detached retirement"
            )

    assert mutations == ["detached retirement retirement"]
    retired = set(tmp_path.glob(".cgl-checkpoint-retired-*.forensic")) - before
    assert len(retired) == 1
    assert next(iter(retired)).read_bytes() == b"retire me\n"
    assert list(live.iterdir()) == []


def test_checkpoint_ambiguous_link_is_fsynced_and_classified_before_raise(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.write_bytes(b"linked bytes\n")
    real_fsync = module.os.fsync
    fsynced = []

    class AmbiguousLinkat:
        argtypes = None
        restype = None

        def __call__(
            self, _source_parent, source_name, target_parent, target_name, _flags
        ):
            os.link(
                os.fsdecode(source_name),
                os.fsdecode(target_name),
                dst_dir_fd=target_parent,
                follow_symlinks=True,
            )
            module.ctypes.set_errno(module.errno.EIO)
            return -1

    class AmbiguousLibc:
        linkat = AmbiguousLinkat()

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced.append(module.profile_identity(profile))
        return result

    monkeypatch.setattr(module.ctypes, "CDLL", lambda *_args, **_kwargs: AmbiguousLibc())
    monkeypatch.setattr(module.os, "fsync", observe_fsync)
    with module.bound_parent_descriptor(source, "ambiguous link") as parent:
        parent_identity = module.profile_identity(os.fstat(parent))
        descriptor = os.open(source.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
        try:
            identity = module.profile_identity(os.fstat(descriptor))
            with pytest.raises(OSError):
                module.link_descriptor_noreplace(
                    descriptor, parent, target.name, identity, "ambiguous link"
                )
        finally:
            os.close(descriptor)

    assert target.read_bytes() == source.read_bytes()
    assert parent_identity in fsynced


@pytest.mark.parametrize("attack", ["detach", "insert", "replace"])
def test_checkpoint_recost_journal_scan_rejects_namespace_change(
    tmp_path, monkeypatch, attack,
):
    module = load_checkpoint_module()
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    journal = transactions / "txn.json"
    write_json(journal, {"transaction_id": "txn", "state": "preparing"})
    detached = tmp_path / "detached-transactions"
    retained_journal = tmp_path / "retained-journal.json"
    paths = {"recost_transactions": transactions}
    real_read = module.read_bound_json_object
    attacked = False

    def attack_after_read(*args, **kwargs):
        nonlocal attacked
        retained = real_read(*args, **kwargs)
        if not attacked and kwargs.get("expected_mode") == 0o644:
            attacked = True
            if attack == "detach":
                transactions.rename(detached)
                transactions.mkdir()
            elif attack == "insert":
                (transactions / "inserted").write_text("concurrent insertion\n")
            else:
                journal.rename(retained_journal)
                write_json(journal, {"transaction_id": "txn", "state": "forged"})
        return retained

    monkeypatch.setattr(module, "read_bound_json_object", attack_after_read)
    expected = "parent path changed during mutation" if attack == "detach" else "changed"
    with pytest.raises(ValueError, match=expected):
        with module.bound_directory_descriptor(
            transactions, "recost transaction directory"
        ) as descriptor:
            module.recost_journal(paths, descriptor)

    assert attacked
    if attack == "detach":
        assert (detached / journal.name).exists()
    elif attack == "replace":
        assert retained_journal.exists()
    else:
        assert journal.exists()


def test_checkpoint_json_predecessor_survives_failed_public_authentication(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "journal.json"
    target.write_text('{"transaction_id": "old"}\n')
    old = target.read_bytes()
    real_authenticate = module.require_bound_entry_security
    failed = False

    def fail_new_public_authentication(parent, name, expected, label, **kwargs):
        nonlocal failed
        result = real_authenticate(parent, name, expected, label, **kwargs)
        if (
            not failed
            and name == target.name
            and label == "JSON atomic write"
            and target.read_bytes() != old
        ):
            failed = True
            raise ValueError("simulated failed public journal authentication")
        return result

    monkeypatch.setattr(
        module, "require_bound_entry_security", fail_new_public_authentication
    )
    with pytest.raises(ValueError):
        module.write_json(target, {"transaction_id": "new"})

    assert failed
    assert json.loads(target.read_text()) == {"transaction_id": "new"}
    assert not list(tmp_path.glob(f".{target.name}.*.tmp"))
    retired = list(tmp_path.parent.glob(".cgl-checkpoint-retired-*.forensic"))
    assert any(path.read_bytes() == old for path in retired)


def test_checkpoint_recost_journal_exactly_recovers_valid_temp_after_forged_public(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    public = transactions / "txn.json"
    public.write_text("forged public bytes\n")
    temporary = transactions / f".{public.name}.123.{'0' * 32}.tmp"
    expected_record = {"transaction_id": "txn", "state": "link-pending"}
    write_json(temporary, expected_record)
    paths = {"recost_transactions": transactions}

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    with module.bound_directory_descriptor(
        transactions, "recost transaction directory"
    ) as descriptor:
        journal, record = module.recost_journal(paths, descriptor)

    assert journal == public
    assert record == expected_record
    assert json.loads(public.read_text()) == expected_record
    assert list(transactions.iterdir()) == [public]
    retired = list(tmp_path.glob(".cgl-checkpoint-retired-*.forensic"))
    assert len(retired) == 1
    assert retired[0].read_bytes() == b"forged public bytes\n"


@pytest.mark.parametrize("phase", ["initial-create", "replacement"])
def test_checkpoint_lustre_recost_journal_recovers_post_link_pre_unlink_state(
    tmp_path, monkeypatch, phase,
):
    module = load_checkpoint_module()
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    public = transactions / "txn.json"
    if phase == "replacement":
        write_json(public, {"transaction_id": "txn", "state": "old"})
    record = {"transaction_id": "txn", "state": phase}
    temporary = leave_lustre_json_post_link_state(
        module, public, record, monkeypatch
    )
    paths = {"recost_transactions": transactions}

    assert temporary.stat().st_ino == public.stat().st_ino
    assert temporary.stat().st_nlink == public.stat().st_nlink == 2
    with module.bound_directory_descriptor(
        transactions, "recost transaction directory"
    ) as descriptor:
        journal, recovered = module.recost_journal(paths, descriptor)

    assert journal == public
    assert recovered == record
    assert not temporary.exists()
    assert public.stat().st_nlink == 1
    assert list(transactions.iterdir()) == [public]


@pytest.mark.parametrize("state", ["public", "retirement-linked", "public-absent"])
def test_checkpoint_lustre_recost_journal_recovers_every_forward_state(
    tmp_path, monkeypatch, state,
):
    module = load_checkpoint_module()
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    public = transactions / "txn.json"
    old_record = {"transaction_id": "txn", "state": "old"}
    new_record = {"transaction_id": "txn", "state": "new"}
    write_json(public, old_record)
    old = public.read_bytes()
    temporary = transactions / f".{public.name}.123.{'0' * 32}.tmp"
    write_json(temporary, new_record)
    paths = {"recost_transactions": transactions}

    with module.bound_directory_descriptor(
        transactions, "fixture recost transaction directory"
    ) as descriptor:
        retired = tmp_path / module.deterministic_retirement_name(
            descriptor, public.name, public.stat()
        )
    if state in {"retirement-linked", "public-absent"}:
        os.link(public, retired)
    if state == "public-absent":
        public.unlink()

    def unsupported_renameat2(*_args, **_kwargs):
        raise OSError(errno.EINVAL, os.strerror(errno.EINVAL))

    monkeypatch.setattr(module, "renameat2", unsupported_renameat2)
    monkeypatch.setattr(module, "renameat2_between", unsupported_renameat2)
    with module.bound_directory_descriptor(
        transactions, "recost transaction directory"
    ) as descriptor:
        journal, record = module.recost_journal(paths, descriptor)

    assert journal == public
    assert record == new_record
    assert json.loads(public.read_text()) == new_record
    assert list(transactions.iterdir()) == [public]
    retained = [
        path
        for path in tmp_path.glob(".cgl-checkpoint-retired-*.forensic")
        if path.read_bytes() == old
    ]
    assert len(retained) == 1


@pytest.mark.parametrize("ambiguous", [False, True])
def test_checkpoint_artifact_review_is_idempotent_for_exact_publication(
    tmp_path, monkeypatch, ambiguous,
):
    module = load_checkpoint_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    target = accounting / "artifact.independent_review.json"
    retained = b'{"review": "exact candidate"}\n'
    expected = hashlib.sha256(retained).hexdigest()
    before_root = set(tmp_path.iterdir())
    if not ambiguous:
        target.write_bytes(retained)
        target.chmod(0o444)

        def forbid_temporary(*_args, **_kwargs):
            raise AssertionError("exact-existing review must not create a temporary")

        monkeypatch.setattr(module, "write_bound_exclusive", forbid_temporary)
    else:
        real_renameat2 = module.renameat2

        def raise_after_publication(*args, **kwargs):
            real_renameat2(*args, **kwargs)
            raise OSError("ambiguous success after review publication")

        monkeypatch.setattr(module, "renameat2", raise_after_publication)
    with module.bound_directory_descriptor(accounting, "accounting directory") as descriptor:
        identity, created = module.create_artifact_review(
            descriptor, target, retained, expected
        )

    assert target.read_bytes() == retained
    assert identity == (target.stat().st_dev, target.stat().st_ino)
    assert created is ambiguous
    assert list(accounting.iterdir()) == [target]
    assert set(tmp_path.iterdir()) == before_root


def test_checkpoint_artifact_review_idempotency_rejects_differing_target(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    target = accounting / "artifact.independent_review.json"
    target.write_bytes(b"different review\n")
    target.chmod(0o444)
    retained = b'{"review": "exact candidate"}\n'
    expected = hashlib.sha256(retained).hexdigest()
    before_root = set(tmp_path.iterdir())

    def forbid_temporary(*_args, **_kwargs):
        raise AssertionError("differing existing review must not create a temporary")

    monkeypatch.setattr(module, "write_bound_exclusive", forbid_temporary)

    with pytest.raises(ValueError, match="already exists or changed"):
        with module.bound_directory_descriptor(
            accounting, "accounting directory"
        ) as descriptor:
            module.create_artifact_review(descriptor, target, retained, expected)

    assert target.read_bytes() == b"different review\n"
    assert list(accounting.iterdir()) == [target]
    assert set(tmp_path.iterdir()) == before_root


def test_checkpoint_artifact_review_authenticates_existing_public_name_before_temporary(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    target = accounting / "artifact.independent_review.json"
    escaped = accounting / "escaped-review.json"
    retained = b'{"review": "exact candidate"}\n'
    expected = hashlib.sha256(retained).hexdigest()
    target.write_bytes(retained)
    target.chmod(0o444)
    target_identity = module.profile_identity(target.stat())
    real_read = module.read_descriptor_bytes
    replaced = False

    def replace_public_name_after_read(descriptor):
        nonlocal replaced
        value = real_read(descriptor)
        if not replaced and module.profile_identity(os.fstat(descriptor)) == target_identity:
            replaced = True
            target.rename(escaped)
            target.write_bytes(value)
            target.chmod(0o444)
        return value

    def forbid_temporary(*_args, **_kwargs):
        raise AssertionError("changed existing review must not create a temporary")

    monkeypatch.setattr(module, "read_descriptor_bytes", replace_public_name_after_read)
    monkeypatch.setattr(module, "write_bound_exclusive", forbid_temporary)
    with pytest.raises(ValueError, match="already exists or changed"):
        with module.bound_directory_descriptor(
            accounting, "accounting directory"
        ) as descriptor:
            module.create_artifact_review(descriptor, target, retained, expected)

    assert replaced
    assert target.read_bytes() == retained
    assert escaped.read_bytes() == retained
    assert sorted(entry.name for entry in accounting.iterdir()) == sorted(
        [escaped.name, target.name]
    )


def test_checkpoint_artifact_review_retry_removes_owner_only_create_remnant_without_forensic(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    target = accounting / "artifact.independent_review.json"
    retained = b'{"review": "exact candidate"}\n'
    expected = hashlib.sha256(retained).hexdigest()
    real_open = module.os.open
    interrupted = False
    before_forensics = set(tmp_path.glob(".cgl-checkpoint-retired-*.forensic"))

    def create_then_raise(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal interrupted
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if (
            not interrupted
            and os.fsdecode(path).startswith(".cgl-checkpoint-review-")
            and flags & os.O_CREAT
        ):
            interrupted = True
            os.close(descriptor)
            raise OSError("simulated exception after review temporary create")
        return descriptor

    monkeypatch.setattr(module.os, "open", create_then_raise)
    with module.bound_directory_descriptor(accounting, "accounting directory") as descriptor:
        with pytest.raises(OSError, match="exception after review temporary create"):
            module.create_artifact_review(descriptor, target, retained, expected)

        temporaries = list(accounting.glob(".cgl-checkpoint-review-*"))
        assert len(temporaries) == 1
        assert temporaries[0].stat().st_size == 0
        assert stat.S_IMODE(temporaries[0].stat().st_mode) == 0o600

        identity, created = module.create_artifact_review(
            descriptor, target, retained, expected
        )
        retry_identity, retry_created = module.create_artifact_review(
            descriptor, target, retained, expected
        )

    assert interrupted
    assert created is True
    assert retry_created is False
    assert identity == retry_identity == module.profile_identity(target.stat())
    assert target.read_bytes() == retained
    assert list(accounting.iterdir()) == [target]
    assert set(tmp_path.glob(".cgl-checkpoint-retired-*.forensic")) == before_forensics


def test_checkpoint_write_bound_exclusive_durably_retains_ambiguous_create(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "ambiguous.tmp"
    real_open = module.os.open
    real_fsync = module.os.fsync
    fsynced = []
    injected = False

    def create_then_raise(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal injected
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if not injected and path == target.name and flags & os.O_CREAT:
            injected = True
            os.close(descriptor)
            raise OSError("simulated exception after exclusive create")
        return descriptor

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced.append(module.profile_identity(profile))
        return result

    with module.bound_parent_descriptor(target, "ambiguous exclusive create") as parent:
        parent_identity = module.profile_identity(os.fstat(parent))
        monkeypatch.setattr(module.os, "open", create_then_raise)
        monkeypatch.setattr(module.os, "fsync", observe_fsync)
        with pytest.raises(OSError, match="simulated exception after exclusive create"):
            module.write_bound_exclusive(
                parent, target.name, b"must not be written\n", 0o644, "ambiguous create"
            )

    assert injected
    assert target.stat().st_size == 0
    assert stat.S_IMODE(target.stat().st_mode) == 0o600
    assert parent_identity in fsynced


def test_checkpoint_mkdir_durably_retains_ambiguous_exact_mode_create(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    target = tmp_path / "created"
    real_mkdir = module.os.mkdir
    real_fsync = module.os.fsync
    fsynced = []
    requested_modes = []
    injected = False

    def create_then_raise(path, mode=0o777, *, dir_fd=None):
        nonlocal injected
        requested_modes.append(mode)
        real_mkdir(path, mode, dir_fd=dir_fd)
        if not injected and path == target.name:
            injected = True
            raise OSError("simulated exception after mkdir")

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced.append(module.profile_identity(profile))
        return result

    ambient = os.umask(0o077)
    try:
        monkeypatch.setattr(module.os, "mkdir", create_then_raise)
        monkeypatch.setattr(module.os, "fsync", observe_fsync)
        with pytest.raises(OSError, match="simulated exception after mkdir"):
            module.mkdir_durable(target)
    finally:
        os.umask(ambient)

    assert injected
    assert requested_modes == [0o755]
    assert stat.S_IMODE(target.stat().st_mode) == 0o755
    assert module.profile_identity(tmp_path.stat()) in fsynced
    assert module.profile_identity(target.stat()) in fsynced
    module.mkdir_durable(target)


def test_checkpoint_promotion_lock_durably_retains_ambiguous_create(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    real_open = module.os.open
    real_fsync = module.os.fsync
    fsynced = []
    injected = False

    def create_then_raise(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal injected
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if not injected and path == lock.name and flags & os.O_CREAT:
            injected = True
            os.close(descriptor)
            raise OSError("simulated exception after lock create")
        return descriptor

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced.append(module.profile_identity(profile))
        return result

    monkeypatch.setattr(module.os, "open", create_then_raise)
    monkeypatch.setattr(module.os, "fsync", observe_fsync)
    with pytest.raises(OSError, match="simulated exception after lock create"):
        with module.promotion_lock({"root": root, "lock": lock}):
            raise AssertionError("ambiguous lock create must fail closed")

    assert injected
    assert lock.read_bytes() == b""
    assert stat.S_IMODE(lock.stat().st_mode) == 0o644
    assert module.profile_identity(root.stat()) in fsynced
    with module.promotion_lock({"root": root, "lock": lock}):
        pass


def test_checkpoint_forensic_copy_durably_retains_ambiguous_create(
    tmp_path, monkeypatch,
):
    module = load_checkpoint_module()
    source = tmp_path / "source"
    target = tmp_path / "forensic"
    source.write_bytes(b"authenticated forensic bytes\n")
    source.chmod(0o444)
    expected = sha256(source)
    temporary = tmp_path / module.forensic_temporary_name(target)
    real_open = module.os.open
    real_fsync = module.os.fsync
    fsynced = []
    injected = False

    def create_then_raise(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal injected
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if not injected and path == temporary.name and flags & os.O_CREAT:
            injected = True
            os.close(descriptor)
            raise OSError("simulated exception after forensic create")
        return descriptor

    def observe_fsync(descriptor):
        profile = os.fstat(descriptor)
        result = real_fsync(descriptor)
        if stat.S_ISDIR(profile.st_mode):
            fsynced.append(module.profile_identity(profile))
        return result

    monkeypatch.setattr(module.os, "open", create_then_raise)
    monkeypatch.setattr(module.os, "fsync", observe_fsync)
    with pytest.raises(OSError, match="simulated exception after forensic create"):
        module.copy_forensic(source, target, expected, expected_mode=0o444)

    assert injected
    assert temporary.stat().st_size == 0
    assert stat.S_IMODE(temporary.stat().st_mode) == 0o600
    assert module.profile_identity(tmp_path.stat()) in fsynced
    module.copy_forensic(source, target, expected, expected_mode=0o444)
    assert target.read_bytes() == source.read_bytes()
    assert not temporary.exists()


@pytest.mark.parametrize(
    ("action_name", "first_preflight"),
    [
        ("promote", "verify_staged_state"),
        ("retire_preparing", "require_empty_queue"),
        ("adopt_legacy_canonical", "require_empty_queue"),
    ],
)
def test_checkpoint_actions_bind_transaction_store_before_first_preflight(
    tmp_path, monkeypatch, action_name, first_preflight,
):
    module = load_checkpoint_module()
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    journal = transactions / "pending.json"
    journal.write_text('{"state": "pending"}\n')
    detached = tmp_path / "detached-transactions"
    paths = {"recost_transactions": transactions}

    @contextmanager
    def no_lock(_paths):
        yield

    def detach_during_first_preflight(*_args, **_kwargs):
        transactions.rename(detached)
        transactions.mkdir()
        raise RuntimeError("stop after transaction-store detachment")

    monkeypatch.setattr(module, "promotion_lock", no_lock)
    monkeypatch.setattr(module, "validate_fixture_options", lambda *_args: None)
    monkeypatch.setattr(module, first_preflight, detach_during_first_preflight)
    action = getattr(module, action_name)
    args = SimpleNamespace(squeue_file=None)

    with pytest.raises(ValueError, match="parent path changed during mutation"):
        action(paths, tmp_path, tmp_path, False, args)

    assert (detached / journal.name).read_bytes() == b'{"state": "pending"}\n'
    assert list(transactions.iterdir()) == []
