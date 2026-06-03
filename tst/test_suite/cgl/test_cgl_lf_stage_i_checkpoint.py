"""Local-fixture regressions for retained Stage I recost publication."""

from __future__ import annotations

import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
STAGE_I = REPOSITORY / "scripts/frontier/cgl_lf_stage_i.py"
CHECKPOINT = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_checkpoint.py"
PRODUCTION_ROOT = "/lustre/orion/ast207/proj-shared/dfielding/CGL"
ARTIFACT = "mks24_stage_i_E03_forcing_policy_R02_t7p25_recost_evidence.json"
EPOCH = "E03-forcing-policy"
EPOCH_SLUG = "E03_forcing_policy"


def sha256(path: Path) -> str:
    """Return one retained fixture digest."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    """Write stable fixture JSON."""

    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def load_checkpoint_module():
    """Load the retained companion without executing its CLI."""

    spec = importlib.util.spec_from_file_location("cgl_lf_stage_i_checkpoint", CHECKPOINT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
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
            sys.executable,
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
    profile = {
        "athena_timeout": "00:55:00",
        "nodes": 1,
        "segment": "R02/s19_rankio_t7p25_t7p5",
        "slurm_walltime": "01:05:00",
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


def checkpoint_command(fixture, action: str, *extra: str,
                       utility_sha256: str | None = None,
                       artifact_sha256: str | None = None,
                       root: Path | None = None,
                       queue_file: Path | None = None,
                       include_queue_fixture: bool = True,
                       generator_relative_path: str | None = None) -> list[str]:
    """Build one fully bound local companion invocation."""

    counts = fixture["counts"]
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
        ARTIFACT,
        "--expected-artifact-sha256",
        artifact_sha256 or fixture["artifact_sha256"],
        "--generator-relative-path",
        generator_relative_path or str(fixture["generator"].relative_to(fixture["root"])),
        "--expected-generator-sha256",
        sha256(fixture["generator"]),
        "--scheduler-relative-path",
        str(fixture["scheduler"].relative_to(fixture["root"])),
        "--expected-scheduler-sha256",
        sha256(fixture["scheduler"]),
        "--authorized-next-segment-profile-json",
        json.dumps(fixture["profile"], sort_keys=True),
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


def expose_legacy_canonical(fixture) -> None:
    """Replace one staged fixture with an unattested canonical-only artifact."""

    fixture["staged"].rename(fixture["canonical"])


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
    assert not fixture["recost_transactions"].exists()


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
    environment = load_checkpoint_module().scheduler_environment()
    assert not any(key.startswith("SLURM_") for key in environment)
    assert environment["CGL_CHECKPOINT_SENTINEL"] == "retained"


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


def test_checkpoint_retires_orphan_journal_temporary(recost_fixture):
    fixture = recost_fixture
    fixture["recost_transactions"].mkdir()
    fixture["recost_forensics"].mkdir(parents=True)
    temporary = (
        fixture["recost_transactions"]
        / f".orphan.json.123.{'0' * 32}.tmp"
    )
    temporary.write_text("partial")
    temporary.chmod(0o600)
    retired = run_checkpoint(fixture, "retire-preparing")
    assert retired.returncode == 0, retired.stderr
    assert fixture["staged"].is_file()
    assert list(fixture["recost_transactions"].iterdir()) == []


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

    def replace_staged_name():
        retained = fixture["staged"].read_bytes()
        fixture["staged"].unlink()
        fixture["staged"].write_bytes(retained)
        fixture["staged"].chmod(0o644)

    with pytest.raises(ValueError, match="links, expected 2"):
        module.unlink_linked_pair_after_empty_queue(
            fixture["root"],
            True,
            str(fixture["queue"]),
            paths,
            args,
            before_unlink=replace_staged_name,
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

    def weaken_accounting_directory():
        fixture["accounting"].chmod(0o777)

    with pytest.raises(ValueError, match="exceeds trusted profile 0755"):
        module.unlink_linked_pair_after_empty_queue(
            fixture["root"],
            True,
            str(fixture["queue"]),
            paths,
            args,
            before_unlink=weaken_accounting_directory,
        )
    assert fixture["staged"].is_file()
    assert fixture["canonical"].is_file()


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


def test_checkpoint_rejects_forged_descriptor_source_profile(recost_fixture, tmp_path):
    fixture = recost_fixture
    root = tmp_path / "forged-source"
    source = root / "scripts/frontier/cgl_lf_stage_i_checkpoint.py"
    source.parent.mkdir(parents=True)
    source.write_bytes(CHECKPOINT.read_bytes())
    source.chmod(0o644)
    descriptor = os.memfd_create("forged-cgl-checkpoint")
    try:
        os.write(descriptor, CHECKPOINT.read_bytes())
        os.lseek(descriptor, 0, os.SEEK_SET)
        environment = dict(os.environ)
        environment["_CGL_LF_RECOST_UTILITY_DESCRIPTOR"] = str(descriptor)
        environment["_CGL_LF_RECOST_UTILITY_SOURCE"] = str(source)
        environment["_CGL_LF_RECOST_REPOSITORY_ROOT"] = str(root)
        completed = subprocess.run(
            [
                sys.executable,
                f"/proc/self/fd/{descriptor}",
                *checkpoint_command(fixture, "verify-staged-recost")[2:],
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
            pass_fds=(descriptor,),
        )
    finally:
        os.close(descriptor)
    assert_rejected(completed, "retained utility mode is 0644, expected 0755")
