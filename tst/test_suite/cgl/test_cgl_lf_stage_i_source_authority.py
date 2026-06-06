"""Adversarial tests for the F-118 current-source-authority publisher."""

from __future__ import annotations

import ast
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
from types import SimpleNamespace

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
PUBLISHER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_source_authority.py"
SPEC = importlib.util.spec_from_file_location("cgl_lf_stage_i_source_authority", PUBLISHER)
assert SPEC is not None and SPEC.loader is not None
authority = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(authority)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retained_tree_snapshot(root: Path) -> dict[str, tuple[object, ...]]:
    """Capture retained bytes and all mutation-relevant inode profile fields."""

    retained = {}
    pending = [root]
    while pending:
        path = pending.pop()
        profile = path.lstat()
        relative = "." if path == root else path.relative_to(root).as_posix()
        payload: bytes | str | None = None
        if stat.S_ISREG(profile.st_mode):
            payload = path.read_bytes()
        elif stat.S_ISLNK(profile.st_mode):
            payload = os.readlink(path)
        elif stat.S_ISDIR(profile.st_mode):
            pending.extend(
                reversed(sorted(path.iterdir(), key=lambda child: child.name))
            )
        retained[relative] = (
            profile.st_mode,
            profile.st_dev,
            profile.st_ino,
            profile.st_uid,
            profile.st_gid,
            profile.st_nlink,
            profile.st_size,
            profile.st_mtime_ns,
            profile.st_ctime_ns,
            payload,
        )
    return retained


def test_cli_reachable_graph_excludes_legacy_namespace_replacement() -> None:
    tree = ast.parse(PUBLISHER.read_text())
    functions = {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    calls = {name: set() for name in functions}
    for name, node in functions.items():
        for child in ast.walk(node):
            if (
                isinstance(child, ast.Call)
                and isinstance(child.func, ast.Name)
                and child.func.id in functions
            ):
                calls[name].add(child.func.id)

    reachable = set()
    pending = ["main"]
    while pending:
        name = pending.pop()
        if name in reachable:
            continue
        reachable.add(name)
        pending.extend(calls.get(name, ()))

    assert reachable.isdisjoint({
        "atomic_write",
        "cleanup_transaction",
        "ensure_catalog",
        "ensure_installed",
        "exchange_bound_entries",
        "publish_bound_file_noreplace",
        "purge_abandoned_staging_transaction",
        "purge_forensic_only_prepublication_transaction_root",
        "purge_retired_transaction",
        "recover_atomic_journal",
        "rename_bound_entry",
        "rename_bound_noreplace",
        "renameat2",
        "renameat2_between",
        "retire_bound_recovery_file",
        "rmdir_bound_entry",
        "rmdir_bound_path",
        "unlink_bound_entry",
        "update_journal",
    })


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def write_bytes(path: Path, payload: bytes, mode: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    path.chmod(mode)


def write_json(path: Path, value: object, mode: int = 0o600) -> str:
    write_bytes(path, canonical_json(value), mode)
    return sha256(path)


def drift_bound_file(parent: int, name: str, drift: str, replacement: bytes,
                     mode: int) -> None:
    if drift == "mode":
        os.chmod(name, 0o666, dir_fd=parent, follow_symlinks=False)
        return
    if drift == "link":
        os.link(
            name,
            f"{name}.extra-link",
            src_dir_fd=parent,
            dst_dir_fd=parent,
            follow_symlinks=False,
        )
        return
    if drift != "content":
        raise AssertionError(f"unknown drift: {drift}")
    os.chmod(name, mode | stat.S_IWUSR, dir_fd=parent, follow_symlinks=False)
    descriptor = os.open(name, os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW, dir_fd=parent)
    try:
        assert os.write(descriptor, replacement) == len(replacement)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.chmod(name, mode, dir_fd=parent, follow_symlinks=False)


def assert_retained_journal_recovery(transaction: Path, payload: bytes) -> Path:
    retained = [
        transaction / name
        for name in authority.JOURNAL_RECOVERY_NAMES
        if (transaction / name).exists()
    ]
    assert len(retained) == 1
    assert retained[0].read_bytes() == payload
    assert stat.S_IMODE(retained[0].stat().st_mode) == 0o600
    return retained[0]


def configure_fake_canonical_namespace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path]:
    boundary = tmp_path / "proj-shared"
    owner = boundary / "dfielding"
    root = owner / "CGL"
    root.mkdir(parents=True)
    boundary.chmod(0o2770)
    owner.chmod(0o2755)
    root.chmod(0o2755)
    profiles = tuple(
        (
            path,
            stat.S_IMODE(path.stat().st_mode),
            path.stat().st_uid,
            path.stat().st_gid,
        )
        for path in (boundary, owner, root)
    )
    monkeypatch.setattr(authority, "CANONICAL_TRUSTED_PROJECT_BOUNDARY", boundary)
    monkeypatch.setattr(authority, "CANONICAL_TRUSTED_PROJECT_GID", boundary.stat().st_gid)
    monkeypatch.setattr(authority, "CANONICAL_OWNER_UID", root.stat().st_uid)
    monkeypatch.setattr(authority, "CANONICAL_PUBLIC_NAMESPACE_PROFILES", profiles)
    monkeypatch.setattr(authority, "DEFAULT_ROOT", root)
    return boundary, root


def git(repository: Path, *arguments: str, check: bool = True) -> subprocess.CompletedProcess:
    completed = subprocess.run(
        ["/usr/bin/git", "-C", str(repository), *arguments],
        check=False,
        capture_output=True,
        text=True,
        env={
            **os.environ,
            "GIT_CONFIG_GLOBAL": "/dev/null",
            "GIT_CONFIG_NOSYSTEM": "1",
            "LC_ALL": "C",
        },
    )
    if check and completed.returncode:
        raise AssertionError(
            f"git {' '.join(arguments)} failed:\n{completed.stdout}\n{completed.stderr}"
        )
    return completed


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def published_binding(path: Path, digest: str, mode: str) -> dict[str, object]:
    return {"path": str(path), "sha256": digest, "mode": mode, "links": 1}


@dataclass
class Campaign:
    repository: Path
    root: Path
    candidates: Path
    publisher: Path
    f115_head: str
    bridge_head: str
    f116_head: str
    final_head: str
    f115_bundle: Path
    bridge_bundle: Path
    f116_bundle: Path
    final_bundle: Path
    evidence_candidate: Path
    provenance_candidate: Path
    plasma_candidate: Path
    audit_candidate: Path
    f115_bindings: dict[str, object]
    f115_digests: dict[str, str]
    f116_bindings: dict[str, object]
    f116_digests: dict[str, str]
    final_verified_revisions: list[str]
    reviewer_ids: dict[str, str]
    authorization: dict[str, object]
    subject_override: str | None
    expected: dict[str, str]

    @classmethod
    def create(cls, temporary: Path) -> "Campaign":
        repository = temporary / "repository"
        root = temporary / "campaign"
        candidates = temporary / "candidates"
        repository.mkdir()
        candidates.mkdir()
        git(repository, "init", "-b", "feature/cgl-landau-fluid")
        git(repository, "config", "user.email", "fixture@example.invalid")
        git(repository, "config", "user.name", "source-authority fixture")

        publisher = repository / authority.PUBLISHER_RELATIVE
        write_bytes(publisher, PUBLISHER.read_bytes(), 0o755)
        for relative, mode_text in authority.REQUIRED_TOOLS.items():
            path = repository / relative
            if path == publisher:
                continue
            write_bytes(path, f"fixture {relative}\n".encode(), int(mode_text, 8))
        write_bytes(repository / "history.txt", b"historical F115 source\n", 0o644)
        git(repository, "add", ".")
        git(repository, "commit", "-m", "fixture F115 source")
        f115_head = git(repository, "rev-parse", "HEAD").stdout.strip()

        accounting = root / "accounting"
        archives = root / "source-archives"
        accounting.mkdir(parents=True)
        archives.mkdir()
        for directory in (root, accounting, archives):
            directory.chmod(0o2755)
        write_bytes(
            root / f".mks24_stage_i_{authority.EXECUTION_EPOCH_SLUG}.lock",
            b"",
            0o644,
        )
        f115_bundle = archives / f"athenak-feature-cgl-through-{f115_head[:9]}.bundle"
        git(repository, "bundle", "create", str(f115_bundle), "HEAD")
        write_bytes(repository / "bridge.txt", b"bridge source\n", 0o644)
        git(repository, "add", "bridge.txt")
        git(repository, "commit", "-m", "fixture bridge source")
        bridge_head = git(repository, "rev-parse", "HEAD").stdout.strip()
        bridge_bundle = archives / f"athenak-feature-cgl-through-{bridge_head[:9]}.bundle"
        git(
            repository,
            "bundle",
            "create",
            str(bridge_bundle),
            "refs/heads/feature/cgl-landau-fluid",
        )
        f115_bundle.chmod(0o644)
        bridge_bundle.chmod(0o644)

        write_bytes(repository / "f116.txt", b"F116 selected source\n", 0o644)
        git(repository, "add", "f116.txt")
        git(repository, "commit", "-m", "fixture F116 selected source")
        f116_head = git(repository, "rev-parse", "HEAD").stdout.strip()
        f116_bundle = archives / f"athenak-feature-cgl-through-{f116_head[:9]}.bundle"
        git(repository, "bundle", "create", str(f116_bundle), "HEAD")
        f116_bundle.chmod(0o644)

        f115_bindings, f115_digests = cls._publish_f115(
            root, repository, f115_bundle, f115_head
        )
        old_readme = b"# Source archives\n\n## AthenaK\n\nHistorical bundles.\n"
        old_sums = f"{sha256(f115_bundle)}  {f115_bundle.name}\n".encode()
        write_bytes(archives / "README.md", old_readme, 0o644)
        write_bytes(archives / "SHA256SUMS", old_sums, 0o644)
        write_bytes(archives / authority.CORRUPT_C7_NAME, b"retained incident evidence\n", 0o644)
        f116_bindings, f116_digests = cls._publish_f116(
            root,
            repository,
            f115_bindings,
            f115_digests,
            f115_head,
            bridge_bundle,
            bridge_head,
            f116_bundle,
            f116_head,
        )

        write_bytes(repository / "final.txt", b"final F118 source\n", 0o644)
        git(repository, "add", "final.txt")
        git(repository, "commit", "-m", "fixture final F118 source")
        final_head = git(repository, "rev-parse", "HEAD").stdout.strip()
        final_bundle = candidates / f"athenak-feature-cgl-through-{final_head[:9]}.bundle"
        git(repository, "bundle", "create", str(final_bundle), "HEAD")
        final_bundle.chmod(0o644)

        campaign = cls(
            repository=repository,
            root=root,
            candidates=candidates,
            publisher=publisher,
            f115_head=f115_head,
            bridge_head=bridge_head,
            f116_head=f116_head,
            final_head=final_head,
            f115_bundle=f115_bundle,
            bridge_bundle=bridge_bundle,
            f116_bundle=f116_bundle,
            final_bundle=final_bundle,
            evidence_candidate=candidates / f"{authority.F118_NAME}.candidate",
            provenance_candidate=(
                candidates / f"{authority.F118_NAME}.provenance_security_review.json.candidate"
            ),
            plasma_candidate=(
                candidates / f"{authority.F118_NAME}.plasma_scientific_review.json.candidate"
            ),
            audit_candidate=(
                candidates / f"{authority.F118_NAME}.publication_audit.json.candidate"
            ),
            f115_bindings=f115_bindings,
            f115_digests=f115_digests,
            f116_bindings=f116_bindings,
            f116_digests=f116_digests,
            final_verified_revisions=[f115_head, bridge_head, f116_head, final_head],
            reviewer_ids={"provenance": "fixture-provenance", "plasma": "fixture-plasma"},
            authorization=dict(authority.AUTHORIZATION),
            subject_override=None,
            expected={},
        )
        campaign.refresh_candidates()
        return campaign

    @staticmethod
    def _publish_f115(
        root: Path, repository: Path, bundle: Path, head: str
    ) -> tuple[dict[str, object], dict[str, str]]:
        del repository
        paths = {key: root / relative for key, relative in authority.F115_PATHS.items()}
        evidence = {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-recovery-supersession-evidence",
            "checkpoint": "F-115",
            "execution_epoch": authority.EXECUTION_EPOCH,
            "implementation": {
                "source_bundle": {
                    "path": bundle.relative_to(root).as_posix(),
                    "sha256": sha256(bundle),
                    "complete_history": True,
                    "head": head,
                    "verified_revisions": [head],
                }
            },
        }
        evidence_sha = write_json(paths["evidence"], evidence, 0o444)
        reviews: dict[str, tuple[dict[str, object], str]] = {}
        for key, kind, decision, agent in (
            ("provenance_review", "provenance-security", "approved-for-publication", "f115-p"),
            ("plasma_review", "plasma-scientific-continuation", "approved", "f115-s"),
        ):
            review = {
                "schema_version": 1,
                "record_type": "stage-i-source-bundle-recovery-supersession-independent-review",
                "checkpoint": "F-115",
                "execution_epoch": authority.EXECUTION_EPOCH,
                "review_kind": kind,
                "decision": decision,
                "published_f115": {"path": str(paths["evidence"]), "sha256": evidence_sha},
                "reviewer": {"agent_id": agent, "identity": f"fixture {agent}"},
            }
            reviews[key] = (review, write_json(paths[key], review, 0o444))
        audit = {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-recovery-supersession-publication-audit",
            "checkpoint": "F-115",
            "execution_epoch": authority.EXECUTION_EPOCH,
            "artifact": published_binding(paths["evidence"], evidence_sha, "0444"),
            "independent_reviews": {
                "reviews_bind_exact_published_f115_sha256": evidence_sha,
                "provenance_security": published_binding(
                    paths["provenance_review"], reviews["provenance_review"][1], "0444"
                ),
                "plasma_scientific_continuation": published_binding(
                    paths["plasma_review"], reviews["plasma_review"][1], "0444"
                ),
            },
            "authority_and_enforcement": {"direct_sbatch_authorized": False},
        }
        audit_sha = write_json(paths["publication_audit"], audit, 0o444)
        digests = {
            "evidence_sha256": evidence_sha,
            "publication_audit_sha256": audit_sha,
            "provenance_review_sha256": reviews["provenance_review"][1],
            "plasma_review_sha256": reviews["plasma_review"][1],
        }
        bindings = {
            key: {
                "path": authority.F115_PATHS[key].as_posix(),
                "sha256": (
                    evidence_sha
                    if key == "evidence"
                    else audit_sha
                    if key == "publication_audit"
                    else reviews[key][1]
                ),
            }
            for key in authority.F115_PATHS
        }
        return bindings, digests

    @staticmethod
    def _publish_f116(
        root: Path,
        repository: Path,
        f115_bindings: dict[str, object],
        f115_digests: dict[str, str],
        f115_head: str,
        bridge_bundle: Path,
        bridge_head: str,
        bundle: Path,
        head: str,
    ) -> tuple[dict[str, object], dict[str, str]]:
        paths = {key: root / relative for key, relative in authority.F116_PATHS.items()}
        archives = root / "source-archives"
        old_readme = (archives / "README.md").read_bytes()
        old_sums = (archives / "SHA256SUMS").read_bytes()
        subject = git(repository, "show", "-s", "--format=%s", head).stdout.strip()
        block = (
            f"`{bundle.name}` records complete history through commit `{head}` "
            f"(`{subject}`). It is the sole current Stage I source-selection bundle "
            "after independently reviewed F-116 publication. It preserves F-115 as "
            "immutable historical R03 s02 authority and does not itself authorize prepare "
            f"or submission. Its SHA-256 is `{sha256(bundle)}`.\n\n"
            f"`{bridge_bundle.name}` records complete history through commit `{bridge_head}`. "
            "It is retained and cataloged as a non-current branch-ref bridge between F-115 "
            "and the final F-116 tooling revision; it is never selected as current source "
            f"authority. Its SHA-256 is `{sha256(bridge_bundle)}`.\n\n"
        ).encode()
        new_readme = old_readme.replace(b"## AthenaK\n\n", b"## AthenaK\n\n" + block, 1)
        new_sums = old_sums + (
            f"{sha256(bridge_bundle)}  {bridge_bundle.name}\n"
            f"{sha256(bundle)}  {bundle.name}\n"
        ).encode()
        before = {
            "readme_sha256": hashlib.sha256(old_readme).hexdigest(),
            "sha256sums_sha256": hashlib.sha256(old_sums).hexdigest(),
            "bridge_listed": False,
            "final_bundle_listed": False,
            "corrupt_c7_listed": False,
        }
        after = {
            "readme_sha256": hashlib.sha256(new_readme).hexdigest(),
            "sha256sums_sha256": hashlib.sha256(new_sums).hexdigest(),
            "bridge_listed_exactly_once": True,
            "final_bundle_listed_exactly_once": True,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
            "sole_current_source_bundle": f"source-archives/{bundle.name}",
        }
        bridge = {
            "path": f"source-archives/{bridge_bundle.name}",
            "sha256": sha256(bridge_bundle),
            "complete_history": True,
            "head": bridge_head,
            "advertised_tip": {
                "revision": bridge_head,
                "name": "refs/heads/feature/cgl-landau-fluid",
            },
            "verified_revisions": sorted([f115_head, bridge_head]),
            "selected_as_current": False,
            "role": "retained-non-current-bridge",
        }
        current = {
            "candidate_path": str(root.parent / "f116-bundle-candidate"),
            "path": f"source-archives/{bundle.name}",
            "sha256": sha256(bundle),
            "complete_history": True,
            "head": head,
            "advertised_tip": {"revision": head, "name": "HEAD"},
            "verified_revisions": sorted([f115_head, bridge_head, head]),
            "selected_as_current": True,
            "subject": subject,
        }
        tools = [
            {
                "path": relative,
                "revision": head,
                "sha256": sha256(repository / relative),
                "mode": mode,
            }
            for relative, mode in sorted(authority.REQUIRED_TOOLS.items())
        ]
        publisher = next(
            item for item in tools if item["path"] == authority.PUBLISHER_RELATIVE.as_posix()
        )
        generated = now()
        evidence = {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-evidence",
            "checkpoint": "F-116",
            "execution_epoch": authority.EXECUTION_EPOCH,
            "generated_utc": generated,
            "scope": {
                "relationship": "current-source-selection-only-supersession",
                "summary": "Select the exact committed final tooling source without execution authority.",
                "preserves": authority.F116_PRESERVES,
                "does_not_authorize": authority.DOES_NOT_AUTHORIZE,
            },
            "predecessor_authorities": {"historical_f115": f115_bindings},
            "implementation": {
                "publisher": publisher,
                "committed_tools": tools,
                "intermediate_36140_bundle": bridge,
                "current_source_bundle": current,
            },
            "source_archive_catalog": {"before": before, "after": after},
            "authorization": authority.AUTHORIZATION,
            "validation": authority.F116_VALIDATION_CLAIMS,
            "publication_requirements": authority.PUBLICATION_REQUIREMENTS,
        }
        evidence_sha = write_json(paths["evidence"], evidence, 0o444)
        verified = {
            "authorization_broadening": False,
            "bridge_selected_as_current": False,
            "corrupt_c7_excluded": True,
            "current_source_selection_only": True,
            "final_bundle_sha256": sha256(bundle),
            "final_head": head,
            "historical_f115_preserved": True,
        }
        reviews: dict[str, tuple[dict[str, object], str]] = {}
        for key, kind, decision, agent in (
            ("provenance_review", "provenance-security", "approved-for-publication", "f116-p"),
            ("plasma_review", "plasma-scientific-continuation", "approved", "f116-s"),
        ):
            review = {
                "schema_version": 1,
                "record_type": (
                    "stage-i-current-source-authority-supersession-independent-review"
                ),
                "checkpoint": "F-116",
                "execution_epoch": authority.EXECUTION_EPOCH,
                "review_kind": kind,
                "decision": decision,
                "reviewed_candidate": {
                    "path": str(root.parent / f"{key}.f116-candidate"),
                    "sha256": evidence_sha,
                },
                "published_f116": {"path": str(paths["evidence"]), "sha256": evidence_sha},
                "reviewer": {"agent_id": agent, "identity": f"fixture {agent}"},
                "reviewed_utc": now(),
                "findings": ["Exact F116 predecessor authority verified."],
                "limitations": [
                    "No execution authority.",
                    authority.INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION,
                ],
                "verified": verified,
            }
            reviews[key] = (review, write_json(paths[key], review, 0o444))
        audit = {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-publication-audit",
            "checkpoint": "F-116",
            "execution_epoch": authority.EXECUTION_EPOCH,
            "published_utc": now(),
            "artifact": published_binding(paths["evidence"], evidence_sha, "0444"),
            "independent_reviews": {
                "reviews_bind_exact_published_f116_sha256": evidence_sha,
                "provenance_security": published_binding(
                    paths["provenance_review"], reviews["provenance_review"][1], "0444"
                ),
                "plasma_scientific_continuation": published_binding(
                    paths["plasma_review"], reviews["plasma_review"][1], "0444"
                ),
            },
            "historical_f115_authority": f115_digests,
            "source_archive_catalog": {
                "readme": published_binding(
                    archives / "README.md", after["readme_sha256"], "0644"
                ),
                "sha256sums": published_binding(
                    archives / "SHA256SUMS", after["sha256sums_sha256"], "0644"
                ),
                "bridge_bundle": {
                    "path": str(bridge_bundle),
                    "sha256": sha256(bridge_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": bridge_head,
                    "role": "retained-non-current-bridge",
                    "selected_as_current": False,
                },
                "current_source_bundle": {
                    "path": str(bundle),
                    "sha256": sha256(bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": head,
                    "selected_as_current": True,
                },
                "corrupt_c7_absent_from_active_checksum_ledger": True,
                "sole_current_source_bundle": str(bundle),
            },
            "authority_and_enforcement": authority.AUTHORIZATION,
            "publication": (
                "recoverable-forward-transaction-with-publication-audit-commit-marker-"
                "under-stage-i-lock"
            ),
        }
        audit_sha = write_json(paths["publication_audit"], audit, 0o444)
        write_bytes(archives / "README.md", new_readme, 0o644)
        write_bytes(archives / "SHA256SUMS", new_sums, 0o644)
        digests = {
            "evidence_sha256": evidence_sha,
            "publication_audit_sha256": audit_sha,
            "provenance_review_sha256": reviews["provenance_review"][1],
            "plasma_review_sha256": reviews["plasma_review"][1],
        }
        bindings = {
            key: {
                "path": authority.F116_PATHS[key].as_posix(),
                "sha256": (
                    evidence_sha
                    if key == "evidence"
                    else audit_sha
                    if key == "publication_audit"
                    else reviews[key][1]
                ),
            }
            for key in authority.F116_PATHS
        }
        return bindings, digests

    @property
    def final_target(self) -> Path:
        return self.root / "source-archives" / self.final_bundle.name

    def committed_tools(self) -> list[dict[str, object]]:
        return [
            {
                "path": relative,
                "revision": self.final_head,
                "sha256": sha256(self.repository / relative),
                "mode": mode,
            }
            for relative, mode in sorted(authority.REQUIRED_TOOLS.items())
        ]

    def refresh_candidates(self) -> None:
        old_readme = (self.root / "source-archives/README.md").read_bytes()
        old_sums = (self.root / "source-archives/SHA256SUMS").read_bytes()
        final_sha = sha256(self.final_bundle)
        bridge_sha = sha256(self.bridge_bundle)
        subject = git(
            self.repository, "show", "-s", "--format=%s", self.final_head
        ).stdout.strip()
        declared_subject = self.subject_override or subject
        new_readme, new_sums = authority.catalog_payloads(
            old_readme,
            old_sums,
            bridge_name=self.bridge_bundle.name,
            predecessor_name=self.f116_bundle.name,
            final_name=self.final_bundle.name,
            final_sha256=final_sha,
            final_revision=self.final_head,
            final_subject=declared_subject,
        )
        before = {
            "readme_sha256": hashlib.sha256(old_readme).hexdigest(),
            "sha256sums_sha256": hashlib.sha256(old_sums).hexdigest(),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed": False,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
        }
        after = {
            "readme_sha256": hashlib.sha256(new_readme).hexdigest(),
            "sha256sums_sha256": hashlib.sha256(new_sums).hexdigest(),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed_exactly_once": True,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
            "historical_f116_preserved": True,
            "all_prior_checksum_entries_preserved": True,
            "sole_current_source_bundle": f"source-archives/{self.final_bundle.name}",
        }
        final = {
            "candidate_path": str(self.final_bundle),
            "path": f"source-archives/{self.final_bundle.name}",
            "sha256": final_sha,
            "complete_history": True,
            "head": self.final_head,
            "advertised_tip": {"revision": self.final_head, "name": "HEAD"},
            "verified_revisions": list(self.final_verified_revisions),
            "selected_as_current": True,
            "subject": declared_subject,
        }
        bridge = {
            "path": f"source-archives/{self.bridge_bundle.name}",
            "sha256": bridge_sha,
            "complete_history": True,
            "head": self.bridge_head,
            "advertised_tip": {
                "revision": self.bridge_head,
                "name": "refs/heads/feature/cgl-landau-fluid",
            },
            "verified_revisions": sorted([self.f115_head, self.bridge_head]),
            "selected_as_current": False,
            "role": "retained-non-current-bridge",
        }
        predecessor = {
            "path": f"source-archives/{self.f116_bundle.name}",
            "sha256": sha256(self.f116_bundle),
            "complete_history": True,
            "head": self.f116_head,
            "advertised_tip": {"revision": self.f116_head, "name": "HEAD"},
            "verified_revisions": sorted(
                [self.f115_head, self.bridge_head, self.f116_head]
            ),
            "selected_as_current": False,
            "role": "retained-non-current-predecessor",
            "subject": git(
                self.repository, "show", "-s", "--format=%s", self.f116_head
            ).stdout.strip(),
        }
        generated = now()
        evidence = {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-evidence",
            "checkpoint": authority.CHECKPOINT,
            "execution_epoch": authority.EXECUTION_EPOCH,
            "generated_utc": generated,
            "scope": {
                "relationship": "current-source-selection-only-supersession",
                "summary": "Select the exact committed final tooling source without execution authority.",
                "preserves": authority.PRESERVES,
                "does_not_authorize": authority.DOES_NOT_AUTHORIZE,
            },
            "predecessor_authorities": {"historical_f116": self.f116_bindings},
            "implementation": {
                "publisher": {
                    "path": authority.PUBLISHER_RELATIVE.as_posix(),
                    "revision": self.final_head,
                    "sha256": sha256(self.publisher),
                    "mode": "0755",
                },
                "committed_tools": self.committed_tools(),
                "intermediate_36140_bundle": bridge,
                "predecessor_current_source_bundle": predecessor,
                "current_source_bundle": final,
            },
            "source_archive_catalog": {"before": before, "after": after},
            "authorization": self.authorization,
            "validation": authority.VALIDATION_CLAIMS,
            "publication_requirements": authority.PUBLICATION_REQUIREMENTS,
        }
        evidence_sha = write_json(self.evidence_candidate, evidence)
        verified = authority.review_verified(final)
        reviews = {}
        for key, kind, decision in (
            ("provenance", "provenance-security", "approved-for-publication"),
            ("plasma", "plasma-scientific-continuation", "approved"),
        ):
            review = {
                "schema_version": 1,
                "record_type": (
                    "stage-i-current-source-authority-supersession-independent-review"
                ),
                "checkpoint": authority.CHECKPOINT,
                "execution_epoch": authority.EXECUTION_EPOCH,
                "review_kind": kind,
                "decision": decision,
                "reviewed_candidate": {
                    "path": str(self.evidence_candidate),
                    "sha256": evidence_sha,
                },
                "published_f118": {
                    "path": str(self.root / authority.F118_PATHS["evidence"]),
                    "sha256": evidence_sha,
                },
                "reviewer": {
                    "agent_id": self.reviewer_ids[key],
                    "identity": f"fixture independent {key} reviewer",
                },
                "reviewed_utc": now(),
                "findings": ["Exact source-only authority and history bindings verified."],
                "limitations": [
                    "No prepare, submit, scheduler, or scientific authority.",
                    authority.INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION,
                ],
                "verified": verified,
            }
            path = self.provenance_candidate if key == "provenance" else self.plasma_candidate
            reviews[key] = write_json(path, review)
        audit = {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-publication-audit",
            "checkpoint": authority.CHECKPOINT,
            "execution_epoch": authority.EXECUTION_EPOCH,
            "published_utc": now(),
            "artifact": published_binding(
                self.root / authority.F118_PATHS["evidence"], evidence_sha, "0444"
            ),
            "independent_reviews": {
                "reviews_bind_exact_published_f118_sha256": evidence_sha,
                "provenance_security": published_binding(
                    self.root / authority.F118_PATHS["provenance_review"],
                    reviews["provenance"],
                    "0444",
                ),
                "plasma_scientific_continuation": published_binding(
                    self.root / authority.F118_PATHS["plasma_review"],
                    reviews["plasma"],
                    "0444",
                ),
            },
            "historical_f116_authority": self.f116_digests,
            "source_archive_catalog": {
                "readme": published_binding(
                    self.root / "source-archives/README.md", after["readme_sha256"], "0644"
                ),
                "sha256sums": published_binding(
                    self.root / "source-archives/SHA256SUMS",
                    after["sha256sums_sha256"],
                    "0644",
                ),
                "bridge_bundle": {
                    "path": str(self.bridge_bundle),
                    "sha256": bridge_sha,
                    "mode": "0644",
                    "links": 1,
                    "head": self.bridge_head,
                    "role": "retained-non-current-bridge",
                    "selected_as_current": False,
                },
                "predecessor_current_source_bundle": {
                    "path": str(self.f116_bundle),
                    "sha256": sha256(self.f116_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": self.f116_head,
                    "role": "retained-non-current-predecessor",
                    "selected_as_current": False,
                },
                "current_source_bundle": {
                    "path": str(self.final_target),
                    "sha256": final_sha,
                    "mode": "0644",
                    "links": 1,
                    "head": self.final_head,
                    "selected_as_current": True,
                },
                "corrupt_c7_absent_from_active_checksum_ledger": True,
                "sole_current_source_bundle": str(self.final_target),
            },
            "authority_and_enforcement": authority.AUTHORIZATION,
            "publication": (
                "recoverable-forward-transaction-with-publication-audit-commit-marker-"
                "under-stage-i-lock"
            ),
        }
        audit_sha = write_json(self.audit_candidate, audit)
        self.expected = {
            "publisher": sha256(self.publisher),
            "bundle": final_sha,
            "evidence": evidence_sha,
            "provenance_review": reviews["provenance"],
            "plasma_review": reviews["plasma"],
            "audit": audit_sha,
        }

    def replace_final_bundle(self, *revisions: str) -> None:
        self.final_bundle.unlink()
        git(self.repository, "bundle", "create", str(self.final_bundle), *revisions)
        self.final_bundle.chmod(0o644)
        self.refresh_candidates()

    def common_command(self) -> list[str]:
        return [
            str(self.publisher),
            "--root",
            str(self.root),
            "--allow-local-root",
            "--expected-publisher-sha256",
            self.expected["publisher"],
        ]

    def draft_evidence_command(self, output: Path, generated_utc: str) -> list[str]:
        return [
            *self.common_command(),
            "draft-evidence",
            "--bundle-candidate",
            str(self.final_bundle),
            "--expected-bundle-sha256",
            sha256(self.final_bundle),
            "--bridge-bundle",
            str(self.bridge_bundle),
            "--expected-bridge-sha256",
            sha256(self.bridge_bundle),
            "--bridge-head",
            self.bridge_head,
            "--evidence-candidate-output",
            str(output),
            "--generated-utc",
            generated_utc,
        ]

    def draft_audit_command(self, evidence: Path, output: Path,
                            published_utc: str) -> list[str]:
        return [
            *self.common_command(),
            "draft-audit",
            "--bundle-candidate",
            str(self.final_bundle),
            "--expected-bundle-sha256",
            sha256(self.final_bundle),
            "--evidence-candidate",
            str(evidence),
            "--expected-evidence-sha256",
            sha256(evidence),
            "--provenance-review-candidate",
            str(self.provenance_candidate),
            "--expected-provenance-review-sha256",
            sha256(self.provenance_candidate),
            "--plasma-review-candidate",
            str(self.plasma_candidate),
            "--expected-plasma-review-sha256",
            sha256(self.plasma_candidate),
            "--audit-candidate-output",
            str(output),
            "--published-utc",
            published_utc,
        ]

    def promote_command(self, interruption: str | None = None) -> list[str]:
        command = [
            *self.common_command(),
            "promote",
            "--bundle-candidate",
            str(self.final_bundle),
            "--expected-bundle-sha256",
            self.expected["bundle"],
            "--evidence-candidate",
            str(self.evidence_candidate),
            "--expected-evidence-sha256",
            self.expected["evidence"],
            "--provenance-review-candidate",
            str(self.provenance_candidate),
            "--expected-provenance-review-sha256",
            self.expected["provenance_review"],
            "--plasma-review-candidate",
            str(self.plasma_candidate),
            "--expected-plasma-review-sha256",
            self.expected["plasma_review"],
            "--audit-candidate",
            str(self.audit_candidate),
            "--expected-audit-sha256",
            self.expected["audit"],
        ]
        if interruption is not None:
            command += ["--simulate-interruption", interruption]
        return command

    def recover_command(self, interruption: str | None = None) -> list[str]:
        command = [
            *self.common_command(),
            "recover",
            "--expected-audit-sha256",
            self.expected["audit"],
        ]
        if interruption is not None:
            command += ["--simulate-interruption", interruption]
        return command

    def verify_command(self) -> list[str]:
        return [
            *self.common_command(),
            "verify",
            "--expected-audit-sha256",
            self.expected["audit"],
        ]

    @staticmethod
    def run(command: list[str]) -> subprocess.CompletedProcess:
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )

    def promote(self, interruption: str | None = None) -> subprocess.CompletedProcess:
        return self.run(self.promote_command(interruption))

    def recover(self, interruption: str | None = None) -> subprocess.CompletedProcess:
        return self.run(self.recover_command(interruption))

    def verify(self) -> subprocess.CompletedProcess:
        return self.run(self.verify_command())


@pytest.fixture
def campaign(tmp_path: Path) -> Campaign:
    return Campaign.create(tmp_path)


def assert_failed(completed: subprocess.CompletedProcess, match: str) -> None:
    assert completed.returncode != 0, completed.stdout
    assert match in completed.stderr, completed.stderr


def f118_public_member(campaign: Campaign, key: str) -> Path:
    members = {
        "bundle": campaign.final_target,
        "evidence": campaign.root / authority.F118_PATHS["evidence"],
        "provenance_review": campaign.root / authority.F118_PATHS["provenance_review"],
        "plasma_review": campaign.root / authority.F118_PATHS["plasma_review"],
        "readme": campaign.root / "source-archives/README.md",
        "sha256sums": campaign.root / "source-archives/SHA256SUMS",
        "audit": campaign.root / authority.F118_PATHS["publication_audit"],
    }
    return members[key]


def hostile_exchange_public_member(path: Path) -> None:
    """Exchange one public authority name for an invalid same-profile inode."""

    replacement = path.with_name(f".{path.name}.hostile-exchange")
    write_bytes(
        replacement,
        b"hostile exchanged public authority\n",
        stat.S_IMODE(path.stat().st_mode),
    )
    os.replace(replacement, path)


def direct_promote(campaign: Campaign) -> None:
    args = SimpleNamespace(
        bundle_candidate=campaign.final_bundle,
        expected_bundle_sha256=campaign.expected["bundle"],
        evidence_candidate=campaign.evidence_candidate,
        expected_evidence_sha256=campaign.expected["evidence"],
        provenance_review_candidate=campaign.provenance_candidate,
        expected_provenance_review_sha256=campaign.expected["provenance_review"],
        plasma_review_candidate=campaign.plasma_candidate,
        expected_plasma_review_sha256=campaign.expected["plasma_review"],
        audit_candidate=campaign.audit_candidate,
        expected_audit_sha256=campaign.expected["audit"],
        simulate_interruption=None,
    )
    layout = authority.root_layout(campaign.root)
    with authority.stage_i_lock(layout):
        authority.promote(
            args,
            campaign.root,
            campaign.repository,
            False,
            campaign.expected["publisher"],
            layout,
        )


def direct_verify(campaign: Campaign) -> None:
    layout = authority.root_layout(campaign.root)
    with authority.stage_i_lock(layout):
        authority.verify_promoted(
            campaign.root,
            campaign.repository,
            False,
            campaign.expected["publisher"],
            layout,
            campaign.expected["audit"],
        )


def bind_external_reviews(campaign: Campaign, evidence_path: Path) -> None:
    evidence = json.loads(evidence_path.read_text())
    evidence_sha256 = sha256(evidence_path)
    final = evidence["implementation"]["current_source_bundle"]
    for path in (campaign.provenance_candidate, campaign.plasma_candidate):
        review = json.loads(path.read_text())
        review["reviewed_candidate"] = {
            "path": str(evidence_path),
            "sha256": evidence_sha256,
        }
        review["published_f118"] = {
            "path": str(campaign.root / authority.F118_PATHS["evidence"]),
            "sha256": evidence_sha256,
        }
        review["reviewed_utc"] = now()
        review["verified"] = authority.review_verified(final)
        write_json(path, review)


def test_promotes_and_verifies_source_selection_only(campaign: Campaign) -> None:
    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0
    sums = (campaign.root / "source-archives/SHA256SUMS").read_text()
    assert sums.count(campaign.f115_bundle.name) == 1
    assert sums.count(campaign.bridge_bundle.name) == 1
    assert sums.count(campaign.f116_bundle.name) == 1
    assert sums.count(campaign.final_bundle.name) == 1
    assert authority.CORRUPT_C7_NAME not in sums
    assert campaign.f115_bundle.exists()
    assert campaign.bridge_bundle.exists()
    assert campaign.f116_bundle.exists()
    assert campaign.final_target.exists()
    audit = json.loads((campaign.root / authority.F118_PATHS["publication_audit"]).read_text())
    assert audit["authority_and_enforcement"] == authority.AUTHORIZATION
    assert audit["authority_and_enforcement"]["prepare_authorized"] is False
    assert audit["authority_and_enforcement"]["submit_authorized"] is False


def test_normal_promote_recover_verify_lifecycle_returns_success(campaign: Campaign) -> None:
    assert campaign.promote().returncode == 0
    assert campaign.recover().returncode == 0
    assert campaign.verify().returncode == 0


@pytest.mark.parametrize(
    "member",
    (
        "bundle",
        "evidence",
        "provenance_review",
        "plasma_review",
        "readme",
        "sha256sums",
    ),
)
def test_pre_audit_public_authority_lease_rejects_immediate_target_exchange(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch, member: str
) -> None:
    real_ensure = authority.ensure_direct_final_file
    audit_target = f118_public_member(campaign, "audit")
    mutated = False

    def exchange_immediately_before_audit_publication(
        payload: bytes, target: Path, expected: str, mode: int, label: str
    ) -> None:
        nonlocal mutated
        if target == audit_target and not mutated:
            hostile_exchange_public_member(f118_public_member(campaign, member))
            mutated = True
        real_ensure(payload, target, expected, mode, label)

    monkeypatch.setattr(
        authority, "ensure_direct_final_file", exchange_immediately_before_audit_publication
    )
    with pytest.raises(ValueError, match="public-authority lease|inode identity changed"):
        direct_promote(campaign)

    assert mutated
    assert not audit_target.exists()


def test_visible_audit_with_post_publication_drift_is_non_authorizing_and_not_repaired(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    real_ensure = authority.ensure_direct_final_file
    audit_target = f118_public_member(campaign, "audit")
    evidence_target = f118_public_member(campaign, "evidence")
    mutated = False

    def drift_immediately_after_audit_publication(
        payload: bytes, target: Path, expected: str, mode: int, label: str
    ) -> None:
        nonlocal mutated
        real_ensure(payload, target, expected, mode, label)
        if target == audit_target and not mutated:
            hostile_exchange_public_member(evidence_target)
            mutated = True

    monkeypatch.setattr(
        authority, "ensure_direct_final_file", drift_immediately_after_audit_publication
    )
    with pytest.raises(ValueError, match="public-authority lease|inode identity changed"):
        direct_promote(campaign)

    assert mutated
    assert sha256(audit_target) == campaign.expected["audit"]
    assert authority.publication_audit_committed(
        authority.root_layout(campaign.root), campaign.expected["audit"]
    )
    drifted_sha256 = sha256(evidence_target)
    assert drifted_sha256 != campaign.expected["evidence"]
    assert campaign.recover().returncode != 0
    assert campaign.verify().returncode != 0
    assert sha256(evidence_target) == drifted_sha256


@pytest.mark.parametrize(
    "member",
    (
        "bundle",
        "evidence",
        "provenance_review",
        "plasma_review",
        "readme",
        "sha256sums",
        "audit",
    ),
)
def test_verify_promoted_public_authority_lease_rejects_final_return_tail_exchange(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch, member: str
) -> None:
    assert campaign.promote().returncode == 0
    real_classify = authority.classify_non_authoritative_recovery_debris
    mutated = False

    def mutate_after_final_validation(layout: dict[str, Path]) -> None:
        nonlocal mutated
        real_classify(layout)
        hostile_exchange_public_member(f118_public_member(campaign, member))
        mutated = True

    monkeypatch.setattr(
        authority, "classify_non_authoritative_recovery_debris", mutate_after_final_validation
    )
    with pytest.raises(ValueError, match="public-authority lease|inode identity changed"):
        direct_verify(campaign)

    assert mutated
    assert f118_public_member(campaign, "audit").exists()
    assert campaign.verify().returncode != 0


@pytest.mark.parametrize("action", ("promote", "recover", "verify"))
def test_cli_lifecycle_final_return_lease_rejects_exchange_after_guards_release(
    campaign: Campaign,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    action: str,
) -> None:
    if action == "recover":
        assert_failed(campaign.promote("after-catalogs"), "simulated interruption")
    elif action == "verify":
        assert campaign.promote().returncode == 0

    command = getattr(campaign, f"{action}_command")()[1:]
    real_active = authority.active_f118_public_authority_lease
    mutated = False

    @contextmanager
    def mutate_at_cli_final_return(lease: authority.F118PublicAuthorityLease):
        nonlocal mutated
        if (
            not mutated
            and authority._ACTIVE_MUTATION_LOCK is None
            and authority._ACTIVE_CANONICAL_PUBLIC_NAMESPACE is None
        ):
            hostile_exchange_public_member(f118_public_member(campaign, "evidence"))
            mutated = True
        with real_active(lease) as retained:
            yield retained

    monkeypatch.setattr(authority, "active_f118_public_authority_lease", mutate_at_cli_final_return)
    monkeypatch.setattr(
        authority,
        "authenticate_self",
        lambda _argv: (
            campaign.publisher,
            campaign.repository,
            campaign.expected["publisher"],
        ),
    )

    with pytest.raises(ValueError, match="public-authority lease|inode identity changed"):
        authority.main(command)

    assert mutated
    assert authority._ACTIVE_MUTATION_LOCK is None
    assert authority._ACTIVE_CANONICAL_PUBLIC_NAMESPACE is None
    assert capsys.readouterr().out == ""
    assert f118_public_member(campaign, "audit").exists()
    assert campaign.verify().returncode != 0


def test_f118_preserves_exact_f116_four_part_authority_and_selected_bundle(
    campaign: Campaign,
) -> None:
    retained = {
        key: retained_tree_snapshot(campaign.root / relative)
        for key, relative in authority.F116_PATHS.items()
    }
    retained["selected_bundle"] = retained_tree_snapshot(campaign.f116_bundle)
    assert all(
        campaign.root / authority.F116_PATHS[key]
        != campaign.root / authority.F118_PATHS[key]
        for key in authority.F116_PATHS
    )

    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0

    assert {
        key: retained_tree_snapshot(campaign.root / relative)
        for key, relative in authority.F116_PATHS.items()
    } == {key: retained[key] for key in authority.F116_PATHS}
    assert retained_tree_snapshot(campaign.f116_bundle) == retained["selected_bundle"]
    assert all((campaign.root / relative).exists() for relative in authority.F118_PATHS.values())


def test_f118_ignores_and_preserves_retained_f116_transaction_namespace(
    campaign: Campaign,
) -> None:
    retained = (
        campaign.root / "accounting" / authority.F116_TRANSACTION_ROOT_NAME
        / f"2026-06-06T000000+0000-{'a' * 32}.staging"
    )
    retained.mkdir(parents=True, mode=0o700)
    write_bytes(retained / "journal.json", b"retained F116 transaction bytes\n", 0o600)
    before = retained_tree_snapshot(retained.parent)

    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0

    assert retained_tree_snapshot(retained.parent) == before
    assert (
        campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    ).exists()


def test_rejects_valid_foreign_archive_and_checksum_entry_after_f116(
    campaign: Campaign,
) -> None:
    archives = campaign.root / "source-archives"
    foreign = archives / "foreign-valid-history.bundle"
    write_bytes(foreign, b"foreign structurally valid history\n", 0o644)
    sums_path = archives / "SHA256SUMS"
    sums_path.write_bytes(
        sums_path.read_bytes() + f"{sha256(foreign)}  {foreign.name}\n".encode()
    )
    assert_failed(
        campaign.promote(),
        "live predecessor source-archive catalog differs from authenticated F116 catalog_after",
    )
    assert not campaign.final_target.exists()


def test_rejects_readme_only_drift_after_f116(campaign: Campaign) -> None:
    readme = campaign.root / "source-archives/README.md"
    readme.write_bytes(readme.read_bytes() + b"\nForeign but structurally valid note.\n")
    assert_failed(
        campaign.promote(),
        "live predecessor source-archive catalog differs from authenticated F116 catalog_after",
    )
    assert not campaign.final_target.exists()


@pytest.mark.parametrize(
    "point",
    [
        "after-staging",
        "after-bundle",
        "after-artifacts",
        "after-readme",
        "after-catalogs",
        "after-audit",
    ],
)
def test_every_publication_interruption_recovers(campaign: Campaign, point: str) -> None:
    assert_failed(campaign.promote(point), f"simulated interruption {point}")
    audit = campaign.root / authority.F118_PATHS["publication_audit"]
    assert audit.exists() is (point == "after-audit")
    if point == "after-audit":
        assert campaign.verify().returncode == 0
    else:
        assert campaign.verify().returncode != 0
    assert campaign.recover().returncode == 0
    assert campaign.verify().returncode == 0
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    retained = list(transactions.iterdir())
    assert any(path.name.endswith(authority.STAGING_TRANSACTION_SUFFIX) for path in retained)


@pytest.mark.parametrize(
    "point",
    [
        "during-staging-after-directory",
        "during-staging-after-first-payload",
        "during-staging-before-journal",
    ],
)
def test_incomplete_staging_interruption_fails_closed_and_remains_bounded(
    campaign: Campaign, point: str
) -> None:
    assert_failed(campaign.promote(point), f"simulated interruption {point}")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    retained = list(transactions.iterdir())
    assert len(retained) == 1
    assert retained[0].name.endswith(authority.STAGING_TRANSACTION_SUFFIX)
    assert not campaign.final_target.exists()
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()
    assert_failed(campaign.promote(), "requires operator disposition")
    assert list(transactions.iterdir()) == retained


def test_complete_staging_journal_is_cleanly_retryable(campaign: Campaign) -> None:
    assert_failed(
        campaign.promote("during-staging-after-journal"),
        "simulated interruption during-staging-after-journal",
    )
    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0


def test_repeated_crash_promotion_cleans_all_private_slots_before_audit_commit(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption after-staging")
    targets = (
        (campaign.final_target, "F118 final source bundle"),
        (campaign.root / authority.F118_PATHS["evidence"], "F118 evidence"),
        (
            campaign.root / authority.F118_PATHS["provenance_review"],
            "F118 provenance review",
        ),
        (
            campaign.root / authority.F118_PATHS["plasma_review"],
            "F118 plasma review",
        ),
        (campaign.root / "source-archives/README.md", "source-archive README"),
        (campaign.root / "source-archives/SHA256SUMS", "source-archive SHA256SUMS"),
        (
            campaign.root / authority.F118_PATHS["publication_audit"],
            "F118 publication audit",
        ),
    )
    private_paths = []
    for target, label in targets:
        names = authority.private_publication_names(target.name, label)
        for name, mode in zip(names, (0o000, 0o600), strict=True):
            private = target.parent / name
            write_bytes(private, b"repeated private writer interruption\n", mode)
            private_paths.append(private)

    assert_failed(campaign.promote("after-catalogs"), "simulated interruption after-catalogs")
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()
    assert any(os.path.lexists(path) for path in private_paths)

    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0
    assert not any(os.path.lexists(path) for path in private_paths)


def test_exact_audit_commit_never_cleans_late_private_slot(campaign: Campaign) -> None:
    assert campaign.promote().returncode == 0
    target = campaign.root / authority.F118_PATHS["evidence"]
    private = target.parent / authority.private_publication_names(
        target.name, "F118 evidence"
    )[0]
    write_bytes(private, b"late non-authoritative debris\n", 0o600)
    before = private.stat()

    assert campaign.promote().returncode == 0
    assert campaign.recover().returncode == 0
    after = private.stat()
    assert (after.st_dev, after.st_ino) == (before.st_dev, before.st_ino)
    assert authority.file_profile_binding(after) == authority.file_profile_binding(before)
    assert private.read_bytes() == b"late non-authoritative debris\n"


def test_precommit_private_cleanup_rejects_hardlinked_incomplete_slot(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-catalogs"), "simulated interruption after-catalogs")
    target = campaign.root / authority.F118_PATHS["evidence"]
    private = target.parent / authority.private_publication_names(
        target.name, "F118 evidence"
    )[0]
    alias = target.parent / "injected-private-alias"
    write_bytes(private, b"hardlinked private debris\n", 0o600)
    os.link(private, alias)

    assert_failed(campaign.promote(), "has 2 links, expected 1")
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()
    assert private.stat().st_nlink == 2
    assert alias.stat().st_nlink == 2


def test_abandoned_staging_cleanup_fails_closed_on_unexpected_entry(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-after-directory"),
        "simulated interruption during-staging-after-directory",
    )
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    staging = next(transactions.iterdir())
    write_bytes(staging / "unreviewed", b"unexpected\n", 0o600)
    assert_failed(campaign.promote(), "staging contains unexpected entries")
    assert not campaign.final_target.exists()
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()


def test_abandoned_staging_cleanup_rejects_visible_f118_target(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-after-directory"),
        "simulated interruption during-staging-after-directory",
    )
    visible = campaign.root / authority.F118_PATHS["evidence"]
    write_bytes(visible, b"uncommitted authority\n", 0o444)
    assert_failed(campaign.promote(), "requires operator disposition")
    assert visible.read_bytes() == b"uncommitted authority\n"
    assert not campaign.final_target.exists()


def test_partial_prejournal_payload_fails_closed_and_remains_bounded(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-after-first-payload"),
        "simulated interruption during-staging-after-first-payload",
    )
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    staging = next(transactions.iterdir())
    payload = staging / authority.TRANSACTION_PAYLOADS["bundle"][0]
    payload.chmod(0o600)
    payload.write_bytes(payload.read_bytes()[:128])
    assert_failed(campaign.promote(), "requires operator disposition")
    assert list(transactions.iterdir()) == [staging]


def test_mode_zero_prejournal_payload_fails_closed_and_remains_bounded(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-after-first-payload"),
        "simulated interruption during-staging-after-first-payload",
    )
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    staging = next(transactions.iterdir())
    payload = staging / authority.TRANSACTION_PAYLOADS["bundle"][0]
    payload.chmod(0o600)
    payload.write_bytes(payload.read_bytes()[:128])
    payload.chmod(0o000)

    assert_failed(campaign.promote(), "requires operator disposition")
    assert list(transactions.iterdir()) == [staging]


def test_incomplete_prepublication_staging_blocks_unbounded_retry(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-after-first-payload"),
        "simulated interruption during-staging-after-first-payload",
    )
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    incomplete = next(transactions.iterdir())
    assert_failed(campaign.recover(), "requires operator disposition")
    assert_failed(campaign.promote(), "only complete authenticated prior-attempt debris")
    assert incomplete.is_dir()
    assert list(transactions.iterdir()) == [incomplete]


def test_mode_zero_partial_journal_staging_fails_closed(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-before-journal"),
        "simulated interruption during-staging-before-journal",
    )
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    incomplete = next(transactions.iterdir())
    write_bytes(incomplete / "journal.json", b'{"partial":', 0o000)

    assert_failed(campaign.promote(), "journal requires operator disposition")
    assert stat.S_IMODE((incomplete / "journal.json").stat().st_mode) == 0o000
    assert list(transactions.iterdir()) == [incomplete]


def test_forensic_only_root_retirement_restores_concurrently_inserted_active_entry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "campaign"
    accounting = root / "accounting"
    transactions = accounting / authority.TRANSACTION_ROOT_NAME
    transactions.mkdir(parents=True)
    forensic = transactions / f".cgl-source-authority-retired-{'a' * 32}.forensic"
    write_bytes(forensic, b"retained forensic bytes\n", 0o600)
    layout = {
        "transactions": transactions,
        **{
            f"f118_{key}": root / relative
            for key, relative in authority.F118_PATHS.items()
        },
    }
    real_renameat2_between = authority.renameat2_between
    inserted = False

    def insert_active_entry_during_retirement(
        parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal inserted
        if (
            not inserted
            and source == transactions.name
            and "pre-publication forensic transaction root retirement" in label
        ):
            inserted = True
            descriptor = os.open(
                source,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent,
            )
            try:
                active = os.open(
                    "active-transaction",
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o600,
                    dir_fd=descriptor,
                )
                os.close(active)
            finally:
                os.close(descriptor)
        real_renameat2_between(parent, source, target_parent, target, flags, label)

    monkeypatch.setattr(
        authority, "renameat2_between", insert_active_entry_during_retirement
    )
    with pytest.raises(ValueError, match="contents changed during atomic retirement; restored"):
        authority.purge_forensic_only_prepublication_transaction_root(layout)

    assert inserted
    assert transactions.is_dir()
    assert (transactions / "active-transaction").is_file()
    assert forensic.read_bytes() == b"retained forensic bytes\n"
    assert not list(root.glob(".cgl-source-authority-retired-directory-*.forensic"))


def test_forensic_root_rollback_revalidates_destination_at_syscall(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    container = tmp_path / "container"
    transaction = container / "transaction"
    container.mkdir()
    transaction.mkdir()
    write_bytes(transaction / "retained", b"retained\n", 0o600)
    with authority.bound_directory(transaction, "transaction contents") as (
        contents,
        _,
    ):
        expected_contents = authority.directory_content_bindings(
            contents, "transaction contents"
        )

    real_renameat2_between = authority.renameat2_between
    inserted = False
    rollback_checked = False

    def insert_then_remove_rollback_authority(
        source_parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal inserted, rollback_checked
        if "contents rollback" in label:
            rollback_checked = True
            container.chmod(0o777)
            try:
                real_renameat2_between(
                    source_parent, source, target_parent, target, flags, label
                )
            finally:
                container.chmod(0o755)
            return
        real_renameat2_between(
            source_parent, source, target_parent, target, flags, label
        )
        if not inserted and label == "transaction retirement":
            inserted = True
            retired = os.open(
                target,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=target_parent,
            )
            try:
                active = os.open(
                    "active",
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o600,
                    dir_fd=retired,
                )
                os.close(active)
            finally:
                os.close(retired)

    monkeypatch.setattr(
        authority, "renameat2_between", insert_then_remove_rollback_authority
    )
    with authority.bound_directory(container, "transaction parent") as (parent, _):
        profile = os.stat(transaction.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="changed during atomic retirement"):
            authority.rmdir_bound_entry(
                parent,
                transaction.name,
                profile,
                "transaction",
                expected_contents=expected_contents,
            )

    assert inserted
    assert rollback_checked
    assert not transaction.exists()
    retired = list(tmp_path.glob(".cgl-source-authority-retired-directory-*.forensic"))
    assert len(retired) == 1
    assert (retired[0] / "retained").read_bytes() == b"retained\n"
    assert (retired[0] / "active").is_file()


def test_bundle_verification_rejects_second_open_substitution(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    valid = campaign.final_bundle.read_bytes()
    corrupt = bytearray(valid)
    corrupt[-1] ^= 0x01
    candidate = campaign.candidates / "corrupt-but-header-valid.bundle"
    substitute = campaign.candidates / "valid-substitute.bundle"
    write_bytes(candidate, bytes(corrupt), 0o644)
    write_bytes(substitute, valid, 0o644)
    real_open = authority.os.open
    candidate_opens = 0

    def substitute_second_open(path: object, flags: int, *args: object,
                               **kwargs: object) -> int:
        nonlocal candidate_opens
        if Path(path) == candidate:
            candidate_opens += 1
            if candidate_opens == 2:
                return real_open(substitute, flags, *args, **kwargs)
        return real_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(authority.os, "open", substitute_second_open)
    with pytest.raises(
        ValueError, match="fails git bundle verify|cannot be reconstructed in isolation"
    ):
        authority.stable_bundle_validation(
            campaign.repository,
            candidate,
            sha256(candidate),
            campaign.final_head,
            "HEAD",
            campaign.final_verified_revisions,
            "adversarial bundle",
        )
    assert candidate_opens == 1


def test_draft_evidence_is_deterministic_and_never_publishes(campaign: Campaign) -> None:
    generated_utc = now()
    first = campaign.candidates / "drafted-evidence-one.json"
    second = campaign.candidates / "drafted-evidence-two.json"
    readme_before = sha256(campaign.root / "source-archives/README.md")
    sums_before = sha256(campaign.root / "source-archives/SHA256SUMS")
    assert campaign.run(campaign.draft_evidence_command(first, generated_utc)).returncode == 0
    assert campaign.run(campaign.draft_evidence_command(second, generated_utc)).returncode == 0
    assert first.read_bytes() == second.read_bytes()
    assert first.stat().st_mode & 0o777 == 0o444
    assert not campaign.final_target.exists()
    assert not (campaign.root / authority.F118_PATHS["evidence"]).exists()
    assert sha256(campaign.root / "source-archives/README.md") == readme_before
    assert sha256(campaign.root / "source-archives/SHA256SUMS") == sums_before


def test_draft_evidence_fails_closed_while_staging_requires_disposition(
    campaign: Campaign,
) -> None:
    assert_failed(
        campaign.promote("during-staging-after-directory"),
        "simulated interruption during-staging-after-directory",
    )
    output = campaign.candidates / "blocked-draft.json"
    assert_failed(
        campaign.run(campaign.draft_evidence_command(output, now())),
        "requires operator disposition",
    )
    assert not output.exists()


def test_candidate_drafting_coexists_with_complete_stale_staging_without_mutation(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    retained_before = retained_tree_snapshot(transactions)
    evidence = campaign.candidates / "coexisting-drafted-evidence.json"
    audit = campaign.candidates / "coexisting-drafted-audit.json"

    assert campaign.run(campaign.draft_evidence_command(evidence, now())).returncode == 0
    assert campaign.run(
        campaign.draft_audit_command(campaign.evidence_candidate, audit, now())
    ).returncode == 0

    assert retained_tree_snapshot(transactions) == retained_before


def test_candidate_drafting_accepts_setgid_complete_stale_staging_without_mutation(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    transactions.chmod(0o2700)
    transaction.chmod(0o2700)
    retained_before = retained_tree_snapshot(transactions)
    output = campaign.candidates / "setgid-coexisting-draft.json"

    assert campaign.run(campaign.draft_evidence_command(output, now())).returncode == 0

    assert retained_tree_snapshot(transactions) == retained_before


@pytest.mark.parametrize(
    "corruption",
    [
        "bare-transaction",
        "retired-transaction",
        "installing-state",
        "committed-state",
        "missing-publisher-revision",
        "missing-publisher-file",
        "non-ancestor-publisher",
        "wrong-publisher-bytes",
        "wrong-publisher-mode",
        "rebound-f118-target",
        "rebound-bundle-target",
        "rebound-catalog-target",
        "mode-zero-journal",
        "mode-zero-payload",
        "malformed-journal",
        "malformed-transaction-path",
        "missing-payload",
        "symlink-payload",
        "unsafe-payload-mode",
        "externally-hardlinked-payload",
        "catalog-digest-rebinding",
        "partial-public-bundle",
        "dangling-public-bundle",
        "partial-public-artifact",
        "catalog-after",
        "catalog-mixed",
        "private-publication-slot",
        "catalog-temporary",
        "legacy-single-link-temporary",
        "catalog-predecessor-recovery",
        "catalog-atomic-recovery",
        "catalog-atomic-recovery-alternate",
        "live-catalog-drift",
    ],
)
def test_drafting_and_new_promotion_reject_non_inert_stale_transaction(
    campaign: Campaign, corruption: str
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    stale_audit = campaign.expected["audit"]
    campaign.reviewer_ids["provenance"] = "fixture-provenance-next-attempt"
    campaign.refresh_candidates()
    assert campaign.expected["audit"] != stale_audit

    journal = transaction / "journal.json"
    payload = transaction / authority.TRANSACTION_PAYLOADS["evidence"][0]
    if corruption == "bare-transaction":
        transaction = transaction.rename(
            transaction.with_name(transaction.name.removesuffix(
                authority.STAGING_TRANSACTION_SUFFIX
            ))
        )
    elif corruption == "retired-transaction":
        transaction = transaction.rename(
            transaction.with_name(
                transaction.name.removesuffix(authority.STAGING_TRANSACTION_SUFFIX)
                + authority.RETIRED_TRANSACTION_SUFFIX
            )
        )
    elif corruption in {"installing-state", "committed-state"}:
        retained_journal = json.loads(journal.read_text())
        retained_journal["state"] = corruption.removesuffix("-state")
        write_json(journal, retained_journal, 0o600)
    elif corruption == "missing-publisher-revision":
        retained_journal = json.loads(journal.read_text())
        retained_journal["publisher"]["revision"] = "0" * 40
        write_json(journal, retained_journal, 0o600)
    elif corruption == "missing-publisher-file":
        git(
            campaign.repository,
            "rm",
            "--cached",
            authority.PUBLISHER_RELATIVE.as_posix(),
        )
        git(campaign.repository, "commit", "-m", "missing publisher file")
        retained_journal = json.loads(journal.read_text())
        retained_journal["publisher"]["revision"] = git(
            campaign.repository, "rev-parse", "HEAD"
        ).stdout.strip()
        write_json(journal, retained_journal, 0o600)
    elif corruption == "non-ancestor-publisher":
        tree = git(campaign.repository, "rev-parse", "HEAD^{tree}").stdout.strip()
        unrelated = git(
            campaign.repository, "commit-tree", tree, "-m", "unrelated publisher"
        ).stdout.strip()
        retained_journal = json.loads(journal.read_text())
        retained_journal["publisher"]["revision"] = unrelated
        write_json(journal, retained_journal, 0o600)
    elif corruption == "wrong-publisher-bytes":
        original = campaign.publisher.read_bytes()
        write_bytes(campaign.publisher, b"committed wrong publisher bytes\n", 0o755)
        git(campaign.repository, "add", authority.PUBLISHER_RELATIVE.as_posix())
        git(campaign.repository, "commit", "-m", "wrong publisher bytes")
        write_bytes(campaign.publisher, original, 0o755)
        retained_journal = json.loads(journal.read_text())
        retained_journal["publisher"]["revision"] = git(
            campaign.repository, "rev-parse", "HEAD"
        ).stdout.strip()
        write_json(journal, retained_journal, 0o600)
    elif corruption == "wrong-publisher-mode":
        git(
            campaign.repository,
            "update-index",
            "--chmod=-x",
            authority.PUBLISHER_RELATIVE.as_posix(),
        )
        git(campaign.repository, "commit", "-m", "wrong publisher mode")
        retained_journal = json.loads(journal.read_text())
        retained_journal["publisher"]["revision"] = git(
            campaign.repository, "rev-parse", "HEAD"
        ).stdout.strip()
        write_json(journal, retained_journal, 0o600)
    elif corruption == "rebound-f118-target":
        retained_journal = json.loads(journal.read_text())
        retained_journal["targets"]["evidence"] = "accounting/rebound-f118-evidence"
        write_json(journal, retained_journal, 0o600)
    elif corruption == "rebound-bundle-target":
        retained_journal = json.loads(journal.read_text())
        retained_journal["targets"][
            "bundle"
        ] = "source-archives/athenak-feature-cgl-through-000000000.bundle"
        write_json(journal, retained_journal, 0o600)
    elif corruption == "rebound-catalog-target":
        retained_journal = json.loads(journal.read_text())
        retained_journal["targets"]["readme"] = "source-archives/REBOUND.md"
        write_json(journal, retained_journal, 0o600)
    elif corruption == "mode-zero-journal":
        journal.chmod(0o000)
    elif corruption == "mode-zero-payload":
        payload.chmod(0o000)
    elif corruption == "malformed-journal":
        write_json(journal, {"malformed": True}, 0o600)
    elif corruption == "malformed-transaction-path":
        transaction = transaction.rename(transactions / "malformed.staging")
    elif corruption == "missing-payload":
        payload.unlink()
    elif corruption == "symlink-payload":
        retained = campaign.candidates / "retained-staged-evidence"
        payload.rename(retained)
        payload.symlink_to(retained)
    elif corruption == "unsafe-payload-mode":
        payload.chmod(0o666)
    elif corruption == "externally-hardlinked-payload":
        os.link(payload, campaign.candidates / "external-staged-evidence-link")
    elif corruption == "catalog-digest-rebinding":
        retained_journal = json.loads(journal.read_text())
        retained_journal["catalog_after"]["readme_sha256"] = "0" * 64
        write_json(journal, retained_journal, 0o600)
    elif corruption == "partial-public-bundle":
        write_bytes(
            campaign.final_target,
            (transaction / authority.TRANSACTION_PAYLOADS["bundle"][0]).read_bytes(),
            0o644,
        )
    elif corruption == "dangling-public-bundle":
        campaign.final_target.symlink_to(campaign.candidates / "missing-public-bundle")
    elif corruption == "partial-public-artifact":
        write_bytes(
            campaign.root / authority.F118_PATHS["evidence"],
            (transaction / authority.TRANSACTION_PAYLOADS["evidence"][0]).read_bytes(),
            0o444,
        )
    elif corruption in {"catalog-after", "catalog-mixed"}:
        write_bytes(
            campaign.root / "source-archives/README.md",
            (transaction / authority.TRANSACTION_PAYLOADS["readme_after"][0]).read_bytes(),
            0o644,
        )
        if corruption == "catalog-after":
            write_bytes(
                campaign.root / "source-archives/SHA256SUMS",
                (
                    transaction / authority.TRANSACTION_PAYLOADS["sha256sums_after"][0]
                ).read_bytes(),
                0o644,
            )
    elif corruption == "private-publication-slot":
        target = campaign.root / authority.F118_PATHS["evidence"]
        write_bytes(
            target.parent / authority.private_publication_names(
                target.name, "stale F118 evidence"
            )[0],
            b"private publication debris\n",
            0o600,
        )
    else:
        transaction_id = authority.logical_transaction_id(transaction)
        target = campaign.root / "source-archives/README.md"
        if corruption == "catalog-temporary":
            name = f".{target.name}.{transaction_id}.tmp"
        elif corruption == "legacy-single-link-temporary":
            name = f".{target.name}.single-link.{transaction_id}.tmp"
        elif corruption == "catalog-predecessor-recovery":
            name = f".{target.name}.predecessor.{transaction_id}.tmp"
        elif corruption == "catalog-atomic-recovery":
            name = authority.atomic_recovery_names(target)[0]
        elif corruption == "catalog-atomic-recovery-alternate":
            name = authority.atomic_recovery_names(target)[1]
        else:
            write_bytes(target, b"unreviewed live catalog drift\n", 0o644)
            name = ""
        if not name:
            output = campaign.candidates / f"blocked-{corruption}-draft.json"
            assert campaign.run(campaign.draft_evidence_command(output, now())).returncode != 0
            assert not output.exists()
            assert campaign.promote().returncode != 0
            assert list(transactions.iterdir()) == [transaction]
            return
        write_bytes(target.parent / name, b"catalog recovery debris\n", 0o600)

    output = campaign.candidates / f"blocked-{corruption}-draft.json"
    assert campaign.run(campaign.draft_evidence_command(output, now())).returncode != 0
    assert not output.exists()
    assert campaign.promote().returncode != 0
    assert list(transactions.iterdir()) == [transaction]


@pytest.mark.parametrize(
    "point",
    ["after-bundle", "after-artifacts", "after-readme", "after-catalogs", "after-audit"],
)
def test_nonmatching_partial_publication_blocks_drafting_and_new_promotion(
    campaign: Campaign, point: str
) -> None:
    assert_failed(campaign.promote(point), f"simulated interruption {point}")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    retained_before = retained_tree_snapshot(transactions)
    stale_audit = campaign.expected["audit"]
    campaign.expected["audit"] = "0" * 64
    assert campaign.expected["audit"] != stale_audit

    output = campaign.candidates / f"blocked-{point}-draft.json"
    assert campaign.run(campaign.draft_evidence_command(output, now())).returncode != 0
    assert not output.exists()
    assert campaign.promote().returncode != 0
    assert retained_tree_snapshot(transactions) == retained_before


@pytest.mark.parametrize(
    ("target_key", "dangling"),
    [
        ("bundle", False),
        ("bundle", True),
        ("evidence", False),
        ("evidence", True),
        ("provenance_review", False),
        ("provenance_review", True),
        ("plasma_review", False),
        ("plasma_review", True),
        ("publication_audit", False),
        ("publication_audit", True),
    ],
)
def test_drafting_rejects_each_stale_public_target_and_dangling_symlink(
    campaign: Campaign, target_key: str, dangling: bool
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    campaign.expected["audit"] = "0" * 64
    if target_key == "bundle":
        target = campaign.final_target
        payload_key = "bundle"
        mode = 0o644
    else:
        target = campaign.root / authority.F118_PATHS[target_key]
        payload_key = "audit" if target_key == "publication_audit" else target_key
        mode = 0o444
    if dangling:
        target.symlink_to(campaign.candidates / f"missing-{target_key}")
    else:
        write_bytes(
            target,
            (transaction / authority.TRANSACTION_PAYLOADS[payload_key][0]).read_bytes(),
            mode,
        )

    output = campaign.candidates / f"blocked-{target_key}-{dangling}-draft.json"
    assert campaign.run(campaign.draft_evidence_command(output, now())).returncode != 0
    assert not output.exists()
    assert campaign.promote().returncode != 0
    assert list(transactions.iterdir()) == [transaction]


@pytest.mark.parametrize(
    ("target_kind", "recovery_kind", "dangling"),
    [
        ("artifact", "temporary", False),
        ("artifact", "temporary", True),
        ("artifact", "single-link", False),
        ("artifact", "single-link", True),
        ("catalog", "temporary", False),
        ("catalog", "temporary", True),
        ("catalog", "single-link", False),
        ("catalog", "single-link", True),
        ("catalog", "predecessor", False),
        ("catalog", "predecessor", True),
    ],
)
def test_drafting_and_new_promotion_reject_foreign_transaction_recovery_names(
    campaign: Campaign, target_kind: str, recovery_kind: str, dangling: bool
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    campaign.reviewer_ids["provenance"] = "fixture-provenance-next-attempt"
    campaign.refresh_candidates()

    target = (
        campaign.root / authority.F118_PATHS["evidence"]
        if target_kind == "artifact"
        else campaign.root / "source-archives/README.md"
    )
    foreign_transaction_id = f"2099-12-31T235959+0000-{'f' * 32}"
    assert foreign_transaction_id != authority.logical_transaction_id(transaction)
    infix = "" if recovery_kind == "temporary" else f"{recovery_kind}."
    recovery = target.parent / f".{target.name}.{infix}{foreign_transaction_id}.tmp"
    if dangling:
        recovery.symlink_to(campaign.candidates / f"missing-{target_kind}-{recovery_kind}")
    else:
        write_bytes(recovery, b"foreign transaction recovery debris\n", 0o600)

    output = campaign.candidates / (
        f"blocked-foreign-{target_kind}-{recovery_kind}-{dangling}-draft.json"
    )
    assert campaign.run(campaign.draft_evidence_command(output, now())).returncode != 0
    assert not output.exists()
    assert campaign.promote().returncode != 0
    assert list(transactions.iterdir()) == [transaction]
    assert os.path.lexists(recovery)


def test_transaction_recovery_scanner_retries_insertion_after_first_identity_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    transaction_id = f"2099-12-31T235959+0000-{'e' * 32}"
    recovery = tmp_path / f".{target.name}.{transaction_id}.tmp"
    real_bindings = authority.direct_child_identity_bindings
    inserted = False

    def insert_after_first_snapshot(
        parent: int, label: str
    ) -> tuple[tuple[str, tuple[int, ...]], ...] | None:
        nonlocal inserted
        bindings = real_bindings(parent, label)
        if not inserted:
            inserted = True
            descriptor = os.open(
                recovery.name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600,
                dir_fd=parent,
            )
            os.close(descriptor)
        return bindings

    monkeypatch.setattr(authority, "direct_child_identity_bindings", insert_after_first_snapshot)
    with pytest.raises(ValueError, match="retains forbidden names"):
        authority.require_transaction_recovery_names_absent(
            tmp_path, target.name, "hostile artifact recovery namespace"
        )

    assert inserted
    assert recovery.is_file()


@pytest.mark.parametrize(
    ("catalog", "name"),
    [
        (False, f".artifact.٢099-12-31T235959+0000-{'a' * 32}.tmp"),
        (False, f".artifact.2099-02-31T235959+0000-{'a' * 32}.tmp"),
        (False, f".artifact.2099-12-31T235959+0000-{'A' * 32}.tmp"),
        (False, f".artifact.2099-12-31T235959+0000-{'a' * 32}.tmp.extra"),
        (False, f".artifact.predecessor.2099-12-31T235959+0000-{'a' * 32}.tmp"),
        (True, f".README.md.predecessors.2099-12-31T235959+0000-{'a' * 32}.tmp"),
    ],
)
def test_transaction_recovery_scanner_accepts_unrelated_near_miss_names(
    tmp_path: Path, catalog: bool, name: str
) -> None:
    target = "README.md" if catalog else "artifact"
    near_miss = tmp_path / name
    near_miss.symlink_to(tmp_path / "missing-near-miss")

    authority.require_transaction_recovery_names_absent(
        tmp_path, target, "near-miss recovery namespace", catalog=catalog
    )

    assert near_miss.is_symlink()


def test_exact_audit_matching_partial_transaction_resumes_under_strict_classifier(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-artifacts"), "simulated interruption after-artifacts")
    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0


def test_drafted_evidence_reviews_and_audit_are_operable_end_to_end(
    campaign: Campaign,
) -> None:
    evidence = campaign.candidates / "publisher-drafted-evidence.json"
    audit_one = campaign.candidates / "publisher-drafted-audit-one.json"
    audit_two = campaign.candidates / "publisher-drafted-audit-two.json"
    assert campaign.run(campaign.draft_evidence_command(evidence, now())).returncode == 0
    bind_external_reviews(campaign, evidence)
    published_utc = now()
    assert campaign.run(
        campaign.draft_audit_command(evidence, audit_one, published_utc)
    ).returncode == 0
    assert campaign.run(
        campaign.draft_audit_command(evidence, audit_two, published_utc)
    ).returncode == 0
    assert audit_one.read_bytes() == audit_two.read_bytes()
    assert audit_one.stat().st_mode & 0o777 == 0o444
    assert not campaign.final_target.exists()
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()

    command = campaign.promote_command()
    replacements = {
        "--evidence-candidate": str(evidence),
        "--expected-evidence-sha256": sha256(evidence),
        "--provenance-review-candidate": str(campaign.provenance_candidate),
        "--expected-provenance-review-sha256": sha256(campaign.provenance_candidate),
        "--plasma-review-candidate": str(campaign.plasma_candidate),
        "--expected-plasma-review-sha256": sha256(campaign.plasma_candidate),
        "--audit-candidate": str(audit_one),
        "--expected-audit-sha256": sha256(audit_one),
    }
    for option, value in replacements.items():
        command[command.index(option) + 1] = value
    assert campaign.run(command).returncode == 0
    verify = campaign.verify_command()
    verify[-1] = sha256(audit_one)
    assert campaign.run(verify).returncode == 0


def test_drafted_candidate_promotion_coexists_with_unchanged_complete_stale_transaction(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    stale = next(transactions.iterdir())
    stale_audit = campaign.expected["audit"]
    stale_before = retained_tree_snapshot(stale)
    evidence = campaign.candidates / "coexisting-end-to-end-evidence.json"
    audit = campaign.candidates / "coexisting-end-to-end-audit.json"

    assert campaign.run(campaign.draft_evidence_command(evidence, now())).returncode == 0
    bind_external_reviews(campaign, evidence)
    assert campaign.run(campaign.draft_audit_command(evidence, audit, now())).returncode == 0
    assert sha256(audit) != stale_audit

    command = campaign.promote_command()
    replacements = {
        "--evidence-candidate": str(evidence),
        "--expected-evidence-sha256": sha256(evidence),
        "--provenance-review-candidate": str(campaign.provenance_candidate),
        "--expected-provenance-review-sha256": sha256(campaign.provenance_candidate),
        "--plasma-review-candidate": str(campaign.plasma_candidate),
        "--expected-plasma-review-sha256": sha256(campaign.plasma_candidate),
        "--audit-candidate": str(audit),
        "--expected-audit-sha256": sha256(audit),
    }
    for option, value in replacements.items():
        command[command.index(option) + 1] = value
    assert campaign.run(command).returncode == 0
    verify = campaign.verify_command()
    verify[-1] = sha256(audit)
    assert campaign.run(verify).returncode == 0

    assert retained_tree_snapshot(stale) == stale_before
    retained_audits = {
        json.loads((transaction / "journal.json").read_text())["expected"]["audit"]
        for transaction in transactions.iterdir()
    }
    assert retained_audits == {stale_audit, sha256(audit)}


def test_draft_audit_rejects_review_for_different_evidence(campaign: Campaign) -> None:
    evidence = campaign.candidates / "publisher-drafted-evidence.json"
    audit = campaign.candidates / "publisher-drafted-audit.json"
    assert campaign.run(campaign.draft_evidence_command(evidence, now())).returncode == 0
    bind_external_reviews(campaign, evidence)
    review = json.loads(campaign.plasma_candidate.read_text())
    review["reviewed_candidate"]["sha256"] = "0" * 64
    write_json(campaign.plasma_candidate, review)
    assert_failed(
        campaign.run(campaign.draft_audit_command(evidence, audit, now())),
        "identity or authorization boundary differs",
    )
    assert not audit.exists()
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()


def test_draft_output_cannot_enter_campaign_root_or_overwrite(
    campaign: Campaign,
) -> None:
    inside = campaign.root / "accounting/draft.json"
    assert_failed(
        campaign.run(campaign.draft_evidence_command(inside, now())),
        "must remain outside the campaign root",
    )
    outside = campaign.candidates / "existing-draft.json"
    write_bytes(outside, b"retain\n", 0o600)
    assert_failed(
        campaign.run(campaign.draft_evidence_command(outside, now())),
        "already exists",
    )
    assert outside.read_bytes() == b"retain\n"


def test_rejects_uncommitted_required_tool(campaign: Campaign) -> None:
    tool = campaign.repository / "scripts/frontier/cgl_lf_stage_i_qualification.py"
    tool.write_text("uncommitted drift\n")
    assert_failed(campaign.promote(), "required source-authority tool is not committed")


def test_rejects_historical_f115_chain_drift(campaign: Campaign) -> None:
    path = campaign.root / authority.F115_PATHS["evidence"]
    path.chmod(0o644)
    path.write_bytes(path.read_bytes() + b" ")
    path.chmod(0o444)
    assert_failed(campaign.promote(), "historical F115 evidence checksum differs")


@pytest.mark.parametrize("key", sorted(authority.F116_PATHS))
def test_rejects_historical_f116_four_part_drift(campaign: Campaign, key: str) -> None:
    path = campaign.root / authority.F116_PATHS[key]
    path.chmod(0o644)
    path.write_bytes(path.read_bytes() + b" ")
    path.chmod(0o444)
    assert_failed(campaign.promote(), f"historical F116 {key} checksum differs")


def test_rejects_historical_f116_selected_bundle_drift(campaign: Campaign) -> None:
    campaign.f116_bundle.write_bytes(campaign.f116_bundle.read_bytes() + b" ")
    assert_failed(campaign.promote(), "historical F116 current source bundle checksum differs")


def test_rejects_catalog_that_drops_f116_predecessor_bundle(campaign: Campaign) -> None:
    sums = campaign.root / "source-archives/SHA256SUMS"
    sums.write_text(
        "".join(
            line
            for line in sums.read_text().splitlines(keepends=True)
            if campaign.f116_bundle.name not in line
        )
    )
    assert_failed(
        campaign.promote(),
        "live predecessor source-archive catalog differs from authenticated F116 catalog_after",
    )


def test_rejects_final_bundle_that_advertises_branch_ref(campaign: Campaign) -> None:
    campaign.replace_final_bundle("refs/heads/feature/cgl-landau-fluid")
    assert_failed(campaign.promote(), "must advertise exactly")


def test_rejects_final_bundle_with_prerequisite_history(campaign: Campaign) -> None:
    campaign.replace_final_bundle("HEAD", f"^{campaign.bridge_head}")
    assert_failed(campaign.promote(), "must be self-contained without prerequisites")


def test_rejects_missing_required_bridge_revision(campaign: Campaign) -> None:
    campaign.final_verified_revisions = [campaign.final_head]
    campaign.refresh_candidates()
    assert_failed(campaign.promote(), "verified revisions omit required history")


def test_rejects_incorrect_final_head_subject(campaign: Campaign) -> None:
    campaign.subject_override = "incorrect final HEAD subject"
    campaign.refresh_candidates()
    assert_failed(campaign.promote(), "final HEAD subject differs")


def test_rejects_nonindependent_reviews(campaign: Campaign) -> None:
    campaign.reviewer_ids["plasma"] = campaign.reviewer_ids["provenance"]
    campaign.refresh_candidates()
    assert_failed(campaign.promote(), "do not have distinct reviewers")


def test_rejects_authorization_broadening(campaign: Campaign) -> None:
    campaign.authorization["submit_authorized"] = True
    campaign.refresh_candidates()
    assert_failed(campaign.promote(), "broadens source-selection-only authority")


def test_rejects_corrupt_c7_reentry_in_active_catalog(campaign: Campaign) -> None:
    corrupt = campaign.root / "source-archives" / authority.CORRUPT_C7_NAME
    sums = campaign.root / "source-archives/SHA256SUMS"
    sums.write_text(sums.read_text() + f"{sha256(corrupt)}  {corrupt.name}\n")
    campaign.refresh_candidates()
    assert_failed(
        campaign.promote(),
        "live predecessor source-archive catalog differs from authenticated F116 catalog_after",
    )


def test_rejects_existing_final_target(campaign: Campaign) -> None:
    write_bytes(campaign.final_target, b"collision\n", 0o644)
    assert_failed(campaign.promote(), "F118 target already exists")


def test_repeated_promotion_resumes_exact_retained_transaction(campaign: Campaign) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    retained = next(transactions.iterdir())
    retained_identity = retained.stat().st_dev, retained.stat().st_ino
    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0
    assert (retained.stat().st_dev, retained.stat().st_ino) == retained_identity


def test_new_promotion_coexists_with_complete_stale_transaction_without_mutation(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    stale = next(transactions.iterdir())
    stale_audit = campaign.expected["audit"]
    stale_before = retained_tree_snapshot(stale)

    campaign.reviewer_ids["provenance"] = "fixture-provenance-next-attempt"
    campaign.refresh_candidates()
    assert campaign.expected["audit"] != stale_audit

    assert campaign.promote().returncode == 0
    assert campaign.verify().returncode == 0
    assert retained_tree_snapshot(stale) == stale_before
    retained_audits = {
        json.loads((transaction / "journal.json").read_text())["expected"]["audit"]
        for transaction in transactions.iterdir()
    }
    assert retained_audits == {stale_audit, campaign.expected["audit"]}


def test_new_promotion_rejects_nonmatching_stale_transaction_payload_tampering(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    stale = next(transactions.iterdir())
    stale_audit = campaign.expected["audit"]
    payload = stale / authority.TRANSACTION_PAYLOADS["evidence"][0]
    payload.chmod(0o644)
    payload.write_bytes(payload.read_bytes() + b" ")
    payload.chmod(0o444)

    campaign.reviewer_ids["provenance"] = "fixture-provenance-next-attempt"
    campaign.refresh_candidates()
    assert campaign.expected["audit"] != stale_audit

    assert_failed(campaign.promote(), "checksum differs")
    assert list(transactions.iterdir()) == [stale]
    assert not campaign.final_target.exists()
    assert not (campaign.root / authority.F118_PATHS["publication_audit"]).exists()


def test_recovery_rejects_transaction_payload_tampering(campaign: Campaign) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    payload = transaction / authority.TRANSACTION_PAYLOADS["evidence"][0]
    payload.chmod(0o644)
    payload.write_bytes(payload.read_bytes() + b" ")
    payload.chmod(0o444)
    assert_failed(campaign.recover(), "checksum differs")


def test_recovery_rejects_journal_catalog_rebinding(campaign: Campaign) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    journal_path = transaction / "journal.json"
    journal = json.loads(journal_path.read_text())
    journal["catalog_before"]["readme_sha256"] = "0" * 64
    write_json(journal_path, journal, 0o600)
    assert_failed(campaign.recover(), "journal catalog bindings differ from payloads")


def test_recovery_rejects_interrupted_legacy_atomic_journal_replace(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    journal = transaction / "journal.json"
    temporary = transaction / f".journal.json.{os.getpid()}.{'a' * 32}.tmp"
    journal.rename(temporary)
    assert_failed(campaign.recover(), "incomplete source-authority staging")
    assert temporary.exists()
    assert not journal.exists()


def test_recovery_rejects_legacy_mode_zero_orphan_atomic_journal_temporary(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    temporary = transaction / f".journal.json.{os.getpid()}.{'a' * 32}.tmp"
    write_bytes(temporary, b"partial journal\n", 0o000)

    assert_failed(campaign.recover(), "mode is 0000, expected 0600")


def test_atomic_journal_exchange_mutation_recovers_exact_predecessor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transaction = tmp_path / "transaction"
    transaction.mkdir(mode=0o700)
    journal = transaction / "journal.json"
    predecessor = canonical_json({"state": "staged"})
    replacement = canonical_json({"state": "installing"})
    write_bytes(journal, predecessor, 0o600)
    real_renameat2 = authority.renameat2
    mutated = False

    def mutate_journal_after_exchange(
        parent: int, source: str, target: str, flags: int, label: str
    ) -> None:
        nonlocal mutated
        real_renameat2(parent, source, target, flags, label)
        if not mutated and flags == authority.RENAME_EXCHANGE:
            mutated = True
            drift_bound_file(
                parent,
                target,
                "content",
                b"X" * len(replacement),
                0o600,
            )

    monkeypatch.setattr(authority, "renameat2", mutate_journal_after_exchange)
    with pytest.raises(ValueError, match="profile changed|checksum differs"):
        authority.atomic_write(journal, replacement, 0o600)

    assert mutated
    assert journal.read_bytes() == b"X" * len(replacement)
    assert (transaction / authority.JOURNAL_RECOVERY_NAME).read_bytes() == predecessor

    monkeypatch.setattr(authority, "renameat2", real_renameat2)
    authority.recover_atomic_journal(transaction)

    assert journal.read_bytes() == predecessor
    assert_retained_journal_recovery(transaction, predecessor)
    assert not any(authority.JOURNAL_TEMP_RE.fullmatch(path.name) for path in transaction.iterdir())


def test_atomic_journal_recovers_after_random_predecessor_retirement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transaction = tmp_path / "transaction"
    transaction.mkdir(mode=0o700)
    journal = transaction / "journal.json"
    predecessor = canonical_json({"state": "staged"})
    replacement = canonical_json({"state": "installing"})
    forged = b"X" * len(replacement)
    write_bytes(journal, predecessor, 0o600)
    real_retire = authority.unlink_bound_entry
    real_fsync = authority.fsync_descriptor
    random_predecessor_retired = False
    mutated = False

    def observe_random_predecessor_retirement(
        parent: int, name: str, expected: os.stat_result, label: str
    ) -> None:
        nonlocal random_predecessor_retired
        real_retire(parent, name, expected, label)
        if authority.JOURNAL_TEMP_RE.fullmatch(name):
            random_predecessor_retired = True

    def mutate_after_random_predecessor_retirement(descriptor: int) -> None:
        nonlocal mutated
        real_fsync(descriptor)
        if (
            random_predecessor_retired
            and not mutated
            and journal.read_bytes() == replacement
        ):
            mutated = True
            drift_bound_file(descriptor, journal.name, "content", forged, 0o600)

    monkeypatch.setattr(authority, "unlink_bound_entry", observe_random_predecessor_retirement)
    monkeypatch.setattr(
        authority, "fsync_descriptor", mutate_after_random_predecessor_retirement
    )
    with pytest.raises(ValueError, match="profile changed|checksum differs"):
        authority.atomic_write(journal, replacement, 0o600)

    assert random_predecessor_retired
    assert mutated
    assert journal.read_bytes() == forged
    assert (transaction / authority.JOURNAL_RECOVERY_NAME).read_bytes() == predecessor
    assert not any(authority.JOURNAL_TEMP_RE.fullmatch(path.name) for path in transaction.iterdir())

    monkeypatch.setattr(authority, "unlink_bound_entry", real_retire)
    monkeypatch.setattr(authority, "fsync_descriptor", real_fsync)
    authority.recover_atomic_journal(transaction)

    assert journal.read_bytes() == predecessor
    assert_retained_journal_recovery(transaction, predecessor)


def test_atomic_journal_recovery_survives_public_forgery_during_orphan_retirement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transaction = tmp_path / "transaction"
    transaction.mkdir(mode=0o700)
    journal = transaction / "journal.json"
    payload = canonical_json({"state": "staged"})
    forged = b"X" * len(payload)
    orphan = transaction / f".journal.json.{os.getpid()}.{'a' * 32}.tmp"
    write_bytes(journal, payload, 0o600)
    write_bytes(orphan, payload, 0o600)
    real_retire = authority.unlink_bound_entry
    mutated = False

    def forge_public_journal_after_orphan_retirement(
        parent: int, name: str, expected: os.stat_result, label: str
    ) -> None:
        nonlocal mutated
        real_retire(parent, name, expected, label)
        if authority.JOURNAL_TEMP_RE.fullmatch(name) and not mutated:
            mutated = True
            drift_bound_file(parent, journal.name, "content", forged, 0o600)

    monkeypatch.setattr(
        authority, "unlink_bound_entry", forge_public_journal_after_orphan_retirement
    )
    with pytest.raises(ValueError, match="profile changed|checksum differs"):
        authority.recover_atomic_journal(transaction)

    assert mutated
    assert journal.read_bytes() == forged
    assert_retained_journal_recovery(transaction, payload)

    monkeypatch.setattr(authority, "unlink_bound_entry", real_retire)
    authority.recover_atomic_journal(transaction)

    assert journal.read_bytes() == payload
    assert_retained_journal_recovery(transaction, payload)


def test_incomplete_sole_atomic_journal_temporary_is_retired_fail_closed(
    tmp_path: Path,
) -> None:
    transaction = tmp_path / "transaction"
    transaction.mkdir(mode=0o700)
    temporary = transaction / f".journal.json.{os.getpid()}.{'a' * 32}.tmp"
    write_bytes(temporary, b"partial journal\n", 0o000)
    retained_before = set(tmp_path.glob(".cgl-source-authority-retired-*.forensic"))

    with pytest.raises(ValueError, match="lacks a complete recoverable journal"):
        authority.recover_atomic_journal(transaction)

    retained_after = set(tmp_path.glob(".cgl-source-authority-retired-*.forensic"))
    assert len(retained_after - retained_before) == 1
    assert not temporary.exists()
    assert not (transaction / "journal.json").exists()


def test_orphan_atomic_journal_temporary_substitution_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transaction = tmp_path / "transaction"
    transaction.mkdir(mode=0o700)
    journal = transaction / "journal.json"
    temporary = transaction / f".journal.json.{os.getpid()}.{'a' * 32}.tmp"
    substitute = transaction / "substitute"
    escaped = transaction / "escaped"
    write_json(journal, {}, 0o600)
    write_bytes(temporary, b"partial journal\n", 0o000)
    write_bytes(substitute, b"raced temporary\n", 0o600)
    real_renameat2_between = authority.renameat2_between
    raced = False

    def substitute_before_retirement(
        parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal raced
        if not raced and source == temporary.name and "retirement" in label:
            raced = True
            os.rename(source, escaped.name, src_dir_fd=parent, dst_dir_fd=parent)
            os.rename(substitute.name, source, src_dir_fd=parent, dst_dir_fd=parent)
        real_renameat2_between(parent, source, target_parent, target, flags, label)

    monkeypatch.setattr(authority, "renameat2_between", substitute_before_retirement)
    with pytest.raises(ValueError, match="changed during atomic retirement"):
        authority.recover_atomic_journal(transaction)

    assert raced
    assert escaped.stat().st_size == len(b"partial journal\n")
    assert not temporary.exists()
    assert b"raced temporary\n" in [
        path.read_bytes()
        for path in tmp_path.glob(".cgl-source-authority-retired-*.forensic")
    ]


def test_recovery_never_repairs_state_beneath_visible_audit_marker(
    campaign: Campaign,
) -> None:
    assert_failed(campaign.promote("after-staging"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    audit_target = campaign.root / authority.F118_PATHS["publication_audit"]
    write_bytes(
        audit_target,
        (transaction / authority.TRANSACTION_PAYLOADS["audit"][0]).read_bytes(),
        0o444,
    )
    readme_before = sha256(campaign.root / "source-archives/README.md")
    assert campaign.recover().returncode != 0
    assert sha256(campaign.root / "source-archives/README.md") == readme_before
    assert not campaign.final_target.exists()


@pytest.mark.parametrize("resume", ["promote", "recover"])
@pytest.mark.parametrize(
    "state",
    ["mode-zero-partial", "owner-only-partial", "owner-only-exact", "linked-private"],
)
def test_high_level_resume_completes_partial_audit_marker(
    campaign: Campaign, resume: str, state: str
) -> None:
    assert_failed(campaign.promote("after-catalogs"), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    audit_target = campaign.root / authority.F118_PATHS["publication_audit"]
    if state == "mode-zero-partial":
        write_bytes(audit_target, b'{"partial":', 0o000)
    elif state == "owner-only-partial":
        write_bytes(audit_target, b'{"partial":', 0o600)
    elif state == "owner-only-exact":
        write_bytes(
            audit_target,
            (transaction / authority.TRANSACTION_PAYLOADS["audit"][0]).read_bytes(),
            0o600,
        )
    else:
        private = audit_target.parent / authority.private_publication_names(
            audit_target.name, "F118 publication audit"
        )[0]
        write_bytes(
            private,
            (transaction / authority.TRANSACTION_PAYLOADS["audit"][0]).read_bytes(),
            0o444,
        )
        os.link(private, audit_target)
        assert audit_target.stat().st_nlink == 2

    assert campaign.verify().returncode != 0
    assert getattr(campaign, resume)().returncode == 0
    assert campaign.verify().returncode == 0
    assert sha256(audit_target) == campaign.expected["audit"]
    assert stat.S_IMODE(audit_target.stat().st_mode) == 0o444
    assert audit_target.stat().st_nlink == 1


@pytest.mark.parametrize(
    ("target_name", "payload_key", "interruption", "resume", "partial_mode"),
    [
        ("README.md", "readme_after", "after-artifacts", "recover", 0o000),
        ("SHA256SUMS", "sha256sums_after", "after-readme", "promote", 0o600),
    ],
)
def test_high_level_resume_completes_owner_only_partial_catalog(
    campaign: Campaign,
    target_name: str,
    payload_key: str,
    interruption: str,
    resume: str,
    partial_mode: int,
) -> None:
    assert_failed(campaign.promote(interruption), "simulated interruption")
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    transaction = next(transactions.iterdir())
    target = campaign.root / "source-archives" / target_name
    expected = (
        transaction / authority.TRANSACTION_PAYLOADS[payload_key][0]
    ).read_bytes()
    write_bytes(target, expected[:31], partial_mode)

    assert getattr(campaign, resume)().returncode == 0
    assert campaign.verify().returncode == 0
    assert target.read_bytes() == expected
    assert stat.S_IMODE(target.stat().st_mode) == 0o644


def test_post_commit_retained_staging_is_non_authoritative_recovery_debris(
    campaign: Campaign,
) -> None:
    assert campaign.promote().returncode == 0
    transactions = campaign.root / "accounting" / authority.TRANSACTION_ROOT_NAME
    retained = next(transactions.iterdir())
    (retained / "journal.json").chmod(0o000)
    (retained / "untrusted-link").symlink_to(campaign.evidence_candidate)

    assert campaign.verify().returncode == 0
    assert campaign.promote().returncode == 0
    assert campaign.recover().returncode == 0


def test_recovery_rejects_forged_final_mode_audit_marker(campaign: Campaign) -> None:
    assert_failed(campaign.promote("after-catalogs"), "simulated interruption")
    audit_target = campaign.root / authority.F118_PATHS["publication_audit"]
    write_bytes(audit_target, b'{"forged": true}\n', 0o444)

    assert_failed(campaign.recover(), "checksum differs")
    assert audit_target.read_bytes() == b'{"forged": true}\n'


def test_rejects_symbolic_link_candidate(campaign: Campaign) -> None:
    retained = campaign.evidence_candidate.with_suffix(".retained")
    campaign.evidence_candidate.rename(retained)
    campaign.evidence_candidate.symlink_to(retained)
    assert_failed(campaign.promote(), "path contains a symbolic link")


def test_rejects_duplicate_security_sensitive_option(campaign: Campaign) -> None:
    command = campaign.promote_command()
    command += ["--evidence-candidate", str(campaign.evidence_candidate)]
    assert_failed(campaign.run(command), "--evidence-candidate must not be supplied more than once")


def test_rejects_equals_form_security_sensitive_option(campaign: Campaign) -> None:
    command = campaign.promote_command()
    index = command.index("--evidence-candidate")
    command[index:index + 2] = [f"--evidence-candidate={campaign.evidence_candidate}"]
    assert_failed(campaign.run(command), "--evidence-candidate must use a separate exact value")


def test_rejects_stage_i_lock_contention(campaign: Campaign) -> None:
    lock = campaign.root / f".mks24_stage_i_{authority.EXECUTION_EPOCH_SLUG}.lock"
    with lock.open("r+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert_failed(campaign.promote(), "another Stage I mutation holds")


@pytest.mark.parametrize("initial", ["absent", "owner-only-partial"])
def test_direct_final_file_publication_never_renames(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, initial: str
) -> None:
    target = tmp_path / "artifact"
    payload = b"immutable reviewed bytes\n"
    digest = hashlib.sha256(payload).hexdigest()
    if initial == "owner-only-partial":
        write_bytes(target, b"partial", 0o600)

    def forbidden(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("direct-final publication must not rename")

    monkeypatch.setattr(authority, "renameat2", forbidden)
    monkeypatch.setattr(authority, "renameat2_between", forbidden)
    monkeypatch.setattr(authority.os, "rename", forbidden)
    authority.ensure_direct_final_file(payload, target, digest, 0o444, "fixture artifact")

    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    assert target.stat().st_nlink == 1


def test_direct_final_hardlink_injection_never_mutates_public_inode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    alias = tmp_path / "injected-alias"
    payload = b"immutable reviewed bytes\n"
    digest = hashlib.sha256(payload).hexdigest()
    real_link = authority.os.link
    real_write = authority.os.write
    real_fchmod = authority.os.fchmod
    published = False

    def inject_alias(source: str, destination: str, *args: object, **kwargs: object) -> None:
        nonlocal published
        real_link(source, destination, *args, **kwargs)
        if destination == target.name:
            parent = kwargs["dst_dir_fd"]
            real_link(
                destination,
                alias.name,
                src_dir_fd=parent,
                dst_dir_fd=parent,
                follow_symlinks=False,
            )
            published = True

    def reject_public_write(descriptor: int, retained: bytes) -> int:
        if published:
            raise AssertionError("public inode was written after publication")
        return real_write(descriptor, retained)

    def reject_public_fchmod(descriptor: int, mode: int) -> None:
        if published:
            raise AssertionError("public inode was chmodded after publication")
        real_fchmod(descriptor, mode)

    monkeypatch.setattr(authority.os, "link", inject_alias)
    monkeypatch.setattr(authority.os, "write", reject_public_write)
    monkeypatch.setattr(authority.os, "fchmod", reject_public_fchmod)
    with pytest.raises(ValueError, match="unsafe link profile"):
        authority.ensure_direct_final_file(payload, target, digest, 0o444, "fixture artifact")

    assert published
    assert target.read_bytes() == payload
    assert alias.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    assert target.stat().st_nlink == 3


def test_direct_final_file_rejects_forged_final_mode_target(tmp_path: Path) -> None:
    target = tmp_path / "artifact"
    payload = b"immutable reviewed bytes\n"
    write_bytes(target, b"forged bytes\n", 0o444)

    with pytest.raises(ValueError, match="cannot recover from mode 0444"):
        authority.ensure_direct_final_file(
            payload,
            target,
            hashlib.sha256(payload).hexdigest(),
            0o444,
            "fixture artifact",
        )
    assert target.read_bytes() == b"forged bytes\n"


def test_direct_final_file_partial_write_is_forward_recoverable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    payload = b"immutable reviewed bytes\n"
    digest = hashlib.sha256(payload).hexdigest()
    real_write = authority.os.write
    interrupted = False

    def partial_then_raise(descriptor: int, retained: bytes) -> int:
        nonlocal interrupted
        if not interrupted:
            interrupted = True
            assert real_write(descriptor, retained[:5]) == 5
            raise RuntimeError("injected direct-final partial write")
        return real_write(descriptor, retained)

    monkeypatch.setattr(authority.os, "write", partial_then_raise)
    with pytest.raises(ValueError, match="retained authenticated recovery entry"):
        authority.ensure_direct_final_file(payload, target, digest, 0o444, "fixture artifact")

    assert interrupted
    private = [
        tmp_path / name
        for name in authority.private_publication_names(target.name, "fixture artifact")
    ]
    assert not target.exists()
    assert stat.S_IMODE(private[0].stat().st_mode) == 0o000
    assert private[0].stat().st_size == 5

    monkeypatch.setattr(authority.os, "write", real_write)
    authority.ensure_direct_final_file(payload, target, digest, 0o444, "fixture artifact")
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    assert [path for path in private if path.exists()] == [private[0]]


def test_direct_final_repeated_private_write_crashes_remain_recoverable_and_bounded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    payload = b"immutable reviewed bytes\n"
    digest = hashlib.sha256(payload).hexdigest()
    real_write = authority.os.write
    private = [
        tmp_path / name
        for name in authority.private_publication_names(target.name, "fixture artifact")
    ]

    for retained_size in (3, 5, 7):
        interrupted = False

        def partial_then_raise(descriptor: int, retained: bytes) -> int:
            nonlocal interrupted
            if not interrupted:
                interrupted = True
                written = min(retained_size, len(retained))
                assert real_write(descriptor, retained[:written]) == written
                raise RuntimeError("injected repeated private write interruption")
            return real_write(descriptor, retained)

        monkeypatch.setattr(authority.os, "write", partial_then_raise)
        with pytest.raises(ValueError, match="retained authenticated recovery entry"):
            authority.ensure_direct_final_file(
                payload, target, digest, 0o444, "fixture artifact"
            )
        assert interrupted
        assert not target.exists()
        assert 1 <= len([path for path in private if path.exists()]) <= 2
        assert all(
            stat.S_IMODE(path.stat().st_mode) == 0o000
            for path in private
            if path.exists()
        )

    monkeypatch.setattr(authority.os, "write", real_write)
    authority.ensure_direct_final_file(payload, target, digest, 0o444, "fixture artifact")
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    assert target.stat().st_nlink == 1
    assert len([path for path in private if path.exists()]) == 1


def test_direct_final_rejects_hardlinked_incomplete_private_slot(
    tmp_path: Path,
) -> None:
    target = tmp_path / "artifact"
    alias = tmp_path / "injected-alias"
    payload = b"immutable reviewed bytes\n"
    private = tmp_path / authority.private_publication_names(
        target.name, "fixture artifact"
    )[0]
    write_bytes(private, b"partial", 0o000)
    os.link(private, alias)

    with pytest.raises(ValueError, match="external hardlink"):
        authority.ensure_direct_final_file(
            payload,
            target,
            hashlib.sha256(payload).hexdigest(),
            0o444,
            "fixture artifact",
        )

    assert not target.exists()
    assert private.stat().st_nlink == 2
    assert alias.stat().st_nlink == 2


def test_direct_final_recycles_bounded_owner_only_private_slots(tmp_path: Path) -> None:
    target = tmp_path / "artifact"
    payload = b"immutable reviewed bytes\n"
    private = [
        tmp_path / name
        for name in authority.private_publication_names(target.name, "fixture artifact")
    ]
    for path in private:
        write_bytes(path, b"partial", 0o600)

    authority.ensure_direct_final_file(
        payload,
        target,
        hashlib.sha256(payload).hexdigest(),
        0o444,
        "fixture artifact",
    )

    assert target.read_bytes() == payload
    assert target.stat().st_nlink == 1
    assert private[0].exists() is False
    assert private[1].read_bytes() == b"partial"


@pytest.mark.parametrize("initial", ["predecessor", "owner-only-partial"])
def test_direct_catalog_completion_never_renames(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, initial: str
) -> None:
    target = tmp_path / "README.md"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(
        target,
        old_payload if initial == "predecessor" else b"partial",
        0o644 if initial == "predecessor" else 0o600,
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("direct catalog completion must not rename")

    monkeypatch.setattr(authority, "renameat2", forbidden)
    monkeypatch.setattr(authority, "renameat2_between", forbidden)
    monkeypatch.setattr(authority.os, "rename", forbidden)
    authority.ensure_direct_catalog(
        target,
        new_payload,
        hashlib.sha256(old_payload).hexdigest(),
        hashlib.sha256(new_payload).hexdigest(),
        "fixture catalog",
    )

    assert target.read_bytes() == new_payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert target.stat().st_nlink == 1


def test_direct_catalog_hardlink_injection_never_mutates_public_inode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "README.md"
    alias = tmp_path / "injected-alias"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(target, old_payload, 0o644)
    real_link = authority.os.link
    real_write = authority.os.write
    real_fchmod = authority.os.fchmod
    published = False

    def inject_alias(source: str, destination: str, *args: object, **kwargs: object) -> None:
        nonlocal published
        real_link(source, destination, *args, **kwargs)
        if destination == target.name:
            parent = kwargs["dst_dir_fd"]
            real_link(
                destination,
                alias.name,
                src_dir_fd=parent,
                dst_dir_fd=parent,
                follow_symlinks=False,
            )
            published = True

    def reject_public_write(descriptor: int, retained: bytes) -> int:
        if published:
            raise AssertionError("public inode was written after publication")
        return real_write(descriptor, retained)

    def reject_public_fchmod(descriptor: int, mode: int) -> None:
        if published:
            raise AssertionError("public inode was chmodded after publication")
        real_fchmod(descriptor, mode)

    monkeypatch.setattr(authority.os, "link", inject_alias)
    monkeypatch.setattr(authority.os, "write", reject_public_write)
    monkeypatch.setattr(authority.os, "fchmod", reject_public_fchmod)
    with pytest.raises(ValueError, match="unsafe link profile"):
        authority.ensure_direct_catalog(
            target,
            new_payload,
            hashlib.sha256(old_payload).hexdigest(),
            hashlib.sha256(new_payload).hexdigest(),
            "fixture catalog",
        )

    assert published
    assert target.read_bytes() == new_payload
    assert alias.read_bytes() == new_payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert target.stat().st_nlink == 3


def test_direct_catalog_partial_write_is_forward_recoverable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "README.md"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(target, old_payload, 0o644)
    real_write = authority.os.write
    interrupted = False

    def partial_then_raise(descriptor: int, retained: bytes) -> int:
        nonlocal interrupted
        if not interrupted:
            interrupted = True
            assert real_write(descriptor, retained[:6]) == 6
            raise RuntimeError("injected direct-catalog partial write")
        return real_write(descriptor, retained)

    monkeypatch.setattr(authority.os, "write", partial_then_raise)
    with pytest.raises(ValueError, match="retained authenticated recovery entry"):
        authority.ensure_direct_catalog(
            target,
            new_payload,
            hashlib.sha256(old_payload).hexdigest(),
            hashlib.sha256(new_payload).hexdigest(),
            "fixture catalog",
        )

    assert interrupted
    private = [
        tmp_path / name
        for name in authority.private_publication_names(target.name, "fixture catalog")
    ]
    assert target.read_bytes() == old_payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert stat.S_IMODE(private[0].stat().st_mode) == 0o000
    assert private[0].stat().st_size == 6

    monkeypatch.setattr(authority.os, "write", real_write)
    authority.ensure_direct_catalog(
        target,
        new_payload,
        hashlib.sha256(old_payload).hexdigest(),
        hashlib.sha256(new_payload).hexdigest(),
        "fixture catalog",
    )
    assert target.read_bytes() == new_payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert [path for path in private if path.exists()] == [private[0]]


def test_direct_catalog_recovers_after_predecessor_unlink_crash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "README.md"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(target, old_payload, 0o644)
    real_unlink = authority.unlink_bound_name_lustre
    interrupted = False

    def unlink_then_raise(
        parent: int, name: str, expected: os.stat_result, label: str
    ) -> None:
        nonlocal interrupted
        real_unlink(parent, name, expected, label)
        if not interrupted and "reviewed predecessor" in label:
            interrupted = True
            raise RuntimeError("injected crash after predecessor unlink")

    monkeypatch.setattr(authority, "unlink_bound_name_lustre", unlink_then_raise)
    with pytest.raises(RuntimeError, match="injected crash after predecessor unlink"):
        authority.ensure_direct_catalog(
            target,
            new_payload,
            hashlib.sha256(old_payload).hexdigest(),
            hashlib.sha256(new_payload).hexdigest(),
            "fixture catalog",
        )

    private = [
        tmp_path / name
        for name in authority.private_publication_names(target.name, "fixture catalog")
    ]
    assert interrupted
    assert not target.exists()
    assert [path.read_bytes() for path in private if path.exists()] == [new_payload]
    assert all(
        stat.S_IMODE(path.stat().st_mode) == 0o644
        for path in private
        if path.exists()
    )

    monkeypatch.setattr(authority, "unlink_bound_name_lustre", real_unlink)
    authority.ensure_direct_catalog(
        target,
        new_payload,
        hashlib.sha256(old_payload).hexdigest(),
        hashlib.sha256(new_payload).hexdigest(),
        "fixture catalog",
    )
    assert target.read_bytes() == new_payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert target.stat().st_nlink == 1
    assert not any(path.exists() for path in private)


def test_hardlink_install_crash_window_is_forward_recoverable(tmp_path: Path) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    write_bytes(temporary, payload, 0o444)
    os.link(temporary, target)
    assert target.stat().st_nlink == 2
    authority.ensure_installed(
        payload,
        target,
        hashlib.sha256(payload).hexdigest(),
        0o444,
        transaction_id,
        "fixture artifact",
    )
    assert target.read_bytes() == payload
    assert target.stat().st_nlink == 1
    assert not temporary.exists()


@pytest.mark.parametrize("remnant_mode", [0o000, 0o444])
def test_partial_artifact_temporary_is_retired_and_recreated(
    tmp_path: Path, remnant_mode: int
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    partial = b"partial"
    write_bytes(temporary, partial, remnant_mode)
    retained_before = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))

    authority.ensure_installed(
        payload,
        target,
        hashlib.sha256(payload).hexdigest(),
        0o444,
        transaction_id,
        "fixture artifact",
    )

    retained_after = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))
    assert len(retained_after - retained_before) == 2
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    assert not temporary.exists()


def test_partial_artifact_temporary_substitution_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped-partial"
    payload = b"immutable reviewed bytes\n"
    write_bytes(temporary, b"partial", 0o000)
    write_bytes(substitute, payload, 0o444)
    real_read_bound_file = authority.read_bound_file
    raced = False

    def substitute_before_temporary_read(
        parent: int, name: str, label: str, **kwargs
    ) -> tuple[bytes, str, os.stat_result]:
        nonlocal raced
        if not raced and name == temporary.name:
            raced = True
            os.rename(name, escaped.name, src_dir_fd=parent, dst_dir_fd=parent)
            os.rename(substitute.name, name, src_dir_fd=parent, dst_dir_fd=parent)
        return real_read_bound_file(parent, name, label, **kwargs)

    monkeypatch.setattr(authority, "read_bound_file", substitute_before_temporary_read)
    with pytest.raises(ValueError, match="inode identity changed before mutation"):
        authority.ensure_installed(
            payload,
            target,
            hashlib.sha256(payload).hexdigest(),
            0o444,
            transaction_id,
            "fixture artifact",
        )

    assert raced
    assert not target.exists()
    assert temporary.read_bytes() == payload
    assert escaped.stat().st_size == len(b"partial")


def test_untrusted_artifact_temporary_is_rejected_without_retirement(
    tmp_path: Path,
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    outside = tmp_path / "outside"
    payload = b"immutable reviewed bytes\n"
    write_bytes(outside, b"untrusted bytes\n", 0o444)
    temporary.symlink_to(outside)
    retained_before = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))

    with pytest.raises(ValueError, match="must be a regular file"):
        authority.ensure_installed(
            payload,
            target,
            hashlib.sha256(payload).hexdigest(),
            0o444,
            transaction_id,
            "fixture artifact",
        )

    assert temporary.is_symlink()
    assert outside.read_bytes() == b"untrusted bytes\n"
    assert set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic")) == retained_before


@pytest.mark.parametrize("drift", ["mode", "link", "content"])
def test_artifact_publication_rejects_in_place_recovery_file_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, drift: str
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    real_publish = authority.publish_bound_file_noreplace
    raced = False

    def drift_before_publication(
        parent: int, source: str, target_name: str,
        expected: os.stat_result, label: str,
        **kwargs: object,
    ) -> os.stat_result:
        nonlocal raced
        if not raced and source == temporary.name:
            raced = True
            drift_bound_file(parent, source, drift, b"X" * len(payload), 0o444)
        return real_publish(parent, source, target_name, expected, label, **kwargs)

    monkeypatch.setattr(authority, "publish_bound_file_noreplace", drift_before_publication)
    with pytest.raises(
        ValueError, match="file security/content profile changed|checksum differs"
    ):
        authority.ensure_installed(
            payload,
            target,
            hashlib.sha256(payload).hexdigest(),
            0o444,
            transaction_id,
            "fixture artifact",
        )

    assert raced
    assert not target.exists()
    assert temporary.exists()


@pytest.mark.parametrize("phase", ["after-rename", "after-fsync"])
def test_artifact_forged_public_bytes_are_retry_recoverable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, phase: str
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    recovery = tmp_path / f".{target.name}.single-link.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    forged = b"X" * len(payload)
    digest = hashlib.sha256(payload).hexdigest()
    real_renameat2 = authority.renameat2
    real_fsync = authority.fsync_descriptor
    mutated = False

    def mutate_after_publication_rename(
        parent: int, source: str, target_name: str, flags: int, label: str
    ) -> None:
        nonlocal mutated
        real_renameat2(parent, source, target_name, flags, label)
        if (
            not mutated
            and phase == "after-rename"
            and flags == authority.RENAME_NOREPLACE
            and target_name == target.name
        ):
            mutated = True
            drift_bound_file(parent, target_name, "content", forged, 0o444)

    def mutate_after_publication_fsync(descriptor: int) -> None:
        nonlocal mutated
        real_fsync(descriptor)
        if (
            not mutated
            and phase == "after-fsync"
            and target.exists()
            and target.read_bytes() == payload
        ):
            mutated = True
            drift_bound_file(descriptor, target.name, "content", forged, 0o444)

    monkeypatch.setattr(authority, "renameat2", mutate_after_publication_rename)
    monkeypatch.setattr(authority, "fsync_descriptor", mutate_after_publication_fsync)
    with pytest.raises(ValueError, match="profile changed|checksum differs"):
        authority.ensure_installed(
            payload, target, digest, 0o444, transaction_id, "fixture artifact"
        )

    assert mutated
    assert target.read_bytes() == forged
    assert recovery.read_bytes() == payload

    monkeypatch.setattr(authority, "renameat2", real_renameat2)
    monkeypatch.setattr(authority, "fsync_descriptor", real_fsync)
    authority.ensure_installed(
        payload, target, digest, 0o444, transaction_id, "fixture artifact"
    )

    assert target.read_bytes() == payload
    assert target.stat().st_nlink == 1
    assert not recovery.exists()


@pytest.mark.parametrize("remnant_mode", [0o000, 0o644])
def test_partial_catalog_temporary_is_retired_and_recreated(
    tmp_path: Path, remnant_mode: int
) -> None:
    target = tmp_path / "README.md"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(target, old_payload, 0o644)
    write_bytes(temporary, b"partial", remnant_mode)
    retained_before = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))

    authority.ensure_catalog(
        target,
        new_payload,
        hashlib.sha256(old_payload).hexdigest(),
        hashlib.sha256(new_payload).hexdigest(),
        transaction_id,
        "fixture catalog",
    )

    retained_after = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))
    assert len(retained_after - retained_before) == 3
    assert target.read_bytes() == new_payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert not temporary.exists()


def test_partial_catalog_temporary_substitution_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "README.md"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped-partial"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(target, old_payload, 0o644)
    write_bytes(temporary, b"partial", 0o000)
    write_bytes(substitute, new_payload, 0o644)
    real_read_bound_file = authority.read_bound_file
    raced = False

    def substitute_before_temporary_read(
        parent: int, name: str, label: str, **kwargs
    ) -> tuple[bytes, str, os.stat_result]:
        nonlocal raced
        if not raced and name == temporary.name:
            raced = True
            os.rename(name, escaped.name, src_dir_fd=parent, dst_dir_fd=parent)
            os.rename(substitute.name, name, src_dir_fd=parent, dst_dir_fd=parent)
        return real_read_bound_file(parent, name, label, **kwargs)

    monkeypatch.setattr(authority, "read_bound_file", substitute_before_temporary_read)
    with pytest.raises(ValueError, match="inode identity changed before mutation"):
        authority.ensure_catalog(
            target,
            new_payload,
            hashlib.sha256(old_payload).hexdigest(),
            hashlib.sha256(new_payload).hexdigest(),
            transaction_id,
            "fixture catalog",
        )

    assert raced
    assert target.read_bytes() == old_payload
    assert temporary.read_bytes() == new_payload
    assert escaped.stat().st_size == len(b"partial")


@pytest.mark.parametrize("drift", ["mode", "link", "content"])
def test_catalog_publication_rejects_in_place_recovery_file_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, drift: str
) -> None:
    target = tmp_path / "README.md"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    write_bytes(target, old_payload, 0o644)
    real_exchange = authority.exchange_bound_entries
    raced = False

    def drift_before_exchange(
        parent: int, source: str, target_name: str,
        source_expected: os.stat_result, target_expected: os.stat_result,
        label: str,
        **kwargs: object,
    ) -> tuple[os.stat_result, os.stat_result]:
        nonlocal raced
        if not raced and source == temporary.name:
            raced = True
            drift_bound_file(parent, source, drift, b"X" * len(new_payload), 0o644)
        return real_exchange(
            parent, source, target_name, source_expected, target_expected, label, **kwargs
        )

    monkeypatch.setattr(authority, "exchange_bound_entries", drift_before_exchange)
    with pytest.raises(
        ValueError, match="file security/content profile changed|checksum differs"
    ):
        authority.ensure_catalog(
            target,
            new_payload,
            hashlib.sha256(old_payload).hexdigest(),
            hashlib.sha256(new_payload).hexdigest(),
            transaction_id,
            "fixture catalog",
        )

    assert raced
    assert target.read_bytes() == old_payload
    assert temporary.exists()


@pytest.mark.parametrize("phase", ["after-exchange", "after-fsync"])
def test_catalog_forged_public_bytes_are_rolled_back_and_retryable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, phase: str
) -> None:
    target = tmp_path / "README.md"
    transaction_id = "fixture"
    old_payload = b"reviewed predecessor\n"
    new_payload = b"reviewed F118 catalog\n"
    forged = b"X" * len(new_payload)
    old_digest = hashlib.sha256(old_payload).hexdigest()
    new_digest = hashlib.sha256(new_payload).hexdigest()
    write_bytes(target, old_payload, 0o644)
    real_renameat2 = authority.renameat2
    real_fsync = authority.fsync_descriptor
    mutated = False

    def mutate_after_catalog_exchange(
        parent: int, source: str, target_name: str, flags: int, label: str
    ) -> None:
        nonlocal mutated
        real_renameat2(parent, source, target_name, flags, label)
        if (
            not mutated
            and phase == "after-exchange"
            and flags == authority.RENAME_EXCHANGE
            and "rollback" not in label
        ):
            mutated = True
            drift_bound_file(parent, target_name, "content", forged, 0o644)

    def mutate_after_catalog_fsync(descriptor: int) -> None:
        nonlocal mutated
        real_fsync(descriptor)
        if (
            not mutated
            and phase == "after-fsync"
            and target.read_bytes() == new_payload
        ):
            mutated = True
            drift_bound_file(descriptor, target.name, "content", forged, 0o644)

    monkeypatch.setattr(authority, "renameat2", mutate_after_catalog_exchange)
    monkeypatch.setattr(authority, "fsync_descriptor", mutate_after_catalog_fsync)
    with pytest.raises(ValueError, match="profile changed|checksum differs"):
        authority.ensure_catalog(
            target,
            new_payload,
            old_digest,
            new_digest,
            transaction_id,
            "fixture catalog",
        )

    assert mutated
    assert target.read_bytes() == old_payload

    monkeypatch.setattr(authority, "renameat2", real_renameat2)
    monkeypatch.setattr(authority, "fsync_descriptor", real_fsync)
    authority.ensure_catalog(
        target,
        new_payload,
        old_digest,
        new_digest,
        transaction_id,
        "fixture catalog",
    )

    assert target.read_bytes() == new_payload
    assert target.stat().st_nlink == 1


def test_completed_single_link_recovery_copy_is_forward_recoverable(tmp_path: Path) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    recovery = tmp_path / f".{target.name}.single-link.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    write_bytes(temporary, payload, 0o444)
    os.link(temporary, target)
    write_bytes(recovery, payload, 0o444)
    assert target.stat().st_nlink == temporary.stat().st_nlink == 2
    assert recovery.stat().st_nlink == 1

    authority.ensure_installed(
        payload,
        target,
        hashlib.sha256(payload).hexdigest(),
        0o444,
        transaction_id,
        "fixture artifact",
    )

    assert target.read_bytes() == payload
    assert target.stat().st_nlink == 1
    assert not temporary.exists()
    assert not recovery.exists()


def test_partial_single_link_recovery_copy_is_retired_and_recreated(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    recovery = tmp_path / f".{target.name}.single-link.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    partial = b"partial"
    write_bytes(temporary, payload, 0o444)
    os.link(temporary, target)
    write_bytes(recovery, partial, 0o000)
    retained_before = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))

    def forbidden_unlink(*args, **kwargs) -> None:
        del args, kwargs
        raise AssertionError("incomplete recovery retirement must never call unlink")

    monkeypatch.setattr(authority.os, "unlink", forbidden_unlink)
    authority.ensure_installed(
        payload,
        target,
        hashlib.sha256(payload).hexdigest(),
        0o444,
        transaction_id,
        "fixture artifact",
    )

    retained_after = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))
    retained_partial = [
        path for path in retained_after - retained_before
        if stat.S_IMODE(path.stat().st_mode) == 0o000 and path.stat().st_size == len(partial)
    ]
    assert len(retained_partial) == 1
    assert target.read_bytes() == payload
    assert target.stat().st_nlink == 1
    assert not temporary.exists()
    assert not recovery.exists()


def test_partial_single_link_recovery_substitution_is_not_accepted_or_retired(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    recovery = tmp_path / f".{target.name}.single-link.{transaction_id}.tmp"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped-partial"
    payload = b"immutable reviewed bytes\n"
    partial = b"partial"
    write_bytes(temporary, payload, 0o444)
    os.link(temporary, target)
    write_bytes(recovery, partial, 0o000)
    write_bytes(substitute, payload, 0o444)
    real_read_bound_file = authority.read_bound_file
    raced = False

    def substitute_before_recovery_read(
        parent: int, name: str, label: str, **kwargs
    ) -> tuple[bytes, str, os.stat_result]:
        nonlocal raced
        if not raced and name == recovery.name:
            raced = True
            os.rename(name, escaped.name, src_dir_fd=parent, dst_dir_fd=parent)
            os.rename(substitute.name, name, src_dir_fd=parent, dst_dir_fd=parent)
        return real_read_bound_file(parent, name, label, **kwargs)

    monkeypatch.setattr(authority, "read_bound_file", substitute_before_recovery_read)
    with pytest.raises(ValueError, match="inode identity changed before mutation"):
        authority.ensure_installed(
            payload,
            target,
            hashlib.sha256(payload).hexdigest(),
            0o444,
            transaction_id,
            "fixture artifact",
        )

    assert raced
    assert target.stat().st_nlink == temporary.stat().st_nlink == 2
    assert recovery.read_bytes() == payload
    assert stat.S_IMODE(recovery.stat().st_mode) == 0o444
    assert escaped.stat().st_size == len(partial)
    assert stat.S_IMODE(escaped.stat().st_mode) == 0o000


def test_hardlink_single_link_recovery_interruption_is_forward_recoverable(
    tmp_path: Path,
) -> None:
    target = tmp_path / "artifact"
    transaction_id = "fixture"
    temporary = tmp_path / f".{target.name}.{transaction_id}.tmp"
    recovery = tmp_path / f".{target.name}.single-link.{transaction_id}.tmp"
    payload = b"immutable reviewed bytes\n"
    write_bytes(target, payload, 0o444)
    write_bytes(temporary, payload, 0o444)
    os.link(temporary, recovery)
    assert target.stat().st_nlink == 1
    assert temporary.stat().st_nlink == recovery.stat().st_nlink == 2

    authority.ensure_installed(
        payload,
        target,
        hashlib.sha256(payload).hexdigest(),
        0o444,
        transaction_id,
        "fixture artifact",
    )

    assert target.read_bytes() == payload
    assert target.stat().st_nlink == 1
    assert not temporary.exists()
    assert not recovery.exists()


def test_verify_rejects_wrong_external_audit_binding(campaign: Campaign) -> None:
    assert campaign.promote().returncode == 0
    command = campaign.verify_command()
    command[-1] = "0" * 64
    assert_failed(campaign.run(command), "publication audit checksum differs")


def test_git_execution_is_descriptor_bound_and_caller_independent(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    for name, value in {
        "PATH": "/tmp/hostile-path",
        "LD_PRELOAD": "/tmp/hostile.so",
        "PYTHONPATH": "/tmp/hostile-python",
        "GIT_CONFIG_PARAMETERS": "hostile",
    }.items():
        monkeypatch.setenv(name, value)
    calls: list[tuple[list[str], dict[str, object]]] = []

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout=b"")

    monkeypatch.setattr(authority.subprocess, "run", fake_run)
    authority.git_run(campaign.repository, ["diff", "--quiet", "--", "final.txt"])

    assert len(calls) == 1
    command, kwargs = calls[0]
    descriptor = int(str(kwargs["executable"]).removeprefix("/proc/self/fd/"))
    assert descriptor in kwargs["pass_fds"]
    assert command[:4] == [
        str(authority.GIT), "--no-replace-objects", "-C", str(campaign.repository)
    ]
    assert "--no-ext-diff" in command
    assert "--no-textconv" in command
    assert f"core.hooksPath={os.devnull}" in command
    assert kwargs["stdin"] == subprocess.DEVNULL
    assert kwargs["env"] == authority.hardened_git_environment()
    assert kwargs["env"]["PATH"] == authority.TRUSTED_SYSTEM_PATH
    assert "LD_PRELOAD" not in kwargs["env"]
    assert "PYTHONPATH" not in kwargs["env"]


def test_publisher_reexec_environment_strips_private_loader_and_interpreter_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for name, value in {
        authority.SELF_DESCRIPTOR_ENV: "99",
        authority.SELF_SOURCE_ENV: "/tmp/hostile-source",
        authority.REPOSITORY_ROOT_ENV: "/tmp/hostile-repository",
        "PATH": "/tmp/hostile-path",
        "LD_PRELOAD": "/tmp/hostile.so",
        "PYTHONPATH": "/tmp/hostile-python",
        "GIT_OBJECT_DIRECTORY": "/tmp/hostile-objects",
    }.items():
        monkeypatch.setenv(name, value)

    environment = authority.reexec_environment(7, 8, PUBLISHER, REPOSITORY)

    assert environment[authority.SELF_DESCRIPTOR_ENV] == "7"
    assert environment[authority.PYTHON_DESCRIPTOR_ENV] == "8"
    assert environment[authority.SELF_SOURCE_ENV] == str(PUBLISHER)
    assert environment[authority.REPOSITORY_ROOT_ENV] == str(REPOSITORY)
    assert environment["PATH"] == authority.TRUSTED_SYSTEM_PATH
    assert environment["PYTHONDONTWRITEBYTECODE"] == "1"
    assert environment["HOME"] == "/nonexistent"
    assert environment["XDG_CONFIG_HOME"] == "/nonexistent"
    assert set(environment) == {
        authority.SELF_DESCRIPTOR_ENV,
        authority.PYTHON_DESCRIPTOR_ENV,
        authority.SELF_SOURCE_ENV,
        authority.REPOSITORY_ROOT_ENV,
        "HOME",
        "LC_ALL",
        "PATH",
        "PYTHONDONTWRITEBYTECODE",
        "XDG_CONFIG_HOME",
    }
    assert "LD_PRELOAD" not in environment
    assert "PYTHONPATH" not in environment
    assert "GIT_OBJECT_DIRECTORY" not in environment


def test_publisher_rejects_orphaned_private_reexec_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv(authority.SELF_DESCRIPTOR_ENV, raising=False)
    monkeypatch.setenv(authority.SELF_SOURCE_ENV, str(PUBLISHER))
    monkeypatch.setenv(authority.REPOSITORY_ROOT_ENV, str(REPOSITORY))

    with pytest.raises(ValueError, match="forbidden without an authenticated descriptor"):
        authority.authenticate_self([
            "--expected-publisher-sha256",
            sha256(PUBLISHER),
        ])


def test_authenticated_publisher_descriptor_controls_private_metadata(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    descriptor = os.open(campaign.publisher, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(), os.O_RDONLY)
    try:
        monkeypatch.setattr(authority, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(authority.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(authority.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(authority.SELF_SOURCE_ENV, str(campaign.publisher))
        monkeypatch.setenv(authority.REPOSITORY_ROOT_ENV, str(campaign.repository))
        monkeypatch.setattr(
            authority.sys,
            "flags",
            type("Flags", (), {"isolated": 1})(),
        )

        source, repository, digest = authority.authenticate_self([
            "--expected-publisher-sha256",
            sha256(campaign.publisher),
        ])
    finally:
        os.close(descriptor)
        os.close(python_descriptor)

    assert source == campaign.publisher
    assert repository == campaign.repository
    assert digest == sha256(campaign.publisher)


def test_authenticated_publisher_rejects_nonisolated_python(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    descriptor = os.open(campaign.publisher, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(), os.O_RDONLY)
    try:
        monkeypatch.setattr(authority, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(authority.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(authority.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(authority.SELF_SOURCE_ENV, str(campaign.publisher))
        monkeypatch.setenv(authority.REPOSITORY_ROOT_ENV, str(campaign.repository))
        monkeypatch.setattr(
            authority.sys,
            "flags",
            type("Flags", (), {"isolated": 0})(),
        )
        with pytest.raises(ValueError, match="interpreter is not isolated"):
            authority.authenticate_self([
                "--expected-publisher-sha256",
                sha256(campaign.publisher),
            ])
    finally:
        os.close(descriptor)
        os.close(python_descriptor)


def test_publisher_rejects_different_root_executable_as_python_descriptor(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    descriptor = os.open(authority.GIT, os.O_RDONLY)
    try:
        monkeypatch.setenv(authority.PYTHON_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setattr(
            authority.sys,
            "flags",
            type("Flags", (), {"isolated": 1})(),
        )
        with pytest.raises(ValueError, match="is not this interpreter"):
            authority.require_authenticated_python_descriptor()
    finally:
        os.close(descriptor)


def test_authenticated_publisher_rejects_named_source_substitution(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    descriptor = os.open(campaign.publisher, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(), os.O_RDONLY)
    retained = campaign.publisher.with_name("retained-publisher.py")
    campaign.publisher.rename(retained)
    write_bytes(campaign.publisher, b"#!/usr/bin/false\n", 0o755)
    try:
        monkeypatch.setattr(authority, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(authority.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(authority.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(authority.SELF_SOURCE_ENV, str(campaign.publisher))
        monkeypatch.setenv(authority.REPOSITORY_ROOT_ENV, str(campaign.repository))
        monkeypatch.setattr(
            authority.sys,
            "flags",
            type("Flags", (), {"isolated": 1})(),
        )
        with pytest.raises(ValueError, match="checksum differs after reexecution"):
            authority.authenticate_self([
                "--expected-publisher-sha256",
                hashlib.sha256(retained.read_bytes()).hexdigest(),
            ])
    finally:
        os.close(descriptor)
        os.close(python_descriptor)


def test_source_authority_lock_replacement_blocks_next_write(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    lock = root / f".mks24_stage_i_{authority.EXECUTION_EPOCH_SLUG}.lock"
    write_bytes(lock, b"", 0o644)
    target = root / "must-not-be-created.json"

    with pytest.raises(ValueError, match="Stage I lock path changed while mutation is active"):
        with authority.stage_i_lock({"lock": lock}):
            lock.rename(root / "retired-lock")
            write_bytes(lock, b"", 0o644)
            authority.write_exclusive(target, b"forbidden\n", 0o600)

    assert not target.exists()


def test_source_authority_exchange_race_does_not_clobber_substitute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped"
    write_bytes(source, b"new\n", 0o600)
    write_bytes(target, b"old\n", 0o600)
    write_bytes(substitute, b"substitute\n", 0o600)
    real_renameat2 = authority.renameat2

    with authority.bound_directory(tmp_path, "exchange parent") as (parent, _):
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        target_profile = os.stat(target.name, dir_fd=parent, follow_symlinks=False)

        def race_after_exchange(
            selected_parent: int, source_name: str, target_name: str,
            flags: int, label: str,
        ) -> None:
            real_renameat2(selected_parent, source_name, target_name, flags, label)
            if flags == authority.RENAME_EXCHANGE and "rollback" not in label:
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

        monkeypatch.setattr(authority, "renameat2", race_after_exchange)
        with pytest.raises(ValueError, match="inode identity changed"):
            authority.exchange_bound_entries(
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


def test_source_authority_requires_noncryptographic_review_limitation(
    campaign: Campaign,
) -> None:
    evidence = json.loads(campaign.evidence_candidate.read_text())
    provenance = json.loads(campaign.provenance_candidate.read_text())
    plasma = json.loads(campaign.plasma_candidate.read_text())
    provenance["limitations"].remove(
        authority.INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION
    )

    with pytest.raises(ValueError, match="non-cryptographic reviewer identity limitation"):
        authority.parse_reviews(
            provenance,
            plasma,
            evidence_candidate_path=campaign.evidence_candidate,
            evidence_sha256=sha256(campaign.evidence_candidate),
            evidence_path=campaign.root / authority.F118_PATHS["evidence"],
            final=evidence["implementation"]["current_source_bundle"],
            generated_utc=evidence["generated_utc"],
        )


def test_source_authority_declared_independence_disclaims_cryptographic_identity() -> None:
    assurance = authority.declared_process_independence_assurance(
        {"provenance-security": "reviewer-a", "plasma-scientific": "reviewer-b"},
        "F118 review",
    )
    assert assurance["cryptographic_identity_verified"] is False
    assert (
        assurance["non_cryptographic_limitation"]
        == authority.INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION
    )


def test_canonical_source_authority_requires_canonical_repository(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="canonical use requires"):
        authority.require_canonical_repository(tmp_path, True)
    authority.require_canonical_repository(authority.CANONICAL_REPOSITORY_ROOT, True)
    authority.require_canonical_repository(tmp_path, False)


def test_canonical_public_namespace_accepts_only_exact_trusted_project_boundary() -> None:
    trusted = SimpleNamespace(
        st_mode=stat.S_IFDIR | 0o2770,
        st_uid=0,
        st_gid=authority.CANONICAL_TRUSTED_PROJECT_GID,
    )
    authority.require_canonical_public_namespace_profile(
        authority.CANONICAL_TRUSTED_PROJECT_BOUNDARY, trusted
    )

    world_writable = SimpleNamespace(
        st_mode=stat.S_IFDIR | 0o2777,
        st_uid=0,
        st_gid=authority.CANONICAL_TRUSTED_PROJECT_GID,
    )
    with pytest.raises(ValueError, match="untrusted"):
        authority.require_canonical_public_namespace_profile(
            authority.CANONICAL_TRUSTED_PROJECT_BOUNDARY, world_writable
        )

    group_writable_descendant = SimpleNamespace(
        st_mode=stat.S_IFDIR | 0o2770,
        st_uid=authority.CANONICAL_OWNER_UID,
        st_gid=authority.CANONICAL_TRUSTED_PROJECT_GID,
    )
    with pytest.raises(ValueError, match="untrusted|group-writable"):
        authority.require_canonical_public_namespace_profile(
            authority.CANONICAL_TRUSTED_PROJECT_BOUNDARY / "dfielding",
            group_writable_descendant,
        )


def test_canonical_boundary_2770_to_2777_blocks_inner_create_before_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    boundary, root = configure_fake_canonical_namespace(tmp_path, monkeypatch)
    target = root / "forbidden"

    with authority.bound_canonical_public_namespace(root, True):
        with authority.bound_directory(root, "canonical mutation parent") as (parent, _):
            boundary.chmod(0o2777)
            try:
                with pytest.raises(ValueError, match="profile is untrusted"):
                    authority.write_bound_exclusive(
                        parent, target.name, b"forbidden\n", 0o600, "canonical create"
                    )
            finally:
                boundary.chmod(0o2770)

    assert not target.exists()


def test_canonical_boundary_2770_to_2777_during_create_is_durable_and_stops_writes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    boundary, root = configure_fake_canonical_namespace(tmp_path, monkeypatch)
    target = root / "ambiguous"
    real_open = authority.os.open
    real_fsync = authority.os.fsync
    created = False
    fsynced: list[int] = []

    def create_then_drift(path: object, flags: int, *args: object, **kwargs: object) -> int:
        nonlocal created
        descriptor = real_open(path, flags, *args, **kwargs)
        if os.fspath(path) == target.name and flags & os.O_CREAT:
            created = True
            boundary.chmod(0o2777)
        return descriptor

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority.os, "open", create_then_drift)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_canonical_public_namespace(root, True):
        with authority.bound_directory(root, "canonical mutation parent") as (parent, _):
            try:
                with pytest.raises(ValueError, match="durably ambiguous after authority loss"):
                    authority.write_bound_exclusive(
                        parent, target.name, b"must not be written\n", 0o600,
                        "canonical create",
                    )
            finally:
                boundary.chmod(0o2770)
            assert parent in fsynced

    assert created
    assert stat.S_IMODE(target.stat().st_mode) == 0o000
    target.chmod(0o600)
    assert target.read_bytes() == b""


def test_write_bound_exclusive_retains_durable_mode_zero_partial_for_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "partial"
    payload = b"reviewed payload\n"
    real_write = authority.os.write
    real_fsync = authority.os.fsync
    failed = False
    fsynced: list[int] = []

    def partial_then_raise(descriptor: int, retained: bytes) -> int:
        nonlocal failed
        if not failed:
            failed = True
            assert real_write(descriptor, retained[:4]) == 4
            raise RuntimeError("injected partial write")
        return real_write(descriptor, retained)

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority.os, "write", partial_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(tmp_path, "create parent") as (parent, _):
        with pytest.raises(ValueError, match="retained authenticated recovery entry"):
            authority.write_bound_exclusive(
                parent, target.name, payload, 0o600, "partial create"
            )
        assert parent in fsynced

    assert failed
    assert stat.S_IMODE(target.stat().st_mode) == 0o000
    target.chmod(0o600)
    assert target.read_bytes() == payload[:4]
    target.chmod(0o000)

    monkeypatch.setattr(authority.os, "write", real_write)
    with authority.bound_directory(tmp_path, "create parent") as (parent, _):
        observed = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        authority.retire_bound_recovery_file(
            parent, target.name, observed, "partial create recovery"
        )
        authority.write_bound_exclusive(
            parent, target.name, payload, 0o600, "retried create"
        )

    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o600


def test_write_bound_exclusive_post_create_exception_fsyncs_mode_zero_recovery(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "post-create"
    real_open = authority.os.open
    real_fsync = authority.os.fsync
    fsynced: list[int] = []
    created = False

    def create_then_raise(path: object, flags: int, *args: object, **kwargs: object) -> int:
        nonlocal created
        descriptor = real_open(path, flags, *args, **kwargs)
        if os.fspath(path) == target.name and flags & os.O_CREAT:
            created = True
            os.close(descriptor)
            raise RuntimeError("injected after create syscall")
        return descriptor

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority.os, "open", create_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(tmp_path, "create parent") as (parent, _):
        with pytest.raises(ValueError, match="retained authenticated recovery entry"):
            authority.write_bound_exclusive(
                parent, target.name, b"must not be written\n", 0o600, "post-create"
            )
        assert parent in fsynced

    assert created
    assert target.stat().st_size == 0
    assert stat.S_IMODE(target.stat().st_mode) == 0o000


@pytest.mark.parametrize(
    ("name", "mode"),
    [("source-authority-transactions", 0o755), ("transaction.staging", 0o700)],
)
def test_mkdir_post_syscall_exception_is_durable_and_recoverable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str, mode: int
) -> None:
    real_mkdir = authority.os.mkdir
    real_fsync = authority.os.fsync
    created = False
    fsynced: list[int] = []

    def mkdir_then_raise(*args, **kwargs) -> None:
        nonlocal created
        real_mkdir(*args, **kwargs)
        created = True
        raise RuntimeError("injected after mkdir syscall")

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority.os, "mkdir", mkdir_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(tmp_path, "mkdir parent") as (parent, _):
        with pytest.raises(
            ValueError, match="retained authenticated recovery directory"
        ):
            authority.mkdir_bound_exclusive(parent, name, mode, "transaction mkdir")
        assert parent in fsynced
        with authority.bound_child_directory(
            parent, name, "retained transaction mkdir", mode=mode
        ):
            pass

    assert created
    assert stat.S_IMODE((tmp_path / name).stat().st_mode) == mode


def test_transaction_directories_retries_insertion_and_never_returns_stale_single(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transactions = tmp_path / "transactions"
    first = transactions / "first"
    second = transactions / "second"
    first.mkdir(parents=True, mode=0o700)
    real_non_forensic = authority.non_forensic_entry_names
    inserted = False

    def insert_second_after_first_scan(parent: int, label: str) -> list[str]:
        nonlocal inserted
        names = real_non_forensic(parent, label)
        if not inserted:
            inserted = True
            os.mkdir(second.name, mode=0o700, dir_fd=parent)
        return names

    monkeypatch.setattr(
        authority, "non_forensic_entry_names", insert_second_after_first_scan
    )
    retained = authority.transaction_directories({"transactions": transactions})

    assert inserted
    assert {path.name for path in retained} == {first.name, second.name}


def test_transaction_directories_retries_same_name_inode_replacement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transactions = tmp_path / "transactions"
    transaction = transactions / "transaction"
    replacement = tmp_path / "replacement"
    escaped = tmp_path / "escaped"
    transaction.mkdir(parents=True, mode=0o700)
    replacement.mkdir(mode=0o700)
    write_bytes(transaction / "old", b"old\n", 0o600)
    write_bytes(replacement / "new", b"new\n", 0o600)
    real_non_forensic = authority.non_forensic_entry_names
    replaced = False

    def replace_after_first_scan(parent: int, label: str) -> list[str]:
        nonlocal replaced
        names = real_non_forensic(parent, label)
        if not replaced:
            replaced = True
            transaction.rename(escaped)
            replacement.rename(transaction)
        return names

    monkeypatch.setattr(authority, "non_forensic_entry_names", replace_after_first_scan)
    retained = authority.transaction_directories({"transactions": transactions})

    assert replaced
    assert retained == [transaction]
    assert (retained[0] / "new").read_bytes() == b"new\n"
    assert (escaped / "old").read_bytes() == b"old\n"


def test_rename_noreplace_post_syscall_exception_fsyncs_namespace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    write_bytes(source, b"source\n", 0o600)
    real_renameat2 = authority.renameat2
    real_fsync = authority.os.fsync
    fsynced: list[int] = []

    def rename_then_raise(*args, **kwargs) -> None:
        real_renameat2(*args, **kwargs)
        raise RuntimeError("injected after namespace syscall")

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority, "renameat2", rename_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(tmp_path, "rename parent") as (parent, _):
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="completed ambiguously.*directory fsync"):
            authority.rename_bound_noreplace(
                parent, source.name, target.name, source_profile, "rename"
            )
        assert parent in fsynced

    assert not source.exists()
    assert target.read_bytes() == b"source\n"


def test_exchange_post_syscall_exception_fsyncs_namespace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    write_bytes(source, b"source\n", 0o600)
    write_bytes(target, b"target\n", 0o600)
    real_renameat2 = authority.renameat2
    real_fsync = authority.os.fsync
    fsynced: list[int] = []

    def exchange_then_raise(*args, **kwargs) -> None:
        real_renameat2(*args, **kwargs)
        raise RuntimeError("injected after namespace syscall")

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority, "renameat2", exchange_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(tmp_path, "exchange parent") as (parent, _):
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        target_profile = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="completed ambiguously.*directory fsync"):
            authority.exchange_bound_entries(
                parent,
                source.name,
                target.name,
                source_profile,
                target_profile,
                "exchange",
            )
        assert parent in fsynced

    assert source.read_bytes() == b"target\n"
    assert target.read_bytes() == b"source\n"


def test_unlink_retirement_post_syscall_exception_fsyncs_both_namespaces(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    container = tmp_path / "container"
    entry = container / "entry"
    container.mkdir()
    write_bytes(entry, b"entry\n", 0o600)
    real_renameat2_between = authority.renameat2_between
    real_fsync = authority.os.fsync
    fsynced: list[int] = []

    def retire_then_raise(*args, **kwargs) -> None:
        real_renameat2_between(*args, **kwargs)
        raise RuntimeError("injected after namespace syscall")

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority, "renameat2_between", retire_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(container, "unlink parent") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="completed ambiguously.*directory fsync"):
            authority.unlink_bound_entry(parent, entry.name, profile, "entry")
        assert len(set(fsynced)) >= 2

    assert not entry.exists()
    assert len(list(tmp_path.glob(".cgl-source-authority-retired-*.forensic"))) == 1


def test_rmdir_retirement_post_syscall_exception_fsyncs_both_namespaces(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    container = tmp_path / "container"
    entry = container / "entry"
    container.mkdir()
    entry.mkdir()
    real_renameat2_between = authority.renameat2_between
    real_fsync = authority.os.fsync
    fsynced: list[int] = []

    def retire_then_raise(*args, **kwargs) -> None:
        real_renameat2_between(*args, **kwargs)
        raise RuntimeError("injected after namespace syscall")

    def count_fsync(descriptor: int) -> None:
        fsynced.append(descriptor)
        real_fsync(descriptor)

    monkeypatch.setattr(authority, "renameat2_between", retire_then_raise)
    monkeypatch.setattr(authority.os, "fsync", count_fsync)
    with authority.bound_directory(container, "rmdir parent") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="completed ambiguously.*directory fsync"):
            authority.rmdir_bound_entry(parent, entry.name, profile, "entry")
        assert len(set(fsynced)) >= 2

    assert not entry.exists()
    assert len(
        list(tmp_path.glob(".cgl-source-authority-retired-directory-*.forensic"))
    ) == 1


def test_recovery_unlink_rejects_entry_inode_substitution(tmp_path: Path) -> None:
    transaction = tmp_path / "transaction"
    transaction.mkdir(mode=0o700)
    entry = transaction / "journal.json"
    replacement = transaction / "replacement"
    write_bytes(entry, b"original\n", 0o600)
    write_bytes(replacement, b"replacement\n", 0o600)

    with authority.bound_directory(transaction, "transaction", mode=0o700) as (
        descriptor,
        _,
    ):
        _, _, profile = authority.read_bound_file(
            descriptor, entry.name, "transaction journal", mode=0o600
        )
        entry.unlink()
        replacement.rename(entry)
        with pytest.raises(ValueError, match="inode identity changed"):
            authority.unlink_bound_entry(
                descriptor, entry.name, profile, "transaction journal"
            )

    assert entry.read_bytes() == b"replacement\n"


@pytest.mark.parametrize("drift", ["mode", "link"])
def test_recovery_retirement_rejects_stale_file_profile(
    tmp_path: Path, drift: str
) -> None:
    entry = tmp_path / "partial"
    write_bytes(entry, b"partial", 0o000)
    retained_before = set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic"))

    with authority.bound_directory(tmp_path, "recovery parent") as (parent, _):
        observed = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        if drift == "mode":
            entry.chmod(0o666)
        else:
            os.link(entry, tmp_path / "partial.extra-link")
        with pytest.raises(ValueError, match="file security/content profile changed"):
            authority.retire_bound_recovery_file(
                parent, entry.name, observed, "deterministic recovery remnant"
            )

    assert entry.exists()
    assert set(tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic")) == retained_before


def test_forensic_namespace_rejects_unsafe_retained_file_profile(tmp_path: Path) -> None:
    forensic = tmp_path / f".cgl-source-authority-retired-{'a' * 32}.forensic"
    write_bytes(forensic, b"retained evidence\n", 0o666)

    with authority.bound_directory(tmp_path, "forensic namespace") as (parent, _):
        with pytest.raises(ValueError, match="must not be group- or world-writable"):
            authority.non_forensic_entry_names(parent, "forensic namespace")

    assert forensic.read_bytes() == b"retained evidence\n"


def test_forensic_namespace_rejects_unbounded_retained_link_profile(tmp_path: Path) -> None:
    forensic = tmp_path / f".cgl-source-authority-retired-{'a' * 32}.forensic"
    write_bytes(forensic, b"retained evidence\n", 0o600)
    os.link(forensic, tmp_path / "retained-link-one")
    os.link(forensic, tmp_path / "retained-link-two")

    with authority.bound_directory(tmp_path, "forensic namespace") as (parent, _):
        with pytest.raises(ValueError, match="invalid link profile"):
            authority.non_forensic_entry_names(parent, "forensic namespace")


def test_forensic_namespace_rejects_external_retained_hardlink(tmp_path: Path) -> None:
    forensic = tmp_path / f".cgl-source-authority-retired-{'a' * 32}.forensic"
    write_bytes(forensic, b"retained evidence\n", 0o600)
    os.link(forensic, tmp_path / "active-external-link")

    with authority.bound_directory(tmp_path, "forensic namespace") as (parent, _):
        with pytest.raises(ValueError, match="external retained hardlink"):
            authority.non_forensic_entry_names(parent, "forensic namespace")


def test_forensic_retirement_never_opens_unbound_dotdot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    container = tmp_path / "container"
    entry = container / "entry"
    container.mkdir()
    write_bytes(entry, b"authenticated entry\n", 0o600)
    real_open = authority.os.open

    def reject_dotdot(path: object, *args: object, **kwargs: object) -> int:
        if os.fspath(path) == "..":
            raise AssertionError("forensic retirement must not open unbound '..'")
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr(authority.os, "open", reject_dotdot)
    with authority.bound_directory(container, "retirement container") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        authority.unlink_bound_entry(parent, entry.name, profile, "authenticated entry")

    assert not entry.exists()
    assert len(list(tmp_path.glob(".cgl-source-authority-retired-*.forensic"))) == 1


def test_moved_retirement_parent_cannot_redirect_forensic_namespace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = tmp_path / "original"
    other = tmp_path / "other"
    container = original / "container"
    replacement = tmp_path / "replacement"
    detached = other / "detached"
    entry = container / "entry"
    original.mkdir()
    other.mkdir()
    container.mkdir()
    replacement.mkdir()
    write_bytes(entry, b"authenticated entry\n", 0o600)
    real_renameat2_between = authority.renameat2_between
    raced = False

    def move_parent_before_retirement(
        source_parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal raced
        if not raced and "retirement" in label:
            raced = True
            container.rename(detached)
            replacement.rename(container)
        real_renameat2_between(
            source_parent, source, target_parent, target, flags, label
        )

    monkeypatch.setattr(authority, "renameat2_between", move_parent_before_retirement)
    with pytest.raises(
        ValueError,
        match=(
            "durably ambiguous after authority loss|"
            "path changed while mutation is active"
        ),
    ):
        with authority.bound_directory(container, "retirement container") as (parent, _):
            profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
            authority.unlink_bound_entry(parent, entry.name, profile, "authenticated entry")

    assert raced
    assert not list(other.glob(".cgl-source-authority-retired-*.forensic"))
    assert not list(original.glob(".cgl-source-authority-retired-*.forensic"))
    assert (detached / entry.name).read_bytes() == b"authenticated entry\n"


def test_forensic_parent_security_drift_during_retirement_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    forensic_parent = tmp_path / "forensic-parent"
    container = forensic_parent / "container"
    entry = container / "entry"
    forensic_parent.mkdir()
    container.mkdir()
    write_bytes(entry, b"authenticated entry\n", 0o600)
    real_renameat2_between = authority.renameat2_between
    raced = False

    def drift_forensic_parent_before_retirement(
        source_parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal raced
        if not raced and "retirement" in label:
            raced = True
            forensic_parent.chmod(0o777)
        real_renameat2_between(
            source_parent, source, target_parent, target, flags, label
        )

    monkeypatch.setattr(
        authority, "renameat2_between", drift_forensic_parent_before_retirement
    )
    with pytest.raises(
        ValueError,
        match=(
            "durably ambiguous after authority loss|"
            "security metadata changed"
        ),
    ):
        with authority.bound_directory(container, "retirement container") as (parent, _):
            profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
            authority.unlink_bound_entry(parent, entry.name, profile, "authenticated entry")
    forensic_parent.chmod(0o755)

    assert raced
    retained = list(
        forensic_parent.glob(".cgl-source-authority-retired-*.forensic")
    )
    assert retained == []
    assert entry.read_bytes() == b"authenticated entry\n"


def test_absent_target_publication_uses_authenticated_source_descriptor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    escaped = tmp_path / "escaped"
    substitute = tmp_path / "substitute"
    target = tmp_path / "target"
    write_bytes(source, b"reviewed source\n", 0o600)
    write_bytes(substitute, b"raced source pathname\n", 0o600)
    real_renameat2 = authority.renameat2
    raced = False

    def race_before_atomic_rename(
        parent: int, source_name: str, target_name: str, flags: int, label: str
    ) -> None:
        nonlocal raced
        if not raced and flags == authority.RENAME_NOREPLACE:
            raced = True
            source.rename(escaped)
            substitute.rename(source)
        real_renameat2(parent, source_name, target_name, flags, label)

    monkeypatch.setattr(authority, "renameat2", race_before_atomic_rename)
    with authority.bound_directory(tmp_path, "publication parent") as (parent, _):
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="inode identity changed"):
            authority.publish_bound_file_noreplace(
                parent, source.name, target.name, source_profile, "reviewed publication"
            )

    assert raced
    assert target.read_bytes() == b"raced source pathname\n"
    assert escaped.read_bytes() == b"reviewed source\n"
    assert not source.exists()


def test_absent_target_publication_never_clobbers_raced_target(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    write_bytes(source, b"reviewed source\n", 0o600)
    real_renameat2 = authority.renameat2
    raced = False

    def race_target_creation(
        parent: int, source_name: str, target_name: str, flags: int, label: str
    ) -> None:
        nonlocal raced
        if not raced and flags == authority.RENAME_NOREPLACE:
            raced = True
            write_bytes(target, b"raced target\n", 0o600)
        real_renameat2(parent, source_name, target_name, flags, label)

    monkeypatch.setattr(authority, "renameat2", race_target_creation)
    with authority.bound_directory(tmp_path, "publication parent") as (parent, _):
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="target already exists"):
            authority.publish_bound_file_noreplace(
                parent, source.name, target.name, source_profile, "reviewed publication"
            )

    assert raced
    assert source.read_bytes() == b"reviewed source\n"
    assert target.read_bytes() == b"raced target\n"


def test_publication_rejects_detached_replaced_parent_after_fsync(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    publication = tmp_path / "publication"
    detached = tmp_path / "detached-publication"
    replacement = tmp_path / "replacement-publication"
    source = publication / "source"
    target = publication / "target"
    lock = tmp_path / ".stage-i.lock"
    publication.mkdir()
    replacement.mkdir()
    write_bytes(source, b"reviewed source\n", 0o600)
    write_bytes(lock, b"", 0o644)
    publication_identity = authority.profile_identity(os.stat(publication))
    real_fsync = authority.fsync_descriptor
    raced = False

    def replace_parent_after_publication_fsync(descriptor: int) -> None:
        nonlocal raced
        real_fsync(descriptor)
        if (
            not raced
            and authority.profile_identity(os.fstat(descriptor)) == publication_identity
            and target.exists()
            and not source.exists()
        ):
            raced = True
            publication.rename(detached)
            replacement.rename(publication)

    monkeypatch.setattr(authority, "fsync_descriptor", replace_parent_after_publication_fsync)
    with authority.stage_i_lock({"lock": lock}):
        with pytest.raises(ValueError, match="path changed while mutation is active"):
            with authority.bound_directory(publication, "publication parent") as (parent, _):
                source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
                authority.publish_bound_file_noreplace(
                    parent,
                    source.name,
                    target.name,
                    source_profile,
                    "reviewed publication",
                )

    assert raced
    assert not target.exists()
    assert list(publication.iterdir()) == []
    assert (detached / target.name).read_bytes() == b"reviewed source\n"


def test_bound_directory_allows_content_and_timestamp_changes(tmp_path: Path) -> None:
    with authority.bound_directory(tmp_path, "publication parent") as (parent, _):
        child = tmp_path / "content-change"
        write_bytes(child, b"content\n", 0o600)
        os.utime(tmp_path, None)
        authority.require_bound_directory_descriptor(parent)
        child.unlink()
        authority.require_bound_directory_descriptor(parent)


@pytest.mark.parametrize("changed_field", ["st_uid", "st_gid"])
def test_bound_directory_rejects_owner_or_group_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, changed_field: str
) -> None:
    with authority.bound_directory(tmp_path, "publication parent") as (parent, _):
        real_fstat = authority.os.fstat

        def drifted_fstat(descriptor: int) -> object:
            profile = real_fstat(descriptor)
            if descriptor != parent:
                return profile
            values = {
                "st_mode": profile.st_mode,
                "st_uid": profile.st_uid,
                "st_gid": profile.st_gid,
                "st_dev": profile.st_dev,
                "st_ino": profile.st_ino,
            }
            values[changed_field] += 1
            return SimpleNamespace(**values)

        with monkeypatch.context() as race:
            race.setattr(authority.os, "fstat", drifted_fstat)
            with pytest.raises(ValueError, match="security metadata changed"):
                authority.require_bound_directory_descriptor(parent)


def test_unlink_race_retires_and_retains_substituted_inode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entry = tmp_path / "entry"
    escaped = tmp_path / "escaped"
    replacement = tmp_path / "replacement"
    write_bytes(entry, b"authenticated entry\n", 0o600)
    write_bytes(replacement, b"raced replacement\n", 0o600)
    real_renameat2_between = authority.renameat2_between
    raced = False

    def race_before_retirement(
        parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal raced
        if not raced and flags == authority.RENAME_NOREPLACE and "retirement" in label:
            raced = True
            os.rename(
                source, escaped.name, src_dir_fd=parent, dst_dir_fd=parent
            )
            os.rename(
                replacement.name, source, src_dir_fd=parent, dst_dir_fd=parent
            )
        real_renameat2_between(parent, source, target_parent, target, flags, label)

    monkeypatch.setattr(authority, "renameat2_between", race_before_retirement)
    with authority.bound_directory(tmp_path, "retirement parent") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="changed during atomic retirement"):
            authority.unlink_bound_entry(parent, entry.name, profile, "authenticated entry")

    assert raced
    assert escaped.read_bytes() == b"authenticated entry\n"
    assert any(
        path.read_bytes() == b"raced replacement\n"
        for path in tmp_path.parent.iterdir()
        if path.is_file()
    )


def test_unlink_post_retirement_fsync_substitution_retains_both_inodes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    container = tmp_path / "container"
    entry = container / "entry"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped"
    container.mkdir()
    write_bytes(entry, b"authenticated entry\n", 0o600)
    write_bytes(substitute, b"raced substitute\n", 0o600)
    real_fsync = authority.fsync_descriptor
    raced = False

    def race_during_final_fsync(descriptor: int) -> None:
        nonlocal raced
        real_fsync(descriptor)
        retired = list(tmp_path.glob(".cgl-source-authority-retired-*.forensic"))
        if not raced and retired:
            raced = True
            retired[0].rename(escaped)
            substitute.rename(retired[0])

    def forbidden_unlink(*args, **kwargs) -> None:
        del args, kwargs
        raise AssertionError("authenticated retirement must never call unlink")

    monkeypatch.setattr(authority, "fsync_descriptor", race_during_final_fsync)
    monkeypatch.setattr(authority.os, "unlink", forbidden_unlink)
    with authority.bound_directory(container, "retirement container") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="changed during atomic retirement"):
            authority.unlink_bound_entry(parent, entry.name, profile, "authenticated entry")

    assert raced
    assert not entry.exists()
    assert escaped.read_bytes() == b"authenticated entry\n"
    retained = list(tmp_path.glob(".cgl-source-authority-retired-*.forensic"))
    assert len(retained) == 1
    assert retained[0].read_bytes() == b"raced substitute\n"


def test_existing_target_preexchange_race_retains_all_inodes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped"
    write_bytes(source, b"reviewed replacement\n", 0o600)
    write_bytes(target, b"reviewed predecessor\n", 0o600)
    write_bytes(substitute, b"raced target\n", 0o600)
    real_renameat2 = authority.renameat2
    raced = False

    def race_before_exchange(
        parent: int, source_name: str, target_name: str, flags: int, label: str
    ) -> None:
        nonlocal raced
        if not raced and flags == authority.RENAME_EXCHANGE:
            raced = True
            os.rename(
                target_name, escaped.name, src_dir_fd=parent, dst_dir_fd=parent
            )
            os.rename(
                substitute.name, target_name, src_dir_fd=parent, dst_dir_fd=parent
            )
        real_renameat2(parent, source_name, target_name, flags, label)

    monkeypatch.setattr(authority, "renameat2", race_before_exchange)
    with authority.bound_directory(tmp_path, "exchange parent") as (parent, _):
        source_profile = os.stat(source.name, dir_fd=parent, follow_symlinks=False)
        target_profile = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="inode identity changed"):
            authority.exchange_bound_entries(
                parent,
                source.name,
                target.name,
                source_profile,
                target_profile,
                "reviewed exchange",
            )

    assert raced
    assert target.read_bytes() == b"reviewed replacement\n"
    assert source.read_bytes() == b"raced target\n"
    assert escaped.read_bytes() == b"reviewed predecessor\n"


def test_atomic_write_final_fsync_target_substitution_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "journal.json"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped"
    write_bytes(target, b"reviewed predecessor\n", 0o600)
    write_bytes(substitute, b"raced target\n", 0o600)
    real_retire = authority.unlink_bound_entry
    real_fsync = authority.fsync_descriptor
    predecessor_retired = False
    raced = False

    def observe_predecessor_retirement(*args, **kwargs) -> None:
        nonlocal predecessor_retired
        real_retire(*args, **kwargs)
        predecessor_retired = True

    def race_during_final_fsync(descriptor: int) -> None:
        nonlocal raced
        real_fsync(descriptor)
        if predecessor_retired and not raced:
            raced = True
            target.rename(escaped)
            substitute.rename(target)

    monkeypatch.setattr(authority, "unlink_bound_entry", observe_predecessor_retirement)
    monkeypatch.setattr(authority, "fsync_descriptor", race_during_final_fsync)
    with pytest.raises(ValueError, match="inode identity changed|checksum differs"):
        authority.atomic_write(target, b"reviewed replacement\n", 0o600)

    assert raced
    assert target.read_bytes() == b"raced target\n"
    assert escaped.read_bytes() == b"reviewed replacement\n"
    assert any(
        path.read_bytes() == b"reviewed predecessor\n"
        for path in tmp_path.parent.glob(".cgl-source-authority-retired-*.forensic")
        if path.is_file()
    )


def test_atomic_write_rejects_detached_replaced_parent_after_fsync(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    accounting = tmp_path / "accounting"
    detached = tmp_path / "detached-accounting"
    replacement = tmp_path / "replacement-accounting"
    target = accounting / "journal.json"
    lock = tmp_path / ".stage-i.lock"
    accounting.mkdir()
    replacement.mkdir()
    write_bytes(target, b"reviewed predecessor\n", 0o600)
    write_bytes(lock, b"", 0o644)
    accounting_identity = authority.profile_identity(os.stat(accounting))
    real_fsync = authority.fsync_descriptor
    raced = False

    def replace_parent_after_atomic_write_fsync(descriptor: int) -> None:
        nonlocal raced
        real_fsync(descriptor)
        if (
            not raced
            and authority.profile_identity(os.fstat(descriptor)) == accounting_identity
            and target.exists()
            and target.read_bytes() == b"reviewed replacement\n"
        ):
            raced = True
            accounting.rename(detached)
            replacement.rename(accounting)

    monkeypatch.setattr(authority, "fsync_descriptor", replace_parent_after_atomic_write_fsync)
    with authority.stage_i_lock({"lock": lock}):
        with pytest.raises(ValueError, match="path changed while mutation is active"):
            authority.atomic_write(target, b"reviewed replacement\n", 0o600)

    assert raced
    assert not target.exists()
    assert list(accounting.iterdir()) == []
    assert (detached / target.name).read_bytes() == b"reviewed replacement\n"
    assert b"reviewed predecessor\n" in [
        path.read_bytes() for path in detached.iterdir() if path.is_file()
    ]


def test_atomic_write_rejects_parent_mode_drift_after_exchange(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    accounting = tmp_path / "accounting"
    target = accounting / "journal.json"
    lock = tmp_path / ".stage-i.lock"
    accounting.mkdir()
    write_bytes(target, b"reviewed predecessor\n", 0o600)
    write_bytes(lock, b"", 0o644)
    accounting_identity = authority.profile_identity(os.stat(accounting))
    real_fsync = authority.fsync_descriptor
    raced = False

    def change_parent_mode_after_exchange(descriptor: int) -> None:
        nonlocal raced
        real_fsync(descriptor)
        if (
            not raced
            and authority.profile_identity(os.fstat(descriptor)) == accounting_identity
            and target.exists()
            and target.read_bytes() == b"reviewed replacement\n"
        ):
            raced = True
            accounting.chmod(0o777)

    monkeypatch.setattr(authority, "fsync_descriptor", change_parent_mode_after_exchange)
    with authority.stage_i_lock({"lock": lock}):
        with pytest.raises(ValueError, match="security metadata changed"):
            authority.atomic_write(target, b"reviewed replacement\n", 0o600)

    assert raced
    assert target.read_bytes() == b"reviewed replacement\n"
    assert accounting.stat().st_mode & 0o777 == 0o777


def test_atomic_write_rollback_retains_substituted_temporary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "journal.json"
    substitute = tmp_path / "substitute"
    write_bytes(target, b"reviewed predecessor\n", 0o600)
    write_bytes(substitute, b"raced temporary\n", 0o600)
    rollback_unlinks: list[str] = []

    def race_then_fail(
        parent: int,
        source_name: str,
        target_name: str,
        source_expected: os.stat_result,
        target_expected: os.stat_result,
        label: str,
    ) -> None:
        del target_name, source_expected, target_expected, label
        escaped = f".{source_name}.escaped"
        os.rename(source_name, escaped, src_dir_fd=parent, dst_dir_fd=parent)
        os.rename(substitute.name, source_name, src_dir_fd=parent, dst_dir_fd=parent)
        raise ValueError("simulated exchange ambiguity")

    def reject_rollback_unlink(path, *args, **kwargs) -> None:
        del args, kwargs
        rollback_unlinks.append(os.fspath(path))
        raise AssertionError("ambiguous rollback attempted unlink")

    monkeypatch.setattr(authority, "exchange_bound_entries", race_then_fail)
    monkeypatch.setattr(authority.os, "unlink", reject_rollback_unlink)
    with pytest.raises(ValueError, match="simulated exchange ambiguity"):
        authority.atomic_write(target, b"reviewed replacement\n", 0o600)

    assert rollback_unlinks == []
    assert target.read_bytes() == b"reviewed predecessor\n"
    retained = [
        path.read_bytes()
        for path in tmp_path.iterdir()
        if path.is_file() and path != target
    ]
    assert b"reviewed replacement\n" in retained
    assert b"raced temporary\n" in retained


def test_recovery_rmdir_rejects_transaction_root_substitution(tmp_path: Path) -> None:
    transactions = tmp_path / "transactions"
    transactions.mkdir()
    retained = tmp_path / "retained-transactions"

    with pytest.raises(ValueError, match="path changed while mutation is active"):
        with authority.bound_directory(
            transactions, "source-authority transaction root"
        ) as (_, profile):
            transactions.rename(retained)
            transactions.mkdir()
            with pytest.raises(ValueError, match="inode identity changed"):
                authority.rmdir_bound_path(
                    transactions, profile, "source-authority transaction root"
                )

    assert transactions.is_dir()
    assert retained.is_dir()


def test_rmdir_race_retires_and_retains_substituted_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entry = tmp_path / "transaction"
    replacement = tmp_path / "replacement"
    escaped = tmp_path / "escaped"
    entry.mkdir()
    replacement.mkdir()
    write_bytes(replacement / "forensic-marker", b"raced directory\n", 0o600)
    real_renameat2_between = authority.renameat2_between
    raced = False

    def race_before_retirement(
        parent: int, source: str, target_parent: int, target: str,
        flags: int, label: str,
    ) -> None:
        nonlocal raced
        if not raced and flags == authority.RENAME_NOREPLACE and "retirement" in label:
            raced = True
            os.rename(
                source, escaped.name, src_dir_fd=parent, dst_dir_fd=parent
            )
            os.rename(
                replacement.name, source, src_dir_fd=parent, dst_dir_fd=parent
            )
        real_renameat2_between(parent, source, target_parent, target, flags, label)

    monkeypatch.setattr(authority, "renameat2_between", race_before_retirement)
    with authority.bound_directory(tmp_path, "transaction parent") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="changed during atomic retirement"):
            authority.rmdir_bound_entry(parent, entry.name, profile, "transaction")

    assert raced
    assert escaped.is_dir()
    assert any(
        (path / "forensic-marker").read_bytes() == b"raced directory\n"
        for path in tmp_path.parent.iterdir()
        if path.is_dir() and (path / "forensic-marker").is_file()
    )


def test_rmdir_post_retirement_fsync_substitution_retains_both_directories(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    container = tmp_path / "container"
    entry = container / "transaction"
    substitute = tmp_path / "substitute"
    escaped = tmp_path / "escaped"
    container.mkdir()
    entry.mkdir()
    substitute.mkdir()
    write_bytes(entry / "authenticated-marker", b"authenticated directory\n", 0o600)
    write_bytes(substitute / "substitute-marker", b"raced directory\n", 0o600)
    real_fsync = authority.fsync_descriptor
    raced = False

    def race_during_final_fsync(descriptor: int) -> None:
        nonlocal raced
        real_fsync(descriptor)
        retired = list(
            tmp_path.glob(".cgl-source-authority-retired-directory-*.forensic")
        )
        if not raced and retired:
            raced = True
            retired[0].rename(escaped)
            substitute.rename(retired[0])

    def forbidden_rmdir(*args, **kwargs) -> None:
        del args, kwargs
        raise AssertionError("authenticated retirement must never call rmdir")

    monkeypatch.setattr(authority, "fsync_descriptor", race_during_final_fsync)
    monkeypatch.setattr(authority.os, "rmdir", forbidden_rmdir)
    with authority.bound_directory(container, "transaction container") as (parent, _):
        profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
        with pytest.raises(ValueError, match="changed during atomic retirement"):
            authority.rmdir_bound_entry(parent, entry.name, profile, "transaction")

    assert raced
    assert not entry.exists()
    assert (escaped / "authenticated-marker").read_bytes() == b"authenticated directory\n"
    retained = list(tmp_path.glob(".cgl-source-authority-retired-directory-*.forensic"))
    assert len(retained) == 1
    assert (retained[0] / "substitute-marker").read_bytes() == b"raced directory\n"
