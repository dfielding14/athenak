"""Focused and adversarial tests for the authenticated Stage I wave planner."""

from __future__ import annotations

import csv
from datetime import datetime, timedelta, timezone
from decimal import Decimal
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import struct
import subprocess

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_wave_plan.py"
NOW = datetime.now(timezone.utc).replace(microsecond=0)


def load_utility():
    """Load the planner without executing its CLI."""

    spec = importlib.util.spec_from_file_location("cgl_lf_stage_i_wave_plan", UTILITY)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


planner = load_utility()


def utc_text(value: datetime = NOW) -> str:
    """Return one canonical UTC timestamp."""

    return value.isoformat()


def touch_now(path: Path, value: datetime = NOW) -> None:
    """Bind one fixture file mtime to the reviewed clock."""

    timestamp = value.timestamp()
    os.utime(path, (timestamp, timestamp))


def write_bytes(path: Path, value: bytes, mode: int = 0o644) -> None:
    """Write one deterministic fixture file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    path.chmod(mode)
    touch_now(path)


def write_text(path: Path, value: str, mode: int = 0o644) -> None:
    """Write one deterministic fixture text file."""

    write_bytes(path, value.encode("utf-8"), mode)


def write_json(path: Path, value: object) -> None:
    """Write stable fixture JSON."""

    if path.exists():
        path.chmod(0o644)
    write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def file_evidence(path: Path) -> dict[str, object]:
    """Return exact planner file evidence."""

    return planner.read_stable_regular_file(path, "fixture evidence")[1]


def sha256(path: Path) -> str:
    """Return one file digest."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(repository: Path, *arguments: str) -> str:
    """Run one successful deterministic fixture Git command."""

    result = subprocess.run(
        ["/usr/bin/git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def segment(index: int, start: str, target: str) -> str:
    """Return one exact Stage I rank-local segment identifier."""

    return f"s{index:02d}_rankio_t{start.replace('.', 'p')}_t{target.replace('.', 'p')}"


def restart_payload(value: Decimal) -> bytes:
    """Return one minimally loadable qualified restart fixture."""

    marker = format(float(value), ".17g")
    parameters = f"<time>\nrestart_time = {marker}\n<par_end>\n".encode()
    header = bytearray(planner.RESTART_MESH_HEADER_SIZE)
    struct.pack_into(planner.RESTART_TIME_FORMAT, header, planner.RESTART_TIME_OFFSET, float(value))
    return parameters + bytes(header) + b"rank-local-restart-payload\n"


def write_batch_script(path: Path) -> bytes:
    """Write one controller-style self-digesting fixture batch script."""

    template = (
        "#!/bin/bash\n"
        f"BATCH_SCRIPT_SHA256={planner.BATCH_SCRIPT_DIGEST_PLACEHOLDER}\n"
        "exit 0\n"
    )
    digest = hashlib.sha256(template.encode()).hexdigest()
    payload = template.replace(planner.BATCH_SCRIPT_DIGEST_PLACEHOLDER, digest).encode()
    write_bytes(path, payload, 0o750)
    return payload


class Campaign:
    """Retain one exact local analogue of the canonical read-only trust boundary."""

    def __init__(self, root: Path, frozen: Path, helper: Path):
        self.root = root
        self.frozen = frozen
        self.helper = helper
        self.manifests: dict[Path, dict[str, object]] = {}
        self.ledger: list[dict[str, str]] = []
        self.reservations: list[dict[str, object]] = []
        self.cumulative = Decimal("0")
        self.next_job = 1000
        self.r03_authorization_path: Path | None = None
        self.r17_ready = False
        self.live_jobs: dict[str, dict[str, object]] = {}
        self.live_scripts: dict[str, bytes] = {}
        self.profile_reviewed_utc = NOW
        self.profile_observed_utc = NOW
        self.estimates = {
            case_id: Decimal("184.888889") if case_id == "R17"
            else Decimal("2.888889") if case_id == "R16"
            else Decimal("23.111111")
            for case_id in planner.ALL_CASES
        }
        self.profiles = {
            case_id: {
                "case_id": case_id,
                "nodes": 8 if case_id == "R17" else 2 if case_id in {"R04", "R12", "R16"} else 1,
                "walltime": "02:00:00",
                "athena_walltime": "01:50:00",
                "next_increment": planner.decimal_text(planner.INITIAL_TARGETS[case_id]),
                "estimated_runtime_seconds": 1000,
                "rationale": f"reviewed fixture profile for {case_id}",
            }
            for case_id in planner.PROFILE_CASES
        }

    def manifest_path(self, case_id: str, segment_id: str) -> Path:
        """Return one exact canonical fixture manifest path."""

        return (
            self.root
            / "runs/mks24-stage-i"
            / planner.EXECUTION_EPOCH
            / case_id
            / segment_id
            / "manifest/prepared_run.json"
        )

    def write_restart_group(
        self, root: Path, count: int, final_time: Decimal, name: str
    ) -> dict[str, object]:
        """Write one complete rank-local restart group."""

        rank_files = []
        payload = restart_payload(final_time)
        for rank in range(count):
            path = root / f"rank_{rank:08d}" / name
            write_bytes(path, payload)
            rank_files.append(
                {"path": str(path), "sha256": sha256(path), "size_bytes": len(payload)}
            )
        return {
            "path": rank_files[0]["path"],
            "rank_files": rank_files,
            "sha256": rank_files[0]["sha256"],
            "size_bytes": len(payload),
            "storage": "per_rank",
        }

    def terminal_restart(
        self, case_id: str, segment_id: str, nodes: int, final_time: Decimal
    ) -> dict[str, object]:
        """Return deterministic retained terminal restart evidence."""

        output = (
            self.root
            / "runs/mks24-stage-i"
            / planner.EXECUTION_EPOCH
            / case_id
            / segment_id
            / "output/rst"
        )
        return self.write_restart_group(
            output,
            nodes * planner.RANKS_PER_NODE,
            final_time,
            f"{case_id}_{segment_id}.00001.rst",
        )

    def parent_record(self, parent_path: Path) -> dict[str, object]:
        """Return exact continuation-parent metadata."""

        parent = self.manifests[parent_path]
        inspection = parent["scientific_inspection"]
        terminal = inspection["terminal_restart"]
        return {
            "case_id": parent["run"]["case_id"],
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "execution_epoch": planner.EXECUTION_EPOCH,
            "final_time": float(inspection["final_time"]),
            "input_sha256": planner.EXPECTED_CASES[parent["run"]["case_id"]][3],
            "manifest": str(parent_path),
            "restart_files": [item["path"] for item in terminal["rank_files"]],
            "restart_sha256": terminal["sha256"],
            "restart_time": float(inspection["final_time"]),
            "result": parent["accounting"]["result"],
            "segment": parent["run"]["segment"],
        }

    def base_manifest(
        self,
        case_id: str,
        segment_id: str,
        state: str,
        nodes: int,
        walltime: str,
        athena_walltime: str,
        parent_path: Path | None,
    ) -> dict[str, object]:
        """Return one exact prepared-manifest fixture."""

        _, _, target = planner.parse_segment(segment_id, "fixture segment")
        requested = planner.walltime_seconds(walltime, "fixture walltime")
        command = {
            "production_utility": {
                "committed": True,
                "path": str(self.helper),
                "revision": planner.CONTROLLER_REVISION,
                "sha256": planner.CONTROLLER_SHA256,
            },
            "source_bundle": {
                "path": str(planner.source_bundle_path()),
                "sha256": planner.SOURCE_BUNDLE_SHA256,
                "verified_revisions": [planner.SOURCE_REVISION, planner.CONTROLLER_REVISION],
            },
            "input_revision": planner.SOURCE_REVISION,
            "input_sha256": planner.EXPECTED_CASES[case_id][3],
            "matrix_sha256": planner.MATRIX_SHA256,
            "executable_revision": planner.SOURCE_REVISION,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "build_manifest": str(planner.build_manifest_path()),
            "overrides": [f"time/tlim={planner.decimal_text(target)}"],
            "time_tlim_target": float(target),
            "athena_walltime": athena_walltime,
            "parent_segment": None,
            "source_restart_file": None,
            "restart_sha256": None,
        }
        if parent_path is not None:
            parent = self.parent_record(parent_path)
            command["parent_segment"] = parent
            command["source_restart_file"] = parent["restart_files"][0]
            command["restart_sha256"] = parent["restart_sha256"]
        batch_path = self.manifest_path(case_id, segment_id).parent / "cgl_lf_stage_i.sbatch"
        batch_payload = write_batch_script(batch_path)
        command["batch_script_sha256"] = planner.normalized_batch_script_sha256(batch_payload)
        return {
            "schema_version": 3,
            "execution_epoch": planner.EXECUTION_EPOCH,
            "state": state,
            "job_id": None,
            "prepared_utc": utc_text(),
            "project_root": str(self.root),
            "run": {
                "case_id": case_id,
                "case_name": planner.EXPECTED_CASES[case_id][0],
                "segment": segment_id,
                "resolution": planner.EXPECTED_CASES[case_id][2],
            },
            "allocation": {
                "nodes": nodes,
                "requested_walltime": walltime,
                "requested_seconds": requested,
                "reserved_node_hours": nodes * requested / 3600,
                "ranks_per_node": planner.RANKS_PER_NODE,
                "cpus_per_task": planner.CPUS_PER_TASK,
            },
            "command": command,
            "paths": {"batch_script": str(batch_path)},
            "scientific_inspection": None,
            "accounting": None,
        }

    def add_recorded(
        self,
        case_id: str,
        segment_id: str,
        *,
        final_time: str | None = None,
        result: str = "accepted",
        parent_path: Path | None = None,
        nodes: int = 1,
        walltime: str = "01:00:00",
        athena_walltime: str = "00:50:00",
        scheduler_state: str = "COMPLETED",
        exit_code: str = "0:0",
    ) -> Path:
        """Retain one recorded fixture segment."""

        path = self.manifest_path(case_id, segment_id)
        manifest = self.base_manifest(
            case_id, segment_id, "recorded", nodes, walltime, athena_walltime, parent_path
        )
        _, start, target = planner.parse_segment(segment_id, "fixture recorded segment")
        final = target if final_time is None else Decimal(final_time)
        elapsed = 600
        actual = Decimal(nodes * elapsed) / Decimal(3600)
        self.cumulative += actual
        job = str(self.next_job)
        self.next_job += 1
        row = {
            "execution_epoch": planner.EXECUTION_EPOCH,
            "job_id": job,
            "submitted_utc": utc_text(),
            "completed_utc": utc_text(),
            "case_id": case_id,
            "case_name": planner.EXPECTED_CASES[case_id][0],
            "segment": segment_id,
            "state": scheduler_state,
            "exit_code": exit_code,
            "nodes": str(nodes),
            "requested_walltime": walltime,
            "elapsed_seconds": str(elapsed),
            "reserved_node_hours": f"{Decimal(nodes * planner.walltime_seconds(walltime, 'fixture')) / Decimal(3600):.6f}",
            "actual_node_hours": f"{actual:.6f}",
            "cumulative_stage_i_node_hours": f"{self.cumulative:.6f}",
            "executable_revision": planner.SOURCE_REVISION,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "input_revision": planner.SOURCE_REVISION,
            "input_file": str(path.parent / "submitted_input.athinput"),
            "output_dir": str(path.parents[1] / "output"),
            "result": result,
            "notes": "reviewed fixture accounting",
        }
        manifest["job_id"] = job
        manifest["accounting"] = row
        manifest["scientific_inspection"] = {
            "required_time": float(target),
            "final_time": float(final),
            "terminal_restart_time": float(final),
            "terminal_restart": self.terminal_restart(case_id, segment_id, nodes, final),
            "clean_for_continuation": result in {"accepted", "clean_partial"},
        }
        self.manifests[path] = manifest
        self.ledger.append(row)
        self.reservations.append(
            {
                "execution_epoch": planner.EXECUTION_EPOCH,
                "manifest": str(path),
                "case_id": case_id,
                "case_name": planner.EXPECTED_CASES[case_id][0],
                "segment": segment_id,
                "nodes": nodes,
                "requested_walltime": walltime,
                "reserved_node_hours": nodes * planner.walltime_seconds(walltime, "fixture") / 3600,
                "state": "recorded",
                "prepared_utc": utc_text(),
                "job_id": job,
                "actual_node_hours": float(actual),
                "result": result,
            }
        )
        assert final > start
        return path

    def add_active(
        self,
        case_id: str,
        segment_id: str,
        *,
        state: str = "submitted",
        parent_path: Path | None = None,
        nodes: int | None = None,
        walltime: str = "02:00:00",
        athena_walltime: str = "01:50:00",
    ) -> Path:
        """Retain one active prepared or submitted fixture segment."""

        nodes = nodes if nodes is not None else self.profiles.get(case_id, {"nodes": 1})["nodes"]
        path = self.manifest_path(case_id, segment_id)
        manifest = self.base_manifest(
            case_id, segment_id, state, nodes, walltime, athena_walltime, parent_path
        )
        job = None
        if state == "submitted":
            job = str(self.next_job)
            self.next_job += 1
            manifest["job_id"] = job
            self.live_jobs[job] = {
                "job_id": job,
                "job_name": planner.expected_job_name(case_id, segment_id),
                "state": "PENDING",
                "nodes": nodes,
                "owner": planner.scheduler_owner(),
                "account": planner.ACCOUNT,
            }
            self.live_scripts[job] = Path(manifest["paths"]["batch_script"]).read_bytes()
        self.manifests[path] = manifest
        reservation = {
            "execution_epoch": planner.EXECUTION_EPOCH,
            "manifest": str(path),
            "case_id": case_id,
            "case_name": planner.EXPECTED_CASES[case_id][0],
            "segment": segment_id,
            "nodes": nodes,
            "requested_walltime": walltime,
            "reserved_node_hours": nodes * planner.walltime_seconds(walltime, "fixture") / 3600,
            "state": state,
            "prepared_utc": utc_text(),
        }
        if job is not None:
            reservation["job_id"] = job
        self.reservations.append(reservation)
        return path

    def add_cancelled(
        self,
        case_id: str,
        segment_id: str,
        *,
        parent_path: Path,
        nodes: int = 1,
        walltime: str = "02:00:00",
        athena_walltime: str = "01:50:00",
    ) -> Path:
        """Retain one cancelled no-start operational identity."""

        path = self.manifest_path(case_id, segment_id)
        manifest = self.base_manifest(
            case_id, segment_id, "cancelled", nodes, walltime, athena_walltime, parent_path
        )
        job = str(self.next_job)
        self.next_job += 1
        manifest["job_id"] = job
        manifest["command"]["source_bundle"] = {
            "path": str(planner.expected_path(planner.R03_F114_SOURCE_BUNDLE)),
            "sha256": planner.R03_F114_SOURCE_BUNDLE_SHA256,
            "verified_revisions": [
                planner.SOURCE_REVISION, planner.R03_F114_CONTROLLER_REVISION
            ],
        }
        self.manifests[path] = manifest
        write_json(path, manifest)
        self.reservations.append(
            {
                "execution_epoch": planner.EXECUTION_EPOCH,
                "manifest": str(path),
                "case_id": case_id,
                "case_name": planner.EXPECTED_CASES[case_id][0],
                "segment": segment_id,
                "nodes": nodes,
                "requested_walltime": walltime,
                "reserved_node_hours": (
                    nodes * planner.walltime_seconds(walltime, "fixture") / 3600
                ),
                "state": "cancelled",
                "prepared_utc": utc_text(),
                "job_id": job,
                "notes": "cancelled before start; operational identity is not reusable",
            }
        )
        return path

    def live_scheduler_job(self, job_id: str) -> dict[str, object]:
        """Return one fixture live Slurm row."""

        if job_id not in self.live_jobs:
            raise ValueError(f"fixture Slurm has no active job {job_id}")
        return dict(self.live_jobs[job_id])

    def live_batch_script(self, job_id: str) -> bytes:
        """Return fixture Slurm-stored script bytes."""

        if job_id not in self.live_scripts:
            raise ValueError(f"fixture Slurm has no stored script for {job_id}")
        return self.live_scripts[job_id]

    def authorize_r03(self, parent_path: Path) -> None:
        """Retain the exact immutable published F-115 authority chain."""

        parent = self.manifests[parent_path]
        endpoint = parent["scientific_inspection"]["final_time"]
        restart = parent["scientific_inspection"]["terminal_restart"]
        cancelled_path = self.manifest_path("R03", planner.R03_F114_SEGMENT)
        cancelled = self.manifests[cancelled_path]
        self.r03_authorization_path = planner.r03_f115_path()
        sole = {
            "athena_walltime": "01:50:00",
            "case_id": "R03",
            "cpus_per_task": planner.CPUS_PER_TASK,
            "executable": str(planner.executable_path()),
            "executable_revision": planner.SOURCE_REVISION,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "matrix": str(planner.frozen_matrix_path()),
            "matrix_sha256": planner.MATRIX_SHA256,
            "nodes": 1,
            "override": "time/tlim=0.5",
            "parent_job_id": parent["job_id"],
            "parent_result": "clean_partial",
            "parent_segment": parent["run"]["segment"],
            "ranks_per_node": planner.RANKS_PER_NODE,
            "restart_file": restart["path"],
            "restart_time": endpoint,
            "segment": planner.R03_F115_SEGMENT,
            "source_bundle": str(planner.source_bundle_path()),
            "source_bundle_sha256": planner.SOURCE_BUNDLE_SHA256,
            "source_dir": str(planner.FROZEN_SOURCE),
            "time_tlim_target": 0.5,
            "walltime": "02:00:00",
        }
        supersession = {
            "segment": {
                "from": planner.R03_F114_SEGMENT,
                "reason": "s01 is an immutable cancelled no-start manifest and directory",
                "to": planner.R03_F115_SEGMENT,
            },
            "source_bundle": {
                "from": str(planner.expected_path(planner.R03_F114_SOURCE_BUNDLE)),
                "to": str(planner.source_bundle_path()),
            },
            "source_bundle_sha256": {
                "from": planner.R03_F114_SOURCE_BUNDLE_SHA256,
                "to": planner.SOURCE_BUNDLE_SHA256,
            },
        }
        cancelled_evidence = file_evidence(cancelled_path)
        write_json(
            self.r03_authorization_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-source-bundle-recovery-supersession-evidence",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "checkpoint": "F-115",
                "implementation": {
                    "source_bundle": {
                        "complete_history": True,
                        "head": planner.CONTROLLER_REVISION,
                        "links": 1,
                        "mode": "0644",
                        "path": planner.source_bundle_path().relative_to(
                            planner.CANONICAL_ROOT
                        ).as_posix(),
                        "sha256": planner.SOURCE_BUNDLE_SHA256,
                        "verified_revisions": list(
                            planner.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS
                        ),
                    },
                    "stage_i_helper": {
                        "links": 1,
                        "mode": "0644",
                        "path": str(self.helper),
                        "sha256": planner.CONTROLLER_SHA256,
                    },
                },
                "authorization": {
                    "sole_next_segment_profile": sole,
                    "supersedes_f114_profile_only_where_explicitly_listed": supersession,
                },
                "cancelled_submission": {
                    "allocated_nodes": 0,
                    "elapsed_seconds": 0,
                    "evidence": {
                        "live_cancelled_manifest": {
                            "path": str(cancelled_path),
                            "sha256": cancelled_evidence["sha256"],
                            "mode": cancelled_evidence["mode"],
                            "links": 1,
                        }
                    },
                    "exit_code": "0:0",
                    "job_id": cancelled["job_id"],
                    "reusable": False,
                    "state": "CANCELLED",
                },
                "publication_requirements": [
                    "Do not reuse cancelled job or the retained s01 directory.",
                    "Prepare only the sole s02 R03 continuation profile.",
                    "Use controller paths; do not invoke direct sbatch.",
                ],
                "scope": {
                    "relationship": (
                        "source-bundle-binding-and-cancelled-segment-identity-supersession"
                    )
                },
            },
        )
        self.r03_authorization_path.chmod(0o444)
        planner.R03_F115_SHA256 = sha256(self.r03_authorization_path)
        published = {
            "path": str(self.r03_authorization_path),
            "sha256": planner.R03_F115_SHA256,
        }

        review_specs = (
            (
                planner.r03_f115_provenance_security_review_path(),
                "provenance-security",
                "approved-for-publication",
                planner.F115_PROVENANCE_REVIEWER,
                "R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256",
            ),
            (
                planner.r03_f115_plasma_scientific_review_path(),
                "plasma-scientific-continuation",
                "approved",
                planner.F115_PLASMA_REVIEWER,
                "R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256",
            ),
        )
        for review_path, kind, decision, reviewer, digest_name in review_specs:
            authority = (
                {"scope": list(planner.F115_PROVENANCE_REVIEW_SCOPE)}
                if kind == "provenance-security"
                else {
                    "verified": {
                        "authorization_limitations": (
                            planner.F115_PLASMA_AUTHORIZATION_LIMITATIONS
                        )
                    }
                }
            )
            write_json(
                review_path,
                {
                    **authority,
                    "schema_version": 1,
                    "record_type": (
                        "stage-i-source-bundle-recovery-supersession-independent-review"
                    ),
                    "execution_epoch": planner.EXECUTION_EPOCH,
                    "checkpoint": "F-115",
                    "decision": decision,
                    "published_f115": published,
                    "review_kind": kind,
                    "reviewed_candidate": {
                        "path": f"/tmp/fixture-{kind}.candidate",
                        "sha256": planner.R03_F115_SHA256,
                    },
                    "reviewed_utc": utc_text(),
                    "reviewer": reviewer,
                },
            )
            review_path.chmod(0o444)
            setattr(planner, digest_name, sha256(review_path))

        audit_path = planner.r03_f115_publication_audit_path()
        write_json(
            audit_path,
            {
                "schema_version": 1,
                "record_type": (
                    "stage-i-source-bundle-recovery-supersession-publication-audit"
                ),
                "execution_epoch": planner.EXECUTION_EPOCH,
                "checkpoint": "F-115",
                "audit_generated_utc": utc_text(),
                "published_utc": utc_text(),
                "publication": (
                    "atomic-write-fsync-rename-fsync-under-canonical-stage-i-lock"
                ),
                "artifact": {
                    "path": str(self.r03_authorization_path),
                    "sha256": planner.R03_F115_SHA256,
                    "mode": "0444",
                    "links": 1,
                },
                "authority_and_enforcement": {
                    "authorization_kind": "procedural-pre-prepare-sole-profile-authority",
                    "direct_sbatch_authorized": False,
                    "enforcement_chain": [
                        "F115 authorizes only the exact embedded profile before prepare.",
                        "Prepared manifest becomes machine-enforced execution intent.",
                    ],
                    "f115_authority": {
                        "path": str(self.r03_authorization_path),
                        "sha256": planner.R03_F115_SHA256,
                        "mode": "0444",
                        "links": 1,
                    },
                    "reuse_cancelled_job_or_s01_authorized": False,
                    "shared_root_acknowledgement_authorized_by_f115": False,
                    (
                        "shared_root_acknowledgement_requires_separate_exact_"
                        "isolation_review_after_prepare"
                    ): True,
                    "sole_next_segment_profile": sole,
                },
                "independent_reviews": {
                    "plasma_scientific_continuation": {
                        "path": str(planner.r03_f115_plasma_scientific_review_path()),
                        "sha256": planner.R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
                        "mode": "0444",
                        "links": 1,
                    },
                    "provenance_security": {
                        "path": str(planner.r03_f115_provenance_security_review_path()),
                        "sha256": planner.R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256,
                        "mode": "0444",
                        "links": 1,
                    },
                    "reviews_bind_exact_published_f115_sha256": planner.R03_F115_SHA256,
                },
                "reproducible_implementation_authority": {
                    "authoritative_source_bundle": {
                        "complete_history": True,
                        "head": planner.CONTROLLER_REVISION,
                        "links": 1,
                        "mode": "0644",
                        "path": str(planner.source_bundle_path()),
                        "sha256": planner.SOURCE_BUNDLE_SHA256,
                    },
                    "committed_stage_i_helper": {
                        "path": "scripts/frontier/cgl_lf_stage_i.py",
                        "sha256": planner.CONTROLLER_SHA256,
                    },
                },
                "source_archive_catalog": {
                    "corrupt_c7_absent_from_active_checksum_ledger": True,
                    "corrupt_c7_retained_as_incident_evidence": True,
                    "full_active_checksum_ledger": "passed",
                    "new_bundle_present_exactly_once": True,
                },
            },
        )
        audit_path.chmod(0o444)
        planner.R03_F115_PUBLICATION_AUDIT_SHA256 = sha256(audit_path)

    def complete_predecessors(self) -> None:
        """Retain authenticated minimal R02-R16 lineages through exact t=10."""

        self.manifests.clear()
        self.ledger.clear()
        self.reservations.clear()
        self.live_jobs.clear()
        self.live_scripts.clear()
        self.cumulative = Decimal("0")
        r02 = self.add_recorded("R02", segment(0, "0", "10"))
        assert r02
        r03_0 = self.add_recorded(
            "R03", segment(0, "0", "0.5"), final_time="0.31282347945569927",
            result="clean_partial",
        )
        self.add_cancelled("R03", planner.R03_F114_SEGMENT, parent_path=r03_0)
        self.authorize_r03(r03_0)
        r03_2 = self.add_recorded(
            "R03", planner.R03_F115_SEGMENT, parent_path=r03_0,
            walltime="02:00:00", athena_walltime="01:50:00",
        )
        self.add_recorded(
            "R03", segment(3, "0.5", "10"), parent_path=r03_2
        )
        self.r03_authorization_path = None
        for case_id in planner.LOWER_CASES:
            first = planner.decimal_text(planner.INITIAL_TARGETS[case_id])
            parent = self.add_recorded(case_id, segment(0, "0", first), nodes=1)
            self.add_recorded(case_id, segment(1, first, "10"), parent_path=parent, nodes=1)

    def write_static_provenance(self) -> None:
        """Retain exact frozen/provenance fixture files and patch known digests."""

        repository = self.helper.parents[2]
        repository.mkdir(parents=True, exist_ok=True)
        git(repository, "init")
        git(repository, "config", "user.name", "Fixture Reviewer")
        git(repository, "config", "user.email", "fixture@example.invalid")
        write_text(self.helper, "qualified source revision\n")
        git(repository, "add", "scripts/frontier/cgl_lf_stage_i.py")
        git(repository, "commit", "-m", "fixture qualified source")
        planner.SOURCE_REVISION = git(repository, "rev-parse", "HEAD")

        write_text(self.helper, "historical F114 controller helper\n")
        git(repository, "add", "scripts/frontier/cgl_lf_stage_i.py")
        git(repository, "commit", "-m", "fixture historical controller")
        planner.R03_F114_CONTROLLER_REVISION = git(repository, "rev-parse", "HEAD")

        write_text(self.helper, "promoted F115 controller helper\n")
        git(repository, "add", "scripts/frontier/cgl_lf_stage_i.py")
        git(repository, "commit", "-m", "fixture promoted controller")
        planner.CONTROLLER_REVISION = git(repository, "rev-parse", "HEAD")
        planner.CONTROLLER_SHA256 = sha256(self.helper)
        planner.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS = (
            planner.SOURCE_REVISION,
            planner.R03_F114_CONTROLLER_REVISION,
            planner.CONTROLLER_REVISION,
        )
        bundle = planner.source_bundle_path()
        bundle.parent.mkdir(parents=True, exist_ok=True)
        git(repository, "bundle", "create", str(bundle), "HEAD")
        bundle.chmod(0o644)
        touch_now(bundle)
        planner.SOURCE_BUNDLE_SHA256 = sha256(bundle)
        executable = planner.executable_path()
        write_text(executable, "qualified executable\n", 0o755)
        planner.EXECUTABLE_SHA256 = sha256(executable)
        build_values = {
            "athena.sha256": f"{planner.EXECUTABLE_SHA256}  {executable}\n",
            "athena-config.txt": "MPI parallelism: ON\n",
            "environment.txt": f"git_revision={planner.SOURCE_REVISION}\n",
            "CMakeCache.txt": "CMAKE_BUILD_TYPE=Release\n",
        }
        build_digests = {}
        for filename, value in build_values.items():
            path = planner.build_manifest_path() / filename
            write_text(path, value)
            build_digests[filename] = sha256(path)
        planner.BUILD_FILE_SHA256 = build_digests

        cases = []
        updated = dict(planner.EXPECTED_CASES)
        for case_id, (name, filename, resolution, _) in planner.EXPECTED_CASES.items():
            path = self.frozen / "inputs/cgl_lf_paper" / filename
            write_text(path, f"frozen input {case_id}\n")
            updated[case_id] = (name, filename, resolution, sha256(path))
            cases.append(
                {
                    "id": case_id,
                    "name": name,
                    "input": f"inputs/cgl_lf_paper/{filename}",
                    "resolution": resolution,
                    "estimated_node_hours": float(self.estimates[case_id]),
                }
            )
        planner.EXPECTED_CASES = updated
        matrix = {
            "schema_version": 1,
            "campaign": planner.CAMPAIGN,
            "cases": cases,
        }
        write_json(planner.frozen_matrix_path(), matrix)
        planner.MATRIX_SHA256 = sha256(planner.frozen_matrix_path())
        write_json(
            planner.qualification_path(),
            {
                "schema_version": 1,
                "execution_epoch": planner.EXECUTION_EPOCH,
                "approved_executable": str(planner.executable_path()),
                "approved_executable_revision": planner.SOURCE_REVISION,
                "approved_executable_sha256": planner.EXECUTABLE_SHA256,
                "build_manifest": str(planner.build_manifest_path()),
                "approved_utc": utc_text(),
                "approved_by": "fixture reviewer",
                "review_notes": "fixture reviewed build",
            },
        )
        for path in planner.transaction_store_paths():
            path.mkdir(parents=True, exist_ok=True)

    def summary_text(self) -> str:
        """Render the exact controller summary consumed by the planner."""

        actual = sum(Decimal(row["actual_node_hours"]) for row in self.ledger)
        active = [
            reservation
            for reservation in self.reservations
            if reservation["state"] in {"prepared", "submitted"}
        ]
        reserved = sum(Decimal(str(item["reserved_node_hours"])) for item in active)
        lines = [
            "# MKS24 Stage I Frontier E03-forcing-policy Budget",
            "",
            f"- Updated UTC: `{utc_text(self.profile_observed_utc)}`",
            f"- Execution epoch: `{planner.EXECUTION_EPOCH}`",
            f"- E03-forcing-policy Stage I actual use: `{actual:.6f}` node-hours",
            f"- Active segment reservations: `{reserved:.6f}` node-hours",
            f"- Unreserved E03 mapped-matrix remainder: `{planner.STAGE_I_BUDGET_NODE_HOURS - actual - reserved:.6f}` node-hours",
            f"- Incremental project remainder after active E03 Stage I use: `{planner.PROJECT_BUDGET_NODE_HOURS - actual - reserved:.6f}` node-hours",
            "",
            "## Recorded Segments",
            "",
        ]
        if self.ledger:
            lines.extend(
                [
                    "| Job | Case/segment | State | Node-hours | Result |",
                    "| --- | --- | --- | ---: | --- |",
                ]
            )
            for row in self.ledger:
                lines.append(
                    f"| `{row['job_id']}` | `{row['case_id']}/{row['segment']}` | "
                    f"{row['state']} | `{Decimal(row['actual_node_hours']):.6f}` | "
                    f"{row['result']} |"
                )
        lines.extend(["", "## Active Reservations", ""])
        active_by_path = sorted(active, key=lambda item: item["manifest"])
        if active_by_path:
            lines.extend(
                [
                    "| Case/segment | Nodes | Walltime | Node-hours | State |",
                    "| --- | ---: | --- | ---: | --- |",
                ]
            )
            for item in active_by_path:
                lines.append(
                    f"| `{item['case_id']}/{item['segment']}` | `{item['nodes']}` | "
                    f"`{item['requested_walltime']}` | "
                    f"`{Decimal(str(item['reserved_node_hours'])):.6f}` | {item['state']} |"
                )
        lines.extend(
            [
                "",
                "Up to 4 distinct R03-R16 cases may be active concurrently, with at most "
                "one prepared submission packet. R17 remains exclusive and last.",
                "",
            ]
        )
        return "\n".join(lines)

    def projection(self, additional: dict[str, Decimal] | None = None) -> dict[str, object]:
        """Reproduce the planner's full remaining-campaign budget fixture."""

        endpoints = {case_id: Decimal("0") for case_id in planner.ALL_CASES}
        for manifest in self.manifests.values():
            if manifest["state"] != "recorded":
                continue
            case_id = manifest["run"]["case_id"]
            final = Decimal(str(manifest["scientific_inspection"]["final_time"]))
            endpoints[case_id] = max(endpoints[case_id], final)
        reserved = {}
        for item in self.reservations:
            if item["state"] in {"prepared", "submitted"}:
                case_id = item["case_id"]
                reserved[case_id] = reserved.get(case_id, Decimal("0")) + Decimal(
                    str(item["reserved_node_hours"])
                )
        for case_id, value in (additional or {}).items():
            reserved[case_id] = reserved.get(case_id, Decimal("0")) + value
        actual = sum((Decimal(row["actual_node_hours"]) for row in self.ledger), Decimal("0"))
        remaining = Decimal("0")
        breakdown = {}
        for case_id in planner.ALL_CASES:
            estimate = self.estimates[case_id]
            progress = min(Decimal("1"), max(Decimal("0"), endpoints[case_id] / Decimal("10")))
            scaled = estimate * (Decimal("1") - progress)
            committed = reserved.get(case_id, Decimal("0"))
            projected = max(scaled, committed)
            remaining += projected
            breakdown[case_id] = {
                "matrix_full_case_node_hours": planner.decimal_text(estimate),
                "authenticated_progress_fraction": planner.decimal_text(progress),
                "scaled_remaining_node_hours": planner.decimal_text(scaled),
                "committed_profile_reserved_node_hours": planner.decimal_text(committed),
                "projected_remaining_node_hours": planner.decimal_text(projected),
            }
        total = actual + remaining
        return {
            "method": (
                "authenticated ledger actuals plus frozen-matrix full-case estimates scaled "
                "by byte-authenticated terminal-lineage progress and floored by committed "
                "active/planned reservations"
            ),
            "actual_stage_i_node_hours": planner.decimal_text(actual),
            "computed_remaining_stage_i_node_hours": planner.decimal_text(remaining),
            "computed_stage_i_total_node_hours": planner.decimal_text(total),
            "promoted_stage_i_envelope_node_hours": planner.decimal_text(
                planner.STAGE_I_BUDGET_NODE_HOURS
            ),
            "project_ceiling_node_hours": planner.decimal_text(planner.PROJECT_BUDGET_NODE_HOURS),
            "computed_stage_i_margin_node_hours": planner.decimal_text(
                planner.STAGE_I_BUDGET_NODE_HOURS - total
            ),
            "case_breakdown": breakdown,
        }

    def write_allocation_authority(self) -> dict[str, dict[str, object]]:
        """Publish an independently reviewed exact allocation authority fixture."""

        profiles = []
        for case_id in sorted(self.profiles):
            value = dict(self.profiles[case_id])
            value["next_increment"] = planner.decimal_text(
                planner.decimal_value(value["next_increment"], "fixture next increment")
            )
            profiles.append(value)
        artifact_path = planner.allocation_authority_path()
        artifact = {
            "schema_version": 1,
            "record_type": "cgl_lf_stage_i_independent_wave_allocation_authority",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "canonical_root": str(self.root),
            "generated_utc": utc_text(),
            "expires_utc": utc_text(NOW + timedelta(hours=1)),
            "matrix_sha256": planner.MATRIX_SHA256,
            "controller_revision": planner.CONTROLLER_REVISION,
            "controller_sha256": planner.CONTROLLER_SHA256,
            "profiles_sha256": planner.sha256_bytes(planner.canonical_json(profiles)),
            "profiles": profiles,
            "projection": self.projection(),
        }
        write_json(artifact_path, artifact)
        audit_path = planner.allocation_authority_audit_path()
        write_json(
            audit_path,
            {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_wave_allocation_authority_publication_audit",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "published_utc": utc_text(),
                "artifact": {
                    "path": str(artifact_path),
                    "sha256": sha256(artifact_path),
                    "mode": "0644",
                },
                "review": {
                    "status": "approved",
                    "reviewed_by": "independent fixture allocation reviewer",
                    "independent_from_profile_author": True,
                },
            },
        )
        return {
            "artifact": file_evidence(artifact_path),
            "publication_audit": file_evidence(audit_path),
        }

    def write_r17_evidence(self) -> dict[str, dict[str, object]]:
        """Retain exact evidence-backed R17 readiness artifacts."""

        if not self.r17_ready:
            return {}
        r17_reserved = Decimal(
            self.profiles["R17"]["nodes"]
            * planner.walltime_seconds(self.profiles["R17"]["walltime"], "fixture R17")
        ) / Decimal(3600)
        readiness_root = (
            self.root
            / "accounting/mks24_stage_i_E03_forcing_policy_R17_64_rank_readiness"
        )
        job_id = "900001"
        elapsed = 600
        scheduler_path = self.root / f"accounting/{job_id}.stage_i.sacct.txt"
        write_text(
            scheduler_path,
            f"{job_id}|cgl_mks24_r17_64_rank_readiness|COMPLETED|0:0|8|{elapsed}|"
            f"{utc_text()}|{utc_text()}\n",
        )
        batch_path = readiness_root / "r17_64_rank_readiness.sbatch"
        write_batch_script(batch_path)
        output_inventory = []
        for rank in range(64):
            path = readiness_root / "output" / f"rank_{rank:08d}" / "readiness.bin"
            write_bytes(path, f"rank {rank} measured output\n".encode())
            output_inventory.append(file_evidence(path))
        final_time = Decimal("0.01")
        terminal = self.write_restart_group(
            readiness_root / "restarts", 64, final_time, "r17_readiness.00001.rst"
        )
        values = {
            "storage": {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_r17_storage_readiness",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "generated_utc": utc_text(),
                "reviewed_utc": utc_text(),
                "reviewed_by": "storage reviewer",
                "decision": "approved",
                "available_bytes": planner.R17_REQUIRED_RETENTION_BYTES * 2,
                "required_retention_bytes": planner.R17_REQUIRED_RETENTION_BYTES,
            },
            "recost": {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_r17_node_hour_recost",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "generated_utc": utc_text(),
                "reviewed_utc": utc_text(),
                "reviewed_by": "recost reviewer",
                "decision": "approved",
                "projection": self.projection({"R17": r17_reserved}),
                "projection_sha256": planner.sha256_bytes(
                    planner.canonical_json(self.projection({"R17": r17_reserved}))
                ),
                "projected_r17_reserved_node_hours": planner.decimal_text(r17_reserved),
            },
            "rank_readiness": {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_r17_64_rank_readiness",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "generated_utc": utc_text(),
                "reviewed_utc": utc_text(),
                "reviewed_by": "64-rank reviewer",
                "decision": "approved",
                "nodes": 8,
                "ranks": 64,
                "job_id": job_id,
                "final_time": planner.decimal_text(final_time),
                "scheduler_evidence": file_evidence(scheduler_path),
                "batch_script": file_evidence(batch_path),
                "output_inventory": output_inventory,
                "terminal_restart": terminal,
                "performance": {
                    "elapsed_seconds": elapsed,
                    "node_hours": planner.decimal_text(Decimal(8 * elapsed) / Decimal(3600)),
                    "simulated_time": planner.decimal_text(final_time),
                    "meshblocks_per_rank": 8,
                },
            },
        }
        paths = {
            "storage": self.root / "accounting/mks24_stage_i_E03_forcing_policy_R17_storage_readiness.json",
            "recost": self.root / "accounting/mks24_stage_i_E03_forcing_policy_R17_node_hour_recost.json",
            "rank_readiness": self.root / "accounting/mks24_stage_i_E03_forcing_policy_R17_64_rank_readiness.json",
        }
        for key, path in paths.items():
            write_json(path, values[key])
        return {key: file_evidence(path) for key, path in paths.items()}

    def refresh(self) -> None:
        """Rewrite exact canonical stores, snapshot, and fresh reviewed profile."""

        for path, manifest in self.manifests.items():
            write_json(path, manifest)
        ledger_buffer = io.StringIO()
        writer = csv.DictWriter(ledger_buffer, fieldnames=planner.LEDGER_COLUMNS)
        writer.writeheader()
        writer.writerows(self.ledger)
        write_text(planner.ledger_path(), ledger_buffer.getvalue())
        write_json(planner.reservations_path(), self.reservations)
        write_text(planner.summary_path(), self.summary_text())

        manifest_paths = sorted(self.manifests)
        manifest_evidence = [file_evidence(path) for path in manifest_paths]
        actual_evidence = {
            "summary": file_evidence(planner.summary_path()),
            "ledger": file_evidence(planner.ledger_path()),
            "reservations": file_evidence(planner.reservations_path()),
            "qualification": file_evidence(planner.qualification_path()),
            "manifests": manifest_evidence,
        }
        active = [
            item for item in self.reservations if item["state"] in {"prepared", "submitted"}
        ]
        reconciliation = {
            "schema_version": 1,
            "record_type": "cgl_lf_stage_i_controller_reconciliation_snapshot",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "generated_utc": utc_text(),
            "observed_utc": utc_text(self.profile_observed_utc),
            "generator": {
                "path": str(self.helper),
                "revision": planner.CONTROLLER_REVISION,
                "sha256": planner.CONTROLLER_SHA256,
            },
            "evidence": actual_evidence,
            "report": {
                "execution_epoch": planner.EXECUTION_EPOCH,
                "root": str(self.root),
                "qualification": {
                    "state": "approved",
                    "path": str(planner.qualification_path()),
                    "sha256": actual_evidence["qualification"]["sha256"],
                    "approved_executable_revision": planner.SOURCE_REVISION,
                    "approved_executable_sha256": planner.EXECUTABLE_SHA256,
                },
                "consistent": True,
                "counts": {
                    "transactions": 0,
                    "reservations": len(self.reservations),
                    "active_reservations": len(active),
                    "ledger_rows": len(self.ledger),
                    "manifests": len(self.manifests),
                },
                "issues": [],
            },
        }
        write_json(planner.reconciliation_path(), reconciliation)

        authority_evidence = self.write_allocation_authority()
        r17_evidence = self.write_r17_evidence()
        input_bindings = {
            case_id: file_evidence(planner.input_path(case_id)) for case_id in planner.ALL_CASES
        }
        build_bindings = {
            filename: file_evidence(planner.build_manifest_path() / filename)
            for filename in planner.BUILD_FILE_SHA256
        }
        profile = {
            "schema_version": 2,
            "record_type": "cgl_lf_stage_i_authenticated_wave_profile",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "canonical_root": str(self.root),
            "review": {
                "decision": "approved",
                "reviewed_by": "fixture wave reviewer",
                "reviewed_utc": utc_text(self.profile_reviewed_utc),
                "state_observed_utc": utc_text(self.profile_observed_utc),
                "notes": "reviewed exact fixture state",
            },
            "provenance": {
                "controller_revision": planner.CONTROLLER_REVISION,
                "controller_sha256": planner.CONTROLLER_SHA256,
                "source_revision": planner.SOURCE_REVISION,
                "matrix_sha256": planner.MATRIX_SHA256,
                "source_bundle_sha256": planner.SOURCE_BUNDLE_SHA256,
                "executable_revision": planner.SOURCE_REVISION,
                "executable_sha256": planner.EXECUTABLE_SHA256,
                "build_manifest": str(planner.build_manifest_path()),
            },
            "evidence": {
                "matrix": file_evidence(planner.frozen_matrix_path()),
                "summary": actual_evidence["summary"],
                "reconciliation": file_evidence(planner.reconciliation_path()),
                "ledger": actual_evidence["ledger"],
                "reservations": actual_evidence["reservations"],
                "qualification": actual_evidence["qualification"],
                "controller_helper": file_evidence(self.helper),
                "source_bundle": file_evidence(planner.source_bundle_path()),
                "executable": file_evidence(planner.executable_path()),
                "build_manifest": build_bindings,
                "inputs": input_bindings,
                "manifests": manifest_evidence,
            },
            "profiles": list(self.profiles.values()),
            "allocation_authority": authority_evidence,
            "r03_authorization": (
                file_evidence(self.r03_authorization_path)
                if self.r03_authorization_path is not None
                else None
            ),
            "r17_readiness": {
                "decision": "approved" if self.r17_ready else "blocked",
                "reviewed_by": "fixture R17 reviewer",
                "reviewed_utc": utc_text(),
                "notes": "ready" if self.r17_ready else "predecessors remain",
                "evidence": r17_evidence,
            },
        }
        write_json(planner.profile_path(), profile)

    def plan(self) -> dict[str, object]:
        """Build one plan through the production path API."""

        return planner.plan_from_paths(
            planner.frozen_matrix_path(),
            planner.profile_path(),
            planner.summary_path(),
            planner.reconciliation_path(),
        )


@pytest.fixture
def campaign(tmp_path, monkeypatch) -> Campaign:
    """Return one default authenticated campaign with incomplete inactive R03."""

    root = tmp_path / "canonical"
    frozen = tmp_path / "frozen"
    helper = tmp_path / "repo/scripts/frontier/cgl_lf_stage_i.py"
    monkeypatch.setattr(planner, "CANONICAL_ROOT", root)
    monkeypatch.setattr(planner, "FROZEN_SOURCE", frozen)
    monkeypatch.setattr(planner, "CONTROLLER_HELPER", helper)
    monkeypatch.setattr(planner, "current_utc", lambda: NOW)
    value = Campaign(root, frozen, helper)
    monkeypatch.setattr(planner, "collect_live_scheduler_job", value.live_scheduler_job)
    monkeypatch.setattr(planner, "collect_live_batch_script", value.live_batch_script)
    monkeypatch.setattr(
        planner,
        "available_storage_bytes",
        lambda path: planner.R17_REQUIRED_RETENTION_BYTES * 2,
    )
    value.write_static_provenance()
    value.add_recorded("R02", segment(0, "0", "10"))
    r03 = value.add_recorded(
        "R03", segment(0, "0", "0.5"), final_time="0.31282347945569927",
        result="clean_partial",
    )
    value.add_cancelled("R03", planner.R03_F114_SEGMENT, parent_path=r03)
    value.authorize_r03(r03)
    value.refresh()
    return value


def rebind_profile_only(campaign: Campaign) -> None:
    """Refresh all expected bindings after an intentional semantic mutation."""

    campaign.refresh()


def republish_f115_fixture_chain(campaign: Campaign) -> None:
    """Rebind fixture review/audit bytes after an intentional F-115 semantic mutation."""

    authority_path = planner.r03_f115_path()
    assert authority_path.is_file()
    authority_path.chmod(0o444)
    planner.R03_F115_SHA256 = sha256(authority_path)
    published = {"path": str(authority_path), "sha256": planner.R03_F115_SHA256}
    review_specs = (
        (
            planner.r03_f115_provenance_security_review_path(),
            "R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256",
            "provenance_security",
        ),
        (
            planner.r03_f115_plasma_scientific_review_path(),
            "R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256",
            "plasma_scientific_continuation",
        ),
    )
    review_digests = {}
    for path, constant, key in review_specs:
        review = json.loads(path.read_text())
        review["published_f115"] = published
        review["reviewed_candidate"]["sha256"] = planner.R03_F115_SHA256
        write_json(path, review)
        path.chmod(0o444)
        digest = sha256(path)
        setattr(planner, constant, digest)
        review_digests[key] = digest

    authority = json.loads(authority_path.read_text())
    audit_path = planner.r03_f115_publication_audit_path()
    audit = json.loads(audit_path.read_text())
    audit["artifact"]["sha256"] = planner.R03_F115_SHA256
    audit["authority_and_enforcement"]["f115_authority"]["sha256"] = planner.R03_F115_SHA256
    audit["authority_and_enforcement"]["sole_next_segment_profile"] = authority[
        "authorization"
    ]["sole_next_segment_profile"]
    audit["independent_reviews"]["reviews_bind_exact_published_f115_sha256"] = (
        planner.R03_F115_SHA256
    )
    for key, digest in review_digests.items():
        audit["independent_reviews"][key]["sha256"] = digest
    write_json(audit_path, audit)
    audit_path.chmod(0o444)
    planner.R03_F115_PUBLICATION_AUDIT_SHA256 = sha256(audit_path)
    campaign.refresh()


def test_deployed_f115_anchors_and_reviewer_authorities_are_exact():
    deployed = load_utility()
    if not deployed.r03_f115_path().is_file():
        pytest.skip("deployed canonical F115 anchors are unavailable")
    authority, authority_evidence = deployed.read_json_file(
        deployed.r03_f115_path(), "deployed F115 authority"
    )
    deployed.require_immutable_publication_evidence(
        authority_evidence, deployed.R03_F115_SHA256, "deployed F115 authority"
    )
    assert authority["authorization"]["sole_next_segment_profile"]["segment"] == (
        deployed.R03_F115_SEGMENT
    )
    assert authority["implementation"]["source_bundle"]["verified_revisions"] == list(
        deployed.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS
    )
    audit, audit_evidence = deployed.read_json_file(
        deployed.r03_f115_publication_audit_path(), "deployed F115 publication audit"
    )
    deployed.require_immutable_publication_evidence(
        audit_evidence,
        deployed.R03_F115_PUBLICATION_AUDIT_SHA256,
        "deployed F115 publication audit",
    )
    assert audit["artifact"]["sha256"] == deployed.R03_F115_SHA256
    specs = (
        (
            deployed.r03_f115_provenance_security_review_path(),
            deployed.R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256,
            "provenance-security",
            "approved-for-publication",
        ),
        (
            deployed.r03_f115_plasma_scientific_review_path(),
            deployed.R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
            "plasma-scientific-continuation",
            "approved",
        ),
    )
    for path, digest, kind, decision in specs:
        review, evidence = deployed.read_json_file(path, f"deployed {kind} review")
        deployed.validate_f115_independent_review(
            review,
            evidence,
            path=path,
            digest=digest,
            kind=kind,
            decision=decision,
            label=f"deployed {kind} review",
        )


def test_deployed_f115_source_bundle_is_independently_usable_and_complete():
    deployed = load_utility()
    if not deployed.source_bundle_path().is_file():
        pytest.skip("deployed canonical F115 source bundle is unavailable")
    validation = deployed.validate_source_bundle_coverage()
    assert validation["git_bundle_verify"] == "passed"
    assert validation["self_contained"] is True
    assert validation["required_revisions"] == list(
        deployed.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS
    )


def test_initial_wave_includes_required_r03_and_exact_sole_profile(campaign):
    plan = campaign.plan()
    assert [packet["case_id"] for packet in plan["wave"]["packets"]] == [
        "R03", "R04", "R12", "R16"
    ]
    r03 = plan["wave"]["packets"][0]
    assert r03["lineage"]["segment"] == planner.R03_F115_SEGMENT
    assert r03["lineage"]["continuation_provenance"]["restart_time"] == "0.31282347945569927"
    assert r03["authorization_state"] == "planning_only_f115_profile_reference_non_authorizing"
    assert isinstance(r03["allocation"]["total_ranks"], int)
    assert "no prepare" in plan["authority"]
    assert "shared-root acknowledgement" in plan["authority"]
    assert plan["record_type"] == "cgl_lf_stage_i_read_only_wave_diagnostic"
    assert plan["disclosures"]


def test_continuation_packet_binds_restart_and_production_provenance(campaign):
    r04 = campaign.add_recorded("R04", segment(0, "0", "0.25"), nodes=1)
    campaign.refresh()
    plan = campaign.plan()
    packet = next(item for item in plan["wave"]["packets"] if item["case_id"] == "R04")
    assert packet["authorization_state"] == "planning_only_independent_profile_non_authorizing"
    assert packet["lineage"]["continuation_provenance"]["parent_manifest"]["path"] == str(r04)
    assert packet["production_provenance"]["controller_helper"]["sha256"] == planner.CONTROLLER_SHA256
    assert packet["production_provenance"]["source_bundle"]["sha256"] == planner.SOURCE_BUNDLE_SHA256
    assert packet["production_provenance"]["frozen_input"]["path"] == str(planner.input_path("R04"))
    assert packet["production_provenance"]["frozen_input"]["sha256"] == planner.EXPECTED_CASES["R04"][3]


def test_fabricated_local_profile_path_is_rejected(campaign, tmp_path):
    fake = tmp_path / "fake-profile.json"
    fake.write_bytes(planner.profile_path().read_bytes())
    fake.chmod(0o644)
    touch_now(fake)
    with pytest.raises(ValueError, match="exact canonical path"):
        planner.plan_from_paths(
            planner.frozen_matrix_path(), fake, planner.summary_path(),
            planner.reconciliation_path(),
        )


def test_stale_timestamp_evidence_is_rejected(campaign):
    campaign.profile_observed_utc = NOW - timedelta(days=1)
    campaign.profile_reviewed_utc = NOW - timedelta(days=1)
    campaign.refresh()
    with pytest.raises(ValueError, match="stale"):
        campaign.plan()


def test_stale_filesystem_mtime_with_current_timestamp_is_rejected(campaign):
    touch_now(planner.profile_path(), NOW - timedelta(days=1))
    with pytest.raises(ValueError, match="filesystem evidence is stale"):
        campaign.plan()


def test_profile_root_and_profile_set_ambiguity_are_rejected(campaign):
    profile = json.loads(planner.profile_path().read_text())
    profile["canonical_root"] = str(campaign.root / "fabricated")
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="identity/root is invalid"):
        campaign.plan()

    campaign.refresh()
    profile = json.loads(planner.profile_path().read_text())
    profile["profiles"].append(profile["profiles"][0])
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="unsupported or duplicated"):
        campaign.plan()

    campaign.refresh()
    profile = json.loads(planner.profile_path().read_text())
    profile["profiles"].pop()
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="cases differ"):
        campaign.plan()


def test_matrix_digest_semantics_and_input_records_are_bound(campaign):
    matrix = json.loads(planner.frozen_matrix_path().read_text())
    matrix["cases"][4]["input"] = matrix["cases"][2]["input"]
    write_json(planner.frozen_matrix_path(), matrix)
    planner.MATRIX_SHA256 = sha256(planner.frozen_matrix_path())
    campaign.refresh()
    with pytest.raises(ValueError, match="R06 semantics"):
        campaign.plan()

    matrix["cases"][4]["input"] = f"inputs/cgl_lf_paper/{planner.EXPECTED_CASES['R06'][1]}"
    write_json(planner.frozen_matrix_path(), matrix)
    planner.MATRIX_SHA256 = sha256(planner.frozen_matrix_path())
    write_text(planner.input_path("R06"), "tampered input\n")
    campaign.refresh()
    with pytest.raises(ValueError, match="R06 frozen input digest"):
        campaign.plan()


@pytest.mark.parametrize("kind", ["helper", "source", "executable", "build"])
def test_committed_helper_source_and_build_provenance_are_bound(campaign, kind):
    paths = {
        "helper": campaign.helper,
        "source": planner.source_bundle_path(),
        "executable": planner.executable_path(),
        "build": planner.build_manifest_path() / "environment.txt",
    }
    write_text(paths[kind], "tampered provenance\n", 0o755 if kind == "executable" else 0o644)
    campaign.refresh()
    with pytest.raises(ValueError, match="digest changed|does not bind"):
        campaign.plan()


def test_fabricated_summary_and_reconciliation_are_rejected(campaign):
    write_text(planner.summary_path(), campaign.summary_text() + "\n# fabricated\n")
    with pytest.raises(ValueError, match="summary evidence changed"):
        campaign.plan()
    campaign.refresh()
    snapshot = json.loads(planner.reconciliation_path().read_text())
    snapshot["generator"]["revision"] = "0" * 40
    write_json(planner.reconciliation_path(), snapshot)
    profile = json.loads(planner.profile_path().read_text())
    profile["evidence"]["reconciliation"] = file_evidence(planner.reconciliation_path())
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="generator is not the promoted controller"):
        campaign.plan()


def test_summary_must_match_canonical_ledger_and_reservations(campaign):
    summary = planner.summary_path().read_text().replace(
        "E03-forcing-policy Stage I actual use: `0.333334`",
        "E03-forcing-policy Stage I actual use: `9.999999`",
    )
    write_text(planner.summary_path(), summary)
    profile = json.loads(planner.profile_path().read_text())
    profile["evidence"]["summary"] = file_evidence(planner.summary_path())
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="summary rows differ|summary .* differs"):
        campaign.plan()


@pytest.mark.parametrize(
    ("case_id", "old_segment", "new_segment"),
    [
        ("R02", segment(0, "0", "10"), segment(99, "9", "10")),
        ("R03", segment(0, "0", "0.5"), segment(99, "0", "0.5")),
    ],
)
def test_r02_and_r03_lineages_require_contiguous_authenticated_parents(
    campaign, case_id, old_segment, new_segment
):
    old_path = campaign.manifest_path(case_id, old_segment)
    campaign.manifests[old_path]["run"]["segment"] = new_segment
    new_path = campaign.manifest_path(case_id, new_segment)
    campaign.manifests[new_path] = campaign.manifests.pop(old_path)
    old_path.unlink()
    reservation = next(item for item in campaign.reservations if item["case_id"] == case_id)
    row = next(item for item in campaign.ledger if item["case_id"] == case_id)
    reservation["manifest"] = str(new_path)
    reservation["segment"] = new_segment
    row["segment"] = new_segment
    campaign.manifests[new_path]["accounting"] = row
    if case_id == "R03":
        campaign.authorize_r03(new_path)
    campaign.refresh()
    with pytest.raises(ValueError, match="indexes are not contiguous|does not start fresh"):
        campaign.plan()


def test_parent_provenance_mismatch_is_rejected(campaign):
    parent = campaign.add_recorded("R04", segment(0, "0", "0.25"))
    child = campaign.add_recorded("R04", segment(1, "0.25", "0.5"), parent_path=parent)
    campaign.manifests[child]["command"]["parent_segment"]["restart_sha256"] = "0" * 64
    campaign.refresh()
    with pytest.raises(ValueError, match="parent provenance is invalid"):
        campaign.plan()


def test_historical_numeric_target_spelling_is_valid_but_ambiguity_is_rejected(campaign):
    path = campaign.manifest_path("R02", segment(0, "0", "10"))
    campaign.manifests[path]["command"]["overrides"] = ["time/tlim=10.0"]
    campaign.refresh()
    campaign.plan()

    campaign.manifests[path]["command"]["overrides"] = ["time/tlim=10", "time/tlim=10.0"]
    campaign.refresh()
    with pytest.raises(ValueError, match="time/tlim override is invalid"):
        campaign.plan()


def test_more_than_one_prepared_packet_is_rejected(campaign):
    campaign.add_active("R04", segment(0, "0", "0.25"), state="prepared")
    campaign.add_active("R12", segment(0, "0", "0.25"), state="prepared")
    campaign.refresh()
    with pytest.raises(ValueError, match="more than one prepared"):
        campaign.plan()


def test_one_prepared_packet_with_submitted_distinct_lanes_is_supported(campaign):
    campaign.add_active("R04", segment(0, "0", "0.25"), state="prepared")
    campaign.add_active("R12", segment(0, "0", "0.25"), state="submitted")
    campaign.refresh()
    plan = campaign.plan()
    assert plan["wave"]["active_lane_count"] == 2
    assert {item["case"] for item in plan["observed_state"]["active_lanes"]} == {"R04", "R12"}


@pytest.mark.parametrize(
    ("state", "exit_code"),
    [("FAILED", "0:0"), ("COMPLETED", "1:0")],
)
def test_recorded_scheduler_state_and_exit_code_must_be_valid(campaign, state, exit_code):
    row = campaign.ledger[0]
    row["state"] = state
    row["exit_code"] = exit_code
    path = campaign.manifest_path("R02", segment(0, "0", "10"))
    campaign.manifests[path]["accounting"] = row
    campaign.refresh()
    with pytest.raises(ValueError, match="scheduler/provenance fields are invalid"):
        campaign.plan()


def test_exact_legacy_scheduler_timestamp_form_is_supported(campaign):
    row = campaign.ledger[0]
    row["submitted_utc"] = "2026-06-03T00:00:00"
    row["completed_utc"] = "2026-06-03T00:10:00"
    path = campaign.manifest_path("R02", segment(0, "0", "10"))
    campaign.manifests[path]["accounting"] = row
    campaign.refresh()
    campaign.plan()

    row["submitted_utc"] = "2026-06-03 00:00:00"
    campaign.manifests[path]["accounting"] = row
    campaign.refresh()
    with pytest.raises(ValueError, match="ambiguous legacy scheduler timestamp"):
        campaign.plan()


@pytest.mark.parametrize(("field", "value"), [("nodes", "+1"), ("elapsed_seconds", "+600")])
def test_recorded_numeric_fields_match_controller_lexical_schema(campaign, field, value):
    row = campaign.ledger[0]
    row[field] = value
    path = campaign.manifest_path("R02", segment(0, "0", "10"))
    campaign.manifests[path]["accounting"] = row
    campaign.refresh()
    with pytest.raises(ValueError, match="nodes/elapsed fields differ"):
        campaign.plan()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("walltime", "03:00:00", "controller maximum"),
        ("athena_walltime", "01:59:00", "shutdown margin"),
        ("nodes", 8, "node policy"),
        ("nodes", "2", "must be an integer"),
    ],
)
def test_profile_resource_policy_is_exact(campaign, field, value, message):
    campaign.profiles["R04"][field] = value
    campaign.refresh()
    with pytest.raises(ValueError, match=message):
        campaign.plan()


def test_profile_rejects_unbounded_continuation_increment(campaign):
    campaign.profiles["R04"]["next_increment"] = "10"
    campaign.refresh()
    with pytest.raises(ValueError, match="exceeds the reviewed packet bound"):
        campaign.plan()


def test_active_resource_and_budget_fields_are_validated(campaign):
    campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    campaign.reservations[-1]["reserved_node_hours"] = 99
    campaign.refresh()
    with pytest.raises(ValueError, match="reservation differs from manifest"):
        campaign.plan()


@pytest.mark.parametrize("mutation", ["future_prepared", "bad_job", "unsupported_field"])
def test_active_reservation_scheduler_fields_are_validated(campaign, mutation):
    path = campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    reservation = campaign.reservations[-1]
    if mutation == "future_prepared":
        future = utc_text(NOW + timedelta(days=1))
        reservation["prepared_utc"] = future
        campaign.manifests[path]["prepared_utc"] = future
        message = "prepared_utc"
    elif mutation == "bad_job":
        reservation["job_id"] = "ambiguous"
        campaign.manifests[path]["job_id"] = "ambiguous"
        message = "job ID"
    else:
        reservation["scheduler_guess"] = "RUNNING"
        message = "unsupported controller fields"
    campaign.refresh()
    with pytest.raises(ValueError, match=message):
        campaign.plan()


def test_active_lane_must_match_fresh_reviewed_profile(campaign):
    campaign.add_active("R04", segment(0, "0", "0.25"), nodes=1)
    campaign.refresh()
    with pytest.raises(ValueError, match="differs from its fresh reviewed profile"):
        campaign.plan()


def test_active_r03_must_match_sole_authorized_profile(campaign):
    parent = campaign.manifest_path("R03", segment(0, "0", "0.5"))
    campaign.add_active(
        "R03", planner.R03_F115_SEGMENT, parent_path=parent,
        walltime="01:30:00", athena_walltime="01:20:00",
    )
    campaign.refresh()
    with pytest.raises(
        ValueError,
        match="resources/restart differ from the sole authorization|"
        "retained R03 F-115 s02 helper/source/resource profile differs",
    ):
        campaign.plan()


def test_duplicate_active_case_lanes_fail_closed(campaign):
    campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    campaign.add_active("R04", segment(1, "0.25", "0.5"), nodes=2)
    campaign.refresh()
    with pytest.raises(ValueError, match="duplicate active case lanes|duplicate manifest"):
        campaign.plan()


def test_four_lane_ceiling_fails_closed(campaign):
    for case_id in ("R04", "R05", "R06", "R07", "R08"):
        target = planner.decimal_text(planner.INITIAL_TARGETS[case_id])
        campaign.add_active(case_id, segment(0, "0", target))
    campaign.refresh()
    with pytest.raises(ValueError, match="four-lane ceiling"):
        campaign.plan()


def test_ten_node_ceiling_fails_closed(campaign):
    for case_id in ("R04", "R05", "R12"):
        campaign.profiles[case_id]["nodes"] = 4
        campaign.add_active(case_id, segment(0, "0", "0.25"), nodes=4)
    campaign.refresh()
    with pytest.raises(ValueError, match="ten-node ceiling"):
        campaign.plan()


def test_planned_wave_respects_ten_node_ceiling(campaign):
    campaign.profiles["R04"]["nodes"] = 4
    campaign.profiles["R12"]["nodes"] = 4
    campaign.refresh()
    plan = campaign.plan()
    assert plan["wave"]["total_nodes"] == planner.MAX_NODES
    assert {"case_id": "R16", "reason": "ten_node_ceiling"} in plan["wave"]["deferred"]


def test_planned_packets_are_charged_against_budget(campaign, monkeypatch):
    monkeypatch.setattr(
        planner, "STAGE_I_BUDGET_NODE_HOURS", campaign.cumulative + Decimal("1")
    )
    campaign.refresh()
    with pytest.raises(ValueError, match="credible remaining-campaign projection exceeds budget"):
        campaign.plan()


def test_r17_requires_bound_measured_reviewed_readiness(campaign):
    campaign.complete_predecessors()
    campaign.r17_ready = False
    campaign.refresh()
    with pytest.raises(ValueError, match="R17 prerequisites are not approved"):
        campaign.plan()
    campaign.r17_ready = True
    campaign.refresh()
    plan = campaign.plan()
    assert plan["wave"]["status"] == "r17_exclusive"
    assert plan["wave"]["total_nodes"] == 8


def test_r17_readiness_digest_measurement_and_review_fail_closed(campaign):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.refresh()
    storage_path = (
        campaign.root
        / "accounting/mks24_stage_i_E03_forcing_policy_R17_storage_readiness.json"
    )
    storage = json.loads(storage_path.read_text())
    storage["available_bytes"] += 1
    write_json(storage_path, storage)
    with pytest.raises(ValueError, match="evidence changed after review"):
        campaign.plan()

    campaign.refresh()
    storage = json.loads(storage_path.read_text())
    storage["available_bytes"] = 1
    write_json(storage_path, storage)
    profile = json.loads(planner.profile_path().read_text())
    profile["r17_readiness"]["evidence"]["storage"] = file_evidence(storage_path)
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="storage evidence does not satisfy"):
        campaign.plan()

    campaign.refresh()
    storage = json.loads(storage_path.read_text())
    storage["reviewed_utc"] = utc_text(NOW + timedelta(days=1))
    write_json(storage_path, storage)
    profile = json.loads(planner.profile_path().read_text())
    profile["r17_readiness"]["evidence"]["storage"] = file_evidence(storage_path)
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="stale or implausibly future"):
        campaign.plan()

    campaign.refresh()
    recost_path = (
        campaign.root
        / "accounting/mks24_stage_i_E03_forcing_policy_R17_node_hour_recost.json"
    )
    recost = json.loads(recost_path.read_text())
    recost["projection"]["computed_stage_i_total_node_hours"] = "0.000001"
    recost["projection_sha256"] = planner.sha256_bytes(
        planner.canonical_json(recost["projection"])
    )
    write_json(recost_path, recost)
    profile = json.loads(planner.profile_path().read_text())
    profile["r17_readiness"]["evidence"]["recost"] = file_evidence(recost_path)
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="differs from reproduced campaign projection"):
        campaign.plan()


def test_r17_last_and_exclusive_are_enforced(campaign):
    campaign.r17_ready = True
    campaign.add_active("R17", segment(0, "0", "0.25"), nodes=8)
    campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    campaign.refresh()
    with pytest.raises(ValueError, match="not exclusive"):
        campaign.plan()


def test_case_policy_classification_is_bound_to_exact_input_semantics(campaign):
    plan = campaign.plan()
    packets = {packet["case_id"]: packet for packet in plan["wave"]["packets"]}
    assert packets["R04"]["physics_policy_codes"] == ["U", "A", "H"]
    assert planner.policy_codes("R06") == ("U", "P", "H")
    assert planner.policy_codes("R14") == ("U", "A", "F")
    campaign.profiles["R04"]["case_id"] = "R06"
    campaign.refresh()
    with pytest.raises(ValueError, match="unsupported or duplicated|cases differ"):
        campaign.plan()


def test_missing_or_wrong_r03_sole_authorization_fails_closed(campaign):
    campaign.r03_authorization_path = None
    campaign.refresh()
    with pytest.raises(ValueError, match="incomplete R03 requires"):
        campaign.plan()
    campaign.authorize_r03(campaign.manifest_path("R03", segment(0, "0", "0.5")))
    authorization = json.loads(campaign.r03_authorization_path.read_text())
    authorization["authorization"]["sole_next_segment_profile"]["restart_time"] = 0.1
    write_json(campaign.r03_authorization_path, authorization)
    campaign.refresh()
    with pytest.raises(ValueError, match="exact immutable published artifact"):
        campaign.plan()


def test_r03_f115_requires_exact_path_digest_and_publication_audit(campaign):
    profile = json.loads(planner.profile_path().read_text())
    lookalike = campaign.root / "accounting/lookalike_F115.json"
    lookalike.write_bytes(campaign.r03_authorization_path.read_bytes())
    lookalike.chmod(0o644)
    touch_now(lookalike)
    profile["r03_authorization"] = file_evidence(lookalike)
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="exact published F-115 artifact"):
        campaign.plan()

    campaign.refresh()
    audit_path = planner.r03_f115_publication_audit_path()
    audit = json.loads(audit_path.read_text())
    audit["authority_and_enforcement"]["sole_next_segment_profile"]["time_tlim_target"] = 1.0
    write_json(audit_path, audit)
    with pytest.raises(ValueError, match="publication audit is not the exact immutable"):
        campaign.plan()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("segment", planner.R03_F114_SEGMENT),
        ("source_bundle", "source-archives/athenak-feature-cgl-through-c7e4fa30e.bundle"),
        ("cpus_per_task", 6),
        ("matrix", "/tmp/unreviewed-matrix.json"),
    ],
)
def test_r03_f115_exact_s02_resource_and_provenance_contract(campaign, field, value):
    authority = json.loads(campaign.r03_authorization_path.read_text())
    authority["authorization"]["sole_next_segment_profile"][field] = value
    write_json(campaign.r03_authorization_path, authority)
    republish_f115_fixture_chain(campaign)
    with pytest.raises(ValueError, match="exact authenticated parent/resources/provenance"):
        campaign.plan()


@pytest.mark.parametrize(
    "field",
    [
        "direct_sbatch_authorized",
        "reuse_cancelled_job_or_s01_authorized",
        "shared_root_acknowledgement_authorized_by_f115",
    ],
)
def test_r03_f115_publication_audit_cannot_over_authorize(campaign, field):
    audit_path = planner.r03_f115_publication_audit_path()
    audit = json.loads(audit_path.read_text())
    audit["authority_and_enforcement"][field] = True
    write_json(audit_path, audit)
    audit_path.chmod(0o444)
    planner.R03_F115_PUBLICATION_AUDIT_SHA256 = sha256(audit_path)
    with pytest.raises(ValueError, match="over-authorizes"):
        campaign.plan()


def test_r03_f115_requires_both_exact_immutable_independent_reviews(campaign):
    review_path = planner.r03_f115_plasma_scientific_review_path()
    review = json.loads(review_path.read_text())
    review["decision"] = "rejected"
    write_json(review_path, review)
    review_path.chmod(0o444)
    with pytest.raises(ValueError, match="exact immutable published artifact"):
        campaign.plan()


def test_r03_f115_reviewer_requires_explicit_campaign_authority(campaign):
    review_path = planner.r03_f115_provenance_security_review_path()
    review = json.loads(review_path.read_text())
    review["reviewer"] = {
        "agent_id": "self-asserted-reviewer",
        "identity": "nonempty but unauthorized reviewer",
    }
    write_json(review_path, review)
    republish_f115_fixture_chain(campaign)
    with pytest.raises(ValueError, match="lacks explicit F-115 authority"):
        campaign.plan()


def test_r03_f115_requires_retained_cancelled_s01_and_plans_only_s02(campaign):
    cancelled_path = campaign.manifest_path("R03", planner.R03_F114_SEGMENT)
    del campaign.manifests[cancelled_path]
    campaign.reservations = [
        item for item in campaign.reservations if item["manifest"] != str(cancelled_path)
    ]
    cancelled_path.unlink()
    campaign.refresh()
    with pytest.raises(ValueError, match="retained cancelled s01"):
        campaign.plan()


def test_r03_f115_cancelled_reservation_job_must_match_manifest_and_f115(campaign):
    reservation = next(
        item for item in campaign.reservations if item["segment"] == planner.R03_F114_SEGMENT
    )
    reservation["job_id"] = "999999"
    campaign.refresh()
    with pytest.raises(ValueError, match="cancelled reservation job/notes differ"):
        campaign.plan()


def test_r03_f115_cancelled_job_cannot_be_reused_by_active_manifest(campaign):
    cancelled_path = campaign.manifest_path("R03", planner.R03_F114_SEGMENT)
    cancelled_job = campaign.manifests[cancelled_path]["job_id"]
    active_path = campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    campaign.manifests[active_path]["job_id"] = cancelled_job
    campaign.reservations[-1]["job_id"] = cancelled_job
    campaign.refresh()
    with pytest.raises(
        ValueError,
        match="non-null job ID binds multiple operational manifest/reservation identities",
    ):
        campaign.plan()


def test_r03_f115_cancelled_job_cannot_be_reused_by_recorded_manifest(campaign):
    cancelled_path = campaign.manifest_path("R03", planner.R03_F114_SEGMENT)
    cancelled_job = campaign.manifests[cancelled_path]["job_id"]
    recorded_path = campaign.add_recorded("R04", segment(0, "0", "0.25"))
    campaign.manifests[recorded_path]["job_id"] = cancelled_job
    campaign.manifests[recorded_path]["accounting"]["job_id"] = cancelled_job
    campaign.reservations[-1]["job_id"] = cancelled_job
    campaign.refresh()
    with pytest.raises(
        ValueError,
        match="non-null job ID binds multiple operational manifest/reservation identities",
    ):
        campaign.plan()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("source", "retained cancelled s01 provenance/resources differ"),
        ("parent", "parent provenance is invalid"),
    ],
)
def test_r03_f115_cancelled_s01_parent_and_historical_source_are_exact(
    campaign, mutation, message
):
    cancelled_path = campaign.manifest_path("R03", planner.R03_F114_SEGMENT)
    command = campaign.manifests[cancelled_path]["command"]
    if mutation == "source":
        command["source_bundle"]["path"] = str(planner.source_bundle_path())
    else:
        command["parent_segment"]["restart_sha256"] = "0" * 64
    campaign.refresh()
    with pytest.raises(ValueError, match=message):
        campaign.plan()


@pytest.mark.parametrize("mutation", ["helper", "source_bundle", "resources"])
def test_recorded_r03_s02_cannot_downgrade_f115_profile(campaign, mutation):
    parent = campaign.manifest_path("R03", segment(0, "0", "0.5"))
    s02 = campaign.add_recorded(
        "R03", planner.R03_F115_SEGMENT, parent_path=parent,
        walltime="02:00:00", athena_walltime="01:50:00",
    )
    campaign.r03_authorization_path = None
    command = campaign.manifests[s02]["command"]
    if mutation == "helper":
        command["production_utility"]["sha256"] = "0" * 64
    elif mutation == "source_bundle":
        command["source_bundle"]["path"] = str(
            planner.expected_path(planner.R03_F114_SOURCE_BUNDLE)
        )
    else:
        campaign.manifests[s02]["allocation"]["cpus_per_task"] = 6
    campaign.refresh()
    with pytest.raises(
        ValueError,
        match="retained R03 F-115 s02 helper/source/resource profile differs|"
        "resource shape differs",
    ):
        campaign.plan()


def test_completed_r03_cannot_bypass_permanent_f115_anchor(campaign):
    campaign.complete_predecessors()
    review_path = planner.r03_f115_plasma_scientific_review_path()
    review = json.loads(review_path.read_text())
    review["reviewer"]["role"] = "self-asserted completed-case reviewer"
    write_json(review_path, review)
    republish_f115_fixture_chain(campaign)
    with pytest.raises(ValueError, match="lacks explicit F-115 authority"):
        campaign.plan()


@pytest.mark.parametrize("mutation", ["helper", "source_bundle", "resources"])
def test_completed_r03_still_rejects_recorded_s02_downgrade(campaign, mutation):
    campaign.complete_predecessors()
    s02 = campaign.manifest_path("R03", planner.R03_F115_SEGMENT)
    if mutation == "helper":
        campaign.manifests[s02]["command"]["production_utility"]["sha256"] = "0" * 64
    elif mutation == "source_bundle":
        campaign.manifests[s02]["command"]["source_bundle"]["sha256"] = (
            planner.R03_F114_SOURCE_BUNDLE_SHA256
        )
    else:
        campaign.manifests[s02]["command"]["athena_walltime"] = "01:40:00"
    campaign.refresh()
    with pytest.raises(ValueError, match="retained R03 F-115 s02 helper/source/resource"):
        campaign.plan()


def test_exact_active_r03_s02_is_valid_after_cancelled_s01(campaign):
    parent = campaign.manifest_path("R03", segment(0, "0", "0.5"))
    campaign.add_active("R03", planner.R03_F115_SEGMENT, parent_path=parent, nodes=1)
    campaign.refresh()
    plan = campaign.plan()
    assert "R03" not in [packet["case_id"] for packet in plan["wave"]["packets"]]
    assert plan["observed_state"]["active_lanes"][0]["segment"] == planner.R03_F115_SEGMENT


def test_planner_discloses_ct_and_provenance_limits_without_claiming_acceptance(campaign):
    plan = campaign.plan()
    serialized = json.dumps(plan)
    assert "maximum normalized CT divergence" not in serialized
    assert "does not authenticate a CT/divergence metric" in serialized
    assert "does not extend or exercise" in serialized
    assert "qualified executable and frozen scientific inputs remain separately bound" in serialized


def test_promoted_source_bundle_must_be_a_usable_git_bundle(campaign):
    bundle = planner.source_bundle_path()
    write_text(bundle, "not a git bundle\n")
    planner.SOURCE_BUNDLE_SHA256 = sha256(bundle)
    campaign.refresh()
    with pytest.raises(ValueError, match="not a Git bundle"):
        campaign.plan()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("duplicate", "inventory is not exact"),
        ("missing", "path component is missing"),
        ("bytes", "bytes differ from retained metadata"),
        ("unloadable", "loadable <par_end>"),
    ],
)
def test_terminal_restart_requires_unique_existing_digest_bound_loadable_files(
    campaign, mutation, message
):
    path = campaign.manifest_path("R03", segment(0, "0", "0.5"))
    terminal = campaign.manifests[path]["scientific_inspection"]["terminal_restart"]
    rank_files = terminal["rank_files"]
    target = Path(rank_files[1]["path"])
    if mutation == "duplicate":
        rank_files[1] = dict(rank_files[0])
        campaign.refresh()
    elif mutation == "missing":
        target.unlink()
    elif mutation == "bytes":
        write_bytes(target, b"changed restart bytes\n")
    else:
        write_bytes(target, b"not a loadable restart\n")
        rank_files[1]["sha256"] = sha256(target)
        rank_files[1]["size_bytes"] = target.stat().st_size
        campaign.refresh()
    with pytest.raises(ValueError, match=message):
        campaign.plan()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("missing_job", "no active job"),
        ("wrong_nodes", "differs from its submitted manifest"),
        ("wrong_owner", "differs from its submitted manifest"),
        ("batch_script", "Slurm-stored batch script differs"),
    ],
)
def test_submitted_jobs_require_live_slurm_and_stored_script_authentication(
    campaign, mutation, message
):
    path = campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    job_id = str(campaign.manifests[path]["job_id"])
    campaign.refresh()
    if mutation == "missing_job":
        del campaign.live_jobs[job_id]
    elif mutation == "wrong_nodes":
        campaign.live_jobs[job_id]["nodes"] = 1
    elif mutation == "wrong_owner":
        campaign.live_jobs[job_id]["owner"] = "fabricated-owner"
    else:
        campaign.live_scripts[job_id] += b"# fabricated\n"
    with pytest.raises(ValueError, match=message):
        campaign.plan()


def test_profile_cannot_self_assert_resources_without_independent_authority(campaign):
    profile = json.loads(planner.profile_path().read_text())
    r04 = next(item for item in profile["profiles"] if item["case_id"] == "R04")
    r04["nodes"] = 4
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="independent published authority"):
        campaign.plan()

    campaign.refresh()
    audit_path = planner.allocation_authority_audit_path()
    audit = json.loads(audit_path.read_text())
    audit["review"]["reviewed_by"] = "fixture wave reviewer"
    write_json(audit_path, audit)
    profile = json.loads(planner.profile_path().read_text())
    profile["allocation_authority"]["publication_audit"] = file_evidence(audit_path)
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="lacks independent approval"):
        campaign.plan()


def test_independent_budget_reproduction_rejects_self_asserted_tiny_projection(campaign):
    authority_path = planner.allocation_authority_path()
    authority = json.loads(authority_path.read_text())
    authority["projection"]["computed_stage_i_total_node_hours"] = "0.000001"
    write_json(authority_path, authority)
    audit_path = planner.allocation_authority_audit_path()
    audit = json.loads(audit_path.read_text())
    audit["artifact"]["sha256"] = sha256(authority_path)
    write_json(audit_path, audit)
    profile = json.loads(planner.profile_path().read_text())
    profile["allocation_authority"] = {
        "artifact": file_evidence(authority_path),
        "publication_audit": file_evidence(audit_path),
    }
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match="projection differs from authenticated campaign state"):
        campaign.plan()


@pytest.mark.parametrize("store_index", [0, 1])
def test_real_transaction_discovery_rejects_stage_i_and_recost_journals(campaign, store_index):
    transaction = planner.transaction_store_paths()[store_index] / "pending.json"
    write_json(transaction, {"pending": True})
    with pytest.raises(ValueError, match="transaction stores are not empty"):
        campaign.plan()


def test_global_toctou_boundary_rejects_late_canonical_mutation(campaign, monkeypatch):
    original = planner.revalidate_global_boundary

    def mutate_then_revalidate(value, transaction_state, submitted_scheduler):
        write_json(planner.reservations_path(), [])
        return original(value, transaction_state, submitted_scheduler)

    monkeypatch.setattr(planner, "revalidate_global_boundary", mutate_then_revalidate)
    with pytest.raises(ValueError, match="global evidence changed before plan completion"):
        campaign.plan()


def test_global_toctou_boundary_rechecks_live_slurm(campaign, monkeypatch):
    campaign.add_active("R04", segment(0, "0", "0.25"), nodes=2)
    campaign.refresh()
    calls = 0

    def changing_live_job(job_id):
        nonlocal calls
        calls += 1
        value = campaign.live_scheduler_job(job_id)
        if calls > 1:
            value["state"] = "RUNNING"
        return value

    monkeypatch.setattr(planner, "collect_live_scheduler_job", changing_live_job)
    with pytest.raises(ValueError, match="live Slurm state changed before plan completion"):
        campaign.plan()


def test_r17_storage_readiness_requires_live_capacity(campaign, monkeypatch):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.refresh()
    monkeypatch.setattr(planner, "available_storage_bytes", lambda path: 1)
    with pytest.raises(ValueError, match="storage evidence does not satisfy"):
        campaign.plan()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("scheduler", "successful measured run"),
        ("output", "output inventory is not exact per rank"),
        ("restart", "loadable <par_end>"),
        ("performance", "performance differs from measured evidence"),
    ],
)
def test_r17_rank_readiness_requires_scheduler_output_restart_and_performance(
    campaign, mutation, message
):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.refresh()
    rank_path = (
        campaign.root
        / "accounting/mks24_stage_i_E03_forcing_policy_R17_64_rank_readiness.json"
    )
    rank = json.loads(rank_path.read_text())
    if mutation == "scheduler":
        scheduler_path = Path(rank["scheduler_evidence"]["path"])
        write_text(
            scheduler_path,
            "900001|fabricated_job|COMPLETED|0:0|8|600|"
            f"{utc_text()}|{utc_text()}\n",
        )
        rank["scheduler_evidence"] = file_evidence(scheduler_path)
    elif mutation == "output":
        rank["output_inventory"][1] = dict(rank["output_inventory"][0])
    elif mutation == "restart":
        item = rank["terminal_restart"]["rank_files"][1]
        restart_path = Path(item["path"])
        write_bytes(restart_path, b"not loadable\n")
        item["sha256"] = sha256(restart_path)
        item["size_bytes"] = restart_path.stat().st_size
    else:
        rank["performance"]["elapsed_seconds"] = 1
    write_json(rank_path, rank)
    profile = json.loads(planner.profile_path().read_text())
    profile["r17_readiness"]["evidence"]["rank_readiness"] = file_evidence(rank_path)
    write_json(planner.profile_path(), profile)
    with pytest.raises(ValueError, match=message):
        campaign.plan()


def test_leaf_and_intermediate_symlink_evidence_is_rejected(campaign, tmp_path):
    link = tmp_path / "summary-link"
    link.symlink_to(planner.summary_path())
    with pytest.raises(ValueError, match="symlink component"):
        planner.read_stable_regular_file(link, "linked summary")

    real = tmp_path / "real"
    write_text(real / "value.json", "{}\n")
    directory_link = tmp_path / "directory-link"
    directory_link.symlink_to(real, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink component"):
        planner.read_stable_regular_file(directory_link / "value.json", "linked directory")


def test_mutation_during_stable_read_is_rejected(tmp_path, monkeypatch):
    path = tmp_path / "changing-evidence.bin"
    write_bytes(path, b"a" * (2 * 1024 * 1024))
    original_read = planner.os.read
    changed = False

    def mutating_read(descriptor, size):
        nonlocal changed
        block = original_read(descriptor, size)
        if block and not changed:
            changed = True
            path.write_bytes(b"fabricated replacement")
        return block

    monkeypatch.setattr(planner.os, "read", mutating_read)
    with pytest.raises(ValueError, match="changed while it was read"):
        planner.read_stable_regular_file(path, "changing evidence")


def test_cli_paths_are_exact_and_output_is_command_free(campaign, capsys):
    result = planner.main(
        [
            "--matrix", str(planner.frozen_matrix_path()),
            "--allocation-profile", str(planner.profile_path()),
            "--summary", str(planner.summary_path()),
            "--reconciliation", str(planner.reconciliation_path()),
        ]
    )
    captured = capsys.readouterr()
    assert result == 0
    assert captured.err == ""
    plan = json.loads(captured.out)
    serialized = json.dumps(plan)
    assert plan["read_only"] is True
    assert "sbatch " not in serialized
    assert "srun " not in serialized
    assert plan["wave"]["packets"][0]["authorization_evidence"]["planner_scope"][
        "direct_sbatch_authorized"
    ] is False


def test_plan_digest_is_deterministic(campaign):
    first = campaign.plan()
    second = campaign.plan()
    assert first == second
    core = {key: value for key, value in first.items() if key != "plan_sha256"}
    assert first["plan_sha256"] == planner.sha256_bytes(planner.canonical_json(core))
