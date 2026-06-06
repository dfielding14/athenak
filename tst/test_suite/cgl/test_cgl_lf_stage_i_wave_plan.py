"""Focused and adversarial tests for the authenticated Stage I wave planner."""

from __future__ import annotations

import csv
from copy import deepcopy
from datetime import datetime, timedelta, timezone
from decimal import Decimal
import ast
import hashlib
import importlib.util
import inspect
import io
import json
import os
from pathlib import Path
import struct
import subprocess
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_wave_plan.py"
CONTROLLER_UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i.py"
CHECKPOINT_UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_checkpoint.py"
QUALIFICATION_TEST_UTILITY = REPOSITORY / "tst/test_suite/cgl/test_cgl_lf_stage_i_qualification.py"
RECOST_UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_recost.py"
SOURCE_AUTHORITY_UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_source_authority.py"
NOW = datetime.now(timezone.utc).replace(microsecond=0)


def load_utility(path=UTILITY, name="cgl_lf_stage_i_wave_plan"):
    """Load one Stage I utility without executing its CLI."""

    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


planner = load_utility()
controller = load_utility(CONTROLLER_UTILITY, "cgl_lf_stage_i_controller_for_wave_tests")


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
    if path.exists():
        path.chmod(0o644)
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


def retained_file(path: Path) -> dict[str, object]:
    """Return the exact retained-file record emitted by the controller."""

    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def write_history(path: Path, columns: dict[str, list[float]]) -> dict[str, list[float]]:
    """Write one finite controller-compatible Athena history."""

    labels = list(columns)
    rows = list(zip(*(columns[label] for label in labels)))
    text = "# " + " ".join(
        f"[{index}]={label}" for index, label in enumerate(labels, start=1)
    ) + "\n"
    text += "\n".join(" ".join(format(value, ".17g") for value in row) for row in rows)
    write_text(path, text + "\n")
    return columns


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
        self.r17_source_authority_override: dict[str, object] | None = None
        self.r17_output_count = 64
        self.r17_duplicate_output_rank = False
        self.r17_decomposition_mutation: str | None = None
        self.r17_account_mutation: str | None = None
        self.r17_physics_mutation: str | None = None
        self.recost_reviewer_agent = "fixture-independent-recost-reviewer"
        self.request_reviewer_agent = "fixture-independent-request-reviewer"
        self.budget_mutator = None
        self.recommended_cases = ["R03", "R04", "R12", "R16"]
        self.r03_next_increment = Decimal("0.25")
        self.recost_path = (
            root / "accounting/mks24_stage_i_E03_forcing_policy_F200_recost_evidence.json"
        )
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
                "nodes": 8 if case_id == "R17" else 4 if case_id in {"R04", "R12"} else 1,
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

    def write_product_group(
        self, root: Path, count: int, name: str, payload: bytes
    ) -> dict[str, object]:
        """Write one exact controller-style rank-local output product."""

        rank_files = []
        for rank in range(count):
            path = root / f"rank_{rank:08d}" / name
            write_bytes(path, payload)
            rank_files.append(retained_file(path))
        return {
            **rank_files[0],
            "storage": "per_rank",
            "rank_files": rank_files,
        }

    def write_continuation_histories(
        self, case_id: str, output: Path, start: Decimal, final: Decimal
    ) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
        """Retain histories and controller-emitted plasma continuation evidence."""

        times = [float(start), float(final)]
        passive = case_id in planner.PASSIVE_CASES
        mhd = {
            "time": times,
            "mass": [1.0, 1.0],
            "tot-E": [10.0, 11.0],
            "lf_nstage": [0.0, 1.0],
            "lf_qface": [0.0, 1.0],
            "lf_qprcap": [0.0, 0.0],
            "lf_qpr10": [0.0, 0.0],
            "lf_qpecap": [0.0, 0.0],
            "lf_qpe10": [0.0, 0.0],
            "lf_qprwrk": [0.0, 0.0],
            "lf_qpewrk": [0.0, 0.0],
            "lf_cpwrk": [0.0, 0.0 if passive else 0.001],
            "lf_cawrk": [0.0, 0.0],
            "lf_hwproj": [0.0, 0.0],
            **{name: [0.0, 0.0] for name in planner.STRICT_LF_FAILURE_COLUMNS},
        }
        user = {
            "time": times,
            "mass": [1.0, 1.0],
            "hard_vol": [0.0, 0.0],
            "force_work": [0.0, 1.0],
            "max_ndiv": [0.0, 1.0e-13],
        }
        mhd_path = output / f"{case_id}.mhd.hst"
        user_path = output / f"{case_id}.user.hst"
        write_history(mhd_path, mhd)
        write_history(user_path, user)
        return (
            retained_file(mhd_path),
            retained_file(user_path),
            controller.continuation_plasma_evidence(case_id, mhd, user),
        )

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
                "revision": self.controller_revision,
                "sha256": self.controller_sha256,
            },
            "source_bundle": {
                "path": str(self.current_bundle),
                "sha256": self.current_bundle_sha256,
                "verified_revisions": self.current_bundle_revisions,
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
        if case_id == "R03" and segment_id == planner.R03_F115_SEGMENT:
            command["production_utility"] = {
                "committed": True,
                "path": str(self.helper),
                "revision": planner.R03_F115_CONTROLLER_REVISION,
                "sha256": planner.R03_F115_CONTROLLER_SHA256,
            }
            command["source_bundle"] = {
                "path": str(planner.f115_source_bundle_path()),
                "sha256": planner.R03_F115_SOURCE_BUNDLE_SHA256,
                "verified_revisions": [
                    planner.SOURCE_REVISION,
                    planner.R03_F115_CONTROLLER_REVISION,
                ],
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
        elapsed_seconds: int = 600,
    ) -> Path:
        """Retain one recorded fixture segment."""

        path = self.manifest_path(case_id, segment_id)
        manifest = self.base_manifest(
            case_id, segment_id, "recorded", nodes, walltime, athena_walltime, parent_path
        )
        _, start, target = planner.parse_segment(segment_id, "fixture recorded segment")
        final = target if final_time is None else Decimal(final_time)
        elapsed = elapsed_seconds
        assert elapsed > 0
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
        accepted = result == "accepted"
        output = path.parents[1] / "output"
        ranks = nodes * planner.RANKS_PER_NODE
        mhd_history, user_history, plasma_evidence = self.write_continuation_histories(
            case_id, output, start, final
        )
        snapshot = self.write_product_group(
            output / "bin",
            ranks,
            f"{case_id}_{segment_id}.00001.bin",
            b"rank-local-snapshot\n",
        )
        terminal = self.terminal_restart(case_id, segment_id, nodes, final)
        inspection = {
            "schema_version": 4,
            "execution_epoch": planner.EXECUTION_EPOCH,
            "inspected_utc": utc_text(),
            "manifest": str(path),
            "job_id": job,
            "case_id": case_id,
            "segment": segment_id,
            "required_time": float(target),
            "final_time": float(final),
            "maximum_strict_failure_counts": {
                key: 0.0 for key in planner.STRICT_LF_FAILURE_COLUMNS
            },
            "checks": {
                "required_time_reached": accepted,
                "strict_lf_failure_counters_zero": True,
                "snapshots_retained": True,
                "terminal_snapshot_retained": True,
                "restart_retained": True,
                "terminal_restart_physical_time_matches_final": True,
                "plasma_continuation_policy": True,
            },
            "accepted": accepted,
            "clean_for_continuation": result in {"accepted", "clean_partial"},
            "mhd_history": mhd_history,
            "user_history": user_history,
            "plasma_continuation_policy": controller.CONTINUATION_PLASMA_POLICY,
            "plasma_continuation_evidence": plasma_evidence,
            "snapshots": [snapshot],
            "snapshot_times": [float(final)],
            "restarts": [terminal],
            "restart_times": [float(final)],
            "restart_time_marker_modes": [
                ["full_precision"] * ranks,
            ],
            "terminal_restart_time": float(final),
            "terminal_restart": terminal,
            "restart_time_marker_bypass": False,
            "final_hardwall_projection_count": 0.0,
        }
        manifest["scientific_inspection"] = inspection
        write_json(path.parent / "segment_inspection.json", inspection)
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
            "source_bundle": str(planner.f115_source_bundle_path()),
            "source_bundle_sha256": planner.R03_F115_SOURCE_BUNDLE_SHA256,
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
                "to": str(planner.f115_source_bundle_path()),
            },
            "source_bundle_sha256": {
                "from": planner.R03_F114_SOURCE_BUNDLE_SHA256,
                "to": planner.R03_F115_SOURCE_BUNDLE_SHA256,
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
                        "head": planner.R03_F115_CONTROLLER_REVISION,
                        "links": 1,
                        "mode": "0644",
                        "path": planner.f115_source_bundle_path().relative_to(
                            planner.CANONICAL_ROOT
                        ).as_posix(),
                        "sha256": planner.R03_F115_SOURCE_BUNDLE_SHA256,
                        "verified_revisions": list(
                            planner.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS
                        ),
                    },
                    "stage_i_helper": {
                        "links": 1,
                        "mode": "0644",
                        "path": str(self.helper),
                        "sha256": planner.R03_F115_CONTROLLER_SHA256,
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
                        "head": planner.R03_F115_CONTROLLER_REVISION,
                        "links": 1,
                        "mode": "0644",
                        "path": str(planner.f115_source_bundle_path()),
                        "sha256": planner.R03_F115_SOURCE_BUNDLE_SHA256,
                    },
                    "committed_stage_i_helper": {
                        "path": "scripts/frontier/cgl_lf_stage_i.py",
                        "sha256": planner.R03_F115_CONTROLLER_SHA256,
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
        for case_id in (case for case in planner.LOWER_CASES if case != "R12"):
            first = planner.decimal_text(planner.INITIAL_TARGETS[case_id])
            parent = self.add_recorded(case_id, segment(0, "0", first), nodes=1)
            self.add_recorded(case_id, segment(1, first, "10"), parent_path=parent, nodes=1)
        next_unallocated_job = max(
            (int(manifest["job_id"]) for manifest in self.manifests.values()), default=0
        ) + 1
        self.next_job = int(planner.R12_HISTORICAL_PARTIAL_JOB_ID)
        self.add_recorded(
            "R12",
            planner.R12_HISTORICAL_PARTIAL_SEGMENT,
            final_time="0.1371931229426507",
            result="clean_partial",
            nodes=planner.R12_FRESH_RERUN_NODES,
            elapsed_seconds=6657,
        )
        self.next_job = max(self.next_job, next_unallocated_job)
        fresh = self.add_recorded(
            "R12",
            planner.R12_FRESH_RERUN_SEGMENT,
            nodes=planner.R12_FRESH_RERUN_NODES,
            walltime=planner.R12_FRESH_RERUN_WALLTIME,
            athena_walltime=planner.R12_FRESH_RERUN_ATHENA_WALLTIME,
        )
        self.add_recorded(
            "R12",
            segment(2, "0.12", "10"),
            parent_path=fresh,
            nodes=planner.R12_FRESH_RERUN_NODES,
        )

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
        planner.R03_F115_CONTROLLER_REVISION = git(repository, "rev-parse", "HEAD")
        planner.R03_F115_CONTROLLER_SHA256 = sha256(self.helper)
        planner.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS = (
            planner.SOURCE_REVISION,
            planner.R03_F114_CONTROLLER_REVISION,
            planner.R03_F115_CONTROLLER_REVISION,
        )
        f115_bundle = planner.f115_source_bundle_path()
        f115_bundle.parent.mkdir(parents=True, exist_ok=True)
        git(repository, "bundle", "create", str(f115_bundle), "HEAD")
        f115_bundle.chmod(0o644)
        touch_now(f115_bundle)
        planner.R03_F115_SOURCE_BUNDLE_SHA256 = sha256(f115_bundle)
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
        repository_matrix = repository / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
        write_json(repository_matrix, matrix)
        write_text(self.helper, "fixture F116 retained bridge controller\n")
        git(repository, "add", "scripts/frontier/cgl_lf_stage_i.py")
        git(repository, "commit", "-m", "fixture F116 retained bridge")
        planner.F116_BRIDGE_REVISION = git(repository, "rev-parse", "HEAD")
        planner.F116_BRIDGE_SOURCE_BUNDLE = (
            f"source-archives/athenak-feature-cgl-through-{planner.F116_BRIDGE_REVISION[:9]}.bundle"
        )
        bridge_bundle = planner.expected_path(planner.F116_BRIDGE_SOURCE_BUNDLE)
        git(repository, "branch", "feature/cgl-landau-fluid")
        git(repository, "bundle", "create", str(bridge_bundle), "feature/cgl-landau-fluid")
        bridge_bundle.chmod(0o644)
        touch_now(bridge_bundle)
        planner.F116_BRIDGE_SHA256 = sha256(bridge_bundle)

        for relative, mode in planner.F116_REQUIRED_TOOLS.items():
            path = repository / relative
            write_text(path, f"fixture promoted F116 tool {relative}\n", int(mode, 8))
        git(repository, "add", "inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
        for relative in planner.F116_REQUIRED_TOOLS:
            git(repository, "add", relative)
        git(repository, "commit", "-m", "fixture current promoted recost/controller bundle")
        git(repository, "branch", "-f", "feature/cgl-landau-fluid", "HEAD")
        self.controller_revision = git(repository, "rev-parse", "HEAD")
        self.controller_sha256 = sha256(self.helper)
        self.generator_sha256 = sha256(
            repository / "scripts/frontier/cgl_lf_stage_i_recost.py"
        )
        self.committed_tools = [
            {
                "path": relative,
                "revision": self.controller_revision,
                "sha256": sha256(repository / relative),
                "mode": mode,
            }
            for relative, mode in sorted(planner.F116_REQUIRED_TOOLS.items())
        ]
        self.current_bundle = planner.promoted_source_bundle_path(self.controller_revision)
        git(repository, "bundle", "create", str(self.current_bundle), "feature/cgl-landau-fluid")
        self.current_bundle.chmod(0o644)
        touch_now(self.current_bundle)
        self.current_bundle_sha256 = sha256(self.current_bundle)
        self.current_bundle_revisions = [
            planner.SOURCE_REVISION,
            planner.R03_F114_CONTROLLER_REVISION,
            planner.R03_F115_CONTROLLER_REVISION,
            planner.F116_BRIDGE_REVISION,
            self.controller_revision,
        ]
        planner.F116_PRODUCTION_REQUIRED_REVISIONS = tuple(self.current_bundle_revisions)
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

    def publish_f116_current_source_authority(self) -> dict[str, object]:
        """Publish the production-schema F116 current-source authority chain."""

        authority_path = planner.f116_current_source_authority_path()
        relative_bundle = self.current_bundle.relative_to(self.root).as_posix()
        bridge_bundle = planner.expected_path(planner.F116_BRIDGE_SOURCE_BUNDLE)
        relative_bridge = bridge_bundle.relative_to(self.root).as_posix()
        readme_path = self.root / "source-archives/README.md"
        sums_path = self.root / "source-archives/SHA256SUMS"
        write_text(
            readme_path,
            "## AthenaK\n\nfixture production F116 current-source authority catalog\n",
        )
        write_text(
            sums_path,
            (
                f"{planner.R03_F115_SOURCE_BUNDLE_SHA256}  "
                f"{Path(planner.R03_F115_SOURCE_BUNDLE).name}\n"
                f"{planner.F116_BRIDGE_SHA256}  {bridge_bundle.name}\n"
                f"{self.current_bundle_sha256}  {self.current_bundle.name}\n"
            ),
        )
        historical = {
            "evidence": {
                "path": planner.r03_f115_path().relative_to(self.root).as_posix(),
                "sha256": planner.R03_F115_SHA256,
            },
            "publication_audit": {
                "path": planner.r03_f115_publication_audit_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": planner.R03_F115_PUBLICATION_AUDIT_SHA256,
            },
            "provenance_review": {
                "path": planner.r03_f115_provenance_security_review_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": planner.R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256,
            },
            "plasma_review": {
                "path": planner.r03_f115_plasma_scientific_review_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": planner.R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
            },
        }
        bridge = {
            "path": relative_bridge,
            "sha256": planner.F116_BRIDGE_SHA256,
            "complete_history": True,
            "head": planner.F116_BRIDGE_REVISION,
            "advertised_tip": {
                "revision": planner.F116_BRIDGE_REVISION,
                "name": "refs/heads/feature/cgl-landau-fluid",
            },
            "verified_revisions": [
                planner.SOURCE_REVISION,
                planner.R03_F114_CONTROLLER_REVISION,
                planner.R03_F115_CONTROLLER_REVISION,
                planner.F116_BRIDGE_REVISION,
            ],
            "selected_as_current": False,
            "role": "retained-non-current-bridge",
        }
        final = {
            "path": relative_bundle,
            "sha256": self.current_bundle_sha256,
            "complete_history": True,
            "head": self.controller_revision,
            "advertised_tip": {
                "revision": self.controller_revision,
                "name": "refs/heads/feature/cgl-landau-fluid",
            },
            "verified_revisions": self.current_bundle_revisions,
            "selected_as_current": True,
            "candidate_path": f"/tmp/{self.current_bundle.name}.candidate",
            "subject": "fixture promoted F116 production source authority",
        }
        after_catalog = {
            "readme_sha256": sha256(readme_path),
            "sha256sums_sha256": sha256(sums_path),
            "bridge_listed_exactly_once": True,
            "final_bundle_listed_exactly_once": True,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
            "sole_current_source_bundle": relative_bundle,
        }
        authority = {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-evidence",
            "checkpoint": "F-116",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "generated_utc": utc_text(),
            "scope": {
                "relationship": "current-source-selection-only-supersession",
                "summary": "Select the exact committed final tooling source without execution authority.",
                "preserves": planner.F116_SCOPE_PRESERVES,
                "does_not_authorize": planner.F116_SCOPE_DOES_NOT_AUTHORIZE,
            },
            "predecessor_authorities": {"historical_f115": historical},
            "implementation": {
                "publisher": next(
                    item
                    for item in self.committed_tools
                    if item["path"] == planner.F116_PUBLISHER_RELATIVE
                ),
                "committed_tools": self.committed_tools,
                "intermediate_36140_bundle": bridge,
                "current_source_bundle": final,
            },
            "source_archive_catalog": {
                "before": {
                    "readme_sha256": "0" * 64,
                    "sha256sums_sha256": "1" * 64,
                    "bridge_listed": False,
                    "final_bundle_listed": False,
                    "corrupt_c7_listed": False,
                },
                "after": after_catalog,
            },
            "authorization": planner.F116_AUTHORIZATION,
            "validation": planner.F116_VALIDATION_CLAIMS,
            "publication_requirements": planner.F116_PUBLICATION_REQUIREMENTS,
        }
        write_json(authority_path, authority)
        authority_path.chmod(0o444)
        verified = {
            "authorization_broadening": False,
            "bridge_selected_as_current": False,
            "corrupt_c7_excluded": True,
            "current_source_selection_only": True,
            "final_bundle_sha256": self.current_bundle_sha256,
            "final_head": self.controller_revision,
            "historical_f115_preserved": True,
        }
        review_specs = (
            (
                planner.f116_provenance_security_review_path(),
                "provenance-security",
                "approved-for-publication",
                "fixture-f116-provenance-reviewer",
                "fixture independent F116 provenance/security reviewer",
            ),
            (
                planner.f116_plasma_scientific_review_path(),
                "plasma-scientific-continuation",
                "approved",
                "fixture-f116-plasma-reviewer",
                "fixture independent F116 plasma/scientific reviewer",
            ),
        )
        for path, kind, decision, agent_id, identity in review_specs:
            write_json(
                path,
                {
                    "schema_version": 1,
                    "record_type": (
                        "stage-i-current-source-authority-supersession-independent-review"
                    ),
                    "checkpoint": "F-116",
                    "execution_epoch": planner.EXECUTION_EPOCH,
                    "review_kind": kind,
                    "decision": decision,
                    "reviewed_candidate": {
                        "path": f"/tmp/{authority_path.name}.candidate",
                        "sha256": sha256(authority_path),
                    },
                    "published_f116": {
                        "path": str(authority_path),
                        "sha256": sha256(authority_path),
                    },
                    "reviewer": {"agent_id": agent_id, "identity": identity},
                    "reviewed_utc": utc_text(),
                    "findings": ["Exact production F116 source-only authority verified."],
                    "limitations": ["No prepare, submit, scheduler, or scientific authority."],
                    "verified": verified,
                },
            )
            path.chmod(0o444)
        audit_path = planner.f116_publication_audit_path()
        write_json(
            audit_path,
            {
                "schema_version": 1,
                "record_type": (
                    "stage-i-current-source-authority-supersession-publication-audit"
                ),
                "checkpoint": "F-116",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "published_utc": utc_text(),
                "artifact": {
                    "path": str(authority_path),
                    "sha256": sha256(authority_path),
                    "mode": "0444",
                    "links": 1,
                },
                "independent_reviews": {
                    "reviews_bind_exact_published_f116_sha256": sha256(authority_path),
                    "provenance_security": {
                        "path": str(planner.f116_provenance_security_review_path()),
                        "sha256": sha256(planner.f116_provenance_security_review_path()),
                        "mode": "0444",
                        "links": 1,
                    },
                    "plasma_scientific_continuation": {
                        "path": str(planner.f116_plasma_scientific_review_path()),
                        "sha256": sha256(planner.f116_plasma_scientific_review_path()),
                        "mode": "0444",
                        "links": 1,
                    },
                },
                "historical_f115_authority": {
                    "evidence_sha256": planner.R03_F115_SHA256,
                    "publication_audit_sha256": planner.R03_F115_PUBLICATION_AUDIT_SHA256,
                    "provenance_review_sha256": (
                        planner.R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256
                    ),
                    "plasma_review_sha256": planner.R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
                },
                "source_archive_catalog": {
                    "readme": {
                        "path": str(readme_path),
                        "sha256": sha256(readme_path),
                        "mode": "0644",
                        "links": 1,
                    },
                    "sha256sums": {
                        "path": str(sums_path),
                        "sha256": sha256(sums_path),
                        "mode": "0644",
                        "links": 1,
                    },
                    "bridge_bundle": {
                        "path": str(bridge_bundle),
                        "sha256": planner.F116_BRIDGE_SHA256,
                        "mode": "0644",
                        "links": 1,
                        "head": planner.F116_BRIDGE_REVISION,
                        "role": "retained-non-current-bridge",
                        "selected_as_current": False,
                    },
                    "current_source_bundle": {
                        "path": str(self.current_bundle),
                        "sha256": self.current_bundle_sha256,
                        "mode": "0644",
                        "links": 1,
                        "head": self.controller_revision,
                        "selected_as_current": True,
                    },
                    "corrupt_c7_absent_from_active_checksum_ledger": True,
                    "sole_current_source_bundle": str(self.current_bundle),
                },
                "authority_and_enforcement": planner.F116_AUTHORIZATION,
                "publication": planner.F116_PUBLICATION_METHOD,
            },
        )
        audit_path.chmod(0o444)
        return {
            "checkpoint": "F-116",
            "evidence": {
                "path": authority_path.relative_to(self.root).as_posix(),
                "sha256": sha256(authority_path),
            },
            "provenance_review": {
                "path": planner.f116_provenance_security_review_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": sha256(planner.f116_provenance_security_review_path()),
            },
            "plasma_review": {
                "path": planner.f116_plasma_scientific_review_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": sha256(planner.f116_plasma_scientific_review_path()),
            },
            "publication_audit": {
                "path": audit_path.relative_to(self.root).as_posix(),
                "sha256": sha256(audit_path),
            },
            "final_source_bundle": {
                "path": relative_bundle,
                "sha256": self.current_bundle_sha256,
                "verified_revisions": self.current_bundle_revisions,
            },
        }

    def publish_f118_current_source_authority(
        self, f116_binding: dict[str, object]
    ) -> dict[str, object]:
        """Publish F118 as current while retaining exact F116 predecessor bytes."""

        f116_evidence = json.loads(planner.f116_current_source_authority_path().read_text())
        implementation = f116_evidence["implementation"]
        bridge = implementation["intermediate_36140_bundle"]
        final = implementation["current_source_bundle"]
        predecessor = dict(final)
        predecessor.pop("candidate_path")
        predecessor["selected_as_current"] = False
        predecessor["role"] = "retained-non-current-predecessor"
        readme = planner.expected_path("source-archives/README.md")
        sums = planner.expected_path("source-archives/SHA256SUMS")
        before = {
            "readme_sha256": sha256(readme),
            "sha256sums_sha256": sha256(sums),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed": False,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
        }
        after = {
            "readme_sha256": sha256(readme),
            "sha256sums_sha256": sha256(sums),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed_exactly_once": True,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
            "historical_f116_preserved": True,
            "all_prior_checksum_entries_preserved": True,
            "sole_current_source_bundle": final["path"],
        }
        authority_path = planner.f118_current_source_authority_path()
        evidence = {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-evidence",
            "checkpoint": "F-118",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "generated_utc": utc_text(),
            "scope": {
                "relationship": "current-source-selection-only-supersession",
                "summary": "Select the exact committed F118 source without execution authority.",
                "preserves": planner.F118_SCOPE_PRESERVES,
                "does_not_authorize": planner.F118_SCOPE_DOES_NOT_AUTHORIZE,
            },
            "predecessor_authorities": {
                "historical_f116": {
                    key: f116_binding[key]
                    for key in (
                        "evidence", "provenance_review", "plasma_review",
                        "publication_audit",
                    )
                }
            },
            "implementation": {
                "publisher": implementation["publisher"],
                "committed_tools": implementation["committed_tools"],
                "intermediate_36140_bundle": bridge,
                "predecessor_current_source_bundle": predecessor,
                "current_source_bundle": final,
            },
            "source_archive_catalog": {"before": before, "after": after},
            "authorization": planner.F118_AUTHORIZATION,
            "validation": planner.F118_VALIDATION_CLAIMS,
            "publication_requirements": planner.F118_PUBLICATION_REQUIREMENTS,
        }
        write_json(authority_path, evidence)
        authority_path.chmod(0o444)
        verified = {
            "authorization_broadening": False,
            "bridge_selected_as_current": False,
            "predecessor_current_source_bundle_selected_as_current": False,
            "corrupt_c7_excluded": True,
            "current_source_selection_only": True,
            "final_bundle_sha256": final["sha256"],
            "final_head": final["head"],
            "historical_f115_preserved": True,
            "historical_f116_preserved": True,
        }
        review_specs = (
            (
                planner.f118_provenance_security_review_path(),
                "provenance-security",
                "approved-for-publication",
                "fixture-f118-provenance-reviewer",
            ),
            (
                planner.f118_plasma_scientific_review_path(),
                "plasma-scientific-continuation",
                "approved",
                "fixture-f118-plasma-reviewer",
            ),
        )
        for path, kind, decision, agent_id in review_specs:
            write_json(
                path,
                {
                    "schema_version": 1,
                    "record_type": (
                        "stage-i-current-source-authority-supersession-independent-review"
                    ),
                    "checkpoint": "F-118",
                    "execution_epoch": planner.EXECUTION_EPOCH,
                    "review_kind": kind,
                    "decision": decision,
                    "reviewed_candidate": {
                        "path": f"/tmp/{authority_path.name}.candidate",
                        "sha256": sha256(authority_path),
                    },
                    "published_f118": {
                        "path": str(authority_path),
                        "sha256": sha256(authority_path),
                    },
                    "reviewer": {
                        "agent_id": agent_id,
                        "identity": f"fixture independent {kind} reviewer",
                    },
                    "reviewed_utc": utc_text(),
                    "findings": ["Exact F118 source-only authority verified."],
                    "limitations": ["No prepare, submit, scheduler, or scientific authority."],
                    "verified": verified,
                },
            )
            path.chmod(0o444)
        f116_digests = {
            "evidence_sha256": f116_binding["evidence"]["sha256"],
            "provenance_review_sha256": f116_binding["provenance_review"]["sha256"],
            "plasma_review_sha256": f116_binding["plasma_review"]["sha256"],
            "publication_audit_sha256": f116_binding["publication_audit"]["sha256"],
        }
        audit_path = planner.f118_publication_audit_path()
        write_json(
            audit_path,
            {
                "schema_version": 1,
                "record_type": (
                    "stage-i-current-source-authority-supersession-publication-audit"
                ),
                "checkpoint": "F-118",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "published_utc": utc_text(),
                "artifact": {
                    "path": str(authority_path), "sha256": sha256(authority_path),
                    "mode": "0444", "links": 1,
                },
                "independent_reviews": {
                    "reviews_bind_exact_published_f118_sha256": sha256(authority_path),
                    "provenance_security": {
                        "path": str(planner.f118_provenance_security_review_path()),
                        "sha256": sha256(planner.f118_provenance_security_review_path()),
                        "mode": "0444", "links": 1,
                    },
                    "plasma_scientific_continuation": {
                        "path": str(planner.f118_plasma_scientific_review_path()),
                        "sha256": sha256(planner.f118_plasma_scientific_review_path()),
                        "mode": "0444", "links": 1,
                    },
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
                        "path": str(planner.expected_path(bridge["path"])),
                        "sha256": bridge["sha256"], "mode": "0644", "links": 1,
                        "head": bridge["head"], "role": "retained-non-current-bridge",
                        "selected_as_current": False,
                    },
                    "predecessor_current_source_bundle": {
                        "path": str(self.current_bundle), "sha256": final["sha256"],
                        "mode": "0644", "links": 1, "head": final["head"],
                        "role": "retained-non-current-predecessor",
                        "selected_as_current": False,
                    },
                    "current_source_bundle": {
                        "path": str(self.current_bundle), "sha256": final["sha256"],
                        "mode": "0644", "links": 1, "head": final["head"],
                        "selected_as_current": True,
                    },
                    "corrupt_c7_absent_from_active_checksum_ledger": True,
                    "sole_current_source_bundle": str(self.current_bundle),
                },
                "authority_and_enforcement": planner.F118_AUTHORIZATION,
                "publication": planner.F118_PUBLICATION_METHOD,
            },
        )
        audit_path.chmod(0o444)
        return {
            "checkpoint": "F-118",
            "evidence": {
                "path": authority_path.relative_to(self.root).as_posix(),
                "sha256": sha256(authority_path),
            },
            "provenance_review": {
                "path": planner.f118_provenance_security_review_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": sha256(planner.f118_provenance_security_review_path()),
            },
            "plasma_review": {
                "path": planner.f118_plasma_scientific_review_path().relative_to(
                    self.root
                ).as_posix(),
                "sha256": sha256(planner.f118_plasma_scientific_review_path()),
            },
            "publication_audit": {
                "path": audit_path.relative_to(self.root).as_posix(),
                "sha256": sha256(audit_path),
            },
            "final_source_bundle": {
                "path": final["path"],
                "sha256": final["sha256"],
                "verified_revisions": final["verified_revisions"],
            },
        }

    def lineage_digest(self) -> str:
        """Return the recost-compatible digest of fixture scientific lineages."""

        value = {}
        for case_id in planner.ALL_CASES:
            records = [
                (path, manifest)
                for path, manifest in self.manifests.items()
                if manifest["state"] == "recorded"
                and manifest["run"]["case_id"] == case_id
                and manifest["accounting"]["result"] in {"accepted", "clean_partial"}
                and not (
                    case_id == "R12"
                    and manifest["run"]["segment"]
                    == planner.R12_HISTORICAL_PARTIAL_SEGMENT
                    and manifest["job_id"] == planner.R12_HISTORICAL_PARTIAL_JOB_ID
                )
            ]
            records.sort(key=lambda item: planner.parse_segment(item[1]["run"]["segment"], "fixture")[0])
            if records:
                value[case_id] = [
                    {
                        "manifest": str(path),
                        "sha256": sha256(path),
                        "job_id": manifest["job_id"],
                        "segment": manifest["run"]["segment"],
                        "result": manifest["accounting"]["result"],
                        "final_time": manifest["scientific_inspection"]["final_time"],
                    }
                    for path, manifest in records
                ]
        return planner.recost_json_sha256(value)

    def next_profile(self, case_id: str) -> dict[str, object]:
        """Return one strict schema-2 recost next-profile fixture."""

        case_manifests = [
            (path, manifest)
            for path, manifest in self.manifests.items()
            if manifest["run"]["case_id"] == case_id
        ]
        recorded = [
            (path, manifest)
            for path, manifest in case_manifests
            if manifest["state"] == "recorded"
            and manifest["accounting"]["result"] in {"accepted", "clean_partial"}
            and not (
                case_id == "R12"
                and manifest["run"]["segment"] == planner.R12_HISTORICAL_PARTIAL_SEGMENT
                and manifest["job_id"] == planner.R12_HISTORICAL_PARTIAL_JOB_ID
            )
        ]
        recorded.sort(key=lambda item: planner.parse_segment(item[1]["run"]["segment"], "fixture")[0])
        r12_fresh_rerun = (
            case_id == "R12"
            and any(
                manifest["run"]["segment"] == planner.R12_HISTORICAL_PARTIAL_SEGMENT
                and manifest["job_id"] == planner.R12_HISTORICAL_PARTIAL_JOB_ID
                for _, manifest in case_manifests
            )
            and not recorded
        )
        next_index = max(
            (
                planner.parse_segment(manifest["run"]["segment"], "fixture")[0]
                for _, manifest in case_manifests
            ),
            default=-1,
        ) + 1
        parent_path = recorded[-1][0] if recorded else None
        parent = recorded[-1][1] if recorded else None
        start = (
            Decimal(str(parent["scientific_inspection"]["final_time"]))
            if parent is not None
            else Decimal("0")
        )
        increment = (
            self.r03_next_increment
            if case_id == "R03"
            else Decimal(str(self.profiles[case_id]["next_increment"]))
        )
        historical_f115_profile = case_id == "R03" and next_index == 2
        target = min(Decimal("10"), start + increment)
        if historical_f115_profile:
            target = Decimal("0.5")
        elif r12_fresh_rerun:
            target = planner.R12_FRESH_RERUN_TARGET
        segment_id = (
            planner.R03_F115_SEGMENT
            if historical_f115_profile
            else planner.R12_FRESH_RERUN_SEGMENT
            if r12_fresh_rerun
            else segment(next_index, planner.decimal_text(start), planner.decimal_text(target))
        )
        nodes = (
            planner.R12_FRESH_RERUN_NODES
            if r12_fresh_rerun
            else 1
            if case_id == "R03"
            else int(self.profiles[case_id]["nodes"])
        )
        walltime = "02:00:00" if case_id == "R03" else self.profiles[case_id]["walltime"]
        athena_walltime = (
            "01:50:00" if case_id == "R03" else self.profiles[case_id]["athena_walltime"]
        )
        terminal = parent["scientific_inspection"]["terminal_restart"] if parent is not None else None
        return {
            "acceptance_criterion": f"fixture reviewed acceptance for {case_id} at {target}",
            "acceptance_policy": "+".join(planner.policy_codes(case_id)),
            "athena_walltime": athena_walltime,
            "build_manifest": str(planner.build_manifest_path()),
            "build_manifest_sha256": planner.build_manifest_inventory_sha256(),
            "case_id": case_id,
            "controller_walltime_max_seconds": planner.MAX_SEGMENT_SECONDS,
            "cpus_per_task": planner.CPUS_PER_TASK,
            "estimated_storage_bytes": 1,
            "executable": str(planner.executable_path()),
            "executable_revision": planner.SOURCE_REVISION,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "input_file": f"inputs/cgl_lf_paper/{planner.EXPECTED_CASES[case_id][1]}",
            "input_revision": planner.SOURCE_REVISION,
            "input_sha256": planner.EXPECTED_CASES[case_id][3],
            "nodes": nodes,
            "output_layout": "rank-local",
            "segment": segment_id,
            "parent_job_id": parent["job_id"] if parent is not None else None,
            "parent_result": parent["accounting"]["result"] if parent is not None else None,
            "parent_segment": parent["run"]["segment"] if parent is not None else None,
            "restart_file": terminal["path"] if terminal is not None else None,
            "restart_file_sha256": terminal["sha256"] if terminal is not None else None,
            "restart_time": float(start) if parent is not None else None,
            "ranks_per_node": planner.RANKS_PER_NODE,
            "recommendation_basis": {"kind": "fixture-reviewed-schema-2-recost"},
            "time_tlim_target": float(target),
            "walltime": walltime,
            "source_bundle": str(
                planner.f115_source_bundle_path()
                if historical_f115_profile
                else self.current_bundle
            ),
            "source_bundle_sha256": (
                planner.R03_F115_SOURCE_BUNDLE_SHA256
                if historical_f115_profile
                else self.current_bundle_sha256
            ),
        }

    def scoped_projection_state(self) -> dict[str, object]:
        """Return the minimal canonical state needed to build a scoped-v2 fixture."""

        lineages = {case_id: ([], None) for case_id in planner.ALL_CASES}
        historical_r12 = []
        for manifest in self.manifests.values():
            accounting = manifest.get("accounting")
            if (
                manifest.get("state") != "recorded"
                or not isinstance(accounting, dict)
                or accounting.get("result") not in {"accepted", "clean_partial"}
            ):
                continue
            case_id = manifest["run"]["case_id"]
            segment_id = manifest["run"]["segment"]
            index, start, target = planner.parse_segment(segment_id, "fixture scoped segment")
            info = {
                "manifest": manifest,
                "case_id": case_id,
                "segment": segment_id,
                "index": index,
                "start": start,
                "target": target,
            }
            if (
                case_id == "R12"
                and segment_id == planner.R12_HISTORICAL_PARTIAL_SEGMENT
                and manifest["job_id"] == planner.R12_HISTORICAL_PARTIAL_JOB_ID
            ):
                historical_r12.append(info)
            else:
                lineages[case_id][0].append(info)
        for lineage, _ in lineages.values():
            lineage.sort(key=lambda item: item["index"])
        return {
            "lineages": lineages,
            "historical_inventory": {"R12": historical_r12},
        }

    def scoped_budget(self, profiles: list[dict[str, object]]) -> dict[str, object]:
        """Build one internally coherent scoped-v2 recost budget fixture."""

        state = self.scoped_projection_state()
        global_basis, r12_basis = planner.expected_scoped_measurement_bases(state)
        actual = sum(
            (Decimal(row["actual_node_hours"]) for row in self.ledger), Decimal("0")
        )
        reserved_by_case = {
            profile["case_id"]: (
                Decimal(
                    profile["nodes"]
                    * planner.walltime_seconds(profile["walltime"], "fixture scoped profile")
                )
                / Decimal(3600)
            )
            for profile in profiles
        }
        reserved = sum(reserved_by_case.values(), Decimal("0"))
        remaining = Decimal("0")
        breakdown = {}
        for case_id in planner.ALL_CASES:
            lineage, _ = state["lineages"][case_id]
            progress = min(
                Decimal("1"),
                max(Decimal("0"), planner.lineage_endpoint(lineage) / Decimal("10")),
            )
            remaining_time = Decimal("10") * (Decimal("1") - progress)
            cells = planner.resolution_cell_count(
                planner.EXPECTED_CASES[case_id][2], f"fixture scoped {case_id} resolution"
            )
            basis = r12_basis if case_id == "R12" else global_basis
            rate = Decimal(
                basis["normalized_node_hours_per_cell_per_simulation_time"]
            )
            observed = rate * Decimal(cells) * remaining_time
            authorized = reserved_by_case.get(case_id, Decimal("0"))
            projected = max(observed, authorized)
            remaining += projected
            breakdown[case_id] = {
                "matrix_full_case_node_hours_reference_only": format(
                    self.estimates[case_id], "f"
                ),
                "authenticated_progress_fraction": format(progress, "f"),
                "remaining_simulation_time": format(remaining_time, "f"),
                "projected_cells": str(cells),
                "projection_measurement_basis": basis,
                "observed_rate_projected_remaining_node_hours": format(observed, "f"),
                "authorized_profile_reserved_node_hours": format(authorized, "f"),
                "projected_remaining_node_hours": format(projected, "f"),
            }
        projected = actual + remaining
        return {
            "method": planner.NODE_HOUR_PROJECTION_METHOD,
            "measurement_basis": global_basis,
            "actual_stage_i_node_hours": format(actual, "f"),
            "authorized_wave_reserved_node_hours": format(reserved, "f"),
            "actual_plus_authorized_wave_node_hours": format(actual + reserved, "f"),
            "computed_remaining_stage_i_node_hours": format(remaining, "f"),
            "computed_stage_i_total_node_hours": format(projected, "f"),
            "promoted_stage_i_envelope_node_hours": format(
                planner.STAGE_I_BUDGET_NODE_HOURS, "f"
            ),
            "project_ceiling_node_hours": format(planner.PROJECT_BUDGET_NODE_HOURS, "f"),
            "computed_stage_i_margin_node_hours": format(
                planner.STAGE_I_BUDGET_NODE_HOURS - projected, "f"
            ),
            "case_breakdown": breakdown,
        }

    def write_strong_r17_readiness(
        self,
        lineage_sha256: str,
        projection_sha256: str,
        storage_sha256: str,
        profile: dict[str, object],
        current_source_authority: dict[str, object],
    ) -> tuple[dict[str, object] | None, str | None]:
        """Publish the recost-compatible R17 readiness/review/audit chain."""

        if not self.r17_ready:
            return None, None

        def binding(path: Path) -> dict[str, str]:
            return {
                "path": path.relative_to(self.root).as_posix(),
                "sha256": sha256(path),
            }

        def write_immutable_json(path: Path, value: object) -> None:
            write_json(path, value)
            path.chmod(0o444)

        logical_locations = [
            [lx1, lx2, lx3, 0]
            for lx1 in range(12)
            for lx2 in range(12)
            for lx3 in range(12)
        ]
        rank_local_outputs = []
        for rank in range(self.r17_output_count):
            path = self.root / f"qualification/output/rank_{rank:08d}/value.bin"
            write_bytes(path, f"R17 output rank {rank}\n".encode())
            rank_local_outputs.append(binding(path))
        if self.r17_duplicate_output_rank:
            duplicate = self.root / "qualification/output/rank_00000000/duplicate.bin"
            write_bytes(duplicate, b"R17 duplicate output rank\n")
            rank_local_outputs[-1] = binding(duplicate)
        rank_local_restarts = []
        for rank in range(64):
            path = self.root / f"qualification/restart/rank_{rank:08d}/value.rst"
            write_bytes(path, f"R17 restart rank {rank}\n".encode())
            rank_local_restarts.append(binding(path))
        output_inventory_sha256 = planner.recost_json_sha256(rank_local_outputs)
        restart_inventory_sha256 = planner.recost_json_sha256(rank_local_restarts)
        rank_inventory = [
            {
                "rank": rank,
                "rank_name": f"rank_{rank:08d}",
                "logical_meshblocks": logical_locations[rank * 27:(rank + 1) * 27],
            }
            for rank in range(64)
        ]
        meshblock_decomposition = {
            "schema_version": 1,
            "record_type": "stage-i-r17-decomposition-evidence",
            "resolution": "384x384x768",
            "mesh_shape": [384, 384, 768],
            "meshblock_shape": [32, 32, 64],
            "logical_meshblock_grid": [12, 12, 12],
            "logical_meshblocks": 1728,
            "ranks": 64,
            "meshblocks_per_rank": 27,
            "complete_block_rank_inventory": rank_inventory,
            "complete_block_rank_inventory_sha256": planner.compact_json_sha256(rank_inventory),
            "terminal_rank_local_output_inventory_sha256": output_inventory_sha256,
            "checks": {
                "exact_resolution": True,
                "exact_rank_count": True,
                "exact_meshblocks_per_rank": True,
                "complete_unique_logical_inventory": True,
            },
        }
        if self.r17_decomposition_mutation == "logical_count":
            meshblock_decomposition["logical_meshblocks"] = 1727
        elif self.r17_decomposition_mutation == "per_rank":
            moved = meshblock_decomposition["complete_block_rank_inventory"][0][
                "logical_meshblocks"
            ].pop()
            meshblock_decomposition["complete_block_rank_inventory"][1][
                "logical_meshblocks"
            ].append(moved)
        elif self.r17_decomposition_mutation == "duplicate":
            meshblock_decomposition["complete_block_rank_inventory"][-1][
                "logical_meshblocks"
            ][-1] = (
                meshblock_decomposition["complete_block_rank_inventory"][0][
                    "logical_meshblocks"
                ][0]
            )
        elif self.r17_decomposition_mutation == "digest":
            meshblock_decomposition["complete_block_rank_inventory_sha256"] = "0" * 64
        elif self.r17_decomposition_mutation == "output_binding":
            meshblock_decomposition["terminal_rank_local_output_inventory_sha256"] = "0" * 64
        if self.r17_decomposition_mutation in {"per_rank", "duplicate"}:
            meshblock_decomposition["complete_block_rank_inventory_sha256"] = (
                planner.compact_json_sha256(
                    meshblock_decomposition["complete_block_rank_inventory"]
                )
            )
        build_inventory = [
            {
                "name": path.name,
                "mode": file_evidence(path)["mode"],
                "sha256": sha256(path),
            }
            for path in sorted(planner.build_manifest_path().iterdir())
        ]
        build_inventory_sha256 = planner.recost_json_sha256(build_inventory)
        assert build_inventory_sha256 == profile["build_manifest_sha256"]

        intent_sha256 = "6" * 64
        execution_contract_sha256 = "7" * 64
        provenance = {
            "source": {"revision": planner.SOURCE_REVISION},
            "source_bundle": {"sha256": self.current_bundle_sha256},
            "matrix": {"sha256": planner.MATRIX_SHA256},
            "executable": {
                "revision": planner.SOURCE_REVISION,
                "sha256": planner.EXECUTABLE_SHA256,
            },
            "build_manifest": {"inventory_sha256": build_inventory_sha256},
        }
        provenance_sha256 = planner.compact_json_sha256(provenance)
        intent = {
            "case_id": "R17",
            "case_name": planner.EXPECTED_CASES["R17"][0],
            "profile_class": "scale_separation_384x384x768",
            "target_time": 0.25,
            "scientific_policy": "active_hardwall",
            "run_basename": "qualification_R17_n08",
            "input": {"sha256": planner.EXPECTED_CASES["R17"][3]},
            "execution_intent_sha256": intent_sha256,
            "execution_contract_sha256": execution_contract_sha256,
        }
        prepared_wave = {
            "provenance": provenance,
            "provenance_sha256": provenance_sha256,
            "waves": [
                {
                    "packets": [
                        {
                            "execution_intent_sha256": intent_sha256,
                            "execution_intent": intent,
                        }
                    ]
                }
            ],
        }
        prepared_wave_path = self.root / "qualification/r17/prepared_wave.json"
        write_json(prepared_wave_path, prepared_wave)
        qualification_evidence_path = (
            self.root / "qualification/r17/R17.qualification_evidence.json"
        )
        write_json(
            qualification_evidence_path,
            {
                "schema_version": 2,
                "record_type": "cgl_lf_stage_i_qualification_evidence",
                "project_root": str(self.root),
                "qualification_root": str(qualification_evidence_path.parent),
                "prepared_wave": {
                    "path": str(prepared_wave_path),
                    "sha256": sha256(prepared_wave_path),
                    "size_bytes": prepared_wave_path.stat().st_size,
                },
                "case_id": "R17",
                "target_time": 0.25,
                "results": [],
            },
        )
        physics_measurements = {
            "finite_rank_outputs": 64,
            "mass_relative_drift_max": "0",
            "mhd_user_mass_mismatch_max": "0",
            "lf_bad_counts_total": 0,
            "normalized_ct_divb_max": "0",
            "normalized_ct_divb_threshold": planner.R17_MAX_NORMALIZED_CT_DIVB_TEXT,
            "normalized_ct_divb_below_threshold": True,
        }
        scientific_checks = dict(planner.R17_SCIENTIFIC_CHECKS)
        if self.r17_physics_mutation == "rank_count":
            physics_measurements["finite_rank_outputs"] = 63
            scientific_checks["complete_rank_inventory"] = False
        elif self.r17_physics_mutation == "mass_drift":
            physics_measurements["mass_relative_drift_max"] = "1.0000000000000002e-12"
            scientific_checks["mass_conserved"] = False
        elif self.r17_physics_mutation == "mass_mismatch":
            physics_measurements["mhd_user_mass_mismatch_max"] = "1.0000000000000002e-12"
            scientific_checks["mass_conserved"] = False
        elif self.r17_physics_mutation == "lf_bad_counts":
            physics_measurements["lf_bad_counts_total"] = 1
            scientific_checks["strict_lf_failure_counters_zero"] = False
        elif self.r17_physics_mutation == "ct_threshold":
            physics_measurements["normalized_ct_divb_threshold"] = "1e-11"
        elif self.r17_physics_mutation == "ct_failure":
            physics_measurements["normalized_ct_divb_max"] = "9.9999999999999998e-13"
            physics_measurements["normalized_ct_divb_below_threshold"] = False
            scientific_checks["normalized_ct_divb_below_threshold"] = False
        scientific_path = self.root / "qualification/r17/R17.scientific_evidence.json"
        scientific = {
            "schema_version": 6,
            "record_type": "cgl_lf_stage_i_qualification_scientific_evidence",
            "case_id": "R17",
            "nodes": 8,
            "execution_intent_sha256": intent_sha256,
            "execution_contract_sha256": execution_contract_sha256,
            "terminal_rank_local_outputs": rank_local_outputs,
            "terminal_rank_local_output_inventory_sha256": output_inventory_sha256,
            "terminal_rank_local_restarts": rank_local_restarts,
            "terminal_rank_local_restart_inventory_sha256": restart_inventory_sha256,
            "r17_decomposition": meshblock_decomposition,
            "physics_measurements": physics_measurements,
            "checks": scientific_checks,
            "accepted_for_operational_qualification": True,
            "accepted_for_profile_selection": False,
        }
        write_json(scientific_path, scientific)

        job_id = "900001"
        submitted = NOW - timedelta(minutes=21)
        started = NOW - timedelta(minutes=20)
        completed = NOW - timedelta(minutes=10)
        measured = NOW - timedelta(minutes=9)
        scheduler_path = self.root / f"accounting/{job_id}.r17_qualification.sacct.txt"
        scheduler_row = "|".join(
            (
                job_id,
                "cglq_R17_n08",
                "COMPLETED",
                "0:0",
                "8",
                "600",
                utc_text(submitted),
                utc_text(completed),
            )
        ) + "\n"
        write_text(scheduler_path, scheduler_row, 0o444)
        account_scheduler_path = self.root / (
            "accounting/fabricated_account_query.txt"
            if self.r17_account_mutation == "path_scheduler"
            else f"accounting/{job_id}.r17_qualification.account.sacct.txt"
        )
        target_row = [
            job_id,
            "cglq_R17_n08",
            "COMPLETED",
            "0:0",
            "8",
            "600",
            utc_text(submitted),
            utc_text(started),
            utc_text(completed),
            "debug" if self.r17_account_mutation == "target_field" else "batch",
            planner.ACCOUNT,
            "fixture-user",
        ]
        account_rows = [target_row]
        if self.r17_account_mutation == "overlap":
            account_rows.append(
                [
                    "900000",
                    "foreign_account_job",
                    "COMPLETED",
                    "0:0",
                    "1",
                    "600",
                    utc_text(submitted),
                    utc_text(started),
                    utc_text(completed),
                    "batch",
                    planner.ACCOUNT,
                    "another-account-user",
                ]
            )
        account_raw = (
            planner.ACCOUNT_SCHEDULER_HEADER
            + "\n"
            + "\n".join("|".join(row) for row in account_rows)
            + "\n"
        )
        write_text(account_scheduler_path, account_raw, 0o444)
        account_jobs = planner.parse_account_scheduler_evidence(account_raw.encode())
        target_scheduler = {
            key: next(record for record in account_jobs if record["job_id"] == job_id)[key]
            for key in (
                "job_id",
                "job_name",
                "state",
                "exit_code",
                "nodes",
                "elapsed_seconds",
                "submit_utc",
                "start_utc",
                "end_utc",
                "partition",
                "account",
            )
        }
        query_contract = planner.r17_account_scheduler_query_contract(started, completed)
        if self.r17_account_mutation == "query_window":
            query_contract["start_utc"] = utc_text(started - timedelta(seconds=2))
        visibility_contract = {
            "private_data": (
                "jobs" if self.r17_account_mutation == "private_data" else "none"
            ),
            "all_users_job_visibility": True,
        }
        account_exclusivity_path = self.root / (
            "accounting/fabricated_account_exclusivity.json"
            if self.r17_account_mutation == "path_exclusivity"
            else "accounting/"
            f"mks24_stage_i_{planner.EXECUTION_EPOCH_SLUG}_R17_account_exclusivity_evidence.json"
        )
        account_exclusivity = {
            "schema_version": 1,
            "record_type": "stage-i-r17-account-exclusivity-evidence",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "measured_utc": utc_text(measured),
            "query_contract": query_contract,
            "visibility_contract": visibility_contract,
            "raw_account_scheduler_sha256": hashlib.sha256(account_raw.encode()).hexdigest(),
            "qualification_job": target_scheduler,
            "qualification_job_sha256": planner.compact_json_sha256(target_scheduler),
            "account_jobs": account_jobs,
            "account_jobs_sha256": planner.compact_json_sha256(account_jobs),
            "overlapping_job_ids": [record["job_id"] for record in account_jobs],
            "exclusive_entire_execution_interval": True,
        }
        write_json(account_exclusivity_path, account_exclusivity)
        account_exclusivity_path.chmod(0o444)
        measurement_author = "fixture qualification measurement author"
        restart_load_path = self.root / (
            "accounting/"
            f"mks24_stage_i_{planner.EXECUTION_EPOCH_SLUG}_R17_restart_load_evidence.json"
        )
        restart_measurements = {
            "rank_local_restart_inventory_sha256": restart_inventory_sha256,
            "loaded_rank_count": 64,
            "load_state": "COMPLETED",
            "load_exit_code": "0:0",
        }
        restart_load_evidence = {
            "schema_version": 1,
            "record_type": "stage-i-r17-restart-load-evidence",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "measured_utc": utc_text(measured),
            "measured_by": measurement_author,
            "job_id": job_id,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "build_manifest_inventory_sha256": build_inventory_sha256,
            "passed": True,
            "measurements": restart_measurements,
        }
        write_immutable_json(restart_load_path, restart_load_evidence)
        physics_validation_path = self.root / (
            "accounting/"
            f"mks24_stage_i_{planner.EXECUTION_EPOCH_SLUG}_R17_physics_validation_evidence.json"
        )
        physics_validation_measurements = {
            "rank_local_output_inventory_sha256": output_inventory_sha256,
            **physics_measurements,
        }
        physics_validation_evidence = {
            "schema_version": 1,
            "record_type": "stage-i-r17-physics-validation-evidence",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "measured_utc": utc_text(measured),
            "measured_by": measurement_author,
            "job_id": job_id,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "build_manifest_inventory_sha256": build_inventory_sha256,
            "passed": True,
            "measurements": physics_validation_measurements,
        }
        write_immutable_json(physics_validation_path, physics_validation_evidence)
        qualification_path = self.root / (
            "accounting/"
            f"mks24_stage_i_{planner.EXECUTION_EPOCH_SLUG}_R17_operational_qualification.json"
        )
        qualification_review_path = qualification_path.with_name(
            f"{qualification_path.name}.independent_review.json"
        )
        parameter_contract = planner.r17_parameter_contract()
        frozen_contract = {
            "case_id": "R17",
            "case_name": planner.EXPECTED_CASES["R17"][0],
            "profile_class": "scale_separation_384x384x768",
            "resolution": "384x384x768",
            "mesh_shape": [384, 384, 768],
            "meshblock_shape": [32, 32, 64],
            "target_time": 0.25,
            "scientific_policy": "active_hardwall",
            "run_basename": "qualification_R17_n08",
            "source_revision": planner.SOURCE_REVISION,
            "source_bundle_sha256": self.current_bundle_sha256,
            "matrix_sha256": planner.MATRIX_SHA256,
            "input_sha256": planner.EXPECTED_CASES["R17"][3],
            "provenance_sha256": provenance_sha256,
            "execution_intent_sha256": intent_sha256,
            "execution_contract_sha256": execution_contract_sha256,
            "parameter_contract": parameter_contract,
            "parameter_contract_sha256": planner.compact_json_sha256(parameter_contract),
            "executable_revision": planner.SOURCE_REVISION,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "build_manifest_inventory_sha256": build_inventory_sha256,
        }
        qualification = {
            "schema_version": 2,
            "record_type": "stage-i-r17-operational-qualification",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "completed_utc": utc_text(completed),
            "measured_utc": utc_text(measured),
            "measured_by": measurement_author,
            "job_id": job_id,
            "state": "COMPLETED",
            "exit_code": "0:0",
            "nodes": 8,
            "ranks": 64,
            "prepared_wave": binding(prepared_wave_path),
            "qualification_evidence": binding(qualification_evidence_path),
            "scientific_evidence": binding(scientific_path),
            "scheduler_evidence": binding(scheduler_path),
            "account_scheduler_evidence": {
                **binding(account_scheduler_path),
                "sha256": (
                    "0" * 64
                    if self.r17_account_mutation == "digest_scheduler"
                    else sha256(account_scheduler_path)
                ),
            },
            "account_exclusivity_evidence": {
                **binding(account_exclusivity_path),
                "sha256": (
                    "0" * 64
                    if self.r17_account_mutation == "digest_exclusivity"
                    else sha256(account_exclusivity_path)
                ),
            },
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "build_manifest_inventory": build_inventory,
            "build_manifest_inventory_sha256": build_inventory_sha256,
            "rank_local_outputs": rank_local_outputs,
            "rank_local_output_inventory_sha256": output_inventory_sha256,
            "rank_local_restarts": rank_local_restarts,
            "rank_local_restart_inventory_sha256": restart_inventory_sha256,
            "decomposition_evidence": meshblock_decomposition,
            "restart_load_evidence": binding(restart_load_path),
            "physics_validation_evidence": binding(physics_validation_path),
            "frozen_science_build_contract": frozen_contract,
            "independent_review_contract": {
                "required": True,
                "path": qualification_review_path.relative_to(self.root).as_posix(),
                "mode": "0444",
                "schema_version": 1,
                "record_type": (
                    "stage-i-r17-operational-qualification-independent-review"
                ),
                "execution_epoch": planner.EXECUTION_EPOCH,
                "decision": "approved",
                "candidate_path": str(qualification_path),
                "candidate_sha256_required": True,
                "reviewer_must_differ_from": [measurement_author],
                "reviewed_after_utc": utc_text(measured),
            },
            "authority": {
                "r17_launch_authorized": False,
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            },
        }
        write_immutable_json(qualification_path, qualification)
        write_immutable_json(
            qualification_review_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-r17-operational-qualification-independent-review",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "reviewed_utc": utc_text(),
                "decision": "approved",
                "reviewer": "fixture qualification reviewer",
                "candidate": {
                    "path": str(qualification_path),
                    "sha256": sha256(qualification_path),
                },
            },
        )
        readiness_path = planner.r17_readiness_path()
        readiness = {
            "schema_version": 1,
            "record_type": "stage-i-r17-readiness",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "root": str(self.root),
            "generated_utc": utc_text(),
            "expires_utc": utc_text(NOW + timedelta(hours=1)),
            "reviewed_by": "fixture readiness reviewer",
            "predecessor_lineages_sha256": lineage_sha256,
            "storage_evidence_sha256": storage_sha256,
            "computed_projection_sha256": projection_sha256,
            "executable_sha256": planner.EXECUTABLE_SHA256,
            "build_manifest_sha256": profile["build_manifest_sha256"],
            "required_retained_bytes": planner.R17_REQUIRED_RETENTION_BYTES,
            "nodes": 8,
            "ranks": 64,
            "storage_ready": True,
            "node_hour_ready": True,
            "rank_64_ready": True,
            "operational_qualification": {
                "path": qualification_path.relative_to(self.root).as_posix(),
                "sha256": sha256(qualification_path),
            },
            "operational_qualification_review": {
                "path": qualification_review_path.relative_to(self.root).as_posix(),
                "sha256": sha256(qualification_review_path),
            },
            "current_source_authority": (
                self.r17_source_authority_override
                if self.r17_source_authority_override is not None
                else current_source_authority
            ),
        }
        write_json(readiness_path, readiness)
        readiness_path.chmod(0o444)
        readiness_review_path = planner.recost_independent_review_path(readiness_path)
        write_json(
            readiness_review_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-r17-readiness-independent-review",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "reviewed_utc": utc_text(),
                "decision": "approved-for-publication",
                "reviewer": "fixture readiness reviewer",
                "candidate": {"path": str(readiness_path), "sha256": sha256(readiness_path)},
            },
        )
        readiness_review_path.chmod(0o444)
        readiness_audit_path = planner.recost_publication_audit_path(readiness_path)
        write_json(
            readiness_audit_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-r17-readiness-publication-audit",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "published_utc": utc_text(),
                "artifact": {
                    "path": str(readiness_path),
                    "sha256": sha256(readiness_path),
                    "mode": "0444",
                    "links": 1,
                },
                "independent_review": {
                    "path": str(readiness_review_path),
                    "sha256": sha256(readiness_review_path),
                    "mode": "0444",
                    "links": 1,
                },
                "authority": {
                    "r17_launch_authorized": False,
                    "scheduler_mutation_authorized": False,
                    "canonical_mutation_authorized": False,
                },
            },
        )
        readiness_audit_path.chmod(0o444)
        expanded_qualification = {
            **qualification,
            "authenticated_rank_local_outputs": qualification["rank_local_outputs"],
            "authenticated_rank_local_restarts": qualification["rank_local_restarts"],
            "authenticated_decomposition_evidence": qualification[
                "decomposition_evidence"
            ],
            "authenticated_account_exclusivity_evidence": account_exclusivity,
            "authenticated_validation_evidence": {
                "restart_load_evidence": {
                    **qualification["restart_load_evidence"],
                    "measured_utc": utc_text(measured),
                    "measured_by": measurement_author,
                    "measurements": restart_measurements,
                },
                "physics_validation_evidence": {
                    **qualification["physics_validation_evidence"],
                    "measured_utc": utc_text(measured),
                    "measured_by": measurement_author,
                    "measurements": physics_validation_measurements,
                },
            },
            "reviewed_by": "fixture qualification reviewer",
        }
        if self.r17_account_mutation == "expanded_drift":
            expanded_qualification["authenticated_account_exclusivity_evidence"] = {
                **account_exclusivity,
                "overlapping_job_ids": [],
            }
        return (
            {
                **readiness,
                "operational_qualification_evidence": expanded_qualification,
                "publication_chain": {
                    "independent_review_sha256": sha256(readiness_review_path),
                    "publication_audit_sha256": sha256(readiness_audit_path),
                    "reviewed_by": "fixture readiness reviewer",
                    "published_utc": utc_text(),
                },
            },
            sha256(readiness_path),
        )

    def refresh(self) -> None:
        """Rewrite canonical stores and the sole promoted schema-2 recost chain."""

        for path, manifest in self.manifests.items():
            write_json(path, manifest)
        ledger_buffer = io.StringIO()
        writer = csv.DictWriter(ledger_buffer, fieldnames=planner.LEDGER_COLUMNS)
        writer.writeheader()
        writer.writerows(self.ledger)
        write_text(planner.ledger_path(), ledger_buffer.getvalue())
        write_json(planner.reservations_path(), self.reservations)
        ledger_evidence = file_evidence(planner.ledger_path())
        reservations_evidence = file_evidence(planner.reservations_path())
        qualification_evidence = file_evidence(planner.qualification_path())
        f116_source_authority = self.publish_f116_current_source_authority()
        current_source_authority = self.publish_f118_current_source_authority(
            f116_source_authority
        )
        active = [
            item for item in self.reservations if item["state"] in {"prepared", "submitted"}
        ]
        counts = {
            "transactions": 0,
            "reservations": len(self.reservations),
            "active_reservations": len(active),
            "ledger_rows": len(self.ledger),
            "manifests": len(self.manifests),
        }
        reconciliation = {
            "execution_epoch": planner.EXECUTION_EPOCH,
            "root": str(self.root),
            "qualification": {
                "state": "approved",
                "path": str(planner.qualification_path()),
                "sha256": qualification_evidence["sha256"],
                "approved_executable_revision": planner.SOURCE_REVISION,
                "approved_executable_sha256": planner.EXECUTABLE_SHA256,
            },
            "consistent": True,
            "counts": counts,
            "issues": [],
        }
        completed_predecessors = all(
            any(
                manifest["state"] == "recorded"
                and manifest["run"]["case_id"] == case_id
                and manifest["accounting"]["result"] == "accepted"
                and Decimal(str(manifest["scientific_inspection"]["final_time"])) == Decimal("10")
                for manifest in self.manifests.values()
            )
            for case_id in tuple(f"R{index:02d}" for index in range(2, 17))
        )
        cases = ["R17"] if completed_predecessors else list(self.recommended_cases)
        profiles = [self.next_profile(case_id) for case_id in cases]
        mode = "sole-next-profile" if len(profiles) == 1 else "bounded-wave"
        budget = self.scoped_budget(profiles)
        if self.budget_mutator is not None:
            self.budget_mutator(budget)
        projection_sha256 = planner.recost_json_sha256(budget)
        lineage_sha256 = self.lineage_digest()
        storage_sha256 = "8" * 64
        r17_readiness, r17_readiness_sha256 = self.write_strong_r17_readiness(
            lineage_sha256,
            projection_sha256,
            storage_sha256,
            profiles[0],
            current_source_authority,
        )
        manifest_bindings = [
            {"path": path.relative_to(self.root).as_posix(), "sha256": sha256(path)}
            for path in sorted(self.manifests)
        ]
        barrier_row = self.ledger[-1]
        generated_utc = utc_text()
        expires_utc = utc_text(NOW + timedelta(hours=1))
        request_path = self.root / "accounting/fixture_F200_recost_request.json"
        write_json(
            request_path,
            {
                "schema_version": 2,
                "record_type": "stage-i-recost-recommendation-request",
                "checkpoint": "F-200",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "generated_utc": generated_utc,
                "expires_utc": expires_utc,
                "requested_by": "fixture recost requester",
                "scope": "fixture reviewed drained-barrier wave",
            },
        )
        request_review_path = request_path.with_name(
            f"{request_path.name}.independent_review.json"
        )
        write_json(
            request_review_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-recost-request-independent-review",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "reviewed_utc": generated_utc,
                "decision": "approved-for-evidence-generation",
                "reviewer": {
                    "agent_id": self.request_reviewer_agent,
                    "role": "independent recost request reviewer",
                    "declared_process_independence": True,
                    "identity_assurance": planner.REQUEST_REVIEW_IDENTITY_ASSURANCE,
                    "identity_assurance_limitation": (
                        planner.REQUEST_REVIEW_IDENTITY_LIMITATION
                    ),
                },
                "candidate": {
                    "path": str(request_path),
                    "sha256": sha256(request_path),
                },
                "scope": {"non_authorizing": True},
            },
        )
        request_review_path.chmod(0o444)
        recost = {
            "schema_version": 2,
            "record_type": "stage-i-recost-recommendation-evidence",
            "checkpoint": "F-200",
            "artifact_name": self.recost_path.name,
            "execution_epoch": planner.EXECUTION_EPOCH,
            "generated_utc": generated_utc,
            "expires_utc": expires_utc,
            "requested_by": "fixture recost requester",
            "scope": "fixture reviewed drained-barrier wave",
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
                "mode": mode,
                "authorizing": False,
                "recommended_next_profiles": profiles,
                "bounded_concurrency": {
                    "max_active_segments": planner.MAX_LANES,
                    "max_wave_nodes": sum(profile["nodes"] for profile in profiles),
                    "r17_exclusive_and_last": True,
                },
                "controller_consumption_state": "reviewed and promoted; controller-mediated consumption",
                "non_authorizing_reason": "planner output and recost evidence do not authorize mutation",
            },
            "barrier": {
                "job_ids": [barrier_row["job_id"]],
                "recorded_segments": [
                    {
                        "job_id": barrier_row["job_id"],
                        "case_id": barrier_row["case_id"],
                        "segment": barrier_row["segment"],
                        "result": barrier_row["result"],
                    }
                ],
                "scheduler_evidence": [],
            },
            "budget": budget,
            "storage": {
                "available_bytes": planner.R17_REQUIRED_RETENTION_BYTES * 3,
                "retained_stage_i_bytes": 1,
                "required_safety_bytes": 1,
                "projected_authorized_wave_growth_bytes": 1,
                "headroom_after_authorized_wave_and_safety_bytes": (
                    planner.R17_REQUIRED_RETENTION_BYTES * 3 - 2
                ),
            },
            "ledger": {
                "rows": len(self.ledger),
                "sha256": ledger_evidence["sha256"],
                "cumulative_stage_i_node_hours": self.ledger[-1][
                    "cumulative_stage_i_node_hours"
                ],
            },
            "reservations": {
                "rows": len(self.reservations),
                "sha256": reservations_evidence["sha256"],
                "active": 0,
            },
            "manifests": {
                "rows": len(manifest_bindings),
                "bindings": manifest_bindings,
                "authenticated_lineages_sha256": lineage_sha256,
            },
            "r17_readiness": r17_readiness,
            "promoted_f113": {},
            "reconcile": reconciliation,
            "provenance": {
                "request_sha256": sha256(request_path),
                "request_independent_review_sha256": sha256(request_review_path),
                "generator_sha256": self.generator_sha256,
                "generator_revision": self.controller_revision,
                "stage_i_helper_sha256": self.controller_sha256,
                "stage_i_helper_revision": self.controller_revision,
                "matrix_sha256": planner.MATRIX_SHA256,
                "matrix_revision": self.controller_revision,
                "source_bundle_sha256": self.current_bundle_sha256,
                "source_bundle_verified_revisions": self.current_bundle_revisions,
                "source_authority": current_source_authority,
                "qualification_approval_sha256": qualification_evidence["sha256"],
                "storage_evidence_sha256": storage_sha256,
                "reconciliation_sha256": "a" * 64,
                "ledger_sha256": ledger_evidence["sha256"],
                "reservations_sha256": reservations_evidence["sha256"],
                "authenticated_lineages_sha256": lineage_sha256,
                "computed_projection_sha256": projection_sha256,
                "r17_readiness_evidence_sha256": r17_readiness_sha256,
            },
        }
        write_json(self.recost_path, recost)
        self.recost_path.chmod(0o444)
        review_path = planner.recost_independent_review_path(self.recost_path)
        write_json(
            review_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-recost-recommendation-independent-review",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "reviewed_utc": utc_text(),
                "decision": "approved-for-publication",
                "reviewer": {
                    "agent_id": self.recost_reviewer_agent,
                    "independent_from_generator": True,
                },
                "candidate": {"path": str(self.recost_path), "sha256": sha256(self.recost_path)},
                "scope": {"non_authorizing": True},
            },
        )
        review_path.chmod(0o444)
        audit_path = planner.recost_publication_audit_path(self.recost_path)
        forensic_path = (
            self.root
            / f"accounting/mks24_stage_i_{planner.EXECUTION_EPOCH_SLUG}_recost_forensics"
            / f"fixture-f200-publish.{self.recost_path.name}.forensic"
        )
        write_bytes(forensic_path, self.recost_path.read_bytes(), 0o444)
        retained_generator = self.root / "accounting/utilities/cgl_lf_stage_i_recost.py"
        write_bytes(
            retained_generator,
            (self.helper.parents[2] / "scripts/frontier/cgl_lf_stage_i_recost.py").read_bytes(),
            0o755,
        )
        review_binding = {
            "path": str(review_path),
            "sha256": sha256(review_path),
            "mode": "0444",
            "links": 1,
        }
        context = {
            "schema_version": 2,
            "artifact_name": self.recost_path.name,
            "checkpoint": "F-200",
            "generated_utc": generated_utc,
            "expires_utc": expires_utc,
            "request": {"path": str(request_path), "sha256": sha256(request_path)},
            "artifact": {
                "basename": self.recost_path.name,
                "sha256": sha256(self.recost_path),
            },
            "independent_review": review_binding,
            "authority": recost["authority"],
            "publication_requirements": recost["publication_requirements"],
            "recommendations": recost["recommendations"],
            "barrier": recost["barrier"],
            "provenance": recost["provenance"],
            "predecessor_recost": recost["predecessor_recost"],
            "reconciliation": recost["reconcile"],
            "budget": recost["budget"],
            "storage": recost["storage"],
            "r17_readiness": recost["r17_readiness"],
            "controller_enforcement": {
                "state": recost["recommendations"]["controller_consumption_state"],
                "launch_authority": False,
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            },
        }
        write_json(
            audit_path,
            {
                "schema_version": 1,
                "record_type": "stage-i-recost-recommendation-publication-audit",
                "execution_epoch": planner.EXECUTION_EPOCH,
                "transaction_id": "fixture-f200-publish",
                "published_utc": utc_text(),
                "artifact": {
                    "path": str(self.recost_path),
                    "sha256": sha256(self.recost_path),
                    "mode": "0444",
                    "links": 1,
                },
                "recost_recommendations": recost["recommendations"],
                "independent_review": review_binding,
                "authority": {
                    "action_authority": False,
                    "scheduler_mutation_authorized": False,
                    "canonical_mutation_authorized": False,
                },
                "counts": counts,
                "generator": {
                    "path": str(retained_generator),
                    "sha256": self.generator_sha256,
                    "mode": "0755",
                    "revision": self.controller_revision,
                },
                "scheduler_evidence": [],
                "source_bundle": {
                    "path": str(self.current_bundle),
                    "sha256": self.current_bundle_sha256,
                    "mode": "0644",
                    "verified_revisions": self.current_bundle_revisions,
                },
                "stage_i_helper": {
                    "path": str(self.helper),
                    "sha256": self.controller_sha256,
                    "committed": True,
                    "reconcile_execution": "descriptor",
                    "revision": self.controller_revision,
                },
                "utility": {
                    "path": str(self.helper),
                    "sha256": self.controller_sha256,
                    "execution": "authenticated-descriptor",
                    "committed": True,
                },
                "forensic_copy": {
                    "path": str(forensic_path),
                    "sha256": sha256(self.recost_path),
                    "mode": "0444",
                    "links": 1,
                    "generalized_vectors": "exact-artifact-payload",
                },
                "publication": (
                    "same-directory-link-fsync-copy-exchange-forensic-retirement-fsync"
                ),
                "generalized_publication_context": context,
            },
        )
        audit_path.chmod(0o444)

    def plan(self) -> dict[str, object]:
        """Build one plan through the production path API."""

        return planner.plan_from_paths(
            self.recost_path,
            planner.recost_independent_review_path(self.recost_path),
            planner.recost_publication_audit_path(self.recost_path),
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
    value.next_job = int(planner.R12_HISTORICAL_PARTIAL_JOB_ID)
    value.add_recorded(
        "R12",
        planner.R12_HISTORICAL_PARTIAL_SEGMENT,
        final_time="0.1371931229426507",
        result="clean_partial",
        nodes=planner.R12_FRESH_RERUN_NODES,
        elapsed_seconds=6657,
    )
    value.refresh()
    return value


def recost_paths(campaign: Campaign) -> tuple[Path, Path, Path]:
    """Return one fixture's exact promoted recost publication chain."""

    return (
        campaign.recost_path,
        planner.recost_independent_review_path(campaign.recost_path),
        planner.recost_publication_audit_path(campaign.recost_path),
    )


def add_credible_above_envelope_measurements(campaign: Campaign) -> None:
    """Add the retained R04/R16 measurements that make scoped-v2 exceed 1400."""

    campaign.add_recorded(
        "R04",
        segment(0, "0", "0.25"),
        nodes=4,
        walltime="02:00:00",
        athena_walltime="01:50:00",
        elapsed_seconds=1360,
    )
    campaign.add_recorded(
        "R16",
        segment(0, "0", "1.5"),
        nodes=1,
        walltime="02:00:00",
        athena_walltime="01:50:00",
        elapsed_seconds=2631,
    )


def current_source_binding(campaign: Campaign) -> dict[str, object]:
    """Return the exact F118 binding retained by the promoted recost."""

    recost = json.loads(campaign.recost_path.read_text())
    return recost["provenance"]["source_authority"]


def republish_changed_f118_review_chain(campaign: Campaign) -> dict[str, object]:
    """Rebind F118 reviews/audit after an intentional semantic mutation."""

    authority_path = planner.f118_current_source_authority_path()
    authority_sha = sha256(authority_path)
    review_paths = (
        planner.f118_provenance_security_review_path(),
        planner.f118_plasma_scientific_review_path(),
    )
    for path in review_paths:
        review = json.loads(path.read_text())
        review["reviewed_candidate"]["sha256"] = authority_sha
        review["published_f118"]["sha256"] = authority_sha
        write_json(path, review)
        path.chmod(0o444)
    audit_path = planner.f118_publication_audit_path()
    audit = json.loads(audit_path.read_text())
    audit["artifact"]["sha256"] = authority_sha
    audit["independent_reviews"]["reviews_bind_exact_published_f118_sha256"] = authority_sha
    audit["independent_reviews"]["provenance_security"]["sha256"] = sha256(review_paths[0])
    audit["independent_reviews"]["plasma_scientific_continuation"]["sha256"] = sha256(
        review_paths[1]
    )
    write_json(audit_path, audit)
    audit_path.chmod(0o444)
    binding = current_source_binding(campaign)
    binding["evidence"]["sha256"] = authority_sha
    binding["provenance_review"]["sha256"] = sha256(review_paths[0])
    binding["plasma_review"]["sha256"] = sha256(review_paths[1])
    binding["publication_audit"]["sha256"] = sha256(audit_path)
    return binding


def test_schema2_recost_and_f116_emit_exact_command_free_wave(campaign):
    plan = campaign.plan()
    assert plan["schema_version"] == 4
    assert [packet["case_id"] for packet in plan["wave"]["packets"]] == [
        "R03", "R04", "R12", "R16"
    ]
    assert plan["wave"]["planned_nodes"] == 10
    r12 = next(packet for packet in plan["wave"]["packets"] if packet["case_id"] == "R12")
    assert r12["lineage"] == {
        "kind": "fresh",
        "segment": planner.R12_FRESH_RERUN_SEGMENT,
        "start_time": "0",
        "target_time": "0.12",
        "restart_required": False,
        "continuation_provenance": None,
    }
    assert r12["allocation"]["nodes"] == planner.R12_FRESH_RERUN_NODES
    assert r12["allocation"]["total_ranks"] == planner.R12_FRESH_RERUN_RANKS
    assert r12["allocation"]["walltime"] == planner.R12_FRESH_RERUN_WALLTIME
    assert (
        r12["allocation"]["athena_walltime"]
        == planner.R12_FRESH_RERUN_ATHENA_WALLTIME
    )
    assert "continuation_eligibility_waiver" not in r12
    assert plan["evidence"]["r12_fresh_rerun_transition"]["waiver_authorized"] is False
    assert plan["observed_state"]["credible_remaining_campaign_projection"] == json.loads(
        campaign.recost_path.read_text()
    )["budget"]
    f116 = plan["wave"]["packets"][1]["production_provenance"][
        "current_source_authority"
    ]
    assert f116 == current_source_binding(campaign)
    assert f116["final_source_bundle"]["sha256"] == campaign.current_bundle_sha256
    assert plan["wave"]["packets"][1]["production_provenance"][
        "current_source_authority_evidence"
    ]["artifact"]["path"] == str(planner.f118_current_source_authority_path())
    assert plan["wave"]["packets"][0]["production_provenance"]["source_bundle"][
        "path"
    ] == str(planner.f115_source_bundle_path())
    assert plan["wave"]["packets"][0]["authorization_state"] == (
        "planning_only_f115_profile_reference_non_authorizing"
    )
    assert plan["wave"]["packets"][0]["authorization_evidence"][
        "sole_next_segment_profile"
    ]["segment"] == planner.R03_F115_SEGMENT
    serialized = json.dumps(plan)
    assert "sbatch " not in serialized
    assert "srun " not in serialized


def test_wave_planner_consumes_scoped_v2_global_and_r12_local_bases(campaign):
    budget = campaign.plan()["observed_state"]["credible_remaining_campaign_projection"]
    global_basis = budget["measurement_basis"]
    r12_basis = budget["case_breakdown"]["R12"]["projection_measurement_basis"]

    assert budget["method"] == planner.NODE_HOUR_PROJECTION_METHOD
    assert global_basis["case_id"] != "R12"
    assert r12_basis["case_id"] == "R12"
    assert r12_basis["job_id"] == planner.R12_HISTORICAL_PARTIAL_JOB_ID
    assert all(
        budget["case_breakdown"][case_id]["projection_measurement_basis"]
        == global_basis
        for case_id in planner.ALL_CASES
        if case_id != "R12"
    )


def test_wave_planner_scoped_v2_schema_matches_recost_producer(campaign, monkeypatch):
    recost = load_utility(RECOST_UTILITY, "cgl_lf_stage_i_recost_for_wave_budget_parity")
    monkeypatch.setattr(
        recost,
        "manifest_identity",
        lambda manifest: (
            manifest["job_id"],
            manifest["run"]["case_id"],
            manifest["run"]["segment"],
            manifest["accounting"]["result"],
        ),
    )
    monkeypatch.setattr(
        recost,
        "final_time",
        lambda manifest: manifest["scientific_inspection"]["final_time"],
    )
    state = campaign.scoped_projection_state()
    lineages = {
        case_id: [info["manifest"] for info in state["lineages"][case_id][0]]
        for case_id in planner.ALL_CASES
    }
    if not lineages["R12"]:
        lineages["R12"] = [
            info["manifest"] for info in state["historical_inventory"]["R12"]
        ]
    matrix = {
        case_id: {
            "_cell_count": planner.resolution_cell_count(
                planner.EXPECTED_CASES[case_id][2], f"fixture parity {case_id} resolution"
            ),
            "_estimated_node_hours": campaign.estimates[case_id],
        }
        for case_id in planner.ALL_CASES
    }
    profiles = [campaign.next_profile(case_id) for case_id in campaign.recommended_cases]
    reserved = sum(
        (
            Decimal(
                profile["nodes"]
                * planner.walltime_seconds(profile["walltime"], "fixture parity profile")
            )
            / Decimal(3600)
            for profile in profiles
        ),
        Decimal("0"),
    )

    produced = recost.calculate_budget(
        campaign.ledger,
        lineages,
        matrix,
        profiles,
        reserved,
        planner.STAGE_I_BUDGET_NODE_HOURS,
        planner.PROJECT_BUDGET_NODE_HOURS,
    )

    def retain_produced_budget(budget):
        budget.clear()
        budget.update(deepcopy(produced))

    campaign.budget_mutator = retain_produced_budget
    campaign.refresh()
    assert (
        campaign.plan()["observed_state"]["credible_remaining_campaign_projection"]
        == produced
    )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("legacy_method", "projection method is not scoped-v2"),
        ("r12_as_global", "global measurement basis differs"),
        ("global_as_r12", "R12 measurement basis differs"),
        ("r12_for_non_r12", "R04 measurement basis differs"),
        ("arithmetic", "R04 projection arithmetic differs"),
    ],
)
def test_wave_planner_rejects_scoped_v2_method_or_basis_confusion(
    campaign, mutation, message
):
    def mutate(budget):
        global_basis = deepcopy(budget["measurement_basis"])
        r12_basis = deepcopy(
            budget["case_breakdown"]["R12"]["projection_measurement_basis"]
        )
        if mutation == "legacy_method":
            budget["method"] = "observed-stage-i-node-hour-rate-v1"
        elif mutation == "r12_as_global":
            budget["measurement_basis"] = r12_basis
        elif mutation == "global_as_r12":
            budget["case_breakdown"]["R12"]["projection_measurement_basis"] = global_basis
        elif mutation == "r12_for_non_r12":
            budget["case_breakdown"]["R04"]["projection_measurement_basis"] = r12_basis
        else:
            budget["case_breakdown"]["R04"][
                "observed_rate_projected_remaining_node_hours"
            ] = "0"

    campaign.budget_mutator = mutate
    campaign.refresh()
    with pytest.raises(ValueError, match=message):
        campaign.plan()


def test_above_envelope_scoped_v2_is_allowed_only_for_mandatory_fresh_r12_wave(
    campaign,
):
    add_credible_above_envelope_measurements(campaign)
    campaign.refresh()
    budget = json.loads(campaign.recost_path.read_text())["budget"]

    assert Decimal(budget["computed_stage_i_total_node_hours"]) > Decimal("1400")
    assert Decimal(budget["actual_plus_authorized_wave_node_hours"]) < Decimal("1400")
    assert Decimal(budget["computed_stage_i_total_node_hours"]) < Decimal("4000")
    plan = campaign.plan()
    assert next(
        packet for packet in plan["wave"]["packets"] if packet["case_id"] == "R12"
    )["lineage"]["kind"] == "fresh"


def test_successor_wave_fails_closed_until_fresh_recost_envelope_transition(
    campaign, monkeypatch
):
    add_credible_above_envelope_measurements(campaign)
    campaign.add_recorded(
        "R12",
        planner.R12_FRESH_RERUN_SEGMENT,
        nodes=planner.R12_FRESH_RERUN_NODES,
        walltime=planner.R12_FRESH_RERUN_WALLTIME,
        athena_walltime=planner.R12_FRESH_RERUN_ATHENA_WALLTIME,
        elapsed_seconds=6657,
    )
    campaign.refresh()
    budget = json.loads(campaign.recost_path.read_text())["budget"]
    assert Decimal(budget["computed_stage_i_total_node_hours"]) > Decimal("1400")
    with pytest.raises(ValueError, match="exact mandatory fresh R12 calibration wave"):
        campaign.plan()

    monkeypatch.setattr(planner, "STAGE_I_BUDGET_NODE_HOURS", Decimal("2000"))
    campaign.refresh()
    transitioned = campaign.plan()
    transitioned_budget = transitioned["observed_state"][
        "credible_remaining_campaign_projection"
    ]
    assert Decimal(transitioned_budget["computed_stage_i_total_node_hours"]) < Decimal(
        transitioned_budget["promoted_stage_i_envelope_node_hours"]
    )
    assert next(
        packet for packet in transitioned["wave"]["packets"] if packet["case_id"] == "R12"
    )["lineage"]["kind"] == "continuation"


def test_fresh_r12_calibration_exception_never_bypasses_project_ceiling(
    campaign, monkeypatch
):
    add_credible_above_envelope_measurements(campaign)
    monkeypatch.setattr(planner, "PROJECT_BUDGET_NODE_HOURS", Decimal("1500"))
    campaign.refresh()
    with pytest.raises(ValueError, match="exceeds project ceiling"):
        campaign.plan()


def test_r03_post_f115_recost_continues_normally_from_accepted_and_clean_partial(campaign):
    r03_s00 = next(
        path
        for path, manifest in campaign.manifests.items()
        if manifest["run"]["case_id"] == "R03"
        and manifest["run"]["segment"] == segment(0, "0", "0.5")
    )
    r03_s02 = campaign.add_recorded(
        "R03",
        planner.R03_F115_SEGMENT,
        parent_path=r03_s00,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    campaign.refresh()

    packet = next(item for item in campaign.plan()["wave"]["packets"] if item["case_id"] == "R03")
    assert packet["lineage"]["segment"] == segment(3, "0.5", "0.75")
    assert packet["lineage"]["continuation_provenance"]["parent_manifest"]["path"] == str(r03_s02)
    assert packet["authorization_state"] == "planning_only_reviewed_recost_profile_non_authorizing"
    assert packet["authorization_evidence"] is None
    assert packet["production_provenance"]["source_bundle"]["path"] == str(campaign.current_bundle)
    assert packet["production_provenance"]["current_source_authority"] == current_source_binding(campaign)

    r03_s03 = campaign.add_recorded(
        "R03",
        segment(3, "0.5", "0.75"),
        final_time="0.625",
        result="clean_partial",
        parent_path=r03_s02,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    campaign.refresh()

    packet = next(item for item in campaign.plan()["wave"]["packets"] if item["case_id"] == "R03")
    assert packet["lineage"]["segment"] == segment(4, "0.625", "0.875")
    assert packet["lineage"]["continuation_provenance"]["parent_manifest"]["path"] == str(r03_s03)
    assert packet["lineage"]["continuation_provenance"]["restart_time"] == "0.625"
    assert packet["authorization_state"] == "planning_only_reviewed_recost_profile_non_authorizing"
    assert packet["authorization_evidence"] is None


def test_f118_branch_tip_is_accepted_and_profile_self_assertion_is_rejected(campaign):
    validation = planner.validate_f118_current_source_authority(current_source_binding(campaign))
    assert validation["source_bundle"]["independent_validation"]["advertised_heads"] == [
        {
            "revision": campaign.controller_revision,
            "name": "refs/heads/feature/cgl-landau-fluid",
        }
    ]
    with pytest.raises(ValueError, match="keys differ"):
        planner.validate_f118_current_source_authority(
            {
                "evidence_sha256": planner.R03_F115_SHA256,
                "publication_audit_sha256": planner.R03_F115_PUBLICATION_AUDIT_SHA256,
                "provenance_review_sha256": (
                    planner.R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256
                ),
                "plasma_review_sha256": planner.R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
            }
        )
    authority_path = planner.f118_current_source_authority_path()
    authority = json.loads(authority_path.read_text())
    authority["authorization"]["prepare_authorized"] = True
    write_json(authority_path, authority)
    authority_path.chmod(0o444)
    rebound = republish_changed_f118_review_chain(campaign)
    with pytest.raises(ValueError, match="source-selection-only scope, or authority differs"):
        planner.validate_f118_current_source_authority(rebound)


def test_f118_requires_distinct_independent_reviews(campaign):
    provenance = json.loads(planner.f118_provenance_security_review_path().read_text())
    plasma_path = planner.f118_plasma_scientific_review_path()
    plasma = json.loads(plasma_path.read_text())
    plasma["reviewer"] = provenance["reviewer"]
    write_json(plasma_path, plasma)
    plasma_path.chmod(0o444)
    rebound = republish_changed_f118_review_chain(campaign)
    with pytest.raises(ValueError, match="distinct reviewers"):
        planner.validate_f118_current_source_authority(rebound)


def test_independent_review_is_declared_process_separation_not_crypto_identity(campaign):
    plan = campaign.plan()
    recost = load_utility(
        RECOST_UTILITY, "cgl_lf_stage_i_recost_review_assurance_invariant"
    )
    assert planner.REQUEST_REVIEW_IDENTITY_ASSURANCE == (
        recost.REQUEST_REVIEW_IDENTITY_ASSURANCE
    )
    assert planner.REQUEST_REVIEW_IDENTITY_LIMITATION == (
        recost.REQUEST_REVIEW_IDENTITY_LIMITATION
    )
    request_reviewer = plan["evidence"]["recost_authority"]["request_review_chain"][
        "independent_review"
    ]["reviewer"]
    assert request_reviewer["declared_process_independence"] is True
    assert request_reviewer["identity_assurance"] == (
        planner.REQUEST_REVIEW_IDENTITY_ASSURANCE
    )
    assert request_reviewer["identity_assurance_limitation"] == (
        planner.REQUEST_REVIEW_IDENTITY_LIMITATION
    )
    f116 = planner.validate_f118_current_source_authority(current_source_binding(campaign))
    assurances = (
        plan["evidence"]["recost_authority"]["independent_review_assurance"],
        plan["evidence"]["recost_authority"]["request_review_chain"][
            "independent_review_assurance"
        ],
        f116["independent_review_assurance"],
    )
    for assurance in assurances:
        assert assurance["basis"] == "declared-process-independence"
        assert assurance["strict_distinct_role_and_agent_declarations"] is True
        assert assurance["cryptographic_identity_verified"] is False
        assert len(set(assurance["declared_role_agents"].values())) == len(
            assurance["declared_role_agents"]
        )
        assert "do not cryptographically authenticate" in assurance[
            "non_cryptographic_limitation"
        ]
    with pytest.raises(ValueError, match="at least two distinct process roles"):
        planner.declared_process_independence_assurance(
            {"reviewer": "fixture reviewer"}, "fixture review"
        )
    campaign.request_reviewer_agent = "fixture recost requester"
    campaign.refresh()
    with pytest.raises(ValueError, match="roles or agents are not strictly distinct"):
        campaign.plan()
    campaign.request_reviewer_agent = "fixture-independent-request-reviewer"
    campaign.recost_reviewer_agent = "fixture recost requester"
    campaign.refresh()
    with pytest.raises(ValueError, match="roles or agents are not strictly distinct"):
        campaign.plan()


def test_f118_reviews_cannot_predate_supersession_evidence(campaign):
    review_path = planner.f118_provenance_security_review_path()
    review = json.loads(review_path.read_text())
    review["reviewed_utc"] = utc_text(NOW - timedelta(seconds=1))
    write_json(review_path, review)
    review_path.chmod(0o444)
    rebound = republish_changed_f118_review_chain(campaign)
    with pytest.raises(ValueError, match="identity, independence, or verification differs"):
        planner.validate_f118_current_source_authority(rebound)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("head", "not live repository HEAD"),
        ("helper", "committed/live tool bytes or mode differ"),
        ("bundle", "bundle digest changed"),
    ],
)
def test_f116_live_head_helper_and_final_bundle_fail_closed(campaign, mutation, message):
    if mutation == "head":
        git(campaign.helper.parents[2], "checkout", "--detach", planner.R03_F115_CONTROLLER_REVISION)
    elif mutation == "helper":
        write_text(campaign.helper, "tampered helper\n")
    else:
        write_text(campaign.current_bundle, "tampered final bundle\n")
    with pytest.raises(ValueError, match=message):
        campaign.plan()


def test_promoted_bundle_rejects_multiple_advertised_tips(campaign, tmp_path):
    bundle = tmp_path / "multiple-tips.bundle"
    git(campaign.helper.parents[2], "bundle", "create", str(bundle), "--all")
    with pytest.raises(ValueError, match="exactly one tip"):
        planner.validate_source_bundle_coverage(
            bundle,
            bundle_sha256=sha256(bundle),
            controller_revision=campaign.controller_revision,
            controller_sha256=campaign.controller_sha256,
            required_revisions=campaign.current_bundle_revisions,
            allow_branch_ref=True,
            label="fixture multiple-tip bundle",
        )


def test_historical_f115_bundle_rejects_branch_tip(campaign, tmp_path):
    repository = campaign.helper.parents[2]
    git(repository, "branch", "fixture-f115-branch", planner.R03_F115_CONTROLLER_REVISION)
    bundle = tmp_path / "f115-branch.bundle"
    git(repository, "bundle", "create", str(bundle), "fixture-f115-branch")
    with pytest.raises(ValueError, match="exact promoted controller tip"):
        planner.validate_source_bundle_coverage(
            bundle,
            bundle_sha256=sha256(bundle),
            controller_revision=planner.R03_F115_CONTROLLER_REVISION,
            controller_sha256=planner.R03_F115_CONTROLLER_SHA256,
            required_revisions=planner.F115_SOURCE_BUNDLE_REQUIRED_REVISIONS,
            allow_branch_ref=False,
            label="fixture historical F115 bundle",
        )


def test_recost_requires_canonical_latest_checkpoint_publication_chain(campaign, tmp_path):
    fake = tmp_path / campaign.recost_path.name
    write_bytes(fake, campaign.recost_path.read_bytes(), 0o444)
    with pytest.raises(ValueError, match="canonical schema"):
        planner.plan_from_paths(
            fake,
            planner.recost_independent_review_path(fake),
            planner.recost_publication_audit_path(fake),
        )
    competing = (
        campaign.root
        / "accounting/mks24_stage_i_E03_forcing_policy_F201_recost_evidence.json.publication_audit.json"
    )
    write_json(
        competing,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "published_utc": utc_text(),
        },
    )
    competing.chmod(0o444)
    with pytest.raises(ValueError, match="uniquely latest"):
        campaign.plan()


def test_minimal_fixture_only_recost_audit_is_rejected(campaign):
    _, review_path, audit_path = recost_paths(campaign)
    write_json(
        audit_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": planner.EXECUTION_EPOCH,
            "published_utc": utc_text(),
            "artifact": {
                "path": str(campaign.recost_path),
                "sha256": sha256(campaign.recost_path),
                "mode": "0444",
                "links": 1,
            },
            "independent_review": {
                "path": str(review_path),
                "sha256": sha256(review_path),
                "mode": "0444",
                "links": 1,
            },
            "authority": {
                "action_authority": False,
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            },
        },
    )
    audit_path.chmod(0o444)
    with pytest.raises(ValueError, match="keys differ"):
        planner.validate_recost_publication_chain(*recost_paths(campaign))


def test_r12_historical_clean_partial_is_inventory_only_and_never_completes(campaign):
    plan = campaign.plan()
    packet = next(item for item in plan["wave"]["packets"] if item["case_id"] == "R12")
    transition = plan["evidence"]["r12_fresh_rerun_transition"]
    assert packet["lineage"]["kind"] == "fresh"
    assert packet["lineage"]["continuation_provenance"] is None
    assert transition["historical_inventory"]["job_id"] == (
        planner.R12_HISTORICAL_PARTIAL_JOB_ID
    )
    assert transition["historical_inventory"]["segment"] == (
        planner.R12_HISTORICAL_PARTIAL_SEGMENT
    )
    assert transition["historical_inventory"]["classification"] == (
        "inventory-only-non-authorizing"
    )
    assert transition["continuation_authorized"] is False
    assert transition["waiver_authorized"] is False
    assert planner.lineage_complete(
        [
            {
                "manifest": {
                    "scientific_inspection": {"final_time": 10},
                    "accounting": {"result": "clean_partial"},
                }
            }
        ]
    ) is False


def test_fresh_r12_preserves_r03_r04_r16_continuations_and_ten_node_wave(campaign):
    r03_s00 = campaign.manifest_path("R03", segment(0, "0", "0.5"))
    r03 = campaign.add_recorded(
        "R03",
        planner.R03_F115_SEGMENT,
        parent_path=r03_s00,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    r04 = campaign.add_recorded(
        "R04",
        segment(0, "0", "0.25"),
        nodes=4,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    r16 = campaign.add_recorded(
        "R16",
        segment(0, "0", "1.5"),
        nodes=1,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    campaign.refresh()
    plan = campaign.plan()
    packets = {packet["case_id"]: packet for packet in plan["wave"]["packets"]}
    assert plan["wave"]["planned_nodes"] == 10
    assert packets["R03"]["lineage"]["continuation_provenance"]["parent_manifest"][
        "path"
    ] == str(r03)
    assert packets["R04"]["lineage"]["continuation_provenance"]["parent_manifest"][
        "path"
    ] == str(r04)
    assert packets["R16"]["lineage"]["continuation_provenance"]["parent_manifest"][
        "path"
    ] == str(r16)
    assert packets["R12"]["lineage"] == {
        "kind": "fresh",
        "segment": planner.R12_FRESH_RERUN_SEGMENT,
        "start_time": "0",
        "target_time": "0.12",
        "restart_required": False,
        "continuation_provenance": None,
    }


def test_r12_continues_normally_after_fresh_s01_is_recorded(campaign):
    fresh = campaign.add_recorded(
        "R12",
        planner.R12_FRESH_RERUN_SEGMENT,
        nodes=planner.R12_FRESH_RERUN_NODES,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    campaign.refresh()
    plan = campaign.plan()
    packet = next(item for item in plan["wave"]["packets"] if item["case_id"] == "R12")
    assert packet["lineage"]["kind"] == "continuation"
    assert packet["lineage"]["segment"] == segment(2, "0.12", "0.24")
    assert packet["lineage"]["continuation_provenance"]["parent_manifest"]["path"] == str(
        fresh
    )
    assert packet["allocation"]["nodes"] == planner.R12_FRESH_RERUN_NODES
    transition = plan["evidence"]["r12_fresh_rerun_transition"]
    assert transition["transition_state"] == "fresh-rerun-established"
    assert transition["fresh_rerun_profile"] is None
    assert transition["fresh_lineage_root"]["path"] == str(fresh)
    assert transition["waiver_authorized"] is False


def test_r12_continuation_increment_must_keep_natural_point_zero_two_alignment(campaign):
    campaign.add_recorded(
        "R12",
        planner.R12_FRESH_RERUN_SEGMENT,
        nodes=planner.R12_FRESH_RERUN_NODES,
        walltime=planner.R12_FRESH_RERUN_WALLTIME,
        athena_walltime=planner.R12_FRESH_RERUN_ATHENA_WALLTIME,
    )
    campaign.profiles["R12"]["next_increment"] = "0.13"
    campaign.refresh()
    with pytest.raises(ValueError, match="continuation recost next profile R12 lineage differs"):
        campaign.plan()


def test_non_r12_continuation_increment_retains_quarter_unit_alignment(campaign):
    campaign.add_recorded(
        "R04",
        segment(0, "0", "0.25"),
        nodes=4,
        walltime="02:00:00",
        athena_walltime="01:50:00",
    )
    campaign.profiles["R04"]["next_increment"] = "0.30"
    campaign.refresh()
    with pytest.raises(ValueError, match="continuation recost next profile R04 lineage differs"):
        campaign.plan()


@pytest.mark.parametrize(
    ("case_id", "segment_id"),
    (
        ("R03", "s00_rankio_t0_t0p5"),
        ("R12", "s00_rankio_t0_t0p25"),
    ),
)
def test_exact_retained_frozen_e03_no_max_ndiv_migration_is_authenticated(
    case_id, segment_id
):
    live = load_utility(
        UTILITY, f"cgl_lf_stage_i_wave_plan_live_migration_{case_id}"
    )
    recost = load_utility(
        RECOST_UTILITY, f"cgl_lf_stage_i_recost_live_migration_{case_id}"
    )
    assert live.FROZEN_E03_NO_MAX_NDIV_MIGRATIONS == (
        recost.FROZEN_E03_NO_MAX_NDIV_MIGRATIONS
    )
    assert live.FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY == (
        recost.FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY
    )
    assert live.FROZEN_E03_CT_DIVERGENCE_REASON == (
        recost.FROZEN_E03_CT_DIVERGENCE_REASON
    )
    path = live.manifest_expected_path(case_id, segment_id)
    if not path.is_file():
        pytest.skip(f"retained {case_id} frozen-E03 migration is unavailable")
    manifest, manifest_evidence = live.read_json_file(path, f"retained {case_id} manifest")
    index, start, target = live.parse_segment(segment_id, f"retained {case_id} segment")
    info = {
        "manifest": manifest,
        "evidence": manifest_evidence,
        "path": str(path),
        "case_id": case_id,
        "segment": segment_id,
        "index": index,
        "start": start,
        "target": target,
        "state": "recorded",
        "nodes": manifest["allocation"]["nodes"],
        "reserved": Decimal(str(manifest["allocation"]["reserved_node_hours"])),
    }
    inspection = manifest["scientific_inspection"]
    final = live.decimal_value(inspection["final_time"], f"retained {case_id} final time")
    retained = live.validate_clean_partial_continuation_evidence(info, inspection, final)
    available = retained["available_diagnostics"]
    plasma = available["plasma_continuation_evidence"]
    mhd = live.parse_controller_history(
        Path(inspection["mhd_history"]["path"]).read_bytes(),
        f"retained {case_id} MHD history",
    )
    user = live.parse_controller_history(
        Path(inspection["user_history"]["path"]).read_bytes(),
        f"retained {case_id} user history",
    )
    assert "max_ndiv" not in user
    with pytest.raises(ValueError, match="irreducibly unavailable.*max_ndiv"):
        live.continuation_plasma_evidence(case_id, mhd, user)
    recost_plasma = recost.require_frozen_e03_no_max_ndiv_migration(
        inspection,
        manifest,
        str(manifest["job_id"]),
        case_id,
        segment_id,
        int(info["nodes"]) * live.RANKS_PER_NODE,
        mhd,
        user,
    )
    live_controller = load_utility(
        CONTROLLER_UTILITY, f"cgl_lf_stage_i_controller_live_migration_{case_id}"
    )
    controller_manifest = json.loads(path.read_text())
    controller_plasma = live_controller.revalidate_continuation_plasma_evidence(
        controller_manifest["scientific_inspection"], controller_manifest
    )["plasma_continuation_evidence"]
    classification_keys = {"authorization_basis", "continuation_eligible"}
    for peer in (recost_plasma, controller_plasma):
        assert {
            key: value for key, value in plasma.items() if key not in classification_keys
        } == {
            key: value for key, value in peer.items() if key not in classification_keys
        }
    assert retained["policy"] == live.FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY
    assert retained["authorizing"] is False
    assert retained["authorization_effect"] == "none"
    assert retained["ct_divergence_claimed"] is False
    assert retained["ct_divergence_reason"] == live.FROZEN_E03_CT_DIVERGENCE_REASON
    assert retained["reason"] == live.FROZEN_E03_CT_DIVERGENCE_REASON
    assert retained["migration_eligibility"]["eligible"] is (case_id != "R12")
    assert retained["migration_eligibility"]["authorizing"] is False
    assert retained["migration_eligibility"]["authorization_effect"] == "none"
    assert retained["migration_eligibility"][
        "requires_exact_recommended_profile_waiver"
    ] is False
    assert available["mhd_history_columns"] == sorted(mhd)
    assert available["user_history_columns"] == sorted(user)
    assert plasma["continuation_authorized"] is False
    assert plasma["continuation_eligible"] is (case_id != "R12")
    assert plasma["eligibility_only"] is True
    assert plasma["schema_version"] == 2
    assert plasma["policy"] == live.FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY
    assert plasma["migration_contract"]["authority"] == {
        "continuation_authorized": False,
        "submission_authorized": False,
        "scheduler_mutation_authorized": False,
        "canonical_mutation_authorized": False,
    }
    assert plasma["normalized_ct_divb_evidence"]["authorizing"] is False
    assert plasma["normalized_ct_divb_evidence"]["reason"] == (
        live.FROZEN_E03_CT_DIVERGENCE_REASON
    )
    assert plasma["measurements"]["normalized_ct_divb_max"] is None
    assert plasma["migration_contract"] == retained["migration_eligibility"][
        "migration_contract"
    ]
    info["clean_partial_continuation_evidence"] = retained
    assert live.terminal_restart(info)["clean_partial_continuation_evidence"] == retained


@pytest.mark.parametrize("mutation", ("job", "input", "rank_count", "inspection", "legacy"))
def test_frozen_e03_no_max_ndiv_migration_rejects_identity_and_legacy_drift(
    mutation,
):
    live = load_utility(
        UTILITY, f"cgl_lf_stage_i_wave_plan_migration_adversary_{mutation}"
    )
    case_id = "R12"
    segment_id = "s00_rankio_t0_t0p25"
    path = live.manifest_expected_path(case_id, segment_id)
    if not path.is_file():
        pytest.skip("retained R12 frozen-E03 migration is unavailable")
    manifest, manifest_evidence = live.read_json_file(path, "retained R12 manifest")
    index, start, target = live.parse_segment(segment_id, "retained R12 segment")
    info = {
        "manifest": deepcopy(manifest),
        "evidence": manifest_evidence,
        "path": str(path),
        "case_id": case_id,
        "segment": segment_id,
        "index": index,
        "start": start,
        "target": target,
        "state": "recorded",
        "nodes": manifest["allocation"]["nodes"],
        "reserved": Decimal(str(manifest["allocation"]["reserved_node_hours"])),
    }
    inspection = deepcopy(manifest["scientific_inspection"])
    mhd = live.parse_controller_history(
        Path(inspection["mhd_history"]["path"]).read_bytes(), "retained R12 MHD history"
    )
    user = live.parse_controller_history(
        Path(inspection["user_history"]["path"]).read_bytes(), "retained R12 user history"
    )
    if mutation == "job":
        info["manifest"]["job_id"] = "9999999"
    elif mutation == "input":
        info["manifest"]["command"]["input_sha256"] = "0" * 64
    elif mutation == "rank_count":
        info["nodes"] = 2
    elif mutation == "inspection":
        inspection["final_hardwall_projection_count"] = 1.0
    else:
        inspection["plasma_continuation_policy"] = live.CONTINUATION_PLASMA_POLICY
    with pytest.raises(ValueError, match="frozen-E03"):
        live.require_frozen_e03_no_max_ndiv_migration(
            info,
            inspection,
            mhd,
            user,
            int(info["nodes"]) * live.RANKS_PER_NODE,
        )


def test_r12_fresh_rerun_transition_rejects_parent_restart_and_profile_drift():
    historical = {
        "segment": planner.R12_HISTORICAL_PARTIAL_SEGMENT,
        "evidence": {"path": "/manifest", "sha256": "5" * 64},
        "manifest": {
            "job_id": planner.R12_HISTORICAL_PARTIAL_JOB_ID,
            "accounting": {"result": "clean_partial"},
        },
    }
    state = {
        "historical_inventory": {"R12": [historical]},
        "lineages": {"R12": ([], None)},
    }
    profile = {
        "case_id": "R12",
        "segment": planner.R12_FRESH_RERUN_SEGMENT,
        "nodes": planner.R12_FRESH_RERUN_NODES,
        "ranks_per_node": planner.RANKS_PER_NODE,
        "walltime": planner.R12_FRESH_RERUN_WALLTIME,
        "athena_walltime": planner.R12_FRESH_RERUN_ATHENA_WALLTIME,
        "parent_job_id": None,
        "parent_result": None,
        "parent_segment": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
    }
    retained = planner.validate_r12_fresh_rerun_transition(state, {"R12": profile})
    assert retained["policy"] == planner.R12_FRESH_RERUN_POLICY
    assert retained["policy"].endswith("-v2")
    assert retained["continuation_authorized"] is False
    assert retained["waiver_authorized"] is False
    for key, value in (
        ("segment", segment(1, "0.1371931229426507", "0.25")),
        ("nodes", 2),
        ("ranks_per_node", 4),
        ("walltime", "01:00:00"),
        ("athena_walltime", "00:50:00"),
        ("parent_job_id", planner.R12_HISTORICAL_PARTIAL_JOB_ID),
        ("restart_file", "/historical/restart"),
    ):
        candidate = deepcopy(profile)
        candidate[key] = value
        with pytest.raises(ValueError, match="R12 fresh-rerun transition profile differs"):
            planner.validate_r12_fresh_rerun_transition(state, {"R12": candidate})


def test_recost_request_review_rejects_any_r12_waiver_scope(campaign):
    recost = json.loads(campaign.recost_path.read_text())
    audit = json.loads(planner.recost_publication_audit_path(campaign.recost_path).read_text())
    request_binding = audit["generalized_publication_context"]["request"]
    request_path = Path(request_binding["path"])
    review_path = request_path.with_name(f"{request_path.name}.independent_review.json")
    review = json.loads(review_path.read_text())
    review["scope"]["r12_frozen_e03_no_ct_profile_waiver"] = {
        "waiver_authorized": True
    }
    write_json(review_path, review)
    review_path.chmod(0o444)
    recost["provenance"]["request_independent_review_sha256"] = sha256(review_path)
    with pytest.raises(ValueError, match="recost request review"):
        planner.validate_recost_request_review_chain(request_binding, recost)


def test_request_review_scope_is_exact_one_key_across_production_tools(campaign):
    recost = load_utility(RECOST_UTILITY, "cgl_lf_stage_i_recost_for_scope_contract")
    checkpoint = load_utility(
        CHECKPOINT_UTILITY, "cgl_lf_stage_i_checkpoint_for_scope_contract"
    )
    scope = campaign.plan()["evidence"]["recost_authority"]["request_review_chain"][
        "independent_review"
    ]["scope"]
    assert scope == {"non_authorizing": True}
    for function in (
        planner.validate_recost_request_review_chain,
        recost.parse_request_independent_review,
        checkpoint.validate_schema2_independent_review,
        controller.latest_published_r17_recost,
    ):
        literals = []
        for node in ast.walk(ast.parse(inspect.getsource(function))):
            if not isinstance(node, ast.Dict):
                continue
            try:
                value = ast.literal_eval(node)
            except (ValueError, TypeError):
                continue
            if isinstance(value, dict) and "non_authorizing" in value:
                literals.append(value)
        assert literals
        assert all(value == {"non_authorizing": True} for value in literals)


def test_wave_planner_contract_contains_no_r12_waiver_planning_path():
    source = UTILITY.read_text()
    assert "r12_frozen_e03_no_ct_profile_waiver" not in source
    assert "continuation_eligibility_waiver" not in source
    assert "validate_r12_frozen_e03_no_ct_profile_waiver" not in source
    assert planner.R12_FRESH_RERUN_SEGMENT == "s01_rankio_t0_t0p12"
    assert planner.R12_FRESH_RERUN_TARGET == Decimal("0.12")
    assert planner.R12_FRESH_RERUN_NODES == 4
    assert planner.R12_FRESH_RERUN_RANKS == 32
    assert planner.R12_FRESH_RERUN_WALLTIME == "02:00:00"
    assert planner.R12_FRESH_RERUN_ATHENA_WALLTIME == "01:50:00"
    assert planner.R12_CONTINUATION_ALIGNMENT == Decimal("0.02")
    assert planner.STANDARD_CONTINUATION_ALIGNMENT == Decimal("0.25")


def test_fresh_r12_identity_is_exact_across_production_tool_contracts():
    recost = load_utility(RECOST_UTILITY, "cgl_lf_stage_i_recost_for_r12_contract")
    checkpoint = load_utility(
        CHECKPOINT_UTILITY, "cgl_lf_stage_i_checkpoint_for_r12_contract"
    )
    assert planner.R12_FRESH_RERUN_SEGMENT == controller.R12_FRESH_RERUN_SEGMENT
    assert planner.R12_FRESH_RERUN_SEGMENT == recost.R12_FRESH_RERUN_SEGMENT
    assert planner.R12_FRESH_RERUN_SEGMENT == checkpoint.FRESH_R12_RERUN["segment"]
    assert planner.R12_FRESH_RERUN_NODES == controller.R12_FRESH_RERUN_NODES
    assert planner.R12_FRESH_RERUN_NODES == checkpoint.FRESH_R12_RERUN["nodes"]
    assert checkpoint.FRESH_R12_RERUN["time_tlim_target"] == 0.12
    assert planner.R12_HISTORICAL_PARTIAL_JOB_ID == recost.R12_HISTORICAL_JOB_ID
    assert planner.R12_HISTORICAL_PARTIAL_JOB_ID == (
        controller.R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID
    )
    assert planner.R12_HISTORICAL_PARTIAL_JOB_ID == (
        checkpoint.HISTORICAL_R12_CLEAN_PARTIAL["job_id"]
    )
    assert planner.R12_HISTORICAL_PARTIAL_SEGMENT == recost.R12_HISTORICAL_SEGMENT
    assert planner.R12_HISTORICAL_PARTIAL_SEGMENT == (
        checkpoint.HISTORICAL_R12_CLEAN_PARTIAL["segment"]
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "legacy",
        "required_reached",
        "strict_failure",
        "accepted",
        "restart_marker_bypass",
        "marker_modes",
        "snapshot_endpoint",
        "plasma_check",
        "plasma_policy",
        "plasma_evidence",
        "history_binding",
        "history_drift",
        "product_binding",
        "retained_drift",
    ],
)
def test_clean_partial_requires_strengthened_physical_continuation_evidence(
    campaign, mutation
):
    parent = campaign.add_recorded(
        "R13",
        segment(0, "0", "0.25"),
        final_time="0.125",
        result="clean_partial",
        nodes=2,
    )
    inspection = campaign.manifests[parent]["scientific_inspection"]
    inspection_path = parent.parent / "segment_inspection.json"
    if mutation == "legacy":
        inspection = {
            key: inspection[key]
            for key in (
                "required_time",
                "final_time",
                "terminal_restart_time",
                "terminal_restart",
                "clean_for_continuation",
            )
        }
        campaign.manifests[parent]["scientific_inspection"] = inspection
    elif mutation == "required_reached":
        inspection["checks"]["required_time_reached"] = True
    elif mutation == "strict_failure":
        inspection["maximum_strict_failure_counts"]["lf_nonfin"] = 1.0
    elif mutation == "accepted":
        inspection["accepted"] = True
    elif mutation == "restart_marker_bypass":
        inspection["restart_time_marker_bypass"] = True
    elif mutation == "marker_modes":
        inspection["restart_time_marker_modes"][0].pop()
    elif mutation == "snapshot_endpoint":
        inspection["snapshot_times"][-1] = 0.124
    elif mutation == "plasma_check":
        inspection["checks"]["plasma_continuation_policy"] = False
    elif mutation == "plasma_policy":
        inspection["plasma_continuation_policy"] = "stage-i-clean-partial-continuation-v1"
    elif mutation == "plasma_evidence":
        inspection["plasma_continuation_evidence"]["continuation_authorized"] = False
    elif mutation == "history_binding":
        inspection["mhd_history"]["sha256"] = "0" * 64
    elif mutation == "product_binding":
        inspection["snapshots"][0]["rank_files"][0]["sha256"] = "0" * 64
    write_json(inspection_path, inspection)
    campaign.refresh()
    if mutation == "history_drift":
        write_text(Path(inspection["mhd_history"]["path"]), "drifted history\n")
    elif mutation == "retained_drift":
        retained = json.loads(inspection_path.read_text())
        retained["clean_for_continuation"] = False
        write_json(inspection_path, retained)
    with pytest.raises(ValueError, match="clean_partial continuation evidence"):
        campaign.plan()


def test_canonical_state_mutation_after_recost_publication_is_rejected(campaign):
    write_json(planner.reservations_path(), [])
    with pytest.raises(ValueError, match="recost reservations binding differs"):
        campaign.plan()


@pytest.mark.parametrize("mutation", ["budget", "profile"])
def test_recost_projection_and_profile_tampering_fail_closed(campaign, mutation):
    recost = json.loads(campaign.recost_path.read_text())
    if mutation == "budget":
        recost["budget"]["computed_stage_i_total_node_hours"] = "0.000001"
    else:
        recost["recommendations"]["recommended_next_profiles"][1]["source_bundle_sha256"] = (
            planner.R03_F115_SOURCE_BUNDLE_SHA256
        )
    write_json(campaign.recost_path, recost)
    campaign.recost_path.chmod(0o444)
    with pytest.raises(
        ValueError,
        match="exact immutable published artifact|binding differs|independent review differs",
    ):
        campaign.plan()


def test_strong_recost_r17_chain_is_sole_exclusive_readiness_path(campaign):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.refresh()
    plan = campaign.plan()
    assert plan["wave"]["status"] == "r17_exclusive"
    assert plan["wave"]["planned_nodes"] == 8
    assert [packet["case_id"] for packet in plan["wave"]["packets"]] == ["R17"]
    assert plan["evidence"]["r17_readiness"]["artifact"]["path"] == str(
        planner.r17_readiness_path()
    )
    recost = json.loads(campaign.recost_path.read_text())
    readiness = json.loads(planner.r17_readiness_path().read_text())
    assert readiness["current_source_authority"] == recost["provenance"]["source_authority"]
    assert set(plan["evidence"]["r17_readiness"]) >= {
        "scheduler_evidence",
        "account_scheduler_evidence",
        "account_exclusivity_evidence",
    }


def test_planner_consumes_real_schema2_r17_qualification_producer_output(
    tmp_path, monkeypatch
):
    suite = load_utility(
        QUALIFICATION_TEST_UTILITY,
        "cgl_lf_stage_i_qualification_suite_for_wave_integration",
    )
    producer = suite.qualification
    fixture = suite.qualification_fixture.__wrapped__(tmp_path, monkeypatch)
    monkeypatch.setattr(producer, "RANKS_PER_NODE", 8)
    _, _, evidence, scheduler = suite.retain_complete_wave(
        fixture, ("R17",), max_nodes=8
    )
    retained = producer.retain_r17_operational_qualification(
        evidence["R17"], "producer measurement agent", scheduler
    )
    root = fixture["root"]
    qualification_path = root / retained["operational_qualification"]["path"]
    qualification = producer.load_json(
        qualification_path, "producer R17 operational qualification"
    )
    measured = datetime.fromisoformat(qualification["measured_utc"])
    reviewed = measured + timedelta(seconds=1)
    review_path = root / qualification["independent_review_contract"]["path"]
    suite.write_json(
        review_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-operational-qualification-independent-review",
            "execution_epoch": producer.EXECUTION_EPOCH,
            "reviewed_utc": reviewed.isoformat(),
            "decision": "approved",
            "reviewer": "independent producer-output reviewer",
            "candidate": {
                "path": str(qualification_path),
                "sha256": producer.sha256(qualification_path),
            },
        },
    )
    review_path.chmod(0o444)

    expected_cases = deepcopy(planner.EXPECTED_CASES)
    r17 = expected_cases["R17"]
    expected_cases["R17"] = (
        r17[0],
        r17[1],
        r17[2],
        producer.CASE_POLICIES["R17"]["input_sha256"],
    )
    monkeypatch.setattr(planner, "CANONICAL_ROOT", root)
    monkeypatch.setattr(planner, "SOURCE_REVISION", producer.FROZEN_SOURCE_REVISION)
    monkeypatch.setattr(planner, "MATRIX_SHA256", producer.FROZEN_MATRIX_SHA256)
    monkeypatch.setattr(planner, "EXPECTED_CASES", expected_cases)
    monkeypatch.setattr(planner, "build_manifest_path", lambda: fixture["build_manifest"])
    monkeypatch.setattr(
        planner,
        "r17_parameter_contract",
        lambda: qualification["frozen_science_build_contract"]["parameter_contract"],
    )
    monkeypatch.setattr(planner, "current_utc", lambda: reviewed)
    qualification_value, qualification_profile = planner.read_json_file(
        qualification_path, "producer R17 qualification"
    )
    review, _ = planner.read_json_file(review_path, "producer R17 review")
    validation, profiles, authors = planner.validate_r17_operational_qualification_contract(
        qualification_value,
        {
            "executable_sha256": qualification["executable_sha256"],
            "build_manifest_sha256": qualification["build_manifest_inventory_sha256"],
        },
        qualification_path,
        qualification_profile,
        review,
        reviewed,
    )
    exclusivity, account_profiles = planner.validate_r17_account_exclusivity_artifacts(
        qualification_value, reviewed
    )
    assert set(validation) == {
        "restart_load_evidence",
        "physics_validation_evidence",
    }
    assert authors == {"producer measurement agent"}
    assert set(profiles) == {
        "prepared_wave",
        "qualification_evidence",
        "scientific_evidence",
        "restart_load_evidence",
        "physics_validation_evidence",
    }
    assert exclusivity["exclusive_entire_execution_interval"] is True
    assert all(profile["mode"] == "0444" for profile in account_profiles.values())
    for key in (
        "prepared_wave",
        "qualification_evidence",
        "scientific_evidence",
        "restart_load_evidence",
        "physics_validation_evidence",
    ):
        candidate = deepcopy(qualification_value)
        candidate[key]["sha256"] = "0" * 64
        with pytest.raises(ValueError):
            planner.validate_r17_operational_qualification_contract(
                candidate,
                {
                    "executable_sha256": qualification["executable_sha256"],
                    "build_manifest_sha256": qualification[
                        "build_manifest_inventory_sha256"
                    ],
                },
                qualification_path,
                qualification_profile,
                review,
                reviewed,
            )
    for key in (
        "build_manifest_inventory_sha256",
        "rank_local_output_inventory_sha256",
        "rank_local_restart_inventory_sha256",
    ):
        candidate = deepcopy(qualification_value)
        candidate[key] = "0" * 64
        with pytest.raises(ValueError):
            planner.validate_r17_operational_qualification_contract(
                candidate,
                {
                    "executable_sha256": qualification["executable_sha256"],
                    "build_manifest_sha256": qualification[
                        "build_manifest_inventory_sha256"
                    ],
                },
                qualification_path,
                qualification_profile,
                review,
                reviewed,
            )
    for mutate in (
        lambda value: value["authority"].__setitem__("r17_launch_authorized", True),
        lambda value: value["frozen_science_build_contract"].__setitem__(
            "input_sha256", "0" * 64
        ),
        lambda value: value["independent_review_contract"].__setitem__(
            "candidate_path", str(qualification_path) + ".drift"
        ),
        lambda value: value.__setitem__("legacy_qualification_fallback", {}),
    ):
        candidate = deepcopy(qualification_value)
        mutate(candidate)
        with pytest.raises(ValueError):
            planner.validate_r17_operational_qualification_contract(
                candidate,
                {
                    "executable_sha256": qualification["executable_sha256"],
                    "build_manifest_sha256": qualification[
                        "build_manifest_inventory_sha256"
                    ],
                },
                qualification_path,
                qualification_profile,
                review,
                reviewed,
            )


@pytest.mark.parametrize("mutation", ["extra", "duplicate_rank"])
def test_r17_readiness_requires_exactly_one_terminal_output_per_rank(campaign, mutation):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    if mutation == "extra":
        campaign.r17_output_count = 65
    else:
        campaign.r17_duplicate_output_rank = True
    campaign.refresh()
    with pytest.raises(ValueError, match="exactly one file for each of 64 ranks"):
        campaign.plan()


@pytest.mark.parametrize(
    "mutation", ["logical_count", "per_rank", "duplicate", "digest", "output_binding"]
)
def test_r17_readiness_rejects_invalid_meshblock_decomposition_proof(campaign, mutation):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.r17_decomposition_mutation = mutation
    campaign.refresh()
    with pytest.raises(ValueError, match="R17 meshblock decomposition proof"):
        campaign.plan()


@pytest.mark.parametrize(
    "mutation",
    (
        "rank_count",
        "mass_drift",
        "mass_mismatch",
        "lf_bad_counts",
        "ct_threshold",
        "ct_failure",
    ),
)
def test_r17_readiness_rejects_coherent_physics_invalid_qualification(campaign, mutation):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.r17_physics_mutation = mutation
    campaign.refresh()
    with pytest.raises(ValueError, match="R17 physics measurements or scientific checks"):
        campaign.plan()


@pytest.mark.parametrize(
    "mutation",
    [
        "overlap",
        "private_data",
        "query_window",
        "path_scheduler",
        "path_exclusivity",
        "digest_scheduler",
        "digest_exclusivity",
        "target_field",
        "expanded_drift",
    ],
)
def test_r17_readiness_authenticates_account_wide_execution_exclusivity(
    campaign, mutation
):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.r17_account_mutation = mutation
    campaign.refresh()
    with pytest.raises(ValueError, match="R17 account"):
        campaign.plan()


def test_r17_readiness_must_bind_authenticated_f116_authority(campaign):
    campaign.complete_predecessors()
    campaign.r17_ready = True
    campaign.r17_source_authority_override = {"checkpoint": "profile-only-self-assertion"}
    campaign.refresh()
    with pytest.raises(ValueError, match="R17 readiness reviewed bindings"):
        campaign.plan()


def test_f116_required_tool_set_is_cross_tool_exact_and_includes_qualification():
    source_authority = load_utility(
        SOURCE_AUTHORITY_UTILITY,
        "cgl_lf_stage_i_source_authority_for_wave_invariant",
    )
    assert planner.F116_REQUIRED_TOOLS == controller.F116_REQUIRED_TOOLS
    assert planner.F116_REQUIRED_TOOLS == source_authority.REQUIRED_TOOLS
    assert planner.F116_REQUIRED_TOOLS[
        "scripts/frontier/cgl_lf_stage_i_qualification.py"
    ] == "0755"


def test_git_and_scheduler_child_environments_strip_caller_execution_controls(
    monkeypatch,
):
    hostile = {
        "PATH": "/tmp/hostile",
        "SBATCH_EXPORT": "ALL",
        "SLURM_CONF": "/tmp/slurm.conf",
        "LD_PRELOAD": "/tmp/preload.so",
        "DYLD_INSERT_LIBRARIES": "/tmp/dyld.so",
        "PYTHONPATH": "/tmp/python",
        "PERL5OPT": "-Mhostile",
        "RUBYOPT": "-rhostile",
        "GIT_CONFIG_GLOBAL": "/tmp/gitconfig",
        "BASH_ENV": "/tmp/bashenv",
        "BASH_FUNC_hostile%%": "() { :; }",
        "_CGL_LF_PRIVATE_DESCRIPTOR": "9",
    }
    for key, value in hostile.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setenv("TZ", "UTC")
    child = planner.hardened_child_environment()
    assert child == {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": planner.TRUSTED_SYSTEM_PATH,
        "XDG_CONFIG_HOME": "/nonexistent",
    }
    git_environment = planner.hardened_git_environment()
    assert git_environment["PATH"] == planner.TRUSTED_SYSTEM_PATH
    assert git_environment["GIT_EXEC_PATH"] == str(planner.GIT_EXEC_PATH)
    assert git_environment["HOME"] == "/nonexistent"
    assert not any(value in git_environment.values() for value in hostile.values())


def test_every_planner_child_process_is_env_hardened_and_not_caller_python():
    tree = ast.parse(UTILITY.read_text())
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "subprocess"
        and node.func.attr == "run"
    ]
    assert calls
    for call in calls:
        assert any(keyword.arg == "env" for keyword in call.keywords)
        assert not any(
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id == "sys"
            and node.attr == "executable"
            for argument in call.args
            for node in ast.walk(argument)
        )
        assert not any(
            isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and "python" in node.value.casefold()
            for argument in call.args
            for node in ast.walk(argument)
        )


def test_git_and_scheduler_subprocesses_receive_only_hardened_environments(
    monkeypatch,
):
    calls = []
    owner = planner.scheduler_owner()

    def runner(argv, **kwargs):
        calls.append((list(argv), kwargs))
        if argv[0] == str(planner.SQUEUE):
            return subprocess.CompletedProcess(
                argv,
                0,
                stdout=f"123|fixture|PENDING|1|{owner}|{planner.ACCOUNT}\n",
                stderr="",
            )
        if argv[0] == str(planner.SCONTROL):
            return subprocess.CompletedProcess(argv, 0, stdout=b"#!/bin/bash\n", stderr=b"")
        return subprocess.CompletedProcess(argv, 0, stdout=b"", stderr=b"")

    monkeypatch.setattr(planner.subprocess, "run", runner)
    planner.collect_live_scheduler_job("123")
    planner.collect_live_batch_script("123")
    planner.git_read_only(Path("/tmp"), ["rev-parse", "HEAD"])
    assert calls[0][1]["env"] == planner.hardened_child_environment()
    assert calls[1][1]["env"] == planner.hardened_child_environment()
    assert calls[2][0][0] == str(planner.GIT)
    assert calls[2][1]["env"] == planner.hardened_git_environment()


def test_cli_accepts_only_recost_chain_and_output_is_command_free(campaign, capsys):
    with pytest.raises(SystemExit):
        planner.parser().parse_args(
            [
                "--matrix", "retired.json",
                "--allocation-profile", "retired.json",
                "--summary", "retired.md",
                "--reconciliation", "retired.json",
            ]
        )
    capsys.readouterr()
    result = planner.main(
        [
            "--recost-artifact", str(campaign.recost_path),
            "--recost-independent-review",
            str(planner.recost_independent_review_path(campaign.recost_path)),
            "--recost-publication-audit",
            str(planner.recost_publication_audit_path(campaign.recost_path)),
        ]
    )
    captured = capsys.readouterr()
    assert result == 0
    assert captured.err == ""
    assert json.loads(captured.out)["read_only"] is True
    assert "sbatch " not in captured.out


def test_global_toctou_and_plan_digest_are_deterministic(campaign, monkeypatch):
    first = campaign.plan()
    second = campaign.plan()
    assert first == second
    core = {key: value for key, value in first.items() if key != "plan_sha256"}
    assert first["plan_sha256"] == planner.sha256_bytes(planner.canonical_json(core))
    original = planner.revalidate_global_boundary

    def mutate_then_revalidate(value, transaction_state, submitted_scheduler):
        write_json(planner.reservations_path(), [])
        return original(value, transaction_state, submitted_scheduler)

    monkeypatch.setattr(planner, "revalidate_global_boundary", mutate_then_revalidate)
    with pytest.raises(ValueError, match="global evidence changed before plan completion"):
        campaign.plan()
