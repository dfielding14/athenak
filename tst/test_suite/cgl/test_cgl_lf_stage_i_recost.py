"""Adversarial local-fixture regressions for the Stage I recost generator."""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timedelta, timezone
from decimal import Decimal
import hashlib
import importlib.util
import json
import math
import os
import fcntl
from pathlib import Path
import re
import stat
import subprocess
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
RECOST = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_recost.py"
STAGE_I = REPOSITORY / "scripts/frontier/cgl_lf_stage_i.py"
QUALIFICATION = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_qualification.py"
QUALIFICATION_TEST = REPOSITORY / "tst/test_suite/cgl/test_cgl_lf_stage_i_qualification.py"
EPOCH = "E03-forcing-policy"
EPOCH_SLUG = "E03_forcing_policy"
F113_NAME = "mks24_stage_i_E03_forcing_policy_F113_controller_transition_evidence.json"
F115_NAME = "mks24_stage_i_E03_forcing_policy_F115_source_bundle_recovery_supersession_evidence.json"
F116_NAME = "mks24_stage_i_E03_forcing_policy_F116_current_source_authority_supersession_evidence.json"
F118_NAME = "mks24_stage_i_E03_forcing_policy_F118_current_source_authority_supersession_evidence.json"
QUALIFICATION_NAME = "mks24_stage_i_E03_forcing_policy_qualification_approval.json"
R17_READINESS_NAME = "mks24_stage_i_E03_forcing_policy_R17_readiness_evidence.json"
AUTHORIZED_REVIEWER = (
    "Codex execution agent with independent Turing, Hubble, Darwin, and Linnaeus audits"
)
F113_REVIEW_AUTHORITY = "campaign-authorized-independent-review"
STORAGE_PROJECTION_METHOD = "observed-stage-i-output-byte-rate-v1"
REQUEST_REVIEW_IDENTITY_ASSURANCE = "declared-process-independence-non-cryptographic"
REQUEST_REVIEW_IDENTITY_LIMITATION = (
    "Reviewer identity and process independence are declared evidence, "
    "not cryptographically proven."
)
F116_REQUIRED_COMMITTED_TOOLS = {
    "scripts/frontier/cgl_lf_stage_i.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_checkpoint.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_qualification.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_recost.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_source_authority.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_validate_segment.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_wave_plan.py": "0644",
}
F116_AUTHORIZATION = {
    "current_source_selection_authorized": True,
    "source_authority_publication_authorized": True,
    "prepare_authorized": False,
    "submit_authorized": False,
    "direct_sbatch_authorized": False,
    "scheduler_mutation_authorized": False,
    "stage_i_execution_state_mutation_authorized": False,
    "scientific_configuration_change_authorized": False,
    "historical_manifest_rebinding_authorized": False,
}
F118_SCOPE_PRESERVES = [
    "The immutable F-116 evidence, reviews, publication audit, selected source bundle, and nested F-115 historical authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
F118_SCOPE_DOES_NOT_AUTHORIZE = [
    "prepare",
    "submit",
    "direct sbatch",
    "scheduler mutation",
    "Stage I execution-state mutation",
    "scientific configuration change",
    "historical manifest rebinding",
]
F118_VALIDATION_CLAIMS = {
    "historical_f116_chain": "passed",
    "historical_f115_chain": "passed",
    "bridge_bundle_complete_history": "passed",
    "predecessor_current_source_bundle_complete_history": "passed",
    "final_bundle_complete_history": "passed",
    "final_bundle_single_head_tip": "passed",
    "final_bundle_required_revisions": "passed",
    "committed_tool_bytes": "passed",
    "corrupt_c7_exclusion_preserved": True,
}
F118_PUBLICATION_REQUIREMENTS = {
    "published_evidence_mode": "0444",
    "published_review_mode": "0444",
    "published_audit_mode": "0444",
    "published_links": 1,
    "publication_audit_is_authority_commit_marker": True,
    "recovery_required_after_interruption": True,
}
F118_PUBLICATION_METHOD = (
    "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
)
F119_REQUESTED_BY = "Codex deterministic F119 draft packet renderer"
F119_SCOPE = (
    "Non-authorizing F119 recost request draft superseding the exact authenticated "
    "failed F117 attempt under current F118 source authority."
)
INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION = (
    "Reviewer roles, agent identifiers, and process separation are retained "
    "declarations; exact artifact digests authenticate reviewed bytes but do not "
    "cryptographically authenticate a human or agent identity."
)
POLICIES = {
    **{case_id: "U+A+H" for case_id in (
        "R02", "R03", "R04", "R05", "R10", "R11", "R12", "R13", "R16", "R17"
    )},
    **{case_id: "U+P+H" for case_id in ("R06", "R07", "R08", "R09")},
    **{case_id: "U+A+F" for case_id in ("R14", "R15")},
}
POLICY_TEXT = {
    "U": (
        "Require the exact endpoint; finite synchronized MHD and user histories; "
        "relative mass drift and MHD/user mass mismatch <= 1e-12; zero lf_dfloor, "
        "lf_pfloor, lf_nonfin, lf_nonpos, lf_hardbd, and hard_vol at every retained "
        "row; positive interval lf_nstage and lf_qface with bounded cap increments; "
        "finite LF heat and pressure-work ledgers; complete ranked products; and one "
        "loadable exact terminal snapshot and restart group. Frozen E03 does not emit "
        "an independently auditable CT-divergence diagnostic, so this recost makes no "
        "CT-divergence bound claim."
    ),
    "A": (
        "For active CGL, require segment and whole-case "
        "|Delta E - Delta force_work| / scale < 1e-8 and finite nonzero applied "
        "pressure work over the developed window."
    ),
    "P": (
        "For passive Delta, do not apply active-CGL total-energy closure; require "
        "lf_cpwrk == lf_cawrk == 0 throughout and retain finite forcing work."
    ),
    "H": (
        "For hardwall cases, require nonnegative monotonic lf_hwproj and review "
        "developed-turbulence activity; zero activity requires physics review rather "
        "than automatic acceptance."
    ),
    "F": (
        "For finite-limiter cases, require lf_hwproj == 0 and finite threshold "
        "occupancy and nu_eff that distinguish the configured limiter rate."
    ),
}
INITIAL_TARGETS = {
    **{f"R{number:02d}": 0.25 for number in range(3, 6)},
    "R06": 0.5,
    **{f"R{number:02d}": 0.25 for number in range(7, 12)},
    "R12": 0.12,
    **{f"R{number:02d}": 0.25 for number in range(13, 16)},
    "R16": 1.5,
    "R17": 0.25,
}
SOLE_COMPATIBLE_KEYS = {
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
}


def sha256(path: Path) -> str:
    """Return one fixture digest."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def mutated_history_copy(source: Path, target: Path, column: str, value: float) -> dict[str, object]:
    """Copy one history while coherently changing a named final-row measurement."""

    lines = source.read_text().splitlines()
    labels = [token.split("=", 1)[1] for token in lines[1].split()[1:]]
    fields = lines[-1].split()
    fields[labels.index(column)] = f"{value:.16e}"
    target.write_text("\n".join([*lines[:-1], " ".join(fields)]) + "\n")
    target.chmod(0o644)
    return {
        "path": str(target),
        "size_bytes": target.stat().st_size,
        "sha256": sha256(target),
    }


def write_json(path: Path, value: object) -> None:
    """Write deterministic fixture JSON."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.chmod(0o644)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    path.chmod(0o644)


def write_immutable_json(path: Path, value: object) -> None:
    """Write deterministic fixture JSON as one immutable publication file."""

    if path.exists():
        path.chmod(0o644)
    write_json(path, value)
    path.chmod(0o444)


def git(repository: Path, *arguments: str) -> str:
    """Run one fixture Git command and return stdout."""

    completed = subprocess.run(
        ["git", *arguments],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def f116_committed_tools(repository: Path, revision: str) -> list[dict[str, object]]:
    """Return the exact source-authority seven-tool fixture contract."""

    return [
        {
            "path": relative,
            "revision": revision,
            "sha256": sha256(repository / relative),
            "mode": mode,
        }
        for relative, mode in sorted(F116_REQUIRED_COMMITTED_TOOLS.items())
    ]


def scheduler_time(value: datetime) -> str:
    """Return the retained naive-UTC sacct timestamp representation."""

    return value.astimezone(timezone.utc).replace(tzinfo=None).isoformat(timespec="seconds")


def acceptance_criterion(case_id: str, target: float) -> str:
    """Return exact case-aware acceptance prose."""

    policy = POLICIES[case_id]
    target_text = str(target).rstrip("0").rstrip(".") if "." in str(target) else str(target)
    clauses = " ".join(f"{code}: {POLICY_TEXT[code]}" for code in policy.split("+"))
    return (
        f"Accept {case_id} only at exact t={target_text} under policies {policy}. "
        f"{clauses} Any endpoint, provenance, product-inventory, scheduler-accounting, "
        "or policy failure blocks acceptance and successor packet planning."
    )


def directory_inventory_sha256(path: Path) -> str:
    """Return the generator's flat build-manifest inventory digest."""

    entries = [
        {
            "name": item.name,
            "mode": f"{stat.S_IMODE(item.stat().st_mode):04o}",
            "sha256": sha256(item),
        }
        for item in sorted(path.iterdir())
    ]
    return hashlib.sha256((json.dumps(entries, sort_keys=True) + "\n").encode()).hexdigest()


def retained_bytes(root: Path) -> int:
    """Return exact regular-file bytes in the retained E03 run tree."""

    run_store = root / "runs/mks24-stage-i" / EPOCH
    return sum(path.stat().st_size for path in run_store.rglob("*") if path.is_file())


def live_available(root: Path) -> int:
    """Return live available filesystem bytes."""

    value = os.statvfs(root)
    return value.f_bavail * value.f_frsize


def terminal_restart(output_dir: Path) -> dict[str, object]:
    """Create and return one complete eight-rank terminal restart group."""

    rank_files = []
    for rank in range(8):
        path = output_dir / "rst" / f"rank_{rank:08d}" / "terminal.rst"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"fixture restart rank {rank}\n".encode())
        path.chmod(0o644)
        rank_files.append(
            {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}
        )
    first = rank_files[0]
    return {
        "path": first["path"],
        "sha256": first["sha256"],
        "size_bytes": first["size_bytes"],
        "storage": "per_rank",
        "rank_files": rank_files,
    }


def terminal_snapshot(output_dir: Path) -> dict[str, object]:
    """Create and return one complete eight-rank terminal snapshot group."""

    rank_files = []
    for rank in range(8):
        path = output_dir / "bin" / f"rank_{rank:08d}" / "terminal.bin"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"fixture snapshot rank {rank}\n".encode())
        path.chmod(0o644)
        rank_files.append(
            {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}
        )
    first = rank_files[0]
    return {
        "path": first["path"],
        "sha256": first["sha256"],
        "size_bytes": first["size_bytes"],
        "storage": "per_rank",
        "rank_files": rank_files,
    }


def history_text(labels: list[str], rows: list[list[float]]) -> str:
    """Render one minimal labeled Athena history file."""

    header = "# " + " ".join(
        f"[{index}]= {label}".replace("= ", "=")
        for index, label in enumerate(labels, start=1)
    )
    return header + "\n" + "\n".join(
        " ".join(format(value, ".17g") for value in row) for row in rows
    ) + "\n"


def retained_file(path: Path) -> dict[str, object]:
    """Return the controller's retained-file schema for a fixture file."""

    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def load_stage_i_module():
    """Load the current production controller for cross-contract fixtures."""

    name = f"_cgl_lf_stage_i_recost_controller_{os.getpid()}_{id(object())}"
    spec = importlib.util.spec_from_file_location(name, STAGE_I)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_qualification_test_module():
    """Load the producer's own integration-fixture module."""

    name = f"_cgl_lf_stage_i_recost_qualification_tests_{os.getpid()}_{id(object())}"
    spec = importlib.util.spec_from_file_location(name, QUALIFICATION_TEST)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


stage_i = load_stage_i_module()


def controller_clean_partial_inspection(
    manifest: Path,
    output_dir: Path,
    *,
    case_id: str,
    segment: str,
    job_id: str,
    final_time: float,
    terminal_restart_value: dict[str, object],
) -> dict[str, object]:
    """Build a full schema-4 inspection using the current controller plasma policy."""

    mhd_labels = [
        "time", "mass", "tot-E", "lf_nstage", "lf_qface", "lf_qprcap",
        "lf_qpr10", "lf_qpecap", "lf_qpe10", "lf_qprwrk", "lf_qpewrk",
        "lf_cpwrk", "lf_cawrk", "lf_hwproj", "lf_dfloor", "lf_pfloor",
        "lf_nonfin", "lf_nonpos", "lf_hardbd",
    ]
    mhd_rows = [
        [0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [final_time, 2.0, 11.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.2,
         1.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]
    user_labels = ["time", "mass", "hard_vol", "force_work", "max_ndiv"]
    user_rows = [
        [0.0, 2.0, 0.0, 0.0, 1.0e-14],
        [final_time, 2.0, 0.0, 1.0, 1.0e-14],
    ]
    mhd_path = output_dir / "fixture.mhd.hst"
    user_path = output_dir / "fixture.user.hst"
    mhd_path.write_text(history_text(mhd_labels, mhd_rows))
    user_path.write_text(history_text(user_labels, user_rows))
    mhd_path.chmod(0o644)
    user_path.chmod(0o644)
    mhd = {
        label: [row[index] for row in mhd_rows]
        for index, label in enumerate(mhd_labels)
    }
    user = {
        label: [row[index] for row in user_rows]
        for index, label in enumerate(user_labels)
    }
    snapshot = terminal_snapshot(output_dir)
    return {
        "schema_version": 4,
        "execution_epoch": EPOCH,
        "inspected_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "manifest": str(manifest),
        "job_id": job_id,
        "case_id": case_id,
        "segment": segment,
        "required_time": float(segment.rsplit("_t", 1)[1].replace("p", ".")),
        "final_time": final_time,
        "maximum_strict_failure_counts": {
            name: max(mhd[name]) for name in stage_i.STRICT_LF_FAILURE_COLUMNS
        },
        "checks": {
            "required_time_reached": False,
            "strict_lf_failure_counters_zero": True,
            "snapshots_retained": True,
            "terminal_snapshot_retained": True,
            "restart_retained": True,
            "terminal_restart_physical_time_matches_final": True,
            "plasma_continuation_policy": True,
        },
        "accepted": False,
        "clean_for_continuation": True,
        "mhd_history": retained_file(mhd_path),
        "user_history": retained_file(user_path),
        "plasma_continuation_policy": stage_i.CONTINUATION_PLASMA_POLICY,
        "plasma_continuation_evidence": stage_i.continuation_plasma_evidence(
            case_id, mhd, user
        ),
        "snapshots": [snapshot],
        "snapshot_times": [final_time],
        "restarts": [terminal_restart_value],
        "restart_times": [final_time],
        "restart_time_marker_modes": [["full_precision"] * 8],
        "terminal_restart": terminal_restart_value,
        "terminal_restart_time": final_time,
        "restart_time_marker_bypass": False,
        "final_hardwall_projection_count": mhd["lf_hwproj"][-1],
    }


def recorded_manifest(
    fixture: dict[str, object],
    *,
    case_id: str,
    segment: str,
    job_id: str,
    nodes: int,
    result: str,
    final_time: float,
    submitted: datetime,
    completed: datetime,
    elapsed_seconds: int,
    terminal: bool,
    scheduler_state: str = "COMPLETED",
    exit_code: str = "0:0",
) -> Path:
    """Create one production-shaped recorded manifest."""

    root = fixture["root"]
    executable = fixture["executable"]
    inputs = fixture["inputs"]
    revision = fixture["revision"]
    assert isinstance(root, Path)
    assert isinstance(executable, Path)
    assert isinstance(inputs, dict)
    assert isinstance(revision, str)
    input_path = inputs[case_id]
    assert isinstance(input_path, Path)
    run_dir = root / "runs/mks24-stage-i" / EPOCH / case_id / segment
    output_dir = run_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = run_dir / "manifest/prepared_run.json"
    terminal_value = terminal_restart(output_dir) if terminal else None
    actual = nodes * elapsed_seconds / 3600
    inspection: dict[str, object] | None = None
    if result in {"accepted", "clean_partial"}:
        inspection = {
            "schema_version": 4,
            "execution_epoch": EPOCH,
            "job_id": job_id,
            "case_id": case_id,
            "segment": segment,
            "accepted": result == "accepted",
            "clean_for_continuation": True,
            "final_time": final_time,
        }
        if result == "clean_partial":
            assert terminal_value is not None
            inspection = controller_clean_partial_inspection(
                manifest,
                output_dir,
                case_id=case_id,
                segment=segment,
                job_id=job_id,
                final_time=final_time,
                terminal_restart_value=terminal_value,
            )
        if terminal_value is not None:
            inspection["terminal_restart"] = terminal_value
            inspection["terminal_restart_time"] = final_time
    write_json(
        manifest,
        {
            "schema_version": 3,
            "execution_epoch": EPOCH,
            "project_root": str(root),
            "state": "recorded",
            "job_id": job_id,
            "run": {
                "case_id": case_id,
                "case_name": f"fixture_case_{case_id}",
                "segment": segment,
            },
            "allocation": {
                "nodes": nodes,
                "ranks_per_node": 8,
                "cpus_per_task": 7,
                "requested_walltime": "02:00:00" if nodes == 1 else "01:00:00",
            },
            "command": {
                "parent_segment": None,
                "executable_revision": revision,
                "executable_sha256": sha256(executable),
                "input_revision": revision,
                "input_file": str(input_path),
                "input_sha256": sha256(input_path),
            },
            "paths": {"output_dir": str(output_dir)},
            "accounting": {
                "job_id": job_id,
                "execution_epoch": EPOCH,
                "case_id": case_id,
                "case_name": f"fixture_case_{case_id}",
                "segment": segment,
                "state": scheduler_state,
                "exit_code": exit_code,
                "nodes": str(nodes),
                "elapsed_seconds": str(elapsed_seconds),
                "actual_node_hours": f"{actual:.6f}",
                "result": result,
                "submitted_utc": scheduler_time(submitted),
                "completed_utc": scheduler_time(completed),
                "executable_revision": revision,
                "executable_sha256": sha256(executable),
                "input_revision": revision,
                "input_file": str(input_path),
                "output_dir": str(output_dir),
            },
            "scientific_inspection": inspection,
        },
    )
    return manifest


def profile(
    fixture: dict[str, object],
    *,
    case_id: str = "R03",
    segment: str | None = None,
    nodes: int = 1,
    parent: bool = True,
    target: float | None = None,
    estimated_storage_bytes: int | None = None,
) -> dict[str, object]:
    """Return one exact authenticated next-segment profile."""

    repository = fixture["repository"]
    root = fixture["root"]
    source_bundle = fixture["source_bundle"]
    executable = fixture["executable"]
    build_manifest = fixture["build_manifest"]
    revision = fixture["revision"]
    inputs = fixture["inputs"]
    restart = fixture["restart"]
    assert all(
        isinstance(item, Path)
        for item in (repository, root, source_bundle, executable, build_manifest, restart)
    )
    assert isinstance(revision, str)
    assert isinstance(inputs, dict)
    input_path = inputs[case_id]
    assert isinstance(input_path, Path)
    if target is None:
        target = 0.5 if parent else INITIAL_TARGETS[case_id]
    if segment is None:
        segment = (
            "s01_rankio_t0p25_t0p5"
            if parent
            else f"s00_rankio_t0_t{str(target).replace('.', 'p')}"
        )
    if estimated_storage_bytes is None:
        manifest = fixture["manifest"]
        assert isinstance(manifest, Path)
        output_dir = Path(json.loads(manifest.read_text())["paths"]["output_dir"])
        basis_bytes = sum(path.stat().st_size for path in output_dir.rglob("*") if path.is_file())
        start = 0.25 if parent else 0.0
        estimated_storage_bytes = max(
            1024, math.ceil(basis_bytes * (float(target) - start) / 0.25)
        )
        if case_id == "R17":
            estimated_storage_bytes = max(958_271_710_272, estimated_storage_bytes)
    if parent:
        manifest = fixture["manifest"]
        assert isinstance(manifest, Path)
        manifest_value = json.loads(manifest.read_text())
        accounting = manifest_value["accounting"]
        observed_interval = Decimal("0.25")
        elapsed = int(accounting["elapsed_seconds"])
        maximum = observed_interval * Decimal(6600) / Decimal(elapsed) * Decimal("0.90")
        recommended = Decimal(str(target)) - Decimal("0.25")
        recommendation_basis = {
            "kind": "measured-production-continuation",
            "job_id": "12345",
            "case_id": case_id,
            "nodes": nodes,
            "observed_result": "clean_partial",
            "observed_simulation_interval": format(observed_interval, "f"),
            "observed_elapsed_seconds": elapsed,
            "athena_walltime_seconds": 6600,
            "runtime_safety_factor": "0.90",
            "maximum_recommended_interval": format(maximum, "f"),
            "recommended_interval": format(recommended, "f"),
        }
    else:
        qualification = fixture["qualification"]
        assert isinstance(qualification, Path)
        recommendation_basis = {
            "kind": "qualified-initial-calibration",
            "qualification_approval_sha256": sha256(qualification),
            "profile_status": "promoted-by-f113-for-measured-calibration",
        }
    return {
        "acceptance_criterion": acceptance_criterion(case_id, target),
        "acceptance_policy": POLICIES[case_id],
        "athena_walltime": "01:50:00",
        "build_manifest": str(build_manifest),
        "build_manifest_sha256": directory_inventory_sha256(build_manifest),
        "case_id": case_id,
        "controller_walltime_max_seconds": 7200,
        "cpus_per_task": 7,
        "estimated_storage_bytes": estimated_storage_bytes,
        "executable": str(executable),
        "executable_revision": revision,
        "executable_sha256": sha256(executable),
        "input_file": str(input_path.relative_to(repository)),
        "input_revision": revision,
        "input_sha256": sha256(input_path),
        "nodes": nodes,
        "output_layout": "rank-local",
        "segment": segment,
        "parent_job_id": "12345" if parent else None,
        "parent_result": "clean_partial" if parent else None,
        "parent_segment": "s00_rankio_t0_t0p5" if parent else None,
        "restart_file": str(restart) if parent else None,
        "restart_file_sha256": sha256(restart) if parent else None,
        "restart_time": 0.25 if parent else None,
        "ranks_per_node": 8,
        "recommendation_basis": recommendation_basis,
        "source_bundle": str(source_bundle),
        "source_bundle_sha256": sha256(source_bundle),
        "time_tlim_target": target,
        "walltime": "02:00:00",
    }


def refresh_request(fixture: dict[str, object], mutate=None) -> None:
    """Refresh every request binding after an intentional fixture mutation."""

    request_path = fixture["request"]
    root = fixture["root"]
    assert isinstance(request_path, Path)
    assert isinstance(root, Path)
    request = json.loads(request_path.read_text())
    if mutate is not None:
        mutate(request)
    for key, fixture_key in (
        ("reconciliation", "reconcile"),
        ("ledger", "ledger"),
        ("reservations", "reservations"),
        ("storage_evidence", "storage"),
        ("ceiling_evidence", "ceiling"),
        ("ceiling_publication_audit", "ceiling_audit"),
        ("qualification_approval", "qualification"),
        ("predecessor_recost", "predecessor"),
        ("predecessor_recost_independent_review", "predecessor_review"),
        ("predecessor_recost_publication_audit", "predecessor_audit"),
    ):
        path = fixture[fixture_key]
        assert isinstance(path, Path)
        request["inputs"][key]["sha256"] = sha256(path)
    for binding in request["inputs"]["scheduler_evidence"]:
        binding["sha256"] = sha256(root / binding["path"])
    for binding in request["inputs"]["manifests"]:
        binding["sha256"] = sha256(root / binding["path"])
    for key, fixture_key in (
        ("evidence", "source_authority"),
        ("publication_audit", "source_authority_audit"),
        ("provenance_review", "source_authority_provenance_review"),
        ("plasma_review", "source_authority_plasma_review"),
    ):
        path = fixture[fixture_key]
        assert isinstance(path, Path)
        request["inputs"]["source_authority"][key]["sha256"] = sha256(path)
    request["inputs"]["source_authority"]["final_source_bundle"] = dict(
        request["inputs"]["source_bundle"]
    )
    readiness = request["inputs"]["r17_readiness_evidence"]
    if readiness is not None:
        readiness["sha256"] = sha256(root / readiness["path"])
    write_json(request_path, request)
    fixture["request_sha256"] = sha256(request_path)
    review = fixture["request_review"]
    assert isinstance(review, Path)
    write_immutable_json(
        review,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-request-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": request["generated_utc"],
            "decision": "approved-for-evidence-generation",
            "reviewer": {
                "agent_id": "fixture-independent-request-reviewer",
                "role": "independent recost request reviewer",
                "declared_process_independence": True,
                "identity_assurance": REQUEST_REVIEW_IDENTITY_ASSURANCE,
                "identity_assurance_limitation": REQUEST_REVIEW_IDENTITY_LIMITATION,
            },
            "candidate": {
                "path": str(request_path),
                "sha256": fixture["request_sha256"],
            },
            "scope": {"non_authorizing": True},
        },
    )


def computed_storage_projections(
    fixture: dict[str, object], profiles: list[dict[str, object]]
) -> tuple[int, list[dict[str, object]]]:
    """Return the generator's independently observed storage projections."""

    module = load_recost_module()
    root = fixture["root"]
    matrix_path = fixture["matrix"]
    assert isinstance(root, Path)
    assert isinstance(matrix_path, Path)
    matrix, _ = module.parse_matrix(json.loads(matrix_path.read_text()))
    manifests = []
    for path in sorted(root.glob(f"runs/mks24-stage-i/{EPOCH}/*/*/manifest/prepared_run.json")):
        value = json.loads(path.read_text())
        value["_manifest_path"] = str(path)
        value["_manifest_sha256"] = sha256(path)
        manifests.append(value)
    projected, evidence, _ = module.observed_storage_projections(
        root, manifests, matrix, profiles
    )
    return projected, evidence


def refresh_storage(
    fixture: dict[str, object],
    *,
    profiles: list[dict[str, object]] | None = None,
    projected_growth: int | None = None,
    available: int | None = None,
) -> None:
    """Refresh independently measured storage evidence."""

    root = fixture["root"]
    storage = fixture["storage"]
    assert isinstance(root, Path)
    assert isinstance(storage, Path)
    value = json.loads(storage.read_text())
    value["retained_stage_i_bytes"] = retained_bytes(root)
    if profiles is not None:
        projected_growth, projections = computed_storage_projections(fixture, profiles)
        value["projection_method"] = STORAGE_PROJECTION_METHOD
        value["profile_projections_sha256"] = hashlib.sha256(
            (json.dumps(projections, sort_keys=True) + "\n").encode()
        ).hexdigest()
    if projected_growth is not None:
        value["projected_authorized_wave_growth_bytes"] = projected_growth
    if available is not None:
        value["available_bytes"] = available
    write_json(storage, value)


def refresh_f118_historical_f116_binding(fixture: dict[str, object]) -> None:
    """Refresh the F118 wrapper after an intentional historical F116 mutation."""

    authority = fixture["source_authority"]
    provenance = fixture["source_authority_provenance_review"]
    plasma = fixture["source_authority_plasma_review"]
    authority_audit = fixture["source_authority_audit"]
    historical = fixture["f116_source_authority"]
    historical_provenance = fixture["f116_source_authority_provenance_review"]
    historical_plasma = fixture["f116_source_authority_plasma_review"]
    historical_audit = fixture["f116_source_authority_audit"]
    assert all(isinstance(item, Path) for item in (authority, provenance, plasma, authority_audit))
    assert all(
        isinstance(item, Path)
        for item in (historical, historical_provenance, historical_plasma, historical_audit)
    )
    authority_value = json.loads(authority.read_text())
    historical_bindings = authority_value["predecessor_authorities"]["historical_f116"]
    for key, path in (
        ("evidence", historical),
        ("publication_audit", historical_audit),
        ("provenance_review", historical_provenance),
        ("plasma_review", historical_plasma),
    ):
        historical_bindings[key]["sha256"] = sha256(path)
    write_immutable_json(authority, authority_value)
    for review_path in (provenance, plasma):
        review_value = json.loads(review_path.read_text())
        review_value["reviewed_candidate"]["sha256"] = sha256(authority)
        review_value["published_f118"]["sha256"] = sha256(authority)
        write_immutable_json(review_path, review_value)
    authority_audit_value = json.loads(authority_audit.read_text())
    authority_audit_value["artifact"]["sha256"] = sha256(authority)
    authority_audit_value["independent_reviews"][
        "reviews_bind_exact_published_f118_sha256"
    ] = sha256(authority)
    authority_audit_value["independent_reviews"]["provenance_security"]["sha256"] = sha256(
        provenance
    )
    authority_audit_value["independent_reviews"]["plasma_scientific_continuation"][
        "sha256"
    ] = sha256(plasma)
    authority_audit_value["historical_f116_authority"] = {
        "evidence_sha256": sha256(historical),
        "publication_audit_sha256": sha256(historical_audit),
        "provenance_review_sha256": sha256(historical_provenance),
        "plasma_review_sha256": sha256(historical_plasma),
    }
    write_immutable_json(authority_audit, authority_audit_value)


def refresh_ceiling(fixture: dict[str, object]) -> None:
    """Refresh the exact F113 publication binding through F115/F116/F118."""

    ceiling = fixture["ceiling"]
    audit = fixture["ceiling_audit"]
    assert isinstance(ceiling, Path)
    assert isinstance(audit, Path)
    value = json.loads(audit.read_text())
    value["artifact"]["sha256"] = sha256(ceiling)
    write_json(audit, value)
    historical = fixture["historical_source_authority"]
    historical_provenance = fixture["historical_source_authority_provenance_review"]
    historical_plasma = fixture["historical_source_authority_plasma_review"]
    historical_audit = fixture["historical_source_authority_audit"]
    assert all(
        isinstance(item, Path)
        for item in (historical, historical_provenance, historical_plasma, historical_audit)
    )
    historical_value = json.loads(historical.read_text())
    historical_value["predecessors"]["f113_controller_transition"]["sha256"] = sha256(ceiling)
    write_immutable_json(historical, historical_value)
    for review_path in (historical_provenance, historical_plasma):
        review_value = json.loads(review_path.read_text())
        review_value["published_f115"]["sha256"] = sha256(historical)
        write_immutable_json(review_path, review_value)
    historical_audit_value = json.loads(historical_audit.read_text())
    historical_audit_value["artifact"]["sha256"] = sha256(historical)
    historical_audit_value["independent_reviews"][
        "reviews_bind_exact_published_f115_sha256"
    ] = sha256(historical)
    historical_audit_value["independent_reviews"]["provenance_security"]["sha256"] = sha256(
        historical_provenance
    )
    historical_audit_value["independent_reviews"]["plasma_scientific_continuation"][
        "sha256"
    ] = sha256(historical_plasma)
    write_immutable_json(historical_audit, historical_audit_value)

    authority = fixture["f116_source_authority"]
    provenance = fixture["f116_source_authority_provenance_review"]
    plasma = fixture["f116_source_authority_plasma_review"]
    authority_audit = fixture["f116_source_authority_audit"]
    assert all(isinstance(item, Path) for item in (authority, provenance, plasma, authority_audit))
    authority_value = json.loads(authority.read_text())
    historical_bindings = authority_value["predecessor_authorities"]["historical_f115"]
    for key, path in (
        ("evidence", historical),
        ("publication_audit", historical_audit),
        ("provenance_review", historical_provenance),
        ("plasma_review", historical_plasma),
    ):
        historical_bindings[key]["sha256"] = sha256(path)
    write_immutable_json(authority, authority_value)
    for review_path in (provenance, plasma):
        review_value = json.loads(review_path.read_text())
        review_value["reviewed_candidate"]["sha256"] = sha256(authority)
        review_value["published_f116"]["sha256"] = sha256(authority)
        write_immutable_json(review_path, review_value)
    authority_audit_value = json.loads(authority_audit.read_text())
    authority_audit_value["artifact"]["sha256"] = sha256(authority)
    authority_audit_value["independent_reviews"][
        "reviews_bind_exact_published_f116_sha256"
    ] = sha256(authority)
    authority_audit_value["independent_reviews"]["provenance_security"]["sha256"] = sha256(
        provenance
    )
    authority_audit_value["independent_reviews"]["plasma_scientific_continuation"][
        "sha256"
    ] = sha256(plasma)
    authority_audit_value["historical_f115_authority"] = {
        "evidence_sha256": sha256(historical),
        "publication_audit_sha256": sha256(historical_audit),
        "provenance_review_sha256": sha256(historical_provenance),
        "plasma_review_sha256": sha256(historical_plasma),
    }
    write_immutable_json(authority_audit, authority_audit_value)
    refresh_f118_historical_f116_binding(fixture)
    refresh_request(fixture)


def configure_f113_audit_supersession(
    fixture: dict[str, object], *, retain_historical_audit: bool
) -> None:
    """Make F116 exactly adopt the retained F113 artifact's missing audit."""

    root = fixture["root"]
    ceiling_audit = fixture["ceiling_audit"]
    source_authority_audit = fixture["source_authority_audit"]
    assert isinstance(root, Path)
    assert isinstance(ceiling_audit, Path)
    assert isinstance(source_authority_audit, Path)
    if not retain_historical_audit:
        ceiling_audit.unlink()
    fixture["ceiling_audit"] = source_authority_audit

    def mutate(request):
        request["inputs"]["ceiling_publication_audit"]["path"] = str(
            source_authority_audit.relative_to(root)
        )

    refresh_request(fixture, mutate)


def refresh_predecessor(fixture: dict[str, object]) -> None:
    """Refresh the selected predecessor publication binding."""

    predecessor = fixture["predecessor"]
    review = fixture["predecessor_review"]
    audit = fixture["predecessor_audit"]
    assert isinstance(predecessor, Path)
    assert isinstance(review, Path)
    assert isinstance(audit, Path)
    predecessor.chmod(0o444)
    review_value = json.loads(review.read_text())
    review_value["candidate"]["sha256"] = sha256(predecessor)
    write_immutable_json(review, review_value)
    value = json.loads(audit.read_text())
    value["artifact"]["sha256"] = sha256(predecessor)
    value["independent_review"]["sha256"] = sha256(review)
    write_immutable_json(audit, value)
    refresh_request(fixture)


def set_profiles(
    fixture: dict[str, object], profiles: list[dict[str, object]], *, mode: str
) -> None:
    """Replace exact profiles and synchronize reviewed storage growth."""

    refresh_storage(fixture, profiles=profiles)

    def mutate(request):
        request["recommendations"] = {
            "mode": mode,
            "max_wave_nodes": sum(int(item["nodes"]) for item in profiles),
            "profiles": profiles,
        }

    refresh_request(fixture, mutate)


def rewrite_scheduler_from_manifest(fixture: dict[str, object]) -> None:
    """Rewrite the primary scheduler row from its bound manifest."""

    manifest = fixture["manifest"]
    scheduler = fixture["scheduler"]
    assert isinstance(manifest, Path)
    assert isinstance(scheduler, Path)
    value = json.loads(manifest.read_text())
    accounting = value["accounting"]
    scheduler.write_text(
        "|".join(
            (
                accounting["job_id"],
                f"cgl_mks24_{EPOCH_SLUG}_{accounting['case_id']}_{accounting['segment']}",
                accounting["state"],
                accounting["exit_code"],
                accounting["nodes"],
                accounting["elapsed_seconds"],
                accounting["submitted_utc"],
                accounting["completed_utc"],
            )
        )
        + "\n"
    )
    scheduler.chmod(0o644)


def add_recorded_wave_job(fixture: dict[str, object]) -> None:
    """Append a second exact recorded job to form a drained-wave barrier."""

    timestamp = fixture["timestamp"]
    root = fixture["root"]
    ledger = fixture["ledger"]
    reservations = fixture["reservations"]
    reconcile = fixture["reconcile"]
    assert isinstance(timestamp, datetime)
    assert all(isinstance(item, Path) for item in (root, ledger, reservations, reconcile))
    manifest = recorded_manifest(
        fixture,
        case_id="R05",
        segment="s00_rankio_t0_t0p25",
        job_id="23456",
        nodes=2,
        result="accepted",
        final_time=0.25,
        submitted=timestamp - timedelta(minutes=50),
        completed=timestamp - timedelta(minutes=20),
        elapsed_seconds=1800,
        terminal=False,
    )
    manifest_accounting = json.loads(manifest.read_text())["accounting"]
    with ledger.open("a") as stream:
        stream.write(
            f"{EPOCH},23456,{scheduler_time(timestamp - timedelta(minutes=50))},"
            f"{scheduler_time(timestamp - timedelta(minutes=20))},R05,fixture_case_R05,"
            "s00_rankio_t0_t0p25,COMPLETED,0:0,2,01:00:00,1800,2.000000,"
            f"1.000000,2.000000,{manifest_accounting['executable_revision']},"
            f"{manifest_accounting['executable_sha256']},{manifest_accounting['input_revision']},"
            f"{manifest_accounting['input_file']},{manifest_accounting['output_dir']},"
            "accepted,fixture\n"
        )
    reservations_value = json.loads(reservations.read_text())
    reservations_value.append(
        {
            "execution_epoch": EPOCH,
            "job_id": "23456",
            "case_id": "R05",
            "case_name": "fixture_case_R05",
            "segment": "s00_rankio_t0_t0p25",
            "nodes": 2,
            "requested_walltime": "01:00:00",
            "state": "recorded",
            "result": "accepted",
            "actual_node_hours": 1.0,
            "manifest": str(manifest),
        }
    )
    write_json(reservations, reservations_value)
    reconcile_value = json.loads(reconcile.read_text())
    reconcile_value["counts"]["reservations"] += 1
    reconcile_value["counts"]["ledger_rows"] += 1
    reconcile_value["counts"]["manifests"] += 1
    write_json(reconcile, reconcile_value)
    scheduler = root / "accounting/23456.stage_i.sacct.txt"
    scheduler.write_text(
        f"23456|cgl_mks24_{EPOCH_SLUG}_R05_s00_rankio_t0_t0p25|COMPLETED|0:0|2|1800|"
        f"{scheduler_time(timestamp - timedelta(minutes=50))}|"
        f"{scheduler_time(timestamp - timedelta(minutes=20))}\n"
    )
    scheduler.chmod(0o644)
    refresh_storage(fixture)

    def mutate(request):
        request["barrier"]["recorded_segments"].append(
            {
                "case_id": "R05",
                "segment": "s00_rankio_t0_t0p25",
                "job_id": "23456",
                "result": "accepted",
            }
        )
        request["inputs"]["manifests"].append(
            {"path": str(manifest.relative_to(root)), "sha256": sha256(manifest)}
        )
        request["inputs"]["scheduler_evidence"].append(
            {"path": str(scheduler.relative_to(root)), "sha256": sha256(scheduler)}
        )

    refresh_request(fixture, mutate)


def add_cancelled_identity(fixture: dict[str, object]) -> Path:
    """Add one authenticated no-start cancelled identity without ledger use."""

    root = fixture["root"]
    reservations = fixture["reservations"]
    reconcile = fixture["reconcile"]
    timestamp = fixture["timestamp"]
    assert all(isinstance(item, Path) for item in (root, reservations, reconcile))
    assert isinstance(timestamp, datetime)
    manifest = (
        root
        / "runs/mks24-stage-i"
        / EPOCH
        / "R04/s00_rankio_t0_t0p25/manifest/prepared_run.json"
    )
    write_json(
        manifest,
        {
            "schema_version": 3,
            "execution_epoch": EPOCH,
            "state": "cancelled",
            "job_id": None,
            "run": {
                "case_id": "R04",
                "case_name": "fixture_case_R04",
                "segment": "s00_rankio_t0_t0p25",
            },
            "allocation": {
                "nodes": 1,
                "ranks_per_node": 8,
                "cpus_per_task": 7,
                "requested_walltime": "02:00:00",
            },
            "accounting": None,
            "scientific_inspection": None,
            "cancellation": {
                "cancelled_utc": timestamp.isoformat(),
                "notes": "Fixture no-start cancellation.",
            },
        },
    )
    reservations_value = json.loads(reservations.read_text())
    reservations_value.append(
        {
            "execution_epoch": EPOCH,
            "case_id": "R04",
            "case_name": "fixture_case_R04",
            "segment": "s00_rankio_t0_t0p25",
            "nodes": 1,
            "requested_walltime": "02:00:00",
            "state": "cancelled",
            "notes": "Fixture no-start cancellation.",
            "manifest": str(manifest),
        }
    )
    write_json(reservations, reservations_value)
    reconcile_value = json.loads(reconcile.read_text())
    reconcile_value["counts"]["reservations"] += 1
    reconcile_value["counts"]["manifests"] += 1
    write_json(reconcile, reconcile_value)
    refresh_storage(fixture)

    def mutate(request):
        request["inputs"]["manifests"].append(
            {"path": str(manifest.relative_to(root)), "sha256": sha256(manifest)}
        )

    refresh_request(fixture, mutate)
    return manifest


def replace_barrier_with_failed_outcome(
    fixture: dict[str, object],
    *,
    result: str = "failed",
    scheduler_state: str = "FAILED",
    exit_code: str = "1:0",
) -> Path:
    """Append one non-scientific terminal job and make it the exact recost barrier."""

    root = fixture["root"]
    ledger = fixture["ledger"]
    reservations = fixture["reservations"]
    reconcile = fixture["reconcile"]
    timestamp = fixture["timestamp"]
    assert all(isinstance(item, Path) for item in (root, ledger, reservations, reconcile))
    assert isinstance(timestamp, datetime)
    manifest = recorded_manifest(
        fixture,
        case_id="R05",
        segment="s00_rankio_t0_t0p25",
        job_id="34567",
        nodes=2,
        result=result,
        final_time=0.0,
        submitted=timestamp - timedelta(minutes=30),
        completed=timestamp - timedelta(minutes=20),
        elapsed_seconds=600,
        terminal=False,
        scheduler_state=scheduler_state,
        exit_code=exit_code,
    )
    accounting = json.loads(manifest.read_text())["accounting"]
    with ledger.open("a") as stream:
        stream.write(
            f"{EPOCH},34567,{accounting['submitted_utc']},{accounting['completed_utc']},"
            f"R05,fixture_case_R05,s00_rankio_t0_t0p25,{scheduler_state},{exit_code},"
            "2,01:00:00,600,"
            "2.000000,0.333333,1.333333,"
            f"{accounting['executable_revision']},{accounting['executable_sha256']},"
            f"{accounting['input_revision']},{accounting['input_file']},"
            f"{accounting['output_dir']},{result},fixture terminal outcome\n"
        )
    reservations_value = json.loads(reservations.read_text())
    reservations_value.append(
        {
            "execution_epoch": EPOCH,
            "job_id": "34567",
            "case_id": "R05",
            "case_name": "fixture_case_R05",
            "segment": "s00_rankio_t0_t0p25",
            "nodes": 2,
            "requested_walltime": "01:00:00",
            "state": "recorded",
            "result": result,
            "actual_node_hours": 0.333333,
            "manifest": str(manifest),
        }
    )
    write_json(reservations, reservations_value)
    reconcile_value = json.loads(reconcile.read_text())
    reconcile_value["counts"]["reservations"] += 1
    reconcile_value["counts"]["ledger_rows"] += 1
    reconcile_value["counts"]["manifests"] += 1
    write_json(reconcile, reconcile_value)
    scheduler = root / "accounting/34567.stage_i.sacct.txt"
    scheduler.write_text(
        f"34567|cgl_mks24_{EPOCH_SLUG}_R05_s00_rankio_t0_t0p25|"
        f"{scheduler_state}|{exit_code}|2|600|"
        f"{accounting['submitted_utc']}|{accounting['completed_utc']}\n"
    )
    scheduler.chmod(0o644)
    refresh_storage(fixture)

    def mutate(request):
        request["barrier"]["recorded_segments"] = [
            {
                "case_id": "R05",
                "segment": "s00_rankio_t0_t0p25",
                "job_id": "34567",
                "result": result,
            }
        ]
        request["inputs"]["manifests"].append(
            {"path": str(manifest.relative_to(root)), "sha256": sha256(manifest)}
        )
        request["inputs"]["scheduler_evidence"] = [
            {"path": str(scheduler.relative_to(root)), "sha256": sha256(scheduler)}
        ]

    refresh_request(fixture, mutate)
    return manifest


def generator_command(
    fixture: dict[str, object],
    *,
    output: Path | None = None,
    generator_sha256: str | None = None,
) -> list[str]:
    """Build one complete local-fixture invocation."""

    generator = fixture["generator"]
    root = fixture["root"]
    request = fixture["request"]
    queue = fixture["queue"]
    default_output = fixture["output"]
    assert all(isinstance(item, Path) for item in (generator, root, request, queue, default_output))
    return [
        sys.executable,
        str(generator),
        "--root",
        str(root),
        "--allow-local-root",
        "--request",
        str(request),
        "--expected-request-sha256",
        str(fixture["request_sha256"]),
        "--output",
        str(output or default_output),
        "--expected-generator-sha256",
        generator_sha256 or sha256(generator),
        "--squeue-file",
        str(queue),
    ]


def run_generator(
    fixture: dict[str, object],
    *,
    output: Path | None = None,
    generator_sha256: str | None = None,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    """Run the retained fixture generator without raising."""

    return subprocess.run(
        generator_command(fixture, output=output, generator_sha256=generator_sha256),
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )


def run_action(
    fixture: dict[str, object],
    action: str,
    *arguments: str,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    """Run one authenticated prerequisite action against the local fixture."""

    generator = fixture["generator"]
    root = fixture["root"]
    queue = fixture["queue"]
    assert all(isinstance(item, Path) for item in (generator, root, queue))
    return subprocess.run(
        [
            sys.executable,
            str(generator),
            "--root",
            str(root),
            "--allow-local-root",
            "--expected-generator-sha256",
            sha256(generator),
            "--squeue-file",
            str(queue),
            action,
            *arguments,
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )


def write_draft_packet(
    fixture: dict[str, object],
    *,
    checkpoint_number: int = 201,
    mutate=None,
) -> Path:
    """Write one explicit non-authorizing request-draft packet."""

    request = fixture["request"]
    accounting = fixture["accounting"]
    assert isinstance(request, Path)
    assert isinstance(accounting, Path)
    packet_value = json.loads(request.read_text())
    packet_value["schema_version"] = 1
    packet_value["record_type"] = "stage-i-recost-request-draft-packet"
    packet_value["checkpoint"] = f"F-{checkpoint_number}"
    packet_value["artifact_name"] = (
        f"mks24_stage_i_{EPOCH_SLUG}_F{checkpoint_number}_recost_evidence.json"
    )
    for key in (
        "reconciliation",
        "ledger",
        "reservations",
        "manifests",
        "scheduler_evidence",
        "storage_evidence",
    ):
        packet_value["inputs"][key] = None
    packet_value["draft_policy"] = {
        "independent_review_created": False,
        "self_approved": False,
        "scheduler_mutation_authorized": False,
        "canonical_mutation_authorized": False,
        "required_storage_safety_bytes": 1024,
    }
    if mutate is not None:
        mutate(packet_value)
    packet = accounting / (
        f"mks24_stage_i_{EPOCH_SLUG}_F{checkpoint_number}_recost_draft_packet.json"
    )
    write_json(packet, packet_value)
    return packet


def write_f119_seed_and_candidate(
    fixture: dict[str, object],
    *,
    seed_mutate=None,
    candidate_mutate=None,
) -> tuple[Path, Path]:
    """Write one exact four-profile F117 seed and derived F119 candidate."""

    predecessor = fixture["predecessor"]
    assert isinstance(predecessor, Path)
    if "_F117_recost_evidence.json" not in predecessor.name:
        replace_fixture_predecessor(fixture, 117)
    profiles = [
        profile(fixture, case_id=case_id, nodes=nodes, parent=case_id == "R03")
        for case_id, nodes in (("R03", 1), ("R04", 4), ("R12", 4), ("R16", 1))
    ]

    def mutate_seed(packet):
        bind_fixture_predecessor(fixture, packet)
        packet["recommendations"] = {
            "mode": "bounded-wave",
            "max_wave_nodes": 10,
            "profiles": profiles,
        }
        packet["draft_policy"]["required_storage_safety_bytes"] = 1024**4
        if seed_mutate is not None:
            seed_mutate(packet)

    seed = write_draft_packet(fixture, checkpoint_number=117, mutate=mutate_seed)
    candidate_value = json.loads(seed.read_text())
    candidate_value.update(
        {
            "checkpoint": "F-119",
            "artifact_name": (
                "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json"
            ),
            "requested_by": F119_REQUESTED_BY,
            "scope": F119_SCOPE,
        }
    )
    candidate_value["inputs"]["ceiling_publication_audit"] = {
        "path": fixture["source_authority_audit"].relative_to(
            fixture["root"]
        ).as_posix(),
        "sha256": sha256(fixture["source_authority_audit"]),
    }
    candidate_value["inputs"]["stage_i_helper"] = {
        "revision": fixture["revision"],
        "sha256": sha256(fixture["helper"]),
    }
    if candidate_mutate is not None:
        candidate_mutate(candidate_value)
    candidate = seed.with_name(
        "mks24_stage_i_E03_forcing_policy_F119_recost_draft_packet.json"
    )
    write_json(candidate, candidate_value)
    return seed, candidate


def f119_render_arguments(
    fixture: dict[str, object],
    seed: Path,
    *,
    generated: datetime | None = None,
    expires: datetime | None = None,
) -> tuple[str, ...]:
    """Return the exact caller-pinned deterministic F119 renderer arguments."""

    timestamp = generated or fixture["timestamp"]
    assert isinstance(timestamp, datetime)
    expiry = expires or timestamp + timedelta(hours=12)
    return (
        "--generated-utc",
        timestamp.isoformat(),
        "--expires-utc",
        expiry.isoformat(),
        "--f117-seed-packet",
        str(seed),
        "--expected-f117-seed-packet-sha256",
        sha256(seed),
        "--f118-source-bundle",
        str(fixture["source_bundle"]),
        "--expected-f118-source-bundle-sha256",
        sha256(fixture["source_bundle"]),
        "--expected-f118-evidence-sha256",
        sha256(fixture["source_authority"]),
        "--expected-f118-publication-audit-sha256",
        sha256(fixture["source_authority_audit"]),
        "--expected-f118-provenance-review-sha256",
        sha256(fixture["source_authority_provenance_review"]),
        "--expected-f118-plasma-review-sha256",
        sha256(fixture["source_authority_plasma_review"]),
    )


def replace_fixture_predecessor(
    fixture: dict[str, object], checkpoint_number: int
) -> tuple[Path, Path, Path]:
    """Replace the fixture predecessor with one exact promoted checkpoint."""

    root = fixture["root"]
    accounting = fixture["accounting"]
    timestamp = fixture["timestamp"]
    old_paths = (
        fixture["predecessor"],
        fixture["predecessor_review"],
        fixture["predecessor_audit"],
    )
    assert isinstance(root, Path)
    assert isinstance(accounting, Path)
    assert isinstance(timestamp, datetime)
    assert all(isinstance(path, Path) for path in old_paths)
    for path in old_paths:
        path.unlink()

    predecessor = (
        accounting
        / f"mks24_stage_i_{EPOCH_SLUG}_F{checkpoint_number}_recost_evidence.json"
    )
    write_immutable_json(
        predecessor,
        {
            "schema_version": 2,
            "record_type": "stage-i-recost-recommendation-evidence",
            "checkpoint": f"F-{checkpoint_number}",
            "artifact_name": predecessor.name,
            "execution_epoch": EPOCH,
            "generated_utc": (timestamp - timedelta(minutes=8)).isoformat(),
            "authority": {
                "authorizing": False,
                "action_authority": "none-until-independent-review-and-publication",
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            },
        },
    )
    review = predecessor.with_name(f"{predecessor.name}.independent_review.json")
    write_immutable_json(
        review,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (timestamp - timedelta(minutes=7)).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": {
                "agent_id": f"fixture-f{checkpoint_number}-reviewer",
                "independent_from_generator": True,
            },
            "candidate": {"path": str(predecessor), "sha256": sha256(predecessor)},
            "scope": {"non_authorizing": True},
        },
    )
    audit = predecessor.with_name(f"{predecessor.name}.publication_audit.json")
    write_immutable_json(
        audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=5)).isoformat(),
            "artifact": {
                "path": str(predecessor),
                "sha256": sha256(predecessor),
                "mode": "0444",
                "links": 1,
            },
            "independent_review": {
                "path": str(review),
                "sha256": sha256(review),
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
    fixture.update(
        {
            "predecessor": predecessor,
            "predecessor_review": review,
            "predecessor_audit": audit,
        }
    )
    return predecessor, review, audit


def bind_fixture_predecessor(
    fixture: dict[str, object], packet: dict[str, object]
) -> None:
    """Bind a draft packet to the fixture's current promoted predecessor."""

    root = fixture["root"]
    assert isinstance(root, Path)
    inputs = packet["inputs"]
    assert isinstance(inputs, dict)
    for key, fixture_key in (
        ("predecessor_recost", "predecessor"),
        ("predecessor_recost_independent_review", "predecessor_review"),
        ("predecessor_recost_publication_audit", "predecessor_audit"),
    ):
        path = fixture[fixture_key]
        assert isinstance(path, Path)
        inputs[key] = {
            "path": path.relative_to(root).as_posix(),
            "sha256": sha256(path),
        }


def write_external_request_review(
    fixture: dict[str, object],
    request: Path,
    *,
    mutate=None,
) -> Path:
    """Retain one immutable external review candidate for managed installation."""

    root = fixture["root"]
    assert isinstance(root, Path)
    request_value = json.loads(request.read_text())
    review = {
        "schema_version": 1,
        "record_type": "stage-i-recost-request-independent-review",
        "execution_epoch": EPOCH,
        "reviewed_utc": request_value["generated_utc"],
        "decision": "approved-for-evidence-generation",
        "reviewer": {
            "agent_id": "external-draft-request-reviewer",
            "role": "independent recost request reviewer",
            "declared_process_independence": True,
            "identity_assurance": REQUEST_REVIEW_IDENTITY_ASSURANCE,
            "identity_assurance_limitation": REQUEST_REVIEW_IDENTITY_LIMITATION,
        },
        "candidate": {"path": str(request), "sha256": sha256(request)},
        "scope": {"non_authorizing": True},
    }
    if mutate is not None:
        mutate(review)
    candidate = root.parent / f"{request.name}.external-independent-review.json"
    write_immutable_json(candidate, review)
    return candidate


def assert_rejected(result: subprocess.CompletedProcess, pattern: str) -> None:
    """Require one fail-closed rejection."""

    assert result.returncode == 1
    assert pattern in result.stderr


def action_report(result: subprocess.CompletedProcess) -> dict[str, object]:
    """Parse the final JSON action report after retained Git query output."""

    return json.loads(result.stdout[result.stdout.index("{"):])


def load_recost_module():
    """Load the retained generator for direct pure-boundary unit probes."""

    name = f"_cgl_lf_stage_i_recost_test_module_{os.getpid()}_{id(object())}"
    spec = importlib.util.spec_from_file_location(name, RECOST)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def recost_fixture(tmp_path):
    """Create one committed source repository and drained authenticated barrier."""

    timestamp = datetime.now(timezone.utc).replace(microsecond=0)
    repository = tmp_path / "repository"
    frontier = repository / "scripts/frontier"
    matrix_dir = repository / "inputs/cgl_lf_paper"
    frontier.mkdir(parents=True)
    matrix_dir.mkdir(parents=True)
    generator = frontier / RECOST.name
    generator.write_bytes(RECOST.read_bytes())
    generator.chmod(0o755)
    helper = frontier / "cgl_lf_stage_i.py"
    helper.write_text(
        "from pathlib import Path\n"
        "import json\n"
        "import sys\n"
        "\n"
        "def reconcile_report(root):\n"
        "    return json.loads(\n"
        "        (root / 'accounting/reviewed_reconciliation.json').read_text()\n"
        "    )\n"
        "\n"
        "if __name__ == '__main__':\n"
        "    root = Path(sys.argv[sys.argv.index('--root') + 1])\n"
        "    print(json.dumps(reconcile_report(root), indent=2, sort_keys=True))\n"
    )
    helper.chmod(0o644)
    for relative, mode in F116_REQUIRED_COMMITTED_TOOLS.items():
        path = repository / relative
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"#!/usr/bin/env python3\n# fixture {path.name}\n")
        path.chmod(int(mode, 8))
    inputs: dict[str, Path] = {}
    cases = []
    for number in range(2, 18):
        case_id = f"R{number:02d}"
        input_path = matrix_dir / f"fixture_{case_id}.athinput"
        input_path.write_text(f"<problem>\ncase = {case_id}\n")
        input_path.chmod(0o644)
        inputs[case_id] = input_path
        cases.append(
            {
                "id": case_id,
                "name": f"fixture_case_{case_id}",
                "input": str(input_path.relative_to(repository)),
                "resolution": "192x192x384",
                "estimated_node_hours": 10.0,
            }
        )
    matrix = matrix_dir / "mks24_stage_i_manifest.json"
    write_json(
        matrix,
        {
            "schema_version": 1,
            "authorization": {
                "project_budget_node_hours": 4000.0,
                "execution_rule": "Execute corrected E03 and complete R17 last.",
            },
            "cases": cases,
        },
    )
    git(repository, "init", "-q")
    git(repository, "add", ".")
    git(
        repository,
        "-c",
        "user.name=CGL fixture",
        "-c",
        "user.email=cgl-fixture@example.invalid",
        "commit",
        "-q",
        "-m",
        "Create recost fixture source",
    )
    f113_revision = git(repository, "rev-parse", "HEAD")
    f113_helper_sha256 = sha256(helper)
    helper.write_text(helper.read_text() + "# current F115 helper revision\n")
    git(repository, "add", str(helper.relative_to(repository)))
    git(
        repository,
        "-c",
        "user.name=CGL fixture",
        "-c",
        "user.email=cgl-fixture@example.invalid",
        "commit",
        "-q",
        "-m",
        "Supersede F113 helper under F115",
    )
    revision = git(repository, "rev-parse", "HEAD")

    root = tmp_path / "root"
    accounting = root / "accounting"
    source_archives = root / "source-archives"
    transaction_store = accounting / f"mks24_stage_i_{EPOCH_SLUG}_transactions"
    recost_transaction_store = accounting / f"mks24_stage_i_{EPOCH_SLUG}_recost_transactions"
    for directory in (root, accounting, source_archives, transaction_store, recost_transaction_store):
        directory.mkdir(parents=True, exist_ok=True)
        directory.chmod(0o755)
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock.write_text("")
    lock.chmod(0o644)

    source_bundle = source_archives / "athenak-fixture.bundle"
    git(repository, "bundle", "create", str(source_bundle), "HEAD")
    source_bundle.chmod(0o644)
    f116_source_bundle = source_archives / "athenak-historical-f116.bundle"
    f116_source_bundle.write_bytes(source_bundle.read_bytes())
    f116_source_bundle.chmod(0o644)
    historical_source_bundle = source_archives / "athenak-historical-f115.bundle"
    git(
        repository,
        "update-ref",
        "refs/heads/feature/cgl-landau-fluid",
        f113_revision,
    )
    git(
        repository,
        "bundle",
        "create",
        str(historical_source_bundle),
        "refs/heads/feature/cgl-landau-fluid",
    )
    historical_source_bundle.chmod(0o644)
    source_archive_readme = source_archives / "README.md"
    source_archive_readme.write_text(
        "Fixture F118 current source catalog retaining "
        f"{historical_source_bundle.name}, {f116_source_bundle.name}, and "
        f"{source_bundle.name}.\n"
    )
    source_archive_readme.chmod(0o644)
    source_archive_sums = source_archives / "SHA256SUMS"
    source_archive_sums.write_text(
        f"{sha256(historical_source_bundle)}  {historical_source_bundle.name}\n"
        f"{sha256(f116_source_bundle)}  {f116_source_bundle.name}\n"
        f"{sha256(source_bundle)}  {source_bundle.name}\n"
    )
    source_archive_sums.chmod(0o644)
    executable = root / "build/frontier-fixture/src/athena"
    executable.parent.mkdir(parents=True)
    executable.write_bytes(b"fixture executable\n")
    executable.chmod(0o755)
    build_manifest = root / "runs/build-manifests/fixture-build"
    build_manifest.mkdir(parents=True)
    (build_manifest / "athena.sha256").write_text(f"{sha256(executable)}  {executable}\n")
    (build_manifest / "environment.txt").write_text(f"git_revision={revision}\n")
    (build_manifest / "athena-config.txt").write_text("MPI parallelism: ON\n")
    for path in build_manifest.iterdir():
        path.chmod(0o644)

    fixture: dict[str, object] = {
        "timestamp": timestamp,
        "repository": repository,
        "generator": generator,
        "helper": helper,
        "matrix": matrix,
        "revision": revision,
        "f113_revision": f113_revision,
        "f113_helper_sha256": f113_helper_sha256,
        "root": root,
        "accounting": accounting,
        "inputs": inputs,
        "executable": executable,
        "build_manifest": build_manifest,
        "source_bundle": source_bundle,
        "f116_source_bundle": f116_source_bundle,
        "historical_source_bundle": historical_source_bundle,
        "transaction_store": transaction_store,
        "recost_transaction_store": recost_transaction_store,
    }
    qualification = accounting / QUALIFICATION_NAME
    write_json(
        qualification,
        {
            "schema_version": 1,
            "execution_epoch": EPOCH,
            "approval_scope": "fixture corrected-build production qualification",
            "approved_by": "fixture independent build reviewer",
            "approved_utc": (timestamp - timedelta(days=1)).isoformat(),
            "approved_executable": str(executable),
            "approved_executable_revision": revision,
            "approved_executable_sha256": sha256(executable),
            "build_manifest": str(build_manifest),
            "review_notes": "Fixture qualification binds exact executable and build manifest.",
        },
    )
    fixture["qualification"] = qualification
    manifest = recorded_manifest(
        fixture,
        case_id="R03",
        segment="s00_rankio_t0_t0p5",
        job_id="12345",
        nodes=1,
        result="clean_partial",
        final_time=0.25,
        submitted=timestamp - timedelta(hours=2),
        completed=timestamp - timedelta(hours=1),
        elapsed_seconds=3600,
        terminal=True,
    )
    manifest_value = json.loads(manifest.read_text())
    restart = Path(manifest_value["scientific_inspection"]["terminal_restart"]["path"])
    manifest_accounting = manifest_value["accounting"]
    fixture["manifest"] = manifest
    fixture["restart"] = restart

    ledger = accounting / f"mks24_stage_i_{EPOCH_SLUG}_node_hours.csv"
    ledger.write_text(
        "execution_epoch,job_id,submitted_utc,completed_utc,case_id,case_name,"
        "segment,state,exit_code,nodes,requested_walltime,elapsed_seconds,"
        "reserved_node_hours,actual_node_hours,cumulative_stage_i_node_hours,"
        "executable_revision,executable_sha256,input_revision,input_file,"
        "output_dir,result,notes\n"
        f"{EPOCH},12345,{scheduler_time(timestamp - timedelta(hours=2))},"
        f"{scheduler_time(timestamp - timedelta(hours=1))},R03,fixture_case_R03,"
        "s00_rankio_t0_t0p5,COMPLETED,0:0,1,02:00:00,3600,2.000000,1.000000,"
        f"1.000000,{manifest_accounting['executable_revision']},"
        f"{manifest_accounting['executable_sha256']},{manifest_accounting['input_revision']},"
        f"{manifest_accounting['input_file']},{manifest_accounting['output_dir']},"
        "clean_partial,fixture\n"
    )
    ledger.chmod(0o644)
    reservations = accounting / f"mks24_stage_i_{EPOCH_SLUG}_reservations.json"
    write_json(
        reservations,
        [
            {
                "execution_epoch": EPOCH,
                "job_id": "12345",
                "case_id": "R03",
                "case_name": "fixture_case_R03",
                "segment": "s00_rankio_t0_t0p5",
                "nodes": 1,
                "requested_walltime": "02:00:00",
                "state": "recorded",
                "result": "clean_partial",
                "actual_node_hours": 1.0,
                "manifest": str(manifest),
            }
        ],
    )
    reconcile = accounting / "reviewed_reconciliation.json"
    write_json(
        reconcile,
        {
            "execution_epoch": EPOCH,
            "root": str(root),
            "consistent": True,
            "issues": [],
            "counts": {
                "transactions": 0,
                "reservations": 1,
                "active_reservations": 0,
                "ledger_rows": 1,
                "manifests": 1,
            },
        },
    )
    scheduler = accounting / "12345.stage_i.sacct.txt"
    scheduler.write_text(
        f"12345|cgl_mks24_{EPOCH_SLUG}_R03_s00_rankio_t0_t0p5|COMPLETED|0:0|1|3600|"
        f"{scheduler_time(timestamp - timedelta(hours=2))}|"
        f"{scheduler_time(timestamp - timedelta(hours=1))}\n"
    )
    scheduler.chmod(0o644)

    ceiling = accounting / F113_NAME
    write_json(
        ceiling,
        {
            "schema_version": 1,
            "record_type": "stage-i-controller-transition-evidence",
            "checkpoint": "F-113",
            "execution_epoch": EPOCH,
            "generated_utc": (timestamp - timedelta(minutes=12)).isoformat(),
            "implementation": {
                "promoted_commit": f113_revision,
                "stage_i_helper": {
                    "path": "scripts/frontier/cgl_lf_stage_i.py",
                    "sha256": f113_helper_sha256,
                },
            },
            "promoted_controls": {
                "campaign_budget_node_hours": 1400.0,
                "project_budget_node_hours": 4000.0,
                "standard_active_lane_limit": 4,
                "prepared_packet_limit": 1,
                "standard_multi_node_profiles": {
                    "R03": [1],
                    "R04-R15": [1, 2, 4],
                    "R16": [1, 2],
                    "R17": [8],
                },
                "r17_policy": "R17 remains exclusive and last.",
            },
        },
    )
    ceiling_audit = accounting / f"{F113_NAME}.publication_audit.json"
    write_json(
        ceiling_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-controller-transition-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=10)).isoformat(),
            "artifact": {"path": str(ceiling), "sha256": sha256(ceiling), "mode": "0644"},
            "review": {
                "status": "approved",
                "reviewed_by": AUTHORIZED_REVIEWER,
                "authority": F113_REVIEW_AUTHORITY,
            },
        },
    )
    historical_source_authority = accounting / F115_NAME
    write_immutable_json(
        historical_source_authority,
        {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-recovery-supersession-evidence",
            "checkpoint": "F-115",
            "execution_epoch": EPOCH,
            "implementation": {
                "source_bundle": {
                    "path": historical_source_bundle.relative_to(root).as_posix(),
                    "sha256": sha256(historical_source_bundle),
                    "complete_history": True,
                }
            },
            "predecessors": {
                "f113_controller_transition": {
                    "path": ceiling.relative_to(root).as_posix(),
                    "sha256": sha256(ceiling),
                }
            },
        },
    )
    historical_provenance_review = accounting / f"{F115_NAME}.provenance_security_review.json"
    historical_plasma_review = accounting / f"{F115_NAME}.plasma_scientific_review.json"
    for path, kind, decision, agent in (
        (
            historical_provenance_review,
            "provenance-security",
            "approved-for-publication",
            "fixture-historical-provenance-reviewer",
        ),
        (
            historical_plasma_review,
            "plasma-scientific-continuation",
            "approved",
            "fixture-historical-plasma-reviewer",
        ),
    ):
        write_immutable_json(
            path,
            {
                "schema_version": 1,
                "record_type": "stage-i-source-bundle-recovery-supersession-independent-review",
                "checkpoint": "F-115",
                "execution_epoch": EPOCH,
                "review_kind": kind,
                "decision": decision,
                "published_f115": {
                    "path": str(historical_source_authority),
                    "sha256": sha256(historical_source_authority),
                },
                "reviewer": {"agent_id": agent},
            },
        )
    historical_source_authority_audit = accounting / f"{F115_NAME}.publication_audit.json"
    write_immutable_json(
        historical_source_authority_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-source-bundle-recovery-supersession-publication-audit",
            "checkpoint": "F-115",
            "execution_epoch": EPOCH,
            "artifact": {
                "path": str(historical_source_authority),
                "sha256": sha256(historical_source_authority),
                "mode": "0444",
                "links": 1,
            },
            "independent_reviews": {
                "reviews_bind_exact_published_f115_sha256": sha256(historical_source_authority),
                "provenance_security": {
                    "path": str(historical_provenance_review),
                    "sha256": sha256(historical_provenance_review),
                    "mode": "0444",
                    "links": 1,
                },
                "plasma_scientific_continuation": {
                    "path": str(historical_plasma_review),
                    "sha256": sha256(historical_plasma_review),
                    "mode": "0444",
                    "links": 1,
                },
            },
            "authority_and_enforcement": {"direct_sbatch_authorized": False},
        },
    )
    source_authority = accounting / F116_NAME
    historical_bindings = {
        "evidence": {
            "path": historical_source_authority.relative_to(root).as_posix(),
            "sha256": sha256(historical_source_authority),
        },
        "publication_audit": {
            "path": historical_source_authority_audit.relative_to(root).as_posix(),
            "sha256": sha256(historical_source_authority_audit),
        },
        "provenance_review": {
            "path": historical_provenance_review.relative_to(root).as_posix(),
            "sha256": sha256(historical_provenance_review),
        },
        "plasma_review": {
            "path": historical_plasma_review.relative_to(root).as_posix(),
            "sha256": sha256(historical_plasma_review),
        },
    }
    current_bundle = {
        "candidate_path": str(f116_source_bundle),
        "path": f116_source_bundle.relative_to(root).as_posix(),
        "sha256": sha256(f116_source_bundle),
        "complete_history": True,
        "head": revision,
        "advertised_tip": {"revision": revision, "name": "HEAD"},
        "verified_revisions": [f113_revision, revision],
        "selected_as_current": True,
        "subject": "Fixture historical F116 source authority",
    }
    bridge_bundle = {
        "path": historical_source_bundle.relative_to(root).as_posix(),
        "sha256": sha256(historical_source_bundle),
        "complete_history": True,
        "head": f113_revision,
        "advertised_tip": {
            "revision": f113_revision,
            "name": "refs/heads/feature/cgl-landau-fluid",
        },
        "verified_revisions": [f113_revision],
        "selected_as_current": False,
        "role": "retained-non-current-bridge",
    }
    f116_catalog_after = {
        "readme_sha256": "1" * 64,
        "sha256sums_sha256": "2" * 64,
        "bridge_listed_exactly_once": True,
        "final_bundle_listed_exactly_once": True,
        "corrupt_c7_listed": False,
        "historical_f115_preserved": True,
        "sole_current_source_bundle": current_bundle["path"],
    }
    write_immutable_json(
        source_authority,
        {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-evidence",
            "checkpoint": "F-116",
            "execution_epoch": EPOCH,
            "generated_utc": (timestamp - timedelta(minutes=9)).isoformat(),
            "scope": {
                "relationship": "current-source-selection-only-supersession",
                "summary": "Select exact current fixture source without execution authority.",
                "preserves": ["Historical F115 authority and all scientific controls."],
                "does_not_authorize": [
                    "prepare",
                    "submit",
                    "direct sbatch",
                    "scheduler mutation",
                ],
            },
            "predecessor_authorities": {"historical_f115": historical_bindings},
            "implementation": {
                "publisher": {"fixture": True},
                "committed_tools": f116_committed_tools(repository, revision),
                "intermediate_36140_bundle": bridge_bundle,
                "current_source_bundle": current_bundle,
            },
            "source_archive_catalog": {
                "before": {"fixture": "before"},
                "after": f116_catalog_after,
            },
            "authorization": F116_AUTHORIZATION,
            "validation": {"fixture": "passed"},
            "publication_requirements": {"fixture": "reviewed publication"},
        },
    )
    f116_verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": sha256(f116_source_bundle),
        "final_head": revision,
        "historical_f115_preserved": True,
    }
    source_authority_provenance_review = accounting / f"{F116_NAME}.provenance_security_review.json"
    source_authority_plasma_review = accounting / f"{F116_NAME}.plasma_scientific_review.json"
    for path, kind, decision, agent in (
        (
            source_authority_provenance_review,
            "provenance-security",
            "approved-for-publication",
            "fixture-current-provenance-reviewer",
        ),
        (
            source_authority_plasma_review,
            "plasma-scientific-continuation",
            "approved",
            "fixture-current-plasma-reviewer",
        ),
    ):
        write_immutable_json(
            path,
            {
                "schema_version": 1,
                "record_type": "stage-i-current-source-authority-supersession-independent-review",
                "checkpoint": "F-116",
                "execution_epoch": EPOCH,
                "review_kind": kind,
                "decision": decision,
                "reviewed_candidate": {
                    "path": str(source_authority.with_suffix(".json.candidate")),
                    "sha256": sha256(source_authority),
                },
                "published_f116": {
                    "path": str(source_authority),
                    "sha256": sha256(source_authority),
                },
                "reviewer": {"agent_id": agent, "identity": f"fixture {kind} reviewer"},
                "reviewed_utc": (timestamp - timedelta(minutes=8)).isoformat(),
                "findings": ["Exact current-source-only authority verified."],
                "limitations": ["No prepare, submit, scheduler, or science authority."],
                "verified": f116_verified,
            },
        )
    source_authority_audit = accounting / f"{F116_NAME}.publication_audit.json"
    write_immutable_json(
        source_authority_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-publication-audit",
            "checkpoint": "F-116",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=7)).isoformat(),
            "artifact": {
                "path": str(source_authority),
                "sha256": sha256(source_authority),
                "mode": "0444",
                "links": 1,
            },
            "independent_reviews": {
                "reviews_bind_exact_published_f116_sha256": sha256(source_authority),
                "provenance_security": {
                    "path": str(source_authority_provenance_review),
                    "sha256": sha256(source_authority_provenance_review),
                    "mode": "0444",
                    "links": 1,
                },
                "plasma_scientific_continuation": {
                    "path": str(source_authority_plasma_review),
                    "sha256": sha256(source_authority_plasma_review),
                    "mode": "0444",
                    "links": 1,
                },
            },
            "historical_f115_authority": {
                "evidence_sha256": sha256(historical_source_authority),
                "publication_audit_sha256": sha256(historical_source_authority_audit),
                "provenance_review_sha256": sha256(historical_provenance_review),
                "plasma_review_sha256": sha256(historical_plasma_review),
            },
            "source_archive_catalog": {
                "readme": {
                    "path": str(source_archive_readme),
                    "sha256": sha256(source_archive_readme),
                    "mode": "0644",
                    "links": 1,
                },
                "sha256sums": {
                    "path": str(source_archive_sums),
                    "sha256": sha256(source_archive_sums),
                    "mode": "0644",
                    "links": 1,
                },
                "bridge_bundle": {
                    "path": str(historical_source_bundle),
                    "sha256": sha256(historical_source_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": f113_revision,
                    "role": "retained-non-current-bridge",
                    "selected_as_current": False,
                },
                "current_source_bundle": {
                    "path": str(f116_source_bundle),
                    "sha256": sha256(f116_source_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": revision,
                    "selected_as_current": True,
                },
                "corrupt_c7_absent_from_active_checksum_ledger": True,
                "sole_current_source_bundle": str(f116_source_bundle),
            },
            "authority_and_enforcement": F116_AUTHORIZATION,
            "publication": (
                "recoverable-forward-transaction-with-publication-audit-commit-marker-"
                "under-stage-i-lock"
            ),
        },
    )
    f116_source_authority = source_authority
    f116_source_authority_provenance_review = source_authority_provenance_review
    f116_source_authority_plasma_review = source_authority_plasma_review
    f116_source_authority_audit = source_authority_audit
    f116_bindings = {
        "evidence": {
            "path": source_authority.relative_to(root).as_posix(),
            "sha256": sha256(source_authority),
        },
        "publication_audit": {
            "path": source_authority_audit.relative_to(root).as_posix(),
            "sha256": sha256(source_authority_audit),
        },
        "provenance_review": {
            "path": source_authority_provenance_review.relative_to(root).as_posix(),
            "sha256": sha256(source_authority_provenance_review),
        },
        "plasma_review": {
            "path": source_authority_plasma_review.relative_to(root).as_posix(),
            "sha256": sha256(source_authority_plasma_review),
        },
    }
    f116_digests = {
        "evidence_sha256": sha256(source_authority),
        "publication_audit_sha256": sha256(source_authority_audit),
        "provenance_review_sha256": sha256(source_authority_provenance_review),
        "plasma_review_sha256": sha256(source_authority_plasma_review),
    }
    f118_predecessor_bundle = dict(current_bundle)
    f118_predecessor_bundle.pop("candidate_path")
    f118_predecessor_bundle["selected_as_current"] = False
    f118_predecessor_bundle["role"] = "retained-non-current-predecessor"
    f118_current_bundle = {
        "candidate_path": str(source_bundle),
        "path": source_bundle.relative_to(root).as_posix(),
        "sha256": sha256(source_bundle),
        "complete_history": True,
        "head": revision,
        "advertised_tip": {"revision": revision, "name": "HEAD"},
        "verified_revisions": [f113_revision, revision],
        "selected_as_current": True,
        "subject": "Fixture current F118 source authority",
    }
    f118_committed_tools = f116_committed_tools(repository, revision)
    f118_publisher = next(
        item
        for item in f118_committed_tools
        if item["path"] == "scripts/frontier/cgl_lf_stage_i_source_authority.py"
    )
    source_authority = accounting / F118_NAME
    write_immutable_json(
        source_authority,
        {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-evidence",
            "checkpoint": "F-118",
            "execution_epoch": EPOCH,
            "generated_utc": (timestamp - timedelta(minutes=6)).isoformat(),
            "scope": {
                "relationship": "current-source-selection-only-supersession",
                "summary": "Select exact F118 fixture source without execution authority.",
                "preserves": F118_SCOPE_PRESERVES,
                "does_not_authorize": F118_SCOPE_DOES_NOT_AUTHORIZE,
            },
            "predecessor_authorities": {"historical_f116": f116_bindings},
            "implementation": {
                "publisher": f118_publisher,
                "committed_tools": f118_committed_tools,
                "intermediate_36140_bundle": bridge_bundle,
                "predecessor_current_source_bundle": f118_predecessor_bundle,
                "current_source_bundle": f118_current_bundle,
            },
            "source_archive_catalog": {
                "before": {
                    "readme_sha256": f116_catalog_after["readme_sha256"],
                    "sha256sums_sha256": f116_catalog_after["sha256sums_sha256"],
                    "bridge_listed_exactly_once": True,
                    "predecessor_current_source_bundle_listed_exactly_once": True,
                    "final_bundle_listed": False,
                    "corrupt_c7_listed": False,
                    "historical_f115_preserved": True,
                },
                "after": {
                    "readme_sha256": sha256(source_archive_readme),
                    "sha256sums_sha256": sha256(source_archive_sums),
                    "bridge_listed_exactly_once": True,
                    "predecessor_current_source_bundle_listed_exactly_once": True,
                    "final_bundle_listed_exactly_once": True,
                    "corrupt_c7_listed": False,
                    "historical_f115_preserved": True,
                    "historical_f116_preserved": True,
                    "all_prior_checksum_entries_preserved": True,
                    "sole_current_source_bundle": source_bundle.relative_to(
                        root
                    ).as_posix(),
                },
            },
            "authorization": F116_AUTHORIZATION,
            "validation": F118_VALIDATION_CLAIMS,
            "publication_requirements": F118_PUBLICATION_REQUIREMENTS,
        },
    )
    f118_verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "predecessor_current_source_bundle_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": sha256(source_bundle),
        "final_head": revision,
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
    }
    source_authority_provenance_review = accounting / f"{F118_NAME}.provenance_security_review.json"
    source_authority_plasma_review = accounting / f"{F118_NAME}.plasma_scientific_review.json"
    for path, kind, decision, agent in (
        (
            source_authority_provenance_review,
            "provenance-security",
            "approved-for-publication",
            "fixture-f118-provenance-reviewer",
        ),
        (
            source_authority_plasma_review,
            "plasma-scientific-continuation",
            "approved",
            "fixture-f118-plasma-reviewer",
        ),
    ):
        write_immutable_json(
            path,
            {
                "schema_version": 1,
                "record_type": "stage-i-current-source-authority-supersession-independent-review",
                "checkpoint": "F-118",
                "execution_epoch": EPOCH,
                "review_kind": kind,
                "decision": decision,
                "reviewed_candidate": {
                    "path": str(source_authority.with_suffix(".json.candidate")),
                    "sha256": sha256(source_authority),
                },
                "published_f118": {
                    "path": str(source_authority),
                    "sha256": sha256(source_authority),
                },
                "reviewer": {"agent_id": agent, "identity": f"fixture {kind} reviewer"},
                "reviewed_utc": (timestamp - timedelta(minutes=5)).isoformat(),
                "findings": ["Exact current F118 source-only authority verified."],
                "limitations": [
                    "No prepare, submit, scheduler, or science authority.",
                    INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION,
                ],
                "verified": f118_verified,
            },
        )
    source_authority_audit = accounting / f"{F118_NAME}.publication_audit.json"
    write_immutable_json(
        source_authority_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-current-source-authority-supersession-publication-audit",
            "checkpoint": "F-118",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=4)).isoformat(),
            "artifact": {
                "path": str(source_authority),
                "sha256": sha256(source_authority),
                "mode": "0444",
                "links": 1,
            },
            "independent_reviews": {
                "reviews_bind_exact_published_f118_sha256": sha256(source_authority),
                "provenance_security": {
                    "path": str(source_authority_provenance_review),
                    "sha256": sha256(source_authority_provenance_review),
                    "mode": "0444",
                    "links": 1,
                },
                "plasma_scientific_continuation": {
                    "path": str(source_authority_plasma_review),
                    "sha256": sha256(source_authority_plasma_review),
                    "mode": "0444",
                    "links": 1,
                },
            },
            "historical_f116_authority": f116_digests,
            "source_archive_catalog": {
                "readme": {
                    "path": str(source_archive_readme),
                    "sha256": sha256(source_archive_readme),
                    "mode": "0644",
                    "links": 1,
                },
                "sha256sums": {
                    "path": str(source_archive_sums),
                    "sha256": sha256(source_archive_sums),
                    "mode": "0644",
                    "links": 1,
                },
                "bridge_bundle": {
                    "path": str(historical_source_bundle),
                    "sha256": sha256(historical_source_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": f113_revision,
                    "role": "retained-non-current-bridge",
                    "selected_as_current": False,
                },
                "predecessor_current_source_bundle": {
                    "path": str(f116_source_bundle),
                    "sha256": sha256(f116_source_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": revision,
                    "role": "retained-non-current-predecessor",
                    "selected_as_current": False,
                },
                "current_source_bundle": {
                    "path": str(source_bundle),
                    "sha256": sha256(source_bundle),
                    "mode": "0644",
                    "links": 1,
                    "head": revision,
                    "selected_as_current": True,
                },
                "corrupt_c7_absent_from_active_checksum_ledger": True,
                "sole_current_source_bundle": str(source_bundle),
            },
            "authority_and_enforcement": F116_AUTHORIZATION,
            "publication": F118_PUBLICATION_METHOD,
        },
    )
    predecessor = accounting / "mks24_stage_i_E03_forcing_policy_F199_recost_evidence.json"
    write_immutable_json(
        predecessor,
        {
            "schema_version": 2,
            "record_type": "stage-i-recost-recommendation-evidence",
            "checkpoint": "F-199",
            "artifact_name": predecessor.name,
            "execution_epoch": EPOCH,
            "generated_utc": (timestamp - timedelta(minutes=8)).isoformat(),
            "authority": {
                "authorizing": False,
                "action_authority": "none-until-independent-review-and-publication",
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            },
        },
    )
    predecessor_review = accounting / f"{predecessor.name}.independent_review.json"
    write_immutable_json(
        predecessor_review,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (timestamp - timedelta(minutes=7)).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": {
                "agent_id": "fixture-predecessor-reviewer",
                "independent_from_generator": True,
            },
            "candidate": {"path": str(predecessor), "sha256": sha256(predecessor)},
            "scope": {"non_authorizing": True},
        },
    )
    predecessor_audit = accounting / f"{predecessor.name}.publication_audit.json"
    write_immutable_json(
        predecessor_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=5)).isoformat(),
            "artifact": {
                "path": str(predecessor),
                "sha256": sha256(predecessor),
                "mode": "0444",
                "links": 1,
            },
            "independent_review": {
                "path": str(predecessor_review),
                "sha256": sha256(predecessor_review),
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
    next_profile = profile(fixture)
    projected_storage, storage_projections = computed_storage_projections(
        fixture, [next_profile]
    )
    storage = accounting / "reviewed_storage_evidence.json"
    write_json(
        storage,
        {
            "schema_version": 1,
            "record_type": "stage-i-storage-evidence",
            "execution_epoch": EPOCH,
            "root": str(root),
            "measured_utc": timestamp.isoformat(),
            "available_bytes": max(20_000, min(live_available(root) // 4, 2_000_000_000_000)),
            "retained_stage_i_bytes": retained_bytes(root),
            "required_safety_bytes": 1024,
            "projected_authorized_wave_growth_bytes": projected_storage,
            "projection_method": STORAGE_PROJECTION_METHOD,
            "profile_projections_sha256": hashlib.sha256(
                (json.dumps(storage_projections, sort_keys=True) + "\n").encode()
            ).hexdigest(),
        },
    )
    queue = tmp_path / "squeue.txt"
    queue.write_text("")
    queue.chmod(0o644)
    output = accounting / "mks24_stage_i_E03_forcing_policy_F200_recost_evidence.json.staged"
    fixture.update(
        {
            "ledger": ledger,
            "reservations": reservations,
            "reconcile": reconcile,
            "scheduler": scheduler,
            "ceiling": ceiling,
            "ceiling_audit": ceiling_audit,
            "source_authority": source_authority,
            "source_authority_audit": source_authority_audit,
            "source_authority_provenance_review": source_authority_provenance_review,
            "source_authority_plasma_review": source_authority_plasma_review,
            "f116_source_authority": f116_source_authority,
            "f116_source_authority_audit": f116_source_authority_audit,
            "f116_source_authority_provenance_review": (
                f116_source_authority_provenance_review
            ),
            "f116_source_authority_plasma_review": f116_source_authority_plasma_review,
            "historical_source_authority": historical_source_authority,
            "historical_source_authority_audit": historical_source_authority_audit,
            "historical_source_authority_provenance_review": historical_provenance_review,
            "historical_source_authority_plasma_review": historical_plasma_review,
            "predecessor": predecessor,
            "predecessor_review": predecessor_review,
            "predecessor_audit": predecessor_audit,
            "storage": storage,
            "queue": queue,
            "output": output,
        }
    )
    request = accounting / "reviewed_recost_request.json"
    request_review = accounting / "reviewed_recost_request.json.independent_review.json"
    fixture["request"] = request
    fixture["request_review"] = request_review
    write_json(
        request,
        {
            "schema_version": 2,
            "record_type": "stage-i-recost-recommendation-request",
            "checkpoint": "F-200",
            "artifact_name": output.name.removesuffix(".staged"),
            "execution_epoch": EPOCH,
            "generated_utc": timestamp.isoformat(),
            "expires_utc": (timestamp + timedelta(hours=12)).isoformat(),
            "requested_by": "fixture-request-author",
            "scope": "Fixture just-recorded clean-partial recost barrier.",
            "barrier": {
                "recorded_segments": [
                    {
                        "case_id": "R03",
                        "segment": "s00_rankio_t0_t0p5",
                        "job_id": "12345",
                        "result": "clean_partial",
                    }
                ]
            },
            "inputs": {
                "reconciliation": {
                    "path": str(reconcile.relative_to(root)),
                    "sha256": sha256(reconcile),
                },
                "ledger": {"path": str(ledger.relative_to(root)), "sha256": sha256(ledger)},
                "reservations": {
                    "path": str(reservations.relative_to(root)),
                    "sha256": sha256(reservations),
                },
                "manifests": [
                    {"path": str(manifest.relative_to(root)), "sha256": sha256(manifest)}
                ],
                "scheduler_evidence": [
                    {"path": str(scheduler.relative_to(root)), "sha256": sha256(scheduler)}
                ],
                "storage_evidence": {
                    "path": str(storage.relative_to(root)),
                    "sha256": sha256(storage),
                },
                "source_bundle": {
                    "path": str(source_bundle.relative_to(root)),
                    "sha256": sha256(source_bundle),
                    "verified_revisions": [f113_revision, revision],
                },
                "matrix": {
                    "path": str(matrix.relative_to(repository)),
                    "revision": revision,
                    "sha256": sha256(matrix),
                },
                "stage_i_helper": {"revision": revision, "sha256": sha256(helper)},
                "ceiling_evidence": {
                    "path": str(ceiling.relative_to(root)),
                    "sha256": sha256(ceiling),
                },
                "ceiling_publication_audit": {
                    "path": str(ceiling_audit.relative_to(root)),
                    "sha256": sha256(ceiling_audit),
                },
                "source_authority": {
                    "checkpoint": "F-118",
                    "evidence": {
                        "path": str(source_authority.relative_to(root)),
                        "sha256": sha256(source_authority),
                    },
                    "publication_audit": {
                        "path": str(source_authority_audit.relative_to(root)),
                        "sha256": sha256(source_authority_audit),
                    },
                    "provenance_review": {
                        "path": str(source_authority_provenance_review.relative_to(root)),
                        "sha256": sha256(source_authority_provenance_review),
                    },
                    "plasma_review": {
                        "path": str(source_authority_plasma_review.relative_to(root)),
                        "sha256": sha256(source_authority_plasma_review),
                    },
                    "final_source_bundle": {
                        "path": str(source_bundle.relative_to(root)),
                        "sha256": sha256(source_bundle),
                        "verified_revisions": [f113_revision, revision],
                    },
                },
                "qualification_approval": {
                    "path": str(qualification.relative_to(root)),
                    "sha256": sha256(qualification),
                },
                "predecessor_recost": {
                    "path": str(predecessor.relative_to(root)),
                    "sha256": sha256(predecessor),
                },
                "predecessor_recost_independent_review": {
                    "path": str(predecessor_review.relative_to(root)),
                    "sha256": sha256(predecessor_review),
                },
                "predecessor_recost_publication_audit": {
                    "path": str(predecessor_audit.relative_to(root)),
                    "sha256": sha256(predecessor_audit),
                },
                "r17_readiness_evidence": None,
                "r17_readiness_independent_review": None,
                "r17_readiness_publication_audit": None,
            },
            "recommendations": {
                "mode": "sole-next-profile",
                "max_wave_nodes": 1,
                "profiles": [next_profile],
            },
        },
    )
    refresh_request(fixture)
    return fixture


def test_generator_emits_deterministic_checkpoint_compatible_sole_artifact(recost_fixture):
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    output = recost_fixture["output"]
    assert isinstance(output, Path)
    first = output.read_bytes()
    artifact = json.loads(first)
    assert artifact["artifact_name"] == output.name.removesuffix(".staged")
    assert artifact["checkpoint"] == "F-200"
    assert artifact["authority"]["authorizing"] is False
    assert artifact["recommendations"]["authorizing"] is False
    assert set(artifact["recommendations"]["sole_next_segment_recommendation"]) == SOLE_COMPATIBLE_KEYS
    assert (
        artifact["recommendations"]["sole_next_segment_recommendation"]
        != artifact["recommendations"]["recommended_next_profiles"][0]
    )
    assert artifact["predecessor_recost"]["checkpoint"] == "F-199"
    assert artifact["barrier"]["job_ids"] == ["12345"]
    assert artifact["budget"]["method"] == "observed-stage-i-scoped-node-hour-rate-v2"
    assert artifact["budget"]["measurement_basis"]["job_id"] == "12345"
    assert float(artifact["budget"]["computed_stage_i_total_node_hours"]) > 500
    assert (
        artifact["budget"]["case_breakdown"]["R02"][
            "matrix_full_case_node_hours_reference_only"
        ]
        == "10.0"
    )
    assert (
        artifact["budget"]["storage_projection_evidence"]["method"]
        == STORAGE_PROJECTION_METHOD
    )
    source_authority = artifact["provenance"]["source_authority"]
    assert set(source_authority) == {
        "checkpoint",
        "evidence",
        "provenance_review",
        "plasma_review",
        "publication_audit",
        "final_source_bundle",
    }
    assert source_authority == json.loads(recost_fixture["request"].read_text())[
        "inputs"
    ]["source_authority"]
    assert source_authority["final_source_bundle"] == json.loads(
        recost_fixture["request"].read_text()
    )["inputs"]["source_bundle"]
    assert artifact["provenance"]["source_bundle_sha256"] != sha256(
        recost_fixture["historical_source_bundle"]
    )
    assert artifact["promoted_f113"]["publication_audit"]["review"]["status"] == "approved"
    assert (
        artifact["provenance"]["f113_historical_helper_revision"]
        == recost_fixture["f113_revision"]
    )
    assert (
        artifact["provenance"]["f113_historical_helper_sha256"]
        == recost_fixture["f113_helper_sha256"]
    )
    assert (
        artifact["provenance"]["f113_historical_helper_revision"]
        != artifact["provenance"]["stage_i_helper_revision"]
    )
    assert (
        artifact["provenance"]["f113_historical_helper_sha256"]
        != artifact["provenance"]["stage_i_helper_sha256"]
    )
    assert artifact["provenance"]["scheduler_sha256"] == sha256(recost_fixture["scheduler"])
    assert output.stat().st_mode & 0o777 == 0o444
    assert output.stat().st_nlink == 1
    verified = run_generator(recost_fixture)
    assert verified.returncode == 0, verified.stderr
    assert output.read_bytes() == first
    output.unlink()
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    assert output.read_bytes() == first


def test_draft_request_action_generates_authenticated_prerequisites_without_self_review(
    recost_fixture,
):
    packet = write_draft_packet(recost_fixture)
    completed = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert completed.returncode == 0, completed.stderr
    report = action_report(completed)
    assert report["independent_review_created"] is False
    assert report["self_approved"] is False
    assert report["managed_independent_review_install"]["command_template"][-1] == (
        "install-request-review"
    )
    assert report["staged_generation_after_external_review"][-1] == str(
        recost_fixture["queue"]
    )

    accounting = recost_fixture["accounting"]
    assert isinstance(accounting, Path)
    prefix = f"mks24_stage_i_{EPOCH_SLUG}_F201"
    request = accounting / f"{prefix}_recost_request.json"
    reconciliation = accounting / f"{prefix}_reconciliation_evidence.json"
    storage = accounting / f"{prefix}_storage_evidence.json"
    review = request.with_name(f"{request.name}.independent_review.json")
    assert all(path.parent == accounting for path in (request, reconciliation, storage))
    assert all(stat.S_IMODE(path.stat().st_mode) == 0o644 for path in (request, reconciliation, storage))
    assert not review.exists()
    storage_value = json.loads(storage.read_text())
    assert storage_value["available_bytes"] == (
        storage_value["required_safety_bytes"]
        + storage_value["projected_authorized_wave_growth_bytes"]
    )

    value = json.loads(request.read_text())
    assert value["schema_version"] == 2
    assert value["record_type"] == "stage-i-recost-recommendation-request"
    assert len(value["inputs"]["manifests"]) == 1
    assert len(value["inputs"]["scheduler_evidence"]) == 1
    assert value["inputs"]["source_authority"]["checkpoint"] == "F-118"
    assert value["inputs"]["matrix"]["path"] == "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
    assert value["inputs"]["stage_i_helper"]["revision"] == recost_fixture["revision"]
    assert value["inputs"]["qualification_approval"]["path"] == (
        f"accounting/{QUALIFICATION_NAME}"
    )
    assert value["inputs"]["predecessor_recost"]["path"].endswith(
        "_F199_recost_evidence.json"
    )

    draft_fixture = dict(recost_fixture)
    draft_fixture["request"] = request
    draft_fixture["request_sha256"] = sha256(request)
    draft_fixture["output"] = accounting / f"{prefix}_recost_evidence.json.staged"
    assert_rejected(
        run_generator(draft_fixture),
        "recost_request.json.independent_review.json",
    )
    external_review = write_external_request_review(recost_fixture, request)
    installed = run_action(
        recost_fixture,
        "install-request-review",
        "--request",
        str(request),
        "--expected-request-sha256",
        sha256(request),
        "--review",
        str(external_review),
        "--expected-review-sha256",
        sha256(external_review),
    )
    assert installed.returncode == 0, installed.stderr
    install_report = action_report(installed)
    assert install_report["created"] is True
    assert install_report["staged_generation"] == report[
        "staged_generation_after_external_review"
    ]
    assert review.read_bytes() == external_review.read_bytes()
    assert stat.S_IMODE(review.stat().st_mode) == 0o444
    generated = run_generator(draft_fixture)
    assert generated.returncode == 0, generated.stderr


def test_install_request_review_rejects_drift_symlinks_and_canonical_sources(
    recost_fixture,
):
    packet = write_draft_packet(recost_fixture, checkpoint_number=205)
    drafted = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert drafted.returncode == 0, drafted.stderr
    accounting = recost_fixture["accounting"]
    root = recost_fixture["root"]
    assert isinstance(accounting, Path)
    assert isinstance(root, Path)
    request = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F205_recost_request.json"
    target = request.with_name(f"{request.name}.independent_review.json")
    external = write_external_request_review(recost_fixture, request)
    arguments = (
        "--request",
        str(request),
        "--expected-request-sha256",
        sha256(request),
        "--review",
        str(external),
        "--expected-review-sha256",
        sha256(external),
    )
    lock = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    with lock.open("r+") as descriptor:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert_rejected(
            run_action(recost_fixture, "install-request-review", *arguments),
            "another Stage I mutation holds",
        )
        fcntl.flock(descriptor, fcntl.LOCK_UN)
    assert_rejected(
        run_action(
            recost_fixture,
            "install-request-review",
            *arguments[:-1],
            "0" * 64,
        ),
        "checksum differs",
    )
    linked = root.parent / "linked-external-review.json"
    linked.symlink_to(external)
    linked_arguments = list(arguments)
    linked_arguments[5] = str(linked)
    assert_rejected(
        run_action(recost_fixture, "install-request-review", *linked_arguments),
        "path contains a symlink",
    )
    canonical_source = accounting / "manual-review-staging/review-candidate.json"
    write_immutable_json(canonical_source, json.loads(external.read_text()))
    canonical_arguments = list(arguments)
    canonical_arguments[5] = str(canonical_source)
    assert_rejected(
        run_action(recost_fixture, "install-request-review", *canonical_arguments),
        "external immutable candidate",
    )
    target.symlink_to(external)
    target_symlink = run_action(recost_fixture, "install-request-review", *arguments)
    assert target_symlink.returncode == 1
    assert target.is_symlink()
    target.unlink()

    installed = run_action(recost_fixture, "install-request-review", *arguments)
    assert installed.returncode == 0, installed.stderr
    assert target.read_bytes() == external.read_bytes()
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    verified = run_action(recost_fixture, "install-request-review", *arguments)
    assert verified.returncode == 0, verified.stderr
    assert action_report(verified)["exact_existing_copy_verified"] is True

    target.chmod(0o644)
    target.write_bytes(target.read_bytes() + b"\n")
    target.chmod(0o444)
    retained = target.read_bytes()
    assert_rejected(
        run_action(recost_fixture, "install-request-review", *arguments),
        "exists with different bytes; refusing to clobber",
    )
    assert target.read_bytes() == retained


def test_draft_request_action_rejects_overauthorization_and_never_clobbers(
    recost_fixture,
):
    overauthorized = write_draft_packet(
        recost_fixture,
        checkpoint_number=202,
        mutate=lambda packet: packet["draft_policy"].update({"self_approved": True}),
    )
    assert_rejected(
        run_action(
            recost_fixture,
            "draft-request",
            "--packet",
            str(overauthorized),
            "--expected-packet-sha256",
            sha256(overauthorized),
        ),
        "over-authorizes",
    )
    accounting = recost_fixture["accounting"]
    assert isinstance(accounting, Path)
    assert not (accounting / f"mks24_stage_i_{EPOCH_SLUG}_F202_recost_request.json").exists()

    packet = write_draft_packet(recost_fixture, checkpoint_number=203)
    first = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert first.returncode == 0, first.stderr
    request = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F203_recost_request.json"
    request.write_text(request.read_text() + "\n")
    retained = request.read_bytes()
    assert_rejected(
        run_action(
            recost_fixture,
            "draft-request",
            "--packet",
            str(packet),
            "--expected-packet-sha256",
            sha256(packet),
        ),
        "exists with different bytes; refusing to clobber",
    )
    assert request.read_bytes() == retained


def test_install_f117_draft_packet_is_managed_no_clobber_or_exact_verify(
    recost_fixture,
):
    packet = write_draft_packet(recost_fixture, checkpoint_number=117)
    source = packet.parent.parent / "reviewed-f117-draft-packet.json"
    packet.replace(source)
    assert not packet.exists()
    packet.write_bytes(b"unrelated owner-only target\n")
    packet.chmod(0o600)
    unrelated = packet.stat()
    rejected = run_action(
        recost_fixture,
        "install-f117-draft-packet",
        "--packet",
        str(source),
        "--expected-packet-sha256",
        sha256(source),
    )
    assert_rejected(rejected, "mode is 0600, expected 0644")
    assert packet.read_bytes() == b"unrelated owner-only target\n"
    assert (packet.stat().st_dev, packet.stat().st_ino) == (
        unrelated.st_dev,
        unrelated.st_ino,
    )
    assert stat.S_IMODE(packet.stat().st_mode) == 0o600
    packet.unlink()

    created = run_action(
        recost_fixture,
        "install-f117-draft-packet",
        "--packet",
        str(source),
        "--expected-packet-sha256",
        sha256(source),
    )
    assert created.returncode == 0, created.stderr
    report = action_report(created)
    assert report["created"] is True
    assert Path(report["path"]) == packet
    assert packet.read_bytes() == source.read_bytes()
    assert stat.S_IMODE(packet.stat().st_mode) == 0o644
    assert report["next_action"][-1] == "draft-request"

    verified = run_action(
        recost_fixture,
        "install-f117-draft-packet",
        "--packet",
        str(source),
        "--expected-packet-sha256",
        sha256(source),
    )
    assert verified.returncode == 0, verified.stderr
    assert action_report(verified)["exact_existing_copy_verified"] is True

    packet.write_bytes(packet.read_bytes() + b"\n")
    packet.chmod(0o644)
    drift = packet.read_bytes()
    assert_rejected(
        run_action(
            recost_fixture,
            "install-f117-draft-packet",
            "--packet",
            str(source),
            "--expected-packet-sha256",
            sha256(source),
        ),
        "exists with different bytes; refusing to clobber",
    )
    assert packet.read_bytes() == drift

    wrong = write_draft_packet(recost_fixture, checkpoint_number=118)
    assert_rejected(
        run_action(
            recost_fixture,
            "install-f117-draft-packet",
            "--packet",
            str(wrong),
            "--expected-packet-sha256",
            sha256(wrong),
        ),
        "restricted to exact F-117",
    )


def test_f119_managed_install_preserves_published_f117(recost_fixture):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    f117_publications = replace_fixture_predecessor(recost_fixture, 117)
    f117_packet, f119_packet = write_f119_seed_and_candidate(
        recost_fixture,
        seed_mutate=lambda packet: bind_fixture_predecessor(recost_fixture, packet),
    )
    external_f119 = root.parent / "externally-reviewed-f119-draft-packet.json"
    f119_packet.replace(external_f119)
    retained_f117 = (f117_packet, *f117_publications)
    f117_snapshot = {
        path: (
            path.read_bytes(),
            path.stat().st_dev,
            path.stat().st_ino,
            stat.S_IMODE(path.stat().st_mode),
        )
        for path in retained_f117
    }

    installed = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external_f119),
        "--expected-packet-sha256",
        sha256(external_f119),
    )
    assert installed.returncode == 0, installed.stderr
    install_report = action_report(installed)
    assert install_report["action"] == "install-f119-draft-packet"
    assert Path(install_report["path"]) == f119_packet
    assert install_report["next_action"][-1] == "draft-request"

    for path, snapshot in f117_snapshot.items():
        assert (
            path.read_bytes(),
            path.stat().st_dev,
            path.stat().st_ino,
            stat.S_IMODE(path.stat().st_mode),
        ) == snapshot

    wrong = write_draft_packet(recost_fixture, checkpoint_number=118)
    assert_rejected(
        run_action(
            recost_fixture,
            "install-f119-draft-packet",
            "--packet",
            str(wrong),
            "--expected-packet-sha256",
            sha256(wrong),
        ),
        "restricted to exact F-119",
    )


def test_f119_candidate_renderer_is_deterministic_read_only_and_verifiable(
    recost_fixture,
):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    arguments = f119_render_arguments(recost_fixture, seed)
    first = run_action(recost_fixture, "render-f119-draft-packet", *arguments)
    second = run_action(recost_fixture, "render-f119-draft-packet", *arguments)
    assert first.returncode == 0, first.stderr
    assert second.returncode == 0, second.stderr
    assert first.stdout == second.stdout
    assert not candidate.exists()
    value = json.loads(first.stdout)
    assert value["checkpoint"] == "F-119"
    assert value["inputs"]["source_authority"]["checkpoint"] == "F-118"
    assert value["inputs"]["ceiling_publication_audit"] == {
        "path": recost_fixture["source_authority_audit"].relative_to(root).as_posix(),
        "sha256": sha256(recost_fixture["source_authority_audit"]),
    }
    assert value["inputs"]["stage_i_helper"] == {
        "revision": recost_fixture["revision"],
        "sha256": sha256(recost_fixture["helper"]),
    }
    assert value["draft_policy"]["required_storage_safety_bytes"] == 1024**4
    assert [
        (profile_value["case_id"], profile_value["nodes"])
        for profile_value in value["recommendations"]["profiles"]
    ] == [("R03", 1), ("R04", 4), ("R12", 4), ("R16", 1)]

    external = root.parent / "rendered-f119-candidate.json"
    external.write_text(first.stdout)
    external.chmod(0o644)
    verified = run_action(
        recost_fixture,
        "verify-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert verified.returncode == 0, verified.stderr
    assert verified.stdout.strip() == sha256(external)


def test_f119_render_verify_install_and_draft_request_end_to_end(
    recost_fixture, monkeypatch
):
    module = load_recost_module()
    root = recost_fixture["root"]
    accounting = recost_fixture["accounting"]
    repository = recost_fixture["repository"]
    generator = recost_fixture["generator"]
    queue = recost_fixture["queue"]
    assert isinstance(root, Path)
    assert isinstance(accounting, Path)
    assert isinstance(repository, Path)
    assert isinstance(generator, Path)
    assert isinstance(queue, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    external = root.parent / "rendered-f119-end-to-end.json"
    external.write_text(rendered.stdout)
    external.chmod(0o644)
    verified = run_action(
        recost_fixture,
        "verify-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert verified.returncode == 0, verified.stderr
    installed = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert installed.returncode == 0, installed.stderr
    real_fstatvfs = module.os.fstatvfs

    def f119_storage_fixture(descriptor):
        retained = list(real_fstatvfs(descriptor))
        retained[3] = max(retained[3], 2 * 1024**4 // retained[1])
        retained[4] = max(retained[4], 2 * 1024**4 // retained[1])
        return os.statvfs_result(retained)

    monkeypatch.setattr(module.os, "fstatvfs", f119_storage_fixture)
    args = module.parse_args(
        [
            "--root",
            str(root),
            "--allow-local-root",
            "--packet",
            str(candidate),
            "--expected-packet-sha256",
            sha256(candidate),
            "--expected-generator-sha256",
            sha256(generator),
            "--squeue-file",
            str(queue),
            "draft-request",
        ]
    )
    with module.stage_i_lock(root) as mutation_lock:
        module.locked_draft_request(
            args,
            root,
            generator,
            repository,
            sha256(generator),
            mutation_lock,
        )
    request = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F119_recost_request.json"
    reconciliation = (
        accounting / f"mks24_stage_i_{EPOCH_SLUG}_F119_reconciliation_evidence.json"
    )
    storage = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F119_storage_evidence.json"
    assert request.read_bytes()
    assert reconciliation.read_bytes()
    assert storage.read_bytes()


@pytest.mark.parametrize(
    ("action", "serialization"),
    (
        ("verify-f119-draft-packet", "compact"),
        ("install-f119-draft-packet", "compact"),
        ("verify-f119-draft-packet", "trailing-space"),
        ("install-f119-draft-packet", "trailing-space"),
    ),
)
def test_f119_verify_and_install_require_byte_exact_stable_json(
    recost_fixture, action, serialization
):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    value = json.loads(rendered.stdout)
    external = root.parent / f"{serialization}-f119-candidate.json"
    if serialization == "compact":
        external.write_text(json.dumps(value, sort_keys=True) + "\n")
    else:
        external.write_text(rendered.stdout.rstrip("\n") + " \n")
    external.chmod(0o644)
    rejected = run_action(
        recost_fixture,
        action,
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert_rejected(rejected, "not byte-exact stable JSON")
    assert not candidate.exists()


@pytest.mark.parametrize(
    "mutation", ("ceiling-publication-audit", "stage-i-helper-revision", "stage-i-helper-sha")
)
def test_f119_rejects_obsolete_or_mutated_f118_runtime_bindings(
    recost_fixture, mutation
):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    value = json.loads(rendered.stdout)
    if mutation == "ceiling-publication-audit":
        value["inputs"]["ceiling_publication_audit"] = {
            "path": recost_fixture["ceiling_audit"].relative_to(root).as_posix(),
            "sha256": sha256(recost_fixture["ceiling_audit"]),
        }
    elif mutation == "stage-i-helper-revision":
        value["inputs"]["stage_i_helper"]["revision"] = recost_fixture["f113_revision"]
    else:
        value["inputs"]["stage_i_helper"]["sha256"] = recost_fixture[
            "f113_helper_sha256"
        ]
    external = root.parent / f"mutated-{mutation}-f119-candidate.json"
    write_json(external, value)
    rejected = run_action(
        recost_fixture,
        "verify-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert_rejected(rejected, "F118-era runtime binding differs")
    assert not candidate.exists()


@pytest.mark.parametrize(
    "mutation",
    (
        "missing-input",
        "extra-input",
        "empty-profiles",
        "wrong-profiles",
        "f116-authority",
        "wrong-identity",
        "wrong-predecessor",
        "wrong-reserve",
        "expired",
    ),
)
def test_f119_install_rejects_semantic_mutations_before_canonical_write(
    recost_fixture, mutation,
):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    value = json.loads(rendered.stdout)
    if mutation == "missing-input":
        value["inputs"].pop("matrix")
    elif mutation == "extra-input":
        value["inputs"]["unexpected"] = None
    elif mutation == "empty-profiles":
        value["recommendations"]["profiles"] = []
    elif mutation == "wrong-profiles":
        value["recommendations"]["profiles"][0]["nodes"] = 2
    elif mutation == "f116-authority":
        value["inputs"]["source_authority"]["checkpoint"] = "F-116"
    elif mutation == "wrong-identity":
        value["checkpoint"] = "F-120"
    elif mutation == "wrong-predecessor":
        value["inputs"]["predecessor_recost"]["sha256"] = "0" * 64
    elif mutation == "wrong-reserve":
        value["draft_policy"]["required_storage_safety_bytes"] = 1024
    else:
        value["generated_utc"] = (
            datetime.now(timezone.utc) - timedelta(hours=2)
        ).isoformat()
        value["expires_utc"] = (
            datetime.now(timezone.utc) - timedelta(hours=1)
        ).isoformat()
    external = root.parent / f"mutated-{mutation}-f119-candidate.json"
    write_json(external, value)
    rejected = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert rejected.returncode == 1
    assert not candidate.exists()
    assert not list(candidate.parent.glob(f".{candidate.name}.retired.sha256-*"))


def test_f119_install_authenticates_and_retires_stale_packet(recost_fixture):
    root = recost_fixture["root"]
    timestamp = recost_fixture["timestamp"]
    assert isinstance(root, Path)
    assert isinstance(timestamp, datetime)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    stale = json.loads(rendered.stdout)
    stale["generated_utc"] = (timestamp - timedelta(minutes=4)).isoformat()
    stale["expires_utc"] = (timestamp - timedelta(seconds=1)).isoformat()
    write_json(candidate, stale)
    stale_payload = candidate.read_bytes()
    external = root.parent / "active-f119-candidate.json"
    external.write_text(rendered.stdout)
    external.chmod(0o644)

    installed = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert installed.returncode == 0, installed.stderr
    report = action_report(installed)
    retired = Path(report["authenticated_retired_predecessor"])
    assert retired.read_bytes() == stale_payload
    assert candidate.read_bytes() == external.read_bytes()


def test_f119_install_never_retires_active_unexpired_valid_packet(recost_fixture):
    root = recost_fixture["root"]
    timestamp = recost_fixture["timestamp"]
    assert isinstance(root, Path)
    assert isinstance(timestamp, datetime)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    first = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(
            recost_fixture, seed, expires=timestamp + timedelta(hours=11)
        ),
    )
    second = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(
            recost_fixture, seed, expires=timestamp + timedelta(hours=12)
        ),
    )
    assert first.returncode == 0, first.stderr
    assert second.returncode == 0, second.stderr
    candidate.write_text(first.stdout)
    candidate.chmod(0o644)
    retained = candidate.read_bytes()
    external = root.parent / "different-active-f119-candidate.json"
    external.write_text(second.stdout)
    external.chmod(0o644)
    rejected = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert_rejected(rejected, "active unexpired valid F119 draft packet cannot be retired")
    assert candidate.read_bytes() == retained
    assert not list(candidate.parent.glob(f".{candidate.name}.retired.sha256-*"))


@pytest.mark.parametrize(
    "namespace", ("request-only", "reconciliation-only", "storage-only", "pre-request-pair")
)
def test_f119_install_rejects_committed_or_partial_request_namespace(
    recost_fixture, namespace
):
    root = recost_fixture["root"]
    accounting = recost_fixture["accounting"]
    assert isinstance(root, Path)
    assert isinstance(accounting, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    external = root.parent / f"{namespace}-f119-candidate.json"
    external.write_text(rendered.stdout)
    external.chmod(0o644)
    prefix = f"mks24_stage_i_{EPOCH_SLUG}_F119"
    paths = {
        "request": accounting / f"{prefix}_recost_request.json",
        "reconciliation": accounting / f"{prefix}_reconciliation_evidence.json",
        "storage": accounting / f"{prefix}_storage_evidence.json",
    }
    selected = {
        "request-only": ("request",),
        "reconciliation-only": ("reconciliation",),
        "storage-only": ("storage",),
        "pre-request-pair": ("reconciliation", "storage"),
    }[namespace]
    for key in selected:
        paths[key].write_bytes(f"hostile {key}\n".encode())
        paths[key].chmod(0o644)
    snapshots = {key: paths[key].read_bytes() for key in selected}
    rejected = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    if namespace == "request-only":
        assert_rejected(rejected, "request is the last commit marker")
    else:
        assert_rejected(rejected, "draft prerequisite namespace is partial")
    assert not candidate.exists()
    assert {key: paths[key].read_bytes() for key in selected} == snapshots


def test_f119_retirement_post_barrier_rejects_concurrent_request_marker(
    recost_fixture, monkeypatch
):
    module = load_recost_module()
    root = recost_fixture["root"]
    repository = recost_fixture["repository"]
    timestamp = recost_fixture["timestamp"]
    accounting = recost_fixture["accounting"]
    assert isinstance(root, Path)
    assert isinstance(repository, Path)
    assert isinstance(timestamp, datetime)
    assert isinstance(accounting, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    stale = json.loads(rendered.stdout)
    stale["generated_utc"] = (timestamp - timedelta(minutes=4)).isoformat()
    stale["expires_utc"] = (timestamp - timedelta(seconds=1)).isoformat()
    write_json(candidate, stale)
    stale_payload = candidate.read_bytes()
    request = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F119_recost_request.json"
    real_move = module.move_bound_name_noreplace
    raced = False

    def create_request_after_retirement(*args, **kwargs):
        nonlocal raced
        moved = real_move(*args, **kwargs)
        if not raced and args[1] == candidate.name:
            request.write_bytes(b"concurrent request commit marker\n")
            request.chmod(0o644)
            raced = True
        return moved

    monkeypatch.setattr(module, "move_bound_name_noreplace", create_request_after_retirement)
    with module.stage_i_lock(root) as mutation_lock:
        with pytest.raises(ValueError, match="request is the last commit marker"):
            module.retire_stale_f119_packet(
                root,
                repository,
                candidate,
                rendered.stdout.encode(),
                mutation_lock,
            )
    assert raced
    assert request.read_bytes() == b"concurrent request commit marker\n"
    assert not candidate.exists()
    retired = list(accounting.glob(f".{candidate.name}.retired.sha256-*"))
    assert len(retired) == 1
    assert retired[0].read_bytes() == stale_payload


@pytest.mark.parametrize("packet_state", ("absent", "exact"))
def test_f119_retirement_early_return_post_barrier_rejects_request_race(
    recost_fixture, monkeypatch, packet_state
):
    module = load_recost_module()
    root = recost_fixture["root"]
    repository = recost_fixture["repository"]
    accounting = recost_fixture["accounting"]
    assert isinstance(root, Path)
    assert isinstance(repository, Path)
    assert isinstance(accounting, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    payload = rendered.stdout.encode()
    if packet_state == "exact":
        candidate.write_bytes(payload)
        candidate.chmod(0o644)
    request = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F119_recost_request.json"
    real_barrier = module.require_f119_install_namespace_barrier
    barriers = 0

    def create_request_at_post_barrier(*args, **kwargs):
        nonlocal barriers
        barriers += 1
        if barriers == 2:
            request.write_bytes(b"concurrent request commit marker\n")
            request.chmod(0o644)
        return real_barrier(*args, **kwargs)

    monkeypatch.setattr(
        module, "require_f119_install_namespace_barrier", create_request_at_post_barrier
    )
    with module.stage_i_lock(root) as mutation_lock:
        with pytest.raises(ValueError, match="request is the last commit marker"):
            module.retire_stale_f119_packet(
                root, repository, candidate, payload, mutation_lock
            )
    assert barriers == 2
    assert request.read_bytes() == b"concurrent request commit marker\n"
    assert candidate.exists() is (packet_state == "exact")


def test_f119_install_post_barrier_rejects_concurrent_request_marker(
    recost_fixture, monkeypatch
):
    module = load_recost_module()
    root = recost_fixture["root"]
    repository = recost_fixture["repository"]
    generator = recost_fixture["generator"]
    queue = recost_fixture["queue"]
    accounting = recost_fixture["accounting"]
    assert isinstance(root, Path)
    assert isinstance(repository, Path)
    assert isinstance(generator, Path)
    assert isinstance(queue, Path)
    assert isinstance(accounting, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    external = root.parent / "concurrent-marker-f119-candidate.json"
    external.write_text(rendered.stdout)
    external.chmod(0o644)
    request = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F119_recost_request.json"
    real_write = module.write_exact_or_verify
    raced = False

    def publish_then_create_request(*args, **kwargs):
        nonlocal raced
        created = real_write(*args, **kwargs)
        if not raced:
            request.write_bytes(b"concurrent request commit marker\n")
            request.chmod(0o644)
            raced = True
        return created

    monkeypatch.setattr(module, "write_exact_or_verify", publish_then_create_request)
    args = module.argparse.Namespace(
        action="install-f119-draft-packet",
        packet=external,
        expected_packet_sha256=sha256(external),
        squeue_file=queue,
    )
    with module.stage_i_lock(root) as mutation_lock:
        with pytest.raises(ValueError, match="request is the last commit marker"):
            module.locked_install_draft_packet(
                args,
                root,
                generator,
                repository,
                sha256(generator),
                mutation_lock,
            )
    assert raced
    assert candidate.read_bytes() == external.read_bytes()
    assert request.read_bytes() == b"concurrent request commit marker\n"


@pytest.mark.parametrize("existing", ("invalid", "request-committed"))
def test_f119_install_preserves_unsupersedable_canonical_state(
    recost_fixture, existing,
):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    seed, candidate = write_f119_seed_and_candidate(recost_fixture)
    candidate.unlink()
    rendered = run_action(
        recost_fixture,
        "render-f119-draft-packet",
        *f119_render_arguments(recost_fixture, seed),
    )
    assert rendered.returncode == 0, rendered.stderr
    if existing == "invalid":
        candidate.write_bytes(b"invalid canonical F119 packet\n")
    else:
        stale = json.loads(rendered.stdout)
        stale["generated_utc"] = (
            datetime.now(timezone.utc) - timedelta(minutes=4)
        ).isoformat()
        stale["expires_utc"] = (
            datetime.now(timezone.utc) - timedelta(seconds=1)
        ).isoformat()
        write_json(candidate, stale)
        request = candidate.with_name(
            "mks24_stage_i_E03_forcing_policy_F119_recost_request.json"
        )
        request.write_bytes(b"last commit marker\n")
    candidate.chmod(0o644)
    retained = candidate.read_bytes()
    external = root.parent / f"{existing}-replacement-f119.json"
    external.write_text(rendered.stdout)
    external.chmod(0o644)
    rejected = run_action(
        recost_fixture,
        "install-f119-draft-packet",
        "--packet",
        str(external),
        "--expected-packet-sha256",
        sha256(external),
    )
    assert rejected.returncode == 1
    assert candidate.read_bytes() == retained


@pytest.mark.parametrize("mutation", ("source-authority", "predecessor"))
def test_f119_draft_request_preflights_authority_before_publishing_prerequisites(
    recost_fixture, mutation
):
    def mutate(packet):
        if mutation == "source-authority":
            packet["inputs"]["source_authority"]["evidence"]["sha256"] = "0" * 64
        else:
            packet["inputs"]["predecessor_recost"]["sha256"] = "0" * 64

    packet = write_draft_packet(
        recost_fixture, checkpoint_number=119, mutate=mutate
    )
    rejected = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert rejected.returncode == 1
    accounting = recost_fixture["accounting"]
    assert isinstance(accounting, Path)
    prefix = f"mks24_stage_i_{EPOCH_SLUG}_F119"
    assert not any(
        path.exists()
        for path in (
            accounting / f"{prefix}_recost_request.json",
            accounting / f"{prefix}_reconciliation_evidence.json",
            accounting / f"{prefix}_storage_evidence.json",
        )
    )


def test_draft_request_validates_qualification_before_publishing_prerequisites(
    recost_fixture,
):
    packet = write_draft_packet(
        recost_fixture,
        checkpoint_number=206,
        mutate=lambda value: value["inputs"]["qualification_approval"].update(
            {"sha256": "0" * 64}
        ),
    )
    rejected = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert_rejected(rejected, "qualification approval checksum differs")
    accounting = recost_fixture["accounting"]
    assert isinstance(accounting, Path)
    prefix = f"mks24_stage_i_{EPOCH_SLUG}_F206"
    assert not any(
        path.exists()
        for path in (
            accounting / f"{prefix}_recost_request.json",
            accounting / f"{prefix}_reconciliation_evidence.json",
            accounting / f"{prefix}_storage_evidence.json",
        )
    )


def test_draft_request_preflights_all_targets_before_first_publication(
    recost_fixture,
):
    packet = write_draft_packet(recost_fixture, checkpoint_number=207)
    accounting = recost_fixture["accounting"]
    assert isinstance(accounting, Path)
    prefix = f"mks24_stage_i_{EPOCH_SLUG}_F207"
    request = accounting / f"{prefix}_recost_request.json"
    reconciliation = accounting / f"{prefix}_reconciliation_evidence.json"
    storage = accounting / f"{prefix}_storage_evidence.json"
    request.write_bytes(b"unrelated retained request namespace\n")
    request.chmod(0o644)
    retained = request.stat()

    rejected = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert_rejected(rejected, "exists with different bytes; refusing to clobber")
    assert request.read_bytes() == b"unrelated retained request namespace\n"
    assert (request.stat().st_dev, request.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )
    assert not reconciliation.exists()
    assert not storage.exists()


def test_draft_prerequisite_trio_rolls_back_failure_before_second_publish(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    publications = (
        (accounting / "F208_reconciliation.json", b"reconciliation\n", "reconciliation"),
        (accounting / "F208_storage.json", b"storage\n", "storage"),
        (accounting / "F208_request.json", b"request\n", "request"),
    )
    real_rename = module.renameat2_noreplace

    def fail_before_second_publish(directory, source, destination, label):
        if destination == publications[1][0].name:
            raise RuntimeError("persistent failure before second trio publish")
        return real_rename(directory, source, destination, label)

    monkeypatch.setattr(module, "renameat2_noreplace", fail_before_second_publish)
    with pytest.raises(
        ValueError, match="transaction-owned canonical links were rolled back"
    ):
        module.publish_draft_prerequisite_trio(publications)
    assert not any(path.exists() for path, _, _ in publications)
    retained_private = publication_private_entries(module, accounting)
    assert len(retained_private) == 3
    assert all(
        stat.S_IMODE(path.stat().st_mode) == 0o644 and path.stat().st_nlink == 1
        for path in retained_private
    )

    monkeypatch.setattr(module, "renameat2_noreplace", real_rename)
    module.publish_draft_prerequisite_trio(publications)
    for path, payload, _ in publications:
        assert path.read_bytes() == payload
        assert stat.S_IMODE(path.stat().st_mode) == 0o644
        assert path.stat().st_nlink == 1
    assert publication_private_entries(module, accounting) == []


def test_draft_prerequisite_trio_preserves_hostile_storage_race_and_recovers(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    publications = (
        (accounting / "F209_reconciliation.json", b"reconciliation\n", "reconciliation"),
        (accounting / "F209_storage.json", b"storage\n", "storage"),
        (accounting / "F209_request.json", b"request\n", "request"),
    )
    storage = publications[1][0]
    real_rename = module.renameat2_noreplace
    raced = False

    def race_storage_after_preflight(directory, source, destination, label):
        nonlocal raced
        if not raced and destination == storage.name:
            storage.write_bytes(b"hostile storage occupant\n")
            storage.chmod(0o644)
            raced = True
        return real_rename(directory, source, destination, label)

    monkeypatch.setattr(module, "renameat2_noreplace", race_storage_after_preflight)
    with pytest.raises(
        ValueError, match="transaction-owned canonical links were rolled back"
    ):
        module.publish_draft_prerequisite_trio(publications)
    assert raced
    assert storage.read_bytes() == b"hostile storage occupant\n"
    assert not publications[0][0].exists()
    assert not publications[2][0].exists()
    assert all(
        not path.exists() or path.read_bytes() != payload
        for path, payload, _ in publications
    )
    retained_private = publication_private_entries(module, accounting)
    assert len(retained_private) == 3
    assert all(
        stat.S_IMODE(path.stat().st_mode) == 0o644 and path.stat().st_nlink == 1
        for path in retained_private
    )

    monkeypatch.setattr(module, "renameat2_noreplace", real_rename)
    storage.unlink()
    module.publish_draft_prerequisite_trio(publications)
    for path, payload, _ in publications:
        assert path.read_bytes() == payload
        assert stat.S_IMODE(path.stat().st_mode) == 0o644
        assert path.stat().st_nlink == 1
    assert publication_private_entries(module, accounting) == []


def test_draft_prerequisite_trio_uses_moves_without_link_or_unlink(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    publications = (
        (accounting / "F210_reconciliation.json", b"reconciliation\n", "reconciliation"),
        (accounting / "F210_storage.json", b"storage\n", "storage"),
        (accounting / "F210_request.json", b"request\n", "request"),
    )
    def forbidden(*_args, **_kwargs):
        raise AssertionError("draft trio publication must not link or unlink")

    monkeypatch.setattr(module.os, "link", forbidden)
    monkeypatch.setattr(module.os, "unlink", forbidden)
    module.publish_draft_prerequisite_trio(publications)
    for path, payload, _ in publications:
        assert path.read_bytes() == payload
        assert stat.S_IMODE(path.stat().st_mode) == 0o644
        assert path.stat().st_nlink == 1
    assert publication_private_entries(module, accounting) == []


def test_bound_quarantine_retains_exact_inode_without_raw_deletion(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    source = parent / "request.json"
    payload = b"request\n"
    source.write_bytes(payload)
    source.chmod(0o644)
    expected = source.stat()

    def forbidden(*_args, **_kwargs):
        raise AssertionError("quarantine retention must not raw-unlink")

    monkeypatch.setattr(module.os, "unlink", forbidden)
    directory = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(ValueError, match="unsafe raw deletion refused"):
            module.unlink_bound_name_via_quarantine(
                directory,
                source.name,
                expected,
                payload,
                0o644,
                "request quarantine",
                expected_links=1,
            )
    finally:
        os.close(directory)
    quarantine = parent / f".{source.name}.delete-{expected.st_dev:x}-{expected.st_ino:x}"
    assert not source.exists()
    assert quarantine.read_bytes() == payload
    assert (quarantine.stat().st_dev, quarantine.stat().st_ino) == (
        expected.st_dev,
        expected.st_ino,
    )


def test_draft_prerequisite_trio_publishes_request_as_last_commit_marker(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    publications = (
        (accounting / "F210_reconciliation.json", b"reconciliation\n", "reconciliation"),
        (accounting / "F210_storage.json", b"storage\n", "storage"),
        (accounting / "F210_request.json", b"request\n", "request"),
    )
    real_rename = module.renameat2_noreplace
    committed = []

    def record_publication_order(directory, source, destination, label):
        result = real_rename(directory, source, destination, label)
        if destination in {path.name for path, _, _ in publications}:
            committed.append(destination)
        return result

    monkeypatch.setattr(module, "renameat2_noreplace", record_publication_order)
    module.publish_draft_prerequisite_trio(publications)
    assert committed == [path.name for path, _, _ in publications]
    assert committed[-1] == publications[-1][0].name
    for path, payload, _ in publications:
        assert path.read_bytes() == payload
        assert path.stat().st_nlink == 1
    assert publication_private_entries(module, accounting) == []


def test_draft_prerequisite_trio_reauthenticates_interrupted_final_durability(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    accounting = tmp_path / "accounting"
    accounting.mkdir()
    publications = (
        (accounting / "F211_reconciliation.json", b"reconciliation\n", "reconciliation"),
        (accounting / "F211_storage.json", b"storage\n", "storage"),
        (accounting / "F211_request.json", b"request\n", "request"),
    )
    real_authenticate = module.durably_authenticate_draft_trio_state
    request_single_attempts = 0

    def interrupt_first_request_single_authentication(*args, **kwargs):
        nonlocal request_single_attempts
        entry = args[2]
        expected_state = args[4]
        if entry.label == "request" and expected_state == "single":
            request_single_attempts += 1
            if request_single_attempts == 1:
                raise RuntimeError("interrupted before final request durability barrier")
        return real_authenticate(*args, **kwargs)

    monkeypatch.setattr(
        module,
        "durably_authenticate_draft_trio_state",
        interrupt_first_request_single_authentication,
    )
    module.publish_draft_prerequisite_trio(publications)

    assert request_single_attempts >= 2
    for path, payload, _ in publications:
        assert path.read_bytes() == payload
        assert stat.S_IMODE(path.stat().st_mode) == 0o644
        assert path.stat().st_nlink == 1
    assert publication_private_entries(module, accounting) == []


def linked_draft_trio_entry(module, parent: Path):
    """Return one exact linked draft-trio entry for direct lifecycle probes."""

    target = parent / "request.json"
    payload = b"exact transaction-owned request\n"
    transaction = module.publication_transaction(target, payload, 0o644, "request")
    private = parent / transaction.private_name(0)
    private.write_bytes(payload)
    private.chmod(0o644)
    os.link(private, target)
    identity = target.stat()
    entry = module.DraftTrioPublicationEntry(
        path=target,
        payload=payload,
        mode=0o644,
        label="request",
        transaction=transaction,
        initial_state="linked",
        private_name=private.name,
        private_identity=identity,
        installed_identity=identity,
    )
    return entry, private


def test_draft_trio_descriptor_bound_rollback_preserves_hostile_swap(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    entry, private = linked_draft_trio_entry(module, parent)
    hostile = parent / "hostile"
    hostile.write_bytes(b"hostile namespace replacement\n")
    hostile.chmod(0o644)
    displaced = parent / "displaced-exact-transaction-owned-link"
    real_rename = module.renameat2_noreplace
    raced = False

    def swap_before_bound_move(directory, source, target, label):
        nonlocal raced
        selected = entry.path.name
        if not raced and source == selected:
            os.rename(source, displaced.name, src_dir_fd=directory, dst_dir_fd=directory)
            os.rename(hostile.name, source, src_dir_fd=directory, dst_dir_fd=directory)
            raced = True
        return real_rename(directory, source, target, label)

    monkeypatch.setattr(module, "renameat2_noreplace", swap_before_bound_move)
    directory = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        parent_profile = os.fstat(directory)
        with pytest.raises(ValueError):
            module.rollback_draft_trio_entries(
                directory, parent_profile, [entry], None
            )
    finally:
        os.close(directory)
    assert raced
    assert entry.path.read_bytes() == b"hostile namespace replacement\n"
    assert displaced.read_bytes() == entry.payload
    assert private.read_bytes() == entry.payload


def test_quarantine_bound_move_preserves_hostile_swap(tmp_path, monkeypatch):
    module = load_recost_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    source = parent / "request.json"
    payload = b"exact transaction-owned request\n"
    source.write_bytes(payload)
    source.chmod(0o644)
    expected = source.stat()
    hostile = parent / "hostile"
    hostile.write_bytes(b"hostile namespace replacement\n")
    hostile.chmod(0o644)
    displaced = parent / "displaced-exact-transaction-owned-request"
    real_rename = module.renameat2_noreplace
    raced = False

    def swap_before_bound_move(directory, selected, target, label):
        nonlocal raced
        if not raced and selected == source.name:
            os.rename(selected, displaced.name, src_dir_fd=directory, dst_dir_fd=directory)
            os.rename(hostile.name, selected, src_dir_fd=directory, dst_dir_fd=directory)
            raced = True
        return real_rename(directory, selected, target, label)

    monkeypatch.setattr(module, "renameat2_noreplace", swap_before_bound_move)
    directory = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(ValueError, match="hostile replacement"):
            module.unlink_bound_name_via_quarantine(
                directory,
                source.name,
                expected,
                payload,
                0o644,
                "request quarantine",
                expected_links=1,
            )
    finally:
        os.close(directory)
    assert raced
    assert source.read_bytes() == b"hostile namespace replacement\n"
    assert displaced.read_bytes() == b"exact transaction-owned request\n"


def test_draft_trio_rejects_ambiguous_exact_private_attempts(tmp_path):
    module = load_recost_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    target = parent / "request.json"
    payload = b"request\n"
    transaction = module.publication_transaction(target, payload, 0o644, "request")
    entry = module.DraftTrioPublicationEntry(
        path=target,
        payload=payload,
        mode=0o644,
        label="request",
        transaction=transaction,
        initial_state="absent",
    )
    for attempt in (0, 1):
        private = parent / transaction.private_name(attempt)
        private.write_bytes(payload)
        private.chmod(0o644)
    directory = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(ValueError, match="ambiguous retained"):
            module.reusable_draft_trio_private(directory, entry)
    finally:
        os.close(directory)


def test_draft_trio_exact_final_rejects_transaction_private_remnants(tmp_path):
    module = load_recost_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    target = parent / "request.json"
    payload = b"request\n"
    target.write_bytes(payload)
    target.chmod(0o644)
    transaction = module.publication_transaction(target, payload, 0o644, "request")
    private = parent / transaction.private_name(0)
    private.write_bytes(payload)
    private.chmod(0o644)
    entry = module.DraftTrioPublicationEntry(
        path=target,
        payload=payload,
        mode=0o644,
        label="request",
        transaction=transaction,
        initial_state="single",
    )
    directory = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(ValueError, match="ambiguous transaction-private remnants"):
            module.classify_draft_trio_target(directory, entry)
    finally:
        os.close(directory)
    assert target.read_bytes() == payload
    assert private.read_bytes() == payload


def test_draft_request_requires_exact_checkpoint_packet_namespace(recost_fixture):
    packet = write_draft_packet(recost_fixture, checkpoint_number=205)
    displaced = packet.with_name("reviewed-draft-packet.json")
    packet.replace(displaced)
    assert_rejected(
        run_action(
            recost_fixture,
            "draft-request",
            "--packet",
            str(displaced),
            "--expected-packet-sha256",
            sha256(displaced),
        ),
        "path differs from its checkpoint namespace",
    )
    assert not (
        recost_fixture["accounting"]
        / f"mks24_stage_i_{EPOCH_SLUG}_F205_recost_request.json"
    ).exists()


def test_f117_managed_workflow_reaches_staged_generation_without_manual_writes(
    recost_fixture,
):
    accounting = recost_fixture["accounting"]
    root = recost_fixture["root"]
    predecessor = recost_fixture["predecessor"]
    predecessor_review = recost_fixture["predecessor_review"]
    predecessor_audit = recost_fixture["predecessor_audit"]
    request_template = recost_fixture["request"]
    assert all(
        isinstance(path, Path)
        for path in (
            accounting,
            root,
            predecessor,
            predecessor_review,
            predecessor_audit,
            request_template,
        )
    )
    request_timestamp = datetime.fromisoformat(
        json.loads(request_template.read_text())["generated_utc"]
    )
    for path in (predecessor_review, predecessor_audit, predecessor):
        path.unlink()

    schema2 = accounting / f"mks24_stage_i_{EPOCH_SLUG}_F115_recost_evidence.json"
    write_immutable_json(
        schema2,
        {
            "schema_version": 2,
            "record_type": "stage-i-recost-recommendation-evidence",
            "checkpoint": "F-115",
            "artifact_name": schema2.name,
            "execution_epoch": EPOCH,
            "generated_utc": (request_timestamp - timedelta(minutes=10)).isoformat(),
            "authority": {
                "authorizing": False,
                "action_authority": "none-until-independent-review-and-publication",
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            },
        },
    )
    schema2_review = schema2.with_name(f"{schema2.name}.independent_review.json")
    write_immutable_json(
        schema2_review,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (request_timestamp - timedelta(minutes=9)).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": {
                "agent_id": "fixture-f115-recost-reviewer",
                "independent_from_generator": True,
            },
            "candidate": {"path": str(schema2), "sha256": sha256(schema2)},
            "scope": {"non_authorizing": True},
        },
    )
    schema2_audit = schema2.with_name(f"{schema2.name}.publication_audit.json")
    write_immutable_json(
        schema2_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": (request_timestamp - timedelta(minutes=8)).isoformat(),
            "artifact": {
                "path": str(schema2),
                "sha256": sha256(schema2),
                "mode": "0444",
                "links": 1,
            },
            "independent_review": {
                "path": str(schema2_review),
                "sha256": sha256(schema2_review),
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

    def bind_schema2_f115(packet):
        packet["inputs"]["predecessor_recost"] = {
            "path": schema2.relative_to(root).as_posix(),
            "sha256": sha256(schema2),
        }
        packet["inputs"]["predecessor_recost_independent_review"] = {
            "path": schema2_review.relative_to(root).as_posix(),
            "sha256": sha256(schema2_review),
        }
        packet["inputs"]["predecessor_recost_publication_audit"] = {
            "path": schema2_audit.relative_to(root).as_posix(),
            "sha256": sha256(schema2_audit),
        }

    packet = write_draft_packet(
        recost_fixture, checkpoint_number=117, mutate=bind_schema2_f115
    )
    external_packet = root.parent / "externally-reviewed-f117-draft-packet.json"
    packet.replace(external_packet)
    installed_packet = run_action(
        recost_fixture,
        "install-f117-draft-packet",
        "--packet",
        str(external_packet),
        "--expected-packet-sha256",
        sha256(external_packet),
    )
    assert installed_packet.returncode == 0, installed_packet.stderr
    packet = Path(action_report(installed_packet)["path"])
    drafted = run_action(
        recost_fixture,
        "draft-request",
        "--packet",
        str(packet),
        "--expected-packet-sha256",
        sha256(packet),
    )
    assert drafted.returncode == 0, drafted.stderr

    prefix = f"mks24_stage_i_{EPOCH_SLUG}_F117"
    request = accounting / f"{prefix}_recost_request.json"
    external_review = write_external_request_review(recost_fixture, request)
    installed_review = run_action(
        recost_fixture,
        "install-request-review",
        "--request",
        str(request),
        "--expected-request-sha256",
        sha256(request),
        "--review",
        str(external_review),
        "--expected-review-sha256",
        sha256(external_review),
    )
    assert installed_review.returncode == 0, installed_review.stderr
    workflow_fixture = dict(recost_fixture)
    workflow_fixture["request"] = request
    workflow_fixture["request_sha256"] = sha256(request)
    workflow_fixture["output"] = accounting / f"{prefix}_recost_evidence.json.staged"
    generated = run_generator(workflow_fixture)
    assert generated.returncode == 0, generated.stderr
    assert stat.S_IMODE(workflow_fixture["output"].stat().st_mode) == 0o444
    artifact = json.loads(workflow_fixture["output"].read_text())
    assert artifact["checkpoint"] == "F-117"
    assert artifact["predecessor_recost"]["checkpoint"] == "F-115"


def test_retain_generator_action_is_no_clobber_or_exact_verify(recost_fixture):
    accounting = recost_fixture["accounting"]
    generator = recost_fixture["generator"]
    assert isinstance(accounting, Path)
    assert isinstance(generator, Path)
    utilities = accounting / "utilities"
    utilities.mkdir()
    legacy_f117 = utilities / "cgl_lf_stage_i_recost.py"
    legacy_f117.write_bytes(b"immutable retained F117 generator\n")
    legacy_f117.chmod(0o755)
    legacy_snapshot = (
        legacy_f117.read_bytes(),
        legacy_f117.stat().st_dev,
        legacy_f117.stat().st_ino,
        stat.S_IMODE(legacy_f117.stat().st_mode),
    )

    created = run_action(recost_fixture, "retain-generator")
    assert created.returncode == 0, created.stderr
    created_report = action_report(created)
    assert created_report["created"] is True
    retained = Path(created_report["path"])
    assert retained.name == f"cgl_lf_stage_i_recost.sha256-{sha256(generator)}.py"
    assert created_report["checkpoint_generator_relative_path"] == (
        retained.relative_to(recost_fixture["root"]).as_posix()
    )
    assert retained.read_bytes() == generator.read_bytes()
    assert stat.S_IMODE(retained.stat().st_mode) == 0o755
    assert (
        legacy_f117.read_bytes(),
        legacy_f117.stat().st_dev,
        legacy_f117.stat().st_ino,
        stat.S_IMODE(legacy_f117.stat().st_mode),
    ) == legacy_snapshot

    verified = run_action(recost_fixture, "retain-generator")
    assert verified.returncode == 0, verified.stderr
    verified_report = action_report(verified)
    assert verified_report["created"] is False
    assert verified_report["exact_existing_copy_verified"] is True

    retained.write_bytes(retained.read_bytes() + b"\n# retained drift\n")
    retained.chmod(0o755)
    drift = retained.read_bytes()
    assert_rejected(
        run_action(recost_fixture, "retain-generator"),
        "exists with different bytes; refusing to clobber",
    )
    assert retained.read_bytes() == drift


def test_retain_generator_keeps_committed_versions_as_content_addressed_siblings(
    recost_fixture,
):
    repository = recost_fixture["repository"]
    generator = recost_fixture["generator"]
    assert isinstance(repository, Path)
    assert isinstance(generator, Path)

    first = run_action(recost_fixture, "retain-generator")
    assert first.returncode == 0, first.stderr
    first_path = Path(action_report(first)["path"])
    first_bytes = first_path.read_bytes()

    generator.write_bytes(generator.read_bytes() + b"\n# committed F119 generator version\n")
    generator.chmod(0o755)
    git(repository, "add", str(generator.relative_to(repository)))
    git(
        repository,
        "-c",
        "user.name=CGL fixture",
        "-c",
        "user.email=cgl-fixture@example.invalid",
        "commit",
        "-q",
        "-m",
        "Commit next recost generator version",
    )
    second = run_action(recost_fixture, "retain-generator")
    assert second.returncode == 0, second.stderr
    second_path = Path(action_report(second)["path"])

    assert second_path != first_path
    assert first_path.read_bytes() == first_bytes
    assert second_path.read_bytes() == generator.read_bytes()
    assert first_path.name == f"cgl_lf_stage_i_recost.sha256-{sha256(first_path)}.py"
    assert second_path.name == f"cgl_lf_stage_i_recost.sha256-{sha256(generator)}.py"
    assert not (first_path.parent / "cgl_lf_stage_i_recost.py").exists()


def test_prerequisite_actions_require_drained_queue_and_empty_transactions(
    recost_fixture,
):
    packet = write_draft_packet(recost_fixture, checkpoint_number=204)
    initial_f117 = write_draft_packet(recost_fixture, checkpoint_number=117)
    initial_f119 = write_draft_packet(recost_fixture, checkpoint_number=119)
    queue = recost_fixture["queue"]
    transaction_store = recost_fixture["transaction_store"]
    request = recost_fixture["request"]
    assert isinstance(queue, Path)
    assert isinstance(transaction_store, Path)
    assert isinstance(request, Path)
    external_review = write_external_request_review(recost_fixture, request)
    review_arguments = (
        "--request",
        str(request),
        "--expected-request-sha256",
        sha256(request),
        "--review",
        str(external_review),
        "--expected-review-sha256",
        sha256(external_review),
    )
    queue.write_text("999|cgl_stage_i_writer|RUNNING\n")
    for action, arguments in (
        (
            "install-f117-draft-packet",
            (
                "--packet",
                str(initial_f117),
                "--expected-packet-sha256",
                sha256(initial_f117),
            ),
        ),
        (
            "install-f119-draft-packet",
            (
                "--packet",
                str(initial_f119),
                "--expected-packet-sha256",
                sha256(initial_f119),
            ),
        ),
        (
            "draft-request",
            (
                "--packet",
                str(packet),
                "--expected-packet-sha256",
                sha256(packet),
            ),
        ),
        ("install-request-review", review_arguments),
        ("retain-generator", ()),
    ):
        assert_rejected(
            run_action(recost_fixture, action, *arguments),
            "active CGL scheduler queue is not drained",
        )
    queue.write_text("")
    write_json(transaction_store / "pending.json", {"state": "pending"})
    for action, arguments in (
        (
            "install-f117-draft-packet",
            (
                "--packet",
                str(initial_f117),
                "--expected-packet-sha256",
                sha256(initial_f117),
            ),
        ),
        (
            "install-f119-draft-packet",
            (
                "--packet",
                str(initial_f119),
                "--expected-packet-sha256",
                sha256(initial_f119),
            ),
        ),
        (
            "draft-request",
            (
                "--packet",
                str(packet),
                "--expected-packet-sha256",
                sha256(packet),
            ),
        ),
        ("install-request-review", review_arguments),
        ("retain-generator", ()),
    ):
        assert_rejected(
            run_action(recost_fixture, action, *arguments),
            "transaction store is not empty",
        )


def test_generator_emits_bounded_wave_as_explicitly_non_authorizing(recost_fixture):
    add_recorded_wave_job(recost_fixture)
    profiles = [
        profile(recost_fixture, case_id="R04", nodes=4, parent=False),
        profile(recost_fixture, case_id="R16", nodes=2, parent=False),
    ]
    set_profiles(recost_fixture, profiles, mode="bounded-wave")
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    recommendations = artifact["recommendations"]
    assert recommendations["mode"] == "bounded-wave"
    assert recommendations["authorizing"] is False
    assert "sole_next_segment_recommendation" not in recommendations
    assert "evidence only" in recommendations["non_authorizing_reason"]
    assert recommendations["bounded_concurrency"]["max_wave_nodes"] == 6
    assert len(artifact["barrier"]["scheduler_evidence"]) == 2
    assert "scheduler_sha256" not in artifact["provenance"]


def test_generator_requires_exact_independent_request_review(recost_fixture):
    review = recost_fixture["request_review"]
    assert isinstance(review, Path)
    review.unlink()
    assert_rejected(run_generator(recost_fixture), "independent_review.json")
    refresh_request(recost_fixture)
    value = json.loads(review.read_text())
    value["reviewer"]["agent_id"] = "fixture-request-author"
    write_immutable_json(review, value)
    assert_rejected(
        run_generator(recost_fixture),
        "lacks the exact declared non-cryptographic process-independence assurance",
    )


def test_request_review_scope_is_exactly_non_authorizing(recost_fixture):
    review = recost_fixture["request_review"]
    assert isinstance(review, Path)
    value = json.loads(review.read_text())
    assert value["scope"] == {"non_authorizing": True}
    value["scope"]["scheduler_mutation_authorized"] = False
    write_immutable_json(review, value)
    assert_rejected(run_generator(recost_fixture), "over-authorizes")


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("declared-independent-false", "lacks the exact declared non-cryptographic"),
        ("cryptographic-identity-claim", "lacks the exact declared non-cryptographic"),
        ("stronger-proof-field", "recost request reviewer schema differs"),
    ),
)
def test_request_review_identity_is_declared_noncryptographic_only(
    recost_fixture, mutation, message,
):
    review = recost_fixture["request_review"]
    assert isinstance(review, Path)
    value = json.loads(review.read_text())
    if mutation == "declared-independent-false":
        value["reviewer"]["declared_process_independence"] = False
    elif mutation == "cryptographic-identity-claim":
        value["reviewer"]["identity_assurance"] = "cryptographically-verified-independent"
    else:
        value["reviewer"]["identity_verified"] = True
    write_immutable_json(review, value)
    assert_rejected(run_generator(recost_fixture), message)


def test_generator_requires_f118_current_source_and_build_qualification_chains(
    recost_fixture,
):
    authority = recost_fixture["source_authority"]
    assert isinstance(authority, Path)
    authority.chmod(0o644)
    assert_rejected(run_generator(recost_fixture), "F118 evidence mode is 0644, expected 0444")

    authority.chmod(0o444)
    qualification = recost_fixture["qualification"]
    assert isinstance(qualification, Path)
    value = json.loads(qualification.read_text())
    value["approved_executable_sha256"] = "0" * 64
    write_json(qualification, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "qualification approval build binding differs")


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("scope-preserves", "F118 scope differs or broadens authority"),
        ("scope-does-not-authorize", "F118 scope differs or broadens authority"),
        ("authorization", "F118 current source authority identity differs"),
        ("validation", "F118 current source authority identity differs"),
        ("publication-requirements", "F118 current source authority identity differs"),
        ("publisher", "F118 publisher binding differs from committed tools"),
        ("current-advertised-tip", "advertised tip is not exact HEAD"),
        ("source-archive-catalog", "source-archive catalog semantics differ"),
        ("review-decision", "F118 provenance_review identity differs"),
        ("review-candidate", "do not bind one exact candidate"),
        ("review-verified", "F118 provenance_review identity differs"),
        ("review-limitation", "lacks the non-cryptographic reviewer identity limitation"),
        ("review-future", "F118 provenance_review review chronology differs"),
        ("audit-publication", "F118 publication audit identity or authority differs"),
        ("audit-bridge", "F118 source-archive catalog authority differs"),
        ("audit-future", "F118 publication audit identity or authority differs"),
    ),
)
def test_generator_requires_exact_f118_contract_review_and_audit_semantics(
    recost_fixture, mutation, message,
):
    authority = recost_fixture["source_authority"]
    provenance = recost_fixture["source_authority_provenance_review"]
    audit = recost_fixture["source_authority_audit"]
    assert all(isinstance(path, Path) for path in (authority, provenance, audit))
    if mutation.startswith("scope-") or mutation in {
        "authorization",
        "validation",
        "publication-requirements",
        "publisher",
        "current-advertised-tip",
        "source-archive-catalog",
    }:
        value = json.loads(authority.read_text())
        if mutation == "scope-preserves":
            value["scope"]["preserves"][0] = "Mutated preservation claim."
        elif mutation == "scope-does-not-authorize":
            value["scope"]["does_not_authorize"].pop()
        elif mutation == "authorization":
            value["authorization"]["prepare_authorized"] = True
        elif mutation == "validation":
            value["validation"]["historical_f116_chain"] = "failed"
        elif mutation == "publication-requirements":
            value["publication_requirements"]["published_audit_mode"] = "0644"
        elif mutation == "publisher":
            value["implementation"]["publisher"]["sha256"] = "0" * 64
        elif mutation == "current-advertised-tip":
            value["implementation"]["current_source_bundle"]["advertised_tip"][
                "name"
            ] = "refs/heads/feature/cgl-landau-fluid"
        else:
            value["source_archive_catalog"]["after"][
                "all_prior_checksum_entries_preserved"
            ] = False
        write_immutable_json(authority, value)
    elif mutation.startswith("review-"):
        value = json.loads(provenance.read_text())
        if mutation == "review-decision":
            value["decision"] = "approved"
        elif mutation == "review-candidate":
            value["reviewed_candidate"]["path"] += ".different"
        elif mutation == "review-verified":
            value["verified"]["current_source_selection_only"] = False
        elif mutation == "review-limitation":
            value["limitations"].remove(INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION)
        else:
            value["reviewed_utc"] = (
                datetime.now(timezone.utc) + timedelta(minutes=10)
            ).isoformat()
        write_immutable_json(provenance, value)
    else:
        value = json.loads(audit.read_text())
        if mutation == "audit-publication":
            value["publication"] = "mutated-publication-method"
        elif mutation == "audit-bridge":
            value["source_archive_catalog"]["bridge_bundle"]["sha256"] = "0" * 64
        else:
            value["published_utc"] = (
                datetime.now(timezone.utc) + timedelta(minutes=10)
            ).isoformat()
        write_immutable_json(audit, value)
    refresh_f118_historical_f116_binding(recost_fixture)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), message)


def test_generator_rejects_legacy_expanded_source_authority_binding(
    recost_fixture,
):
    request = recost_fixture["request"]
    review = recost_fixture["request_review"]
    assert isinstance(request, Path)
    assert isinstance(review, Path)
    value = json.loads(request.read_text())
    authority = value["inputs"]["source_authority"]
    value["inputs"]["source_authority"] = {
        "checkpoint": "F-116",
        "evidence_sha256": authority["evidence"]["sha256"],
        "publication_audit_sha256": authority["publication_audit"]["sha256"],
        "historical_f115": {"legacy_digest_summary": True},
        "current_source_bundle": authority["final_source_bundle"],
    }
    write_json(request, value)
    recost_fixture["request_sha256"] = sha256(request)
    review_value = json.loads(review.read_text())
    review_value["candidate"]["sha256"] = recost_fixture["request_sha256"]
    write_immutable_json(review, review_value)
    assert_rejected(
        run_generator(recost_fixture),
        "F118 current source authority bindings schema differs",
    )


def test_generator_rejects_legacy_source_authority_publisher_schema(recost_fixture):
    authority = recost_fixture["source_authority"]
    assert isinstance(authority, Path)
    value = json.loads(authority.read_text())
    value["predecessors"] = value.pop("predecessor_authorities")
    value["implementation"] = {
        "source_bundle": value["implementation"]["current_source_bundle"]
    }
    write_immutable_json(authority, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "F118 evidence schema differs")


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("missing", "must contain exactly seven tools"),
        ("duplicate", "identity, order, revision, or mode differs"),
        ("mode", "identity, order, revision, or mode differs"),
        ("digest", "historical revision bytes differ"),
    ),
)
def test_f118_consumer_requires_exact_seven_committed_tools(
    recost_fixture, mutation, message,
):
    module = load_recost_module()
    repository = recost_fixture["repository"]
    revision = recost_fixture["revision"]
    assert isinstance(repository, Path)
    assert isinstance(revision, str)
    tools = f116_committed_tools(repository, revision)
    assert module.authenticate_f118_committed_tools(tools, repository, revision) == tools
    invalid = deepcopy(tools)
    if mutation == "missing":
        invalid.pop()
    elif mutation == "duplicate":
        invalid[1]["path"] = invalid[0]["path"]
    elif mutation == "mode":
        invalid[0]["mode"] = "0755" if invalid[0]["mode"] == "0644" else "0644"
    else:
        invalid[0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match=message):
        module.authenticate_f118_committed_tools(invalid, repository, revision)


def test_generator_rejects_f118_six_part_binding_to_historical_f115_bundle(
    recost_fixture,
):
    request = recost_fixture["request"]
    review = recost_fixture["request_review"]
    historical = recost_fixture["historical_source_bundle"]
    assert isinstance(request, Path)
    assert isinstance(review, Path)
    assert isinstance(historical, Path)
    value = json.loads(request.read_text())
    value["inputs"]["source_authority"]["final_source_bundle"] = {
        "path": historical.relative_to(recost_fixture["root"]).as_posix(),
        "sha256": sha256(historical),
        "verified_revisions": [recost_fixture["f113_revision"]],
    }
    write_json(request, value)
    recost_fixture["request_sha256"] = sha256(request)
    review_value = json.loads(review.read_text())
    review_value["candidate"]["sha256"] = recost_fixture["request_sha256"]
    write_immutable_json(review, review_value)
    assert_rejected(
        run_generator(recost_fixture),
        "F118 six-part binding does not select the current source bundle",
    )


@pytest.mark.parametrize(
    ("terminal_result", "scheduler_state", "exit_code"),
    (
        ("failed", "FAILED", "1:0"),
        ("aborted", "CANCELLED", "0:15"),
        ("rejected", "COMPLETED", "0:0"),
    ),
)
def test_generator_accounts_cancelled_and_terminal_without_fabricating_science(
    recost_fixture,
    terminal_result,
    scheduler_state,
    exit_code,
):
    add_cancelled_identity(recost_fixture)
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    assert artifact["ledger"]["rows"] == 1
    assert artifact["reservations"]["rows"] == 2
    recost_fixture["output"].unlink()
    retry = profile(
        recost_fixture,
        case_id="R04",
        parent=False,
        segment="s01_rankio_t0_t0p25",
    )
    set_profiles(recost_fixture, [retry], mode="sole-next-profile")
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    recost_fixture["output"].unlink()

    replace_barrier_with_failed_outcome(
        recost_fixture,
        result=terminal_result,
        scheduler_state=scheduler_state,
        exit_code=exit_code,
    )
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    assert artifact["barrier"]["recorded_segments"][0]["result"] == terminal_result
    assert artifact["barrier"]["scheduler_evidence"][0]["state"] == scheduler_state
    assert artifact["ledger"]["rows"] == 2


@pytest.mark.parametrize(
    ("case_id", "nodes"),
    tuple((f"R{number:02d}", 4) for number in range(4, 16)) + (("R16", 2),),
)
def test_generator_accepts_qualification_backed_initial_profiles_r04_through_r16(
    recost_fixture,
    case_id,
    nodes,
):
    set_profiles(
        recost_fixture,
        [profile(recost_fixture, case_id=case_id, nodes=nodes, parent=False)],
        mode="sole-next-profile",
    )
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    retained = artifact["recommendations"]["recommended_next_profiles"][0]
    assert retained["case_id"] == case_id
    assert retained["nodes"] == nodes
    assert retained["time_tlim_target"] == INITIAL_TARGETS[case_id]
    assert retained["recommendation_basis"]["kind"] == "qualified-initial-calibration"


def test_generator_recommendations_are_measured_and_ct_claim_is_truthful(recost_fixture):
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    basis = artifact["recommendations"]["recommended_next_profiles"][0][
        "recommendation_basis"
    ]
    assert basis["kind"] == "measured-production-continuation"
    serialized = json.dumps(artifact)
    assert "maximum normalized CT divergence" not in serialized
    assert "makes no CT-divergence bound claim" in serialized

    recost_fixture["output"].unlink()

    def mutate(request):
        request["recommendations"]["profiles"][0]["recommendation_basis"][
            "maximum_recommended_interval"
        ] = "9"

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "recommendation basis differs")


def test_generator_rejects_cross_name_replay_and_checkpoint_mismatch(recost_fixture):
    other = recost_fixture["accounting"] / "mks24_stage_i_E03_forcing_policy_F201_recost_evidence.json.staged"
    assert isinstance(other, Path)
    assert_rejected(run_generator(recost_fixture, output=other), "artifact name differs")

    def mutate(request):
        request["checkpoint"] = "F-201"

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "artifact name and checkpoint ID differ")


def test_generator_rejects_active_queue_reservation_and_transactions(recost_fixture):
    queue = recost_fixture["queue"]
    assert isinstance(queue, Path)
    queue.write_text("999|unrelated-job|RUNNING\n")
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    recost_fixture["output"].unlink()
    queue.write_text("999|cgl_other_stage_i_writer|RUNNING\n")
    assert_rejected(run_generator(recost_fixture), "active CGL scheduler queue is not drained")
    queue.write_text("")
    reservations = recost_fixture["reservations"]
    assert isinstance(reservations, Path)
    value = json.loads(reservations.read_text())
    value[0]["state"] = "submitted"
    write_json(reservations, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "remains active")


def test_generator_rejects_reservation_numeric_type_replay(recost_fixture):
    reservations = recost_fixture["reservations"]
    assert isinstance(reservations, Path)
    value = json.loads(reservations.read_text())
    value[0]["nodes"] = True
    write_json(reservations, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "reservation 12345 nodes must be an integer")


@pytest.mark.parametrize("store_key", ["transaction_store", "recost_transaction_store"])
def test_generator_rejects_transaction_journal(recost_fixture, store_key):
    store = recost_fixture[store_key]
    assert isinstance(store, Path)
    write_json(store / "pending.json", {"state": "pending"})
    assert_rejected(run_generator(recost_fixture), "transaction store is not empty")


def test_generator_rejects_stale_hash_and_self_asserted_projection(recost_fixture):
    ledger = recost_fixture["ledger"]
    assert isinstance(ledger, Path)
    ledger.write_text(ledger.read_text() + "\n")
    assert_rejected(run_generator(recost_fixture), "Stage I ledger checksum differs")

    fixture = recost_fixture
    refresh_request(fixture)

    def mutate(request):
        request["budget_projection"] = {"projected_stage_i_total_node_hours": "1"}

    refresh_request(fixture, mutate)
    assert_rejected(run_generator(fixture), "recost request schema differs")


def test_generator_requires_exact_f113_path_digest_and_publication_audit(recost_fixture):
    root = recost_fixture["root"]
    ceiling = recost_fixture["ceiling"]
    assert isinstance(root, Path)
    assert isinstance(ceiling, Path)
    lookalike = root / "accounting/lookalike_F113.json"
    lookalike.write_bytes(ceiling.read_bytes())
    lookalike.chmod(0o644)

    def mutate(request):
        request["inputs"]["ceiling_evidence"] = {
            "path": str(lookalike.relative_to(root)),
            "sha256": sha256(lookalike),
        }

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "not the exact promoted F113 artifact")


def test_generator_rejects_f113_audit_replay_and_fabricated_lane(recost_fixture):
    ceiling = recost_fixture["ceiling"]
    assert isinstance(ceiling, Path)
    value = json.loads(ceiling.read_text())
    value["promoted_controls"]["campaign_budget_node_hours"] = 1300.0
    write_json(ceiling, value)
    refresh_request(recost_fixture)
    assert_rejected(
        run_generator(recost_fixture),
        "F115 source authority F113 predecessor binding differs",
    )

    value["promoted_controls"]["standard_active_lane_limit"] = 5
    write_json(ceiling, value)
    refresh_ceiling(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "lane limit must remain exactly four")


def test_generator_rejects_self_asserted_f113_reviewer_authority(recost_fixture):
    audit = recost_fixture["ceiling_audit"]
    assert isinstance(audit, Path)
    value = json.loads(audit.read_text())
    value["review"]["reviewed_by"] = "self-asserted reviewer"
    write_json(audit, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "lacks committed campaign authority")


def test_generator_accepts_exact_f116_supersession_of_missing_f113_audit(recost_fixture):
    configure_f113_audit_supersession(
        recost_fixture, retain_historical_audit=False
    )
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    assert artifact["promoted_f113"]["publication_audit"]["supersession"]["status"] == (
        "controlled-exact-transitive-F115-F116-F118-supersession"
    )


def test_generator_accepts_exact_f118_audit_supersession_while_f113_audit_is_retained(
    recost_fixture,
):
    historical = recost_fixture["ceiling_audit"]
    assert isinstance(historical, Path)
    snapshot = historical.read_bytes()
    configure_f113_audit_supersession(
        recost_fixture, retain_historical_audit=True
    )
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr
    artifact = json.loads(recost_fixture["output"].read_text())
    supersession = artifact["promoted_f113"]["publication_audit"]["supersession"]
    assert supersession["status"] == (
        "exact-retained-F118-current-source-authority-supersession"
    )
    assert supersession["historical_publication_audit_absent"] is False
    assert historical.read_bytes() == snapshot


def test_generator_rejects_f116_supersession_without_exact_historical_f115(recost_fixture):
    audit = recost_fixture["f116_source_authority_audit"]
    assert isinstance(audit, Path)
    audit_value = json.loads(audit.read_text())
    audit_value["historical_f115_authority"]["evidence_sha256"] = "0" * 64
    write_immutable_json(audit, audit_value)
    refresh_f118_historical_f116_binding(recost_fixture)
    refresh_request(recost_fixture)
    assert_rejected(
        run_generator(recost_fixture),
        "historical F116 publication audit identity differs",
    )


def test_generator_rejects_f113_identity_and_publication_chronology(recost_fixture):
    ceiling = recost_fixture["ceiling"]
    timestamp = recost_fixture["timestamp"]
    assert isinstance(ceiling, Path)
    assert isinstance(timestamp, datetime)
    value = json.loads(ceiling.read_text())
    value["checkpoint"] = "F-112"
    write_json(ceiling, value)
    refresh_ceiling(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "checkpoint differs from F-113")

    value["checkpoint"] = "F-113"
    value["generated_utc"] = timestamp.isoformat()
    write_json(ceiling, value)
    refresh_ceiling(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "publication predates")


def test_generator_authenticates_historical_f113_helper_bytes(recost_fixture):
    ceiling = recost_fixture["ceiling"]
    assert isinstance(ceiling, Path)
    value = json.loads(ceiling.read_text())
    value["implementation"]["stage_i_helper"]["sha256"] = "0" * 64
    write_json(ceiling, value)
    refresh_ceiling(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "historical revision bytes differ")


def test_generator_rejects_scheduler_chronology_and_filename_replay(recost_fixture):
    manifest = recost_fixture["manifest"]
    ledger = recost_fixture["ledger"]
    assert isinstance(manifest, Path)
    assert isinstance(ledger, Path)
    value = json.loads(manifest.read_text())
    completed = value["accounting"]["completed_utc"]
    value["accounting"]["completed_utc"] = value["accounting"]["submitted_utc"]
    write_json(manifest, value)
    ledger.write_text(ledger.read_text().replace(completed, value["accounting"]["submitted_utc"], 1))
    ledger.chmod(0o644)
    rewrite_scheduler_from_manifest(recost_fixture)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "chronology is invalid")


def test_generator_rejects_scheduler_completion_after_request(recost_fixture):
    manifest = recost_fixture["manifest"]
    ledger = recost_fixture["ledger"]
    timestamp = recost_fixture["timestamp"]
    assert isinstance(manifest, Path)
    assert isinstance(ledger, Path)
    assert isinstance(timestamp, datetime)
    value = json.loads(manifest.read_text())
    original = value["accounting"]["completed_utc"]
    replayed = scheduler_time(timestamp + timedelta(minutes=10))
    value["accounting"]["completed_utc"] = replayed
    write_json(manifest, value)
    ledger.write_text(ledger.read_text().replace(original, replayed, 1))
    ledger.chmod(0o644)
    rewrite_scheduler_from_manifest(recost_fixture)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "completes after the request")


def test_generator_rejects_internally_replayed_ledger_provenance(recost_fixture):
    ledger = recost_fixture["ledger"]
    timestamp = recost_fixture["timestamp"]
    assert isinstance(ledger, Path)
    assert isinstance(timestamp, datetime)
    original = scheduler_time(timestamp - timedelta(hours=2))
    replayed = scheduler_time(timestamp - timedelta(hours=3))
    ledger.write_text(ledger.read_text().replace(original, replayed, 1))
    ledger.chmod(0o644)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "accounting submitted_utc differs")


def test_generator_rejects_ledger_allocation_arithmetic_replay(recost_fixture):
    ledger = recost_fixture["ledger"]
    assert isinstance(ledger, Path)
    ledger.write_text(ledger.read_text().replace("3600,2.000000,1.000000", "3600,3.000000,1.000000"))
    ledger.chmod(0o644)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "reserved use differs from allocation arithmetic")


@pytest.mark.parametrize(
    ("retained", "replayed", "message"),
    (
        (
            "2.000000,1.000000,1.000000",
            "2.000001,1.000000,1.000000",
            "reserved use differs from allocation arithmetic",
        ),
        (
            "2.000000,1.000000,1.000000",
            "2.000000,1.000001,1.000000",
            "actual use differs from scheduler arithmetic",
        ),
        (
            "2.000000,1.000000,1.000000",
            "2.000000,1.000000,1.000001",
            "cumulative use differs from ledger arithmetic",
        ),
        (
            "2.000000,1.000000,1.000000",
            "2.0,1.000000,1.000000",
            "reserved use must be a canonical six-decimal node-hour string",
        ),
    ),
)
def test_generator_rejects_last_decimal_and_noncanonical_ledger_replays(
    recost_fixture, retained, replayed, message
):
    ledger = recost_fixture["ledger"]
    assert isinstance(ledger, Path)
    ledger.write_text(ledger.read_text().replace(retained, replayed, 1))
    ledger.chmod(0o644)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), message)


def test_generator_rejects_last_decimal_reservation_replay(recost_fixture):
    reservations = recost_fixture["reservations"]
    assert isinstance(reservations, Path)
    value = json.loads(reservations.read_text())
    value[0]["actual_node_hours"] = 1.000001
    write_json(reservations, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "reservation actual use differs for job 12345")


def test_generator_accepts_raw_scheduler_reservation_after_exact_quantization(
    recost_fixture,
):
    replace_barrier_with_failed_outcome(recost_fixture)
    reservations = recost_fixture["reservations"]
    assert isinstance(reservations, Path)
    value = json.loads(reservations.read_text())
    value[-1]["actual_node_hours"] = 2 * 600 / 3600.0
    write_json(reservations, value)
    refresh_request(recost_fixture)
    completed = run_generator(recost_fixture)
    assert completed.returncode == 0, completed.stderr


def test_generator_rejects_last_decimal_cumulative_drift_across_rows(recost_fixture):
    add_recorded_wave_job(recost_fixture)
    ledger = recost_fixture["ledger"]
    assert isinstance(ledger, Path)
    ledger.write_text(
        ledger.read_text().replace(
            "2.000000,1.000000,2.000000,",
            "2.000000,1.000000,2.000001,",
            1,
        )
    )
    ledger.chmod(0o644)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "cumulative use differs from ledger arithmetic")


def test_generator_rejects_reordered_barrier_suffix(recost_fixture):
    add_recorded_wave_job(recost_fixture)
    ledger = recost_fixture["ledger"]
    assert isinstance(ledger, Path)

    add_fixture = recost_fixture

    def reverse_scheduler(request):
        request["inputs"]["scheduler_evidence"].reverse()

    refresh_request(add_fixture, reverse_scheduler)
    assert_rejected(run_generator(add_fixture), "scheduler evidence order differs")

    def reverse_barrier(request):
        request["inputs"]["scheduler_evidence"].reverse()
        request["barrier"]["recorded_segments"].reverse()

    refresh_request(add_fixture, reverse_barrier)
    assert_rejected(run_generator(add_fixture), "not the exact recorded ledger suffix")


@pytest.mark.parametrize(
    ("field", "value", "pattern"),
    [
        ("ranks_per_node", 7, "ranks per node differ"),
        ("cpus_per_task", 8, "CPUs per task differ"),
        ("controller_walltime_max_seconds", 7000, "controller walltime maximum differs"),
        ("athena_walltime", "01:55:01", "ten-minute shutdown margin"),
        ("output_layout", "shared", "output layout must be rank-local"),
        ("acceptance_policy", "U+A+F", "acceptance policy differs"),
        ("acceptance_criterion", "x" * 200, "acceptance criterion differs"),
    ],
)
def test_generator_rejects_inexact_execution_profile_fields(
    recost_fixture, field, value, pattern
):
    def mutate(request):
        request["recommendations"]["profiles"][0][field] = value

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), pattern)


def test_generator_rejects_executable_build_and_input_replays(recost_fixture):
    def mutate_executable(request):
        request["recommendations"]["profiles"][0]["executable_sha256"] = "0" * 64

    refresh_request(recost_fixture, mutate_executable)
    assert_rejected(run_generator(recost_fixture), "executable checksum differs")


def test_generator_rejects_build_environment_revision_replay(recost_fixture):
    environment = recost_fixture["build_manifest"] / "environment.txt"
    assert isinstance(environment, Path)
    environment.write_text("git_revision=" + "0" * 40 + "\n")
    environment.chmod(0o644)

    def mutate(request):
        request["recommendations"]["profiles"][0]["build_manifest_sha256"] = (
            directory_inventory_sha256(recost_fixture["build_manifest"])
        )

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "build environment revision differs")


def test_generator_rejects_input_and_absolute_path_replays(recost_fixture):
    def mutate_input(request):
        request["recommendations"]["profiles"][0]["input_file"] = (
            "inputs/cgl_lf_paper/fixture_R04.athinput"
        )

    refresh_request(recost_fixture, mutate_input)
    assert_rejected(run_generator(recost_fixture), "input differs from the matrix")

    root = recost_fixture["root"]
    assert isinstance(root, Path)

    def mutate_path(request):
        request["recommendations"]["profiles"][0]["executable"] = str(
            root / "build/frontier-fixture/../frontier-fixture/src/athena"
        )

    refresh_request(recost_fixture, mutate_path)
    assert_rejected(run_generator(recost_fixture), "must be an absolute normalized path")


def test_generator_rejects_unbounded_plan_target_and_segment_lineage(recost_fixture):
    continuation = profile(
        recost_fixture,
        target=1.0,
        segment="s01_rankio_t0p25_t1",
    )
    set_profiles(recost_fixture, [continuation], mode="sole-next-profile")
    assert_rejected(run_generator(recost_fixture), "increment exceeds the plan")

    fresh = profile(
        recost_fixture,
        case_id="R04",
        parent=False,
        target=0.5,
        segment="s00_rankio_t0_t0p5",
    )
    set_profiles(recost_fixture, [fresh], mode="sole-next-profile")
    assert_rejected(run_generator(recost_fixture), "target differs from the plan")


def test_generator_authenticates_complete_terminal_restart_group(recost_fixture):
    manifest = recost_fixture["manifest"]
    assert isinstance(manifest, Path)
    value = json.loads(manifest.read_text())
    second = Path(value["scientific_inspection"]["terminal_restart"]["rank_files"][1]["path"])
    second.write_bytes(b"replayed restart bytes\n")
    second.chmod(0o644)
    assert_rejected(run_generator(recost_fixture), "terminal restart rank 1 checksum differs")


def test_generator_rejects_self_asserted_storage_estimate(recost_fixture):
    def mutate(request):
        request["recommendations"]["profiles"][0]["estimated_storage_bytes"] += 1

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "differs from observed-rate projection")


def test_source_bundle_is_complete_and_pathname_races_fail_closed(
    recost_fixture, monkeypatch
):
    module = load_recost_module()
    repository = recost_fixture["repository"]
    bundle = recost_fixture["source_bundle"]
    revision = recost_fixture["revision"]
    assert isinstance(repository, Path)
    assert isinstance(bundle, Path)
    assert isinstance(revision, str)

    original_run = module.subprocess.run
    replacement = bundle.with_name("replacement.bundle")
    replacement.write_bytes(b"not a source bundle\n")
    replacement.chmod(0o644)
    retained = bundle.with_name("retained-original.bundle")
    expected = sha256(bundle)
    raced = False

    def race(command, *args, **kwargs):
        nonlocal raced
        if not raced and "bundle" in command and "verify" in command:
            bundle.rename(retained)
            replacement.rename(bundle)
            raced = True
        return original_run(command, *args, **kwargs)

    monkeypatch.setattr(module.subprocess, "run", race)
    with pytest.raises(ValueError, match="pathname changed during verification"):
        module.require_valid_git_bundle(repository, bundle, expected, [revision])
    assert raced


def test_source_bundle_rejects_prerequisite_only_history(recost_fixture):
    module = load_recost_module()
    repository = recost_fixture["repository"]
    root = recost_fixture["root"]
    revision = recost_fixture["revision"]
    assert isinstance(repository, Path)
    assert isinstance(root, Path)
    assert isinstance(revision, str)
    new_source = repository / "new-source.txt"
    new_source.write_text("new source history\n")
    git(repository, "add", str(new_source.relative_to(repository)))
    git(
        repository,
        "-c",
        "user.name=CGL fixture",
        "-c",
        "user.email=cgl-fixture@example.invalid",
        "commit",
        "-q",
        "-m",
        "Create prerequisite-only bundle fixture",
    )
    head = git(repository, "rev-parse", "HEAD")
    bundle = root / "source-archives/prerequisite.bundle"
    git(repository, "bundle", "create", str(bundle), "HEAD", f"^{revision}")
    bundle.chmod(0o644)
    with pytest.raises(ValueError, match="self-contained without prerequisites"):
        module.require_valid_git_bundle(repository, bundle, sha256(bundle), [head])


def test_source_bundle_requires_exact_advertised_head_selection(recost_fixture):
    module = load_recost_module()
    repository = recost_fixture["repository"]
    bundle = recost_fixture["source_bundle"]
    revision = recost_fixture["revision"]
    assert isinstance(repository, Path)
    assert isinstance(bundle, Path)
    assert isinstance(revision, str)
    module.require_valid_git_bundle(
        repository,
        bundle,
        sha256(bundle),
        [revision],
        expected_advertised_tip=(revision, "HEAD"),
    )

    master_only = bundle.with_name("master-only.bundle")
    git(repository, "bundle", "create", str(master_only), "refs/heads/master")
    master_only.chmod(0o644)
    with pytest.raises(ValueError, match="advertised tip differs"):
        module.require_valid_git_bundle(
            repository,
            master_only,
            sha256(master_only),
            [revision],
            expected_advertised_tip=(revision, "HEAD"),
        )

    git(repository, "commit", "--allow-empty", "-q", "-m", "unreviewed descendant")
    descendant = git(repository, "rev-parse", "HEAD")
    descendant_bundle = bundle.with_name("unreviewed-descendant.bundle")
    git(repository, "bundle", "create", str(descendant_bundle), "HEAD")
    descendant_bundle.chmod(0o644)
    with pytest.raises(ValueError, match="advertised tip differs"):
        module.require_valid_git_bundle(
            repository,
            descendant_bundle,
            sha256(descendant_bundle),
            [revision],
            expected_advertised_tip=(revision, "HEAD"),
        )
    assert descendant != revision


@pytest.mark.parametrize(
    "fixture_key", ("historical_source_bundle", "f116_source_bundle")
)
def test_f118_bridge_and_predecessor_bundles_bind_exact_live_bytes(
    recost_fixture, fixture_key,
):
    bundle = recost_fixture[fixture_key]
    assert isinstance(bundle, Path)
    bundle.write_bytes(bundle.read_bytes() + b"\nmutated live bundle bytes\n")
    bundle.chmod(0o644)
    assert_rejected(run_generator(recost_fixture), "checksum differs")


def test_generator_rejects_replayed_terminal_rank_inventory(recost_fixture):
    manifest = recost_fixture["manifest"]
    assert isinstance(manifest, Path)
    value = json.loads(manifest.read_text())
    rank_file = value["scientific_inspection"]["terminal_restart"]["rank_files"][1]
    original = Path(rank_file["path"])
    replay = original.parent.parent / "rank_99999999" / original.name
    replay.parent.mkdir()
    replay.write_bytes(original.read_bytes())
    replay.chmod(0o644)
    rank_file.update(
        {"path": str(replay), "sha256": sha256(replay), "size_bytes": replay.stat().st_size}
    )
    write_json(manifest, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "terminal restart rank inventory differs")


def test_generator_rejects_broken_recorded_parent_lineage(recost_fixture):
    manifest = recost_fixture["manifest"]
    root = recost_fixture["root"]
    assert isinstance(manifest, Path)
    assert isinstance(root, Path)
    value = json.loads(manifest.read_text())
    value["command"]["parent_segment"] = {
        "case_id": "R03",
        "manifest": str(root / "runs/mks24-stage-i" / EPOCH / "missing.json"),
        "segment": "s99_rankio_t0_t0p1",
        "result": "clean_partial",
        "final_time": 0.1,
    }
    write_json(manifest, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "recorded lineage parent is missing")


def test_authenticated_lineage_rejects_endpoint_replay_and_orphan_root(recost_fixture):
    module = load_recost_module()
    root = recost_fixture["root"]
    manifest = recost_fixture["manifest"]
    assert isinstance(root, Path)
    assert isinstance(manifest, Path)

    replayed = json.loads(manifest.read_text())
    replayed["_manifest_path"] = str(manifest)
    replayed["_manifest_sha256"] = sha256(manifest)
    replayed["scientific_inspection"]["final_time"] = 0.75
    with pytest.raises(ValueError, match="final time is outside its segment interval"):
        module.authenticated_case_lineages(root, [replayed])

    orphan = json.loads(manifest.read_text())
    orphan["_manifest_path"] = str(manifest)
    orphan["_manifest_sha256"] = sha256(manifest)
    orphan_segment = "s01_rankio_t0p25_t0p5"
    orphan["run"]["segment"] = orphan_segment
    orphan["accounting"]["segment"] = orphan_segment
    orphan["accounting"]["result"] = "accepted"
    orphan["scientific_inspection"]["segment"] = orphan_segment
    orphan["scientific_inspection"]["final_time"] = 0.5
    orphan["scientific_inspection"]["accepted"] = True
    orphan["scientific_inspection"]["terminal_restart_time"] = 0.5
    with pytest.raises(ValueError, match="recorded lineage root is not s00 from t=0"):
        module.authenticated_case_lineages(root, [orphan])


def test_generator_rejects_nonlatest_malformed_and_nonpreceding_predecessor(recost_fixture):
    accounting = recost_fixture["accounting"]
    timestamp = recost_fixture["timestamp"]
    assert isinstance(accounting, Path)
    assert isinstance(timestamp, datetime)
    newer = accounting / "mks24_stage_i_E03_forcing_policy_F198_recost_evidence.json"
    predecessor_value = json.loads(recost_fixture["predecessor"].read_text())
    predecessor_value["checkpoint"] = "F-198"
    predecessor_value["artifact_name"] = newer.name
    write_immutable_json(newer, predecessor_value)
    newer_review = accounting / f"{newer.name}.independent_review.json"
    review_value = json.loads(recost_fixture["predecessor_review"].read_text())
    review_value["reviewed_utc"] = (timestamp - timedelta(minutes=2)).isoformat()
    review_value["candidate"] = {"path": str(newer), "sha256": sha256(newer)}
    write_immutable_json(newer_review, review_value)
    newer_audit = accounting / f"{newer.name}.publication_audit.json"
    write_immutable_json(
        newer_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=1)).isoformat(),
            "artifact": {
                "path": str(newer),
                "sha256": sha256(newer),
                "mode": "0444",
                "links": 1,
            },
            "independent_review": {
                "path": str(newer_review),
                "sha256": sha256(newer_review),
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
    assert_rejected(run_generator(recost_fixture), "predecessor recost is not the latest")

    newer.unlink()
    newer_review.unlink()
    newer_audit.unlink()
    malformed = accounting / "mks24_stage_i_E03_forcing_policy_F197_recost_evidence.json"
    write_json(malformed, {"record_type": "fixture-recost"})
    malformed_audit = accounting / f"{malformed.name}.publication_audit.json"
    malformed_audit.write_text("{")
    malformed_audit.chmod(0o644)
    assert_rejected(run_generator(recost_fixture), "is not valid JSON")

    malformed.unlink()
    malformed_audit.unlink()
    predecessor = recost_fixture["predecessor"]
    assert isinstance(predecessor, Path)
    value = json.loads(predecessor.read_text())
    value["generated_utc"] = timestamp.isoformat()
    write_json(predecessor, value)
    refresh_predecessor(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "publication predates its artifact")

    value["generated_utc"] = (timestamp - timedelta(minutes=8)).isoformat()
    value["checkpoint"] = "F-200"
    write_json(predecessor, value)
    refresh_predecessor(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "artifact name and checkpoint ID differ")


def test_generator_rejects_legacy_predecessor_without_strict_publication_chain(recost_fixture):
    old = recost_fixture["predecessor"]
    old_audit = recost_fixture["predecessor_audit"]
    accounting = recost_fixture["accounting"]
    timestamp = recost_fixture["timestamp"]
    root = recost_fixture["root"]
    assert all(isinstance(item, Path) for item in (old, old_audit, accounting, root))
    assert isinstance(timestamp, datetime)
    old.unlink()
    old_audit.unlink()
    predecessor = accounting / "mks24_stage_i_E03_forcing_policy_legacy_recost_evidence.json"
    write_json(
        predecessor,
        {
            "record_type": "stage-i-legacy-recost-checkpoint",
            "checkpoint": "F-200",
            "execution_epoch": EPOCH,
            "generated_utc": (timestamp - timedelta(minutes=8)).isoformat(),
        },
    )
    audit = accounting / f"{predecessor.name}.publication_audit.json"
    write_json(
        audit,
        {
            "record_type": "observed-publication",
            "execution_epoch": EPOCH,
            "published_utc": (timestamp - timedelta(minutes=5)).isoformat(),
            "artifact": {
                "path": str(predecessor),
                "sha256": sha256(predecessor),
                "mode": "0644",
                "links": 1,
            },
        },
    )
    recost_fixture["predecessor"] = predecessor
    recost_fixture["predecessor_audit"] = audit

    def mutate(request):
        request["inputs"]["predecessor_recost"]["path"] = str(predecessor.relative_to(root))
        request["inputs"]["predecessor_recost_publication_audit"]["path"] = str(
            audit.relative_to(root)
        )

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "predecessor recost mode is 0644, expected 0444")


@pytest.mark.parametrize(
    ("mutation", "pattern"),
    [
        (lambda r, now: r.__setitem__("generated_utc", (now - timedelta(days=2)).isoformat()), "stale"),
        (lambda r, now: r.__setitem__("expires_utc", r["generated_utc"]), "expiry must follow"),
        (
            lambda r, now: r.__setitem__(
                "expires_utc", (now + timedelta(hours=25)).isoformat()
            ),
            "lifetime exceeds 24 hours",
        ),
        (
            lambda r, now: r.update(
                {
                    "generated_utc": (now - timedelta(hours=2)).isoformat(),
                    "expires_utc": (now - timedelta(hours=1)).isoformat(),
                }
            ),
            "expired",
        ),
    ],
)
def test_generator_rejects_invalid_generation_and_expiry(recost_fixture, mutation, pattern):
    now = recost_fixture["timestamp"]
    assert isinstance(now, datetime)
    refresh_request(recost_fixture, lambda request: mutation(request, now))
    assert_rejected(run_generator(recost_fixture), pattern)


def test_generator_rejects_storage_replays_and_exhaustion(recost_fixture):
    storage = recost_fixture["storage"]
    root = recost_fixture["root"]
    assert isinstance(storage, Path)
    assert isinstance(root, Path)
    value = json.loads(storage.read_text())
    value["retained_stage_i_bytes"] += 1
    write_json(storage, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "retained-byte count differs")

    refresh_storage(recost_fixture, available=live_available(root) + 1_000_000_000_000)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "exceeds live available bytes")

    refresh_storage(recost_fixture, available=2047)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "storage headroom is exhausted")


def test_generator_independently_rejects_budget_exhaustion(recost_fixture):
    ceiling = recost_fixture["ceiling"]
    assert isinstance(ceiling, Path)
    value = json.loads(ceiling.read_text())
    value["promoted_controls"]["campaign_budget_node_hours"] = 1.0
    write_json(ceiling, value)
    refresh_ceiling(recost_fixture)
    assert_rejected(run_generator(recost_fixture), "actual use plus authorized wave exceeds")


def test_generator_enforces_four_lanes_ten_nodes_and_r17_exclusive_last(recost_fixture):
    five = [
        profile(recost_fixture, case_id=case_id, parent=False)
        for case_id in ("R04", "R05", "R06", "R07", "R08")
    ]
    set_profiles(recost_fixture, five, mode="bounded-wave")
    assert_rejected(run_generator(recost_fixture), "promoted lane limit")

    twelve = [
        profile(recost_fixture, case_id=case_id, parent=False, nodes=4)
        for case_id in ("R04", "R05", "R06")
    ]
    set_profiles(recost_fixture, twelve, mode="bounded-wave")
    assert_rejected(run_generator(recost_fixture), "10-node wave ceiling")

    r17 = profile(
        recost_fixture,
        case_id="R17",
        parent=False,
        nodes=8,
        estimated_storage_bytes=958_271_710_272,
    )
    r04 = profile(recost_fixture, case_id="R04", parent=False)
    set_profiles(recost_fixture, [r17, r04], mode="bounded-wave")
    assert_rejected(run_generator(recost_fixture), "R17 authorization must be exclusive")
    set_profiles(recost_fixture, [r17], mode="sole-next-profile")
    assert_rejected(run_generator(recost_fixture), "R17 must remain last")


def test_generator_rejects_non_r17_and_accepts_reviewed_r17_operational_readiness(
    recost_fixture,
):
    root = recost_fixture["root"]
    assert isinstance(root, Path)
    readiness = root / "accounting" / R17_READINESS_NAME
    write_json(readiness, {"fixture": True})

    def mutate(request):
        request["inputs"]["r17_readiness_evidence"] = {
            "path": str(readiness.relative_to(root)),
            "sha256": sha256(readiness),
        }

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "non-R17 recost must not retain")

    pytest.skip("superseded by the real producer-output integration test")

    module = load_recost_module()
    now = datetime.now(timezone.utc).replace(microsecond=0)
    build_manifest = recost_fixture["build_manifest"]
    assert isinstance(build_manifest, Path)
    build_inventory = [
        {
            "name": path.name,
            "mode": f"{stat.S_IMODE(path.stat().st_mode):04o}",
            "sha256": sha256(path),
        }
        for path in sorted(build_manifest.iterdir())
    ]
    build_inventory_sha256 = hashlib.sha256(
        (json.dumps(build_inventory, sort_keys=True) + "\n").encode()
    ).hexdigest()
    scheduler = root / "accounting/98765.r17_qualification.sacct.txt"
    submitted_text = (now - timedelta(minutes=20)).isoformat()
    started_text = submitted_text
    completed_text = (now - timedelta(minutes=10)).isoformat()
    scheduler.write_text(
        "98765|r17_qualification|COMPLETED|0:0|8|600|"
        f"{submitted_text}|{completed_text}\n"
    )
    scheduler.chmod(0o644)
    outputs = []
    restarts = []
    for rank in range(64):
        rank_dir = root / "runs/r17-fixture" / f"rank_{rank:08d}"
        rank_dir.mkdir(parents=True)
        output_file = rank_dir / "output.bin"
        restart_file = rank_dir / "restart.rst"
        output_file.write_bytes(f"output {rank}\n".encode())
        restart_file.write_bytes(f"restart {rank}\n".encode())
        output_file.chmod(0o644)
        restart_file.chmod(0o644)
        outputs.append(
            {"path": output_file.relative_to(root).as_posix(), "sha256": sha256(output_file)}
        )
        restarts.append(
            {"path": restart_file.relative_to(root).as_posix(), "sha256": sha256(restart_file)}
        )
    restart_load = root / "accounting/r17_restart_load_evidence.json"
    physics = root / "accounting/r17_physics_validation_evidence.json"
    account_scheduler = root / "accounting/98765.r17_qualification.account.sacct.txt"
    account_exclusivity = root / (
        f"accounting/mks24_stage_i_{EPOCH_SLUG}_R17_account_exclusivity_evidence.json"
    )
    output_inventory_sha256 = hashlib.sha256(
        (json.dumps(outputs, sort_keys=True) + "\n").encode()
    ).hexdigest()
    restart_inventory_sha256 = hashlib.sha256(
        (json.dumps(restarts, sort_keys=True) + "\n").encode()
    ).hexdigest()
    write_immutable_json(
        restart_load,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-restart-load-evidence",
            "execution_epoch": EPOCH,
            "measured_utc": (now - timedelta(minutes=8)).isoformat(),
            "measured_by": "fixture restart-load measurement agent",
            "job_id": "98765",
            "executable_sha256": "4" * 64,
            "build_manifest_inventory_sha256": build_inventory_sha256,
            "passed": True,
            "measurements": {
                "rank_local_restart_inventory_sha256": restart_inventory_sha256,
                "loaded_rank_count": 64,
                "load_state": "COMPLETED",
                "load_exit_code": "0:0",
            },
        },
    )
    logical_locations = [
        [lx1, lx2, lx3, 0]
        for lx1 in range(12)
        for lx2 in range(12)
        for lx3 in range(12)
    ]
    block_rank_inventory = [
        {
            "rank": rank,
            "rank_name": f"rank_{rank:08d}",
            "logical_meshblocks": logical_locations[rank * 27:(rank + 1) * 27],
        }
        for rank in range(64)
    ]
    decomposition_value = {
        "schema_version": 1,
        "record_type": "stage-i-r17-decomposition-evidence",
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "logical_meshblock_grid": [12, 12, 12],
        "logical_meshblocks": 1728,
        "ranks": 64,
        "meshblocks_per_rank": 27,
        "complete_block_rank_inventory": block_rank_inventory,
        "complete_block_rank_inventory_sha256": hashlib.sha256(
            json.dumps(
                block_rank_inventory, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest(),
        "terminal_rank_local_output_inventory_sha256": output_inventory_sha256,
        "checks": {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        },
    }
    physics_value = {
        "schema_version": 1,
        "record_type": "stage-i-r17-physics-validation-evidence",
        "execution_epoch": EPOCH,
        "measured_utc": (now - timedelta(minutes=8)).isoformat(),
        "measured_by": "fixture plasma measurement agent",
        "job_id": "98765",
        "executable_sha256": "4" * 64,
        "build_manifest_inventory_sha256": build_inventory_sha256,
        "passed": True,
        "measurements": {
            "rank_local_output_inventory_sha256": output_inventory_sha256,
            "finite_rank_outputs": 64,
            "mass_relative_drift_max": "0",
            "mhd_user_mass_mismatch_max": "0",
            "lf_bad_counts_total": 0,
            "normalized_ct_divb_max": "1e-14",
            "normalized_ct_divb_threshold": "1e-12",
            "normalized_ct_divb_below_threshold": True,
        },
    }
    write_immutable_json(physics, physics_value)
    account_scheduler.write_text(
        module.ACCOUNT_SCHEDULER_HEADER
        + "\n"
        + "|".join(
            [
                "98765",
                "r17_qualification",
                "COMPLETED",
                "0:0",
                "8",
                "600",
                submitted_text,
                started_text,
                completed_text,
                "batch",
                "AST207",
                "fixture",
            ]
        )
        + "\n"
    )
    account_scheduler.chmod(0o644)
    account_jobs = module.parse_account_scheduler_text(account_scheduler.read_bytes())
    qualification_job = {
        key: account_jobs[0][key]
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
    account_exclusivity_value = {
        "schema_version": 1,
        "record_type": "stage-i-r17-account-exclusivity-evidence",
        "execution_epoch": EPOCH,
        "measured_utc": (now - timedelta(minutes=8)).isoformat(),
        "query_contract": module.expected_account_query_contract(
            now - timedelta(minutes=20), now - timedelta(minutes=10)
        ),
        "visibility_contract": {
            "private_data": "none",
            "all_users_job_visibility": True,
        },
        "raw_account_scheduler_sha256": sha256(account_scheduler),
        "qualification_job": qualification_job,
        "qualification_job_sha256": module.compact_json_sha256(qualification_job),
        "account_jobs": account_jobs,
        "account_jobs_sha256": module.compact_json_sha256(account_jobs),
        "overlapping_job_ids": ["98765"],
        "exclusive_entire_execution_interval": True,
    }
    write_immutable_json(account_exclusivity, account_exclusivity_value)
    qualification = root / "accounting/r17_operational_qualification.json"
    write_json(
        qualification,
        {
            "schema_version": 2,
            "record_type": "stage-i-r17-operational-qualification",
            "execution_epoch": EPOCH,
            "completed_utc": completed_text,
            "job_id": "98765",
            "scheduler_evidence": {
                "path": str(scheduler.relative_to(root)),
                "sha256": sha256(scheduler),
            },
            "account_scheduler_evidence": {
                "path": account_scheduler.relative_to(root).as_posix(),
                "sha256": sha256(account_scheduler),
            },
            "account_exclusivity_evidence": {
                "path": account_exclusivity.relative_to(root).as_posix(),
                "sha256": sha256(account_exclusivity),
            },
            "state": "COMPLETED",
            "exit_code": "0:0",
            "nodes": 8,
            "ranks": 64,
            "executable_sha256": "4" * 64,
            "build_manifest_inventory": build_inventory,
            "rank_local_outputs": outputs,
            "rank_local_restarts": restarts,
            "decomposition_evidence": decomposition_value,
            "restart_load_evidence": {
                "path": restart_load.relative_to(root).as_posix(),
                "sha256": sha256(restart_load),
            },
            "physics_validation_evidence": {
                "path": physics.relative_to(root).as_posix(),
                "sha256": sha256(physics),
            },
        },
    )
    qualification_review = root / "accounting/r17_operational_qualification.json.independent_review.json"
    write_immutable_json(
        qualification_review,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-operational-qualification-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (now - timedelta(minutes=5)).isoformat(),
            "decision": "approved",
            "reviewer": "independent operational reviewer",
            "candidate": {"path": str(qualification), "sha256": sha256(qualification)},
        },
    )
    value = {
        "schema_version": 1,
        "record_type": "stage-i-r17-readiness",
        "execution_epoch": EPOCH,
        "root": str(root),
        "generated_utc": (now - timedelta(minutes=1)).isoformat(),
        "expires_utc": (now + timedelta(hours=2)).isoformat(),
        "reviewed_by": "fixture readiness reviewer",
        "predecessor_lineages_sha256": "1" * 64,
        "storage_evidence_sha256": "2" * 64,
        "computed_projection_sha256": "3" * 64,
        "executable_sha256": "4" * 64,
        "build_manifest_sha256": build_inventory_sha256,
        "required_retained_bytes": 958_271_710_272,
        "nodes": 8,
        "ranks": 64,
        "storage_ready": True,
        "node_hour_ready": True,
        "rank_64_ready": True,
        "operational_qualification": {
            "path": str(qualification.relative_to(root)),
            "sha256": sha256(qualification),
        },
        "operational_qualification_review": {
            "path": str(qualification_review.relative_to(root)),
            "sha256": sha256(qualification_review),
        },
    }
    tracker = module.InputTracker()
    parsed = module.parse_r17_readiness(
        value,
        root,
        now,
        now + timedelta(hours=1),
        "1" * 64,
        "2" * 64,
        "3" * 64,
        {
            "executable_sha256": "4" * 64,
            "build_manifest": str(build_manifest),
            "build_manifest_sha256": build_inventory_sha256,
        },
        tracker,
    )
    assert parsed["operational_qualification_evidence"]["job_id"] == "98765"
    assert (
        parsed["operational_qualification_evidence"][
            "build_manifest_inventory_sha256"
        ]
        == build_inventory_sha256
    )
    assert (
        parsed["operational_qualification_evidence"][
            "authenticated_decomposition_evidence"
        ]["logical_meshblocks"]
        == 1728
    )
    assert (
        parsed["operational_qualification_evidence"][
            "authenticated_account_exclusivity_evidence"
        ]["overlapping_job_ids"]
        == ["98765"]
    )
    invalid_outputs = json.loads(qualification.read_text())
    invalid_outputs["rank_local_outputs"].append(dict(outputs[0]))
    with pytest.raises(ValueError, match="must retain exactly 64 files"):
        module.parse_r17_operational_qualification(
            invalid_outputs,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    invalid_qualification = json.loads(qualification.read_text())
    inline = invalid_qualification["decomposition_evidence"]
    inline["complete_block_rank_inventory"][0]["logical_meshblocks"][0] = list(
        inline["complete_block_rank_inventory"][0]["logical_meshblocks"][1]
    )
    inline["complete_block_rank_inventory"][0]["logical_meshblocks"].sort()
    inline["complete_block_rank_inventory_sha256"] = module.compact_json_sha256(
        inline["complete_block_rank_inventory"]
    )
    with pytest.raises(ValueError, match="does not retain 27 unique blocks per rank"):
        module.parse_r17_operational_qualification(
            invalid_qualification,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    legacy_decomposition = json.loads(qualification.read_text())
    legacy_decomposition["meshblock_decomposition_evidence"] = {
        "path": "accounting/legacy-self-asserted-decomposition.json",
        "sha256": "0" * 64,
    }
    del legacy_decomposition["decomposition_evidence"]
    with pytest.raises(ValueError, match="R17 operational qualification schema differs"):
        module.parse_r17_operational_qualification(
            legacy_decomposition,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    invalid_account = json.loads(account_exclusivity.read_text())
    invalid_account["overlapping_job_ids"] = []
    write_immutable_json(account_exclusivity, invalid_account)
    invalid_qualification = json.loads(qualification.read_text())
    invalid_qualification["account_exclusivity_evidence"]["sha256"] = sha256(
        account_exclusivity
    )
    with pytest.raises(ValueError, match="account exclusivity evidence differs"):
        module.parse_r17_operational_qualification(
            invalid_qualification,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    write_immutable_json(account_exclusivity, account_exclusivity_value)
    invalid_inventory = json.loads(qualification.read_text())
    invalid_inventory["build_manifest_inventory"] = build_inventory[:-1]
    with pytest.raises(ValueError, match="build binding differs"):
        module.parse_r17_operational_qualification(
            invalid_inventory,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    invalid_ct = json.loads(physics.read_text())
    invalid_ct["measurements"]["normalized_ct_divb_below_threshold"] = False
    write_immutable_json(physics, invalid_ct)
    invalid_qualification = json.loads(qualification.read_text())
    invalid_qualification["physics_validation_evidence"]["sha256"] = sha256(physics)
    with pytest.raises(ValueError, match="physics-validation evidence measurements differ"):
        module.parse_r17_operational_qualification(
            invalid_qualification,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    write_immutable_json(physics, physics_value)
    invalid_physics = json.loads(physics.read_text())
    invalid_physics["measurements"]["mass_relative_drift_max"] = "0.01"
    write_immutable_json(physics, invalid_physics)
    invalid_qualification = json.loads(qualification.read_text())
    invalid_qualification["physics_validation_evidence"]["sha256"] = sha256(physics)
    with pytest.raises(ValueError, match="physics-validation evidence measurements differ"):
        module.parse_r17_operational_qualification(
            invalid_qualification,
            root,
            now,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            qualification,
            sha256(qualification),
            {
                "path": qualification_review.relative_to(root).as_posix(),
                "sha256": sha256(qualification_review),
            },
            module.InputTracker(),
        )
    write_immutable_json(physics, physics_value)
    first_restart = root / restarts[0]["path"]
    first_restart.write_bytes(b"tampered restart\n")
    with pytest.raises(ValueError, match="checksum differs"):
        module.parse_r17_readiness(
            value,
            root,
            now,
            now + timedelta(hours=1),
            "1" * 64,
            "2" * 64,
            "3" * 64,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            module.InputTracker(),
        )
    first_restart.write_bytes(b"restart 0\n")
    first_restart.chmod(0o644)

    write_immutable_json(readiness, value)
    readiness_review = readiness.with_name(f"{readiness.name}.independent_review.json")
    write_immutable_json(
        readiness_review,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-readiness-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": (now - timedelta(seconds=30)).isoformat(),
            "decision": "approved-for-publication",
            "reviewer": "fixture readiness reviewer",
            "candidate": {"path": str(readiness), "sha256": sha256(readiness)},
        },
    )
    readiness_audit = readiness.with_name(f"{readiness.name}.publication_audit.json")
    write_immutable_json(
        readiness_audit,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-readiness-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": now.isoformat(),
            "artifact": {
                "path": str(readiness),
                "sha256": sha256(readiness),
                "mode": "0444",
                "links": 1,
            },
            "independent_review": {
                "path": str(readiness_review),
                "sha256": sha256(readiness_review),
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
    chain = module.parse_r17_readiness_publication_chain(
        root,
        readiness,
        sha256(readiness),
        value,
        {"path": readiness_review.relative_to(root).as_posix(), "sha256": sha256(readiness_review)},
        {"path": readiness_audit.relative_to(root).as_posix(), "sha256": sha256(readiness_audit)},
        now,
        module.InputTracker(),
    )
    assert chain["reviewed_by"] == "fixture readiness reviewer"
    value["predecessor_lineages_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="reviewed bindings differ"):
        module.parse_r17_readiness(
            value,
            root,
            now,
            now + timedelta(hours=1),
            "1" * 64,
            "2" * 64,
            "3" * 64,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            module.InputTracker(),
        )
    value["predecessor_lineages_sha256"] = "1" * 64
    value["storage_ready"] = 1
    with pytest.raises(ValueError, match="flags must be true booleans"):
        module.parse_r17_readiness(
            value,
            root,
            now,
            now + timedelta(hours=1),
            "1" * 64,
            "2" * 64,
            "3" * 64,
            {
                "executable_sha256": "4" * 64,
                "build_manifest": str(build_manifest),
                "build_manifest_sha256": build_inventory_sha256,
            },
            module.InputTracker(),
        )


def test_r17_recost_consumes_real_schema2_producer_output(tmp_path, monkeypatch):
    """Feed the qualification producer's retained schema-2 output into recost."""

    producer_tests = load_qualification_test_module()
    producer = producer_tests.qualification
    producer_fixture = producer_tests.qualification_fixture.__wrapped__(
        tmp_path, monkeypatch
    )
    monkeypatch.setattr(producer, "RANKS_PER_NODE", 8)
    wave, _, evidence, scheduler = producer_tests.retain_complete_wave(
        producer_fixture, ("R17",), max_nodes=8
    )
    retention = producer.retain_r17_operational_qualification(
        evidence["R17"], "producer measurement agent", scheduler
    )
    root = producer_fixture["root"]
    qualification_path = root / retention["operational_qualification"]["path"]
    operational = json.loads(qualification_path.read_text())
    assert stat.S_IMODE(qualification_path.stat().st_mode) == 0o444
    for key in (
        "scheduler_evidence",
        "account_scheduler_evidence",
        "account_exclusivity_evidence",
        "restart_load_evidence",
        "physics_validation_evidence",
    ):
        assert stat.S_IMODE((root / operational[key]["path"]).stat().st_mode) == 0o444

    measured = datetime.fromisoformat(operational["measured_utc"].replace("Z", "+00:00"))
    reviewed = measured + timedelta(seconds=1)
    request_time = reviewed + timedelta(seconds=1)
    review_path = root / operational["independent_review_contract"]["path"]
    write_immutable_json(
        review_path,
        {
            "schema_version": 1,
            "record_type": "stage-i-r17-operational-qualification-independent-review",
            "execution_epoch": EPOCH,
            "reviewed_utc": reviewed.isoformat(),
            "decision": "approved",
            "reviewer": "independent producer-output reviewer",
            "candidate": {
                "path": str(qualification_path),
                "sha256": sha256(qualification_path),
            },
        },
    )
    packet = producer_tests.packet_map(wave)[("R17", 8)]
    provenance = wave["provenance"]
    profile = {
        "source_bundle": provenance["source_bundle"]["path"],
        "source_bundle_sha256": provenance["source_bundle"]["sha256"],
        "executable": provenance["executable"]["path"],
        "executable_revision": provenance["executable"]["revision"],
        "executable_sha256": provenance["executable"]["sha256"],
        "input_revision": provenance["source"]["revision"],
        "input_sha256": packet["execution_intent"]["input"]["sha256"],
        "build_manifest": provenance["build_manifest"]["path"],
        "build_manifest_sha256": provenance["build_manifest"]["inventory_sha256"],
        "nodes": 8,
        "ranks_per_node": 8,
        "cpus_per_task": 1,
        "time_tlim_target": 0.25,
    }
    module = load_recost_module()
    review_binding = {
        "path": review_path.relative_to(root).as_posix(),
        "sha256": sha256(review_path),
    }
    parsed = module.parse_r17_operational_qualification(
        operational,
        root,
        request_time,
        profile,
        provenance["matrix"]["sha256"],
        qualification_path,
        sha256(qualification_path),
        review_binding,
        module.InputTracker(),
    )
    assert parsed["authenticated_frozen_science_build_contract"] == (
        operational["frozen_science_build_contract"]
    )
    assert parsed["authenticated_decomposition_evidence"]["logical_meshblocks"] == 1728
    assert parsed["authenticated_account_exclusivity_evidence"][
        "exclusive_entire_execution_interval"
    ] is True

    reduced = {
        key: operational[key]
        for key in (
            "schema_version", "record_type", "execution_epoch", "completed_utc",
            "job_id", "state", "exit_code", "nodes", "ranks",
        )
    }
    with pytest.raises(ValueError, match="schema differs"):
        module.parse_r17_operational_qualification(
            reduced,
            root,
            request_time,
            profile,
            provenance["matrix"]["sha256"],
            qualification_path,
            sha256(qualification_path),
            review_binding,
            module.InputTracker(),
        )

    invalid = deepcopy(operational)
    invalid["frozen_science_build_contract"]["parameter_contract_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="frozen science/build binding differs"):
        module.parse_r17_operational_qualification(
            invalid,
            root,
            request_time,
            profile,
            provenance["matrix"]["sha256"],
            qualification_path,
            sha256(qualification_path),
            review_binding,
            module.InputTracker(),
        )

    scientific_path = root / operational["scientific_evidence"]["path"]
    scientific = json.loads(scientific_path.read_text())
    invalid_ct = deepcopy(scientific)
    invalid_ct["user_history"] = mutated_history_copy(
        Path(scientific["user_history"]["path"]),
        scientific_path.with_name("coherent-invalid-ct.user.hst"),
        "max_ndiv",
        1.0e-12,
    )
    with pytest.raises(ValueError, match="mass, LF, hard-bound, or CT gate failed"):
        module.authenticate_r17_physics_contract(
            invalid_ct,
            root,
            operational["rank_local_outputs"],
            operational["rank_local_restarts"],
            module.InputTracker(),
        )
    invalid_mass = deepcopy(scientific)
    invalid_mass["mhd_history"] = mutated_history_copy(
        Path(scientific["mhd_history"]["path"]),
        scientific_path.with_name("coherent-invalid-mass.mhd.hst"),
        "mass",
        2.01,
    )
    with pytest.raises(ValueError, match="mass, LF, hard-bound, or CT gate failed"):
        module.authenticate_r17_physics_contract(
            invalid_mass,
            root,
            operational["rank_local_outputs"],
            operational["rank_local_restarts"],
            module.InputTracker(),
        )

    scheduler_path = root / operational["scheduler_evidence"]["path"]
    scheduler_path.chmod(0o644)
    with pytest.raises(ValueError, match="mode is 0644, expected 0444"):
        module.parse_r17_operational_qualification(
            operational,
            root,
            request_time,
            profile,
            provenance["matrix"]["sha256"],
            qualification_path,
            sha256(qualification_path),
            review_binding,
            module.InputTracker(),
        )


def test_final_storage_boundary_rejects_run_tree_race(recost_fixture):
    module = load_recost_module()
    root = recost_fixture["root"]
    storage = recost_fixture["storage"]
    assert isinstance(root, Path)
    assert isinstance(storage, Path)
    reviewed = json.loads(storage.read_text())
    race = root / "runs/mks24-stage-i" / EPOCH / "race.bin"
    race.write_bytes(b"late storage mutation\n")
    race.chmod(0o644)
    with pytest.raises(ValueError, match="changed before staged output creation"):
        module.require_live_storage_boundary(
            root,
            reviewed["available_bytes"],
            reviewed["retained_stage_i_bytes"],
            reviewed["required_safety_bytes"],
            reviewed["projected_authorized_wave_growth_bytes"],
        )


def test_storage_measurement_validates_and_excludes_intentional_r02_bundle_symlink(
    tmp_path,
):
    module = load_recost_module()
    run_store = tmp_path / "runs/mks24-stage-i/E03-forcing-policy"
    target = run_store / "R02/s00/output/value.bin"
    target.parent.mkdir(parents=True)
    target.write_bytes(b"12345678")
    target.chmod(0o644)
    link = run_store / "bundles/R02/cases/case/value.bin"
    link.parent.mkdir(parents=True)
    link.symlink_to(target)
    assert module.directory_tree_regular_bytes(run_store, "fixture Stage I storage") == 8


def test_storage_measurement_rejects_noninventory_and_escaping_symlinks(tmp_path):
    module = load_recost_module()
    run_store = tmp_path / "runs/mks24-stage-i/E03-forcing-policy"
    target = run_store / "R02/s00/output/value.bin"
    target.parent.mkdir(parents=True)
    target.write_bytes(b"12345678")
    target.chmod(0o644)
    unauthorized = run_store / "R02/s00/output/link.bin"
    unauthorized.symlink_to(target)
    with pytest.raises(ValueError, match="unauthorized symlink"):
        module.directory_tree_regular_bytes(run_store, "fixture Stage I storage")
    unauthorized.unlink()
    outside = tmp_path / "outside.bin"
    outside.write_bytes(b"outside")
    outside.chmod(0o644)
    escaping = run_store / "bundles/R02/cases/case/value.bin"
    escaping.parent.mkdir(parents=True)
    escaping.symlink_to(outside)
    with pytest.raises(ValueError, match="unsafe R02 bundle symlink"):
        module.directory_tree_regular_bytes(run_store, "fixture Stage I storage")


def test_locked_main_rechecks_mutable_boundaries_after_storage(monkeypatch, tmp_path):
    module = load_recost_module()
    events = []

    class Tracker:
        def reauthenticate_all(self):
            events.append("inputs")

    build = module.BuildResult(
        payload=b"{}\n",
        tracker=Tracker(),
        reconcile={},
        helper_path=tmp_path / "helper",
        helper_sha256="1" * 64,
        helper_revision="2" * 40,
        matrix_path=tmp_path / "matrix",
        matrix_sha256="3" * 64,
        matrix_revision="4" * 40,
        storage_available_bytes=10_000,
        storage_retained_stage_i_bytes=1_000,
        storage_required_safety_bytes=1_000,
        projected_storage_bytes=1_000,
        storage_measurements=(),
    )
    monkeypatch.setattr(module, "build_payload", lambda *args, **kwargs: build)
    monkeypatch.setattr(
        module, "require_empty_transaction_stores", lambda *args: events.append("transactions")
    )
    monkeypatch.setattr(module, "require_drained_queue", lambda *args: events.append("queue"))
    monkeypatch.setattr(
        module, "require_output_namespace_empty", lambda *args: events.append("namespace")
    )
    monkeypatch.setattr(module, "require_committed_file", lambda *args, **kwargs: "revision")
    monkeypatch.setattr(module, "run_authenticated_reconcile", lambda *args, **kwargs: {})
    monkeypatch.setattr(
        module, "require_live_storage_boundary", lambda *args: events.append("storage")
    )
    monkeypatch.setattr(module, "write_staged_output", lambda *args: events.append("write"))

    module.locked_main(
        object(),
        tmp_path,
        tmp_path / "artifact.json.staged",
        tmp_path / "generator",
        tmp_path / "repository",
        "5" * 64,
    )
    assert events[-5:] == ["storage", "transactions", "queue", "namespace", "write"]


def test_generator_rejects_barrier_suffix_replay(recost_fixture):
    add_recorded_wave_job(recost_fixture)

    def mutate(request):
        request["barrier"]["recorded_segments"] = request["barrier"]["recorded_segments"][:1]
        request["inputs"]["scheduler_evidence"] = request["inputs"]["scheduler_evidence"][:1]

    refresh_request(recost_fixture, mutate)
    assert_rejected(run_generator(recost_fixture), "not the exact recorded ledger suffix")


def test_generator_rejects_writable_directory_intermediate_symlink_and_output_symlink(
    recost_fixture,
):
    accounting = recost_fixture["accounting"]
    build_manifest = recost_fixture["build_manifest"]
    assert isinstance(accounting, Path)
    assert isinstance(build_manifest, Path)
    accounting.chmod(0o775)
    assert_rejected(run_generator(recost_fixture), "exceeds trusted profile 0755")
    accounting.chmod(0o755)

    source_archives = recost_fixture["root"] / "source-archives"
    source_archives.chmod(0o775)
    assert_rejected(run_generator(recost_fixture), "source archive directory mode")
    source_archives.chmod(0o755)

    target = build_manifest.with_name("fixture-build-target")
    build_manifest.rename(target)
    build_manifest.symlink_to(target, target_is_directory=True)
    assert_rejected(run_generator(recost_fixture), "path contains a symlink")

    build_manifest.unlink()
    target.rename(build_manifest)
    other = accounting / "mks24_stage_i_E03_forcing_policy_F200_recost_evidence.json.staged"
    symlink_target = accounting / "do-not-overwrite"
    symlink_target.write_text("retained\n")
    other.symlink_to(symlink_target)
    assert_rejected(run_generator(recost_fixture), "output namespace is not empty")
    assert symlink_target.read_text() == "retained\n"
    assert symlink_target.read_text() == "retained\n"


def test_generator_requires_itself_committed(recost_fixture):
    generator = recost_fixture["generator"]
    assert isinstance(generator, Path)
    generator.write_text(generator.read_text() + "\n# uncommitted fixture drift\n")
    generator.chmod(0o755)
    assert_rejected(
        run_generator(recost_fixture, generator_sha256=sha256(generator)),
        "must be committed before recost generation",
    )


def test_generator_rejects_preseeded_private_reexec_paths_without_valid_descriptor(
    recost_fixture,
):
    environment = dict(os.environ)
    environment.update(
        {
            "_CGL_LF_STAGE_I_RECOST_PYTHON_DESCRIPTOR": "99",
            "_CGL_LF_STAGE_I_RECOST_SOURCE": "/attacker/source.py",
            "_CGL_LF_STAGE_I_RECOST_REPOSITORY_ROOT": "/attacker/repository",
        }
    )
    assert_rejected(
        run_generator(recost_fixture, environment=environment),
        "private generator reexecution path is forbidden without an authenticated descriptor",
    )
    environment["_CGL_LF_STAGE_I_RECOST_DESCRIPTOR"] = "0"
    assert_rejected(
        run_generator(recost_fixture, environment=environment),
        "generator descriptor marker is not attached to this execution",
    )


def test_authenticated_reexec_descriptor_controls_private_source_metadata(
    recost_fixture, monkeypatch,
):
    module = load_recost_module()
    generator = recost_fixture["generator"]
    repository = recost_fixture["repository"]
    assert isinstance(generator, Path)
    assert isinstance(repository, Path)
    descriptor = os.open(generator, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(strict=True), os.O_RDONLY)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(generator))
        monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, str(repository))
        monkeypatch.setattr(module, "require_isolated_python_reexec", lambda value: None)
        assert module.authenticate_self(
            ["--expected-generator-sha256", sha256(generator)]
        ) == (generator, repository, sha256(generator))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(repository / "attacker.py"))
        with pytest.raises(ValueError, match="source/repository relationship differs"):
            module.authenticate_self(
                ["--expected-generator-sha256", sha256(generator)]
            )
    finally:
        os.close(descriptor)
        os.close(python_descriptor)


def test_authenticated_reexec_rejects_named_source_inode_substitution(
    recost_fixture, monkeypatch,
):
    module = load_recost_module()
    generator = recost_fixture["generator"]
    repository = recost_fixture["repository"]
    assert isinstance(generator, Path)
    assert isinstance(repository, Path)
    descriptor = os.open(generator, os.O_RDONLY)
    python_descriptor = os.open(Path("/proc/self/exe").resolve(strict=True), os.O_RDONLY)
    displaced = generator.with_name(f"{generator.name}.displaced")
    generator.rename(displaced)
    generator.write_bytes(displaced.read_bytes())
    generator.chmod(0o755)
    try:
        monkeypatch.setattr(module, "__file__", f"/proc/self/fd/{descriptor}")
        monkeypatch.setenv(module.SELF_DESCRIPTOR_ENV, str(descriptor))
        monkeypatch.setenv(module.PYTHON_DESCRIPTOR_ENV, str(python_descriptor))
        monkeypatch.setenv(module.SELF_SOURCE_ENV, str(generator))
        monkeypatch.setenv(module.REPOSITORY_ROOT_ENV, str(repository))
        monkeypatch.setattr(module, "require_isolated_python_reexec", lambda value: None)
        with pytest.raises(ValueError, match="source pathname is not the inherited descriptor"):
            module.authenticate_self(
                ["--expected-generator-sha256", sha256(generator)]
            )
    finally:
        os.close(descriptor)
        os.close(python_descriptor)


def test_authenticated_python_reexec_requires_isolated_running_interpreter():
    if sys.flags.isolated:
        pytest.skip("test runner is already isolated")
    module = load_recost_module()
    descriptor = os.open(Path("/proc/self/exe").resolve(strict=True), os.O_RDONLY)
    try:
        with pytest.raises(ValueError, match="Python reexecution is not isolated"):
            module.require_isolated_python_reexec(descriptor)
    finally:
        os.close(descriptor)


def test_reexec_and_git_environments_strip_caller_path_config_and_private_state(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    for key, value in {
        "PATH": str(tmp_path / "hostile-path"),
        "LD_PRELOAD": str(tmp_path / "hostile.so"),
        "PYTHONPATH": str(tmp_path / "hostile-python"),
        "GIT_CONFIG_GLOBAL": str(tmp_path / "hostile-global"),
        "GIT_CONFIG_SYSTEM": str(tmp_path / "hostile-system"),
        "GIT_CONFIG_PARAMETERS": "'core.repositoryformatversion=99'",
        module.SELF_DESCRIPTOR_ENV: "99",
        module.SELF_SOURCE_ENV: "/hostile/source",
        module.REPOSITORY_ROOT_ENV: "/hostile/repository",
    }.items():
        monkeypatch.setenv(key, value)
    reexec = module.reexec_environment(7, 8, RECOST, REPOSITORY)
    assert reexec[module.SELF_DESCRIPTOR_ENV] == "7"
    assert reexec[module.PYTHON_DESCRIPTOR_ENV] == "8"
    assert reexec[module.SELF_SOURCE_ENV] == str(RECOST)
    assert reexec[module.REPOSITORY_ROOT_ENV] == str(REPOSITORY)
    assert "LD_PRELOAD" not in reexec
    assert "PYTHONPATH" not in reexec
    assert "GIT_CONFIG_PARAMETERS" not in reexec
    assert reexec["XDG_CONFIG_HOME"] == "/nonexistent"

    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    module.git_run(tmp_path, ["status"], capture_output=True)
    assert len(calls) == 1
    command, kwargs = calls[0]
    assert command[:2] == [str(module.GIT), "--no-replace-objects"]
    assert kwargs["executable"].startswith("/proc/self/fd/")
    assert kwargs["pass_fds"]
    environment = kwargs["env"]
    assert environment == module.hardened_git_environment()
    assert environment["PATH"] == module.TRUSTED_SYSTEM_PATH
    assert environment["GIT_CONFIG_GLOBAL"] == os.devnull
    assert environment["GIT_CONFIG_SYSTEM"] == os.devnull
    assert environment["GIT_CONFIG_NOSYSTEM"] == "1"
    assert "LD_PRELOAD" not in environment
    assert "GIT_CONFIG_PARAMETERS" not in environment


def test_mutation_lock_pathname_replacement_blocks_managed_write(tmp_path):
    module = load_recost_module()
    root = tmp_path / "root"
    root.mkdir()
    lock_path = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock_path.write_bytes(b"")
    lock_path.chmod(0o644)
    displaced = root / "displaced.lock"
    target = root / "managed.json"
    with module.stage_i_lock(root) as mutation_lock:
        lock_path.rename(displaced)
        lock_path.write_bytes(b"")
        lock_path.chmod(0o644)
        with pytest.raises(ValueError, match="lock pathname changed"):
            module.write_exact_or_verify(
                target,
                b"{}\n",
                mode=0o644,
                label="managed fixture",
                mutation_lock=mutation_lock,
            )
        assert not target.exists()
        lock_path.unlink()
        displaced.rename(lock_path)


def publication_private_entries(module, parent: Path) -> list[Path]:
    """Return exact deterministic recost private-attempt entries."""

    return sorted(
        path
        for path in parent.iterdir()
        if module.PUBLICATION_PRIVATE_NAME_PATTERN.fullmatch(path.name) is not None
    )


def test_finalized_private_inode_no_replace_move_publication(tmp_path, monkeypatch):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    real_rename = module.renameat2_noreplace
    events = []

    def forbidden_link_or_unlink(*_args, **_kwargs):
        raise AssertionError("direct-final publication must not link or unlink")

    def tracked_rename(directory, source, destination, label):
        private = os.stat(source, dir_fd=directory, follow_symlinks=False)
        assert destination == target.name
        assert stat.S_IMODE(private.st_mode) == 0o644
        assert private.st_nlink == 1
        result = real_rename(directory, source, destination, label)
        public = os.stat(destination, dir_fd=directory, follow_symlinks=False)
        assert (public.st_dev, public.st_ino) == (private.st_dev, private.st_ino)
        assert public.st_nlink == 1
        events.append(("move", source))
        return result

    monkeypatch.setattr(module, "renameat2_noreplace", tracked_rename)
    monkeypatch.setattr(module.os, "link", forbidden_link_or_unlink)
    monkeypatch.setattr(module.os, "unlink", forbidden_link_or_unlink)
    assert module.write_exact_or_verify(
        target, payload, mode=0o644, label="managed fixture"
    )

    assert [event for event, _name in events] == ["move"]
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert target.stat().st_nlink == 1
    assert list(parent.iterdir()) == [target]


def test_publication_transaction_binds_target_payload_mode_producer_and_revision(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    target = tmp_path / "managed.json"
    base = module.publication_transaction(target, b"managed\n", 0o644, "managed fixture")
    assert base.target == str(target)
    assert base.payload_sha256 == hashlib.sha256(b"managed\n").hexdigest()
    assert base.payload_size == len(b"managed\n")
    assert base.final_mode == 0o644
    assert base.producer == "scripts/frontier/cgl_lf_stage_i_recost.py"
    assert re.fullmatch(r"[0-9a-f]{64}", base.producer_revision)
    assert module.PUBLICATION_PRIVATE_NAME_PATTERN.fullmatch(base.private_name(0))

    variants = {
        module.publication_transaction(
            tmp_path / "other.json", b"managed\n", 0o644, "managed fixture"
        ).transaction_id,
        module.publication_transaction(target, b"different\n", 0o644, "managed fixture").transaction_id,
        module.publication_transaction(target, b"managed\n", 0o444, "managed fixture").transaction_id,
        module.publication_transaction(target, b"managed\n", 0o644, "other fixture").transaction_id,
    }
    monkeypatch.setattr(module, "publication_producer_revision", lambda: "0" * 64)
    variants.add(
        module.publication_transaction(
            target, b"managed\n", 0o644, "managed fixture"
        ).transaction_id
    )
    assert base.transaction_id not in variants
    assert len(variants) == 5


@pytest.mark.parametrize("mode", [0o400, 0o600, 0o666])
def test_direct_final_rejects_unmanaged_final_mode_before_creation(tmp_path, mode):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"

    with pytest.raises(ValueError, match="not a managed publication mode"):
        module.write_exact_or_verify(
            target, b"managed\n", mode=mode, label="managed fixture"
        )
    assert not target.exists()


def test_direct_final_exact_existing_copy_is_verification_only(tmp_path, monkeypatch):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    target.write_bytes(b"managed\n")
    target.chmod(0o644)
    retained = target.stat()

    def forbidden_mutation(*_args, **_kwargs):
        raise AssertionError("exact occupied target verification must not mutate")

    monkeypatch.setattr(module.os, "write", forbidden_mutation)
    monkeypatch.setattr(module.os, "fchmod", forbidden_mutation)
    monkeypatch.setattr(module.os, "rename", forbidden_mutation)
    monkeypatch.setattr(module.os, "replace", forbidden_mutation)
    monkeypatch.setattr(module.os, "link", forbidden_mutation)
    monkeypatch.setattr(module.os, "unlink", forbidden_mutation)
    assert not module.write_exact_or_verify(
        target, b"managed\n", mode=0o644, label="managed fixture"
    )
    assert target.read_bytes() == b"managed\n"
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )


@pytest.mark.parametrize(
    ("occupied_state", "mode", "payload"),
    [
        ("unrelated-owner-only", 0o600, b"unrelated owner-only target\n"),
        ("unrelated-owner-only-hardlink", 0o600, b"unrelated owner-only target\n"),
        ("different-final-bytes", 0o644, b"different\n"),
        ("unbound-final-hardlink", 0o644, b"managed\n"),
    ],
)
def test_occupied_public_target_is_never_mutated(
    tmp_path, monkeypatch, occupied_state, mode, payload,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    target.write_bytes(payload)
    target.chmod(mode)
    extra = parent / "extra-link"
    if "hardlink" in occupied_state:
        os.link(target, extra)
    retained = target.stat()

    def forbidden_mutation(*_args, **_kwargs):
        raise AssertionError("occupied public target must never authorize mutation")

    monkeypatch.setattr(module.os, "write", forbidden_mutation)
    monkeypatch.setattr(module.os, "fchmod", forbidden_mutation)
    monkeypatch.setattr(module.os, "link", forbidden_mutation)
    monkeypatch.setattr(module.os, "unlink", forbidden_mutation)
    with pytest.raises(ValueError):
        module.write_exact_or_verify(
            target, b"managed\n", mode=0o644, label="managed fixture"
        )
    assert target.read_bytes() == payload
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )
    assert stat.S_IMODE(target.stat().st_mode) == mode
    assert target.stat().st_nlink == retained.st_nlink


def test_occupied_transaction_private_name_is_preserved_and_rejected(tmp_path):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    transaction = module.publication_transaction(target, payload, 0o644, "managed fixture")
    collision = parent / transaction.private_name(0)
    collision.write_bytes(b"unrelated private collision\n")
    collision.chmod(0o600)
    retained = collision.stat()

    with pytest.raises(ValueError, match="deterministic public target is absent"):
        module.write_exact_or_verify(
            target, payload, mode=0o644, label="managed fixture"
        )
    assert collision.read_bytes() == b"unrelated private collision\n"
    assert (collision.stat().st_dev, collision.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )
    assert not target.exists()
    assert publication_private_entries(module, parent) == [collision]


def test_linked_recovery_rejects_transaction_named_different_inode_without_unlink(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    target.write_bytes(payload)
    target.chmod(0o644)
    unrelated_link = parent / "unrelated-public-link"
    os.link(target, unrelated_link)
    transaction = module.publication_transaction(target, payload, 0o644, "managed fixture")
    decoy = parent / transaction.private_name(0)
    decoy.write_bytes(payload)
    decoy.chmod(0o644)
    target_profile = target.stat()
    decoy_profile = decoy.stat()

    def forbidden_unlink(*_args, **_kwargs):
        raise AssertionError("different-inode transaction decoy must never be unlinked")

    monkeypatch.setattr(module.os, "unlink", forbidden_unlink)
    with pytest.raises(ValueError, match="deterministic public target is occupied"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    assert (target.stat().st_dev, target.stat().st_ino) == (
        target_profile.st_dev,
        target_profile.st_ino,
    )
    assert (decoy.stat().st_dev, decoy.stat().st_ino) == (
        decoy_profile.st_dev,
        decoy_profile.st_ino,
    )
    assert target.stat().st_nlink == 2
    assert decoy.stat().st_nlink == 1


def test_public_target_race_fails_closed_and_preserves_finalized_private_inode(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    real_rename = module.renameat2_noreplace
    raced = False

    def create_target_before_move(directory, source, destination, label):
        nonlocal raced
        if not raced and destination == target.name:
            target.write_bytes(b"raced replacement\n")
            target.chmod(0o644)
            raced = True
        return real_rename(directory, source, destination, label)

    monkeypatch.setattr(module, "renameat2_noreplace", create_target_before_move)
    with pytest.raises(ValueError, match="deterministic public target is occupied.*no rollback"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    assert raced
    assert target.read_bytes() == b"raced replacement\n"
    retained_private = publication_private_entries(module, parent)
    assert len(retained_private) == 1
    assert retained_private[0].read_bytes() == payload
    assert stat.S_IMODE(retained_private[0].stat().st_mode) == 0o644
    with pytest.raises(ValueError, match="ambiguous transaction-private remnants"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    assert target.read_bytes() == b"raced replacement\n"


def test_post_private_create_exception_preserves_unknown_attempt_and_retry_recovers(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    real_open = module.os.open
    injected = False

    def create_then_raise(name, flags, mode=0o777, *, dir_fd=None):
        nonlocal injected
        if (
            not injected
            and module.PUBLICATION_PRIVATE_NAME_PATTERN.fullmatch(name) is not None
            and flags & os.O_EXCL
        ):
            injected = True
            descriptor = real_open(name, flags, mode, dir_fd=dir_fd)
            os.close(descriptor)
            raise RuntimeError("reported open failure after private creation")
        return real_open(name, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(module.os, "open", create_then_raise)
    retained_umask = os.umask(0o777)
    try:
        with pytest.raises(ValueError, match="deterministic public target is absent"):
            module.write_exact_or_verify(
                target, b"managed\n", mode=0o644, label="managed fixture"
            )
    finally:
        os.umask(retained_umask)
    assert injected
    assert not target.exists()
    retained_private = publication_private_entries(module, parent)
    assert len(retained_private) == 1
    assert retained_private[0].read_bytes() == b""
    assert stat.S_IMODE(retained_private[0].stat().st_mode) == 0o600
    retained = retained_private[0].stat()
    monkeypatch.setattr(module.os, "open", real_open)
    assert module.write_exact_or_verify(
        target, b"managed\n", mode=0o644, label="managed fixture"
    )
    assert target.read_bytes() == b"managed\n"
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert not retained_private[0].exists()
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )


def test_direct_final_partial_write_failure_is_forward_recoverable(tmp_path, monkeypatch):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed payload\n"
    real_write = module.os.write
    injected = False

    def write_partial_then_raise(descriptor, retained):
        nonlocal injected
        if not injected:
            injected = True
            real_write(descriptor, retained[:4])
            raise RuntimeError("reported write failure after partial direct-final write")
        return real_write(descriptor, retained)

    monkeypatch.setattr(module.os, "write", write_partial_then_raise)
    with pytest.raises(ValueError, match="deterministic public target is absent.*no rollback"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    assert injected
    assert not target.exists()
    retained_private = publication_private_entries(module, parent)
    assert len(retained_private) == 1
    assert stat.S_IMODE(retained_private[0].stat().st_mode) == 0o600
    assert retained_private[0].read_bytes() == payload[:4]
    retained = retained_private[0].stat()
    monkeypatch.setattr(module.os, "write", real_write)
    assert module.write_exact_or_verify(
        target, payload, mode=0o644, label="managed fixture"
    )
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert not retained_private[0].exists()
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )


@pytest.mark.parametrize("operation", ["private-fsync", "private-fchmod"])
def test_private_finalize_crash_points_leave_public_absent_and_retry(
    tmp_path, monkeypatch, operation,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed payload\n"
    real_fsync = module.os.fsync
    real_fchmod = module.os.fchmod
    injected = False

    def fsync_then_raise(descriptor):
        nonlocal injected
        result = real_fsync(descriptor)
        if operation == "private-fsync" and not injected and stat.S_ISREG(os.fstat(descriptor).st_mode):
            injected = True
            raise RuntimeError("reported private fsync failure after success")
        return result

    def fchmod_then_raise(descriptor, mode):
        nonlocal injected
        result = real_fchmod(descriptor, mode)
        if operation == "private-fchmod" and not injected:
            injected = True
            raise RuntimeError("reported private fchmod failure after success")
        return result

    monkeypatch.setattr(module.os, "fsync", fsync_then_raise)
    monkeypatch.setattr(module.os, "fchmod", fchmod_then_raise)
    with pytest.raises(ValueError, match="deterministic public target is absent"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    assert injected
    assert not target.exists()
    retained_private = publication_private_entries(module, parent)
    assert len(retained_private) == 1
    retained = retained_private[0].stat()

    monkeypatch.setattr(module.os, "fsync", real_fsync)
    monkeypatch.setattr(module.os, "fchmod", real_fchmod)
    assert module.write_exact_or_verify(
        target, payload, mode=0o644, label="managed fixture"
    )
    assert target.read_bytes() == payload
    assert not retained_private[0].exists()
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )


def test_direct_final_lock_loss_before_commit_is_forward_recoverable(tmp_path, monkeypatch):
    module = load_recost_module()
    root = tmp_path / "root"
    root.mkdir()
    lock_path = root / f".mks24_stage_i_{EPOCH_SLUG}.lock"
    lock_path.write_bytes(b"")
    lock_path.chmod(0o644)
    displaced = root / "displaced.lock"
    target = root / "managed.json"
    payload = b"managed\n"
    real_write = module.os.write
    real_fsync = module.os.fsync
    replaced = False
    fsyncs_after_authority_loss = 0

    def write_then_replace_lock(descriptor, retained):
        nonlocal replaced
        result = real_write(descriptor, retained)
        if not replaced:
            replaced = True
            lock_path.rename(displaced)
            lock_path.write_bytes(b"")
            lock_path.chmod(0o644)
        return result

    def count_fsync_after_authority_loss(descriptor):
        nonlocal fsyncs_after_authority_loss
        if replaced and displaced.exists():
            fsyncs_after_authority_loss += 1
        return real_fsync(descriptor)

    monkeypatch.setattr(module.os, "write", write_then_replace_lock)
    monkeypatch.setattr(module.os, "fsync", count_fsync_after_authority_loss)
    with module.stage_i_lock(root) as mutation_lock:
        try:
            with pytest.raises(
                ValueError, match="authority revalidation failed and no rollback"
            ):
                module.write_exact_or_verify(
                    target,
                    payload,
                    mode=0o644,
                    label="managed fixture",
                    mutation_lock=mutation_lock,
                )
            assert not target.exists()
            retained_private = publication_private_entries(module, root)
            assert len(retained_private) == 1
            retained = retained_private[0].stat()
            assert fsyncs_after_authority_loss == 0
        finally:
            lock_path.unlink()
            displaced.rename(lock_path)
        assert module.write_exact_or_verify(
            target,
            payload,
            mode=0o644,
            label="managed fixture",
            mutation_lock=mutation_lock,
        )
    assert replaced
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    assert not retained_private[0].exists()
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )
    assert fsyncs_after_authority_loss == 0


def test_linked_transaction_inode_is_retained_without_unsafe_unlink(
    tmp_path, monkeypatch
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    transaction = module.publication_transaction(target, payload, 0o644, "managed fixture")
    private = parent / transaction.private_name(0)
    private.write_bytes(payload)
    private.chmod(0o644)
    os.link(private, target)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("legacy linked state must not raw-unlink")

    monkeypatch.setattr(module.os, "unlink", forbidden)
    with pytest.raises(ValueError, match="deterministic public target is occupied.*links 2"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    retained_private = publication_private_entries(module, parent)
    assert len(retained_private) == 1
    assert target.stat().st_nlink == 2
    assert (target.stat().st_dev, target.stat().st_ino) == (
        retained_private[0].stat().st_dev,
        retained_private[0].stat().st_ino,
    )


def test_ambiguous_successful_no_replace_move_is_classified_and_completed(
    tmp_path, monkeypatch
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    real_rename = module.renameat2_noreplace
    injected = False

    def rename_then_raise(directory, source, destination, label):
        nonlocal injected
        result = real_rename(directory, source, destination, label)
        if not injected:
            injected = True
            raise RuntimeError("reported rename failure after success")
        return result

    monkeypatch.setattr(module, "renameat2_noreplace", rename_then_raise)
    assert module.write_exact_or_verify(
        target, b"managed\n", mode=0o644, label="managed fixture"
    )
    assert injected
    assert target.read_bytes() == b"managed\n"
    assert target.stat().st_nlink == 1
    assert publication_private_entries(module, parent) == []


def test_exact_final_rejects_transaction_private_remnants(tmp_path):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    target.write_bytes(payload)
    target.chmod(0o644)
    transaction = module.publication_transaction(target, payload, 0o644, "managed fixture")
    private = parent / transaction.private_name(0)
    private.write_bytes(payload)
    private.chmod(0o644)
    target_snapshot = target.stat()
    private_snapshot = private.stat()

    with pytest.raises(ValueError, match="ambiguous transaction-private remnants"):
        module.write_exact_or_verify(
            target, payload, mode=0o644, label="managed fixture"
        )
    with pytest.raises(ValueError, match="ambiguous transaction-private remnants"):
        module.preflight_exact_or_absent(
            target, payload, mode=0o644, label="managed fixture"
        )
    assert (target.stat().st_dev, target.stat().st_ino) == (
        target_snapshot.st_dev,
        target_snapshot.st_ino,
    )
    assert (private.stat().st_dev, private.stat().st_ino) == (
        private_snapshot.st_dev,
        private_snapshot.st_ino,
    )


def test_staged_output_rejects_unknown_0600_and_accepts_exact_final_commit_marker(tmp_path):
    module = load_recost_module()
    parent = tmp_path / "accounting"
    parent.mkdir()
    output = parent / "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json.staged"
    output.write_bytes(b"unrelated owner-only staged target\n")
    output.chmod(0o600)
    retained = output.stat()

    with pytest.raises(ValueError, match="output namespace is not empty"):
        module.require_output_namespace_empty(output)
    with pytest.raises(ValueError, match="mode is 0600, expected 0444"):
        module.write_staged_output(output, b'{"checkpoint":"F-117"}\n')
    assert output.read_bytes() == b"unrelated owner-only staged target\n"
    assert (output.stat().st_dev, output.stat().st_ino) == (
        retained.st_dev,
        retained.st_ino,
    )

    output.unlink()
    payload = b'{"checkpoint":"F-117"}\n'
    module.write_staged_output(output, payload)
    module.require_output_namespace_empty(output)
    module.write_staged_output(output, payload)
    assert output.read_bytes() == payload
    assert stat.S_IMODE(output.stat().st_mode) == 0o444

def test_direct_final_post_publication_fsync_failure_is_forward_recoverable(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"managed\n"
    real_fsync = module.os.fsync
    injected = False

    def fsync_then_raise_once(descriptor):
        nonlocal injected
        result = real_fsync(descriptor)
        if (
            not injected
            and stat.S_ISDIR(os.fstat(descriptor).st_mode)
            and target.exists()
            and stat.S_IMODE(target.stat().st_mode) == 0o644
        ):
            injected = True
            raise RuntimeError("reported directory fsync failure after publication")
        return result

    monkeypatch.setattr(module.os, "fsync", fsync_then_raise_once)
    assert module.write_exact_or_verify(
        target, payload, mode=0o644, label="managed fixture"
    )
    assert injected
    assert target.read_bytes() == payload
    assert stat.S_IMODE(target.stat().st_mode) == 0o644
    monkeypatch.setattr(module.os, "fsync", real_fsync)
    assert not module.write_exact_or_verify(
        target, payload, mode=0o644, label="managed fixture"
    )


def test_exact_existing_copy_verification_rejects_parent_path_substitution(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent = tmp_path / "parent"
    parent.mkdir()
    target = parent / "managed.json"
    payload = b"exact\n"
    target.write_bytes(payload)
    target.chmod(0o644)
    displaced_parent = tmp_path / "displaced-parent"
    original_hash = module.sha256_descriptor
    swapped = False

    def swap_parent_after_hash(descriptor):
        nonlocal swapped
        digest = original_hash(descriptor)
        if not swapped:
            swapped = True
            parent.rename(displaced_parent)
            parent.mkdir()
            replacement = parent / target.name
            replacement.write_bytes(payload)
            replacement.chmod(0o644)
        return digest

    monkeypatch.setattr(module, "sha256_descriptor", swap_parent_after_hash)
    with pytest.raises(ValueError, match="parent pathname changed during publication"):
        module.write_exact_or_verify(target, payload, mode=0o644, label="managed fixture")
    assert target.read_bytes() == payload
    assert (displaced_parent / target.name).read_bytes() == payload


def test_generator_ignores_hostile_git_path_and_global_system_configuration(
    recost_fixture, tmp_path,
):
    hostile_bin = tmp_path / "bin"
    hostile_bin.mkdir()
    marker = tmp_path / "caller-git-ran"
    hostile_git = hostile_bin / "git"
    hostile_git.write_text(f"#!/bin/sh\n: > {marker}\nexit 99\n")
    hostile_git.chmod(0o755)
    hostile_config = tmp_path / "hostile.gitconfig"
    hostile_config.write_text("[core]\n\trepositoryformatversion = 99\n")
    environment = dict(os.environ)
    environment.update(
        {
            "PATH": str(hostile_bin),
            "GIT_CONFIG_GLOBAL": str(hostile_config),
            "GIT_CONFIG_SYSTEM": str(hostile_config),
            "GIT_CONFIG_NOSYSTEM": "0",
            "GIT_CONFIG_PARAMETERS": "'core.repositoryformatversion=99'",
            "GIT_EXEC_PATH": str(hostile_bin),
        }
    )
    completed = run_generator(recost_fixture, environment=environment)
    assert completed.returncode == 0, completed.stderr
    assert not marker.exists()


def test_git_queries_disable_repository_local_exec_and_worktree_redirect(tmp_path):
    module = load_recost_module()
    repository = tmp_path / "repository"
    repository.mkdir()
    git(repository, "init", "-q")
    tracked = repository / "tracked"
    tracked.write_text("committed\n")
    git(repository, "add", "tracked")
    git(
        repository,
        "-c",
        "user.name=CGL fixture",
        "-c",
        "user.email=cgl-fixture@example.invalid",
        "commit",
        "-q",
        "-m",
        "Create local-config fixture",
    )
    marker = tmp_path / "local-config-executed"
    executable = tmp_path / "hostile-local-config"
    executable.write_text(f"#!/bin/sh\n: > {marker}\nexit 0\n")
    executable.chmod(0o755)
    redirected = tmp_path / "redirected-worktree"
    redirected.mkdir()
    (redirected / "tracked").write_text("committed\n")
    git(repository, "config", "core.fsmonitor", str(executable))
    git(repository, "config", "diff.external", str(executable))
    git(repository, "config", "core.worktree", str(redirected))
    tracked.write_text("changed\n")

    completed = module.git_run(repository, ["diff", "--quiet", "--", "tracked"])
    assert completed.returncode == 1
    assert not marker.exists()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("schema", "lacks controller schema-4 evidence"),
        ("failure-count", "continuation physics evidence differs"),
        ("legacy-checks", "clean-partial checks schema differs"),
        ("legacy-schema", "not eligible for frozen-E03 no-max_ndiv migration"),
        ("plasma-drift", "continuation physics evidence differs"),
    ),
)
def test_clean_partial_continuation_requires_controller_schema4_physics_evidence(
    recost_fixture, mutation, message,
):
    manifest = recost_fixture["manifest"]
    assert isinstance(manifest, Path)
    value = json.loads(manifest.read_text())
    inspection = value["scientific_inspection"]
    if mutation == "schema":
        inspection["schema_version"] = 3
    elif mutation == "failure-count":
        inspection["maximum_strict_failure_counts"]["lf_nonfin"] = 1
    elif mutation == "legacy-checks":
        del inspection["checks"]["plasma_continuation_policy"]
    elif mutation == "legacy-schema":
        del inspection["checks"]["plasma_continuation_policy"]
        del inspection["plasma_continuation_policy"]
        del inspection["plasma_continuation_evidence"]
    else:
        inspection["plasma_continuation_evidence"]["measurements"][
            "normalized_ct_divb_max"
        ] = 0.0
    write_json(manifest, value)
    refresh_request(recost_fixture)
    assert_rejected(run_generator(recost_fixture), message)


@pytest.mark.parametrize("case_id", ("R03", "R12"))
def test_clean_partial_accepts_exact_frozen_e03_no_max_ndiv_live_migration(
    case_id, monkeypatch,
):
    """Use exact retained no-max_ndiv records without mutating canonical state."""

    segments = {
        "R03": "s00_rankio_t0_t0p5",
        "R12": "s00_rankio_t0_t0p25",
    }
    manifest_path = Path(
        "/lustre/orion/ast207/proj-shared/dfielding/CGL/runs/mks24-stage-i/"
        f"E03-forcing-policy/{case_id}/{segments[case_id]}/manifest/prepared_run.json"
    )
    if not manifest_path.is_file():
        pytest.skip(f"retained live-compatible {case_id} path is unavailable")
    module = load_recost_module()
    manifest = json.loads(manifest_path.read_text())
    manifest["_manifest_path"] = str(manifest_path)
    assert module.manifest_identity(manifest) == (
        manifest["job_id"],
        case_id,
        segments[case_id],
        "clean_partial",
    )
    inspection = manifest["scientific_inspection"]
    module.require_clean_partial_controller_evidence(
        inspection,
        manifest,
        manifest["job_id"],
        inspection["final_time"],
        manifest["allocation"]["nodes"] * 8,
    )
    mhd = module.parse_controller_history(
        Path(inspection["mhd_history"]["path"]).read_bytes(), f"live {case_id} MHD history"
    )
    user = module.parse_controller_history(
        Path(inspection["user_history"]["path"]).read_bytes(), f"live {case_id} user history"
    )
    recost_evidence = module.require_frozen_e03_no_max_ndiv_migration(
        inspection,
        manifest,
        manifest["job_id"],
        case_id,
        segments[case_id],
        manifest["allocation"]["nodes"] * 8,
        mhd,
        user,
    )
    controller_manifest = json.loads(manifest_path.read_text())
    controller_evidence = stage_i.revalidate_continuation_plasma_evidence(
        controller_manifest["scientific_inspection"], controller_manifest
    )["plasma_continuation_evidence"]
    assert recost_evidence["authorization_basis"] == (
        "none; historical frozen-E03 evidence is inventory-only and cannot "
        "authorize continuation"
    )
    assert recost_evidence == controller_evidence
    drifted = deepcopy(manifest)
    drifted["command"]["input_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="manifest or inspection bytes differ"):
        module.require_clean_partial_controller_evidence(
            drifted["scientific_inspection"],
            drifted,
            drifted["job_id"],
            drifted["scientific_inspection"]["final_time"],
            drifted["allocation"]["nodes"] * 8,
        )
    contracts = deepcopy(module.FROZEN_E03_NO_MAX_NDIV_MIGRATIONS)
    contracts[(case_id, manifest["job_id"])]["independent_validation"]["sha256"] = "0" * 64
    monkeypatch.setattr(module, "FROZEN_E03_NO_MAX_NDIV_MIGRATIONS", contracts)
    with pytest.raises(ValueError, match="controller migration binding"):
        module.require_frozen_e03_no_max_ndiv_migration(
            inspection,
            manifest,
            manifest["job_id"],
            case_id,
            segments[case_id],
            manifest["allocation"]["nodes"] * 8,
            mhd,
            user,
        )


def test_r12_historical_clean_partial_requires_exact_fresh_four_node_rerun():
    """Keep retained s00 evidence non-authorizing and require the exact fresh s01 run."""

    module = load_recost_module()
    terminal = ("4766856", "R12", "s00_rankio_t0_t0p25", "clean_partial")
    null_parent = (None, None, None, None, None, None)
    assert module.PLAN_INITIAL_TARGETS["R12"] == stage_i.R12_FRESH_RERUN_TARGET == 0.12
    assert module.R12_FRESH_RERUN_SEGMENT == stage_i.R12_FRESH_RERUN_SEGMENT
    assert module.require_authoritative_r12_fresh_rerun(
        "R12",
        terminal,
        stage_i.R12_FRESH_RERUN_SEGMENT,
        1,
        0.0,
        stage_i.R12_FRESH_RERUN_TARGET,
        stage_i.R12_FRESH_RERUN_NODES,
        stage_i.R12_FRESH_RERUN_RANKS_PER_NODE,
        stage_i.R12_FRESH_RERUN_WALLTIME,
        stage_i.R12_FRESH_RERUN_ATHENA_WALLTIME,
        null_parent,
        1,
    )
    assert not module.require_authoritative_r12_fresh_rerun(
        "R03",
        ("4766828", "R03", "s02_rankio_t0p5_t0p75", "accepted"),
        "s03_rankio_t0p75_t1",
        3,
        0.75,
        1.0,
        1,
        8,
        "02:00:00",
        "01:50:00",
        null_parent,
        3,
    )
    for mutation in (
        {"segment": "s01_rankio_t0p25_t0p5"},
        {"nodes": 1},
        {"ranks_per_node": 4},
        {"walltime": "01:59:59"},
        {"athena_walltime": "01:49:59"},
        {"parent_fields": ("4766856", "clean_partial", "s00_rankio_t0_t0p25", None, None, 0.25)},
        {"target": 0.25},
    ):
        values = {
            "segment": "s01_rankio_t0_t0p12",
            "nodes": 4,
            "ranks_per_node": 8,
            "walltime": "02:00:00",
            "athena_walltime": "01:50:00",
            "parent_fields": null_parent,
            "target": 0.12,
        }
        values.update(mutation)
        with pytest.raises(ValueError, match="inventory-only"):
            module.require_authoritative_r12_fresh_rerun(
                "R12",
                terminal,
                values["segment"],
                1,
                0.0,
                values["target"],
                values["nodes"],
                values["ranks_per_node"],
                values["walltime"],
                values["athena_walltime"],
                values["parent_fields"],
                1,
            )


def test_r12_fresh_profile_identity_matches_controller_contract(recost_fixture):
    """Bind the recost fresh recommendation to the controller's exact launch contract."""

    module = load_recost_module()
    fresh = profile(
        recost_fixture,
        case_id="R12",
        segment=module.R12_FRESH_RERUN_SEGMENT,
        nodes=module.R12_FRESH_RERUN_NODES,
        parent=False,
        target=module.R12_FRESH_RERUN_TARGET,
    )
    assert {
        "segment": fresh["segment"],
        "nodes": fresh["nodes"],
        "ranks_per_node": fresh["ranks_per_node"],
        "walltime": fresh["walltime"],
        "athena_walltime": fresh["athena_walltime"],
        "time_tlim_target": fresh["time_tlim_target"],
    } == {
        "segment": stage_i.R12_FRESH_RERUN_SEGMENT,
        "nodes": stage_i.R12_FRESH_RERUN_NODES,
        "ranks_per_node": stage_i.R12_FRESH_RERUN_RANKS_PER_NODE,
        "walltime": stage_i.R12_FRESH_RERUN_WALLTIME,
        "athena_walltime": stage_i.R12_FRESH_RERUN_ATHENA_WALLTIME,
        "time_tlim_target": stage_i.R12_FRESH_RERUN_TARGET,
    }
    assert fresh["nodes"] * fresh["ranks_per_node"] == stage_i.R12_FRESH_RERUN_TOTAL_RANKS
    assert all(
        fresh[field] is None
        for field in (
            "parent_job_id",
            "parent_result",
            "parent_segment",
            "restart_file",
            "restart_file_sha256",
            "restart_time",
        )
    )
    stage_i.require_fresh_r12_profile_contract(
        fresh, None, None, fresh["time_tlim_target"]
    )


def test_continuation_target_alignment_is_r12_specific_and_exact():
    module = load_recost_module()
    for target in (0.12, 0.14, 0.24, 10.0):
        module.require_continuation_target_alignment("R12", target)
    for target in (0.121, 0.25, 0.31):
        with pytest.raises(ValueError, match="exact 0.02 fresh-R12 output cadence"):
            module.require_continuation_target_alignment("R12", target)

    for target in (0.25, 0.5, 10.0):
        module.require_continuation_target_alignment("R03", target)
    for target in (0.24, 0.3, 0.55):
        with pytest.raises(ValueError, match="exact quarter-time output cadence"):
            module.require_continuation_target_alignment("R03", target)


def test_generator_retains_quarter_alignment_and_measured_90_percent_cap(recost_fixture):
    nonquarter = profile(
        recost_fixture,
        target=0.55,
        segment="s01_rankio_t0p25_t0p55",
    )
    set_profiles(recost_fixture, [nonquarter], mode="sole-next-profile")
    assert_rejected(run_generator(recost_fixture), "exact quarter-time output cadence")

    measured_overrun = profile(
        recost_fixture,
        target=0.75,
        segment="s01_rankio_t0p25_t0p75",
    )
    set_profiles(recost_fixture, [measured_overrun], mode="sole-next-profile")
    assert_rejected(
        run_generator(recost_fixture), "interval exceeds measured runtime recommendation"
    )


def test_r12_fresh_rerun_supersedes_historical_inventory_in_authoritative_lineage(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    historical = {
        "_manifest_path": str(tmp_path / "historical.json"),
        "identity": ("4766856", "R12", "s00_rankio_t0_t0p25", "clean_partial"),
        "final_time": 0.1371931229426507,
        "parent": None,
    }
    fresh = {
        "_manifest_path": str(tmp_path / "fresh.json"),
        "identity": ("5000000", "R12", "s01_rankio_t0_t0p12", "clean_partial"),
        "final_time": 0.12,
        "parent": None,
    }
    monkeypatch.setattr(module, "manifest_identity", lambda manifest: manifest["identity"])
    monkeypatch.setattr(module, "final_time", lambda manifest: manifest["final_time"])
    monkeypatch.setattr(module, "manifest_parent_path", lambda manifest: manifest["parent"])
    lineages = module.authenticated_case_lineages(tmp_path, [historical, fresh])
    assert lineages == {"R12": [fresh]}


def test_budget_treats_historical_r12_as_inventory_only_but_preserves_accounting(
    monkeypatch,
):
    module = load_recost_module()
    r03 = {
        "identity": ("5000000", "R03", "s00_rankio_t0_t0p25", "accepted"),
        "final_time": 0.25,
    }
    historical_r12 = {
        "identity": ("4766856", "R12", "s00_rankio_t0_t0p25", "clean_partial"),
        "final_time": 0.1371931229426507,
    }
    monkeypatch.setattr(module, "manifest_identity", lambda manifest: manifest["identity"])
    monkeypatch.setattr(module, "final_time", lambda manifest: manifest["final_time"])

    budget = module.calculate_budget(
        [
            {"job_id": "5000000", "actual_node_hours": "1"},
            {"job_id": "4766856", "actual_node_hours": "100"},
        ],
        {"R03": [r03], "R12": [historical_r12]},
        {
            "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            "R12": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
        },
        [],
        Decimal("0"),
        Decimal("10000"),
        Decimal("10000"),
    )

    assert budget["actual_stage_i_node_hours"] == "101"
    assert budget["measurement_basis"]["job_id"] == "5000000"
    assert (
        budget["case_breakdown"]["R12"]["projection_measurement_basis"]["job_id"]
        == "4766856"
    )
    assert budget["case_breakdown"]["R12"]["authenticated_progress_fraction"] == "0"
    assert budget["case_breakdown"]["R12"]["remaining_simulation_time"] == "10.0"


def test_budget_uses_fresh_r12_for_throughput_and_progress(monkeypatch):
    module = load_recost_module()
    fresh_r12 = {
        "identity": ("5000001", "R12", "s01_rankio_t0_t0p12", "clean_partial"),
        "final_time": 0.12,
    }
    monkeypatch.setattr(module, "manifest_identity", lambda manifest: manifest["identity"])
    monkeypatch.setattr(module, "final_time", lambda manifest: manifest["final_time"])

    budget = module.calculate_budget(
        [
            {"job_id": "4766856", "actual_node_hours": "100"},
            {"job_id": "5000001", "actual_node_hours": "1"},
            {"job_id": "5000002", "actual_node_hours": "1"},
        ],
        {
            "R03": [{
                "identity": ("5000002", "R03", "s00_rankio_t0_t0p25", "accepted"),
                "final_time": 0.25,
            }],
            "R12": [fresh_r12],
        },
        {
            "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            "R12": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
        },
        [],
        Decimal("0"),
        Decimal("1000"),
        Decimal("1000"),
    )

    assert budget["actual_stage_i_node_hours"] == "102"
    assert budget["measurement_basis"]["job_id"] == "5000002"
    assert (
        budget["case_breakdown"]["R12"]["projection_measurement_basis"]["job_id"]
        == "5000001"
    )
    assert budget["case_breakdown"]["R12"]["authenticated_progress_fraction"] == "0.012"
    assert budget["case_breakdown"]["R12"]["remaining_simulation_time"] == "9.8800"


def test_budget_does_not_exclude_different_r12_job_lookalike(monkeypatch):
    module = load_recost_module()
    identity = ("4766857", "R12", "s00_rankio_t0_t0p25", "clean_partial")
    manifest = {"identity": identity, "final_time": 0.25}
    monkeypatch.setattr(module, "manifest_identity", lambda value: value["identity"])
    monkeypatch.setattr(module, "final_time", lambda value: value["final_time"])

    budget = module.calculate_budget(
        [
            {"job_id": "5000000", "actual_node_hours": "0.1"},
            {"job_id": identity[0], "actual_node_hours": "1"},
        ],
        {
            "R03": [{
                "identity": ("5000000", "R03", "s00_rankio_t0_t0p25", "accepted"),
                "final_time": 0.25,
            }],
            "R12": [manifest],
        },
        {
            "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            "R12": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
        },
        [],
        Decimal("0"),
        Decimal("100"),
        Decimal("100"),
    )

    assert budget["measurement_basis"]["job_id"] == "5000000"
    assert (
        budget["case_breakdown"]["R12"]["projection_measurement_basis"]["job_id"]
        == identity[0]
    )
    assert budget["case_breakdown"]["R12"]["authenticated_progress_fraction"] == "0.025"


@pytest.mark.parametrize(
    "identity",
    (
        ("4766856", "R11", "s00_rankio_t0_t0p25", "clean_partial"),
        ("4766856", "R12", "s00_rankio_t0_t0p25", "accepted"),
        ("4766856", "R12", "s00_rankio_t0_t0p5", "clean_partial"),
    ),
)
def test_budget_rejects_partial_historical_r12_identity_collisions(monkeypatch, identity):
    module = load_recost_module()
    monkeypatch.setattr(module, "manifest_identity", lambda value: value["identity"])
    monkeypatch.setattr(module, "final_time", lambda value: value["final_time"])

    with pytest.raises(ValueError, match="historical R12 inventory identity partially collides"):
        module.calculate_budget(
            [
                {"job_id": "5000000", "actual_node_hours": "0.1"},
                {"job_id": "4766856", "actual_node_hours": "1"},
            ],
            {
                "R03": [{
                    "identity": ("5000000", "R03", "s00_rankio_t0_t0p25", "accepted"),
                    "final_time": 0.25,
                }],
                "R12": [{"identity": identity, "final_time": 0.25}],
            },
            {
                "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
                "R12": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            },
            [],
            Decimal("0"),
            Decimal("100"),
            Decimal("100"),
        )


@pytest.mark.parametrize(
    ("lineage_case", "identity"),
    (
        ("R11", ("4766856", "R12", "s00_rankio_t0_t0p25", "clean_partial")),
        ("R12", ("5000001", "R03", "s00_rankio_t0_t0p25", "accepted")),
    ),
)
def test_budget_rejects_lineage_key_identity_mismatch(monkeypatch, lineage_case, identity):
    module = load_recost_module()
    monkeypatch.setattr(module, "manifest_identity", lambda value: value["identity"])
    monkeypatch.setattr(module, "final_time", lambda value: value["final_time"])

    with pytest.raises(ValueError, match="budget lineage identity differs"):
        module.calculate_budget(
            [
                {"job_id": "5000000", "actual_node_hours": "0.1"},
                {"job_id": identity[0], "actual_node_hours": "1"},
            ],
            {
                "R03": [{
                    "identity": ("5000000", "R03", "s00_rankio_t0_t0p25", "accepted"),
                    "final_time": 0.25,
                }],
                lineage_case: [{"identity": identity, "final_time": 0.25}],
            },
            {
                "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
                lineage_case: {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            },
            [],
            Decimal("0"),
            Decimal("100"),
            Decimal("100"),
        )


@pytest.mark.parametrize("result", ("cancelled", "failed", "accepted"))
def test_protected_historical_r12_job_rejects_noninventory_results(result):
    module = load_recost_module()
    with pytest.raises(ValueError, match="historical R12 inventory identity partially collides"):
        module.validated_manifest_identity(
            "4766856", "R12", "s00_rankio_t0_t0p25", result
        )


def test_budget_live_projection_allows_exact_fresh_r12_calibration_above_envelope(
    monkeypatch,
):
    module = load_recost_module()
    standard_cells = 192 * 192 * 384
    rows = [
        {"job_id": "R02-total", "actual_node_hours": "23.996113"},
        {"job_id": "R03-total", "actual_node_hours": "2.942500"},
        {"job_id": "4766847", "actual_node_hours": "1.511111"},
        {"job_id": "4766856", "actual_node_hours": "7.396667"},
        {"job_id": "4766866", "actual_node_hours": "0.730833"},
    ]
    manifests = {
        "R02": [{"identity": ("R02-total", "R02", "s00_rankio_t0_t10", "accepted"),
                 "final_time": 10.0}],
        "R03": [{"identity": ("R03-total", "R03", "s00_rankio_t0_t0p5", "accepted"),
                 "final_time": 0.5}],
        "R04": [{"identity": ("4766847", "R04", "s01_rankio_t0_t0p25", "accepted"),
                 "final_time": 0.25}],
        "R12": [{"identity": module.R12_HISTORICAL_INVENTORY_IDENTITY,
                 "final_time": 0.1371931229426507}],
        "R16": [{"identity": ("4766866", "R16", "s00_rankio_t0_t1p5", "accepted"),
                 "final_time": 1.5}],
    }
    matrix = {
        f"R{case:02d}": {
            "_cell_count": (
                96 * 96 * 192 if case == 16
                else 384 * 384 * 768 if case == 17
                else standard_cells
            ),
            "_estimated_node_hours": Decimal("1"),
        }
        for case in range(2, 18)
    }
    profiles = [
        {"case_id": "R03", "nodes": 1, "walltime": "02:00:00"},
        {"case_id": "R04", "nodes": 4, "walltime": "02:00:00"},
        {
            "case_id": "R12",
            "segment": "s01_rankio_t0_t0p12",
            "nodes": 4,
            "ranks_per_node": 8,
            "walltime": "02:00:00",
            "athena_walltime": "01:50:00",
            "time_tlim_target": 0.12,
        },
        {"case_id": "R16", "nodes": 1, "walltime": "02:00:00"},
    ]
    monkeypatch.setattr(module, "manifest_identity", lambda value: value["identity"])
    monkeypatch.setattr(module, "final_time", lambda value: value["final_time"])

    budget = module.calculate_budget(
        rows,
        manifests,
        matrix,
        profiles,
        Decimal("20"),
        Decimal("1400"),
        Decimal("4000"),
    )

    assert budget["measurement_basis"]["job_id"] == "4766847"
    assert (
        budget["case_breakdown"]["R12"]["projection_measurement_basis"]["job_id"]
        == "4766856"
    )
    assert budget["actual_stage_i_node_hours"] == "36.577224"
    assert budget["case_breakdown"]["R12"]["authenticated_progress_fraction"] == "0"
    assert Decimal(budget["computed_remaining_stage_i_node_hours"]) == Decimal(
        "1749.920383467427316981173717"
    )
    assert Decimal(budget["computed_stage_i_total_node_hours"]) == Decimal(
        "1786.497607467427316981173717"
    )
    assert Decimal(budget["computed_stage_i_margin_node_hours"]) == Decimal(
        "-386.497607467427316981173717"
    )


def test_budget_above_envelope_fails_without_exact_fresh_r12_calibration(monkeypatch):
    module = load_recost_module()
    monkeypatch.setattr(module, "manifest_identity", lambda value: value["identity"])
    monkeypatch.setattr(module, "final_time", lambda value: value["final_time"])

    with pytest.raises(ValueError, match="projection exceeds the promoted envelope"):
        module.calculate_budget(
            [
                {"job_id": "5000000", "actual_node_hours": "1"},
                {"job_id": "4766856", "actual_node_hours": "10"},
            ],
            {
                "R03": [{
                    "identity": ("5000000", "R03", "s00_rankio_t0_t0p25", "accepted"),
                    "final_time": 0.25,
                }],
                "R12": [{
                    "identity": module.R12_HISTORICAL_INVENTORY_IDENTITY,
                    "final_time": 0.1,
                }],
            },
            {
                "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
                "R12": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            },
            [],
            Decimal("0"),
            Decimal("100"),
            Decimal("2000"),
        )


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("parent_job_id", "4766856"),
        ("parent_result", "clean_partial"),
        ("parent_segment", "s00_rankio_t0_t0p25"),
        ("restart_file", "/tmp/restart"),
        ("restart_file_sha256", "a" * 64),
        ("restart_time", 0.1371931229426507),
    ),
)
def test_budget_above_envelope_rejects_nonfresh_r12_parent_fields(monkeypatch, field, value):
    module = load_recost_module()
    monkeypatch.setattr(module, "manifest_identity", lambda item: item["identity"])
    monkeypatch.setattr(module, "final_time", lambda item: item["final_time"])
    profile = {
        "case_id": "R12",
        "segment": "s01_rankio_t0_t0p12",
        "nodes": 4,
        "ranks_per_node": 8,
        "walltime": "02:00:00",
        "athena_walltime": "01:50:00",
        "time_tlim_target": 0.12,
        field: value,
    }

    with pytest.raises(ValueError, match="historical clean partial is inventory-only"):
        module.calculate_budget(
            [
                {"job_id": "5000000", "actual_node_hours": "1"},
                {"job_id": "4766856", "actual_node_hours": "10"},
            ],
            {
                "R03": [{
                    "identity": ("5000000", "R03", "s00_rankio_t0_t0p25", "accepted"),
                    "final_time": 0.25,
                }],
                "R12": [{
                    "identity": module.R12_HISTORICAL_INVENTORY_IDENTITY,
                    "final_time": 0.1,
                }],
            },
            {
                "R03": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
                "R12": {"_cell_count": 1, "_estimated_node_hours": Decimal("10")},
            },
            [profile],
            Decimal("8"),
            Decimal("100"),
            Decimal("2000"),
        )


def test_authoritative_lineage_allows_only_authenticated_non_scientific_index_gaps(
    tmp_path, monkeypatch,
):
    module = load_recost_module()
    parent_path = tmp_path / "s00.json"
    parent = {
        "_manifest_path": str(parent_path),
        "identity": ("4762472", "R03", "s00_rankio_t0_t0p5", "clean_partial"),
        "final_time": 0.312823,
        "parent": None,
        "command": {"executable_sha256": "e", "input_sha256": "i"},
    }
    cancelled = {
        "_manifest_path": str(tmp_path / "s01.json"),
        "identity": ("4766485", "R03", "s01_rankio_t0p312823_t0p5", "cancelled"),
        "final_time": 0.0,
        "parent": parent_path,
    }
    accepted = {
        "_manifest_path": str(tmp_path / "s02.json"),
        "identity": ("4766828", "R03", "s02_rankio_t0p312823_t0p5", "accepted"),
        "final_time": 0.5,
        "parent": parent_path,
        "command": {
            "parent_segment": {
                "execution_epoch": EPOCH,
                "segment": "s00_rankio_t0_t0p5",
                "result": "clean_partial",
                "job_id": "4762472",
                "final_time": 0.312823,
                "restart_time": 0.312823,
                "executable_sha256": "e",
                "input_sha256": "i",
                "restart_sha256": "r",
                "restart_files": ["rank0.rst"],
            }
        },
    }
    monkeypatch.setattr(module, "manifest_identity", lambda manifest: manifest["identity"])
    monkeypatch.setattr(module, "final_time", lambda manifest: manifest["final_time"])
    monkeypatch.setattr(module, "manifest_parent_path", lambda manifest: manifest["parent"])
    monkeypatch.setattr(
        module,
        "terminal_restart_binding",
        lambda root, manifest, label: {
            "sha256": "r",
            "rank_files": [{"path": "rank0.rst"}],
        },
    )
    assert module.authenticated_case_lineages(tmp_path, [parent, cancelled, accepted]) == {
        "R03": [parent, accepted]
    }
    with pytest.raises(ValueError, match="parent identity differs"):
        module.authenticated_case_lineages(tmp_path, [parent, accepted])


def test_locked_reconcile_uses_authenticated_in_process_report(tmp_path):
    module = load_recost_module()
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
        "root": str(root),
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
    with lock.open("r+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        report = module.run_authenticated_reconcile(
            helper, sha256(helper), root, stage_i_lock_held=True
        )
    assert report["consistent"] is True
    assert report["root"] == str(root)


def write_failed_f117_attempt(module, root: Path) -> dict[str, Path]:
    """Write one internally bound failed F117 draft state and pin its digests."""

    accounting = root / "accounting"
    accounting.mkdir(parents=True, exist_ok=True)
    paths = {
        key: root / relative
        for key, relative in module.F117_FAILED_ATTEMPT_RELATIVES.items()
    }
    generated = (datetime.now(timezone.utc) - timedelta(minutes=5)).isoformat()
    reconciliation = {
        "execution_epoch": EPOCH,
        "root": str(root),
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
    storage = {
        "schema_version": 1,
        "record_type": "stage-i-storage-evidence",
        "execution_epoch": EPOCH,
        "root": str(root),
        "measured_utc": generated,
        "available_bytes": 3,
        "retained_stage_i_bytes": 0,
        "required_safety_bytes": 1,
        "projected_authorized_wave_growth_bytes": 2,
        "projection_method": STORAGE_PROJECTION_METHOD,
        "profile_projections_sha256": "a" * 64,
    }
    write_json(paths["reconciliation"], reconciliation)
    write_json(paths["storage"], storage)
    legacy = {
        "path": f"accounting/{module.LEGACY_F114_RECOST_NAME}",
        "sha256": module.LEGACY_F114_RECOST_SHA256,
    }
    legacy_audit = {
        "path": f"accounting/{module.LEGACY_F114_RECOST_NAME}.publication_audit.json",
        "sha256": module.LEGACY_F114_PUBLICATION_AUDIT_SHA256,
    }
    inputs = {
        "reconciliation": None,
        "ledger": None,
        "reservations": None,
        "manifests": None,
        "scheduler_evidence": None,
        "storage_evidence": None,
        "source_bundle": {"fixture": True},
        "matrix": {"fixture": True},
        "stage_i_helper": {"fixture": True},
        "ceiling_evidence": {"fixture": True},
        "ceiling_publication_audit": {"fixture": True},
        "source_authority": {"checkpoint": "F-116"},
        "qualification_approval": {"fixture": True},
        "predecessor_recost": legacy,
        "predecessor_recost_independent_review": None,
        "predecessor_recost_publication_audit": legacy_audit,
        "r17_readiness_evidence": None,
        "r17_readiness_independent_review": None,
        "r17_readiness_publication_audit": None,
    }
    packet = {
        "schema_version": 1,
        "record_type": "stage-i-recost-request-draft-packet",
        "checkpoint": "F-117",
        "artifact_name": "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json",
        "execution_epoch": EPOCH,
        "generated_utc": generated,
        "expires_utc": (datetime.now(timezone.utc) + timedelta(hours=1)).isoformat(),
        "requested_by": "fixture failed F117 author",
        "scope": "Exact failed F117 draft state.",
        "barrier": {"recorded_segments": []},
        "inputs": inputs,
        "recommendations": {"mode": "bounded-wave", "max_wave_nodes": 1, "profiles": []},
        "draft_policy": {
            "independent_review_created": False,
            "self_approved": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
            "required_storage_safety_bytes": 1,
        },
    }
    write_json(paths["packet"], packet)
    request = deepcopy(packet)
    request.pop("draft_policy")
    request["schema_version"] = 2
    request["record_type"] = "stage-i-recost-recommendation-request"
    request["inputs"].update(
        {
            "reconciliation": {
                "path": module.F117_FAILED_ATTEMPT_RELATIVES["reconciliation"].as_posix(),
                "sha256": sha256(paths["reconciliation"]),
            },
            "ledger": {"fixture": True},
            "reservations": {"fixture": True},
            "manifests": [],
            "scheduler_evidence": [],
            "storage_evidence": {
                "path": module.F117_FAILED_ATTEMPT_RELATIVES["storage"].as_posix(),
                "sha256": sha256(paths["storage"]),
            },
        }
    )
    write_json(paths["request"], request)
    module.F117_FAILED_ATTEMPT_SHA256 = {
        key: sha256(path) for key, path in paths.items()
    }
    return paths


def legacy_f114_arguments(module, root: Path, monkeypatch):
    """Return one exact local legacy F114 predecessor argument tuple."""

    accounting = root / "accounting"
    accounting.mkdir(parents=True, exist_ok=True)
    artifact = accounting / module.LEGACY_F114_RECOST_NAME
    generated = datetime.now(timezone.utc) - timedelta(minutes=10)
    published = generated + timedelta(minutes=1)
    current = generated + timedelta(minutes=5)
    recost = {
        "schema_version": 1,
        "record_type": "stage-i-clean-partial-recost-checkpoint",
        "execution_epoch": EPOCH,
        "checkpoint": "F-114",
        "generated_utc": generated.isoformat(),
    }
    write_json(artifact, recost)
    audit = artifact.with_name(f"{artifact.name}.publication_audit.json")
    audit_value = {
        "schema_version": 1,
        "record_type": "observed-publication",
        "execution_epoch": EPOCH,
        "published_utc": published.isoformat(),
        "artifact": {
            "path": str(artifact),
            "sha256": sha256(artifact),
            "mode": "0644",
            "links": 1,
        },
    }
    write_json(audit, audit_value)
    monkeypatch.setattr(module, "LEGACY_F114_RECOST_SHA256", sha256(artifact))
    monkeypatch.setattr(module, "LEGACY_F114_PUBLICATION_AUDIT_SHA256", sha256(audit))
    return artifact, recost, audit, audit_value, current


def test_legacy_f114_bootstrap_is_exact_and_consumable_once(monkeypatch, tmp_path):
    module = load_recost_module()
    root = tmp_path / "root"
    accounting = root / "accounting"
    accounting.mkdir(parents=True)
    artifact = accounting / module.LEGACY_F114_RECOST_NAME
    generated = datetime.now(timezone.utc) - timedelta(minutes=10)
    published = generated + timedelta(minutes=1)
    current = generated + timedelta(minutes=5)
    recost = {
        "schema_version": 1,
        "record_type": "stage-i-clean-partial-recost-checkpoint",
        "execution_epoch": EPOCH,
        "checkpoint": "F-114",
        "generated_utc": generated.isoformat(),
    }
    write_json(artifact, recost)
    audit = artifact.with_name(f"{artifact.name}.publication_audit.json")
    audit_value = {
        "schema_version": 1,
        "record_type": "observed-publication",
        "execution_epoch": EPOCH,
        "published_utc": published.isoformat(),
        "artifact": {
            "path": str(artifact),
            "sha256": sha256(artifact),
            "mode": "0644",
            "links": 1,
        },
    }
    write_json(audit, audit_value)
    monkeypatch.setattr(module, "LEGACY_F114_RECOST_SHA256", sha256(artifact))
    monkeypatch.setattr(module, "LEGACY_F114_PUBLICATION_AUDIT_SHA256", sha256(audit))
    parsed = module.parse_predecessor_recost(
        root,
        artifact,
        sha256(artifact),
        recost,
        None,
        None,
        None,
        audit,
        sha256(audit),
        audit_value,
        current,
        "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json",
        "F-117",
        module.InputTracker(),
    )
    assert parsed["bootstrap"] == "exact-retained-legacy-F114-once"
    assert parsed["independent_review_sha256"] is None
    with pytest.raises(ValueError, match="reserved exactly for first schema-2 F-117"):
        module.parse_predecessor_recost(
            root,
            artifact,
            sha256(artifact),
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current,
            "mks24_stage_i_E03_forcing_policy_F118_recost_evidence.json",
            "F-118",
            module.InputTracker(),
        )
    with pytest.raises(ValueError, match="independent-review path differs"):
        module.parse_predecessor_recost(
            root,
            artifact,
            "0" * 64,
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current,
            "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json",
            "F-117",
            module.InputTracker(),
        )

    schema2 = accounting / "mks24_stage_i_E03_forcing_policy_F199_recost_evidence.json"
    write_immutable_json(schema2, {"schema_version": 2})
    schema2_audit = schema2.with_name(f"{schema2.name}.publication_audit.json")
    write_immutable_json(
        schema2_audit,
        {
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": (published + timedelta(minutes=1)).isoformat(),
            "artifact": {
                "path": str(schema2),
                "sha256": sha256(schema2),
                "mode": "0444",
                "links": 1,
            },
        },
    )
    with pytest.raises(ValueError, match="already consumed"):
        module.parse_predecessor_recost(
            root,
            artifact,
            sha256(artifact),
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current,
            "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json",
            "F-117",
            module.InputTracker(),
        )


def test_f119_exactly_supersedes_authenticated_failed_f117_attempt(monkeypatch, tmp_path):
    module = load_recost_module()
    root = tmp_path / "root"
    artifact, recost, audit, audit_value, current = legacy_f114_arguments(
        module, root, monkeypatch
    )
    paths = write_failed_f117_attempt(module, root)

    parsed = module.parse_predecessor_recost(
        root,
        artifact,
        sha256(artifact),
        recost,
        None,
        None,
        None,
        audit,
        sha256(audit),
        audit_value,
        current,
        "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json",
        "F-119",
        module.InputTracker(),
    )

    assert parsed["bootstrap"] == (
        "exact-retained-legacy-F114-after-authenticated-F117-failed-attempt"
    )
    assert parsed["superseded_failed_attempt"] == {
        "checkpoint": "F-117",
        **{
            key: {
                "path": module.F117_FAILED_ATTEMPT_RELATIVES[key].as_posix(),
                "sha256": sha256(path),
            }
            for key, path in paths.items()
        },
        "status": "authenticated-unpromoted-failed-attempt",
    }


@pytest.mark.parametrize("operation", ("missing", "mutated"))
@pytest.mark.parametrize("failed_key", ("packet", "request", "reconciliation", "storage"))
def test_f119_legacy_supersession_rejects_failed_f117_drift(
    monkeypatch, tmp_path, failed_key, operation
):
    module = load_recost_module()
    root = tmp_path / "root"
    artifact, recost, audit, audit_value, current = legacy_f114_arguments(
        module, root, monkeypatch
    )
    paths = write_failed_f117_attempt(module, root)
    if operation == "missing":
        paths[failed_key].unlink()
    else:
        paths[failed_key].write_bytes(paths[failed_key].read_bytes() + b"\n")
        paths[failed_key].chmod(0o644)

    with pytest.raises((OSError, ValueError)):
        module.parse_predecessor_recost(
            root,
            artifact,
            sha256(artifact),
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current,
            "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json",
            "F-119",
            module.InputTracker(),
        )


@pytest.mark.parametrize("relative_index", range(5))
def test_f119_legacy_supersession_rejects_any_f117_promotion_namespace(
    monkeypatch, tmp_path, relative_index
):
    module = load_recost_module()
    root = tmp_path / "root"
    artifact, recost, audit, audit_value, current = legacy_f114_arguments(
        module, root, monkeypatch
    )
    write_failed_f117_attempt(module, root)
    promoted = root / module.F117_FORBIDDEN_PROMOTION_RELATIVES[relative_index]
    write_json(promoted, {"forbidden": True})

    with pytest.raises(ValueError, match="requires no F117 artifact"):
        module.parse_predecessor_recost(
            root,
            artifact,
            sha256(artifact),
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current,
            "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json",
            "F-119",
            module.InputTracker(),
        )


@pytest.mark.parametrize("checkpoint_number", (116, 118, 120, 200))
def test_failed_f117_legacy_supersession_is_f119_only(
    monkeypatch, tmp_path, checkpoint_number
):
    module = load_recost_module()
    root = tmp_path / "root"
    artifact, recost, audit, audit_value, current = legacy_f114_arguments(
        module, root, monkeypatch
    )
    write_failed_f117_attempt(module, root)

    with pytest.raises(ValueError, match="authenticated failed-F117 supersession by F-119"):
        module.parse_predecessor_recost(
            root,
            artifact,
            sha256(artifact),
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current,
            f"mks24_stage_i_E03_forcing_policy_F{checkpoint_number}_recost_evidence.json",
            f"F-{checkpoint_number}",
            module.InputTracker(),
        )


def test_f119_legacy_supersession_rejects_any_schema2_predecessor(
    monkeypatch, tmp_path
):
    module = load_recost_module()
    root = tmp_path / "root"
    artifact, recost, audit, audit_value, current = legacy_f114_arguments(
        module, root, monkeypatch
    )
    write_failed_f117_attempt(module, root)
    schema2 = root / "accounting/mks24_stage_i_E03_forcing_policy_F118_recost_evidence.json"
    write_immutable_json(schema2, {"schema_version": 2})
    schema2_audit = schema2.with_name(f"{schema2.name}.publication_audit.json")
    write_immutable_json(
        schema2_audit,
        {
            "record_type": "stage-i-recost-recommendation-publication-audit",
            "execution_epoch": EPOCH,
            "published_utc": current.isoformat(),
            "artifact": {
                "path": str(schema2),
                "sha256": sha256(schema2),
                "mode": "0444",
                "links": 1,
            },
        },
    )

    with pytest.raises(ValueError, match="already consumed by schema-2 publication"):
        module.parse_predecessor_recost(
            root,
            artifact,
            sha256(artifact),
            recost,
            None,
            None,
            None,
            audit,
            sha256(audit),
            audit_value,
            current + timedelta(minutes=1),
            "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json",
            "F-119",
            module.InputTracker(),
        )
