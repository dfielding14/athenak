#!/usr/bin/env python3
"""Prepare and account Frontier debug-QOS CGL-LF validation jobs.

This utility intentionally never submits a Slurm job.  It writes an
inspectable debug-QOS batch script, reserves conservative node hours, checks
that the debug queue is empty immediately before manual submission, and
records completed allocation use from ``sacct``.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


ROOT_DIR = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
ACCOUNT = "AST207"
QOS = "debug"
PARTITION = "batch"
BUDGET_NODE_HOURS = 1000.0
MAX_DEBUG_SECONDS = 2 * 60 * 60
LEDGER_COLUMNS = (
    "job_id",
    "submitted_utc",
    "completed_utc",
    "campaign",
    "run_name",
    "purpose",
    "state",
    "exit_code",
    "nodes",
    "requested_walltime",
    "elapsed_seconds",
    "reserved_node_hours",
    "actual_node_hours",
    "cumulative_actual_node_hours",
    "git_revision",
    "executable_sha256",
    "input_file",
    "output_dir",
    "result",
    "notes",
)
NONTERMINAL_STATES = {
    "PENDING",
    "RUNNING",
    "CONFIGURING",
    "COMPLETING",
    "SUSPENDED",
}


def utc_now() -> str:
    """Return a UTC ISO timestamp without subsecond variability."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_json(path: Path, value: object) -> None:
    """Write deterministic JSON used as retained campaign metadata."""

    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_root(root: Path, allow_local_root: bool) -> Path:
    """Require the prescribed project root except during an explicit local test."""

    resolved = root.expanduser().resolve()
    if resolved != DEFAULT_ROOT and not allow_local_root:
        raise ValueError(
            f"Frontier run root must be {DEFAULT_ROOT}; "
            "use --allow-local-root only for offline validation"
        )
    return resolved


def require_beneath_root(path: Path, root: Path, label: str,
                         allow_local_root: bool) -> None:
    """Keep actual Frontier sources and inputs beneath the project run root."""

    if allow_local_root:
        return
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(
            f"{label} must be beneath the Frontier CGL root {root}: {path}"
        ) from error


def parse_walltime(value: str) -> int:
    """Parse an HH:MM:SS walltime and return seconds."""

    match = re.fullmatch(r"(\d{2}):([0-5]\d):([0-5]\d)", value)
    if match is None:
        raise ValueError(f"walltime must use HH:MM:SS: {value}")
    hours, minutes, seconds = (int(item) for item in match.groups())
    return hours * 3600 + minutes * 60 + seconds


def node_hours(nodes: int, seconds: int) -> float:
    """Compute allocated node hours."""

    return nodes * seconds / 3600.0


def layout(root: Path) -> dict[str, Path]:
    """Return prescribed project storage locations."""

    accounting = root / "accounting"
    return {
        "root": root,
        "accounting": accounting,
        "ledger": accounting / "frontier_node_hours.csv",
        "reservations": accounting / "frontier_reservations.json",
        "summary": accounting / "frontier_budget_summary.md",
        "runs": root / "runs",
        "logs_build": root / "logs" / "build",
        "logs_slurm": root / "logs" / "slurm",
        "build": root / "build",
        "inputs": root / "inputs" / "archived",
        "scripts_build": root / "scripts" / "build",
        "scripts_submit": root / "scripts" / "submit",
        "scripts_accounting": root / "scripts" / "accounting",
    }


def initialize(root: Path) -> dict[str, Path]:
    """Initialize retained layout and accounting metadata if absent."""

    paths = layout(root)
    for name, path in paths.items():
        if name not in {"root", "ledger", "reservations", "summary"}:
            path.mkdir(parents=True, exist_ok=True)
    paths["accounting"].mkdir(parents=True, exist_ok=True)
    if not paths["ledger"].exists():
        with paths["ledger"].open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerow(LEDGER_COLUMNS)
    if not paths["reservations"].exists():
        write_json(paths["reservations"], [])
    refresh_summary(paths)
    return paths


def read_ledger(paths: dict[str, Path]) -> list[dict[str, str]]:
    """Read recorded top-level allocations."""

    with paths["ledger"].open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_reservations(paths: dict[str, Path]) -> list[dict[str, object]]:
    """Read prepared and completed reservation records."""

    data = json.loads(paths["reservations"].read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"reservation store is not a list: {paths['reservations']}")
    return data


def active_reservations(reservations: list[dict[str, object]]
                        ) -> list[dict[str, object]]:
    """Return reservations that conservatively consume unrecorded budget."""

    return [
        item for item in reservations
        if item.get("state") in {"prepared", "submitted"}
    ]


def refresh_summary(paths: dict[str, Path]) -> None:
    """Regenerate human-readable project budget accounting."""

    ledger = read_ledger(paths)
    reservations = read_reservations(paths)
    consumed = sum(float(row["actual_node_hours"]) for row in ledger)
    active = active_reservations(reservations)
    reserved = sum(float(item["reserved_node_hours"]) for item in active)
    remaining = BUDGET_NODE_HOURS - consumed - reserved
    lines = [
        "# Frontier CGL-LF Testing Budget",
        "",
        f"- Updated UTC: `{utc_now()}`",
        f"- Authorized test budget: `{BUDGET_NODE_HOURS:.6f}` node-hours",
        f"- Recorded actual use: `{consumed:.6f}` node-hours",
        f"- Active conservative reservations: `{reserved:.6f}` node-hours",
        f"- Uncommitted remaining budget: `{remaining:.6f}` node-hours",
        "",
        "## Recorded Allocations",
        "",
    ]
    if ledger:
        lines.extend([
            "| Job | Run | State | Actual node-hours | Cumulative | Result |",
            "| --- | --- | --- | ---: | ---: | --- |",
        ])
        for row in ledger:
            lines.append(
                "| `{job_id}` | `{campaign}/{run_name}` | {state} | "
                "`{actual_node_hours}` | `{cumulative_actual_node_hours}` | "
                "{result} |".format(**row)
            )
    else:
        lines.append("No Frontier allocations have been recorded.")
    lines.extend(["", "## Active Reservations", ""])
    if active:
        lines.extend([
            "| Run | Purpose | Nodes | Walltime | Reserved node-hours | State |",
            "| --- | --- | ---: | --- | ---: | --- |",
        ])
        for item in active:
            lines.append(
                "| `{campaign}/{run_name}` | {purpose} | `{nodes}` | "
                "`{requested_walltime}` | `{reserved_node_hours:.6f}` | "
                "{state} |".format(**item)
            )
    else:
        lines.append("No prepared or submitted debug validation job is reserved.")
    lines.extend([
        "",
        "Only manually inspected, sequential `debug`-QOS validation jobs may "
        "use this testing ledger. Paper-production execution is excluded.",
        "",
    ])
    paths["summary"].write_text("\n".join(lines), encoding="utf-8")


def git_revision(source_dir: Path, allow_dirty: bool,
                 test_revision: str | None) -> str:
    """Read an immutable source revision for a proposed Frontier executable."""

    if test_revision is not None:
        if not allow_dirty:
            raise ValueError("--test-git-revision is restricted to offline tests")
        return test_revision
    revision = subprocess.run(
        ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "-C", str(source_dir), "status", "--short"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if status and not allow_dirty:
        raise ValueError(
            "refusing Frontier preparation from a dirty source tree; "
            "commit and rebuild an immutable executable first"
        )
    return revision


def validate_debug_input(path: Path) -> None:
    """Reject paper-production input decks from the debug-only utility."""

    lower_name = path.name.lower()
    text = path.read_text(encoding="utf-8").lower()
    if "paper_standard" in lower_name or "paper_nulim" in lower_name:
        raise ValueError(
            "paper-standard and limiter-production inputs may not be "
            "prepared for Frontier debug execution"
        )
    if re.search(r"(?m)^\s*paper_grade\s*=\s*(true|1)\b", text):
        raise ValueError(
            "paper-grade production inputs may not be prepared for "
            "Frontier debug execution"
        )


def quote(value: str | Path) -> str:
    """Quote a shell literal in a generated submission script."""

    return shlex.quote(str(value))


def generated_batch_script(manifest: dict[str, object], script_path: Path) -> str:
    """Generate one explicit debug-only Slurm script for manual submission."""

    run = manifest["run"]
    allocation = manifest["allocation"]
    paths = manifest["paths"]
    command = manifest["command"]
    execution_target = str(command.get("execution_target", "gpu"))
    ranks_per_node = int(allocation["ranks_per_node"])
    cpus_per_task = int(allocation["cpus_per_task"])
    nodes = int(allocation["nodes"])
    overrides = [quote(value) for value in command["overrides"]]
    override_text = " ".join(overrides)
    if override_text:
        override_text = " " + override_text
    restart_file = command.get("restart_file")
    restart_value = quote(str(restart_file)) if restart_file else "''"
    job_name = re.sub(r"[^A-Za-z0-9_]+", "_", str(run["run_name"]))[:60]
    if execution_target == "gpu":
        gpu_request = "#SBATCH --gpus-per-node=8\n"
        accelerator_modules = (
            "module load craype-accel-amd-gfx90a\n"
            "module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2\n"
        )
        target_environment = """export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export HSA_XNACK=1
"""
        target_launch = "  --gpus-per-task=1 --gpu-bind=closest \\\n"
    elif execution_target == "cpu":
        gpu_request = ""
        accelerator_modules = "module load cpe/25.09 cray-mpich/9.0.1\n"
        target_environment = """export MPICH_GPU_SUPPORT_ENABLED=0
unset MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED MPICH_OFI_NIC_POLICY
unset MPICH_GPU_IPC_CACHE_MAX_SIZE HSA_XNACK
"""
        target_launch = ""
    else:
        raise ValueError(f"unsupported execution target: {execution_target}")
    return f"""#!/bin/bash
#SBATCH -J cgl_{job_name}
#SBATCH -A {ACCOUNT}
#SBATCH -o {paths["slurm_log"]}
#SBATCH -p {PARTITION}
#SBATCH -q {QOS}
#SBATCH -t {allocation["requested_walltime"]}
#SBATCH -N {nodes}
{gpu_request}#SBATCH --threads-per-core=1

set -euo pipefail

CGL_ROOT={quote(manifest["project_root"])}
RUN_MANIFEST={quote(script_path.parent / "prepared_run.json")}
ATHENA={quote(command["executable"])}
INPUT={quote(command["input_file"])}
RESTART={restart_value}
OUT_DIR={quote(paths["output_dir"])}
RESTART_DIR={quote(paths["restart_dir"])}
ENV_LOG={quote(paths["environment_log"])}
RANKS_PER_NODE={ranks_per_node}
CPUS_PER_TASK={cpus_per_task}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * RANKS_PER_NODE))"

test -f "${{RUN_MANIFEST}}"
test -x "${{ATHENA}}"
test -f "${{INPUT}}"
if [[ -n "${{RESTART}}" ]]; then
  test -f "${{RESTART}}"
fi
mkdir -p "${{OUT_DIR}}" "${{RESTART_DIR}}"

module restore
module load PrgEnv-cray
{accelerator_modules}module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
{target_environment}export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export OMP_NUM_THREADS=1

{{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "slurm_job_id=${{SLURM_JOB_ID:?Missing SLURM_JOB_ID}}"
  echo "prepared_manifest=${{RUN_MANIFEST}}"
  echo "nodes=${{NNODES}}"
  echo "ranks=${{NRANKS}}"
  module -t list 2>&1
  env | LC_ALL=C sort | grep -E '^(MPICH_|FI_|HSA_|OMP_|ROCR_|HIP_|CRAY_)' || true
}} > "${{ENV_LOG}}"

RUN_ARGS=(-i "${{INPUT}}")
if [[ -n "${{RESTART}}" ]]; then
  RUN_ARGS=(-r "${{RESTART}}")
fi

srun -N "${{NNODES}}" -n "${{NRANKS}}" --ntasks-per-node="${{RANKS_PER_NODE}}" \
  -c "${{CPUS_PER_TASK}}" --threads-per-core=1 --cpu-bind=threads \
{target_launch}  "${{ATHENA}}" "${{RUN_ARGS[@]}}" -d "${{OUT_DIR}}" \
  -t {quote(command["athena_walltime"])} \
  job/basename={quote(run["run_name"])}{override_text}

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "${{ENV_LOG}}"
"""


def reservation_usage(paths: dict[str, Path]) -> tuple[float, float]:
    """Return recorded use and currently reserved conservative use."""

    consumed = sum(
        float(row["actual_node_hours"]) for row in read_ledger(paths)
    )
    reserved = sum(
        float(item["reserved_node_hours"])
        for item in active_reservations(read_reservations(paths))
    )
    return consumed, reserved


def prepare(args: argparse.Namespace) -> Path:
    """Create a retained run manifest and a manual debug submission script."""

    root = require_root(Path(args.root), args.allow_local_root)
    paths = initialize(root)
    if args.allow_dirty_source and not args.allow_local_root:
        raise ValueError("--allow-dirty-source is restricted to offline validation")
    executable = Path(args.executable).expanduser().resolve()
    input_path = Path(args.input_file).expanduser().resolve()
    source_dir = Path(args.source_dir).expanduser().resolve()
    restart_path = (
        Path(args.restart_file).expanduser().resolve()
        if args.restart_file is not None else None
    )
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"executable is missing or not executable: {executable}")
    if not input_path.is_file():
        raise ValueError(f"input deck is missing: {input_path}")
    if restart_path is not None and not restart_path.is_file():
        raise ValueError(f"restart file is missing: {restart_path}")
    require_beneath_root(executable, root, "executable", args.allow_local_root)
    require_beneath_root(input_path, root, "input deck", args.allow_local_root)
    require_beneath_root(source_dir, root, "source directory", args.allow_local_root)
    if restart_path is not None:
        require_beneath_root(
            restart_path, root, "restart file", args.allow_local_root
        )
    validate_debug_input(input_path)
    nodes = args.nodes
    if nodes < 1:
        raise ValueError("--nodes must be positive")
    requested_seconds = parse_walltime(args.walltime)
    if requested_seconds > MAX_DEBUG_SECONDS:
        raise ValueError("Frontier debug walltime must not exceed 02:00:00")
    athena_seconds = parse_walltime(args.athena_walltime)
    margin = 600 if requested_seconds > 3600 else 300
    if athena_seconds > requested_seconds - margin:
        raise ValueError(
            "Athena walltime must leave at least five minutes before a "
            "one-hour-or-shorter job, or ten minutes for a longer job"
        )
    revision = git_revision(
        source_dir, args.allow_dirty_source,
        args.test_git_revision,
    )
    reserved_node_hours = node_hours(nodes, requested_seconds)
    consumed, already_reserved = reservation_usage(paths)
    if consumed + already_reserved + reserved_node_hours > BUDGET_NODE_HOURS:
        raise ValueError("proposed reservation exceeds the 1000 node-hour budget")
    run_dir = paths["runs"] / args.campaign / args.run_name
    manifest_dir = run_dir / "manifest"
    manifest_path = manifest_dir / "prepared_run.json"
    if manifest_path.exists():
        raise ValueError(f"run has already been prepared: {manifest_path}")
    output_dir = run_dir / "output"
    restart_dir = run_dir / "restart"
    analysis_dir = run_dir / "analysis"
    for path in (manifest_dir, output_dir, restart_dir, analysis_dir):
        path.mkdir(parents=True, exist_ok=True)
    archived_input = manifest_dir / "submitted_input.athinput"
    shutil.copy2(input_path, archived_input)
    archived_restart = None
    if restart_path is not None:
        archived_restart = manifest_dir / "submitted_restart.rst"
        shutil.copy2(restart_path, archived_restart)
    batch_script = manifest_dir / "cgl_lf_debug.sbatch"
    manifest: dict[str, object] = {
        "schema_version": 1,
        "state": "prepared",
        "prepared_utc": utc_now(),
        "project_root": str(root),
        "policy": {
            "account": ACCOUNT,
            "partition": PARTITION,
            "qos": QOS,
            "production_execution_allowed": False,
            "authorized_testing_budget_node_hours": BUDGET_NODE_HOURS,
            "sequential_manual_submission_required": True,
        },
        "run": {
            "campaign": args.campaign,
            "run_name": args.run_name,
            "purpose": args.purpose,
            "acceptance_criterion": args.acceptance_criterion,
        },
        "allocation": {
            "nodes": nodes,
            "requested_walltime": args.walltime,
            "requested_seconds": requested_seconds,
            "reserved_node_hours": reserved_node_hours,
            "ranks_per_node": args.ranks_per_node,
            "cpus_per_task": args.cpus_per_task,
        },
        "command": {
            "source_dir": str(source_dir),
            "git_revision": revision,
            "executable": str(executable),
            "executable_sha256": sha256(executable),
            "source_input_file": str(input_path),
            "input_file": str(archived_input),
            "input_sha256": sha256(archived_input),
            "source_restart_file": (
                str(restart_path) if restart_path is not None else None
            ),
            "restart_file": (
                str(archived_restart) if archived_restart is not None else None
            ),
            "restart_sha256": (
                sha256(archived_restart) if archived_restart is not None else None
            ),
            "execution_target": args.execution_target,
            "athena_walltime": args.athena_walltime,
            "overrides": args.override,
        },
        "paths": {
            "run_dir": str(run_dir),
            "output_dir": str(output_dir),
            "restart_dir": str(restart_dir),
            "analysis_dir": str(analysis_dir),
            "environment_log": str(manifest_dir / "run_environment.txt"),
            "batch_script": str(batch_script),
            "slurm_log": str(paths["logs_slurm"] / "%x.%j.log"),
        },
    }
    batch_script.write_text(
        generated_batch_script(manifest, batch_script), encoding="utf-8"
    )
    batch_script.chmod(0o750)
    write_json(manifest_path, manifest)
    reservations = read_reservations(paths)
    reservations.append({
        "manifest": str(manifest_path),
        "campaign": args.campaign,
        "run_name": args.run_name,
        "purpose": args.purpose,
        "nodes": nodes,
        "requested_walltime": args.walltime,
        "reserved_node_hours": reserved_node_hours,
        "state": "prepared",
        "prepared_utc": manifest["prepared_utc"],
    })
    write_json(paths["reservations"], reservations)
    refresh_summary(paths)
    print(f"Prepared Frontier debug validation bundle: {manifest_path}")
    print(f"Reserved node-hours: {reserved_node_hours:.6f}")
    print("Run check-submit immediately before manually issuing sbatch.")
    return manifest_path


def read_manifest(path: Path) -> dict[str, object]:
    """Read one prepared-run manifest."""

    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"invalid run manifest: {path}")
    return value


def reservation_for_manifest(reservations: list[dict[str, object]],
                             manifest_path: Path) -> dict[str, object]:
    """Find exactly one reservation for a prepared run."""

    matches = [
        item for item in reservations
        if Path(str(item.get("manifest", ""))).resolve() == manifest_path.resolve()
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one reservation for {manifest_path}, found {len(matches)}"
        )
    return matches[0]


def debug_queue_output(args: argparse.Namespace) -> str:
    """Read or query current debug-QOS jobs for the user."""

    if args.squeue_file is not None:
        return Path(args.squeue_file).read_text(encoding="utf-8")
    user = os.environ.get("USER")
    if not user:
        raise ValueError("USER is not set; cannot check existing debug jobs")
    try:
        return subprocess.run(
            [
                "squeue", "-h", "-u", user, "-p", PARTITION, "-q", QOS,
                "-o", "%i|%T|%j",
            ],
            check=True, capture_output=True, text=True,
        ).stdout
    except FileNotFoundError as error:
        raise ValueError("squeue is unavailable; refusing submission check") from error


def check_submit(args: argparse.Namespace) -> int:
    """Fail closed unless a prepared job is budgeted and debug queue is empty."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "prepared":
        raise ValueError("only a prepared, not-yet-submitted job can be checked")
    reservation = reservation_for_manifest(read_reservations(paths), manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is not available for manual submission")
    consumed, reserved = reservation_usage(paths)
    if consumed + reserved > BUDGET_NODE_HOURS:
        raise ValueError("active reservations exceed authorized budget")
    queued = [line for line in debug_queue_output(args).splitlines() if line.strip()]
    if queued:
        raise ValueError(
            "another debug-QOS job is present; sequential submission required: "
            + "; ".join(queued)
        )
    script = manifest["paths"]["batch_script"]
    print("Submission preflight passed. Submit manually, then record its job ID:")
    print(f"  sbatch {shlex.quote(str(script))}")
    print(
        "  python3 scripts/frontier/cgl_lf_frontier.py mark-submitted "
        f"--manifest {shlex.quote(str(manifest_path))} --job-id <jobid>"
    )
    return 0


def mark_submitted(args: argparse.Namespace) -> int:
    """Associate a manually submitted Slurm allocation with its reservation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "prepared":
        raise ValueError("manifest is not in prepared state")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is not in prepared state")
    manifest["state"] = "submitted"
    manifest["job_id"] = args.job_id
    manifest["submitted_recorded_utc"] = utc_now()
    reservation["state"] = "submitted"
    reservation["job_id"] = args.job_id
    write_json(manifest_path, manifest)
    write_json(paths["reservations"], reservations)
    refresh_summary(paths)
    print(f"Marked job {args.job_id} submitted; record allocation use after completion.")
    return 0


def sacct_output(args: argparse.Namespace, paths: dict[str, Path]) -> str:
    """Read or query a pipe-delimited top-level Slurm allocation record."""

    if args.sacct_file is not None:
        output = Path(args.sacct_file).read_text(encoding="utf-8")
    else:
        try:
            output = subprocess.run(
                [
                    "sacct", "-X", "-j", args.job_id,
                    "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                    "ElapsedRaw,Submit,End", "-n", "-P",
                ],
                check=True, capture_output=True, text=True,
            ).stdout
        except FileNotFoundError as error:
            raise ValueError("sacct is unavailable; cannot record allocation") from error
    saved = paths["accounting"] / f"{args.job_id}.sacct.txt"
    saved.write_text(output, encoding="utf-8")
    return output


def parse_sacct(output: str, job_id: str) -> dict[str, str]:
    """Select exactly one top-level allocation row for a completed job."""

    records = []
    for row in csv.reader(output.splitlines(), delimiter="|"):
        if row and row[0] == job_id:
            while row and row[-1] == "":
                row.pop()
            records.append(row)
    if len(records) != 1 or len(records[0]) != 8:
        raise ValueError(
            f"expected exactly one eight-column allocation row for {job_id}"
        )
    keys = (
        "job_id", "job_name", "state", "exit_code", "nodes",
        "elapsed_seconds", "submitted_utc", "completed_utc",
    )
    record = dict(zip(keys, records[0]))
    state = record["state"].split()[0].split("+")[0]
    if state in NONTERMINAL_STATES:
        raise ValueError(f"job {job_id} is not complete: {record['state']}")
    record["state"] = state
    return record


def record(args: argparse.Namespace) -> int:
    """Record one completed allocation and resolve its reservation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "submitted":
        raise ValueError("only a submitted manifest can be accounted")
    if str(manifest.get("job_id")) != args.job_id:
        raise ValueError("job ID does not match submitted manifest")
    ledger = read_ledger(paths)
    if any(row["job_id"] == args.job_id for row in ledger):
        raise ValueError(f"job {args.job_id} is already recorded")
    sacct = parse_sacct(sacct_output(args, paths), args.job_id)
    requested_nodes = int(manifest["allocation"]["nodes"])
    allocated_nodes = int(sacct["nodes"])
    if allocated_nodes != requested_nodes:
        raise ValueError("allocated nodes differ from the prepared reservation")
    actual = node_hours(allocated_nodes, int(sacct["elapsed_seconds"]))
    cumulative = sum(float(row["actual_node_hours"]) for row in ledger) + actual
    if cumulative > BUDGET_NODE_HOURS:
        raise ValueError("recorded actual use exceeds the authorized budget")
    allocation = manifest["allocation"]
    command = manifest["command"]
    run = manifest["run"]
    output_dir = manifest["paths"]["output_dir"]
    row = {
        "job_id": args.job_id,
        "submitted_utc": sacct["submitted_utc"],
        "completed_utc": sacct["completed_utc"],
        "campaign": run["campaign"],
        "run_name": run["run_name"],
        "purpose": run["purpose"],
        "state": sacct["state"],
        "exit_code": sacct["exit_code"],
        "nodes": str(allocated_nodes),
        "requested_walltime": allocation["requested_walltime"],
        "elapsed_seconds": sacct["elapsed_seconds"],
        "reserved_node_hours": f"{float(allocation['reserved_node_hours']):.6f}",
        "actual_node_hours": f"{actual:.6f}",
        "cumulative_actual_node_hours": f"{cumulative:.6f}",
        "git_revision": command["git_revision"],
        "executable_sha256": command["executable_sha256"],
        "input_file": command["input_file"],
        "output_dir": output_dir,
        "result": args.result,
        "notes": args.notes,
    }
    with paths["ledger"].open("a", newline="", encoding="utf-8") as stream:
        csv.DictWriter(stream, fieldnames=LEDGER_COLUMNS).writerow(row)
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    reservation["state"] = "recorded"
    reservation["actual_node_hours"] = actual
    reservation["result"] = args.result
    manifest["state"] = "recorded"
    manifest["accounting"] = row
    write_json(paths["reservations"], reservations)
    write_json(manifest_path, manifest)
    refresh_summary(paths)
    print(
        f"Recorded {actual:.6f} node-hours for job {args.job_id}; "
        f"cumulative={cumulative:.6f}; "
        f"remaining={BUDGET_NODE_HOURS - cumulative:.6f}"
    )
    return 0


def cancel(args: argparse.Namespace) -> int:
    """Release a prepared reservation that has never been submitted."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "prepared":
        raise ValueError("only an unsubmitted prepared reservation may be cancelled")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation has already progressed beyond preparation")
    reservation["state"] = "cancelled"
    reservation["notes"] = args.notes
    manifest["state"] = "cancelled"
    manifest["cancellation_notes"] = args.notes
    write_json(paths["reservations"], reservations)
    write_json(manifest_path, manifest)
    refresh_summary(paths)
    print(f"Cancelled unsubmitted reservation for {manifest_path}")
    return 0


def self_test() -> int:
    """Exercise policy and accounting logic without Frontier or Slurm."""

    with tempfile.TemporaryDirectory(prefix="cgl_lf_frontier_test_") as directory:
        root = Path(directory) / "CGL"
        executable = Path(directory) / "athena"
        executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        executable.chmod(0o750)
        smoke_input = Path(directory) / "cgl_lf_paper_smoke_active_beta10.athinput"
        smoke_input.write_text("<problem>\npaper_grade = false\n", encoding="utf-8")
        empty_queue = Path(directory) / "squeue.txt"
        empty_queue.write_text("", encoding="utf-8")
        arguments = argparse.Namespace(
            root=str(root),
            allow_local_root=True,
            allow_dirty_source=True,
            executable=str(executable),
            input_file=str(smoke_input),
            restart_file=None,
            source_dir=str(ROOT_DIR),
            test_git_revision="self-test-revision",
            nodes=1,
            walltime="00:30:00",
            athena_walltime="00:25:00",
            ranks_per_node=8,
            cpus_per_task=7,
            execution_target="gpu",
            campaign="self-test",
            run_name="valid",
            purpose="offline accounting validation",
            acceptance_criterion="all self-test assertions pass",
            override=[],
        )
        manifest_path = prepare(arguments)
        script_text = Path(str(
            read_manifest(manifest_path)["paths"]["batch_script"]
        )).read_text(encoding="utf-8")
        expected_capture = "grep -E '^(MPICH_|FI_|HSA_|OMP_|ROCR_|HIP_|CRAY_)'"
        if expected_capture not in script_text:
            raise ValueError("self-test failed to retain exported runtime settings")
        if "MPICH_GPU_SUPPORT_ENABLED=1" not in script_text:
            raise ValueError("self-test failed to configure the GPU target")
        check_submit(argparse.Namespace(
            manifest=str(manifest_path), allow_local_root=True,
            squeue_file=str(empty_queue),
        ))
        mark_submitted(argparse.Namespace(
            manifest=str(manifest_path), allow_local_root=True, job_id="12345",
        ))
        sacct = Path(directory) / "sacct.txt"
        sacct.write_text(
            "12345|cgl_self_test|COMPLETED|0:0|1|600|"
            "2026-05-25T00:00:00|2026-05-25T00:10:00|\n",
            encoding="utf-8",
        )
        record(argparse.Namespace(
            manifest=str(manifest_path), allow_local_root=True, job_id="12345",
            sacct_file=str(sacct), result="passed", notes="offline self-test",
        ))
        ledger = read_ledger(layout(root))
        if len(ledger) != 1 or ledger[0]["actual_node_hours"] != "0.166667":
            raise ValueError("self-test failed to calculate actual node hours")
        duplicate_rejected = False
        try:
            record(argparse.Namespace(
                manifest=str(manifest_path), allow_local_root=True, job_id="12345",
                sacct_file=str(sacct), result="passed", notes="duplicate",
            ))
        except ValueError:
            duplicate_rejected = True
        if not duplicate_rejected:
            raise ValueError("self-test failed to reject duplicate accounting")
        restart = Path(directory) / "checkpoint.rst"
        restart.write_text("restart-test-data\n", encoding="utf-8")
        arguments.run_name = "valid_restart"
        arguments.restart_file = str(restart)
        restart_manifest_path = prepare(arguments)
        restart_manifest = read_manifest(restart_manifest_path)
        archived_restart = Path(str(restart_manifest["command"]["restart_file"]))
        if (not archived_restart.is_file()
                or restart_manifest["command"]["restart_sha256"] != sha256(restart)):
            raise ValueError("self-test failed to archive restart provenance")
        cancel(argparse.Namespace(
            manifest=str(restart_manifest_path), allow_local_root=True,
            notes="offline restart preparation check complete",
        ))
        arguments.run_name = "valid_cpu"
        arguments.execution_target = "cpu"
        cpu_manifest_path = prepare(arguments)
        cpu_script_text = Path(str(
            read_manifest(cpu_manifest_path)["paths"]["batch_script"]
        )).read_text(encoding="utf-8")
        if ("MPICH_GPU_SUPPORT_ENABLED=0" not in cpu_script_text
                or "--gpus-per-task" in cpu_script_text
                or "#SBATCH --gpus-per-node" in cpu_script_text
                or "rocm/6.4.2" in cpu_script_text):
            raise ValueError("self-test failed to generate a CPU MPI launch")
        cancel(argparse.Namespace(
            manifest=str(cpu_manifest_path), allow_local_root=True,
            notes="offline CPU MPI preparation check complete",
        ))
        arguments.restart_file = None
        arguments.execution_target = "gpu"
        arguments.run_name = "bad_walltime"
        arguments.walltime = "02:00:01"
        rejected_walltime = False
        try:
            prepare(arguments)
        except ValueError:
            rejected_walltime = True
        if not rejected_walltime:
            raise ValueError("self-test failed to reject long debug walltime")
        production = Path(directory) / "cgl_lf_paper_standard_beta10.athinput"
        production.write_text("<problem>\npaper_grade = true\n", encoding="utf-8")
        arguments.walltime = "00:30:00"
        arguments.run_name = "bad_production"
        arguments.input_file = str(production)
        rejected_production = False
        try:
            prepare(arguments)
        except ValueError:
            rejected_production = True
        if not rejected_production:
            raise ValueError("self-test failed to reject production input")
    print("Frontier CGL-LF campaign utility self-test passed.")
    return 0


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--root", default=str(DEFAULT_ROOT))
    command.add_argument(
        "--allow-local-root", action="store_true",
        help="Permit a non-Frontier root only for offline utility validation.",
    )
    subparsers = command.add_subparsers(dest="action", required=True)
    subparsers.add_parser("init", help="Initialize prescribed campaign layout.")
    prepare_parser = subparsers.add_parser(
        "prepare", help="Reserve budget and create a manual debug submission bundle."
    )
    prepare_parser.add_argument("--campaign", required=True)
    prepare_parser.add_argument("--run-name", required=True)
    prepare_parser.add_argument("--purpose", required=True)
    prepare_parser.add_argument("--acceptance-criterion", required=True)
    prepare_parser.add_argument("--executable", required=True)
    prepare_parser.add_argument("--input-file", required=True)
    prepare_parser.add_argument(
        "--restart-file",
        help="Archive and resume from a retained restart file beneath the run root.",
    )
    prepare_parser.add_argument("--source-dir", default=str(ROOT_DIR))
    prepare_parser.add_argument("--nodes", type=int, required=True)
    prepare_parser.add_argument("--walltime", required=True)
    prepare_parser.add_argument("--athena-walltime", required=True)
    prepare_parser.add_argument("--ranks-per-node", type=int, default=8)
    prepare_parser.add_argument("--cpus-per-task", type=int, default=7)
    prepare_parser.add_argument(
        "--execution-target", choices=("gpu", "cpu"), default="gpu",
        help="Select GPU-aware or Kokkos-Serial CPU/MPI launch configuration.",
    )
    prepare_parser.add_argument("--override", action="append", default=[])
    prepare_parser.add_argument(
        "--allow-dirty-source", action="store_true",
        help="Offline testing only; Frontier preparation requires a clean revision.",
    )
    prepare_parser.add_argument(
        "--test-git-revision",
        help="Offline testing only; bypass git revision lookup.",
    )
    submit_parser = subparsers.add_parser(
        "check-submit", help="Check queue and reservation immediately before sbatch."
    )
    submit_parser.add_argument("--manifest", required=True)
    submit_parser.add_argument(
        "--squeue-file", help="Offline test input in place of querying squeue."
    )
    submitted_parser = subparsers.add_parser(
        "mark-submitted", help="Attach a manually returned Slurm job ID."
    )
    submitted_parser.add_argument("--manifest", required=True)
    submitted_parser.add_argument("--job-id", required=True)
    record_parser = subparsers.add_parser(
        "record", help="Record completed allocation use from sacct."
    )
    record_parser.add_argument("--manifest", required=True)
    record_parser.add_argument("--job-id", required=True)
    record_parser.add_argument("--result", required=True)
    record_parser.add_argument("--notes", default="")
    record_parser.add_argument(
        "--sacct-file", help="Offline test input in place of querying sacct."
    )
    cancel_parser = subparsers.add_parser(
        "cancel", help="Release an unsubmitted prepared reservation."
    )
    cancel_parser.add_argument("--manifest", required=True)
    cancel_parser.add_argument("--notes", required=True)
    subparsers.add_parser("summary", help="Regenerate the accounting summary.")
    subparsers.add_parser("self-test", help="Run offline policy/accounting checks.")
    return command


def main() -> int:
    """Command-line entry point."""

    args = parser().parse_args()
    try:
        if args.action == "self-test":
            return self_test()
        root = require_root(Path(args.root), args.allow_local_root)
        if args.action == "init":
            initialize(root)
            print(f"Initialized Frontier CGL-LF campaign root: {root}")
            return 0
        if args.action == "prepare":
            prepare(args)
            return 0
        if args.action == "check-submit":
            return check_submit(args)
        if args.action == "mark-submitted":
            return mark_submitted(args)
        if args.action == "record":
            return record(args)
        if args.action == "cancel":
            return cancel(args)
        if args.action == "summary":
            paths = initialize(root)
            refresh_summary(paths)
            print(f"Wrote {paths['summary']}")
            return 0
        raise ValueError(f"unsupported action: {args.action}")
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"Frontier CGL-LF campaign utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
