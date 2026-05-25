#!/usr/bin/env python3
"""Prepare and account Frontier MKS24 Stage I production segments.

This utility is separate from ``cgl_lf_frontier.py`` because that launcher is
restricted to non-production ``debug`` QOS work.  Stage I jobs use Frontier's
``batch`` partition with its default production QOS, are prepared and
submitted one segment at a time, and must belong to the tracked mapped-case
manifest.
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


ROOT_DIR = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
DEFAULT_MATRIX = ROOT_DIR / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
ACCOUNT = "AST207"
PARTITION = "batch"
PRODUCTION_QOS = "normal (Frontier default; no -q directive)"
PROJECT_BUDGET_NODE_HOURS = 4000.0
PRIOR_QUALIFICATION_NODE_HOURS = 0.851670
STAGE_I_RESERVED_NODE_HOURS = 1068.888889
MAX_SEGMENT_SECONDS = 24 * 60 * 60
LEDGER_COLUMNS = (
    "job_id",
    "submitted_utc",
    "completed_utc",
    "case_id",
    "case_name",
    "segment",
    "state",
    "exit_code",
    "nodes",
    "requested_walltime",
    "elapsed_seconds",
    "reserved_node_hours",
    "actual_node_hours",
    "cumulative_stage_i_node_hours",
    "executable_revision",
    "executable_sha256",
    "input_revision",
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
    """Return a deterministic UTC timestamp."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_json(path: Path, value: object) -> None:
    """Write stable JSON metadata."""

    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_walltime(value: str) -> int:
    """Parse an HH:MM:SS allocation walltime."""

    match = re.fullmatch(r"(\d{2}):([0-5]\d):([0-5]\d)", value)
    if match is None:
        raise ValueError(f"walltime must use HH:MM:SS: {value}")
    hours, minutes, seconds = (int(item) for item in match.groups())
    return hours * 3600 + minutes * 60 + seconds


def node_hours(nodes: int, seconds: int) -> float:
    """Compute allocated node hours."""

    return nodes * seconds / 3600.0


def require_root(root: Path, allow_local_root: bool) -> Path:
    """Require the declared project run root outside offline validation."""

    resolved = root.expanduser().resolve()
    if resolved != DEFAULT_ROOT and not allow_local_root:
        raise ValueError(
            f"Stage I root must be {DEFAULT_ROOT}; "
            "use --allow-local-root only for offline validation"
        )
    return resolved


def require_beneath_root(path: Path, root: Path, label: str,
                         allow_local_root: bool) -> None:
    """Keep submitted products in the project filesystem."""

    if allow_local_root:
        return
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath {root}: {path}") from error


def layout(root: Path) -> dict[str, Path]:
    """Return retained Stage I production locations."""

    accounting = root / "accounting"
    return {
        "root": root,
        "accounting": accounting,
        "ledger": accounting / "mks24_stage_i_node_hours.csv",
        "reservations": accounting / "mks24_stage_i_reservations.json",
        "summary": accounting / "mks24_stage_i_budget_summary.md",
        "runs": root / "runs" / "mks24-stage-i",
        "logs_slurm": root / "logs" / "slurm",
    }


def initialize(root: Path) -> dict[str, Path]:
    """Create the production layout and accounting stores."""

    paths = layout(root)
    for key in ("accounting", "runs", "logs_slurm"):
        paths[key].mkdir(parents=True, exist_ok=True)
    if not paths["ledger"].exists():
        with paths["ledger"].open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerow(LEDGER_COLUMNS)
    if not paths["reservations"].exists():
        write_json(paths["reservations"], [])
    refresh_summary(paths)
    return paths


def read_ledger(paths: dict[str, Path]) -> list[dict[str, str]]:
    """Read retained Stage I allocation records."""

    with paths["ledger"].open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_reservations(paths: dict[str, Path]) -> list[dict[str, object]]:
    """Read all Stage I preparation records."""

    value = json.loads(paths["reservations"].read_text(encoding="utf-8"))
    if not isinstance(value, list):
        raise ValueError("Stage I reservation store must contain a list")
    return value


def active_reservations(reservations: list[dict[str, object]]
                        ) -> list[dict[str, object]]:
    """Return unaccounted prepared or submitted segment reservations."""

    return [
        item for item in reservations
        if item.get("state") in {"prepared", "submitted"}
    ]


def reservation_usage(paths: dict[str, Path]) -> tuple[float, float]:
    """Return actual and actively reserved Stage I node-hours."""

    actual = sum(float(row["actual_node_hours"]) for row in read_ledger(paths))
    reserved = sum(
        float(item["reserved_node_hours"])
        for item in active_reservations(read_reservations(paths))
    )
    return actual, reserved


def refresh_summary(paths: dict[str, Path]) -> None:
    """Regenerate the human-readable production accounting summary."""

    ledger = read_ledger(paths)
    reservations = read_reservations(paths)
    actual = sum(float(row["actual_node_hours"]) for row in ledger)
    active = active_reservations(reservations)
    reserved = sum(float(item["reserved_node_hours"]) for item in active)
    stage_remaining = STAGE_I_RESERVED_NODE_HOURS - actual - reserved
    project_remaining = (
        PROJECT_BUDGET_NODE_HOURS - PRIOR_QUALIFICATION_NODE_HOURS
        - actual - reserved
    )
    lines = [
        "# MKS24 Stage I Frontier Production Budget",
        "",
        f"- Updated UTC: `{utc_now()}`",
        f"- Project ceiling: `{PROJECT_BUDGET_NODE_HOURS:.6f}` node-hours",
        f"- Previously recorded qualification use: "
        f"`{PRIOR_QUALIFICATION_NODE_HOURS:.6f}` node-hours",
        f"- Stage I production reservation: "
        f"`{STAGE_I_RESERVED_NODE_HOURS:.6f}` node-hours",
        f"- Stage I actual use: `{actual:.6f}` node-hours",
        f"- Active segment reservations: `{reserved:.6f}` node-hours",
        f"- Unreserved Stage I remainder: `{stage_remaining:.6f}` node-hours",
        f"- Project remainder after active Stage I use: "
        f"`{project_remaining:.6f}` node-hours",
        "",
        "## Recorded Segments",
        "",
    ]
    if ledger:
        lines.extend([
            "| Job | Case/segment | State | Node-hours | Result |",
            "| --- | --- | --- | ---: | --- |",
        ])
        for row in ledger:
            lines.append(
                "| `{job_id}` | `{case_id}/{segment}` | {state} | "
                "`{actual_node_hours}` | {result} |".format(**row)
            )
    else:
        lines.append("No Stage I production allocation has been recorded.")
    lines.extend(["", "## Active Reservations", ""])
    if active:
        lines.extend([
            "| Case/segment | Nodes | Walltime | Node-hours | State |",
            "| --- | ---: | --- | ---: | --- |",
        ])
        for item in active:
            lines.append(
                "| `{case_id}/{segment}` | `{nodes}` | `{requested_walltime}` | "
                "`{reserved_node_hours:.6f}` | {state} |".format(**item)
            )
    else:
        lines.append("No prepared or submitted Stage I segment is reserved.")
    lines.extend([
        "",
        "Only one Stage I segment may be prepared or submitted at a time. "
        "Jobs use the `batch` partition with Frontier's default production "
        "`normal` QOS; the `debug` QOS is not used for paper production.",
        "",
    ])
    paths["summary"].write_text("\n".join(lines), encoding="utf-8")


def load_matrix(path: Path) -> dict[str, object]:
    """Read the mapped Stage I matrix."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid matrix manifest: {path}")
    return value


def parse_input_parameters(path: Path) -> dict[str, str]:
    """Read block-qualified input values while ignoring comments."""

    block = ""
    values: dict[str, str] = {}
    for original in path.read_text(encoding="utf-8").splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if "=" in line and block:
            parameter, value = line.split("=", 1)
            values[f"{block}/{parameter.strip()}"] = value.strip()
    return values


def validate_matrix(matrix_path: Path, source_dir: Path) -> dict[str, object]:
    """Validate unique case mapping and admitted file-level aliases."""

    matrix = load_matrix(matrix_path)
    cases = matrix.get("cases")
    if not isinstance(cases, list):
        raise ValueError("matrix cases must be a list")
    expected_ids = [f"R{number:02d}" for number in range(2, 18)]
    ids = [str(case.get("id")) for case in cases]
    if ids != expected_ids:
        raise ValueError(f"Stage I case identifiers must be {expected_ids}: {ids}")
    inputs = [str(case.get("input")) for case in cases]
    if len(set(inputs)) != len(inputs):
        raise ValueError("canonical Stage I inputs are not unique")
    if matrix.get("authorization", {}).get("mapped_unique_runs") != len(cases):
        raise ValueError("mapped run count disagrees with matrix cases")
    for path_text in inputs:
        if not (source_dir / path_text).is_file():
            raise ValueError(f"matrix input does not exist: {path_text}")
    excluded = {
        str(item.get("input"))
        for item in matrix.get("excluded_definitions", [])
    }
    if excluded.intersection(inputs):
        raise ValueError("an excluded input appears in the canonical run matrix")
    for alias in matrix.get("aliases", []):
        equivalent = alias.get("equivalent_input")
        if equivalent is None:
            continue
        canonical = source_dir / str(alias["canonical_input"])
        comparison = source_dir / str(equivalent)
        if not canonical.is_file() or not comparison.is_file():
            raise ValueError(f"alias inputs are missing for {alias['canonical_case']}")
        ignored = set(alias.get("comparison_ignore_parameters", []))
        canonical_values = parse_input_parameters(canonical)
        comparison_values = parse_input_parameters(comparison)
        for key in ignored:
            canonical_values.pop(key, None)
            comparison_values.pop(key, None)
        if canonical_values != comparison_values:
            raise ValueError(
                f"alias is not physically identical for {alias['canonical_case']}"
            )
    return matrix


def case_for_id(matrix: dict[str, object], case_id: str) -> dict[str, object]:
    """Select one unique canonical Stage I case."""

    matches = [
        case for case in matrix["cases"]
        if str(case.get("id")) == case_id
    ]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate Stage I case ID: {case_id}")
    return matches[0]


def git_revision_for_input(source_dir: Path, input_path: Path,
                           matrix_path: Path) -> str:
    """Require submitted input and matrix to match committed source content."""

    revision = subprocess.run(
        ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    relative_paths = [
        str(input_path.relative_to(source_dir)),
        str(matrix_path.relative_to(source_dir)),
    ]
    for diff_args in (["diff", "--quiet", "--"], ["diff", "--cached", "--quiet", "--"]):
        result = subprocess.run(
            ["git", "-C", str(source_dir), *diff_args, *relative_paths],
            check=False,
        )
        if result.returncode != 0:
            raise ValueError(
                "submitted input and matrix must be committed before production"
            )
    return revision


def read_build_provenance(executable: Path,
                          build_manifest: Path) -> dict[str, str]:
    """Verify an archived executable against its immutable build manifest."""

    digest_file = build_manifest / "athena.sha256"
    environment_file = build_manifest / "environment.txt"
    if not digest_file.is_file() or not environment_file.is_file():
        raise ValueError(f"incomplete build manifest directory: {build_manifest}")
    recorded_sha = digest_file.read_text(encoding="utf-8").split()[0]
    actual_sha = sha256(executable)
    if recorded_sha != actual_sha:
        raise ValueError("executable digest does not match build manifest")
    revision_match = re.search(
        r"(?m)^git_revision=([0-9a-f]{40})\s*$",
        environment_file.read_text(encoding="utf-8"),
    )
    if revision_match is None:
        raise ValueError("build manifest does not record a full git revision")
    return {
        "revision": revision_match.group(1),
        "sha256": actual_sha,
        "manifest_dir": str(build_manifest),
    }


def quote(value: str | Path) -> str:
    """Quote a shell literal in the generated Slurm script."""

    return shlex.quote(str(value))


def generated_batch_script(manifest: dict[str, object],
                           manifest_path: Path) -> str:
    """Generate a manually submitted normal-QOS production segment."""

    run = manifest["run"]
    allocation = manifest["allocation"]
    command = manifest["command"]
    paths = manifest["paths"]
    overrides = " ".join(quote(value) for value in command["overrides"])
    if overrides:
        overrides = " " + overrides
    restart = command.get("restart_file")
    restart_literal = quote(str(restart)) if restart else "''"
    job_name = f"cgl_mks24_{run['case_id']}_{run['segment']}"
    job_name = re.sub(r"[^A-Za-z0-9_]+", "_", job_name)[:60]
    return f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -A {ACCOUNT}
#SBATCH -o {paths["slurm_log"]}
#SBATCH -p {PARTITION}
# Frontier default normal QOS is required for production; do not add -q debug.
#SBATCH -t {allocation["requested_walltime"]}
#SBATCH -N {allocation["nodes"]}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -euo pipefail

RUN_MANIFEST={quote(manifest_path)}
ATHENA={quote(command["executable"])}
INPUT={quote(command["input_file"])}
RESTART={restart_literal}
OUT_DIR={quote(paths["output_dir"])}
ENV_LOG={quote(paths["environment_log"])}
RANKS_PER_NODE={allocation["ranks_per_node"]}
CPUS_PER_TASK={allocation["cpus_per_task"]}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * RANKS_PER_NODE))"

test -f "${{RUN_MANIFEST}}"
test -x "${{ATHENA}}"
test -f "${{INPUT}}"
if [[ -n "${{RESTART}}" ]]; then
  test -f "${{RESTART}}"
fi
mkdir -p "${{OUT_DIR}}"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
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
  --gpus-per-task=1 --gpu-bind=closest \
  "${{ATHENA}}" "${{RUN_ARGS[@]}}" -d "${{OUT_DIR}}" \
  -t {quote(command["athena_walltime"])} \
  job/basename={quote(run["run_basename"])}{overrides}

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "${{ENV_LOG}}"
"""


def prepare(args: argparse.Namespace) -> Path:
    """Create one retained, sequentially submitted production segment."""

    root = require_root(Path(args.root), args.allow_local_root)
    paths = initialize(root)
    if active_reservations(read_reservations(paths)):
        raise ValueError(
            "another Stage I segment is prepared or submitted; "
            "record or cancel it before preparing a new segment"
        )
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    executable = Path(args.executable).expanduser().resolve()
    build_manifest = Path(args.build_manifest).expanduser().resolve()
    restart = (
        Path(args.restart_file).expanduser().resolve()
        if args.restart_file else None
    )
    require_beneath_root(source_dir, root, "source directory", args.allow_local_root)
    require_beneath_root(executable, root, "executable", args.allow_local_root)
    require_beneath_root(
        build_manifest, root, "build manifest", args.allow_local_root
    )
    if restart is not None:
        require_beneath_root(restart, root, "restart file", args.allow_local_root)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"executable is unavailable: {executable}")
    matrix = validate_matrix(matrix_path, source_dir)
    case = case_for_id(matrix, args.case_id)
    input_path = source_dir / str(case["input"])
    input_revision = git_revision_for_input(source_dir, input_path, matrix_path)
    provenance = read_build_provenance(executable, build_manifest)
    if restart is not None and not restart.is_file():
        raise ValueError(f"restart file is unavailable: {restart}")
    if args.nodes < 1:
        raise ValueError("--nodes must be positive")
    requested_seconds = parse_walltime(args.walltime)
    if requested_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError("a Stage I segment may not request more than 24 hours")
    athena_seconds = parse_walltime(args.athena_walltime)
    if athena_seconds > requested_seconds - 600:
        raise ValueError("Athena walltime must leave ten minutes for shutdown")
    segment_hours = node_hours(args.nodes, requested_seconds)
    actual, reserved = reservation_usage(paths)
    if actual + reserved + segment_hours > STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("proposed segment exceeds the Stage I reservation")
    if (
        PRIOR_QUALIFICATION_NODE_HOURS + actual + reserved + segment_hours
        > PROJECT_BUDGET_NODE_HOURS
    ):
        raise ValueError("proposed segment exceeds the project ceiling")
    run_dir = paths["runs"] / args.case_id / args.segment
    manifest_dir = run_dir / "manifest"
    manifest_path = manifest_dir / "prepared_run.json"
    if manifest_path.exists():
        raise ValueError(f"segment is already prepared: {manifest_path}")
    output_dir = run_dir / "output"
    for path in (manifest_dir, output_dir):
        path.mkdir(parents=True, exist_ok=True)
    archived_input = manifest_dir / "submitted_input.athinput"
    archived_matrix = manifest_dir / "mks24_stage_i_manifest.json"
    shutil.copy2(input_path, archived_input)
    shutil.copy2(matrix_path, archived_matrix)
    archived_restart = None
    if restart is not None:
        archived_restart = manifest_dir / "submitted_restart.rst"
        shutil.copy2(restart, archived_restart)
    batch_script = manifest_dir / "cgl_lf_stage_i.sbatch"
    manifest: dict[str, object] = {
        "schema_version": 1,
        "state": "prepared",
        "prepared_utc": utc_now(),
        "project_root": str(root),
        "policy": {
            "account": ACCOUNT,
            "partition": PARTITION,
            "qos": PRODUCTION_QOS,
            "project_budget_node_hours": PROJECT_BUDGET_NODE_HOURS,
            "prior_qualification_node_hours": PRIOR_QUALIFICATION_NODE_HOURS,
            "stage_i_reserved_node_hours": STAGE_I_RESERVED_NODE_HOURS,
            "sequential_manual_submission_required": True,
        },
        "run": {
            "case_id": args.case_id,
            "case_name": case["name"],
            "segment": args.segment,
            "run_basename": f"{case['name']}_{args.segment}",
            "resolution": case["resolution"],
            "figure_roles": case["figure_roles"],
            "acceptance_criterion": args.acceptance_criterion,
        },
        "allocation": {
            "nodes": args.nodes,
            "requested_walltime": args.walltime,
            "requested_seconds": requested_seconds,
            "reserved_node_hours": segment_hours,
            "ranks_per_node": args.ranks_per_node,
            "cpus_per_task": args.cpus_per_task,
        },
        "command": {
            "source_dir": str(source_dir),
            "input_revision": input_revision,
            "source_input_file": str(input_path),
            "input_file": str(archived_input),
            "input_sha256": sha256(archived_input),
            "matrix_file": str(archived_matrix),
            "matrix_sha256": sha256(archived_matrix),
            "executable": str(executable),
            "executable_revision": provenance["revision"],
            "executable_sha256": provenance["sha256"],
            "build_manifest": provenance["manifest_dir"],
            "source_restart_file": str(restart) if restart else None,
            "restart_file": str(archived_restart) if archived_restart else None,
            "restart_sha256": sha256(archived_restart) if archived_restart else None,
            "athena_walltime": args.athena_walltime,
            "overrides": args.override,
        },
        "paths": {
            "run_dir": str(run_dir),
            "output_dir": str(output_dir),
            "environment_log": str(manifest_dir / "run_environment.txt"),
            "batch_script": str(batch_script),
            "slurm_log": str(paths["logs_slurm"] / "%x.%j.log"),
        },
    }
    batch_script.write_text(
        generated_batch_script(manifest, manifest_path), encoding="utf-8"
    )
    batch_script.chmod(0o750)
    write_json(manifest_path, manifest)
    reservations = read_reservations(paths)
    reservations.append({
        "manifest": str(manifest_path),
        "case_id": args.case_id,
        "case_name": case["name"],
        "segment": args.segment,
        "nodes": args.nodes,
        "requested_walltime": args.walltime,
        "reserved_node_hours": segment_hours,
        "state": "prepared",
        "prepared_utc": manifest["prepared_utc"],
    })
    write_json(paths["reservations"], reservations)
    refresh_summary(paths)
    print(f"Prepared Stage I segment: {manifest_path}")
    print(f"Reserved node-hours: {segment_hours:.6f}")
    return manifest_path


def read_manifest(path: Path) -> dict[str, object]:
    """Read one segment preparation manifest."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid run manifest: {path}")
    return value


def reservation_for_manifest(reservations: list[dict[str, object]],
                             manifest_path: Path) -> dict[str, object]:
    """Return exactly one reservation corresponding to a run manifest."""

    matches = [
        item for item in reservations
        if Path(str(item.get("manifest", ""))).resolve() == manifest_path.resolve()
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one reservation for {manifest_path}, found {len(matches)}"
        )
    return matches[0]


def production_queue_output(args: argparse.Namespace) -> str:
    """Query this user's currently queued Stage I production jobs."""

    if args.squeue_file:
        return Path(args.squeue_file).read_text(encoding="utf-8")
    user = os.environ.get("USER")
    if not user:
        raise ValueError("USER is not set; cannot check production queue")
    return subprocess.run(
        ["squeue", "-h", "-u", user, "-p", PARTITION, "-o", "%i|%T|%j"],
        check=True, capture_output=True, text=True,
    ).stdout


def check_submit(args: argparse.Namespace) -> int:
    """Fail closed unless a prepared segment may be manually submitted."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "prepared":
        raise ValueError("only a prepared segment can be submission-checked")
    reservation = reservation_for_manifest(read_reservations(paths), manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is unavailable for submission")
    actual, reserved = reservation_usage(paths)
    if actual + reserved > STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("active Stage I reservation exceeds its ceiling")
    lines = [
        line for line in production_queue_output(args).splitlines()
        if line.strip() and "cgl_mks24_" in line
    ]
    if lines:
        raise ValueError(
            "another Stage I production segment is queued: " + "; ".join(lines)
        )
    script = str(manifest["paths"]["batch_script"])
    if not args.allow_local_root and not args.skip_slurm_test:
        result = subprocess.run(
            ["sbatch", "--test-only", script],
            check=True, capture_output=True, text=True,
        )
        print(result.stdout.strip())
    print("Submission preflight passed. Submit manually, then record its job ID:")
    print(f"  sbatch {shlex.quote(script)}")
    print(
        "  python3 scripts/frontier/cgl_lf_stage_i.py mark-submitted "
        f"--manifest {shlex.quote(str(manifest_path))} --job-id <jobid>"
    )
    return 0


def mark_submitted(args: argparse.Namespace) -> int:
    """Attach a Slurm job ID to the reserved Stage I segment."""

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
    print(f"Marked Stage I job {args.job_id} submitted.")
    return 0


def sacct_output(args: argparse.Namespace, paths: dict[str, Path]) -> str:
    """Read or query a top-level Slurm allocation record."""

    if args.sacct_file:
        output = Path(args.sacct_file).read_text(encoding="utf-8")
    else:
        output = subprocess.run(
            [
                "sacct", "-X", "-j", args.job_id,
                "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                "ElapsedRaw,Submit,End", "-n", "-P",
            ],
            check=True, capture_output=True, text=True,
        ).stdout
    (paths["accounting"] / f"{args.job_id}.stage_i.sacct.txt").write_text(
        output, encoding="utf-8"
    )
    return output


def parse_sacct(output: str, job_id: str) -> dict[str, str]:
    """Select a completed top-level allocation row."""

    records: list[list[str]] = []
    for row in csv.reader(output.splitlines(), delimiter="|"):
        if row and row[0] == job_id:
            while row and not row[-1]:
                row.pop()
            records.append(row)
    if len(records) != 1 or len(records[0]) != 8:
        raise ValueError(f"expected one eight-column allocation row for {job_id}")
    keys = (
        "job_id", "job_name", "state", "exit_code", "nodes",
        "elapsed_seconds", "submitted_utc", "completed_utc",
    )
    result = dict(zip(keys, records[0]))
    state = result["state"].split()[0].split("+")[0]
    if state in NONTERMINAL_STATES:
        raise ValueError(f"job {job_id} is not complete: {result['state']}")
    result["state"] = state
    return result


def record(args: argparse.Namespace) -> int:
    """Account a completed segment and release its reservation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "submitted":
        raise ValueError("only a submitted segment can be accounted")
    if str(manifest.get("job_id")) != args.job_id:
        raise ValueError("job ID does not match the submitted segment")
    ledger = read_ledger(paths)
    if any(row["job_id"] == args.job_id for row in ledger):
        raise ValueError(f"job {args.job_id} is already accounted")
    sacct = parse_sacct(sacct_output(args, paths), args.job_id)
    nodes = int(sacct["nodes"])
    if nodes != int(manifest["allocation"]["nodes"]):
        raise ValueError("allocated nodes differ from the prepared reservation")
    actual = node_hours(nodes, int(sacct["elapsed_seconds"]))
    cumulative = sum(float(row["actual_node_hours"]) for row in ledger) + actual
    if cumulative > STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("actual use exceeds the Stage I reservation")
    command = manifest["command"]
    run = manifest["run"]
    allocation = manifest["allocation"]
    row = {
        "job_id": args.job_id,
        "submitted_utc": sacct["submitted_utc"],
        "completed_utc": sacct["completed_utc"],
        "case_id": run["case_id"],
        "case_name": run["case_name"],
        "segment": run["segment"],
        "state": sacct["state"],
        "exit_code": sacct["exit_code"],
        "nodes": str(nodes),
        "requested_walltime": allocation["requested_walltime"],
        "elapsed_seconds": sacct["elapsed_seconds"],
        "reserved_node_hours": f"{float(allocation['reserved_node_hours']):.6f}",
        "actual_node_hours": f"{actual:.6f}",
        "cumulative_stage_i_node_hours": f"{cumulative:.6f}",
        "executable_revision": command["executable_revision"],
        "executable_sha256": command["executable_sha256"],
        "input_revision": command["input_revision"],
        "input_file": command["input_file"],
        "output_dir": manifest["paths"]["output_dir"],
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
        f"Recorded {actual:.6f} node-hours for {run['case_id']}/{run['segment']}; "
        f"Stage I cumulative={cumulative:.6f}."
    )
    return 0


def cancel(args: argparse.Namespace) -> int:
    """Release a segment that was prepared but never submitted."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    paths = initialize(root)
    if manifest.get("state") != "prepared":
        raise ValueError("only an unsubmitted segment may be cancelled")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation has progressed beyond preparation")
    reservation["state"] = "cancelled"
    reservation["notes"] = args.notes
    manifest["state"] = "cancelled"
    manifest["cancellation_notes"] = args.notes
    write_json(paths["reservations"], reservations)
    write_json(manifest_path, manifest)
    refresh_summary(paths)
    print(f"Cancelled Stage I reservation: {manifest_path}")
    return 0


def parser() -> argparse.ArgumentParser:
    """Build command-line parsing."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--root", default=str(DEFAULT_ROOT))
    command.add_argument(
        "--allow-local-root", action="store_true",
        help="Permit a non-project root only for offline validation.",
    )
    actions = command.add_subparsers(dest="action", required=True)
    validate = actions.add_parser("validate-matrix")
    validate.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    validate.add_argument("--source-dir", default=str(ROOT_DIR))
    actions.add_parser("init")
    prepare_parser = actions.add_parser("prepare")
    prepare_parser.add_argument("--case-id", required=True)
    prepare_parser.add_argument("--segment", required=True)
    prepare_parser.add_argument("--acceptance-criterion", required=True)
    prepare_parser.add_argument("--executable", required=True)
    prepare_parser.add_argument("--build-manifest", required=True)
    prepare_parser.add_argument("--source-dir", default=str(ROOT_DIR))
    prepare_parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    prepare_parser.add_argument("--restart-file")
    prepare_parser.add_argument("--nodes", type=int, required=True)
    prepare_parser.add_argument("--walltime", required=True)
    prepare_parser.add_argument("--athena-walltime", required=True)
    prepare_parser.add_argument("--ranks-per-node", type=int, default=8)
    prepare_parser.add_argument("--cpus-per-task", type=int, default=7)
    prepare_parser.add_argument("--override", action="append", default=[])
    submit_parser = actions.add_parser("check-submit")
    submit_parser.add_argument("--manifest", required=True)
    submit_parser.add_argument("--squeue-file")
    submit_parser.add_argument("--skip-slurm-test", action="store_true")
    submitted = actions.add_parser("mark-submitted")
    submitted.add_argument("--manifest", required=True)
    submitted.add_argument("--job-id", required=True)
    recorded = actions.add_parser("record")
    recorded.add_argument("--manifest", required=True)
    recorded.add_argument("--job-id", required=True)
    recorded.add_argument("--result", required=True)
    recorded.add_argument("--notes", default="")
    recorded.add_argument("--sacct-file")
    cancelled = actions.add_parser("cancel")
    cancelled.add_argument("--manifest", required=True)
    cancelled.add_argument("--notes", required=True)
    actions.add_parser("summary")
    return command


def main() -> int:
    """Command-line entry point."""

    args = parser().parse_args()
    try:
        if args.action == "validate-matrix":
            matrix = validate_matrix(
                Path(args.matrix).expanduser().resolve(),
                Path(args.source_dir).expanduser().resolve(),
            )
            print(f"Validated {len(matrix['cases'])} mapped Stage I cases.")
            return 0
        root = require_root(Path(args.root), args.allow_local_root)
        if args.action == "init":
            initialize(root)
            print(f"Initialized Stage I production accounting beneath {root}.")
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
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"Stage I production utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
