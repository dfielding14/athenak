#!/usr/bin/env python3
"""Run repeatable CPU/MPI CR tracer particle benchmark cases."""

from __future__ import annotations

import argparse
import json
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import Sequence

import cr_tracer_inspect


ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "inputs" / "particles"

CPU_TIME_RE = re.compile(r"cpu time used\s*=\s*([0-9.eE+-]+)")
PUPS_RE = re.compile(r"particle-updates/cpu_second\s*=\s*([0-9.eE+-]+)")
MBC_RE = re.compile(r"MeshBlock-cycles\s*=\s*([0-9.eE+-]+)")
AMR_RE = re.compile(r"(\d+) MeshBlocks created, (\d+) deleted by AMR")
LIVE_MB_RE = re.compile(r"Current number of MeshBlocks\s*=\s*(\d+)")
PERF_RE = re.compile(r"Particle performance ([^:]+):\s*(.*)")


BASE_CASES: Dict[str, Dict[str, object]] = {
    "push": {
        "input": "random_particle_drift.athinput",
        "description": "drift push throughput without AMR",
        "overrides": [
            "job/basename=bench_push",
            "mesh/nx1=32", "mesh/nx2=32", "mesh/nx3=32",
            "meshblock/nx1=16", "meshblock/nx2=16", "meshblock/nx3=16",
            "particles/ppc=0.25",
            "particles/nspecies=2",
            "particles/log_performance=true",
            "particles/check_consistency_mode=none",
            "problem/particle_position=random",
            "problem/particle_velocity=random",
            "time/nlim=8",
            "time/tlim=0.2",
            "output1/dt=1000.0",
            "output2/dt=1000.0",
        ],
    },
    "amr_remap": {
        "input": "cr_tracer_drift_amr_perf.athinput",
        "description": "zero-velocity drift AMR remap throughput",
        "overrides": [
            "job/basename=bench_amr_remap",
            "particles/log_performance=true",
            "particles/check_consistency_mode=none",
        ],
    },
    "exchange": {
        "input": "random_particle_drift.athinput",
        "description": "high-velocity subcycled drift exchange throughput",
        "requires_mpi": True,
        "overrides": [
            "job/basename=bench_exchange",
            "mesh/nx1=64", "mesh/nx2=8", "mesh/nx3=8",
            "meshblock/nx1=8", "meshblock/nx2=8", "meshblock/nx3=8",
            "particles/ppc=0.125",
            "particles/nspecies=2",
            "particles/log_performance=true",
            "particles/check_consistency_mode=counts",
            "particles/subcycle=true",
            "particles/subcycle_max_steps=128",
            "problem/particle_position=random",
            "problem/particle_velocity=uniform",
            "problem/v0x=20.0",
            "problem/v0y=0.0",
            "problem/v0z=0.0",
            "time/nlim=4",
            "time/tlim=0.2",
            "output1/dt=1000.0",
            "output2/dt=1000.0",
        ],
    },
    "histogram": {
        "input": "cr_tracer_boris_amr_stress.athinput",
        "description": (
            "reduced CR histogram, spectrum, joint-spectrum, and moment output"),
        "overrides": [
            "job/basename=bench_histogram",
            "particles/log_performance=true",
            "particles/check_consistency_mode=none",
            "time/nlim=4",
            "time/tlim=0.04",
        ],
    },
    "histogram_per_rank": {
        "input": "cr_tracer_boris_amr_stress.athinput",
        "description": (
            "per-rank CR histogram, spectrum, joint-spectrum, and moment output"),
        "overrides": [
            "job/basename=bench_histogram_rank",
            "particles/log_performance=true",
            "particles/check_consistency_mode=none",
            "time/nlim=4",
            "time/tlim=0.04",
            "output3/df_single_file_per_rank=1",
            "output4/dxhist_single_file_per_rank=1",
            "output5/pmom_single_file_per_rank=1",
            "output6/drh_single_file_per_rank=1",
            "output7/dparh_single_file_per_rank=1",
            "output8/pspec_single_file_per_rank=1",
            "output9/pspec2_single_file_per_rank=1",
        ],
    },
    "restart": {
        "input": "random_particle_drift.athinput",
        "description": "particle restart write and inspector read throughput",
        "inspect_restart": True,
        "overrides": [
            "job/basename=bench_restart",
            "mesh/nx1=32", "mesh/nx2=32", "mesh/nx3=32",
            "meshblock/nx1=16", "meshblock/nx2=16", "meshblock/nx3=16",
            "particles/ppc=0.25",
            "particles/nspecies=2",
            "particles/log_performance=true",
            "particles/check_consistency_mode=none",
            "problem/particle_position=random",
            "problem/particle_velocity=random",
            "time/nlim=2",
            "time/tlim=0.05",
            "output1/file_type=prst",
            "output1/dt=0.01",
            "output2/dt=1000.0",
        ],
    },
    "boris_amr": {
        "input": "cr_tracer_boris_amr_perf.athinput",
        "description": "end-to-end Boris, MHD, AMR, and particle migration",
        "overrides": [
            "job/basename=bench_boris_amr",
            "particles/log_performance=true",
            "particles/check_consistency_mode=none",
        ],
    },
}


MODE_OVERRIDES: Dict[str, List[str]] = {
    "default": [
        "particles/amr_remap=device_table",
        "particles/exchange_mode=alltoall_counts",
    ],
    "host_tree": [
        "particles/amr_remap=host_tree",
        "particles/exchange_mode=alltoall_counts",
    ],
    "allgather": [
        "particles/amr_remap=device_table",
        "particles/exchange_mode=allgather",
    ],
    "fallback": [
        "particles/amr_remap=host_tree",
        "particles/exchange_mode=allgather",
    ],
}


def _split_csv(text: str) -> List[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def _parse_int_csv(text: str) -> List[int]:
    return [int(item) for item in _split_csv(text)]


def _default_athena() -> Path:
    candidates = [
        ROOT / "build-cr-robust-mpi" / "src" / "athena",
        ROOT / "build-cr-perf-mpi" / "src" / "athena",
        ROOT / "build" / "src" / "athena",
        ROOT / "tst" / "build" / "src" / "athena",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _metadata(args: argparse.Namespace) -> Dict[str, object]:
    metadata: Dict[str, object] = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": sys.version.split()[0],
        "athena": str(args.athena),
        "mpiexec": args.mpiexec,
    }
    for key, command in {
        "git_branch": ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        "git_commit": ["git", "rev-parse", "HEAD"],
        "mpi_version": [args.mpiexec, "--version"],
    }.items():
        try:
            proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                  timeout=10)
            if proc.returncode == 0:
                metadata[key] = (proc.stdout or proc.stderr).splitlines()[0]
        except (OSError, subprocess.SubprocessError):
            pass
    return metadata


def _parse_perf_fields(text: str) -> Dict[str, int]:
    fields: Dict[str, int] = {}
    for item in text.split():
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        try:
            fields[key] = int(value)
        except ValueError:
            continue
    return fields


def _parse_output(text: str) -> Dict[str, object]:
    metrics: Dict[str, object] = {"performance": {}}
    if match := CPU_TIME_RE.search(text):
        metrics["cpu_time"] = float(match.group(1))
    if match := PUPS_RE.search(text):
        metrics["particle_updates_per_cpu_second"] = float(match.group(1))
    if match := MBC_RE.search(text):
        metrics["meshblock_cycles"] = float(match.group(1))
    if match := AMR_RE.search(text):
        metrics["meshblocks_created"] = int(match.group(1))
        metrics["meshblocks_deleted"] = int(match.group(2))
    if match := LIVE_MB_RE.search(text):
        metrics["meshblocks_live"] = int(match.group(1))

    perf: Dict[str, Dict[str, int]] = {}
    for match in PERF_RE.finditer(text):
        label = match.group(1)
        values = _parse_perf_fields(match.group(2))
        current = perf.setdefault(label, {})
        for key, value in values.items():
            current[key] = current.get(key, 0) + value
    metrics["performance"] = perf
    return metrics


def _median(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    return float(statistics.median(values))


def _summarize_runs(runs: Sequence[Dict[str, object]]) -> Dict[str, object]:
    numeric_keys = [
        "wall_time",
        "cpu_time",
        "particle_updates_per_cpu_second",
        "meshblock_cycles",
        "meshblocks_created",
        "meshblocks_deleted",
        "meshblocks_live",
        "restart_read_wall_time",
        "restart_total",
    ]
    summary: Dict[str, object] = {}
    for key in numeric_keys:
        values = [run[key] for run in runs if isinstance(run.get(key), (int, float))]
        value = _median([float(v) for v in values])
        if value is not None:
            summary[key + "_median"] = value

    labels = sorted({
        label for run in runs
        for label in run.get("performance", {}).keys()
    })
    perf_summary: Dict[str, Dict[str, float]] = {}
    for label in labels:
        field_names = sorted({
            field for run in runs
            for field in run.get("performance", {}).get(label, {}).keys()
        })
        perf_summary[label] = {}
        for field in field_names:
            values = [
                run.get("performance", {}).get(label, {}).get(field)
                for run in runs
                if field in run.get("performance", {}).get(label, {})
            ]
            median = _median([float(v) for v in values if isinstance(v, (int, float))])
            if median is not None:
                perf_summary[label][field + "_median"] = median
    summary["performance"] = perf_summary

    sent = 0.0
    bytes_sent = 0.0
    messages = 0.0
    for label_data in perf_summary.values():
        sent += float(label_data.get("sent_median", 0.0))
        bytes_sent += float(label_data.get("bytes_median", 0.0))
        messages += float(label_data.get("messages_median", 0.0))
    if sent > 0.0:
        summary["bytes_per_sent_particle_median"] = bytes_sent/sent
        summary["messages_per_exchange_median"] = messages
    return summary


def _run_one(args: argparse.Namespace, case_name: str, case: Dict[str, object],
             mode_name: str, rank: int, repeat: int, run_root: Path) -> Dict[str, object]:
    run_dir = run_root / f"{case_name}_{mode_name}_n{rank}_r{repeat}"
    run_dir.mkdir(parents=True, exist_ok=True)
    command = [str(args.athena), "-i", str(INPUTS / str(case["input"])), "-d",
               str(run_dir)]
    if rank > 1 or args.mpi_n1:
        command = [args.mpiexec, "-n", str(rank)] + command
    command += list(case["overrides"])
    command += MODE_OVERRIDES[mode_name]

    start = time.perf_counter()
    proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    wall = time.perf_counter() - start
    output = proc.stdout + proc.stderr
    if proc.returncode != 0:
        raise RuntimeError(
            f"Benchmark case {case_name}/{mode_name}/n{rank}/r{repeat} failed\n"
            + output)

    result = _parse_output(output)
    result.update({
        "case": case_name,
        "mode": mode_name,
        "rank": rank,
        "repeat": repeat,
        "wall_time": wall,
        "run_dir": str(run_dir),
    })
    if case.get("inspect_restart"):
        read_start = time.perf_counter()
        restart = cr_tracer_inspect.summarize_restart(run_dir)
        result["restart_read_wall_time"] = time.perf_counter() - read_start
        result["restart_total"] = restart["total"]
        result["restart_rank_counts"] = restart["rank_counts"]
        result["restart_species_counts"] = restart["species_counts"]
    if not args.keep_runs:
        shutil.rmtree(run_dir, ignore_errors=True)
    return result


def _write_markdown(path: Path, data: Dict[str, object]) -> None:
    lines = [
        "# CR Tracer Benchmark Results",
        "",
        "## Metadata",
        "",
    ]
    for key, value in data["metadata"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += ["", "## Medians", ""]
    lines.append(
        "| Case | Mode | Ranks | CPU s | Wall s | PUPS | Sent | Bytes | Messages |")
    lines.append(
        "|------|------|-------|-------|--------|------|------|-------|----------|")
    for row in data["summaries"]:
        summary = row["summary"]
        perf = summary.get("performance", {})
        exchange = perf.get("particle_mpi_exchange", {})
        sent = exchange.get("sent_median", "")
        bytes_sent = exchange.get("bytes_median", "")
        messages = exchange.get("messages_median", "")
        lines.append(
            f"| {row['case']} | {row['mode']} | {row['rank']} | "
            f"{summary.get('cpu_time_median', '')} | "
            f"{summary.get('wall_time_median', '')} | "
            f"{summary.get('particle_updates_per_cpu_second_median', '')} | "
            f"{sent} | {bytes_sent} | {messages} |")
    path.write_text("\n".join(lines) + "\n")


def run_benchmarks(args: argparse.Namespace) -> Dict[str, object]:
    selected_cases = _split_csv(args.cases)
    selected_modes = _split_csv(args.modes)
    ranks = _parse_int_csv(args.ranks)
    unknown_cases = sorted(set(selected_cases) - set(BASE_CASES))
    unknown_modes = sorted(set(selected_modes) - set(MODE_OVERRIDES))
    if unknown_cases:
        raise ValueError(f"Unknown benchmark cases: {', '.join(unknown_cases)}")
    if unknown_modes:
        raise ValueError(f"Unknown benchmark modes: {', '.join(unknown_modes)}")

    run_root = Path(args.run_root) if args.run_root else Path(
        tempfile.mkdtemp(prefix="cr_tracer_bench_", dir="/tmp"))
    run_root.mkdir(parents=True, exist_ok=True)
    runs: List[Dict[str, object]] = []
    try:
        for case_name in selected_cases:
            case = BASE_CASES[case_name]
            for mode_name in selected_modes:
                for rank in ranks:
                    if case.get("requires_mpi") and rank < 2:
                        continue
                    for repeat in range(args.repeats):
                        print(f"running case={case_name} mode={mode_name} "
                              f"rank={rank} repeat={repeat + 1}/{args.repeats}",
                              flush=True)
                        runs.append(_run_one(args, case_name, case, mode_name, rank,
                                             repeat, run_root))
    finally:
        if not args.keep_runs and args.run_root is None:
            shutil.rmtree(run_root, ignore_errors=True)

    summaries = []
    for case_name in selected_cases:
        for mode_name in selected_modes:
            for rank in ranks:
                grouped = [
                    run for run in runs
                    if run["case"] == case_name and run["mode"] == mode_name
                    and run["rank"] == rank
                ]
                if grouped:
                    summaries.append({
                        "case": case_name,
                        "mode": mode_name,
                        "rank": rank,
                        "summary": _summarize_runs(grouped),
                    })
    return {
        "metadata": _metadata(args),
        "cases": {name: BASE_CASES[name]["description"] for name in selected_cases},
        "runs": runs,
        "summaries": summaries,
    }


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run CR tracer CPU/MPI benchmark cases and summarize medians.")
    parser.add_argument("--athena", type=Path, default=_default_athena())
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--ranks", default="1,2,4")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--cases",
        default="push,amr_remap,exchange,histogram,histogram_per_rank,restart")
    parser.add_argument("--modes", default="default")
    parser.add_argument("--mpi-n1", action="store_true",
                        help="Launch one-rank cases through mpiexec -n 1")
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--keep-runs", action="store_true")
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.repeats < 1:
        raise ValueError("--repeats must be >= 1")
    if not args.athena.exists():
        raise FileNotFoundError(f"Athena executable not found: {args.athena}")

    data = run_benchmarks(args)
    text = json.dumps(data, indent=2, sort_keys=True)
    if args.output_json:
        args.output_json.write_text(text + "\n")
    else:
        print(text)
    if args.output_md:
        _write_markdown(args.output_md, data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
