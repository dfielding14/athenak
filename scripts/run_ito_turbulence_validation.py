#!/usr/bin/env python3
"""Run paired old/corrected Ito-2 turbulence validation experiments."""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import os
import platform
import re
import resource
import shlex
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from analyze_ito_turbulence_validation import BIN_CONVERT, discover_output


REPO_ROOT = Path(__file__).resolve().parents[1]
FORMAL_RK2_CFL_LIMIT = 1.0
DEFAULT_ITO_PROBABILITY_TARGET = 0.99
INITIAL_3D_CFL_GUARD = DEFAULT_ITO_PROBABILITY_TARGET / 3.0
FAILURE_PATTERN = re.compile(r"FATAL ERROR|Segmentation fault|\bnan\b", re.IGNORECASE)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(*args: str) -> str | None:
    process = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return process.stdout.strip() if process.returncode == 0 else None


@functools.lru_cache(maxsize=None)
def executable_metadata(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    stat = resolved.stat()
    configuration = subprocess.run(
        [str(resolved), "-c"],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "configuration_return_code": configuration.returncode,
        "configuration": configuration.stdout,
    }


def capture_repo_state(output_root: Path) -> dict[str, Any]:
    status = subprocess.run(
        ["git", "status", "--short", "--branch"],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    ).stdout
    patch = subprocess.run(
        ["git", "diff", "--binary", "HEAD"],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    ).stdout
    status_path = output_root / "source_status.txt"
    patch_path = output_root / "source_state.patch"
    status_path.write_text(status, encoding="utf-8")
    patch_path.write_bytes(patch)
    return {
        "head": git_value("rev-parse", "HEAD"),
        "status_file": str(status_path),
        "patch_file": str(patch_path),
        "patch_sha256": hashlib.sha256(patch).hexdigest(),
    }


def command_for_run(
    executable: Path,
    input_file: Path,
    run_dir: Path,
    launcher: str,
    overrides: list[str],
) -> list[str]:
    return [
        *shlex.split(launcher),
        str(executable.expanduser().resolve()),
        "-i",
        str(input_file.expanduser().resolve()),
        "-d",
        str(run_dir.expanduser().resolve()),
        *overrides,
    ]


def run_one(
    label: str,
    executable: Path,
    input_file: Path,
    run_dir: Path,
    launcher: str,
    overrides: list[str],
    dry_run: bool,
) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=True)
    command = command_for_run(executable, input_file, run_dir, launcher, overrides)
    metadata: dict[str, Any] = {
        "label": label,
        "command": command,
        "command_shell": shlex.join(command),
        "input": {
            "path": str(input_file.resolve()),
            "sha256": sha256_file(input_file),
        },
        "executable": executable_metadata(executable),
        "launcher": launcher,
        "overrides": overrides,
        "repo": {
            "root": str(REPO_ROOT),
            "head": git_value("rev-parse", "HEAD"),
            "branch": git_value("branch", "--show-current"),
            "describe": git_value("describe", "--always", "--dirty"),
        },
        "host": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "python": sys.version,
            "cpu_count": os.cpu_count(),
        },
        "environment": {
            key: os.environ.get(key)
            for key in (
                "OMP_NUM_THREADS",
                "KOKKOS_NUM_THREADS",
                "CUDA_VISIBLE_DEVICES",
                "SLURM_JOB_ID",
                "SLURM_NTASKS",
            )
            if os.environ.get(key) is not None
        },
        "start_utc": utc_now(),
        "status": "dry_run" if dry_run else "running",
    }
    metadata_path = run_dir / "run_metadata.json"

    if dry_run:
        metadata["end_utc"] = utc_now()
        metadata["elapsed_seconds"] = 0.0
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        print(metadata["command_shell"])
        return metadata

    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    log_path = run_dir / "athena.log"
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.run(
            command,
            cwd=REPO_ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed = time.perf_counter() - start
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    detected_failure = FAILURE_PATTERN.search(log_text)

    metadata.update(
        {
            "end_utc": utc_now(),
            "elapsed_seconds": elapsed,
            "return_code": process.returncode,
            "child_user_seconds": usage_after.ru_utime - usage_before.ru_utime,
            "child_system_seconds": usage_after.ru_stime - usage_before.ru_stime,
            "child_max_rss": None,
            "child_max_rss_note": (
                "Not reported: RUSAGE_CHILDREN peak RSS is cumulative across prior "
                "children and cannot be attributed to one paired lane."
            ),
            "log": str(log_path),
            "failure_text": detected_failure.group(0) if detected_failure else None,
            "status": (
                "passed"
                if process.returncode == 0 and detected_failure is None
                else "failed"
            ),
        }
    )
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def inspect_final_mesh_outputs(run_dir: Path) -> dict[str, Any]:
    gas_files = discover_output(run_dir, "bin", "gas", "bin")
    tracer_files = discover_output(run_dir, "bin", "tracer_ngp", "bin")
    if not gas_files or not tracer_files:
        raise FileNotFoundError("missing gas or tracer binary output")
    gas = BIN_CONVERT.read_binary_as_athdf(str(gas_files[-1]), quantities=["dens"])
    tracer = BIN_CONVERT.read_binary_as_athdf(str(tracer_files[-1]), quantities=["pdens"])
    density = np.asarray(gas["dens"], dtype=np.float64)
    particle_density = np.asarray(tracer["pdens"], dtype=np.float64)
    if not np.all(np.isfinite(density)) or np.min(density) <= 0.0:
        raise ValueError("gas density is non-finite or non-positive")
    if not np.all(np.isfinite(particle_density)) or np.min(particle_density) < 0.0:
        raise ValueError("particle density is non-finite or negative")
    if not np.isclose(gas["Time"], tracer["Time"], rtol=0.0, atol=1.0e-12):
        raise ValueError("gas and particle outputs are at different times")
    return {
        "gas": density,
        "tracer": particle_density,
        "time": float(gas["Time"]),
        "cycle": int(gas["NumCycles"]),
        "gas_file": str(gas_files[-1]),
        "tracer_file": str(tracer_files[-1]),
    }


def compare_gas(old: np.ndarray, corrected: np.ndarray) -> dict[str, float]:
    difference = corrected - old
    denominator = float(np.sqrt(np.mean(old**2)))
    return {
        "relative_l2": float(np.sqrt(np.mean(difference**2)) / denominator),
        "max_abs": float(np.max(np.abs(difference))),
    }


def preflight(args: argparse.Namespace) -> int:
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    cfl_values = [float(value) for value in args.cfl_values.split(",")]
    report: dict[str, Any] = {
        "formal_rk2_cfl_limit": FORMAL_RK2_CFL_LIMIT,
        "initial_3d_cfl_guard": INITIAL_3D_CFL_GUARD,
        "gas_match_tolerance": args.gas_match_tolerance,
        "runs": [],
        "selected_cfl": None,
        "source_state": capture_repo_state(output_root),
    }

    for cfl in cfl_values:
        cfl_name = f"cfl_{cfl:.6f}".replace(".", "p")
        cfl_result: dict[str, Any] = {"cfl": cfl, "old": {}, "corrected": {}}
        overrides = [f"time/cfl_number={cfl:.16g}", *args.override]
        for label, executable in (
            ("old", args.old_exe),
            ("corrected", args.corrected_exe),
        ):
            run_dir = output_root / cfl_name / label
            metadata = run_one(
                label,
                executable,
                args.input,
                run_dir,
                args.launcher,
                overrides,
                args.dry_run,
            )
            cfl_result[label]["runtime"] = metadata
            if not args.dry_run and metadata["status"] == "passed":
                try:
                    health = inspect_final_mesh_outputs(run_dir)
                    cfl_result[label]["health"] = {
                        key: value
                        for key, value in health.items()
                        if key not in ("gas", "tracer")
                    }
                    cfl_result[label]["_gas"] = health["gas"]
                except Exception as error:
                    cfl_result[label]["health_error"] = str(error)

        if not args.dry_run:
            old_gas = cfl_result["old"].pop("_gas", None)
            corrected_gas = cfl_result["corrected"].pop("_gas", None)
            if old_gas is not None and corrected_gas is not None:
                gas_comparison = compare_gas(old_gas, corrected_gas)
                cfl_result["gas_comparison"] = gas_comparison
                cfl_result["passed"] = (
                    gas_comparison["relative_l2"] <= args.gas_match_tolerance
                )
            else:
                cfl_result["passed"] = False
        else:
            cfl_result["passed"] = None
        report["runs"].append(cfl_result)

    eligible = [
        result["cfl"]
        for result in report["runs"]
        if result["passed"]
        and (
            result["cfl"] <= INITIAL_3D_CFL_GUARD + 1.0e-12
            or args.allow_above_initial_guard
        )
    ]
    if eligible:
        report["selected_cfl"] = max(eligible)
    report_path = output_root / "preflight_summary.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for result in report["runs"]:
        old_status = result["old"]["runtime"]["status"]
        corrected_status = result["corrected"]["runtime"]["status"]
        gas_error = result.get("gas_comparison", {}).get("relative_l2", math.nan)
        print(
            f"CFL={result['cfl']:.12g} old={old_status} corrected={corrected_status} "
            f"gas_rel_l2={gas_error:.3e} accepted={result['passed']}"
        )
    print(f"selected_cfl={report['selected_cfl']}")
    print(f"report={report_path}")
    if args.dry_run:
        return 0
    return 0 if report["selected_cfl"] is not None else 1


def run_pair(args: argparse.Namespace) -> int:
    if (
        args.cfl > INITIAL_3D_CFL_GUARD + 1.0e-12
        and not args.allow_above_initial_guard
    ):
        raise SystemExit(
            f"CFL {args.cfl} exceeds the initial three-dimensional guard "
            f"{INITIAL_3D_CFL_GUARD:.12g}; run a full-duration qualification and pass "
            "--allow-above-initial-guard only with explicit justification."
        )

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    overrides = [f"time/cfl_number={args.cfl:.16g}", *args.override]
    runs: dict[str, Any] = {}
    for label, executable in (("old", args.old_exe), ("corrected", args.corrected_exe)):
        run_dir = output_root / label
        runs[label] = run_one(
            label,
            executable,
            args.input,
            run_dir,
            args.launcher,
            overrides,
            args.dry_run,
        )

    manifest = {
        "input": str(args.input.resolve()),
        "cfl": args.cfl,
        "old_run": str((output_root / "old").resolve()),
        "corrected_run": str((output_root / "corrected").resolve()),
        "runs": runs,
        "source_state": capture_repo_state(output_root),
    }
    (output_root / "pair_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if args.dry_run:
        return 0
    if any(metadata["status"] != "passed" for metadata in runs.values()):
        return 1

    old_health = inspect_final_mesh_outputs(output_root / "old")
    corrected_health = inspect_final_mesh_outputs(output_root / "corrected")
    gas_comparison = compare_gas(old_health["gas"], corrected_health["gas"])
    manifest["gas_comparison"] = gas_comparison
    (output_root / "pair_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if gas_comparison["relative_l2"] > args.gas_match_tolerance:
        print(
            f"paired gas fields differ: relative L2={gas_comparison['relative_l2']:.3e}",
            file=sys.stderr,
        )
        return 1

    if not args.skip_analysis:
        analysis_command = [
            sys.executable,
            str(REPO_ROOT / "scripts/analyze_ito_turbulence_validation.py"),
            "--old-run",
            str(output_root / "old"),
            "--corrected-run",
            str(output_root / "corrected"),
            "--output-dir",
            str(output_root / "analysis"),
            "--deposition",
            args.deposition,
        ]
        print(shlex.join(analysis_command), flush=True)
        result = subprocess.run(analysis_command, cwd=REPO_ROOT, check=False)
        if result.returncode != 0:
            return result.returncode
    return 0


def add_common_arguments(parser: argparse.ArgumentParser):
    parser.add_argument("--old-exe", type=Path, required=True)
    parser.add_argument("--corrected-exe", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--launcher",
        default="",
        help='Optional launcher prefix, for example "mpirun -np 8".',
    )
    parser.add_argument("--gas-match-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--allow-above-initial-guard", action="store_true")
    parser.add_argument(
        "--override",
        action="append",
        default=[],
        help="Additional AthenaK block/name=value override; may be repeated.",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser("preflight")
    add_common_arguments(preflight_parser)
    preflight_parser.add_argument(
        "--cfl-values",
        default="0.25,0.30,0.32,0.329,0.33,0.331,0.34",
    )
    preflight_parser.set_defaults(function=preflight)

    pair_parser = subparsers.add_parser("run-pair")
    add_common_arguments(pair_parser)
    pair_parser.add_argument("--cfl", type=float, default=INITIAL_3D_CFL_GUARD)
    pair_parser.add_argument("--skip-analysis", action="store_true")
    pair_parser.add_argument(
        "--deposition", choices=("auto", "cic", "ngp"), default="auto"
    )
    pair_parser.set_defaults(function=run_pair)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for attribute in ("old_exe", "corrected_exe", "input"):
        path = getattr(args, attribute)
        if not path.expanduser().exists():
            raise SystemExit(f"{attribute.replace('_', '-')} does not exist: {path}")
        setattr(args, attribute, path.expanduser().resolve())
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
