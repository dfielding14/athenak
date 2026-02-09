#!/usr/bin/env python3
"""Run publication-priority PIC/MHD-PIC cases and archive artifacts."""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Sequence

from pic_publication_manifest import PublicationCase
from pic_publication_manifest import select_cases


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _run_command(
    cmd: Sequence[str],
    cwd: Path,
    log_path: Path,
    dry_run: bool,
    env: dict[str, str] | None = None,
) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("cwd: " + str(cwd) + "\n")
        log.write("cmd: " + " ".join(cmd) + "\n\n")
        if dry_run:
            log.write("dry-run: command not executed\n")
            return 0
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            env=env,
        )
        return int(proc.returncode)


def _cleanup_patterns(repo_root: Path, patterns: Sequence[str]) -> None:
    for pattern in patterns:
        for path in repo_root.glob(pattern):
            if path.is_file():
                path.unlink()


def _archive_outputs(
    repo_root: Path,
    case_dir: Path,
    patterns: Sequence[str],
) -> list[str]:
    files_dir = case_dir / "files"
    copied: list[str] = []
    for pattern in patterns:
        for path in sorted(repo_root.glob(pattern)):
            if not path.is_file():
                continue
            rel = path.relative_to(repo_root)
            dst = files_dir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dst)
            copied.append(str(rel))
    return copied


def _module_command(module_name: str, python_exec: str) -> list[str]:
    code = (
        "import importlib\n"
        f"m = importlib.import_module('{module_name}')\n"
        "m.run()\n"
        "ok = m.analyze()\n"
        "import sys\n"
        "sys.exit(0 if ok else 2)\n"
    )
    return [python_exec, "-c", code]


def _deck_command(
    repo_root: Path,
    case: PublicationCase,
    mpiexec: str,
    nproc: int,
) -> tuple[list[str], Path]:
    athena_cwd = repo_root / "tst" / "build" / "src"
    input_path = repo_root / case.input_deck
    cmd = ["./athena", "-i", str(input_path)] + list(case.athena_args)
    if nproc > 1:
        cmd = [mpiexec, "-n", str(nproc)] + cmd
    return cmd, athena_cwd


def _case_nproc(case: PublicationCase, tier: str) -> int:
    if tier == "hpc":
        return case.ranks_hpc
    if tier == "local":
        return case.ranks_local
    return case.ranks_hpc if case.profile == "hpc" else case.ranks_local


def _build_athena(repo_root: Path, cmake_args: Sequence[str], dry_run: bool) -> int:
    build_log = repo_root / "tst" / ".codex" / "pic_publication_build.log"
    cfg_cmd = ["cmake", "-S", "..", "-B", "build"] + list(cmake_args)
    rc = _run_command(cfg_cmd, repo_root / "tst", build_log, dry_run)
    if rc != 0:
        return rc
    bld_cmd = ["cmake", "--build", "build", "-j8"]
    return _run_command(bld_cmd, repo_root / "tst", build_log, dry_run)


def _prepare_tst_build_dir(repo_root: Path, dry_run: bool) -> None:
    """Ensure `tst/build` is a real build directory, not a symlink."""
    tst_build = repo_root / "tst" / "build"
    if tst_build.is_symlink():
        if dry_run:
            return
        tst_build.unlink()


def _athena_exists(repo_root: Path) -> bool:
    return (repo_root / "tst" / "build" / "src" / "athena").exists()


def _athena_config(repo_root: Path, dry_run: bool) -> dict[str, object]:
    if dry_run or not _athena_exists(repo_root):
        return {"available": False, "mpi_enabled": False, "raw": ""}

    proc = subprocess.run(
        ["./athena", "-c"],
        cwd=str(repo_root / "tst" / "build" / "src"),
        capture_output=True,
        text=True,
        check=False,
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    mpi_on = "MPI parallelism:            ON" in output
    return {
        "available": proc.returncode == 0,
        "mpi_enabled": mpi_on,
        "raw": output,
    }


def _parse_csv(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run publication-priority PIC/MHD-PIC suite and archive artifacts."
    )
    parser.add_argument(
        "--groups",
        default=(
            "entity_core_proxy,benchmarks_proxy,"
            "entity_core_publication,benchmarks_publication,shock_story"
        ),
        help="Comma-separated case groups.",
    )
    parser.add_argument(
        "--tier",
        choices=("local", "hpc", "all"),
        default="local",
        help="Execution tier filter.",
    )
    parser.add_argument(
        "--cases",
        default="",
        help="Optional comma-separated case_id filter.",
    )
    parser.add_argument(
        "--run-root",
        default="tst/.codex/pic_publication_runs",
        help="Artifact root directory.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Run identifier (default: timestamp).",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip configure/build step.",
    )
    parser.add_argument(
        "--cmake-arg",
        action="append",
        default=[],
        help="Additional cmake configure arg; repeat for multiple args.",
    )
    parser.add_argument(
        "--enable-mpi",
        action="store_true",
        help="Append -DAthena_ENABLE_MPI=ON during configure.",
    )
    parser.add_argument(
        "--disable-mpi",
        action="store_true",
        help="Append -DAthena_ENABLE_MPI=OFF during configure.",
    )
    parser.add_argument(
        "--mpiexec",
        default=os.environ.get("MPIEXEC", "mpiexec"),
        help="MPI launcher for nproc > 1.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used for module-run cases.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List selected cases and exit.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue after failed cases.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print and log commands without executing.",
    )
    args = parser.parse_args()
    if args.enable_mpi and args.disable_mpi:
        parser.error("Use at most one of --enable-mpi or --disable-mpi.")

    repo_root = Path(__file__).resolve().parents[2]
    run_id = args.run_id if args.run_id else _timestamp()
    run_root = repo_root / args.run_root
    run_dir = run_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    groups = _parse_csv(args.groups)
    selected = select_cases(groups, args.tier)

    case_filter = set(_parse_csv(args.cases))
    if case_filter:
        selected = [case for case in selected if case.case_id in case_filter]

    if args.list:
        for case in selected:
            print(f"{case.case_id:36s} group={case.group:12s} runner={case.runner}")
        return 0

    summary = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "repo_root": str(repo_root),
        "tier": args.tier,
        "groups": groups,
        "cases": [],
    }

    _prepare_tst_build_dir(repo_root, args.dry_run)
    cmake_args = list(args.cmake_arg)
    if args.enable_mpi:
        cmake_args.append("-DAthena_ENABLE_MPI=ON")
    if args.disable_mpi:
        cmake_args.append("-DAthena_ENABLE_MPI=OFF")

    if not args.skip_build:
        rc = _build_athena(repo_root, cmake_args, args.dry_run)
        if rc != 0:
            summary["build_status"] = "failed"
            summary["build_return_code"] = rc
            (run_dir / "summary.json").write_text(
                json.dumps(summary, indent=2), encoding="utf-8"
            )
            return rc
    elif not args.dry_run and not _athena_exists(repo_root):
        summary["build_status"] = "missing"
        summary["error"] = (
            "skip-build requested but tst/build/src/athena is missing; "
            "rerun without --skip-build"
        )
        (run_dir / "summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )
        return 3
    summary["build_status"] = "skipped" if args.skip_build else "ok"
    summary["athena_config"] = _athena_config(repo_root, args.dry_run)

    overall_ok = True
    for case in selected:
        case_dir = run_dir / "cases" / case.case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        case_log = case_dir / "run.log"

        _cleanup_patterns(repo_root, case.output_globs)

        if case.runner == "module":
            cmd = _module_command(case.module, args.python)
            cwd = repo_root / "tst"
            nproc = None
            case_env = dict(os.environ)
            case_env["MPIEXEC"] = args.mpiexec
        elif case.runner == "deck":
            nproc = _case_nproc(case, args.tier)
            cmd, cwd = _deck_command(repo_root, case, args.mpiexec, nproc)
            case_env = None
        else:
            raise RuntimeError(f"Unsupported runner '{case.runner}'")

        rc = _run_command(cmd, cwd, case_log, args.dry_run, env=case_env)
        copied = _archive_outputs(repo_root, case_dir, case.output_globs)

        case_record = {
            "case_id": case.case_id,
            "group": case.group,
            "runner": case.runner,
            "profile": case.profile,
            "note": case.note,
            "cwd": str(cwd),
            "command": cmd,
            "nproc": nproc,
            "return_code": rc,
            "status": "ok" if rc == 0 else "failed",
            "copied_outputs": copied,
            "log": str(case_log.relative_to(repo_root)),
        }
        summary["cases"].append(case_record)
        (case_dir / "case.json").write_text(
            json.dumps(case_record, indent=2), encoding="utf-8"
        )

        if rc != 0:
            overall_ok = False
            if not args.keep_going:
                break

    summary["overall_status"] = "ok" if overall_ok else "failed"
    (run_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print("run_dir:", run_dir)
    print("overall_status:", summary["overall_status"])
    return 0 if overall_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
