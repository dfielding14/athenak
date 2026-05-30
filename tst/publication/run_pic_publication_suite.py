#!/usr/bin/env python3
"""Run exploratory PIC/MHD-PIC engineering cases and archive artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Sequence

from pic_publication_manifest import ArtifactCase
from pic_publication_manifest import LEGACY_GROUP_ALIASES
from pic_publication_manifest import canonicalize_groups
from pic_publication_manifest import select_cases
from artifact_lineage import write_companion_manifest


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _frontier_runtime_detected(path: Path) -> bool:
    host = os.uname().nodename.lower()
    return (
        str(path).startswith("/lustre/orion/")
        or "frontier" in host
        or os.environ.get("LMOD_SYSTEM_NAME", "").lower() == "frontier"
    )


def _in_slurm_environment() -> bool:
    return bool(os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fobj:
        for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
) -> list[dict[str, str]]:
    files_dir = case_dir / "files"
    copied: list[dict[str, str]] = []
    for pattern in patterns:
        for path in sorted(repo_root.glob(pattern)):
            if not path.is_file():
                continue
            rel = path.relative_to(repo_root)
            dst = files_dir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dst)
            copied.append(
                {
                    "source_path": str(rel),
                    "archive_path": str(dst.relative_to(case_dir)),
                    "sha256": _sha256(dst),
                }
            )
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
    case: ArtifactCase,
    mpiexec: str,
    nproc: int,
) -> tuple[list[str], Path]:
    athena_cwd = repo_root / "tst" / "build" / "src"
    input_path = repo_root / case.input_deck
    cmd = ["./athena", "-i", str(input_path)] + list(case.athena_args)
    if nproc > 1:
        cmd = [mpiexec, "-n", str(nproc)] + cmd
    return cmd, athena_cwd


def _case_nproc(case: ArtifactCase, tier: str) -> int:
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


def _display_path(path: Path, repo_root: Path) -> str:
    try:
        return str(path.relative_to(repo_root))
    except ValueError:
        return str(path)


def _git_output(repo_root: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        check=False,
    )
    return (proc.stdout or "") + (proc.stderr or "")


def _file_record(path: Path, root: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(root)),
        "sha256": _sha256(path),
    }


def _write_artifact_manifest(run_dir: Path, repo_root: Path) -> None:
    manifest_path = run_dir / "artifact_manifest.json"
    files = [
        _file_record(path, run_dir)
        for path in sorted(run_dir.rglob("*"))
        if path.is_file() and path != manifest_path
    ]
    athena = repo_root / "tst" / "build" / "src" / "athena"
    manifest = {
        "schema_version": 1,
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_exploratory_artifacts",
        "git_commit": _git_output(repo_root, "rev-parse", "HEAD").strip(),
        "git_status": _git_output(
            repo_root, "status", "--short", "--branch"
        ).splitlines(),
        "executable": str(athena) if athena.is_file() else "",
        "executable_sha256": _sha256(athena) if athena.is_file() else "",
        "files": files,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    lineage_path = run_dir / "lineage_manifest.json"
    lineage_outputs = [
        path
        for path in sorted(run_dir.rglob("*"))
        if path.is_file() and path != lineage_path
    ]
    write_companion_manifest(
        lineage_path,
        generator="run_pic_publication_suite.py",
        physical_model="exploratory_pic_artifact_bundle",
        inputs=(),
        outputs=lineage_outputs,
        executables=(athena,) if athena.is_file() else (),
    )


def _write_summary(run_dir: Path, repo_root: Path, summary: dict[str, object]) -> None:
    (run_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    _write_artifact_manifest(run_dir, repo_root)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run exploratory PIC/MHD-PIC engineering suite and archive artifacts."
    )
    parser.add_argument(
        "--groups",
        default=(
            "entity_core_engineering,benchmark_engineering,"
            "extended_entity_engineering,extended_benchmark_engineering,shock_engineering"
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
    parser.add_argument(
        "--local-exploratory-run",
        action="store_true",
        help="Acknowledge that direct execution is limited to local engineering use.",
    )
    args = parser.parse_args()
    if args.enable_mpi and args.disable_mpi:
        parser.error("Use at most one of --enable-mpi or --disable-mpi.")

    repo_root = Path(__file__).resolve().parents[2]
    if not args.list and not args.dry_run:
        if args.tier != "local":
            parser.error(
                "direct HPC/all-tier execution is prohibited; use immutable "
                "snapshots and the approved external scheduler wrapper"
            )
        if not args.local_exploratory_run:
            parser.error("local execution requires --local-exploratory-run")
        if _in_slurm_environment() or _frontier_runtime_detected(repo_root):
            parser.error(
                "direct execution is local-only; use the approved external "
                "manifest, validator, ledger, and scheduler wrapper"
            )
    run_id = args.run_id if args.run_id else _timestamp()
    run_root = repo_root / args.run_root
    run_dir = run_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    groups = _parse_csv(args.groups)
    canonical_groups = sorted(canonicalize_groups(groups))
    selected = select_cases(groups, args.tier)

    case_filter = set(_parse_csv(args.cases))
    if case_filter:
        selected = [case for case in selected if case.case_id in case_filter]

    if args.list:
        for group in groups:
            if group in LEGACY_GROUP_ALIASES:
                print(
                    f"warning: legacy group '{group}' maps to "
                    f"'{LEGACY_GROUP_ALIASES[group]}'"
                )
        for case in selected:
            print(
                f"{case.case_id:36s} group={case.group:28s} "
                f"runner={case.runner:6s} evidence_class={case.evidence_class}"
            )
        return 0

    summary = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "repo_root": str(repo_root),
        "tier": args.tier,
        "requested_groups": groups,
        "canonical_groups": canonical_groups,
        "legacy_group_aliases_used": {
            group: LEGACY_GROUP_ALIASES[group]
            for group in groups
            if group in LEGACY_GROUP_ALIASES
        },
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_exploratory_artifacts",
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
            _write_summary(run_dir, repo_root, summary)
            return rc
    elif not args.dry_run and not _athena_exists(repo_root):
        summary["build_status"] = "missing"
        summary["error"] = (
            "skip-build requested but tst/build/src/athena is missing; "
            "rerun without --skip-build"
        )
        _write_summary(run_dir, repo_root, summary)
        return 3
    summary["build_status"] = "skipped" if args.skip_build else "ok"
    summary["athena_config"] = _athena_config(repo_root, args.dry_run)

    overall_ok = True
    for case in selected:
        case_dir = run_dir / "cases" / case.case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        case_log = case_dir / "run.log"

        if not args.dry_run:
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
        copied = (
            []
            if args.dry_run
            else _archive_outputs(repo_root, case_dir, case.output_globs)
        )
        module_source = (
            repo_root / "tst" / (case.module.replace(".", "/") + ".py")
            if case.module
            else None
        )

        case_record = {
            "case_id": case.case_id,
            "group": case.group,
            "runner": case.runner,
            "profile": case.profile,
            "evidence_class": case.evidence_class,
            "not_sun_bai_reproduction": True,
            "qualification_status": "unqualified_exploratory_artifact",
            "physical_model": case.physical_model,
            "legacy_identifier": case.legacy_identifier,
            "input_deck": case.input_deck,
            "input_deck_sha256": (
                _sha256(repo_root / case.input_deck) if case.input_deck else ""
            ),
            "module_source": "" if module_source is None else str(module_source),
            "module_source_sha256": (
                _sha256(module_source)
                if module_source is not None and module_source.is_file()
                else ""
            ),
            "note": case.note,
            "cwd": str(cwd),
            "command": cmd,
            "nproc": nproc,
            "return_code": rc,
            "status": "ok" if rc == 0 else "failed",
            "copied_outputs": copied,
            "log": _display_path(case_log, repo_root),
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
    _write_summary(run_dir, repo_root, summary)
    print("run_dir:", run_dir)
    print("overall_status:", summary["overall_status"])
    return 0 if overall_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
