#!/usr/bin/env python3
"""Create and run tiered CR-shock scan matrices with machine-readable summaries."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import shlex
import subprocess
import time


@dataclass
class ScanCase:
    case_id: str
    mach: float
    ppc: float
    nx: int
    nx3: int
    mb1: int
    mb3: int
    command: list[str]
    status: str = "pending"
    return_code: int | None = None
    runtime_s: float = 0.0


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _parse_float_list(raw: str) -> list[float]:
    out = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        out.append(float(item))
    if not out:
        raise ValueError("list cannot be empty")
    return out


def _parse_int_list(raw: str) -> list[int]:
    return [int(v) for v in _parse_float_list(raw)]


def _resolve_default_lists(tier: str) -> tuple[str, str, str]:
    if tier == "hpc":
        return ("1.0,1.2,1.4", "8,16,24", "256,384,512")
    return ("0.8,1.0,1.2", "2,4,8", "32,48,64")


def _resolve_nproc(tier: str, nproc: int | None) -> int:
    if nproc is not None:
        return int(nproc)
    if tier == "hpc":
        return 4096
    return 8


def _resolve_deck(repo_root: Path, tier: str, deck_arg: str) -> Path:
    if deck_arg:
        if deck_arg in {"local", "hpc"}:
            tag = deck_arg
            return repo_root / f"inputs/tests/pic_amr_shock_lb_publication_{tag}.athinput"
        return (repo_root / deck_arg).resolve()
    return repo_root / f"inputs/tests/pic_amr_shock_lb_publication_{tier}.athinput"


def _command_str(parts: list[str]) -> str:
    return " ".join(shlex.quote(p) for p in parts)


def _launch_prefix(
    template: str, mpiexec: str, nproc: int, gpus_per_rank: int
) -> list[str]:
    if template:
        rendered = template.format(nproc=nproc, gpus_per_rank=gpus_per_rank)
        return shlex.split(rendered)
    return [mpiexec, "-n", str(nproc)] if nproc > 1 else []


def _build_cases(
    tier: str,
    deck: Path,
    mach_values: list[float],
    ppc_values: list[float],
    res_values: list[int],
    launch_prefix: list[str],
    extra_overrides: list[str],
) -> list[ScanCase]:
    cases: list[ScanCase] = []
    for mach in mach_values:
        for ppc in ppc_values:
            for nx in res_values:
                nx3 = max(8, nx // 4)
                mb1 = max(4, nx // 8)
                mb3 = max(4, nx3 // 2)
                case_id = (
                    f"pic_shock_scan_{tier}_m{mach:.2f}_ppc{ppc:.2f}_n{nx}"
                    .replace(".", "p")
                    .replace("-", "m")
                )

                overrides = [
                    f"job/basename={case_id}",
                    f"problem/ot_mach={mach:.6f}",
                    f"particles/ppc={ppc:.6f}",
                    f"mesh/nx1={nx}",
                    f"mesh/nx2={nx}",
                    f"mesh/nx3={nx3}",
                    f"meshblock/nx1={mb1}",
                    f"meshblock/nx2={mb1}",
                    f"meshblock/nx3={mb3}",
                ]
                overrides.extend(extra_overrides)
                cmd = ["./athena", "-i", str(deck)] + overrides
                full_cmd = launch_prefix + cmd
                cases.append(
                    ScanCase(
                        case_id=case_id,
                        mach=mach,
                        ppc=ppc,
                        nx=nx,
                        nx3=nx3,
                        mb1=mb1,
                        mb3=mb3,
                        command=full_cmd,
                    )
                )
    return cases


def _write_command_script(path: Path, cwd: Path, cases: list[ScanCase]) -> None:
    with path.open("w", encoding="utf-8") as fobj:
        fobj.write("#!/usr/bin/env bash\nset -euo pipefail\n\n")
        fobj.write(f"cd {shlex.quote(str(cwd))}\n\n")
        for case in cases:
            fobj.write(_command_str(case.command) + "\n")
    path.chmod(0o755)


def _write_summary_csv(path: Path, cases: list[ScanCase]) -> None:
    rows = []
    for case in cases:
        rows.append(
            {
                "case_id": case.case_id,
                "mach": case.mach,
                "ppc": case.ppc,
                "nx": case.nx,
                "nx3": case.nx3,
                "mb1": case.mb1,
                "mb3": case.mb3,
                "status": case.status,
                "return_code": case.return_code,
                "runtime_s": f"{case.runtime_s:.3f}",
                "command": _command_str(case.command),
            }
        )

    with path.open("w", encoding="utf-8", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create and run CR-shock scan matrices for local/HPC tiers."
    )
    parser.add_argument("--tier", choices=("local", "hpc"), default="local")
    parser.add_argument(
        "--deck",
        default="",
        help="Optional deck selector: local|hpc|<path>. Defaults to tier deck.",
    )
    parser.add_argument("--mach-list", default="")
    parser.add_argument("--ppc-list", default="")
    parser.add_argument("--res-list", default="")
    parser.add_argument(
        "--nproc",
        type=int,
        default=None,
        help="MPI ranks to launch (default: 8 local, 4096 hpc).",
    )
    parser.add_argument(
        "--launch-template",
        default="",
        help=(
            "Optional launcher prefix template with {nproc} and {gpus_per_rank}, "
            "for example: 'srun -n {nproc} --gpus-per-task={gpus_per_rank}'."
        ),
    )
    parser.add_argument("--gpus-per-rank", type=int, default=0)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument(
        "--extra-override",
        action="append",
        default=[],
        help="Additional Athena override (block/param=value). Repeatable.",
    )
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument(
        "--output-dir",
        default="tst/.codex/pic_publication_runs",
        help="Directory for command scripts and scan summaries.",
    )
    parser.add_argument("--run-id", default="")
    parser.add_argument("--athena-cwd", default="tst/build/src")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    run_id = args.run_id if args.run_id else f"pic_shock_scan_{_timestamp()}"
    output_root = (repo_root / args.output_dir).resolve()
    scan_dir = output_root / run_id / "reviews" / "shock_step4"
    scan_dir.mkdir(parents=True, exist_ok=True)

    default_mach, default_ppc, default_res = _resolve_default_lists(args.tier)
    mach_values = _parse_float_list(args.mach_list or default_mach)
    ppc_values = _parse_float_list(args.ppc_list or default_ppc)
    res_values = _parse_int_list(args.res_list or default_res)

    if len(mach_values) < 3 or len(ppc_values) < 3 or len(res_values) < 2:
        raise ValueError(
            "scan matrix must include >=3 Mach, >=3 CR loading, >=2 resolution"
        )

    nproc = _resolve_nproc(args.tier, args.nproc)
    deck = _resolve_deck(repo_root, args.tier, args.deck)
    if not deck.exists():
        raise FileNotFoundError(f"deck not found: {deck}")

    launch_prefix = _launch_prefix(
        args.launch_template, args.mpiexec, nproc, args.gpus_per_rank
    )
    cases = _build_cases(
        args.tier,
        deck,
        mach_values,
        ppc_values,
        res_values,
        launch_prefix,
        args.extra_override,
    )
    if args.max_cases > 0:
        cases = cases[: args.max_cases]

    athena_cwd = (repo_root / args.athena_cwd).resolve()
    script_path = scan_dir / "scan_commands.sh"
    summary_json = scan_dir / "scan_summary.json"
    summary_csv = scan_dir / "scan_summary.csv"

    _write_command_script(script_path, athena_cwd, cases)

    if args.run:
        for case in cases:
            t0 = time.time()
            proc = subprocess.run(case.command, cwd=str(athena_cwd), check=False)
            case.runtime_s = float(time.time() - t0)
            case.return_code = int(proc.returncode)
            case.status = "passed" if proc.returncode == 0 else "failed"
            if proc.returncode != 0 and not args.keep_going:
                break
    else:
        for case in cases:
            case.status = "pending"
            case.return_code = None

    cases_dict = [asdict(case) for case in cases]
    summary = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "run_id": run_id,
        "tier": args.tier,
        "deck": str(deck),
        "athena_cwd": str(athena_cwd),
        "nproc": nproc,
        "launch_prefix": launch_prefix,
        "matrix": {
            "mach_values": mach_values,
            "ppc_values": ppc_values,
            "res_values": res_values,
            "n_cases": len(cases),
        },
        "extra_overrides": args.extra_override,
        "cases": cases_dict,
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    if cases:
        _write_summary_csv(summary_csv, cases)

    print("scan_dir:", scan_dir)
    print("command_script:", script_path)
    print("summary_json:", summary_json)
    print("summary_csv:", summary_csv)
    print("n_cases:", len(cases))

    if not args.run:
        print("mode: emit-only (use --run to execute)")

    failed = [case for case in cases if case.status == "failed"]
    return 0 if not failed else 2


if __name__ == "__main__":
    raise SystemExit(main())
