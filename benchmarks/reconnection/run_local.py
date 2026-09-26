#!/usr/bin/env python3
"""Run one reproducible Harris-sheet case using at most four single-threaded MPI ranks."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
THREAD_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="new case directory")
    parser.add_argument("--cells-per-di", type=int, choices=(2, 4, 8))
    parser.add_argument("--model", choices=("current_limited", "constant", "ideal"))
    parser.add_argument("--tlim", type=positive, required=True)
    parser.add_argument("--cfl", type=positive, default=0.4)
    parser.add_argument("--sts-ratio", type=positive, default=32)
    parser.add_argument("--ranks", type=int, choices=range(1, 5), default=3)
    parser.add_argument("--output-dt", type=positive, default=0.1)
    parser.add_argument("--history-dt", type=positive, default=0.01)
    parser.add_argument("--restart-dt", type=positive, default=1)
    parser.add_argument("--nlim", type=int, default=-1, help="absolute cycle limit for pilots")
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--wall-limit", help="Athena wall limit HH:MM:SS; writes final outputs")
    args = parser.parse_args()
    if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
        parser.error("--binary must name an executable")
    if args.nlim < -1 or args.cfl > 1:
        parser.error("require nlim >= -1 and CFL <= 1")
    if args.wall_limit and not re.fullmatch(r"\d+:[0-5]\d:[0-5]\d", args.wall_limit):
        parser.error("--wall-limit must use HH:MM:SS")
    if args.output.exists():
        parser.error("--output already exists; choose a new directory")
    if args.restart and not args.restart.is_file():
        parser.error("restart checkpoint does not exist")
    if not args.restart and args.cells_per_di is None:
        parser.error("new runs require --cells-per-di")

    binary, output = args.binary.resolve(), args.output.resolve()
    model = args.model or "current_limited"
    physics = {"inherited_from": str(args.restart.resolve())} if args.restart else {
        "model": model, "cells_per_di": args.cells_per_di, "nx1": 200*args.cells_per_di,
        "nx2": 100*args.cells_per_di, "meshblock": [40, 40],
        "ohmic_resistivity": {"current_limited": 1e-6, "constant": 5e-6, "ideal": 0}[model],
        "resistivity_integrator": "sts" if model == "current_limited" else "explicit"}
    text = "# Restart overrides; mesh and physics are inherited.\n"
    settings = {"time/tlim": args.tlim, "time/nlim": args.nlim,
                "time/cfl_number": args.cfl, "time/sts_max_dt_ratio": args.sts_ratio}
    if args.restart:
        for parent in args.restart.resolve().parents:
            source = parent / "run.json"
            if source.is_file():
                saved = json.loads(source.read_text())["physics"]
                for key in ("model", "cells_per_di"):
                    requested = getattr(args, key)
                    if requested is not None and key in saved and requested != saved[key]:
                        parser.error(f"restart {key} differs from {source}")
                physics = dict(saved, inherited_from=str(args.restart.resolve()))
                break
    else:
        text = (ROOT / "inputs/mhd/resistive_harris.athinput").read_text()
        settings.update({"mesh/nx1": physics["nx1"], "mesh/nx2": physics["nx2"],
                         "meshblock/nx1": 40, "meshblock/nx2": 40,
                         "mhd/resistivity_model": "current_limited" if model == "current_limited"
                         else "constant", "mhd/ohmic_resistivity": physics["ohmic_resistivity"],
                         "mhd/resistivity_integrator": physics["resistivity_integrator"],
                         "time/sts_integrator": "rkl2" if model == "current_limited" else "none"})
        if model != "current_limited":
            settings.update({"output5/variable": "mhd_bmag", "output5/id": "bmag"})
    for index in (1, 2, 4, 5):
        settings[f"output{index}/dt"] = args.output_dt
    settings.update({"output3/dt": args.history_dt, "output6/file_type": "rst",
                     "output6/dt": args.restart_dt})
    for key, value in settings.items():
        block, parameter = key.split("/")
        text += f"\n<{block}>\n{parameter} = {value}\n"
    output.mkdir(parents=True)
    (output / "input.athinput").write_text(text)
    command = shlex.split(os.environ.get("MPIEXEC", "mpirun"))
    command += ["-np", str(args.ranks), str(binary), "-i", str(output / "input.athinput")]
    if args.restart:
        command += ["-r", str(args.restart.resolve())]
    if args.wall_limit:
        command += ["-t", args.wall_limit]
    environment = dict(os.environ, **{key: "1" for key in THREAD_KEYS})
    manifest = {"status": "running", "physics": physics, "command": command,
                "arguments": {key: str(value) if isinstance(value, Path) else value
                              for key, value in vars(args).items()},
                "thread_environment": {key: environment[key] for key in THREAD_KEYS},
                "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                                    text=True).strip(),
                "git_status": subprocess.check_output(["git", "status", "--porcelain"],
                                                      cwd=ROOT, text=True),
                "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                "started_unix": time.time()}
    destination = output / "run.json"
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    started = time.monotonic()
    try:
        with (output / "stdout.log").open("w") as log:
            result = subprocess.run(command, cwd=output, env=environment, stdout=log,
                                    stderr=subprocess.STDOUT, check=False)
        manifest["returncode"] = result.returncode
    except (OSError, KeyboardInterrupt) as error:
        manifest.update(returncode=1, error=f"{type(error).__name__}: {error}")
    manifest["wall_seconds"] = time.monotonic() - started
    log = (output / "stdout.log").read_text(errors="replace")
    final = re.findall(r"^time=(\S+) cycle=(\d+)", log, re.MULTILINE)
    termination = re.findall(r"Terminating on (.+)", log)
    counts = re.search(r"STS sweeps = (\d+) STS stages = (\d+)", log)
    manifest.update(final_time=float(final[-1][0]) if final else None,
                    final_cycle=int(final[-1][1]) if final else None,
                    termination=termination[-1] if termination else None,
                    sts_sweeps=int(counts[1]) if counts else 0,
                    sts_stages=int(counts[2]) if counts else 0)
    manifest["reached_tlim"] = bool(
        final and (manifest["termination"] == "time limit"
                   or manifest["final_time"] >= args.tlim))
    manifest["status"] = ("failed" if manifest["returncode"] or not final else
                          "complete" if manifest["reached_tlim"] else "partial")
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{manifest['status']}: {output}")
    return 1 if manifest["status"] == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
