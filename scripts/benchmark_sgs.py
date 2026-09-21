#!/usr/bin/env python3
"""Measure SGS output overhead on the target machine and filesystem.

Example (inside a GPU allocation):
  python3 scripts/benchmark_sgs.py build/src/athena case.athinput \
      --launcher 'srun -n 4' --output-dir /scratch/sgs-benchmark time/nlim=1000

All other outputs remain enabled. Wall time includes launch/setup/initial output;
Athena's timer includes evolution/final output. Logs, inputs, and files are retained.
Run long enough to produce several scheduled SGS snapshots beyond initial/final output.
Use the Kokkos profiling regions to separate SGS calculation from cbin file writing.
"""

import argparse
import json
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import time


def prepare_inputs(source, overrides):
    """Apply Athena overrides and remove every SGS block from the baseline copy."""
    source = re.split(r"(?m)^\s*<par_end>", source, maxsplit=1)[0]
    for override in overrides:
        match = re.fullmatch(r"([^/\s]+)/([^=\s]+)=([^\n]+)", override)
        if not match:
            raise ValueError(f"Expected block/parameter=value: {override}")
        block, key, value = match.groups()
        source += f"\n<{block}>\n{key} = {value}\n"
    parts = re.split(r"(?m)^\s*<([^>]+)>[^\n]*(?:\n|$)", source)
    blocks = {}
    for name, body in zip(parts[1::2], parts[2::2]):
        fields = re.findall(r"(?m)^\s*(\w+)\s*=\s*([^#\n]*)", body)
        blocks.setdefault(name, {}).update((key, value.strip()) for key, value in fields)
    sgs = {name for name, fields in blocks.items()
           if name.startswith("output") and fields.get("variable") == "hydro_sgs_2d"}
    if not sgs:
        raise ValueError("Input must contain hydro_sgs_2d output blocks")
    baseline = parts[0] + "".join(f"<{name}>\n{body}"
        for name, body in zip(parts[1::2], parts[2::2]) if name not in sgs)
    return {"off": baseline, "on": source}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("executable", type=Path)
    parser.add_argument("input", type=Path)
    parser.add_argument("overrides", nargs="*")
    parser.add_argument("--launcher", default="", help="e.g. 'srun -n 4'")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_intermixed_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    inputs = prepare_inputs(args.input.read_text(), args.overrides)
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    samples = []
    endpoint = None
    for repeat in range(args.repeats):
        for mode in (("off", "on") if repeat % 2 == 0 else ("on", "off")):
            directory = root / f"{repeat:02d}-{mode}"
            directory.mkdir()  # Refuse to overwrite an earlier benchmark.
            input_file = directory / "case.athinput"
            input_file.write_text(inputs[mode])
            command = [*shlex.split(args.launcher), str(args.executable.resolve()),
                       "-i", str(input_file), "-d", str(directory)]
            with (directory / "run.log").open("w") as log:
                start = time.perf_counter()
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
                wall = time.perf_counter() - start
            if result.returncode:
                raise RuntimeError(f"Run failed; see {directory / 'run.log'}")
            log = (directory / "run.log").read_text()
            final = re.findall(r"(?m)^time=(\S+) cycle=(\d+)$", log)[-1]
            if endpoint is not None and final != endpoint:
                raise RuntimeError("Runs reached different times/cycles; compare run logs")
            endpoint = final
            solver = float(re.findall(r"cpu time used\s*=\s*(\S+)", log)[-1])
            files = list(directory.rglob("*.cbin"))
            if mode == "on" and not files:
                raise RuntimeError("SGS-on run wrote no cbin files")
            sample = dict(mode=mode, repeat=repeat, wall_s=wall, solver_s=solver,
                          cbin_files=len(files),
                          cbin_bytes=sum(p.stat().st_size for p in files))
            samples.append(sample)
            (root / "results.json").write_text(json.dumps(samples, indent=2) + "\n")
            print(json.dumps(sample), flush=True)
    medians = {mode: {key: statistics.median(s[key] for s in samples if s["mode"] == mode)
                     for key in ("wall_s", "solver_s", "cbin_files", "cbin_bytes")}
               for mode in ("off", "on")}
    if medians["on"]["cbin_files"] <= medians["off"]["cbin_files"]:
        raise RuntimeError("SGS-on run produced no additional cbin files; check cadence")
    summary = {"medians": medians, "fractional_overhead": {
        key: medians["on"][key] / medians["off"][key] - 1
        for key in ("wall_s", "solver_s")}}
    report = json.dumps({**summary, "samples": samples}, indent=2) + "\n"
    (root / "results.json").write_text(report)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
