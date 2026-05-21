#!/usr/bin/env python3
"""Time short drag-particle runs for deposition/interpolation profiling."""

import argparse
import csv
import subprocess
import sys
import tempfile
import time


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", default="build/src/athena")
    parser.add_argument(
        "--input", default="inputs/particles/hello_drag_particles.athinput"
    )
    parser.add_argument("--ppc", type=float, nargs="+", default=[1.0, 4.0, 16.0])
    parser.add_argument(
        "--coupling",
        nargs="+",
        default=["host_cell", "cloud_in_cell"],
        choices=["host_cell", "cloud_in_cell"],
    )
    parser.add_argument("--tlim", type=float, default=0.05)
    parser.add_argument("--output", default="-")
    return parser.parse_args()


def run_case(args, ppc, coupling):
    with tempfile.TemporaryDirectory(prefix="drag-profile-") as run_dir:
        command = [
            args.athena,
            "-i",
            args.input,
            "-d",
            run_dir,
            f"particles/ppc={ppc}",
            f"drag_particles/interpolation={coupling}",
            f"drag_particles/deposition={coupling}",
            f"time/tlim={args.tlim}",
        ]
        start = time.perf_counter()
        process = subprocess.run(command, check=False, capture_output=True, text=True)
        elapsed = time.perf_counter() - start
    if process.returncode != 0:
        raise RuntimeError(
            f"case failed for ppc={ppc}, coupling={coupling}\n"
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}"
        )
    return elapsed


def main():
    args = parse_args()
    rows = []
    for ppc in args.ppc:
        for coupling in args.coupling:
            elapsed = run_case(args, ppc, coupling)
            rows.append({"ppc": ppc, "coupling": coupling, "seconds": elapsed})

    stream = None if args.output == "-" else open(args.output, "w", newline="")
    try:
        writer = csv.DictWriter(stream or sys.stdout, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if stream is not None:
            stream.close()


if __name__ == "__main__":
    main()
