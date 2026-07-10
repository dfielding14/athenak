#!/usr/bin/env python3
"""Stream-compare two reducer output trees without loading PDF payloads into RAM."""

import argparse
import array
import math
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-9)
    parser.add_argument("--global-rtol", type=float, default=1.0e-12)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    parser.add_argument("--chunk-doubles", type=int, default=1_048_576)
    return parser.parse_args()


def relative_files(root, pattern):
    return {path.relative_to(root) for path in root.glob(pattern) if path.is_file()}


def compare_payload(
    reference,
    candidate,
    rtol,
    global_rtol,
    atol,
    chunk_doubles,
):
    if reference.stat().st_size != candidate.stat().st_size:
        raise RuntimeError(
            f"size mismatch for {reference.name}: "
            f"{reference.stat().st_size} != {candidate.stat().st_size}"
        )
    if reference.stat().st_size % 8:
        raise RuntimeError(f"payload size is not a multiple of double precision: {reference}")

    max_absolute = 0.0
    max_relative = 0.0
    file_scale = 0.0
    values = 0
    chunk_bytes = chunk_doubles * 8

    with reference.open("rb") as left, candidate.open("rb") as right:
        while True:
            left_bytes = left.read(chunk_bytes)
            right_bytes = right.read(chunk_bytes)
            if not left_bytes and not right_bytes:
                break
            if len(left_bytes) != len(right_bytes):
                raise RuntimeError(f"short read mismatch for {reference.name}")
            left_values = array.array("d")
            right_values = array.array("d")
            left_values.frombytes(left_bytes)
            right_values.frombytes(right_bytes)
            for index, (a, b) in enumerate(zip(left_values, right_values)):
                if not math.isfinite(a) or not math.isfinite(b):
                    raise RuntimeError(
                        f"non-finite value in {reference.name} at double {values + index}: "
                        f"{a!r}, {b!r}"
                    )
                file_scale = max(file_scale, abs(a), abs(b))
            values += len(left_values)

    compared_values = 0
    with reference.open("rb") as left, candidate.open("rb") as right:
        while True:
            left_bytes = left.read(chunk_bytes)
            right_bytes = right.read(chunk_bytes)
            if not left_bytes and not right_bytes:
                break
            if len(left_bytes) != len(right_bytes):
                raise RuntimeError(f"short read mismatch for {reference.name}")
            left_values = array.array("d")
            right_values = array.array("d")
            left_values.frombytes(left_bytes)
            right_values.frombytes(right_bytes)
            for index, (a, b) in enumerate(zip(left_values, right_values)):
                absolute = abs(a - b)
                scale = max(abs(a), abs(b))
                relative = absolute / scale if scale else 0.0
                max_absolute = max(max_absolute, absolute)
                max_relative = max(max_relative, relative)
                tolerance = atol + rtol * scale + global_rtol * file_scale
                if absolute > tolerance:
                    raise RuntimeError(
                        f"value mismatch in {reference.name} at double {compared_values + index}: "
                        f"{a:.17g} != {b:.17g}; abs={absolute:.6g}, rel={relative:.6g}, "
                        f"tolerance={tolerance:.6g}"
                    )
            compared_values += len(left_values)
    if compared_values != values:
        raise RuntimeError(f"payload changed while comparing: {reference}")
    return max_absolute, max_relative, values, file_scale


def main():
    args = parse_args()
    if not args.reference.is_dir() or not args.candidate.is_dir():
        raise SystemExit("reference and candidate must both be output directories")
    if args.chunk_doubles <= 0:
        raise SystemExit("--chunk-doubles must be positive")
    if args.rtol < 0.0 or args.global_rtol < 0.0 or args.atol < 0.0:
        raise SystemExit("comparison tolerances must be non-negative")

    headers = relative_files(args.reference, "*/gotham.header.pdf")
    candidate_headers = relative_files(args.candidate, "*/gotham.header.pdf")
    if headers != candidate_headers:
        raise SystemExit("header file sets differ")
    for relative in sorted(headers):
        if (args.reference / relative).read_bytes() != (args.candidate / relative).read_bytes():
            raise SystemExit(f"header mismatch: {relative}")

    payloads = relative_files(args.reference, "*/gotham.*.pdf") - headers
    candidate_payloads = relative_files(args.candidate, "*/gotham.*.pdf") - candidate_headers
    if payloads != candidate_payloads:
        raise SystemExit("payload file sets differ")
    if not payloads:
        raise SystemExit("no payloads found")

    worst_absolute = 0.0
    worst_relative = 0.0
    total_values = 0
    for relative in sorted(payloads):
        absolute, relative_error, values, file_scale = compare_payload(
            args.reference / relative,
            args.candidate / relative,
            rtol=args.rtol,
            global_rtol=args.global_rtol,
            atol=args.atol,
            chunk_doubles=args.chunk_doubles,
        )
        worst_absolute = max(worst_absolute, absolute)
        worst_relative = max(worst_relative, relative_error)
        total_values += values
        print(
            f"matched {relative}: values={values} "
            f"scale={file_scale:.6g} max_abs={absolute:.6g} max_rel={relative_error:.6g}"
        )

    print(
        f"matched {len(payloads)} payloads and {total_values} doubles; "
        f"global_max_abs={worst_absolute:.6g} global_max_rel={worst_relative:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
