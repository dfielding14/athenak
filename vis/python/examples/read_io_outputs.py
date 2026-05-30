#!/usr/bin/env python3
"""Print a compact summary of AthenaK binary, PDF, or spherical-slice output."""

from pathlib import Path
import argparse
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import bin_convert  # noqa: E402
from io_reader_common import (  # noqa: E402
    add_reader_limit_arguments,
    reader_limits_from_args,
)
from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("format", choices=("bin", "cbin", "pdf", "sphslice"))
    parser.add_argument("path", type=Path)
    parser.add_argument(
        "--assemble-shards",
        action="store_true",
        help="discover sibling rank/node shards for binary outputs",
    )
    add_reader_limit_arguments(parser)
    args = parser.parse_args()
    limits = reader_limits_from_args(args)

    if args.format == "bin":
        data = bin_convert.read_binary(
            str(args.path), args.assemble_shards, limits=limits
        )
        print(f"binary meshblocks={data['n_mbs']} variables={data['var_names']}")
    elif args.format == "cbin":
        data = bin_convert.read_coarsened_binary(
            str(args.path), args.assemble_shards, limits=limits
        )
        print(f"coarsened meshblocks={data['n_mbs']} variables={data['var_names']}")
    elif args.format == "pdf":
        data = read_pdf(str(args.path), limits=limits)
        dims = [entry["variable"] for entry in data["header"]["dimensions"]]
        print(f"pdf shape={data['pdf'].shape} axes={dims} time={data['time']}")
    else:
        data = read_sphslice(str(args.path), limits=limits)
        print(
            f"sphslice shape={data['data'].shape} variables={data['variables']} "
            f"radius={data['radius']} time={data['time']}"
        )


if __name__ == "__main__":
    main()
