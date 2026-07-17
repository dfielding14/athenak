#!/usr/bin/env python3
"""Read a complete AthenaK MHD snapshot and merged tracks with MPI."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


PYTHON_VIS = Path(__file__).resolve().parents[1]
if str(PYTHON_VIS) not in sys.path:
    sys.path.insert(0, str(PYTHON_VIS))

from cr_visualization.cr_data import (  # noqa: E402
    discover_rank_files,
    read_merged_track_partition,
    read_rank_meshblocks,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhd-rank0", required=True, type=Path)
    parser.add_argument("--merged-tracks", required=True, type=Path)
    parser.add_argument("--quantities", nargs="+")
    parser.add_argument("--track-fields", nargs="+")
    parser.add_argument("--time-min", type=float)
    parser.add_argument("--time-max", type=float)
    parser.add_argument("--time-stride", type=int, default=1)
    parser.add_argument("--particle-batch", type=int, default=4)
    return parser.parse_args()


def read_local_data(args: argparse.Namespace, comm) -> tuple[list[dict], dict]:
    """Read this MPI process's MHD files and contiguous particle-row range."""

    rank = comm.Get_rank()
    size = comm.Get_size()
    rank_files = discover_rank_files(args.mhd_rank0) if rank == 0 else None
    rank_files = [Path(path) for path in comm.bcast(rank_files, root=0)]

    local_meshblocks = []
    for filename in rank_files[rank::size]:
        local_meshblocks.extend(
            read_rank_meshblocks(filename, quantities=args.quantities)
        )
    local_tracks = read_merged_track_partition(
        args.merged_tracks,
        partition=rank,
        num_partitions=size,
        time_min=args.time_min,
        time_max=args.time_max,
        time_stride=args.time_stride,
        fields=args.track_fields,
        particle_batch=args.particle_batch,
    )
    return local_meshblocks, local_tracks


def array_bytes(meshblocks: list[dict]) -> int:
    total = 0
    for block in meshblocks:
        total += sum(block[name].nbytes for name in block["VariableNames"])
    return total


def main() -> None:
    args = parse_args()
    try:
        from mpi4py import MPI
    except ImportError as error:
        raise SystemExit("read_full_dataset_mpi.py requires mpi4py") from error

    comm = MPI.COMM_WORLD
    meshblocks, tracks = read_local_data(args, comm)
    local_counts = np.array(
        [len(meshblocks), tracks["values"].shape[0], array_bytes(meshblocks),
         tracks["values"].nbytes],
        dtype=np.int64,
    )
    total_counts = np.zeros_like(local_counts)
    comm.Reduce(local_counts, total_counts, op=MPI.SUM, root=0)
    if comm.Get_rank() == 0:
        print(
            "distributed read complete: "
            f"MeshBlocks={total_counts[0]} particles={total_counts[1]} "
            f"MHD_GiB={total_counts[2] / 2**30:.3f} "
            f"tracks_GiB={total_counts[3] / 2**30:.3f}"
        )

    # Pass meshblocks and tracks to the distributed visualization pipeline here.


if __name__ == "__main__":
    main()
