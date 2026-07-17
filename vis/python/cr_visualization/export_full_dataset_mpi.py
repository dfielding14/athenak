#!/usr/bin/env python3
"""Export a distributed AthenaK MHD snapshot and tracks as XDMF/HDF5."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


PYTHON_VIS = Path(__file__).resolve().parents[1]
if str(PYTHON_VIS) not in sys.path:
    sys.path.insert(0, str(PYTHON_VIS))

from cr_visualization.read_full_dataset_mpi import read_local_data  # noqa: E402
from cr_visualization.xdmf_export import (  # noqa: E402
    piece_path,
    visualization_track_fields,
    write_visualization_piece,
    write_xdmf_collection,
    xmf_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhd-rank0", required=True, type=Path)
    parser.add_argument("--merged-tracks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--quantities", nargs="+")
    parser.add_argument("--track-fields", nargs="+")
    parser.add_argument("--time-min", type=float)
    parser.add_argument("--time-max", type=float)
    parser.add_argument("--time-stride", type=int, default=1)
    parser.add_argument("--particle-batch", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        from mpi4py import MPI
    except ImportError as error:
        raise SystemExit("export_full_dataset_mpi.py requires mpi4py") from error

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    output = xmf_path(args.output).resolve()
    payload = piece_path(output, rank)
    local_exists = payload.exists() or (rank == 0 and output.exists())
    if comm.allreduce(local_exists, op=MPI.LOR):
        raise SystemExit(f"refusing to replace an existing export for {output}")

    args.track_fields = visualization_track_fields(args.track_fields)
    meshblocks, tracks = read_local_data(args, comm)
    geometry = None
    if rank == 0 and meshblocks:
        geometry = (
            meshblocks[0].get("DomainBounds"),
            meshblocks[0].get("PeriodicAxes", ()),
        )
    domain_bounds, periodic_axes = comm.bcast(geometry, root=0)
    metadata = write_visualization_piece(
        payload,
        meshblocks,
        tracks,
        track_geometry="complete",
        particle_batch=args.particle_batch,
        domain_bounds=domain_bounds,
        periodic_axes=periodic_axes,
    )
    pieces = comm.gather(metadata, root=0)
    if rank == 0:
        write_xdmf_collection(output, pieces)

    local_counts = np.array(
        [
            len(meshblocks),
            tracks["values"].shape[0],
            tracks["values"].shape[0] * tracks["values"].shape[1],
        ],
        dtype=np.int64,
    )
    total_counts = np.zeros_like(local_counts)
    comm.Reduce(local_counts, total_counts, op=MPI.SUM, root=0)
    if rank == 0:
        print(
            f"{output}: MeshBlocks={total_counts[0]} "
            f"particles={total_counts[1]} track_points={total_counts[2]} "
            f"pieces={comm.Get_size()}"
        )


if __name__ == "__main__":
    main()
