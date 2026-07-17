#!/usr/bin/env python3
"""Read one AthenaK MeshBlock and the trajectories that pass through it."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PYTHON_VIS = Path(__file__).resolve().parents[1]
if str(PYTHON_VIS) not in sys.path:
    sys.path.insert(0, str(PYTHON_VIS))

from cr_visualization.cr_data import (  # noqa: E402
    find_track_visits,
    read_meshblock,
    read_merged_track_subset,
    write_meshblock_bundle,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhd-bin", required=True, type=Path)
    parser.add_argument("--raw-trk", required=True, nargs="+", type=Path)
    parser.add_argument("--merged-tracks", required=True, type=Path)
    parser.add_argument("--meshblock", type=int, default=0)
    parser.add_argument("--quantities", nargs="+")
    parser.add_argument("--track-fields", nargs="+")
    parser.add_argument("--species", nargs="+", type=int)
    parser.add_argument("--max-particles", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--time-min", type=float)
    parser.add_argument("--time-max", type=float)
    parser.add_argument("--time-stride", type=int, default=1)
    parser.add_argument(
        "--inside-only",
        action="store_true",
        help="replace samples outside the MeshBlock with NaN",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional combined MeshBlock and trajectory HDF5 file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    meshblock = read_meshblock(
        args.mhd_bin,
        meshblock_index=args.meshblock,
        quantities=args.quantities,
    )
    visits = find_track_visits(
        args.raw_trk,
        meshblock["Bounds"],
        time_min=args.time_min,
        time_max=args.time_max,
    )
    tracks = read_merged_track_subset(
        args.merged_tracks,
        visits.output_tags,
        bounds=meshblock["Bounds"],
        species=args.species,
        max_particles=args.max_particles,
        seed=args.seed,
        time_min=args.time_min,
        time_max=args.time_max,
        time_stride=args.time_stride,
        fields=args.track_fields,
        inside_only=args.inside_only,
    )

    shape = meshblock[meshblock["VariableNames"][0]].shape
    print(
        f"MHD MeshBlock {args.meshblock}: shape[z,y,x]={shape} "
        f"bounds={meshblock['Bounds'].tolist()}"
    )
    print(
        f"track shard: frames={visits.frames_seen} records={visits.records_seen} "
        f"inside={visits.records_inside} particles={visits.output_tags.size}"
    )
    print(
        f"loaded trajectories: particles={tracks['values'].shape[0]} "
        f"times={tracks['values'].shape[1]} fields={tracks['fields']}"
    )

    # A bundle is convenient for tools that should not reopen the source files.
    if args.output is not None:
        write_meshblock_bundle(args.output, meshblock, tracks)
        print(args.output)


if __name__ == "__main__":
    main()
