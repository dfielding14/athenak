#!/usr/bin/env python3
"""Read AthenaK lagrangian_mc particle thermodynamic history files."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np


FILE_HEADER = struct.Struct("32siiii")
BLOCK_PREFIX = struct.Struct("16siiiiii")


def read_history(path: str | Path) -> dict[str, np.ndarray]:
    path = Path(path)
    rows_i: list[np.ndarray] = []
    rows_r: list[np.ndarray] = []
    nscalars = None
    real_dtype = None

    with path.open("rb") as handle:
        header = handle.read(FILE_HEADER.size)
        if len(header) != FILE_HEADER.size:
            raise ValueError(f"{path} is too small to be a thermo-history file")
        magic, version, real_size, file_nscalars, _ = FILE_HEADER.unpack(header)
        if not magic.rstrip(b"\0").startswith(b"ATHK_PRTCL_THERMO_HISTORY"):
            raise ValueError(f"{path} has an unrecognized magic string")
        if version != 1:
            raise ValueError(f"unsupported thermo-history version {version}")
        real_dtype = np.float64 if real_size == 8 else np.float32
        nscalars = file_nscalars

        while True:
            prefix = handle.read(BLOCK_PREFIX.size)
            if not prefix:
                break
            if len(prefix) != BLOCK_PREFIX.size:
                raise ValueError("truncated block header")
            magic, block_version, nrecords, block_nscalars, int_per, real_per, cycle = (
                BLOCK_PREFIX.unpack(prefix)
            )
            time = struct.unpack("d" if real_size == 8 else "f", handle.read(real_size))[0]
            if not magic.rstrip(b"\0").startswith(b"ATHKTHPBLK"):
                raise ValueError("unrecognized block magic")
            if block_version != 1 or block_nscalars != nscalars:
                raise ValueError("incompatible block metadata")
            ints = np.fromfile(handle, dtype=np.int32, count=nrecords * int_per)
            reals = np.fromfile(handle, dtype=real_dtype, count=nrecords * real_per)
            if ints.size != nrecords * int_per or reals.size != nrecords * real_per:
                raise ValueError("truncated block payload")
            rows_i.append(ints.reshape(nrecords, int_per))
            rows_r.append(reals.reshape(nrecords, real_per))

    if rows_i:
        ints = np.vstack(rows_i)
        reals = np.vstack(rows_r)
    else:
        ints = np.empty((0, 4), dtype=np.int32)
        reals = np.empty((0, 12 + (nscalars or 0)), dtype=real_dtype or np.float64)

    columns = {
        "time": reals[:, 0],
        "cycle": ints[:, 0],
        "tag": ints[:, 1],
        "seed_id": ints[:, 2],
        "x1": reals[:, 1],
        "x2": reals[:, 2],
        "x3": reals[:, 3],
        "gid": ints[:, 3],
        "rho": reals[:, 4],
        "pressure": reals[:, 5],
        "temperature": reals[:, 6],
        "specific_entropy": reals[:, 7],
        "internal_energy": reals[:, 8],
        "v1": reals[:, 9],
        "v2": reals[:, 10],
        "v3": reals[:, 11],
    }
    for n in range(nscalars or 0):
        columns[f"scalar{n}"] = reals[:, 12 + n]
    return columns


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="Path to *.thp thermo-history file")
    parser.add_argument("--npz", help="Optional output .npz file")
    args = parser.parse_args()

    data = read_history(args.path)
    print(f"records: {len(data['time'])}")
    print("columns:", " ".join(data))
    if args.npz:
        np.savez(args.npz, **data)


if __name__ == "__main__":
    main()
