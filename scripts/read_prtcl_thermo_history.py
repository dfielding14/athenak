#!/usr/bin/env python3
"""Read AthenaK lagrangian_mc particle thermodynamic history files."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np


FILE_HEADER = struct.Struct("32siiii")
BLOCK_PREFIX_V1 = struct.Struct("16siiiiii")
BLOCK_PREFIX_V2 = struct.Struct("16siiiii")


def _read_exact(handle, size: int, label: str) -> bytes:
    data = handle.read(size)
    if len(data) != size:
        raise ValueError(f"truncated {label}")
    return data


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
        magic, version, real_size, count_or_nfields, reserved_or_names_size = (
            FILE_HEADER.unpack(header)
        )
        if not magic.rstrip(b"\0").startswith(b"ATHK_PRTCL_THERMO_HISTORY"):
            raise ValueError(f"{path} has an unrecognized magic string")
        if version not in (1, 2):
            raise ValueError(f"unsupported thermo-history version {version}")
        if real_size not in (4, 8):
            raise ValueError(f"unsupported real size {real_size}")
        real_dtype = np.float64 if real_size == 8 else np.float32
        if version == 1:
            nscalars = count_or_nfields
            field_names = [
                "rho",
                "pressure",
                "temperature",
                "specific_entropy",
                "internal_energy",
                "v1",
                "v2",
                "v3",
            ] + [f"scalar{n}" for n in range(nscalars)]
        else:
            nfields = count_or_nfields
            names_size = reserved_or_names_size
            names_blob = _read_exact(handle, names_size, "field schema").decode()
            field_names = names_blob.split("\n") if names_blob else []
            if len(field_names) != nfields:
                raise ValueError("field schema length does not match header")
            nscalars = 0

        while True:
            prefix_size = BLOCK_PREFIX_V1.size if version == 1 else BLOCK_PREFIX_V2.size
            prefix = handle.read(prefix_size)
            if not prefix:
                break
            if len(prefix) != prefix_size:
                raise ValueError("truncated block header")
            if version == 1:
                magic, block_version, nrecords, block_nscalars, int_per, real_per, cycle = (
                    BLOCK_PREFIX_V1.unpack(prefix)
                )
                if block_nscalars != nscalars:
                    raise ValueError("incompatible block scalar metadata")
            else:
                magic, block_version, nrecords, int_per, real_per, cycle = (
                    BLOCK_PREFIX_V2.unpack(prefix)
                )
            time = struct.unpack("d" if real_size == 8 else "f",
                                 _read_exact(handle, real_size, "block time"))[0]
            if not magic.rstrip(b"\0").startswith(b"ATHKTHPBLK"):
                raise ValueError("unrecognized block magic")
            if block_version != version:
                raise ValueError("incompatible block metadata")
            expected_real_per = (
                12 + (nscalars or 0) if version == 1 else 4 + len(field_names)
            )
            if int_per != 4 or real_per != expected_real_per:
                raise ValueError("incompatible block record size")
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
        if version == 1:
            real_per = 12 + (nscalars or 0)
        else:
            real_per = 4 + len(field_names)
        reals = np.empty((0, real_per), dtype=real_dtype or np.float64)

    columns = {
        "time": reals[:, 0],
        "cycle": ints[:, 0],
        "tag": ints[:, 1],
        "seed_id": ints[:, 2],
        "x1": reals[:, 1],
        "x2": reals[:, 2],
        "x3": reals[:, 3],
        "gid": ints[:, 3],
    }
    for n, name in enumerate(field_names):
        columns[name] = reals[:, 4 + n]
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
