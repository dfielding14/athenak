#!/usr/bin/env python3
"""Regression test for the lagrangian_mc thermo-history reader."""

from __future__ import annotations

import struct
import tempfile
from pathlib import Path

import numpy as np

from read_prtcl_thermo_history import read_history


FILE_HEADER = struct.Struct("32siiii")
BLOCK_PREFIX_V2 = struct.Struct("16siiiii")


def _padded(value: bytes, size: int) -> bytes:
    return value + b"\0" * (size - len(value))


def _write_block(
    handle, time: float, cycle: int,
    records: list[tuple[int, int, int, int, list[float]]],
) -> None:
    int_per_record = 4
    real_per_record = 8
    handle.write(BLOCK_PREFIX_V2.pack(_padded(b"ATHKTHPBLK", 16), 2, len(records),
                                      int_per_record, real_per_record, cycle))
    handle.write(struct.pack("d", time))
    ints = []
    reals = []
    for cycle_i, tag, seed_id, gid, real_values in records:
        ints.extend([cycle_i, tag, seed_id, gid])
        reals.extend(real_values)
    np.asarray(ints, dtype=np.int32).tofile(handle)
    np.asarray(reals, dtype=np.float64).tofile(handle)


def test_v2_named_columns() -> None:
    fields = ["density", "temperature", "mach", "scalar0"]
    names_blob = "\n".join(fields).encode()

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "synthetic.thp"
        with path.open("wb") as handle:
            handle.write(FILE_HEADER.pack(_padded(b"ATHK_PRTCL_THERMO_HISTORY", 32),
                                          2, 8, len(fields), len(names_blob)))
            handle.write(names_blob)
            _write_block(handle, 0.0, 0, [
                (0, 10, 1, 4, [0.0, 0.1, 0.2, 0.3, 1.0, 0.6, 0.01, 0.25]),
                (0, 11, 1, 4, [0.0, 0.4, 0.5, 0.6, 2.0, 0.7, 0.02, 0.25]),
            ])
            _write_block(handle, 0.1, 5, [
                (5, 10, 1, 4, [0.1, 0.2, 0.3, 0.4, 1.5, 0.65, 0.03, 0.25]),
            ])

        data = read_history(path)
        assert list(data) == [
            "time", "cycle", "tag", "seed_id", "x1", "x2", "x3", "gid", *fields,
        ]
        assert len(data["time"]) == 3
        np.testing.assert_array_equal(data["tag"], np.array([10, 11, 10]))
        np.testing.assert_allclose(data["time"], np.array([0.0, 0.0, 0.1]))
        np.testing.assert_allclose(data["density"], np.array([1.0, 2.0, 1.5]))
        np.testing.assert_allclose(data["mach"], np.array([0.01, 0.02, 0.03]))
        np.testing.assert_allclose(data["scalar0"], np.array([0.25, 0.25, 0.25]))


def main() -> None:
    test_v2_named_columns()
    print("prtcl_thermo_history reader test passed")


if __name__ == "__main__":
    main()
