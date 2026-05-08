#!/usr/bin/env python3
"""Read AthenaK fixed-column tracked-particle `.trk` files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import struct


HEADER_RE = re.compile(
    r"time=\s*(?P<time>[-+0-9.eE]+)\s+"
    r"nranks=\s*(?P<nranks>\d+)\s+"
    r"cycle=\s*(?P<cycle>\d+)\s+"
    r"ntracked_prtcls=\s*(?P<ntrack>\d+)\s+"
    r"nvalues=\s*(?P<nvalues>\d+)\s+"
    r"variables=\s*(?P<variables>.*)$"
)


@dataclass(frozen=True)
class TrackRecord:
    time: float
    cycle: int
    nranks: int
    ntrack: int
    variables: tuple[str, ...]
    rows: list[dict[str, float]]


def _parse_header(line: bytes) -> tuple[float, int, int, int, tuple[str, ...]]:
    text = line.decode("ascii", errors="strict").strip()
    match = HEADER_RE.search(text)
    if match is None:
        raise ValueError(f"unrecognized tracked-particle header: {text!r}")
    variables = tuple(match.group("variables").split())
    return (
        float(match.group("time")),
        int(match.group("cycle")),
        int(match.group("nranks")),
        int(match.group("ntrack")),
        variables,
    )


def read_trk(path: str | Path) -> list[TrackRecord]:
    path = Path(path)
    data = path.read_bytes()
    records: list[TrackRecord] = []
    pos = 0
    marker = b"# AthenaK tracked particle data"
    while True:
        header_start = data.find(marker, pos)
        if header_start < 0:
            break
        header_end = data.find(b"\n", header_start)
        if header_end < 0:
            raise ValueError(f"unterminated header in {path}")
        time, cycle, nranks, ntrack, variables = _parse_header(
            data[header_start:header_end]
        )
        nvalues = len(variables)
        payload_start = header_end + 1
        payload_nbytes = ntrack * nvalues * 4
        payload_end = payload_start + payload_nbytes
        if payload_end > len(data):
            raise ValueError(
                f"truncated payload in {path}: need {payload_nbytes} bytes after header"
            )
        values = struct.unpack(
            f"<{ntrack * nvalues}f", data[payload_start:payload_end]
        )
        rows: list[dict[str, float]] = []
        for tag in range(ntrack):
            offset = tag * nvalues
            row = {name: values[offset + i] for i, name in enumerate(variables)}
            row.setdefault("tag", float(tag))
            row["track_index"] = float(tag)
            rows.append(row)
        records.append(
            TrackRecord(
                time=time,
                cycle=cycle,
                nranks=nranks,
                ntrack=ntrack,
                variables=variables,
                rows=rows,
            )
        )
        pos = payload_end
    return records


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trk_file", type=Path)
    args = parser.parse_args()
    for record in read_trk(args.trk_file):
        active = sum(1 for row in record.rows if int(round(row["active"])) == 1)
        print(
            f"time={record.time:g} cycle={record.cycle} "
            f"ntrack={record.ntrack} active={active} "
            f"variables={' '.join(record.variables)}"
        )


if __name__ == "__main__":
    main()
