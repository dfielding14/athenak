#!/usr/bin/env python3
"""Inspect and validate AthenaK cosmic-ray tracer particle outputs."""

from __future__ import annotations

import argparse
import math
import re
import struct
from collections import Counter
from pathlib import Path
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import Sequence


RESTART_FIELDS = 17
PPD_FIELDS = 4


def _rank_from_path(path: Path) -> int:
    match = re.search(r"rank_(\d+)", str(path))
    return int(match.group(1)) if match else 0


def _file_number(path: Path) -> int:
    match = re.search(r"\.(\d+)\.[^.]+$", path.name)
    return int(match.group(1)) if match else -1


def _restart_files(path: Path, latest: bool = True) -> List[Path]:
    if path.is_file():
        return [path]

    files = sorted(path.rglob("*.prst"))
    if not latest or not files:
        return files

    max_number = max(_file_number(file_path) for file_path in files)
    return [file_path for file_path in files if _file_number(file_path) == max_number]


def _ppd_file(path: Path, latest: bool = True) -> Optional[Path]:
    if path.is_file() and path.suffix == ".ppd":
        return path

    ppd_root = path / "ppd" if (path / "ppd").is_dir() else path
    files = sorted(ppd_root.glob("*.ppd"))
    if not files:
        return None
    if latest:
        return max(files, key=_file_number)
    return files[-1]


def _try_restart_decode(blob: bytes, real_size: int) -> Optional[Dict]:
    fmt = "d" if real_size == 8 else "f"
    if len(blob) < 3 * real_size:
        return None

    header = struct.unpack("<" + 3 * fmt, blob[: 3 * real_size])
    count_value = header[2]
    count = int(round(count_value))
    if count < 0 or abs(count_value - count) > 1.0e-6:
        return None

    expected_size = real_size * (3 + RESTART_FIELDS * count)
    if len(blob) != expected_size:
        return None

    data_fmt = "<" + str(RESTART_FIELDS * count) + fmt
    data = struct.unpack(data_fmt, blob[3 * real_size:]) if count else ()
    return {
        "time": header[0],
        "dt": header[1],
        "count": count,
        "data": data,
        "real_size": real_size,
    }


def read_prst_file(path: Path) -> Dict:
    """Read one per-rank particle restart file and validate its length."""
    blob = path.read_bytes()
    decoded = _try_restart_decode(blob, 8)
    if decoded is None:
        decoded = _try_restart_decode(blob, 4)
    if decoded is None:
        raise ValueError(f"{path} is not a well-formed particle restart file")

    particles = []
    data = decoded["data"]
    for offset in range(0, len(data), RESTART_FIELDS):
        record = data[offset: offset + RESTART_FIELDS]
        pgid = int(round(record[0]))
        tag = int(round(record[1]))
        species = int(round(record[2]))
        reals = tuple(float(value) for value in record[3:])
        if not all(math.isfinite(value) for value in reals):
            raise ValueError(f"{path} contains non-finite particle data")
        particles.append({
            "pgid": pgid,
            "tag": tag,
            "species": species,
            "reals": reals,
        })

    return {
        "path": path,
        "rank": _rank_from_path(path),
        "time": decoded["time"],
        "dt": decoded["dt"],
        "count": decoded["count"],
        "real_size": decoded["real_size"],
        "particles": particles,
    }


def read_ppd_file(path: Path) -> Dict:
    """Read one shared compact particle position dump."""
    blob = path.read_bytes()
    header_end = blob.find(b"\n")
    if header_end < 0:
        raise ValueError(f"{path} has no text header")

    header = blob[:header_end].decode("ascii", errors="replace")
    payload = blob[header_end + 1:]
    record_size = PPD_FIELDS * struct.calcsize("<f")
    if len(payload) % record_size != 0:
        raise ValueError(f"{path} has malformed ppd payload length")

    particles = []
    for record in struct.iter_unpack("<ffff", payload):
        species = int(round(record[0]))
        xyz = tuple(float(value) for value in record[1:])
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"{path} contains non-finite positions")
        particles.append({"species": species, "xyz": xyz})

    return {
        "path": path,
        "header": header,
        "count": len(particles),
        "species_counts": Counter(p["species"] for p in particles),
        "particles": particles,
    }


def _histogram_record(path: Path, marker: bytes, count: int) -> List[int]:
    blob = path.read_bytes()
    index = blob.rfind(marker)
    if index < 0:
        raise ValueError(f"{path} does not contain a recognizable histogram header")

    line_end = blob.find(b"\n", index)
    if line_end < 0:
        raise ValueError(f"{path} has a malformed histogram header")

    start = line_end + 1
    while start < len(blob) and blob[start:start + 1] in (b"\n", b"\r"):
        start += 1

    byte_count = count * struct.calcsize("<i")
    end = start + byte_count
    if end > len(blob):
        raise ValueError(f"{path} has an incomplete histogram record")
    return list(struct.unpack("<" + str(count) + "i", blob[start:end]))


def read_df_file(path: Path, nspecies: int, nbin: int) -> List[List[int]]:
    values = _histogram_record(
        path, b"# AthenaK particle distribution function", nspecies * nbin)
    return [values[sp * nbin: (sp + 1) * nbin] for sp in range(nspecies)]


def read_dxh_file(path: Path, nspecies: int, nbin: int) -> List[List[List[int]]]:
    values = _histogram_record(
        path, b"# AthenaK particle displacement histogram", 3 * nspecies * nbin)
    result = []
    for sp in range(nspecies):
        species_values = []
        base = sp * 3 * nbin
        for component in range(3):
            start = base + component * nbin
            species_values.append(values[start: start + nbin])
        result.append(species_values)
    return result


def summarize_restart(path: Path, latest: bool = True) -> Dict:
    files = _restart_files(path, latest=latest)
    if not files:
        raise FileNotFoundError(f"No particle restart files found below {path}")

    restarts = [read_prst_file(file_path) for file_path in files]
    total = sum(restart["count"] for restart in restarts)
    rank_counts = Counter()
    species_counts = Counter()
    tags = set()
    duplicates = []

    for restart in restarts:
        rank_counts[restart["rank"]] += restart["count"]
        for particle in restart["particles"]:
            key = (particle["species"], particle["tag"])
            if key in tags:
                duplicates.append(key)
            tags.add(key)
            species_counts[particle["species"]] += 1

    if duplicates:
        raise ValueError(f"Duplicate particle tags within species: {duplicates[:5]}")

    return {
        "files": files,
        "total": total,
        "rank_counts": dict(sorted(rank_counts.items())),
        "species_counts": dict(sorted(species_counts.items())),
        "restart": restarts,
    }


def summarize_ppd(path: Path, latest: bool = True) -> Optional[Dict]:
    ppd_path = _ppd_file(path, latest=latest)
    if ppd_path is None:
        return None
    return read_ppd_file(ppd_path)


def validate_expected_counts(summary: Dict, expected_total: Optional[int],
                             expected_species: Optional[Sequence[int]]) -> None:
    if expected_total is not None and summary["total"] != expected_total:
        raise ValueError(
            f"Expected {expected_total} particles, found {summary['total']}")

    if expected_species is None:
        return

    for species, expected in enumerate(expected_species):
        actual = summary["species_counts"].get(species, 0)
        if actual != expected:
            raise ValueError(
                f"Expected species {species} count {expected}, found {actual}")


def validate_histograms(path: Path, expected_species: Sequence[int],
                        nbin: int) -> Dict:
    nspecies = len(expected_species)
    df_path = path / "df"
    dxh_path = path / "dxh"
    if df_path.is_dir():
        df_file = max(df_path.glob("*.df"), key=_file_number)
    else:
        df_file = df_path
    if dxh_path.is_dir():
        dxh_file = max(dxh_path.glob("*.dxh"), key=_file_number)
    else:
        dxh_file = dxh_path

    df = read_df_file(df_file, nspecies, nbin)
    dxh = read_dxh_file(dxh_file, nspecies, nbin)
    for species, expected in enumerate(expected_species):
        df_total = sum(df[species])
        if df_total != expected:
            raise ValueError(
                f"DF species {species} sums to {df_total}, expected {expected}")
        for component in range(3):
            dxh_total = sum(dxh[species][component])
            if dxh_total != expected:
                raise ValueError(
                    f"DXH species {species} component {component} sums to "
                    f"{dxh_total}, expected {expected}")
    return {"df": df, "dxh": dxh, "df_file": df_file, "dxh_file": dxh_file}


def _format_counts(counts: Dict[int, int]) -> str:
    return " ".join(f"{key}={value}" for key, value in sorted(counts.items()))


def _parse_counts(text: Optional[str]) -> Optional[List[int]]:
    if text is None:
        return None
    return [int(value) for value in text.split(",") if value]


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Inspect AthenaK cosmic-ray tracer particle outputs.")
    parser.add_argument("path", type=Path, help="Run directory, prst dir, or file")
    parser.add_argument("--all", action="store_true",
                        help="Inspect all restart files instead of latest output")
    parser.add_argument("--expected-total", type=int)
    parser.add_argument("--expected-species-counts",
                        help="Comma-separated expected counts, e.g. 512,512")
    parser.add_argument("--histogram-nbin", type=int,
                        help="Also validate df/dxh histogram sums")
    args = parser.parse_args(list(argv) if argv is not None else None)

    expected_species = _parse_counts(args.expected_species_counts)
    summary = summarize_restart(args.path, latest=not args.all)
    validate_expected_counts(summary, args.expected_total, expected_species)

    print("restart:")
    print(f"  files: {len(summary['files'])}")
    print(f"  total: {summary['total']}")
    print(f"  ranks: {_format_counts(summary['rank_counts'])}")
    print(f"  species: {_format_counts(summary['species_counts'])}")

    ppd_summary = summarize_ppd(args.path, latest=not args.all)
    if ppd_summary is not None:
        print("ppd:")
        print(f"  file: {ppd_summary['path']}")
        print(f"  total: {ppd_summary['count']}")
        print(f"  species: {_format_counts(ppd_summary['species_counts'])}")

    if args.histogram_nbin is not None:
        if expected_species is None:
            raise ValueError("--histogram-nbin requires --expected-species-counts")
        hist = validate_histograms(args.path, expected_species,
                                   args.histogram_nbin)
        print("histograms:")
        print(f"  df: {hist['df_file']}")
        print(f"  dxh: {hist['dxh_file']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
