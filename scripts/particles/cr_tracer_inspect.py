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


def _psamp_files(path: Path) -> List[Path]:
    if path.is_file() and path.suffix == ".psamp":
        return [path]

    sample_root = path / "psamp" if (path / "psamp").is_dir() else path
    return sorted(sample_root.rglob("*.psamp"))


def sample_hash(species: int, tag: int) -> int:
    """Return the deterministic hash used by AthenaK psamp selection."""
    value = ((species & 0xFFFFFFFF) << 32) ^ (tag & 0xFFFFFFFF)
    value = (value + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (value ^ (value >> 31)) & 0xFFFFFFFFFFFFFFFF


def sample_selected(species: int, tag: int, sample_species: int,
                    sample_stride: int, sample_offset: int) -> bool:
    if sample_species >= 0 and species != sample_species:
        return False
    return sample_hash(species, tag) % sample_stride == sample_offset


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


def _read_particle_metadata(path: Path) -> Optional[Dict[str, str]]:
    meta_path = path.with_suffix(path.suffix + ".pmeta")
    if not meta_path.exists():
        return None
    metadata: Dict[str, str] = {"path": str(meta_path)}
    for line in meta_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(maxsplit=1)
        if len(parts) == 2:
            metadata[parts[0]] = parts[1]
    return metadata


def _fnv1a64(blob: bytes) -> int:
    value = 1469598103934665603
    for byte in blob:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def _validate_particle_metadata(path: Path, blob: bytes, decoded: Dict,
                                metadata: Optional[Dict[str, str]]) -> str:
    if metadata is None:
        return "legacy_prst_no_metadata"

    if metadata.get("magic") != "AKPRST":
        raise ValueError(f"{path} metadata has invalid magic {metadata.get('magic')}")
    if metadata.get("binary_payload") != "legacy_prst":
        raise ValueError(f"{path} metadata has unsupported binary payload")

    checks = {
        "real_size": decoded["real_size"],
        "particle_count": decoded["count"],
        "byte_count": len(blob),
        "record_real_count": RESTART_FIELDS,
    }
    for key, expected in checks.items():
        actual = int(metadata.get(key, -1))
        if actual != expected:
            raise ValueError(
                f"{path} metadata {key}={actual}, expected {expected}")

    if "checksum_fnv1a64" in metadata:
        expected_checksum = int(metadata["checksum_fnv1a64"])
        actual_checksum = _fnv1a64(blob)
        if actual_checksum != expected_checksum:
            raise ValueError(
                f"{path} checksum mismatch: metadata={expected_checksum} "
                f"actual={actual_checksum}")

    return f"AKPRST-v{metadata.get('version', 'unknown')}"


def read_prst_file(path: Path) -> Dict:
    """Read one per-rank particle restart file and validate its length."""
    blob = path.read_bytes()
    decoded = _try_restart_decode(blob, 8)
    if decoded is None:
        decoded = _try_restart_decode(blob, 4)
    if decoded is None:
        raise ValueError(
            f"{path} is not a well-formed particle restart file; "
            f"byte_count={len(blob)} does not match a float or double legacy "
            "particle restart header/data length")

    metadata = _read_particle_metadata(path)
    format_version = _validate_particle_metadata(path, blob, decoded, metadata)

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
        "format_version": format_version,
        "metadata": metadata,
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


def _histogram_record(path: Path, marker: bytes, count: int,
                      skip_comment_lines: bool = False) -> List[int]:
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
    if skip_comment_lines:
        while start < len(blob) and blob[start:start + 1] == b"#":
            line_end = blob.find(b"\n", start)
            if line_end < 0:
                raise ValueError(f"{path} has a malformed histogram metadata line")
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


def read_scalar_histogram(path: Path, marker: bytes, nspecies: int,
                          nbin: int) -> List[List[int]]:
    values = _histogram_record(path, marker, nspecies * nbin)
    return [values[sp * nbin: (sp + 1) * nbin] for sp in range(nspecies)]


def read_pspec_file(path: Path, nspecies: int, nbin: int) -> List[List[int]]:
    values = _histogram_record(
        path, b"# AthenaK particle spectrum", nspecies * nbin,
        skip_comment_lines=True)
    return [values[sp * nbin: (sp + 1) * nbin] for sp in range(nspecies)]


def read_pspec2_file(path: Path, nspecies: int,
                     nbin1: int, nbin2: int) -> List[List[List[int]]]:
    values = _histogram_record(
        path, b"# AthenaK particle joint spectrum", nspecies * nbin1 * nbin2,
        skip_comment_lines=True)
    result = []
    for species in range(nspecies):
        species_values = []
        base = species * nbin1 * nbin2
        for bin1 in range(nbin1):
            start = base + bin1 * nbin2
            species_values.append(values[start: start + nbin2])
        result.append(species_values)
    return result


def read_pmom_file(path: Path, nspecies: int) -> List[Dict[str, float]]:
    """Read the latest per-species particle moment block from a pmom file."""
    if path.is_dir():
        path = max(path.glob("*.pmom"), key=_file_number)
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values = stripped.split()
        if len(values) not in (14, 26):
            raise ValueError(f"{path} has malformed particle-moment row: {line}")
        rows.append(values)
    if len(rows) < nspecies:
        raise ValueError(f"{path} does not contain {nspecies} particle-moment rows")

    latest = rows[-nspecies:]
    names = [
        "species",
        "count",
        "mean_mu",
        "mean_mu2",
        "anisotropy",
        "mean_dx",
        "mean_dy",
        "mean_dz",
        "mean_dx2",
        "mean_dy2",
        "mean_dz2",
        "mean_dparallel",
        "mean_dparallel2",
        "mean_speed2",
        "mean_vx",
        "mean_vy",
        "mean_vz",
        "mean_vx2",
        "mean_vy2",
        "mean_vz2",
        "mean_vxvy",
        "mean_vxvz",
        "mean_vyvz",
        "mean_dxdy",
        "mean_dxdz",
        "mean_dydz",
    ]
    result = []
    for row in latest:
        parsed = {name: float(value) for name, value in zip(names, row)}
        parsed["species"] = int(round(parsed["species"]))
        parsed["count"] = int(round(parsed["count"]))
        if not all(math.isfinite(value) for value in parsed.values()):
            raise ValueError(f"{path} contains non-finite particle moments")
        result.append(parsed)
    return result


def _parse_metadata_line(line: str) -> Dict[str, str]:
    result = {}
    for token in line.split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        result[key] = value
    return result


def read_psamp_file(path: Path, latest: bool = True) -> Dict:
    """Read one per-rank deterministic particle sample file."""
    blocks = []
    current: Optional[Dict] = None
    columns: Optional[List[str]] = None
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("# AthenaK particle sample"):
            if current is not None:
                blocks.append(current)
            current = {
                "path": path,
                "rank": _rank_from_path(path),
                "metadata": {},
                "columns": [],
                "rows": [],
            }
            columns = None
            continue
        if stripped.startswith("# metadata:"):
            if current is None:
                raise ValueError(f"{path}:{line_number} metadata before header")
            current["metadata"].update(
                _parse_metadata_line(stripped.split(":", 1)[1]))
            continue
        if stripped.startswith("# columns:"):
            if current is None:
                raise ValueError(f"{path}:{line_number} columns before header")
            columns = stripped.split(":", 1)[1].split()
            if len(columns) < 4 or columns[:4] != ["rank", "pgid", "tag", "species"]:
                raise ValueError(f"{path}:{line_number} has malformed columns")
            current["columns"] = columns
            continue
        if stripped.startswith("#"):
            continue
        if current is None or columns is None:
            raise ValueError(f"{path}:{line_number} sample row before columns")

        values = stripped.split()
        if len(values) != len(columns):
            raise ValueError(
                f"{path}:{line_number} has {len(values)} values for "
                f"{len(columns)} columns")
        row = {
            "rank": int(values[0]),
            "pgid": int(values[1]),
            "tag": int(values[2]),
            "species": int(values[3]),
            "fields": {},
        }
        for name, value in zip(columns[4:], values[4:]):
            parsed = float(value)
            if not math.isfinite(parsed):
                raise ValueError(f"{path}:{line_number} has non-finite {name}")
            row["fields"][name] = parsed
        current["rows"].append(row)

    if current is not None:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"{path} does not contain a particle sample block")
    return blocks[-1] if latest else {"path": path, "blocks": blocks}


def summarize_samples(path: Path, latest: bool = True) -> Optional[Dict]:
    files = _psamp_files(path)
    if not files:
        return None

    sample_blocks = [read_psamp_file(file_path, latest=latest) for file_path in files]
    rows = []
    rank_counts = Counter()
    species_counts = Counter()
    tags = set()
    duplicates = []
    field_names = set()
    for block in sample_blocks:
        block_rows = block["rows"] if latest else [
            row for nested in block["blocks"] for row in nested["rows"]]
        for row in block_rows:
            rows.append(row)
            rank_counts[row["rank"]] += 1
            species_counts[row["species"]] += 1
            key = (row["species"], row["tag"])
            if key in tags:
                duplicates.append(key)
            tags.add(key)
            field_names.update(row["fields"].keys())

    if latest and duplicates:
        raise ValueError(f"Duplicate sampled particle tags: {duplicates[:5]}")

    return {
        "files": files,
        "total": len(rows),
        "rank_counts": dict(sorted(rank_counts.items())),
        "species_counts": dict(sorted(species_counts.items())),
        "fields": sorted(field_names),
        "rows": rows,
        "blocks": sample_blocks,
    }


def expected_sample_counts_from_restart(summary: Dict, sample_species: int,
                                        sample_stride: int,
                                        sample_offset: int) -> Dict:
    species_counts = Counter()
    rank_counts = Counter()
    total = 0
    for restart in summary["restart"]:
        for particle in restart["particles"]:
            species = particle["species"]
            tag = particle["tag"]
            if sample_selected(species, tag, sample_species,
                               sample_stride, sample_offset):
                total += 1
                species_counts[species] += 1
                rank_counts[restart["rank"]] += 1
    return {
        "total": total,
        "species_counts": dict(sorted(species_counts.items())),
        "rank_counts": dict(sorted(rank_counts.items())),
    }


def summarize_restart(path: Path, latest: bool = True) -> Dict:
    files = _restart_files(path, latest=latest)
    if not files:
        raise FileNotFoundError(f"No particle restart files found below {path}")

    restarts = [read_prst_file(file_path) for file_path in files]
    total = sum(restart["count"] for restart in restarts)
    rank_counts = Counter()
    species_counts = Counter()
    pgid_counts = Counter()
    real_sizes = Counter()
    format_versions = Counter()
    tags = set()
    duplicates = []

    for restart in restarts:
        rank_counts[restart["rank"]] += restart["count"]
        real_sizes[restart["real_size"]] += 1
        format_versions[restart["format_version"]] += 1
        for particle in restart["particles"]:
            key = (particle["species"], particle["tag"])
            if key in tags:
                duplicates.append(key)
            tags.add(key)
            species_counts[particle["species"]] += 1
            pgid_counts[particle["pgid"]] += 1

    if duplicates:
        raise ValueError(f"Duplicate particle tags within species: {duplicates[:5]}")

    return {
        "files": files,
        "total": total,
        "rank_counts": dict(sorted(rank_counts.items())),
        "species_counts": dict(sorted(species_counts.items())),
        "pgid_counts": dict(sorted(pgid_counts.items())),
        "real_sizes": dict(sorted(real_sizes.items())),
        "format_versions": dict(sorted(format_versions.items())),
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
    optional = {}
    scalar_specs = {
        "drh": b"# AthenaK particle scalar displacement histogram",
        "dparh": b"# AthenaK particle parallel displacement histogram",
    }
    for name, marker in scalar_specs.items():
        hist_path = path / name
        if hist_path.is_dir():
            hist_file = max(hist_path.glob(f"*.{name}"), key=_file_number)
        else:
            hist_file = hist_path
        if hist_file.exists():
            optional[name] = {
                "file": hist_file,
                "histogram": read_scalar_histogram(hist_file, marker, nspecies, nbin),
            }
    pspec_path = path / "pspec"
    if pspec_path.is_dir():
        pspec_file = max(pspec_path.glob("*.pspec"), key=_file_number)
    else:
        pspec_file = pspec_path
    if pspec_file.exists():
        optional["pspec"] = {
            "file": pspec_file,
            "histogram": read_pspec_file(pspec_file, nspecies, nbin),
        }
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
        for name, data in optional.items():
            hist_total = sum(data["histogram"][species])
            if hist_total != expected:
                raise ValueError(
                    f"{name} species {species} sums to {hist_total}, "
                    f"expected {expected}")
    result = {"df": df, "dxh": dxh, "df_file": df_file, "dxh_file": dxh_file}
    result.update(optional)
    return result


def validate_spectra(path: Path, expected_species: Sequence[int],
                     nbin: int) -> Dict:
    nspecies = len(expected_species)
    pspec_path = path / "pspec"
    if pspec_path.is_dir():
        pspec_file = max(pspec_path.glob("*.pspec"), key=_file_number)
    else:
        pspec_file = pspec_path
    spectrum = read_pspec_file(pspec_file, nspecies, nbin)
    for species, expected in enumerate(expected_species):
        total = sum(spectrum[species])
        if total != expected:
            raise ValueError(
                f"PSPEC species {species} sums to {total}, expected {expected}")
    return {"file": pspec_file, "histogram": spectrum}


def validate_joint_spectra(path: Path, expected_species: Sequence[int],
                           nbin1: int, nbin2: int) -> Dict:
    nspecies = len(expected_species)
    pspec2_path = path / "pspec2"
    if pspec2_path.is_dir():
        pspec2_file = max(pspec2_path.glob("*.pspec2"), key=_file_number)
    else:
        pspec2_file = pspec2_path
    spectrum = read_pspec2_file(pspec2_file, nspecies, nbin1, nbin2)
    for species, expected in enumerate(expected_species):
        total = sum(sum(row) for row in spectrum[species])
        if total != expected:
            raise ValueError(
                f"PSPEC2 species {species} sums to {total}, expected {expected}")
    return {"file": pspec2_file, "histogram": spectrum}


def validate_moments(path: Path, expected_species: Sequence[int]) -> List[Dict]:
    nspecies = len(expected_species)
    pmom_path = path / "pmom"
    if pmom_path.is_dir():
        pmom_file = max(pmom_path.glob("*.pmom"), key=_file_number)
    else:
        pmom_file = pmom_path
    moments = read_pmom_file(pmom_file, nspecies)
    for species, expected in enumerate(expected_species):
        if moments[species]["species"] != species:
            raise ValueError(
                f"PMOM row {species} has species {moments[species]['species']}")
        if moments[species]["count"] != expected:
            raise ValueError(
                f"PMOM species {species} count {moments[species]['count']}, "
                f"expected {expected}")
        if abs(moments[species]["mean_mu"]) > 1.0 + 1.0e-12:
            raise ValueError(f"PMOM species {species} has |mean_mu| > 1")
        if moments[species]["mean_mu2"] < -1.0e-12:
            raise ValueError(f"PMOM species {species} has negative mean_mu2")
    return moments


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
    parser.add_argument("--spectrum-nbin", type=int,
                        help="Also validate pspec spectrum sums")
    parser.add_argument("--joint-spectrum-nbin",
                        help="Also validate pspec2 sums as nbin1,nbin2")
    args = parser.parse_args(list(argv) if argv is not None else None)

    expected_species = _parse_counts(args.expected_species_counts)
    summary = summarize_restart(args.path, latest=not args.all)
    validate_expected_counts(summary, args.expected_total, expected_species)

    print("restart:")
    print(f"  files: {len(summary['files'])}")
    print(f"  total: {summary['total']}")
    print(f"  ranks: {_format_counts(summary['rank_counts'])}")
    print(f"  species: {_format_counts(summary['species_counts'])}")
    print(f"  meshblocks: {_format_counts(summary['pgid_counts'])}")
    print(f"  real_sizes: {_format_counts(summary['real_sizes'])}")
    print(f"  formats: {_format_counts(summary['format_versions'])}")

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
        for name in ("drh", "dparh", "pspec"):
            if name in hist:
                print(f"  {name}: {hist[name]['file']}")

    if args.spectrum_nbin is not None:
        if expected_species is None:
            raise ValueError("--spectrum-nbin requires --expected-species-counts")
        spectrum = validate_spectra(args.path, expected_species,
                                    args.spectrum_nbin)
        print("spectra:")
        print(f"  pspec: {spectrum['file']}")

    if args.joint_spectrum_nbin is not None:
        if expected_species is None:
            raise ValueError(
                "--joint-spectrum-nbin requires --expected-species-counts")
        nbin_values = _parse_counts(args.joint_spectrum_nbin)
        if nbin_values is None or len(nbin_values) != 2:
            raise ValueError("--joint-spectrum-nbin must be nbin1,nbin2")
        joint = validate_joint_spectra(
            args.path, expected_species, nbin_values[0], nbin_values[1])
        print("joint_spectra:")
        print(f"  pspec2: {joint['file']}")

    pmom_path = args.path / "pmom"
    if expected_species is not None and pmom_path.exists():
        moments = validate_moments(args.path, expected_species)
        print("moments:")
        for row in moments:
            print(
                "  species {species}: count={count} mean_mu={mean_mu:.6e} "
                "mean_mu2={mean_mu2:.6e} anisotropy={anisotropy:.6e}".format(
                    **row))

    sample_summary = summarize_samples(args.path, latest=not args.all)
    if sample_summary is not None:
        print("samples:")
        print(f"  files: {len(sample_summary['files'])}")
        print(f"  total: {sample_summary['total']}")
        print(f"  ranks: {_format_counts(sample_summary['rank_counts'])}")
        print(f"  species: {_format_counts(sample_summary['species_counts'])}")
        print(f"  fields: {','.join(sample_summary['fields'])}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
