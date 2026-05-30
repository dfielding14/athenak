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
RELATIVISTIC_OUTPUT_SCHEMA = "akcr_particle_output_v1"
DF_HISTOGRAM_LAYOUT = "species-major int32 histogram; mu=dot(v,B)/(|v||B|)"
DXH_HISTOGRAM_LAYOUT = "species-major int32 histograms ordered dx, dy, dz"
SCALAR_HISTOGRAM_LAYOUT = "species-major int32 histogram"
JOINT_HISTOGRAM_LAYOUT = "species-major int32 histogram, bin1 then bin2"
TYPED_V2_MAGIC = b"AKPRST2\0"
TYPED_V2_HEADER_BYTES = 144
TYPED_V2_HEADER_CHECKSUM_OFFSET = 120
TYPED_V2_INTEGER_FIELDS = 3
TYPED_V2_REAL_FIELDS = 22
TYPED_V2_RECORD = struct.Struct("<iii22d")
TYPED_V2_RECORD_BYTES = TYPED_V2_RECORD.size
TYPED_V2_HEADER = struct.Struct("<8sHH9IQQqdddQQQQIQ4x")
SIGNED_INT_MAX = (1 << 31) - 1


def _rank_from_path(path: Path) -> int:
    match = re.search(r"rank_(\d+)", str(path))
    return int(match.group(1)) if match else 0


def _file_number(path: Path) -> int:
    match = re.search(r"\.(\d+)\.[^.]+$", path.name)
    return int(match.group(1)) if match else -1


def _typed_v2_file_number(path: Path) -> int:
    match = re.search(r"\.(\d{5})\.[^.]+$", path.name)
    if match is None:
        raise ValueError(
            f"{path} typed-v2 restart filename is missing five-digit checkpoint number")
    return int(match.group(1))


def _typed_v2_checkpoint_id(cycle: int, mesh_time: float, mesh_dt: float,
                            file_number: int, saved_nranks: int) -> int:
    if saved_nranks <= 0 or saved_nranks > (1 << 31) - 1:
        raise ValueError("typed-v2 restart saved rank count exceeds runtime int range")
    return _fnv1a64(struct.pack(
        "<qddii", cycle, mesh_time, mesh_dt, file_number, saved_nranks))


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


def _typed_v2_manifest_path(path: Path) -> Path:
    match = re.fullmatch(r"rank_(\d+)", path.parent.name)
    if match is None:
        raise ValueError(
            f"{path} typed-v2 restart shard is not below a rank directory")
    return path.parent.parent / f"{path.name}.manifest"


def _is_typed_v2_restart(path: Path, blob: bytes) -> bool:
    if blob.startswith(TYPED_V2_MAGIC) or blob.startswith(b"AKPRST2"):
        return True
    if blob and TYPED_V2_MAGIC.startswith(blob):
        return True
    try:
        return _typed_v2_manifest_path(path).exists()
    except ValueError:
        return False


def _typed_v2_uint(text: str, label: str, maximum: int = (1 << 64) - 1,
                   context: str = "manifest") -> int:
    if re.fullmatch(r"[0-9]+", text) is None:
        raise ValueError(f"typed-v2 restart {context} has invalid {label}={text}")
    try:
        value = int(text)
    except ValueError as error:
        raise ValueError(
            f"typed-v2 restart {context} has invalid {label}={text}") from error
    if value < 0 or value > maximum:
        raise ValueError(f"typed-v2 restart {context} has invalid {label}={text}")
    return value


def _typed_v2_int(text: str, label: str, context: str = "manifest") -> int:
    if re.fullmatch(r"-?[0-9]+", text) is None:
        raise ValueError(
            f"typed-v2 restart {context} has invalid {label}={text}")
    try:
        return int(text)
    except ValueError as error:
        raise ValueError(
            f"typed-v2 restart {context} has invalid {label}={text}") from error


def _typed_v2_float(text: str, label: str, context: str = "manifest") -> float:
    if re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?",
                    text) is None:
        raise ValueError(
            f"typed-v2 restart {context} has invalid {label}={text}")
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError(
            f"typed-v2 restart {context} has invalid {label}={text}") from error
    if not math.isfinite(value):
        raise ValueError(f"typed-v2 restart {context} has non-finite {label}={text}")
    return value


def _typed_v2_safe_relative_path(path: str) -> bool:
    return (
        bool(path) and not path.startswith("/") and ".." not in path
        and "\\" not in path
    )


def _read_typed_v2_manifest(path: Path, shard_name: str) -> Dict:
    try:
        lines = path.read_text().splitlines()
    except FileNotFoundError as error:
        raise ValueError(f"typed-v2 restart requires manifest {path}") from error
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(f"typed-v2 restart manifest {path} could not be read") from error

    values: Dict[str, str] = {}
    shard_rows = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        key = fields[0]
        if key == "shard":
            if len(fields) != 7:
                raise ValueError(
                    f"typed-v2 restart manifest {path}:{line_number} "
                    "has malformed shard row")
            shard_rows.append({
                "rank": _typed_v2_uint(fields[1], "shard rank", (1 << 32) - 1),
                "relative_path": fields[2],
                "byte_count": _typed_v2_uint(fields[3], "shard byte_count"),
                "local_count": _typed_v2_uint(fields[4], "shard local_count"),
                "payload_checksum": _typed_v2_uint(
                    fields[5], "shard payload_checksum"),
                "header_checksum": _typed_v2_uint(
                    fields[6], "shard header_checksum"),
            })
            continue
        if len(fields) != 2:
            raise ValueError(
                f"typed-v2 restart manifest {path}:{line_number} "
                f"has malformed key {key!r}")
        if key in values:
            raise ValueError(f"typed-v2 restart manifest {path} repeats key {key!r}")
        values[key] = fields[1]

    def required(key: str) -> str:
        try:
            return values[key]
        except KeyError as error:
            raise ValueError(
                f"typed-v2 restart manifest {path} is missing {key!r}") from error

    if required("magic") != "AKPRST-MANIFEST":
        raise ValueError(f"typed-v2 restart manifest {path} has invalid magic")
    if required("version") != "2.0":
        raise ValueError(f"typed-v2 restart manifest {path} has unsupported version")
    if required("topology_policy") != "reject_rank_count_change":
        raise ValueError(
            f"typed-v2 restart manifest {path} has unsupported topology policy")
    if required("paired_mesh_checkpoint") != "required":
        raise ValueError(
            f"typed-v2 restart manifest {path} does not require a paired checkpoint")

    manifest = {
        "path": path,
        "values": values,
        "checkpoint_id": _typed_v2_uint(required("checkpoint_id"), "checkpoint_id"),
        "saved_nranks": _typed_v2_uint(
            required("saved_nranks"), "saved_nranks", (1 << 32) - 1),
        "global_count": _typed_v2_uint(required("global_count"), "global_count"),
        "mesh_cycle": _typed_v2_int(required("mesh_cycle"), "mesh_cycle"),
        "mesh_time": _typed_v2_float(required("mesh_time"), "mesh_time"),
        "mesh_dt": _typed_v2_float(required("mesh_dt"), "mesh_dt"),
        "particle_dtnew": _typed_v2_float(
            required("particle_dtnew"), "particle_dtnew"),
        "config_fingerprint": _typed_v2_uint(
            required("config_fingerprint"), "config_fingerprint"),
        "mesh_byte_count": _typed_v2_uint(
            required("mesh_byte_count"), "mesh_byte_count"),
        "mesh_checksum": _typed_v2_uint(required("mesh_checksum"), "mesh_checksum"),
        "mesh_topology_hash": _typed_v2_uint(
            required("mesh_topology_hash"), "mesh_topology_hash"),
        "mesh_restart": required("mesh_restart"),
        "shards": {},
    }
    if manifest["saved_nranks"] == 0:
        raise ValueError(f"typed-v2 restart manifest {path} has zero saved ranks")
    if len(shard_rows) != manifest["saved_nranks"]:
        raise ValueError(
            f"typed-v2 restart manifest {path} shard coverage count mismatch: "
            f"found {len(shard_rows)}, expected {manifest['saved_nranks']}")
    if not _typed_v2_safe_relative_path(manifest["mesh_restart"]):
        raise ValueError(
            f"typed-v2 restart manifest {path} has unsafe mesh restart path "
            f"{manifest['mesh_restart']!r}")

    manifest_total = 0
    for row in shard_rows:
        rank = row["rank"]
        relative_path = row["relative_path"]
        if rank >= manifest["saved_nranks"]:
            raise ValueError(
                f"typed-v2 restart manifest {path} shard rank {rank} is out of range")
        if not _typed_v2_safe_relative_path(relative_path):
            raise ValueError(
                f"typed-v2 restart manifest {path} has unsafe shard path "
                f"{relative_path!r}")
        expected_path = f"rank_{rank:08d}/{shard_name}"
        if relative_path != expected_path:
            raise ValueError(
                f"typed-v2 restart manifest {path} topology mismatch for rank {rank}: "
                f"path={relative_path!r}, expected={expected_path!r}")
        if rank in manifest["shards"]:
            raise ValueError(
                f"typed-v2 restart manifest {path} repeats shard rank {rank}")
        manifest["shards"][rank] = row
        if row["local_count"] > (1 << 64) - 1 - manifest_total:
            raise ValueError(
                f"typed-v2 restart manifest {path} shard count sum overflows uint64")
        manifest_total += row["local_count"]

    expected_ranks = set(range(manifest["saved_nranks"]))
    if set(manifest["shards"]) != expected_ranks:
        raise ValueError(f"typed-v2 restart manifest {path} shard coverage is incomplete")
    if manifest_total != manifest["global_count"]:
        raise ValueError(
            f"typed-v2 restart manifest {path} particle count mismatch: "
            f"shards sum to {manifest_total}, global_count={manifest['global_count']}")
    return manifest


def _read_typed_v2_mesh_witness(manifest_path: Path, manifest: Dict) -> Dict:
    run_root = manifest_path.parent.parent
    mesh_path = run_root / manifest["mesh_restart"]
    witness_path = Path(str(mesh_path) + ".rmeta")
    try:
        lines = witness_path.read_text().splitlines()
    except FileNotFoundError as error:
        raise ValueError(
            f"typed-v2 restart requires mesh witness {witness_path}") from error
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(
            f"typed-v2 restart mesh witness {witness_path} could not be read") from error

    values: Dict[str, str] = {}
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) != 2:
            raise ValueError(
                f"typed-v2 restart mesh witness {witness_path}:{line_number} "
                "has malformed row")
        key, value = fields
        if key == "shard":
            raise ValueError(
                f"typed-v2 restart mesh witness {witness_path} "
                "contains unexpected shard row")
        if key in values:
            raise ValueError(
                f"typed-v2 restart mesh witness {witness_path} "
                f"repeats key {key!r}")
        values[key] = value

    def required(key: str) -> str:
        try:
            return values[key]
        except KeyError as error:
            raise ValueError(
                f"typed-v2 restart mesh witness {witness_path} "
                f"is missing {key!r}") from error

    if required("magic") != "AKRST-WITNESS":
        raise ValueError(
            f"typed-v2 restart mesh witness {witness_path} has invalid magic")
    if required("version") != "1":
        raise ValueError(
            f"typed-v2 restart mesh witness {witness_path} has unsupported version")

    context = "mesh witness"
    witness = {
        "path": witness_path,
        "mesh_path": mesh_path,
        "checkpoint_id": _typed_v2_uint(
            required("checkpoint_id"), "checkpoint_id", context=context),
        "saved_nranks": _typed_v2_uint(
            required("saved_nranks"), "saved_nranks", (1 << 32) - 1,
            context=context),
        "mesh_cycle": _typed_v2_int(
            required("mesh_cycle"), "mesh_cycle", context=context),
        "mesh_time": _typed_v2_float(
            required("mesh_time"), "mesh_time", context=context),
        "mesh_dt": _typed_v2_float(required("mesh_dt"), "mesh_dt", context=context),
        "mesh_byte_count": _typed_v2_uint(
            required("mesh_byte_count"), "mesh_byte_count", context=context),
        "mesh_checksum": _typed_v2_uint(
            required("mesh_checksum"), "mesh_checksum", context=context),
        "mesh_topology_hash": _typed_v2_uint(
            required("mesh_topology_hash"), "mesh_topology_hash", context=context),
    }
    for key in (
            "checkpoint_id", "saved_nranks", "mesh_cycle", "mesh_time", "mesh_dt",
            "mesh_byte_count", "mesh_checksum", "mesh_topology_hash"):
        if witness[key] != manifest[key]:
            raise ValueError(
                f"typed-v2 restart mesh witness {witness_path} {key} mismatch: "
                f"witness={witness[key]}, manifest={manifest[key]}")

    try:
        mesh_blob = mesh_path.read_bytes()
    except OSError as error:
        raise ValueError(
            f"typed-v2 restart mesh checkpoint {mesh_path} could not be read") from error
    if len(mesh_blob) != witness["mesh_byte_count"]:
        raise ValueError(
            f"typed-v2 restart mesh checkpoint {mesh_path} byte count mismatch: "
            f"actual={len(mesh_blob)}, witness={witness['mesh_byte_count']}")
    actual_checksum = _fnv1a64(mesh_blob)
    if actual_checksum != witness["mesh_checksum"]:
        raise ValueError(
            f"typed-v2 restart mesh checkpoint {mesh_path} checksum mismatch: "
            f"actual={actual_checksum}, witness={witness['mesh_checksum']}")
    mesh_file_number = _typed_v2_file_number(mesh_path)
    expected_checkpoint_id = _typed_v2_checkpoint_id(
        witness["mesh_cycle"], witness["mesh_time"], witness["mesh_dt"],
        mesh_file_number, witness["saved_nranks"])
    if witness["checkpoint_id"] != expected_checkpoint_id:
        raise ValueError(
            f"typed-v2 restart mesh witness {witness_path} checkpoint ID is inconsistent")
    witness["file_number"] = mesh_file_number
    return witness


def _decode_typed_v2_header(path: Path, blob: bytes) -> Dict:
    if len(blob) < TYPED_V2_HEADER_BYTES:
        raise ValueError(
            f"{path} typed-v2 restart has truncated header: "
            f"byte_count={len(blob)}, expected at least {TYPED_V2_HEADER_BYTES}")

    fields = TYPED_V2_HEADER.unpack_from(blob)
    (magic, schema_major, schema_minor, header_bytes, flags, endian_marker,
     integer_fields, real_fields, record_bytes, saved_nranks, saved_rank,
     topology_policy, local_count, global_count, mesh_cycle, mesh_time, mesh_dt,
     particle_dtnew, checkpoint_id, payload_bytes, payload_checksum,
     header_checksum, requires_manifest, config_fingerprint) = fields
    if magic != TYPED_V2_MAGIC:
        raise ValueError(f"{path} typed-v2 restart has corrupt magic")
    if (schema_major, schema_minor) != (2, 0):
        raise ValueError(
            f"{path} typed-v2 restart schema version mismatch: "
            f"found {schema_major}.{schema_minor}, expected 2.0")
    if header_bytes != TYPED_V2_HEADER_BYTES:
        raise ValueError(
            f"{path} typed-v2 restart header width mismatch: "
            f"found {header_bytes}, expected {TYPED_V2_HEADER_BYTES}")
    if flags != 1:
        raise ValueError(f"{path} typed-v2 restart mode flags mismatch: found {flags}")
    if endian_marker != 0x01020304:
        raise ValueError(f"{path} typed-v2 restart endian marker mismatch")
    if (integer_fields, real_fields, record_bytes) != (
            TYPED_V2_INTEGER_FIELDS, TYPED_V2_REAL_FIELDS, TYPED_V2_RECORD_BYTES):
        raise ValueError(
            f"{path} typed-v2 restart record type mismatch: "
            f"found int32_fields={integer_fields}, f64_fields={real_fields}, "
            f"record_bytes={record_bytes}; expected int32_fields=3, "
            f"f64_fields=22, record_bytes={TYPED_V2_RECORD_BYTES}")
    if topology_policy != 1:
        raise ValueError(f"{path} typed-v2 restart topology policy mismatch")
    if requires_manifest != 1:
        raise ValueError(f"{path} typed-v2 restart does not require a manifest")

    checksum_bytes = bytearray(blob[:TYPED_V2_HEADER_BYTES])
    checksum_bytes[
        TYPED_V2_HEADER_CHECKSUM_OFFSET:TYPED_V2_HEADER_CHECKSUM_OFFSET + 8
    ] = b"\0" * 8
    actual_header_checksum = _fnv1a64(checksum_bytes)
    if actual_header_checksum != header_checksum:
        raise ValueError(
            f"{path} typed-v2 restart header checksum mismatch: "
            f"stored={header_checksum}, actual={actual_header_checksum}")
    if not all(math.isfinite(value) for value in (mesh_time, mesh_dt, particle_dtnew)):
        raise ValueError(f"{path} typed-v2 restart header has non-finite timestep data")
    if saved_nranks == 0 or saved_rank >= saved_nranks:
        raise ValueError(f"{path} typed-v2 restart has invalid saved rank topology")
    if local_count > (1 << 31) - 1:
        raise ValueError(f"{path} typed-v2 restart local count exceeds runtime int range")
    expected_payload_bytes = local_count * TYPED_V2_RECORD_BYTES
    if payload_bytes != expected_payload_bytes:
        raise ValueError(
            f"{path} typed-v2 restart payload count mismatch: "
            f"payload_bytes={payload_bytes}, expected {expected_payload_bytes} "
            f"for local_count={local_count}")
    expected_bytes = TYPED_V2_HEADER_BYTES + payload_bytes
    if len(blob) < expected_bytes:
        raise ValueError(
            f"{path} typed-v2 restart has truncated payload: "
            f"byte_count={len(blob)}, expected {expected_bytes}")
    if len(blob) > expected_bytes:
        raise ValueError(
            f"{path} typed-v2 restart has corrupt trailing data: "
            f"byte_count={len(blob)}, expected {expected_bytes}")
    actual_payload_checksum = _fnv1a64(blob[TYPED_V2_HEADER_BYTES:])
    if actual_payload_checksum != payload_checksum:
        raise ValueError(
            f"{path} typed-v2 restart payload checksum mismatch: "
            f"stored={payload_checksum}, actual={actual_payload_checksum}")
    return {
        "saved_nranks": saved_nranks,
        "saved_rank": saved_rank,
        "local_count": local_count,
        "global_count": global_count,
        "mesh_cycle": mesh_cycle,
        "mesh_time": mesh_time,
        "mesh_dt": mesh_dt,
        "particle_dtnew": particle_dtnew,
        "checkpoint_id": checkpoint_id,
        "payload_bytes": payload_bytes,
        "payload_checksum": payload_checksum,
        "header_checksum": header_checksum,
        "config_fingerprint": config_fingerprint,
    }


def _validate_typed_v2_header_manifest(path: Path, header: Dict,
                                       manifest: Dict, shard: Dict) -> None:
    checks = {
        "saved_nranks": manifest["saved_nranks"],
        "saved_rank": shard["rank"],
        "local_count": shard["local_count"],
        "global_count": manifest["global_count"],
        "mesh_cycle": manifest["mesh_cycle"],
        "mesh_time": manifest["mesh_time"],
        "mesh_dt": manifest["mesh_dt"],
        "particle_dtnew": manifest["particle_dtnew"],
        "checkpoint_id": manifest["checkpoint_id"],
        "config_fingerprint": manifest["config_fingerprint"],
        "payload_checksum": shard["payload_checksum"],
        "header_checksum": shard["header_checksum"],
    }
    for key, expected in checks.items():
        actual = header[key]
        if actual != expected:
            raise ValueError(
                f"{path} typed-v2 restart {key} mismatch: "
                f"header={actual}, manifest={expected}")


def _read_typed_v2_checkpoint(path: Path, selected_blob: bytes) -> Dict:
    particle_file_number = _typed_v2_file_number(path)
    manifest_path = _typed_v2_manifest_path(path)
    manifest = _read_typed_v2_manifest(manifest_path, path.name)
    mesh_witness = _read_typed_v2_mesh_witness(manifest_path, manifest)
    if particle_file_number != mesh_witness["file_number"]:
        raise ValueError(
            f"{path} typed-v2 restart particle and mesh checkpoint numbers do not match")
    shards = {}
    selected_path = path.resolve()
    for shard in manifest["shards"].values():
        shard_path = manifest_path.parent / shard["relative_path"]
        if _typed_v2_file_number(shard_path) != mesh_witness["file_number"]:
            raise ValueError(
                f"{shard_path} typed-v2 restart particle and mesh checkpoint "
                "numbers do not match")
        try:
            blob = (selected_blob if shard_path.resolve() == selected_path
                    else shard_path.read_bytes())
        except OSError as error:
            raise ValueError(
                f"typed-v2 restart manifest {manifest_path} references missing "
                f"shard {shard_path}") from error
        if len(blob) != shard["byte_count"]:
            condition = "truncated" if len(blob) < shard["byte_count"] else "corrupt"
            raise ValueError(
                f"{shard_path} typed-v2 restart has {condition} shard byte count: "
                f"actual={len(blob)}, manifest={shard['byte_count']}")
        header = _decode_typed_v2_header(shard_path, blob)
        _validate_typed_v2_header_manifest(shard_path, header, manifest, shard)
        shards[shard_path.resolve()] = {"blob": blob, "header": header}
    if selected_path not in shards:
        raise ValueError(
            f"{path} typed-v2 restart shard is not covered by manifest {manifest_path}")
    return {"manifest": manifest, "mesh_witness": mesh_witness, "shards": shards}


def _decode_typed_v2_restart(path: Path, checkpoint: Dict) -> Dict:
    shard = checkpoint["shards"][path.resolve()]
    header = shard["header"]
    particles = []
    payload = shard["blob"][TYPED_V2_HEADER_BYTES:]
    for record in TYPED_V2_RECORD.iter_unpack(payload):
        pgid, tag, species = record[:TYPED_V2_INTEGER_FIELDS]
        reals = tuple(float(value) for value in record[TYPED_V2_INTEGER_FIELDS:])
        if not all(math.isfinite(value) for value in reals):
            raise ValueError(f"{path} typed-v2 restart contains non-finite particle data")
        particles.append({
            "pgid": pgid,
            "tag": tag,
            "species": species,
            "reals": reals,
        })
    if len(particles) != header["local_count"]:
        raise ValueError(
            f"{path} typed-v2 restart particle count mismatch: "
            f"decoded={len(particles)}, header={header['local_count']}")

    manifest = checkpoint["manifest"]
    metadata = dict(manifest["values"])
    metadata["path"] = str(manifest["path"])
    return {
        "path": path,
        "rank": header["saved_rank"],
        "time": header["mesh_time"],
        "dt": header["mesh_dt"],
        "count": header["local_count"],
        "real_size": struct.calcsize("<d"),
        "format_version": "AKPRST-v2.0",
        "metadata": metadata,
        "particles": particles,
    }


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


def read_prst_file(path: Path,
                   _typed_v2_checkpoints: Optional[Dict[Path, Dict]] = None) -> Dict:
    """Read one per-rank particle restart file and validate its length."""
    blob = path.read_bytes()
    if _is_typed_v2_restart(path, blob):
        manifest_path = _typed_v2_manifest_path(path).resolve()
        checkpoint = None
        if _typed_v2_checkpoints is not None:
            checkpoint = _typed_v2_checkpoints.get(manifest_path)
        if checkpoint is None:
            checkpoint = _read_typed_v2_checkpoint(path, blob)
            if _typed_v2_checkpoints is not None:
                _typed_v2_checkpoints[manifest_path] = checkpoint
        return _decode_typed_v2_restart(path, checkpoint)

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
    header_count = _header_nonnegative_int(path, header, "particles")
    metadata: Dict[str, str] = {}
    start = header_end + 1
    while start < len(blob) and blob[start:start + 1] == b"#":
        line_end = blob.find(b"\n", start)
        if line_end < 0:
            raise ValueError(f"{path} has a malformed ppd metadata line")
        try:
            comment = blob[start:line_end].decode("ascii")
        except UnicodeDecodeError as error:
            raise ValueError(f"{path} has a non-ASCII ppd metadata line") from error
        if comment.startswith("# metadata:"):
            parsed = _parse_metadata_line(comment.split(":", 1)[1])
            duplicates = metadata.keys() & parsed.keys()
            if duplicates:
                raise ValueError(
                    f"{path} repeats ppd metadata keys {sorted(duplicates)}")
            metadata.update(parsed)
        start = line_end + 1
    _validate_relativistic_metadata(
        path,
        metadata,
        "ppd",
        ("columns", "position_units", "payload"),
        {
            "columns": "species_x_y_z",
            "position_units": "code_length",
            "payload": "float32",
        },
    )
    payload = blob[start:]
    record_size = PPD_FIELDS * struct.calcsize("<f")
    if len(payload) % record_size != 0:
        raise ValueError(f"{path} has malformed ppd payload length")

    particles = []
    for record in struct.iter_unpack("<ffff", payload):
        species_value = float(record[0])
        if (not math.isfinite(species_value) or species_value < 0.0
                or species_value > SIGNED_INT_MAX
                or not species_value.is_integer()):
            raise ValueError(
                f"{path} contains non-integral, negative, or out-of-range "
                "ppd species")
        species = int(species_value)
        xyz = tuple(float(value) for value in record[1:])
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"{path} contains non-finite positions")
        particles.append({"species": species, "xyz": xyz})
    if len(particles) != header_count:
        raise ValueError(
            f"{path} ppd header particles={header_count}, "
            f"payload count={len(particles)}")

    return {
        "path": path,
        "header": header,
        "metadata": metadata,
        "count": len(particles),
        "species_counts": Counter(p["species"] for p in particles),
        "particles": particles,
    }


def _histogram_record_details(path: Path, marker: bytes, count: int,
                              skip_comment_lines: bool = False) -> Dict:
    blob = path.read_bytes()
    index = blob.rfind(marker)
    if index < 0:
        raise ValueError(f"{path} does not contain a recognizable histogram header")

    line_end = blob.find(b"\n", index)
    if line_end < 0:
        raise ValueError(f"{path} has a malformed histogram header")

    header = blob[index:line_end].decode("ascii", errors="replace")
    metadata: Dict[str, str] = {}
    layout: Optional[str] = None
    byte_count = count * struct.calcsize("<i")
    start = len(blob) - byte_count
    if start < line_end + 1:
        raise ValueError(f"{path} has an incomplete histogram record")
    if skip_comment_lines:
        try:
            prefix = blob[line_end + 1:start].decode("ascii")
        except UnicodeDecodeError as error:
            raise ValueError(
                f"{path} has a non-ASCII histogram metadata prefix") from error
        for comment in prefix.splitlines():
            if not comment:
                continue
            if comment.startswith("# metadata:"):
                parsed = _parse_metadata_line(comment.split(":", 1)[1])
                duplicates = metadata.keys() & parsed.keys()
                if duplicates:
                    raise ValueError(
                        f"{path} repeats histogram metadata keys "
                        f"{sorted(duplicates)}")
                metadata.update(parsed)
            elif comment.startswith("# layout:"):
                if layout is not None:
                    raise ValueError(f"{path} repeats histogram layout metadata")
                layout = comment.split(":", 1)[1].strip()
            else:
                raise ValueError(
                    f"{path} has an unrecognized histogram metadata prefix line "
                    f"{comment!r}")
    end = start + byte_count
    values = list(struct.unpack("<" + str(count) + "i", blob[start:end]))
    if any(value < 0 for value in values):
        raise ValueError(f"{path} contains negative histogram bin count")
    return {
        "header": header,
        "metadata": metadata,
        "layout": layout,
        "values": values,
    }


def _histogram_record(path: Path, marker: bytes, count: int,
                      skip_comment_lines: bool = False) -> List[int]:
    return _histogram_record_details(
        path, marker, count, skip_comment_lines=skip_comment_lines)["values"]


def _validate_scalar_histogram_metadata(path: Path, record: Dict,
                                        nspecies: int, nbin: int) -> None:
    for name, expected in (("nspecies", nspecies), ("nbin", nbin)):
        if name in record["metadata"]:
            actual = _metadata_nonnegative_int(path, record["metadata"], name)
            if actual != expected:
                raise ValueError(
                    f"{path} metadata {name}={actual}, expected {expected}")


def _validate_histogram_layout(path: Path, record: Dict,
                               expected: str) -> None:
    layout = record["layout"]
    if layout is not None and layout != expected:
        raise ValueError(
            f"{path} histogram layout {layout!r}, expected {expected!r}")
    if record["metadata"].get("mode") == "relativistic_hc" and layout is None:
        raise ValueError(
            f"{path} relativistic_hc histogram is missing layout metadata")


def read_df_record(path: Path, nspecies: int, nbin: int) -> Dict:
    record = _histogram_record_details(
        path, b"# AthenaK particle distribution function", nspecies * nbin,
        skip_comment_lines=True)
    _validate_histogram_layout(path, record, DF_HISTOGRAM_LAYOUT)
    _validate_scalar_histogram_metadata(path, record, nspecies, nbin)
    _validate_relativistic_scalar_histogram_metadata(
        path, record["metadata"], "df", "mu", "dimensionless", nspecies, nbin)
    values = record.pop("values")
    record["histogram"] = [
        values[sp * nbin: (sp + 1) * nbin] for sp in range(nspecies)]
    return record


def read_df_file(path: Path, nspecies: int, nbin: int) -> List[List[int]]:
    return read_df_record(path, nspecies, nbin)["histogram"]


def read_dxh_record(path: Path, nspecies: int, nbin: int) -> Dict:
    record = _histogram_record_details(
        path, b"# AthenaK particle displacement histogram", 3 * nspecies * nbin,
        skip_comment_lines=True)
    _validate_histogram_layout(path, record, DXH_HISTOGRAM_LAYOUT)
    _validate_scalar_histogram_metadata(path, record, nspecies, nbin)
    _validate_relativistic_scalar_histogram_metadata(
        path, record["metadata"], "dxh", "dx_dy_dz", "code_length",
        nspecies, nbin)
    values = record.pop("values")
    result = []
    for sp in range(nspecies):
        species_values = []
        base = sp * 3 * nbin
        for component in range(3):
            start = base + component * nbin
            species_values.append(values[start: start + nbin])
        result.append(species_values)
    record["histogram"] = result
    return record


def read_dxh_file(path: Path, nspecies: int, nbin: int) -> List[List[List[int]]]:
    return read_dxh_record(path, nspecies, nbin)["histogram"]


def read_scalar_histogram_record(path: Path, marker: bytes, nspecies: int,
                                 nbin: int) -> Dict:
    record = _histogram_record_details(
        path, marker, nspecies * nbin, skip_comment_lines=True)
    _validate_histogram_layout(path, record, SCALAR_HISTOGRAM_LAYOUT)
    _validate_scalar_histogram_metadata(path, record, nspecies, nbin)
    if marker == b"# AthenaK particle scalar displacement histogram":
        _validate_relativistic_scalar_histogram_metadata(
            path, record["metadata"], "drh", "displacement_norm", "code_length",
            nspecies, nbin)
    elif marker == b"# AthenaK particle parallel displacement histogram":
        _validate_relativistic_scalar_histogram_metadata(
            path,
            record["metadata"],
            "dparh",
            "dparallel",
            "code_length",
            nspecies,
            nbin,
            {"definition": "accumulated_midpoint_sum_dx_dot_Bhat"},
        )
    values = record.pop("values")
    record["histogram"] = [
        values[sp * nbin: (sp + 1) * nbin] for sp in range(nspecies)]
    return record


def read_scalar_histogram(path: Path, marker: bytes, nspecies: int,
                          nbin: int) -> List[List[int]]:
    return read_scalar_histogram_record(path, marker, nspecies, nbin)["histogram"]


def read_pspec_record(path: Path, nspecies: int, nbin: int) -> Dict:
    record = _histogram_record_details(
        path, b"# AthenaK particle spectrum", nspecies * nbin,
        skip_comment_lines=True)
    _validate_histogram_layout(path, record, SCALAR_HISTOGRAM_LAYOUT)
    for name, expected in (("nspecies", nspecies), ("nbin", nbin)):
        if name in record["metadata"]:
            actual = _metadata_nonnegative_int(path, record["metadata"], name)
            if actual != expected:
                raise ValueError(
                    f"{path} metadata {name}={actual}, expected {expected}")
    if _declares_relativistic_metadata(record["metadata"]):
        quantity = record["metadata"].get("quantity")
        if quantity not in RELATIVISTIC_PSPEC_QUANTITIES:
            raise ValueError(
                f"{path} relativistic_hc spectrum rejects ambiguous or "
                f"missing quantity {quantity!r}")
        required = []
        if quantity == "log10_kinetic_energy_model":
            required.extend(("log10_floor", "log10_floor_units"))
        _validate_relativistic_histogram_metadata(
            path,
            record["metadata"],
            "pspec",
            ("quantity", "quantity_units", "histogram_units", "nspecies", "nbin",
             "vmin", "vmax", "reduce", *required),
            (("nspecies", nspecies), ("nbin", nbin)),
            (("vmin", "vmax"),),
            {
                "quantity_units": RELATIVISTIC_PSPEC_UNITS[quantity],
                "histogram_units": "particle_count",
                **({
                    "log10_floor_units": "code_velocity_squared",
                } if quantity == "log10_kinetic_energy_model" else {}),
            },
        )
        if quantity == "log10_kinetic_energy_model":
            log10_floor = _metadata_finite_float(
                path, record["metadata"], "log10_floor")
            if log10_floor <= 0.0:
                raise ValueError(f"{path} metadata requires log10_floor > 0")
    values = record.pop("values")
    record["histogram"] = [
        values[sp * nbin: (sp + 1) * nbin] for sp in range(nspecies)]
    return record


def read_pspec_file(path: Path, nspecies: int, nbin: int) -> List[List[int]]:
    return read_pspec_record(path, nspecies, nbin)["histogram"]


def read_pspec2_record(path: Path, nspecies: int,
                       nbin1: int, nbin2: int) -> Dict:
    record = _histogram_record_details(
        path, b"# AthenaK particle joint spectrum", nspecies * nbin1 * nbin2,
        skip_comment_lines=True)
    _validate_histogram_layout(path, record, JOINT_HISTOGRAM_LAYOUT)
    for name, expected in (
            ("nspecies", nspecies), ("nbin1", nbin1), ("nbin2", nbin2)):
        if name in record["metadata"]:
            actual = _metadata_nonnegative_int(path, record["metadata"], name)
            if actual != expected:
                raise ValueError(
                    f"{path} metadata {name}={actual}, expected {expected}")
    if _declares_relativistic_metadata(record["metadata"]):
        quantity = record["metadata"].get("quantity")
        if quantity not in RELATIVISTIC_PSPEC2_QUANTITIES:
            raise ValueError(
                f"{path} relativistic_hc joint spectrum rejects ambiguous or "
                f"missing quantity {quantity!r}")
        axis1_units, axis2_units = RELATIVISTIC_PSPEC2_UNITS[quantity]
        _validate_relativistic_histogram_metadata(
            path,
            record["metadata"],
            "pspec2",
            ("quantity", "axis1_units", "axis2_units", "histogram_units",
             "nspecies", "nbin1", "nbin2", "vmin1", "vmax1", "vmin2",
             "vmax2", "reduce"),
            (("nspecies", nspecies), ("nbin1", nbin1), ("nbin2", nbin2)),
            (("vmin1", "vmax1"), ("vmin2", "vmax2")),
            {
                "axis1_units": axis1_units,
                "axis2_units": axis2_units,
                "histogram_units": "particle_count",
            },
        )
    values = record.pop("values")
    result = []
    for species in range(nspecies):
        species_values = []
        base = species * nbin1 * nbin2
        for bin1 in range(nbin1):
            start = base + bin1 * nbin2
            species_values.append(values[start: start + nbin2])
        result.append(species_values)
    record["histogram"] = result
    return record


def read_pspec2_file(path: Path, nspecies: int,
                     nbin1: int, nbin2: int) -> List[List[List[int]]]:
    return read_pspec2_record(path, nspecies, nbin1, nbin2)["histogram"]


def read_pmom_record(path: Path, nspecies: int) -> Dict:
    """Read the latest per-species particle moment block and its metadata."""
    if path.is_dir():
        path = max(path.glob("*.pmom"), key=_file_number)
    lines = path.read_text().splitlines()
    starts = [
        index for index, line in enumerate(lines)
        if line.startswith("# AthenaK particle moments")
    ]
    if not starts:
        raise ValueError(f"{path} is missing particle-moment header")
    start = starts[-1]
    metadata: Dict[str, str] = {}
    rows = []
    for line in lines[start + 1:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("# metadata:"):
            parsed = _parse_metadata_line(stripped.split(":", 1)[1])
            duplicates = metadata.keys() & parsed.keys()
            if duplicates:
                raise ValueError(
                    f"{path} repeats pmom metadata keys {sorted(duplicates)}")
            metadata.update(parsed)
            continue
        if stripped.startswith("#"):
            continue
        values = stripped.split()
        if len(values) not in (14, 26):
            raise ValueError(f"{path} has malformed particle-moment row: {line}")
        rows.append(values)
    if len(rows) != nspecies:
        raise ValueError(
            f"{path} latest block contains {len(rows)} particle-moment rows, "
            f"expected {nspecies}")

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
    for row in rows:
        parsed = {name: float(value) for name, value in zip(names, row)}
        if not all(math.isfinite(value) for value in parsed.values()):
            raise ValueError(f"{path} contains non-finite particle moments")
        for name in ("species", "count"):
            if parsed[name] < 0.0 or not parsed[name].is_integer():
                raise ValueError(
                    f"{path} contains non-integral or negative {name}")
            parsed[name] = int(parsed[name])
        result.append(parsed)
    modern = _validate_relativistic_metadata(
        path,
        metadata,
        "pmom",
        (
            "quantity_basis",
            "mu_units",
            "displacement_units",
            "displacement_second_moment_units",
            "velocity_units",
            "velocity_second_moment_units",
            "dparallel_definition",
        ),
        {
            "quantity_basis": "physical_velocity_shadow",
            "mu_units": "dimensionless",
            "displacement_units": "code_length",
            "displacement_second_moment_units": "code_length_squared",
            "velocity_units": "code_velocity",
            "velocity_second_moment_units": "code_velocity_squared",
            "dparallel_definition": "accumulated_midpoint_sum_dx_dot_Bhat",
        },
    )
    if modern and any(len(row) != 26 for row in rows):
        raise ValueError(f"{path} relativistic_hc pmom requires 26-column rows")
    if modern:
        for expected_species, moment in enumerate(result):
            if moment["species"] != expected_species:
                raise ValueError(
                    f"{path} relativistic_hc pmom row {expected_species} has "
                    f"species {moment['species']}")
    return {"path": path, "metadata": metadata, "moments": result}


def read_pmom_file(path: Path, nspecies: int) -> List[Dict[str, float]]:
    """Read the latest per-species particle moment block."""
    return read_pmom_record(path, nspecies)["moments"]


def _parse_metadata_line(line: str) -> Dict[str, str]:
    result = {}
    for token in line.split():
        if "=" not in token:
            raise ValueError(f"Malformed metadata token {token!r}")
        key, value = token.split("=", 1)
        if not key or not value or key in result:
            raise ValueError(f"Malformed or duplicate metadata key {key!r}")
        result[key] = value
    return result


def _metadata_nonnegative_int(path: Path, metadata: Dict[str, str], key: str) -> int:
    value = metadata[key]
    if re.fullmatch(r"[0-9]+", value) is None:
        raise ValueError(f"{path} metadata has invalid {key}={value!r}")
    return int(value)


def _header_nonnegative_int(path: Path, header: str, key: str,
                            minimum: int = 0,
                            maximum: int = SIGNED_INT_MAX) -> int:
    values = re.findall(
        rf"(?:^|\s){re.escape(key)}\s*=\s*(\S+)", header)
    if len(values) != 1:
        raise ValueError(
            f"{path} header requires exactly one {key} integer token")
    value = values[0]
    if re.fullmatch(r"[0-9]+", value) is None:
        raise ValueError(f"{path} header has invalid {key}={value!r}")
    result = int(value)
    if result < minimum or result > maximum:
        raise ValueError(f"{path} header has out-of-range {key}={value!r}")
    return result


def _metadata_signed_int(path: Path, metadata: Dict[str, str], key: str,
                         minimum: int = -SIGNED_INT_MAX - 1,
                         maximum: int = SIGNED_INT_MAX) -> int:
    value = metadata[key]
    if re.fullmatch(r"-?[0-9]+", value) is None:
        raise ValueError(f"{path} metadata has invalid {key}={value!r}")
    result = int(value)
    if result < minimum or result > maximum:
        raise ValueError(f"{path} metadata has out-of-range {key}={value!r}")
    return result


def _metadata_finite_float(path: Path, metadata: Dict[str, str], key: str) -> float:
    try:
        value = float(metadata[key])
    except ValueError as error:
        raise ValueError(
            f"{path} metadata has invalid {key}={metadata[key]!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"{path} metadata has non-finite {key}={metadata[key]!r}")
    return value


def _validate_relativistic_metadata(
    path: Path,
    metadata: Dict[str, str],
    family: str,
    required: Sequence[str] = (),
    expected: Optional[Dict[str, str]] = None,
) -> bool:
    modern = _declares_relativistic_metadata(metadata)
    if not modern:
        return False
    expected_values = {
        "schema": RELATIVISTIC_OUTPUT_SCHEMA,
        "mode": "relativistic_hc",
        "units": "code_model",
        **(expected or {}),
    }
    missing = sorted({
        key for key in (*expected_values, "c_model", "alpha_s", *required)
        if key not in metadata
    })
    if missing:
        raise ValueError(
            f"{path} relativistic_hc {family} metadata is missing {missing}")
    for key, value in expected_values.items():
        if metadata[key] != value:
            raise ValueError(
                f"{path} relativistic_hc {family} metadata {key}="
                f"{metadata[key]!r}, expected {value!r}")
    c_model = _metadata_finite_float(path, metadata, "c_model")
    alpha_s = _metadata_finite_float(path, metadata, "alpha_s")
    if c_model <= 0.0 or alpha_s == 0.0:
        raise ValueError(
            f"{path} relativistic_hc {family} metadata requires c_model > 0 "
            "and alpha_s != 0")
    return True


def _declares_relativistic_metadata(metadata: Dict[str, str]) -> bool:
    return "mode" in metadata or "schema" in metadata


def _validate_relativistic_histogram_metadata(
    path: Path,
    metadata: Dict[str, str],
    family: str,
    required: Sequence[str],
    integer_fields: Sequence[tuple[str, int]],
    range_fields: Sequence[tuple[str, str]],
    expected: Optional[Dict[str, str]] = None,
) -> None:
    if not _validate_relativistic_metadata(
            path, metadata, family, required, expected):
        return
    for key, value in integer_fields:
        actual = _metadata_nonnegative_int(path, metadata, key)
        if actual != value:
            raise ValueError(f"{path} metadata {key}={actual}, expected {value}")
    reduce = _metadata_nonnegative_int(path, metadata, "reduce")
    if reduce not in (0, 1):
        raise ValueError(f"{path} metadata reduce={reduce}, expected 0 or 1")
    for lower, upper in range_fields:
        lower_value = _metadata_finite_float(path, metadata, lower)
        upper_value = _metadata_finite_float(path, metadata, upper)
        if lower_value >= upper_value:
            raise ValueError(
                f"{path} metadata requires {lower} < {upper}")


def _validate_relativistic_scalar_histogram_metadata(
    path: Path,
    metadata: Dict[str, str],
    family: str,
    quantity: str,
    quantity_units: str,
    nspecies: int,
    nbin: int,
    additions: Optional[Dict[str, str]] = None,
) -> None:
    _validate_relativistic_histogram_metadata(
        path,
        metadata,
        family,
        ("quantity", "quantity_units", "histogram_units", "nspecies", "nbin",
         "vmin", "vmax", "reduce", *(additions or {})),
        (("nspecies", nspecies), ("nbin", nbin)),
        (("vmin", "vmax"),),
        {
            "quantity": quantity,
            "quantity_units": quantity_units,
            "histogram_units": "particle_count",
            **(additions or {}),
        },
    )


PSAMP_COLUMN_ALIASES = {
    "rank": "rank",
    "pgid": "pgid",
    "gid": "pgid",
    "tag": "tag",
    "species": "species",
    "x": "x",
    "x1": "x",
    "y": "y",
    "x2": "y",
    "z": "z",
    "x3": "z",
    "vx": "vx",
    "v1": "vx",
    "vy": "vy",
    "v2": "vy",
    "vz": "vz",
    "v3": "vz",
    "m": "mass",
    "mass": "mass",
    "bx": "bx",
    "b1": "bx",
    "by": "by",
    "b2": "by",
    "bz": "bz",
    "b3": "bz",
    "dx": "dx",
    "dy": "dy",
    "dz": "dz",
    "dpar": "dpar",
    "dparallel": "dpar",
    "speed": "speed",
    "p": "speed",
    "energy": "energy",
    "e": "energy",
    "bmag": "bmag",
    "mu": "mu",
    "vpar": "vpar",
    "vparallel": "vpar",
    "vperp": "vperp",
    "magnetic_moment": "magnetic_moment",
    "velocity_magnetic_moment_proxy": "velocity_magnetic_moment_proxy",
    "wx": "wx",
    "w1": "wx",
    "wy": "wy",
    "w2": "wy",
    "wz": "wz",
    "w3": "wz",
    "wmag": "wmag",
    "gamma": "gamma",
    "kinetic_energy_model": "kinetic_energy_model",
    "cex": "cex",
    "ce1": "cex",
    "cey": "cey",
    "ce2": "cey",
    "cez": "cez",
    "ce3": "cez",
    "work": "work",
    "alpha": "alpha_s",
    "alpha_s": "alpha_s",
    "r_larmor_over_dx_min": "r_larmor_over_dx_min",
}

RELATIVISTIC_PSAMP_REJECTED_ALIASES = {
    "p", "energy", "e", "mass", "m", "magnetic_moment",
}
RELATIVISTIC_PSPEC_QUANTITIES = {
    "speed", "wmag", "kinetic_energy_model", "log10_kinetic_energy_model",
    "mu", "velocity_magnetic_moment_proxy",
}
RELATIVISTIC_PSPEC_UNITS = {
    "speed": "code_velocity",
    "wmag": "code_velocity",
    "kinetic_energy_model": "code_velocity_squared",
    "log10_kinetic_energy_model": "log10_code_velocity_squared",
    "mu": "dimensionless",
    "velocity_magnetic_moment_proxy": "code_velocity_squared_per_code_B",
}
RELATIVISTIC_PSPEC2_QUANTITIES = {
    "mu_speed", "mu_wmag", "mu_kinetic_energy_model", "vpar_vperp",
}
RELATIVISTIC_PSPEC2_UNITS = {
    "mu_speed": ("dimensionless", "code_velocity"),
    "mu_wmag": ("dimensionless", "code_velocity"),
    "mu_kinetic_energy_model": ("dimensionless", "code_velocity_squared"),
    "vpar_vperp": ("code_velocity", "code_velocity"),
}


def _canonical_psamp_columns(
        path: Path, line_number: int, columns: Sequence[str],
        metadata: Dict[str, str]) -> List[str]:
    canonical = []
    for name in columns:
        normalized = name.strip().lower()
        if (metadata.get("mode") == "relativistic_hc" and
                normalized in RELATIVISTIC_PSAMP_REJECTED_ALIASES):
            raise ValueError(
                f"{path}:{line_number} relativistic_hc sample rejects "
                f"ambiguous column {name!r}")
        if normalized not in PSAMP_COLUMN_ALIASES:
            raise ValueError(
                f"{path}:{line_number} has unknown sample column {name!r}")
        canonical.append(PSAMP_COLUMN_ALIASES[normalized])
    duplicates = sorted(
        name for name, count in Counter(canonical).items() if count > 1)
    if duplicates:
        raise ValueError(
            f"{path}:{line_number} has duplicate canonical sample columns "
            f"{duplicates}")
    if len(canonical) < 4 or canonical[:4] != ["rank", "pgid", "tag", "species"]:
        raise ValueError(f"{path}:{line_number} has malformed columns")
    return canonical


def _modern_psamp_rank(path: Path) -> int:
    match = re.fullmatch(r"rank_([0-9]{8})", path.parent.name)
    if match is None or path.parent.parent.name != "psamp":
        raise ValueError(
            f"{path} relativistic_hc sample requires canonical "
            "psamp/rank_######## shard path")
    return int(match.group(1))


def _psamp_header_nranks(path: Path, header: str) -> Optional[int]:
    if re.search(r"(?:^|\s)nranks\s*=", header) is None:
        return None
    return _header_nonnegative_int(path, header, "nranks", minimum=1)


def _validate_psamp_block(path: Path, block: Dict) -> None:
    columns = block["columns"]
    if not columns:
        raise ValueError(f"{path} particle sample block is missing columns")
    metadata = block["metadata"]
    if "field_count" in metadata:
        field_count = _metadata_nonnegative_int(path, metadata, "field_count")
        if field_count != len(columns) - 4:
            raise ValueError(
                f"{path} metadata field_count={field_count}, expected "
                f"{len(columns) - 4}")
    if "sample_count" in metadata:
        sample_count = _metadata_nonnegative_int(path, metadata, "sample_count")
        if sample_count != len(block["rows"]):
            raise ValueError(
                f"{path} metadata sample_count={sample_count}, expected "
                f"{len(block['rows'])}")
    modern = _validate_relativistic_metadata(
        path,
        metadata,
        "psamp",
        (
            "selected_species",
            "sample_stride",
            "sample_offset",
            "sample_count",
            "field_count",
            "position_units",
            "velocity_units",
            "w_units",
            "cE_units",
            "gamma_units",
            "kinetic_energy_model_units",
            "work_units",
            "alpha_s_units",
            "r_larmor_over_dx_min_units",
        ),
        {
            "position_units": "code_length",
            "velocity_units": "code_velocity",
            "w_units": "code_velocity",
            "cE_units": "code_velocity_times_code_B",
            "gamma_units": "dimensionless",
            "kinetic_energy_model_units": "code_velocity_squared",
            "work_units": "code_velocity_squared",
            "alpha_s_units": "inverse_code_time_per_code_B",
            "r_larmor_over_dx_min_units": "dimensionless",
        },
    )
    if modern:
        rank = _modern_psamp_rank(path)
        nranks = block["nranks"]
        if nranks is None:
            raise ValueError(
                f"{path} relativistic_hc sample header is missing nranks")
        if rank >= nranks:
            raise ValueError(
                f"{path} relativistic_hc sample shard rank {rank} is outside "
                f"nranks={nranks}")
        selected_species = _metadata_signed_int(
            path, metadata, "selected_species", minimum=-1)
        sample_stride = _metadata_signed_int(
            path, metadata, "sample_stride", minimum=1)
        sample_offset = _metadata_signed_int(
            path, metadata, "sample_offset", minimum=0)
        if sample_offset >= sample_stride:
            raise ValueError(
                f"{path} relativistic_hc sample requires "
                "0 <= sample_offset < sample_stride")
        for row in block["rows"]:
            for name in ("rank", "pgid", "tag", "species"):
                if row[name] < 0 or row[name] > SIGNED_INT_MAX:
                    raise ValueError(
                        f"{path} relativistic_hc sample has out-of-range {name}")
            if row["rank"] != rank:
                raise ValueError(
                    f"{path} relativistic_hc sample row rank {row['rank']} "
                    f"does not match shard rank {rank}")
            if not sample_selected(
                    row["species"], row["tag"], selected_species,
                    sample_stride, sample_offset):
                raise ValueError(
                    f"{path} relativistic_hc sample row does not satisfy "
                    "selection metadata")


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
                _validate_psamp_block(path, current)
                blocks.append(current)
            current = {
                "path": path,
                "rank": _rank_from_path(path),
                "nranks": _psamp_header_nranks(path, stripped),
                "metadata": {},
                "columns": [],
                "rows": [],
            }
            columns = None
            continue
        if stripped.startswith("# metadata:"):
            if current is None:
                raise ValueError(f"{path}:{line_number} metadata before header")
            if columns is not None:
                raise ValueError(f"{path}:{line_number} metadata after columns")
            parsed = _parse_metadata_line(stripped.split(":", 1)[1])
            duplicates = current["metadata"].keys() & parsed.keys()
            if duplicates:
                raise ValueError(
                    f"{path}:{line_number} repeats metadata keys "
                    f"{sorted(duplicates)}")
            current["metadata"].update(parsed)
            continue
        if stripped.startswith("# columns:"):
            if current is None:
                raise ValueError(f"{path}:{line_number} columns before header")
            if columns is not None:
                raise ValueError(f"{path}:{line_number} repeats sample columns")
            columns = _canonical_psamp_columns(
                path, line_number, stripped.split(":", 1)[1].split(),
                current["metadata"])
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
        _validate_psamp_block(path, current)
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

    typed_v2_checkpoints: Dict[Path, Dict] = {}
    restarts = [
        read_prst_file(file_path, _typed_v2_checkpoints=typed_v2_checkpoints)
        for file_path in files
    ]
    total = sum(restart["count"] for restart in restarts)
    rank_counts = Counter()
    species_counts = Counter()
    pgid_counts = Counter()
    real_sizes = Counter()
    format_versions = Counter()
    tags = set()
    duplicates = []

    for restart in restarts:
        checkpoint_number = _file_number(restart["path"])
        rank_counts[restart["rank"]] += restart["count"]
        real_sizes[restart["real_size"]] += 1
        format_versions[restart["format_version"]] += 1
        for particle in restart["particles"]:
            key = (checkpoint_number, particle["species"], particle["tag"])
            if key in tags:
                duplicates.append(key)
            tags.add(key)
            species_counts[particle["species"]] += 1
            pgid_counts[particle["pgid"]] += 1

    if duplicates:
        raise ValueError(
            f"Duplicate particle tags within checkpoint and species: {duplicates[:5]}")

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

    df_record = read_df_record(df_file, nspecies, nbin)
    dxh_record = read_dxh_record(dxh_file, nspecies, nbin)
    df = df_record["histogram"]
    dxh = dxh_record["histogram"]
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
            optional[name] = {"file": hist_file, **read_scalar_histogram_record(
                hist_file, marker, nspecies, nbin)}
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
    result = {
        "df": df,
        "dxh": dxh,
        "df_file": df_file,
        "dxh_file": dxh_file,
        "df_record": df_record,
        "dxh_record": dxh_record,
    }
    result.update(optional)
    pspec_path = path / "pspec"
    if pspec_path.exists():
        result["pspec"] = validate_spectra(path, expected_species, nbin)
    return result


def validate_spectra(path: Path, expected_species: Sequence[int],
                     nbin: int) -> Dict:
    nspecies = len(expected_species)
    pspec_path = path / "pspec"
    if pspec_path.is_dir():
        pspec_file = max(pspec_path.glob("*.pspec"), key=_file_number)
    else:
        pspec_file = pspec_path
    record = read_pspec_record(pspec_file, nspecies, nbin)
    spectrum = record["histogram"]
    for species, expected in enumerate(expected_species):
        total = sum(spectrum[species])
        if total != expected:
            raise ValueError(
                f"PSPEC species {species} sums to {total}, expected {expected}")
    return {"file": pspec_file, **record}


def validate_joint_spectra(path: Path, expected_species: Sequence[int],
                           nbin1: int, nbin2: int) -> Dict:
    nspecies = len(expected_species)
    pspec2_path = path / "pspec2"
    if pspec2_path.is_dir():
        pspec2_file = max(pspec2_path.glob("*.pspec2"), key=_file_number)
    else:
        pspec2_file = pspec2_path
    record = read_pspec2_record(pspec2_file, nspecies, nbin1, nbin2)
    spectrum = record["histogram"]
    for species, expected in enumerate(expected_species):
        total = sum(sum(row) for row in spectrum[species])
        if total != expected:
            raise ValueError(
                f"PSPEC2 species {species} sums to {total}, expected {expected}")
    return {"file": pspec2_file, **record}


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
        for name in ("drh", "dparh"):
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
