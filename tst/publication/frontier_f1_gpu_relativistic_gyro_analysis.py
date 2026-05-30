#!/usr/bin/env python3
"""Validate the registered one-rank Frontier GPU relativistic gyro oracle."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import struct


LIGHT_SPEED = 3.0
INITIAL_MOMENTUM = (1.0, 0.0, 0.0)
MAGNETIC_FIELD = (0.0, 0.0, 1.0)
VELOCITY_TOLERANCE = 2.0e-6
MAPPING_PATTERN = re.compile(
    r"^rank=(?P<rank>\d+) host=(?P<host>\S+) "
    r"ROCR_VISIBLE_DEVICES=(?P<rocr>\S+) "
    r"GPU_DEVICE_ORDINAL=(?P<ordinal>\S+) "
    r"HIP_VISIBLE_DEVICES=(?P<hip>\S+)$"
)


def require_file(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"Missing F1 artifact: {path}")
    return path.read_text(encoding="utf-8")


def require_nonempty(path: Path) -> str:
    text = require_file(path)
    if not text.strip():
        raise ValueError(f"Empty F1 artifact: {path}")
    return text


def require_linked_library(ldd_text: str, library: str) -> None:
    pattern = re.compile(
        rf"^\s*{re.escape(library)}[^ ]*\s+=>\s+(?!not found\b)\S+",
        re.MULTILINE,
    )
    if pattern.search(ldd_text) is None:
        raise ValueError(f"Athena executable is not linked against {library}")


def cross(left: tuple[float, float, float],
          right: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def add(left: tuple[float, float, float],
        right: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(left[index] + right[index] for index in range(3))


def scale(value: float,
          vector: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(value * component for component in vector)


def norm2(vector: tuple[float, float, float]) -> float:
    return sum(component * component for component in vector)


def expected_velocity(cycle: int, time: float) -> tuple[float, float, float]:
    state = INITIAL_MOMENTUM
    dt = time / cycle
    for _ in range(cycle):
        gamma = math.sqrt(1.0 + norm2(state) / (LIGHT_SPEED * LIGHT_SPEED))
        rotation = scale(0.5 * dt / gamma, MAGNETIC_FIELD)
        correction = scale(2.0 / (1.0 + norm2(rotation)), rotation)
        state_prime = add(state, cross(state, rotation))
        state = add(state, cross(state_prime, correction))
    gamma = math.sqrt(1.0 + norm2(state) / (LIGHT_SPEED * LIGHT_SPEED))
    return scale(1.0 / gamma, state)


def parse_particle_vtk(path: Path) -> dict[str, object]:
    contents = path.read_bytes()
    header = re.search(
        rb"# AthenaK particle data at time=\s*([^ ]+)\s+nranks=.*cycle=([0-9]+)",
        contents,
    )
    if header is None:
        raise ValueError(f"Missing particle VTK header: {path}")
    time = float(header.group(1))
    if not math.isfinite(time):
        raise ValueError(f"Non-finite particle VTK time: {path}")
    cycle = int(header.group(2))

    points = re.search(rb"\nPOINTS\s+([0-9]+)\s+float\n", contents)
    if points is None:
        raise ValueError(f"Missing POINTS marker: {path}")
    count = int(points.group(1))
    offset = points.end() + 12 * count
    for name in ("gid", "ptag", "species"):
        marker = re.search(
            rb"\nSCALARS " + name.encode("ascii") + rb" int\nLOOKUP_TABLE default\n",
            contents[offset:],
        )
        if marker is None:
            raise ValueError(f"Missing {name} marker: {path}")
        offset += marker.end() + 4 * count
    for name in ("deltaf_f0", "deltaf_weight"):
        marker = re.search(
            rb"\nSCALARS " + name.encode("ascii") + rb" float\nLOOKUP_TABLE default\n",
            contents[offset:],
        )
        if marker is None:
            raise ValueError(f"Missing {name} marker: {path}")
        offset += marker.end() + 4 * count
    marker = re.search(rb"\nVECTORS vel float\n", contents[offset:])
    if marker is None:
        raise ValueError(f"Missing velocity marker: {path}")
    offset += marker.end()
    payload = contents[offset:offset + 12 * count]
    if len(payload) != 12 * count:
        raise ValueError(f"Truncated velocity payload: {path}")
    velocities = [
        struct.unpack(">fff", payload[12 * index:12 * (index + 1)])
        for index in range(count)
    ]
    if any(
        not math.isfinite(component)
        for velocity in velocities
        for component in velocity
    ):
        raise ValueError(f"Non-finite particle velocity: {path}")
    return {"time": time, "cycle": cycle, "count": count, "velocities": velocities}


def analyze(artifact_dir: Path) -> dict[str, object]:
    snapshot_verification = require_nonempty(
        artifact_dir / "snapshot_verification.txt"
    )
    if snapshot_verification.strip() != "immutable snapshot verified":
        raise ValueError("Unexpected snapshot verification artifact")
    environment_text = require_nonempty(artifact_dir / "environment.allowlist.txt")
    require_nonempty(artifact_dir / "modules.txt")
    rocm_smi_text = require_nonempty(artifact_dir / "rocm_smi.txt")
    if (
        "ROCm System Management Interface" not in rocm_smi_text
        or "AMD INSTINCT MI200" not in rocm_smi_text
    ):
        raise ValueError("Missing expected Frontier ROCm device identity")
    ldd_text = require_nonempty(artifact_dir / "athena_ldd.txt")
    for library in ("libamdhip64", "libmpi_amd", "libmpi_gtl_hsa"):
        require_linked_library(ldd_text, library)

    mapping_text = require_nonempty(artifact_dir / "gpu_mapping.txt").strip()
    mapping_match = MAPPING_PATTERN.fullmatch(mapping_text)
    if mapping_match is None:
        raise ValueError(f"Malformed GPU mapping line: {mapping_text}")
    mapping = mapping_match.groupdict()
    if mapping["rank"] != "0":
        raise ValueError(f"Expected GPU mapping rank 0, found {mapping['rank']}")
    if re.fullmatch(r"\d+", mapping["rocr"]) is None:
        raise ValueError("ROCR_VISIBLE_DEVICES must be one numeric GPU identifier")

    environment = dict(
        line.split("=", 1)
        for line in environment_text.splitlines()
        if "=" in line
    )
    if environment.get("MPICH_GPU_SUPPORT_ENABLED") != "1":
        raise ValueError("MPICH_GPU_SUPPORT_ENABLED must be 1")
    if environment.get("SLURM_EXPORT_ENV") != "ALL":
        raise ValueError("SLURM_EXPORT_ENV must be ALL for compute steps")

    stdout = require_nonempty(artifact_dir / "athena_stdout.txt")
    stderr = require_file(artifact_dir / "athena_stderr.txt")
    if stderr:
        raise ValueError("Athena stderr is not empty")
    if "physical_mode=paper_test_particle" not in stdout:
        raise ValueError("Missing paper_test_particle runtime identity")
    paths = sorted((artifact_dir / "output" / "pvtk").glob(
        "f1_gpu_relativistic_gyro.prtcl_all.*.part.vtk"
    ))
    if not paths:
        raise ValueError("Missing particle VTK output")
    snapshot = parse_particle_vtk(paths[-1])
    if snapshot["cycle"] != 2:
        raise ValueError(f"Expected cycle 2 output, found {snapshot['cycle']}")
    if snapshot["count"] != 64:
        raise ValueError(f"Expected 64 particles, found {snapshot['count']}")
    expected = expected_velocity(int(snapshot["cycle"]), float(snapshot["time"]))
    if any(not math.isfinite(component) for component in expected):
        raise ValueError("Non-finite expected Boris velocity")
    max_error = max(
        abs(measured[index] - expected[index])
        for measured in snapshot["velocities"]
        for index in range(3)
    )
    if not math.isfinite(max_error):
        raise ValueError("Non-finite Boris velocity error")
    if max_error > VELOCITY_TOLERANCE:
        raise ValueError(f"Boris velocity error exceeds tolerance: {max_error}")
    result = {
        "schema_version": 1,
        "status": "pass",
        "physical_mode": "paper_test_particle",
        "execution_context": {
            "gpu_mapping_rank": int(mapping["rank"]),
            "rocr_visible_devices": mapping["rocr"],
            "hip_linked": True,
            "gpu_aware_mpi_linked": True,
            "mpich_gpu_support_enabled": environment["MPICH_GPU_SUPPORT_ENABLED"],
            "slurm_export_env": environment["SLURM_EXPORT_ENV"],
            "snapshot_verification": "pass",
        },
        "cycle": snapshot["cycle"],
        "particle_count": snapshot["count"],
        "expected_velocity": expected,
        "max_abs_velocity_error": max_error,
        "velocity_tolerance": VELOCITY_TOLERANCE,
    }
    (artifact_dir / "analysis.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.artifact_dir)


if __name__ == "__main__":
    main()
