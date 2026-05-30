#!/usr/bin/env python3
"""Validate the registered one-rank Frontier GPU paper-mode coupling oracle."""

from __future__ import annotations

import argparse
import io
import json
import math
from pathlib import Path
import re
import struct
import sys


LIGHT_SPEED = 3.0
EXPECTED_PARTICLE_COUNT = 64
EXPECTED_MHD_CELL_COUNT = 64
MOMENTUM_TOLERANCE = 2.0e-6
ENERGY_TOLERANCE = 3.0e-5
PARTICLE_MOMENTUM_LIVENESS_MINIMUM = 1.0e-3
COEFFICIENT_INVARIANCE_TOLERANCE = 1.0e-12
RUNTIME_TOKENS = (
    "physical_mode=paper_mhd_pic",
    "state=momentum_p_over_m",
    "C=3 ",
    "background=coupled",
    "feedback=coupled",
    "induction=ideal_mhd_only",
)
MAPPING_PATTERN = re.compile(
    r"^rank=(?P<rank>\d+) host=(?P<host>\S+) "
    r"ROCR_VISIBLE_DEVICES=(?P<rocr>\S+) "
    r"GPU_DEVICE_ORDINAL=(?P<ordinal>\S+) "
    r"HIP_VISIBLE_DEVICES=(?P<hip>\S+)$"
)


def require_file(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"Missing F1 coupling artifact: {path}")
    return path.read_text(encoding="utf-8")


def require_nonempty(path: Path) -> str:
    text = require_file(path)
    if not text.strip():
        raise ValueError(f"Empty F1 coupling artifact: {path}")
    return text


def require_linked_library(ldd_text: str, library: str) -> None:
    pattern = re.compile(
        rf"^\s*{re.escape(library)}[^ ]*\s+=>\s+(?!not found\b)\S+",
        re.MULTILINE,
    )
    if pattern.search(ldd_text) is None:
        raise ValueError(f"Athena executable is not linked against {library}")


def require_finite(values: tuple[float, ...] | list[float], label: str) -> None:
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"Non-finite {label}")


def read_exact(stream: io.BytesIO, size: int, *, path: Path, label: str) -> bytes:
    payload = stream.read(size)
    if len(payload) != size:
        raise ValueError(f"Truncated {label} in Athena binary output: {path}")
    return payload


def parse_assignment(line: bytes, *, path: Path, label: str) -> str:
    if b"=" not in line:
        raise ValueError(f"Malformed {label} in Athena binary output: {path}")
    return line.split(b"=", 1)[1].strip().decode("ascii")


def parse_athena_binary(path: Path) -> dict[str, object]:
    """Parse the bounded native-endian Athena binary output version 1.1 subset."""
    contents = path.read_bytes()
    stream = io.BytesIO(contents)
    if stream.readline().strip() != b"Athena binary output version=1.1":
        raise ValueError(f"Unsupported Athena binary output header: {path}")

    preheader_lines = int(parse_assignment(
        stream.readline(), path=path, label="preheader size"
    ))
    metadata = {}
    for _ in range(preheader_lines - 1):
        line = stream.readline()
        if b"=" not in line:
            raise ValueError(f"Malformed Athena binary preheader: {path}")
        key, value = line.decode("ascii").split("=", 1)
        metadata[key.strip()] = value.strip()

    nvars = int(parse_assignment(
        stream.readline(), path=path, label="number of variables"
    ))
    variables_line = stream.readline()
    if b":" not in variables_line:
        raise ValueError(f"Malformed variable list in Athena binary output: {path}")
    variables = [
        token.decode("ascii") for token in variables_line.split(b":", 1)[1].split()
    ]
    if nvars != len(variables) or not variables:
        raise ValueError(f"Inconsistent variable count in Athena binary output: {path}")

    ascii_header_size = int(parse_assignment(
        stream.readline(), path=path, label="simulation header size"
    ))
    read_exact(stream, ascii_header_size, path=path, label="simulation header")

    location_size = int(metadata["size of location"])
    variable_size = int(metadata["size of variable"])
    if location_size not in (4, 8) or variable_size not in (4, 8):
        raise ValueError(f"Unsupported Athena binary float size: {path}")
    location_code = "f" if location_size == 4 else "d"
    variable_code = "f" if variable_size == 4 else "d"

    integrals = {name: 0.0 for name in variables}
    values = {name: [] for name in variables}
    geometries = []
    blocks = 0
    while stream.tell() < len(contents):
        bounds = struct.unpack(
            "=6i", read_exact(stream, 24, path=path, label="meshblock bounds")
        )
        read_exact(stream, 16, path=path, label="meshblock logical location")
        geometry = struct.unpack(
            "=" + 6 * location_code,
            read_exact(stream, 6 * location_size, path=path, label="meshblock geometry"),
        )
        nx1 = bounds[1] - bounds[0] + 1
        nx2 = bounds[3] - bounds[2] + 1
        nx3 = bounds[5] - bounds[4] + 1
        cells = nx1 * nx2 * nx3
        if cells <= 0:
            raise ValueError(f"Invalid meshblock bounds in Athena binary output: {path}")
        volume = (
            (geometry[1] - geometry[0])
            * (geometry[3] - geometry[2])
            * (geometry[5] - geometry[4])
        )
        require_finite(list(geometry), f"meshblock geometry in {path}")
        if not math.isfinite(volume) or volume <= 0.0:
            raise ValueError(f"Invalid meshblock volume in Athena binary output: {path}")
        geometries.append(geometry)
        for name in variables:
            payload = read_exact(
                stream, cells * variable_size, path=path, label=f"{name} payload"
            )
            block_values = struct.unpack("=" + str(cells) + variable_code, payload)
            require_finite(list(block_values), f"{name} payload in {path}")
            integrals[name] += math.fsum(block_values) * volume / cells
            values[name].extend(block_values)
        blocks += 1
    if blocks == 0:
        raise ValueError(f"Athena binary output contains no meshblocks: {path}")
    time = float(metadata["time"])
    if not math.isfinite(time):
        raise ValueError(f"Non-finite Athena binary time: {path}")
    require_finite(list(integrals.values()), f"integrated payload in {path}")
    return {
        "time": time,
        "cycle": int(metadata["cycle"]),
        "variables": variables,
        "integrals": integrals,
        "values": values,
        "geometries": geometries,
        "meshblock_count": blocks,
    }


def find_marker(contents: bytes, pattern: bytes, offset: int, *, path: Path,
                label: str) -> tuple[int, re.Match[bytes]]:
    match = re.search(pattern, contents[offset:])
    if match is None:
        raise ValueError(f"Missing {label} in particle VTK output: {path}")
    return offset + match.end(), match


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

    offset, points = find_marker(
        contents, rb"\nPOINTS\s+([0-9]+)\s+float\n", 0, path=path, label="POINTS marker"
    )
    count = int(points.group(1))
    offset += 12 * count
    for name in ("gid", "ptag", "species"):
        pattern = (
            rb"\nSCALARS " + name.encode("ascii") + rb" int\nLOOKUP_TABLE default\n"
        )
        offset, _ = find_marker(contents, pattern, offset, path=path, label=name)
        offset += 4 * count
    for name in ("deltaf_f0", "deltaf_weight"):
        pattern = (
            rb"\nSCALARS " + name.encode("ascii") + rb" float\nLOOKUP_TABLE default\n"
        )
        offset, _ = find_marker(contents, pattern, offset, path=path, label=name)
        offset += 4 * count
    offset, _ = find_marker(
        contents, rb"\nVECTORS vel float\n", offset, path=path, label="velocity marker"
    )
    payload = contents[offset:offset + 12 * count]
    if len(payload) != 12 * count:
        raise ValueError(f"Truncated velocity payload in particle VTK output: {path}")

    momentum = [0.0, 0.0, 0.0]
    energy = 0.0
    for index in range(count):
        velocity = struct.unpack(">fff", payload[12 * index:12 * (index + 1)])
        require_finite(list(velocity), f"particle velocity in {path}")
        velocity_squared = math.fsum(component * component for component in velocity)
        if velocity_squared >= LIGHT_SPEED * LIGHT_SPEED:
            raise ValueError(f"Particle speed reaches artificial light speed: {path}")
        gamma = 1.0 / math.sqrt(1.0 - velocity_squared / (LIGHT_SPEED * LIGHT_SPEED))
        for component in range(3):
            momentum[component] += gamma * velocity[component]
        energy += (gamma - 1.0) * LIGHT_SPEED * LIGHT_SPEED
    return {
        "time": time,
        "cycle": cycle,
        "count": count,
        "momentum": tuple(momentum),
        "energy": energy,
    }


def subtract(left: tuple[float, ...], right: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(left[index] - right[index] for index in range(len(left)))


def max_abs(values: tuple[float, ...] | list[float]) -> float:
    require_finite(values, "comparison values")
    return max(abs(value) for value in values)


def require_runtime_identity(stdout: str, label: str) -> str:
    traces = [line for line in stdout.splitlines() if "PIC runtime model:" in line]
    if len(traces) != 1:
        raise ValueError(f"Expected one PIC runtime identity trace for {label}")
    trace = traces[0]
    for token in RUNTIME_TOKENS:
        if token not in trace:
            raise ValueError(f"Missing runtime identity token for {label}: {token}")
    return trace


def require_binary(path: Path, quantity: str, expected_cycle: int) -> dict[str, object]:
    snapshot = parse_athena_binary(path)
    if snapshot["cycle"] != expected_cycle:
        raise ValueError(f"Expected cycle {expected_cycle} Athena binary output: {path}")
    if quantity not in snapshot["variables"]:
        raise ValueError(f"Missing {quantity} in Athena binary output: {path}")
    if snapshot["meshblock_count"] != 1:
        raise ValueError(f"Expected one meshblock in Athena binary output: {path}")
    if len(snapshot["values"][quantity]) != EXPECTED_MHD_CELL_COUNT:
        raise ValueError(f"Expected {EXPECTED_MHD_CELL_COUNT} MHD cells: {path}")
    return snapshot


def parse_case(artifact_dir: Path, label: str) -> dict[str, object]:
    basename = "f1_gpu_paper_coupling_" + label
    output_dir = artifact_dir / "output" / label
    stdout = require_nonempty(artifact_dir / f"athena_{label}_stdout.txt")
    stderr = require_file(artifact_dir / f"athena_{label}_stderr.txt")
    if stderr:
        raise ValueError(f"Athena stderr is not empty for {label}")
    runtime_identity = require_runtime_identity(stdout, label)

    particles_initial = parse_particle_vtk(
        output_dir / "pvtk" / f"{basename}.prtcl_all.00000.part.vtk"
    )
    particles_final = parse_particle_vtk(
        output_dir / "pvtk" / f"{basename}.prtcl_all.00002.part.vtk"
    )
    if particles_initial["cycle"] != 0 or particles_final["cycle"] != 2:
        raise ValueError(f"Unexpected particle VTK cycles for {label}")
    if (
        particles_initial["count"] != EXPECTED_PARTICLE_COUNT
        or particles_final["count"] != EXPECTED_PARTICLE_COUNT
    ):
        raise ValueError(f"Expected {EXPECTED_PARTICLE_COUNT} particles for {label}")

    mhd_delta = []
    for file_id, quantity in (
        ("mhd_u_m1", "mom1"),
        ("mhd_u_m2", "mom2"),
        ("mhd_u_m3", "mom3"),
    ):
        initial = require_binary(
            output_dir / "bin" / f"{basename}.{file_id}.00000.bin", quantity, 0
        )
        final = require_binary(
            output_dir / "bin" / f"{basename}.{file_id}.00002.bin", quantity, 2
        )
        if initial["geometries"] != final["geometries"]:
            raise ValueError(f"Initial and final {quantity} geometries differ for {label}")
        mhd_delta.append(final["integrals"][quantity] - initial["integrals"][quantity])
    energy_initial = require_binary(
        output_dir / "bin" / f"{basename}.mhd_u_e.00000.bin", "ener", 0
    )
    energy_final = require_binary(
        output_dir / "bin" / f"{basename}.mhd_u_e.00002.bin", "ener", 2
    )
    if energy_initial["geometries"] != energy_final["geometries"]:
        raise ValueError(f"Initial and final ener geometries differ for {label}")
    bcc_final = require_binary(
        output_dir / "bin" / f"{basename}.mhd_bcc.00002.bin", "bcc1", 2
    )
    for quantity in ("bcc1", "bcc2", "bcc3"):
        if quantity not in bcc_final["variables"]:
            raise ValueError(f"Missing {quantity} in final MHD magnetic output for {label}")

    particle_momentum_delta = subtract(
        particles_final["momentum"], particles_initial["momentum"]
    )
    total_momentum_delta = tuple(
        particle_momentum_delta[index] + mhd_delta[index] for index in range(3)
    )
    total_energy_delta = (
        particles_final["energy"]
        - particles_initial["energy"]
        + energy_final["integrals"]["ener"]
        - energy_initial["integrals"]["ener"]
    )
    momentum_error = max_abs(total_momentum_delta)
    energy_error = abs(total_energy_delta)
    momentum_liveness = max_abs(particle_momentum_delta)
    if momentum_error > MOMENTUM_TOLERANCE:
        raise ValueError(f"Momentum conservation error exceeds tolerance for {label}")
    if energy_error > ENERGY_TOLERANCE:
        raise ValueError(f"Energy conservation error exceeds tolerance for {label}")
    if momentum_liveness <= PARTICLE_MOMENTUM_LIVENESS_MINIMUM:
        raise ValueError(f"Particle momentum coupling is inactive for {label}")

    return {
        "runtime_identity": runtime_identity,
        "particle_count": particles_final["count"],
        "particle_momentum_initial": particles_initial["momentum"],
        "particle_momentum_final": particles_final["momentum"],
        "particle_energy_initial": particles_initial["energy"],
        "particle_energy_final": particles_final["energy"],
        "mhd_momentum_delta": mhd_delta,
        "mhd_energy_delta": (
            energy_final["integrals"]["ener"] - energy_initial["integrals"]["ener"]
        ),
        "total_momentum_delta": total_momentum_delta,
        "total_energy_delta": total_energy_delta,
        "max_abs_momentum_conservation_error": momentum_error,
        "abs_energy_conservation_error": energy_error,
        "max_abs_particle_momentum_delta": momentum_liveness,
        "bcc_final": bcc_final["values"],
        "bcc_geometries": bcc_final["geometries"],
    }


def require_execution_context(artifact_dir: Path) -> dict[str, object]:
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
        line.split("=", 1) for line in environment_text.splitlines() if "=" in line
    )
    if environment.get("MPICH_GPU_SUPPORT_ENABLED") != "1":
        raise ValueError("MPICH_GPU_SUPPORT_ENABLED must be 1")
    if environment.get("SLURM_EXPORT_ENV") != "ALL":
        raise ValueError("SLURM_EXPORT_ENV must be ALL for compute steps")
    return {
        "gpu_mapping_rank": int(mapping["rank"]),
        "rocr_visible_devices": mapping["rocr"],
        "hip_linked": True,
        "gpu_aware_mpi_linked": True,
        "mpich_gpu_support_enabled": environment["MPICH_GPU_SUPPORT_ENABLED"],
        "slurm_export_env": environment["SLURM_EXPORT_ENV"],
        "snapshot_verification": "pass",
    }


def analyze(artifact_dir: Path) -> dict[str, object]:
    execution_context = require_execution_context(artifact_dir)
    cases = {label: parse_case(artifact_dir, label) for label in ("coeff0", "coeff7")}
    coeff0 = cases["coeff0"]
    coeff7 = cases["coeff7"]

    bcc_errors = {}
    if coeff0["bcc_geometries"] != coeff7["bcc_geometries"]:
        raise ValueError("Coefficient comparison magnetic geometries differ")
    for quantity in ("bcc1", "bcc2", "bcc3"):
        left = coeff0["bcc_final"][quantity]
        right = coeff7["bcc_final"][quantity]
        if len(left) != len(right):
            raise ValueError(f"Coefficient comparison {quantity} payload sizes differ")
        error = max_abs([left[index] - right[index] for index in range(len(left))])
        if error > COEFFICIENT_INVARIANCE_TOLERANCE:
            raise ValueError(f"Coefficient-invariant induction check failed for {quantity}")
        bcc_errors[quantity] = error
    momentum_difference = subtract(
        coeff0["particle_momentum_final"], coeff7["particle_momentum_final"]
    )
    momentum_invariance_error = max_abs(momentum_difference)
    if momentum_invariance_error > COEFFICIENT_INVARIANCE_TOLERANCE:
        raise ValueError("Coefficient-invariant particle momentum check failed")

    for case in cases.values():
        del case["bcc_final"]
        del case["bcc_geometries"]
    result = {
        "schema_version": 1,
        "status": "pass",
        "physical_mode": "paper_mhd_pic",
        "execution_context": execution_context,
        "cases": cases,
        "coefficient_invariance": {
            "max_abs_bcc_error": bcc_errors,
            "max_abs_particle_momentum_error": momentum_invariance_error,
        },
        "thresholds": {
            "expected_particle_count": EXPECTED_PARTICLE_COUNT,
            "expected_mhd_cell_count": EXPECTED_MHD_CELL_COUNT,
            "momentum_conservation_tolerance": MOMENTUM_TOLERANCE,
            "energy_conservation_tolerance": ENERGY_TOLERANCE,
            "particle_momentum_liveness_minimum": PARTICLE_MOMENTUM_LIVENESS_MINIMUM,
            "coefficient_invariance_tolerance": COEFFICIENT_INVARIANCE_TOLERANCE,
        },
        "athena_binary_parser": {
            "format_version": "1.1",
            "payload_byteorder": sys.byteorder,
            "dependencies": "python_stdlib_only",
        },
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
