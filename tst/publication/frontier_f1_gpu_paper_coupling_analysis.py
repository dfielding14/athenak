#!/usr/bin/env python3
"""Validate the registered one-rank Frontier GPU paper-mode coupling oracle."""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import stat
import struct
import sys


def _load_artifact_helpers() -> object:
    inherited_fd = os.environ.pop("PIC_F1_ANALYSIS_HELPER_FD", None)
    path = (
        Path("/proc/self/fd") / inherited_fd
        if inherited_fd is not None
        else Path(__file__).with_name("frontier_f1_structured_artifacts.py")
    )
    loader = importlib.machinery.SourceFileLoader(
        "_frontier_f1_structured_artifacts", str(path)
    )
    spec = importlib.util.spec_from_loader("_frontier_f1_structured_artifacts", loader)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load structured F1 artifact helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_ARTIFACT_HELPERS = _load_artifact_helpers()
StructuredArtifactTree = _ARTIFACT_HELPERS.StructuredArtifactTree
load_inventory = _ARTIFACT_HELPERS.load_inventory
offline_analysis_receipt = _ARTIFACT_HELPERS.offline_analysis_receipt
read_inventory_bytes = _ARTIFACT_HELPERS.read_inventory_bytes
require_inventory_sha256 = _ARTIFACT_HELPERS.require_inventory_sha256
validate_frontier_mpich_diagnostic_stderr = (
    _ARTIFACT_HELPERS.validate_frontier_mpich_diagnostic_stderr
)
write_result_exclusive = _ARTIFACT_HELPERS.write_result_exclusive


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
RUNTIME_ALLOWLIST_KEYS = [
    "PIC_FRONTIER_PROFILE",
    "HSA_XNACK",
    "MPICH_ENV_DISPLAY",
    "MPICH_VERSION_DISPLAY",
    "MPICH_GPU_SUPPORT_ENABLED",
    "MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED",
    "MPICH_OFI_NIC_POLICY",
    "MPICH_GPU_IPC_CACHE_MAX_SIZE",
    "MPICH_MPIIO_HINTS",
    "MPICH_OFI_NUM_CQ_ENTRIES",
    "FI_MR_CACHE_MONITOR",
    "FI_CXI_RX_MATCH_MODE",
    "OMP_NUM_THREADS",
    "SLURM_EXPORT_ENV",
    "ROCM_PATH",
    "LOADEDMODULES",
    "_LMFILES_",
    "MODULEPATH",
]
EXPECTED_BASELINE = {
    "PIC_FRONTIER_PROFILE": "frontier_minimum_supported",
    "HSA_XNACK": "0",
    "MPICH_ENV_DISPLAY": "1",
    "MPICH_VERSION_DISPLAY": "1",
    "MPICH_GPU_SUPPORT_ENABLED": "1",
    "MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED": "<unset>",
    "MPICH_OFI_NIC_POLICY": "<unset>",
    "MPICH_GPU_IPC_CACHE_MAX_SIZE": "<unset>",
    "MPICH_MPIIO_HINTS": "<unset>",
    "MPICH_OFI_NUM_CQ_ENTRIES": "<unset>",
    "FI_MR_CACHE_MONITOR": "<unset>",
    "FI_CXI_RX_MATCH_MODE": "<unset>",
    "OMP_NUM_THREADS": "<unset>",
    "SLURM_EXPORT_ENV": "ALL",
    "ROCM_PATH": "/opt/rocm-6.2.4",
}
EXPECTED_PROVENANCE_SHA256 = {
    "LOADEDMODULES": "0fa38c7f44ada1f16e61007e577d996267ebbda53b6434a5bc27f97451019eea",
    "_LMFILES_": "dfd37544d54b83574d67af2dd3835b7bc5100f7424bdd2e0a009bb67b9a5174b",
    "MODULEPATH": "825366ca3c91985b4a4ff7ff84ec507bed7d2e433d460a34ba98ee3708b2c9c4",
}
GPU_LAUNCH_PATTERN = re.compile(
    r"^PIC trusted GPU launch: rank=(?P<rank>[0-9]+) host=(?P<host>\S+) "
    r"ROCR_VISIBLE_DEVICES=(?P<rocr>[0-9]+) "
    r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa$",
    re.MULTILINE,
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
        rf"^\s*{re.escape(library)}[.]so(?:[.][0-9]+)*\s+=>\s+(?!not found\b)\S+",
        re.MULTILINE,
    )
    if pattern.search(ldd_text) is None:
        raise ValueError(f"Athena executable is not linked against {library}")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _parse_structured_allowlist(data: bytes) -> dict[str, str]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("Structured runtime allowlist is not UTF-8") from error
    if "\r" in text or not text.endswith("\n"):
        raise ValueError("Structured runtime allowlist must use terminated LF records")
    values = {}
    for record in text.splitlines():
        key, separator, value = record.partition("=")
        if not separator or not key or not value or key in values:
            raise ValueError("Structured runtime allowlist contains a malformed record")
        values[key] = value
    if list(values) != RUNTIME_ALLOWLIST_KEYS:
        raise ValueError("Structured runtime allowlist key set or order is not authorized")
    for key, expected in EXPECTED_BASELINE.items():
        if values[key] != expected:
            raise ValueError(f"Structured runtime allowlist {key} differs from the baseline")
    for key, expected in EXPECTED_PROVENANCE_SHA256.items():
        if _sha256(values[key].encode("utf-8")) != expected:
            raise ValueError(
                f"Structured runtime allowlist {key} differs from reviewed provenance"
            )
    return values


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


def parse_athena_binary_bytes(contents: bytes, *, label: str) -> dict[str, object]:
    """Parse the bounded native-endian Athena binary output version 1.1 subset."""
    stream = io.BytesIO(contents)
    if stream.readline().strip() != b"Athena binary output version=1.1":
        raise ValueError(f"Unsupported Athena binary output header: {label}")

    preheader_lines = int(parse_assignment(
        stream.readline(), path=Path(label), label="preheader size"
    ))
    metadata = {}
    for _ in range(preheader_lines - 1):
        line = stream.readline()
        if b"=" not in line:
            raise ValueError(f"Malformed Athena binary preheader: {label}")
        key, value = line.decode("ascii").split("=", 1)
        metadata[key.strip()] = value.strip()

    nvars = int(parse_assignment(
        stream.readline(), path=Path(label), label="number of variables"
    ))
    variables_line = stream.readline()
    if b":" not in variables_line:
        raise ValueError(f"Malformed variable list in Athena binary output: {label}")
    variables = [
        token.decode("ascii") for token in variables_line.split(b":", 1)[1].split()
    ]
    if nvars != len(variables) or not variables:
        raise ValueError(f"Inconsistent variable count in Athena binary output: {label}")

    ascii_header_size = int(parse_assignment(
        stream.readline(), path=Path(label), label="simulation header size"
    ))
    read_exact(stream, ascii_header_size, path=Path(label), label="simulation header")

    location_size = int(metadata["size of location"])
    variable_size = int(metadata["size of variable"])
    if location_size not in (4, 8) or variable_size not in (4, 8):
        raise ValueError(f"Unsupported Athena binary float size: {label}")
    location_code = "f" if location_size == 4 else "d"
    variable_code = "f" if variable_size == 4 else "d"

    integrals = {name: 0.0 for name in variables}
    values = {name: [] for name in variables}
    geometries = []
    blocks = 0
    while stream.tell() < len(contents):
        bounds = struct.unpack(
            "=6i", read_exact(stream, 24, path=Path(label), label="meshblock bounds")
        )
        read_exact(stream, 16, path=Path(label), label="meshblock logical location")
        geometry = struct.unpack(
            "=" + 6 * location_code,
            read_exact(stream, 6 * location_size, path=Path(label), label="meshblock geometry"),
        )
        nx1 = bounds[1] - bounds[0] + 1
        nx2 = bounds[3] - bounds[2] + 1
        nx3 = bounds[5] - bounds[4] + 1
        cells = nx1 * nx2 * nx3
        if cells <= 0:
            raise ValueError(f"Invalid meshblock bounds in Athena binary output: {label}")
        volume = (
            (geometry[1] - geometry[0])
            * (geometry[3] - geometry[2])
            * (geometry[5] - geometry[4])
        )
        require_finite(list(geometry), f"meshblock geometry in {label}")
        if not math.isfinite(volume) or volume <= 0.0:
            raise ValueError(f"Invalid meshblock volume in Athena binary output: {label}")
        geometries.append(geometry)
        for name in variables:
            payload = read_exact(
                stream, cells * variable_size, path=Path(label), label=f"{name} payload"
            )
            block_values = struct.unpack("=" + str(cells) + variable_code, payload)
            require_finite(list(block_values), f"{name} payload in {label}")
            integrals[name] += math.fsum(block_values) * volume / cells
            values[name].extend(block_values)
        blocks += 1
    if blocks == 0:
        raise ValueError(f"Athena binary output contains no meshblocks: {label}")
    time = float(metadata["time"])
    if not math.isfinite(time):
        raise ValueError(f"Non-finite Athena binary time: {label}")
    require_finite(list(integrals.values()), f"integrated payload in {label}")
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
    match = re.match(pattern, contents[offset:])
    if match is None:
        raise ValueError(f"Missing {label} in particle VTK output: {path}")
    return offset + match.end(), match


def parse_particle_vtk_bytes(contents: bytes, *, label: str) -> dict[str, object]:
    header = re.search(
        rb"# AthenaK particle data at time=\s*([^ ]+)\s+nranks=.*cycle=([0-9]+)",
        contents,
    )
    if header is None:
        raise ValueError(f"Missing particle VTK header: {label}")
    time = float(header.group(1))
    if not math.isfinite(time):
        raise ValueError(f"Non-finite particle VTK time: {label}")
    cycle = int(header.group(2))

    offset, points = find_marker(
        contents,
        rb"\nPOINTS\s+([0-9]+)\s+float\n",
        header.end(),
        path=Path(label),
        label="POINTS marker",
    )
    count = int(points.group(1))
    offset += 12 * count
    for name in ("gid", "ptag", "species"):
        pattern = (
            rb"\nSCALARS " + name.encode("ascii") + rb" int\nLOOKUP_TABLE default\n"
        )
        offset, _ = find_marker(contents, pattern, offset, path=Path(label), label=name)
        offset += 4 * count
    for name in ("deltaf_f0", "deltaf_weight"):
        pattern = (
            rb"\nSCALARS " + name.encode("ascii") + rb" float\nLOOKUP_TABLE default\n"
        )
        offset, _ = find_marker(contents, pattern, offset, path=Path(label), label=name)
        offset += 4 * count
    offset, _ = find_marker(
        contents, rb"\nVECTORS vel float\n", offset, path=Path(label), label="velocity marker"
    )
    payload = contents[offset:offset + 12 * count]
    if len(payload) != 12 * count:
        raise ValueError(f"Truncated velocity payload in particle VTK output: {label}")
    if offset + len(payload) != len(contents):
        raise ValueError(f"Unexpected particle VTK suffix: {label}")

    momentum = [0.0, 0.0, 0.0]
    energy = 0.0
    for index in range(count):
        velocity = struct.unpack(">fff", payload[12 * index:12 * (index + 1)])
        require_finite(list(velocity), f"particle velocity in {label}")
        velocity_squared = math.fsum(component * component for component in velocity)
        if velocity_squared >= LIGHT_SPEED * LIGHT_SPEED:
            raise ValueError(f"Particle speed reaches artificial light speed: {label}")
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


def require_binary(
    artifact_tree: object,
    inventory: dict[str, dict[str, object]],
    relative: str,
    quantity: str,
    expected_cycle: int,
) -> dict[str, object]:
    snapshot = parse_athena_binary_bytes(
        read_inventory_bytes(artifact_tree, inventory, relative),
        label=relative,
    )
    if snapshot["cycle"] != expected_cycle:
        raise ValueError(f"Expected cycle {expected_cycle} Athena binary output: {relative}")
    if quantity not in snapshot["variables"]:
        raise ValueError(f"Missing {quantity} in Athena binary output: {relative}")
    if snapshot["meshblock_count"] != 1:
        raise ValueError(f"Expected one meshblock in Athena binary output: {relative}")
    if len(snapshot["values"][quantity]) != EXPECTED_MHD_CELL_COUNT:
        raise ValueError(f"Expected {EXPECTED_MHD_CELL_COUNT} MHD cells: {relative}")
    return snapshot


def require_trusted_gpu_launch(stdout: str, label: str) -> dict[str, object]:
    matches = list(GPU_LAUNCH_PATTERN.finditer(stdout))
    if len(matches) != 1 or matches[0].group("rank") != "0":
        raise ValueError(f"Missing exact rank-zero trusted GPU launch preflight for {label}")
    return {
        "rank": 0,
        "host": matches[0].group("host"),
        "rocr_visible_devices": matches[0].group("rocr"),
        "linkage": ["libamdhip64", "libmpi_amd", "libmpi_gtl_hsa"],
    }


def registered_case_output_paths(
    inventory: dict[str, dict[str, object]], label: str
) -> list[str]:
    if label not in {"coeff0", "coeff7"}:
        raise ValueError(f"Unregistered coupling case: {label}")
    basename = "f1_gpu_paper_coupling_" + label
    output_dir = f"output/{label}"
    expected = [
        f"{output_dir}/pvtk/{basename}.prtcl_all.{cycle:05d}.part.vtk"
        for cycle in range(3)
    ]
    expected.extend(
        f"{output_dir}/bin/{basename}.{file_id}.{cycle:05d}.bin"
        for file_id in ("mhd_bcc", "mhd_u_e", "mhd_u_m1", "mhd_u_m2", "mhd_u_m3")
        for cycle in range(3)
    )
    registered_paths = []
    for registered_label in ("coeff0", "coeff7"):
        registered_basename = "f1_gpu_paper_coupling_" + registered_label
        registered_dir = f"output/{registered_label}"
        registered_paths.extend(
            f"{registered_dir}/pvtk/{registered_basename}.prtcl_all.{cycle:05d}.part.vtk"
            for cycle in range(3)
        )
        registered_paths.extend(
            f"{registered_dir}/bin/{registered_basename}.{file_id}.{cycle:05d}.bin"
            for file_id in ("mhd_bcc", "mhd_u_e", "mhd_u_m1", "mhd_u_m2", "mhd_u_m3")
            for cycle in range(3)
        )
    paths = sorted(path for path in inventory if path.startswith("output/"))
    if paths != sorted(registered_paths):
        raise ValueError("Coupling output set differs from registered artifacts")
    return expected


def parse_case(
    artifact_tree: object,
    inventory: dict[str, dict[str, object]],
    label: str,
) -> dict[str, object]:
    basename = "f1_gpu_paper_coupling_" + label
    output_dir = f"output/{label}"
    registered_case_output_paths(inventory, label)
    stdout = read_inventory_bytes(
        artifact_tree, inventory, f"athena_{label}_stdout.txt"
    ).decode("utf-8")
    stderr = read_inventory_bytes(
        artifact_tree, inventory, f"athena_{label}_stderr.txt"
    ).decode("utf-8")
    if not stdout.strip():
        raise ValueError(f"Athena stdout is empty for {label}")
    validate_frontier_mpich_diagnostic_stderr(stderr)
    runtime_identity = require_runtime_identity(stdout, label)
    gpu_launch = require_trusted_gpu_launch(stdout, label)

    particles_initial_relative = f"{output_dir}/pvtk/{basename}.prtcl_all.00000.part.vtk"
    particles_final_relative = f"{output_dir}/pvtk/{basename}.prtcl_all.00002.part.vtk"
    particles_initial = parse_particle_vtk_bytes(
        read_inventory_bytes(artifact_tree, inventory, particles_initial_relative),
        label=particles_initial_relative,
    )
    particles_final = parse_particle_vtk_bytes(
        read_inventory_bytes(artifact_tree, inventory, particles_final_relative),
        label=particles_final_relative,
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
            artifact_tree,
            inventory,
            f"{output_dir}/bin/{basename}.{file_id}.00000.bin",
            quantity,
            0,
        )
        final = require_binary(
            artifact_tree,
            inventory,
            f"{output_dir}/bin/{basename}.{file_id}.00002.bin",
            quantity,
            2,
        )
        if initial["geometries"] != final["geometries"]:
            raise ValueError(f"Initial and final {quantity} geometries differ for {label}")
        mhd_delta.append(final["integrals"][quantity] - initial["integrals"][quantity])
    energy_initial = require_binary(
        artifact_tree,
        inventory,
        f"{output_dir}/bin/{basename}.mhd_u_e.00000.bin",
        "ener",
        0,
    )
    energy_final = require_binary(
        artifact_tree,
        inventory,
        f"{output_dir}/bin/{basename}.mhd_u_e.00002.bin",
        "ener",
        2,
    )
    if energy_initial["geometries"] != energy_final["geometries"]:
        raise ValueError(f"Initial and final ener geometries differ for {label}")
    bcc_final = require_binary(
        artifact_tree,
        inventory,
        f"{output_dir}/bin/{basename}.mhd_bcc.00002.bin",
        "bcc1",
        2,
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
        "gpu_launch": gpu_launch,
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


def require_execution_context(
    artifact_tree: object,
    inventory: dict[str, dict[str, object]],
) -> dict[str, object]:
    structured_paths = [
        "f1-coupling-coeff0.environment.allowlist.txt",
        "f1-coupling-coeff7.environment.allowlist.txt",
    ]
    if not all(path in inventory for path in structured_paths):
        raise ValueError("Structured coupling runtime allowlist set is incomplete")
    payloads = [
        read_inventory_bytes(artifact_tree, inventory, path) for path in structured_paths
    ]
    parsed = [_parse_structured_allowlist(data) for data in payloads]
    if parsed[0] != parsed[1]:
        raise ValueError("Structured coupling runtime allowlists differ")
    values = parsed[0]
    return {
        "launch_contract": "trusted_trampoline_athena_argv_v1",
        "runtime_profile": values["PIC_FRONTIER_PROFILE"],
        "runtime_allowlist_sha256": {
            path: _sha256(data)
            for path, data in zip(structured_paths, payloads)
        },
        "mpich_gpu_support_enabled": values["MPICH_GPU_SUPPORT_ENABLED"],
        "slurm_export_env": values["SLURM_EXPORT_ENV"],
    }


def _analyze_tree(artifact_tree: object) -> dict[str, object]:
    inventory = load_inventory(artifact_tree)
    execution_context = require_execution_context(artifact_tree, inventory)
    cases = {
        label: parse_case(artifact_tree, inventory, label)
        for label in ("coeff0", "coeff7")
    }
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
    return result


def analyze(artifact_dir: Path) -> dict[str, object]:
    with StructuredArtifactTree(artifact_dir) as artifact_tree:
        return _analyze_tree(artifact_tree)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    parser.add_argument("--artifact-dir-fd", type=int)
    parser.add_argument("--verify-artifact-inventory-sha256")
    parser.add_argument("--verify-result-sha256")
    args = parser.parse_args()
    verification_values = (
        args.artifact_dir_fd,
        args.verify_artifact_inventory_sha256,
        args.verify_result_sha256,
    )
    if any(value is not None for value in verification_values) and not all(
        value is not None for value in verification_values
    ):
        raise ValueError("Offline verification requires the complete parent binding")
    if args.artifact_dir_fd is not None and args.artifact_dir_fd < 0:
        raise ValueError("Offline verification artifact descriptor is malformed")
    if (
        sys.executable != _ARTIFACT_HELPERS.TRUSTED_PYTHON
        or not sys.flags.isolated
        or not sys.dont_write_bytecode
    ):
        raise ValueError("Offline analysis requires the trusted Python -I -B runner")
    analyzer_sha256 = ""
    support_module_sha256 = ""
    if args.verify_result_sha256 is None:
        analyzer_sha256 = _ARTIFACT_HELPERS.read_only_file_sha256(Path(__file__))
        support_module_sha256 = _ARTIFACT_HELPERS.read_only_file_sha256(
            Path(_ARTIFACT_HELPERS.__file__)
        )
    with StructuredArtifactTree(
        args.artifact_dir, inherited_root_fd=args.artifact_dir_fd
    ) as artifact_tree:
        if args.verify_artifact_inventory_sha256 is not None:
            require_inventory_sha256(
                artifact_tree, args.verify_artifact_inventory_sha256
            )
        result = _analyze_tree(artifact_tree)
        if args.verify_result_sha256 is not None:
            require_inventory_sha256(
                artifact_tree, args.verify_artifact_inventory_sha256
            )
            expected = args.verify_result_sha256
            if (
                re.fullmatch(r"[0-9a-f]{64}", expected) is None
                or hashlib.sha256(
                    _ARTIFACT_HELPERS.canonical_json_bytes(result)
                ).hexdigest()
                != expected
            ):
                raise ValueError("Offline analysis recomputation differs from bound result")
            return
        write_result_exclusive(
            artifact_tree, "analysis/analysis.json", result
        )
        write_result_exclusive(
            artifact_tree,
            "analysis/offline_analysis_receipt.json",
            offline_analysis_receipt(
                artifact_tree,
                analyzer_path=Path(__file__),
                analyzer_sha256=analyzer_sha256,
                support_module_sha256=support_module_sha256,
                result=result,
            ),
        )


if __name__ == "__main__":
    main()
