#!/usr/bin/env python3
"""Validate the registered one-rank Frontier GPU relativistic gyro oracle."""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import importlib.util
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
write_result_exclusive = _ARTIFACT_HELPERS.write_result_exclusive


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
        raise ValueError(f"Missing F1 artifact: {path}")
    return path.read_text(encoding="utf-8")


def require_nonempty(path: Path) -> str:
    text = require_file(path)
    if not text.strip():
        raise ValueError(f"Empty F1 artifact: {path}")
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


def _structured_execution_context(
    artifact_tree: object,
    inventory: dict[str, dict[str, object]],
) -> dict[str, object]:
    data = read_inventory_bytes(
        artifact_tree, inventory, "f1-gyro.environment.allowlist.txt"
    )
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
    return {
        "launch_contract": "trusted_trampoline_athena_argv_v1",
        "runtime_profile": values["PIC_FRONTIER_PROFILE"],
        "runtime_allowlist_sha256": _sha256(data),
        "mpich_gpu_support_enabled": values["MPICH_GPU_SUPPORT_ENABLED"],
        "slurm_export_env": values["SLURM_EXPORT_ENV"],
    }


def _legacy_execution_context(artifact_dir: Path) -> dict[str, object]:
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
    return {
        "gpu_mapping_rank": int(mapping["rank"]),
        "rocr_visible_devices": mapping["rocr"],
        "hip_linked": True,
        "gpu_aware_mpi_linked": True,
        "mpich_gpu_support_enabled": environment["MPICH_GPU_SUPPORT_ENABLED"],
        "slurm_export_env": environment["SLURM_EXPORT_ENV"],
        "snapshot_verification": "pass",
    }


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

    points = re.search(rb"\nPOINTS\s+([0-9]+)\s+float\n", contents)
    if points is None:
        raise ValueError(f"Missing POINTS marker: {label}")
    count = int(points.group(1))
    offset = points.end() + 12 * count
    for name in ("gid", "ptag", "species"):
        marker = re.match(
            rb"\nSCALARS " + name.encode("ascii") + rb" int\nLOOKUP_TABLE default\n",
            contents[offset:],
        )
        if marker is None:
            raise ValueError(f"Missing {name} marker: {label}")
        offset += marker.end() + 4 * count
    for name in ("deltaf_f0", "deltaf_weight"):
        marker = re.match(
            rb"\nSCALARS " + name.encode("ascii") + rb" float\nLOOKUP_TABLE default\n",
            contents[offset:],
        )
        if marker is None:
            raise ValueError(f"Missing {name} marker: {label}")
        offset += marker.end() + 4 * count
    marker = re.match(rb"\nVECTORS vel float\n", contents[offset:])
    if marker is None:
        raise ValueError(f"Missing velocity marker: {label}")
    offset += marker.end()
    payload = contents[offset:offset + 12 * count]
    if len(payload) != 12 * count:
        raise ValueError(f"Truncated velocity payload: {label}")
    if offset + len(payload) != len(contents):
        raise ValueError(f"Unexpected particle VTK suffix: {label}")
    velocities = [
        struct.unpack(">fff", payload[12 * index:12 * (index + 1)])
        for index in range(count)
    ]
    if any(
        not math.isfinite(component)
        for velocity in velocities
        for component in velocity
    ):
        raise ValueError(f"Non-finite particle velocity: {label}")
    return {"time": time, "cycle": cycle, "count": count, "velocities": velocities}


def require_trusted_gpu_launch(stdout: str) -> dict[str, object]:
    matches = list(GPU_LAUNCH_PATTERN.finditer(stdout))
    if len(matches) != 1 or matches[0].group("rank") != "0":
        raise ValueError("Missing exact rank-zero trusted GPU launch preflight")
    return {
        "rank": 0,
        "host": matches[0].group("host"),
        "rocr_visible_devices": matches[0].group("rocr"),
        "linkage": ["libamdhip64", "libmpi_amd", "libmpi_gtl_hsa"],
    }


def registered_particle_output(inventory: dict[str, dict[str, object]]) -> str:
    expected = "output/pvtk/f1_gpu_relativistic_gyro.prtcl_all.00002.part.vtk"
    paths = sorted(
        path for path in inventory
        if path.startswith("output/pvtk/f1_gpu_relativistic_gyro.prtcl_all.")
        and path.endswith(".part.vtk")
    )
    if paths != [expected]:
        raise ValueError("Gyro particle VTK output set differs from registered cycle 2 artifact")
    return expected


def _analyze_tree(artifact_tree: object) -> dict[str, object]:
    inventory = load_inventory(artifact_tree)
    execution_context = _structured_execution_context(artifact_tree, inventory)

    stdout = read_inventory_bytes(artifact_tree, inventory, "athena_stdout.txt").decode(
        "utf-8"
    )
    stderr = read_inventory_bytes(artifact_tree, inventory, "athena_stderr.txt").decode(
        "utf-8"
    )
    if stderr:
        raise ValueError("Athena stderr is not empty")
    if "physical_mode=paper_test_particle" not in stdout:
        raise ValueError("Missing paper_test_particle runtime identity")
    execution_context["gpu_launch"] = require_trusted_gpu_launch(stdout)
    expected_path = registered_particle_output(inventory)
    snapshot = parse_particle_vtk_bytes(
        read_inventory_bytes(artifact_tree, inventory, expected_path),
        label=expected_path,
    )
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
        "execution_context": execution_context,
        "cycle": snapshot["cycle"],
        "particle_count": snapshot["count"],
        "expected_velocity": expected,
        "max_abs_velocity_error": max_error,
        "velocity_tolerance": VELOCITY_TOLERANCE,
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
