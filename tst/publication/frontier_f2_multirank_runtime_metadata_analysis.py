#!/usr/bin/env python3
"""Validate the bounded Frontier F2 eight-rank runtime-metadata slice."""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import importlib.util
import os
import re
from pathlib import Path
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


EXPECTED_RANKS = set(range(8))
EXPECTED_DEVICES = set(range(8))
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
GPU_PREFLIGHT = re.compile(
    r"^PIC trusted GPU launch: rank=(?P<rank>[0-9]+) "
    r"host=(?P<host>\S+) ROCR_VISIBLE_DEVICES=(?P<device>[0-9]+) "
    r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa$"
)
RUNTIME_MODEL_PREFIX = "PIC runtime model:"
EXPECTED_RUNTIME_MODEL = (
    "PIC runtime model: physical_mode=extended_mhd_pic state=momentum_p_over_m "
    "C=3 background=coupled feedback=coupled induction=ideal_mhd_only "
    "deposition=tsc deltaf=off deltaf_adapt=off deltaf_adapt_interval=0 "
    "expanding_box=off expansion_law=linear wave_damping=off nu_in=0 "
    "lb_cost_per_particle=0 max_cell_cross=2 theta_max=0.3 restart_schema=7"
)
EXPECTED_RANK_COUNT = "Number of parallel ranks = 8"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _decode(
    artifact_tree: object,
    inventory: dict[str, dict[str, object]],
    relative: str,
    *,
    require_nonempty: bool,
) -> tuple[bytes, str]:
    data = read_inventory_bytes(artifact_tree, inventory, relative)
    if require_nonempty and not data:
        raise ValueError(f"F2 runtime artifact is empty: {relative}")
    try:
        return data, data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"F2 runtime artifact is not UTF-8: {relative}") from error


def _parse_allowlist(
    artifact_tree: object,
    inventory: dict[str, dict[str, object]],
) -> tuple[bytes, dict[str, str]]:
    data, text = _decode(
        artifact_tree,
        inventory,
        "f2-runtime-metadata.environment.allowlist.txt",
        require_nonempty=True,
    )
    if "\r" in text or not text.endswith("\n"):
        raise ValueError("Runtime allowlist must use terminated LF records")
    values = {}
    for record in text.splitlines():
        key, separator, value = record.partition("=")
        if not separator or not key or not value:
            raise ValueError("Runtime allowlist contains a malformed record")
        if key in values:
            raise ValueError(f"Runtime allowlist contains duplicate key: {key}")
        values[key] = value
    if list(values) != RUNTIME_ALLOWLIST_KEYS:
        raise ValueError("Runtime allowlist key set or order is not authorized")
    for key, expected in EXPECTED_BASELINE.items():
        if values[key] != expected:
            raise ValueError(f"Runtime allowlist {key} differs from the baseline")
    for key, expected in EXPECTED_PROVENANCE_SHA256.items():
        if _sha256(values[key].encode("utf-8")) != expected:
            raise ValueError(f"Runtime allowlist {key} differs from reviewed provenance")
    return data, values


def _rank_gpu_bindings(stdout: str) -> list[dict[str, object]]:
    bindings = []
    for line in stdout.splitlines():
        match = GPU_PREFLIGHT.fullmatch(line)
        if match is not None:
            bindings.append(
                {
                    "rank": int(match.group("rank")),
                    "host": match.group("host"),
                    "rocr_visible_device": int(match.group("device")),
                }
            )
    ranks = {binding["rank"] for binding in bindings}
    devices = {binding["rocr_visible_device"] for binding in bindings}
    hosts = {binding["host"] for binding in bindings}
    if len(bindings) != len(EXPECTED_RANKS) or ranks != EXPECTED_RANKS:
        raise ValueError("F2 trusted GPU preflight does not contain exactly ranks 0 through 7")
    if devices != EXPECTED_DEVICES:
        raise ValueError("F2 trusted GPU preflight does not bind exactly devices 0 through 7")
    if len(hosts) != 1:
        raise ValueError("F2 trusted GPU preflight does not remain on one Frontier node")
    return sorted(bindings, key=lambda binding: int(binding["rank"]))


def _require_runtime_model(stdout: str) -> str:
    models = [
        line for line in stdout.splitlines() if line.startswith(RUNTIME_MODEL_PREFIX)
    ]
    if len(models) != 1:
        raise ValueError("F2 Athena stdout does not contain exactly one runtime-model line")
    if models[0] != EXPECTED_RUNTIME_MODEL:
        raise ValueError("F2 Athena runtime model differs from the reviewed identity")
    return models[0]


def _analyze_tree(artifact_tree: object) -> dict[str, object]:
    inventory = load_inventory(artifact_tree)
    allowlist_data, allowlist = _parse_allowlist(artifact_tree, inventory)
    stdout_data, stdout = _decode(
        artifact_tree, inventory, "athena_stdout.txt", require_nonempty=True
    )
    stderr_data, stderr = _decode(
        artifact_tree, inventory, "athena_stderr.txt", require_nonempty=True
    )
    validate_frontier_mpich_diagnostic_stderr(stderr)
    if stdout.splitlines().count(EXPECTED_RANK_COUNT) != 1:
        raise ValueError("F2 Athena stdout does not report exactly one eight-rank record")
    bindings = _rank_gpu_bindings(stdout)
    runtime_model = _require_runtime_model(stdout)
    return {
        "schema_version": 1,
        "status": "pass",
        "evidence_class": "frontier_f2_clean_candidate_multirank_runtime_metadata",
        "launch_contract": "trusted_trampoline_athena_argv_v1",
        "parser_runtime_smoke": "pass",
        "runtime_profile": allowlist["PIC_FRONTIER_PROFILE"],
        "parallel_ranks": 8,
        "hosts": sorted({str(binding["host"]) for binding in bindings}),
        "rank_gpu_bindings": bindings,
        "runtime_model": runtime_model,
        "runtime_artifacts": {
            "f2-runtime-metadata.environment.allowlist.txt": _sha256(allowlist_data),
            "athena_stdout.txt": _sha256(stdout_data),
            "athena_stderr.txt": _sha256(stderr_data),
        },
    }


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
        write_result_exclusive(artifact_tree, "analysis/analysis.json", result)
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
