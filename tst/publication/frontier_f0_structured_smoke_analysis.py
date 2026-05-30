#!/usr/bin/env python3
"""Validate artifacts from the structured Frontier F0 parser-contract smoke."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat


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


def _read_regular(path: Path, *, require_nonempty: bool) -> bytes:
    metadata = path.lstat()
    if not stat.S_ISREG(metadata.st_mode):
        raise ValueError(f"F0 artifact is not a regular file: {path}")
    data = path.read_bytes()
    if require_nonempty and not data:
        raise ValueError(f"F0 artifact is empty: {path}")
    return data


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _parse_allowlist(path: Path) -> tuple[bytes, dict[str, str]]:
    data = _read_regular(path, require_nonempty=True)
    if path.stat().st_mode & 0o222:
        raise ValueError("Runtime allowlist must be read-only")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("Runtime allowlist is not UTF-8") from error
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


def analyze(artifact_dir: Path) -> dict[str, object]:
    allowlist_data, allowlist = _parse_allowlist(
        artifact_dir / "f0-parser.environment.allowlist.txt"
    )
    stdout = _read_regular(artifact_dir / "athena_stdout.txt", require_nonempty=True)
    stderr = _read_regular(artifact_dir / "athena_stderr.txt", require_nonempty=False)
    return {
        "schema_version": 1,
        "status": "pass",
        "evidence_class": "frontier_f0_admission_smoke_candidate",
        "launch_contract": "trusted_trampoline_athena_argv_v1",
        "parser_smoke": "pass",
        "runtime_profile": allowlist["PIC_FRONTIER_PROFILE"],
        "runtime_artifacts": {
            "f0-parser.environment.allowlist.txt": _sha256(allowlist_data),
            "athena_stdout.txt": _sha256(stdout),
            "athena_stderr.txt": _sha256(stderr),
        },
    }


def _write_result(path: Path, result: dict[str, object]) -> None:
    data = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode("utf-8")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(data)
            stream.flush()
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    args = parser.parse_args()
    _write_result(args.artifact_dir / "analysis.json", analyze(args.artifact_dir))


if __name__ == "__main__":
    main()
