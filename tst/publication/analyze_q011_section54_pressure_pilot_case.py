#!/usr/bin/env python3
"""Verify one frozen trusted-trampoline Q-011 pressure-pilot raw tree.

The emitted descriptor binds the complete raw-tree inventory and the minimal
member set that may enter the aggregate pressure-pilot bundle.  Runtime
metadata remains in the frozen trampoline tree; the aggregate publisher must
not copy it into the analyzer-facing raw bundle.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
from typing import Any, Mapping, Sequence


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
canonical_json_bytes = _ARTIFACT_HELPERS.canonical_json_bytes
load_inventory = _ARTIFACT_HELPERS.load_inventory
read_inventory_bytes = _ARTIFACT_HELPERS.read_inventory_bytes
require_inventory_sha256 = _ARTIFACT_HELPERS.require_inventory_sha256
validate_frontier_mpich_diagnostic_stderr = (
    _ARTIFACT_HELPERS.validate_frontier_mpich_diagnostic_stderr
)
write_result_exclusive = _ARTIFACT_HELPERS.write_result_exclusive


CASE_DESCRIPTOR_PATH = "analysis/analysis.json"
CASE_IDS = ("ps_p0_1p00", "ps_p0_0p05", "ps_p0_0p10", "ps_p0_0p20")
_CASES = (
    ("ps_p0_1p00", 1.0, "1.0"),
    ("ps_p0_0p05", 0.05, "0.05"),
    ("ps_p0_0p10", 0.1, "0.10"),
    ("ps_p0_0p20", 0.2, "0.20"),
)
_CASE_BY_ID = {case_id: (ps_p0, argv) for case_id, ps_p0, argv in _CASES}
_TIMES = (0.0, 15.0, 30.0, 45.0, 60.0)
_PRODUCTS = ("mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all")
_FIXED_OVERRIDES = (
    "mesh/nx1=100",
    "mesh/x1max=1200",
    "mesh/nx2=20",
    "mesh/x2max=240",
    "mesh_refinement/refinement=none",
    "mesh_refinement/num_levels=1",
    "time/tlim=60",
    "time/nlim=4096",
    "time/ndiag=50",
    "problem/ps_enable_curvature_amr=false",
    "problem/ps_feedback_diag_dcycle=50",
    "output1/variable=mhd_w_bcc",
    "output1/id=mhd_w_bcc",
    "output1/dt=15",
    "output2/dt=15",
    "output3/dt=15",
    "output4/dt=15",
    "output5/dt=15",
    "output6/dt=15",
)
EVIDENCE_CLASS = "engineering_calibration_only"
QUALIFICATION_EFFECT = "none_no_sun_bai_claim_no_execution_authorization"
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_FNV1A64_PATTERN = re.compile(r"[0-9a-f]{16}")
_RANK_DIRECTORY_PATTERN = re.compile(r"rank_([0-9]{8})")
_GPU_PREFLIGHT = re.compile(
    r"^PIC trusted GPU launch: rank=(?P<rank>[0-9]+) "
    r"host=(?P<host>\S+) ROCR_VISIBLE_DEVICES=(?P<device>[0-9]+) "
    r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa$"
)
_RUNTIME_ALLOWLIST_KEYS = [
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
_EXPECTED_BASELINE = {
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
_EXPECTED_PROVENANCE_SHA256 = {
    "LOADEDMODULES": "0fa38c7f44ada1f16e61007e577d996267ebbda53b6434a5bc27f97451019eea",
    "_LMFILES_": "dfd37544d54b83574d67af2dd3835b7bc5100f7424bdd2e0a009bb67b9a5174b",
    "MODULEPATH": "825366ca3c91985b4a4ff7ff84ec507bed7d2e433d460a34ba98ee3708b2c9c4",
}


class PressurePilotCaseError(ValueError):
    """Raised when one frozen pressure-pilot case fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PressurePilotCaseError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fnv1a64(payload: bytes) -> str:
    value = 14695981039346656037
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def _reject_constant(value: str) -> None:
    raise PressurePilotCaseError(f"JSON constant is forbidden: {value}")


def _decode_json(payload: bytes, label: str) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PressurePilotCaseError(f"{label}: JSON is not UTF-8") from error
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=_reject_constant,
        )
    except json.JSONDecodeError as error:
        raise PressurePilotCaseError(f"{label}: malformed JSON") from error


def _object(value: object, expected: set[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    mapping = value
    _require(set(mapping) == expected, f"{label}: keys drifted")
    return mapping


def _list(value: object, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected list")
    return value


def _text(value: object, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected nonempty text")
    return value


def _integer(value: object, label: str) -> int:
    _require(type(value) is int, f"{label}: expected integer")
    return value


def _sha256_text(value: object, label: str) -> str:
    text = _text(value, label)
    _require(_SHA256_PATTERN.fullmatch(text) is not None, f"{label}: malformed SHA-256")
    return text


def _relative_path(value: object, label: str) -> str:
    text = _text(value, label)
    path = PurePosixPath(text)
    _require(
        not path.is_absolute()
        and path.as_posix() == text
        and text != "."
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe relative path",
    )
    return text


def _strict_equal(actual: object, expected: object, label: str) -> None:
    _require(type(actual) is type(expected), f"{label}: primitive type drifted")
    if type(expected) is dict:
        _require(set(actual) == set(expected), f"{label}: object keys drifted")
        for key in sorted(expected):
            _strict_equal(actual[key], expected[key], f"{label}/{key}")
    elif type(expected) is list:
        _require(len(actual) == len(expected), f"{label}: list length drifted")
        for index, (left, right) in enumerate(zip(actual, expected)):
            _strict_equal(left, right, f"{label}[{index}]")
    else:
        _require(actual == expected, f"{label}: value drifted")


def _read(
    tree: object,
    inventory: dict[str, dict[str, object]],
    relative: str,
    *,
    nonempty: bool = True,
) -> bytes:
    payload = read_inventory_bytes(tree, inventory, relative)
    if nonempty:
        _require(bool(payload), f"raw pressure-pilot artifact is empty: {relative}")
    return payload


def _binding(path: str, payload: bytes) -> dict[str, str]:
    return {"path": path, "sha256": _sha256(payload)}


def _member(source_path: str, target_path: str, payload: bytes) -> dict[str, object]:
    return {
        "source_path": source_path,
        "path": target_path,
        "size": len(payload),
        "sha256": _sha256(payload),
    }


def _snapshot_source_paths(case_id: str, index: int) -> dict[str, str]:
    suffix = f"{index:05d}"
    return {
        "mhd_w_bcc": f"output/bin/{case_id}.mhd_w_bcc.{suffix}.bin",
        "bmag": f"output/bin/{case_id}.bmag.{suffix}.bin",
        "prtcl_jx": f"output/bin/{case_id}.prtcl_jx.{suffix}.bin",
        "j2": f"output/bin/{case_id}.j2.{suffix}.bin",
        "prtcl_all": f"output/pvtk/{case_id}.prtcl_all.{suffix}.part.vtk",
    }


def _snapshot_target_paths(case_id: str, index: int) -> dict[str, str]:
    return {
        name: f"cases/{case_id}/{source.removeprefix('output/')}"
        for name, source in _snapshot_source_paths(case_id, index).items()
    }


def _action_id(case_id: str) -> str:
    _require(case_id in _CASE_BY_ID, f"unregistered pressure-pilot case: {case_id}")
    return f"q011-pressure-{case_id.replace('_', '-')}"


def _allowlist(
    tree: object,
    inventory: dict[str, dict[str, object]],
    case_id: str,
) -> tuple[str, bytes, dict[str, str]]:
    relative = f"{_action_id(case_id)}.environment.allowlist.txt"
    payload = _read(tree, inventory, relative)
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PressurePilotCaseError("runtime allowlist is not UTF-8") from error
    _require("\r" not in text and text.endswith("\n"), "runtime allowlist must use terminated LF records")
    values = {}
    for record in text.splitlines():
        key, separator, value = record.partition("=")
        _require(bool(separator and key and value), "runtime allowlist contains a malformed record")
        _require(key not in values, f"runtime allowlist contains duplicate key: {key}")
        values[key] = value
    _require(list(values) == _RUNTIME_ALLOWLIST_KEYS, "runtime allowlist key set or order is not authorized")
    for key, expected in _EXPECTED_BASELINE.items():
        _require(values[key] == expected, f"runtime allowlist {key} differs from the baseline")
    for key, expected in _EXPECTED_PROVENANCE_SHA256.items():
        _require(
            _sha256(values[key].encode("utf-8")) == expected,
            f"runtime allowlist {key} differs from reviewed provenance",
        )
    return relative, payload, values


def _trusted_gpu_preflight(stdout: str) -> list[dict[str, object]]:
    bindings = []
    for line in stdout.splitlines():
        match = _GPU_PREFLIGHT.fullmatch(line)
        if match is not None:
            bindings.append(
                {
                    "rank": int(match.group("rank")),
                    "host": match.group("host"),
                    "rocr_visible_device": int(match.group("device")),
                }
            )
    ranks = sorted(int(binding["rank"]) for binding in bindings)
    devices = sorted(int(binding["rocr_visible_device"]) for binding in bindings)
    _require(bool(bindings), "trusted GPU launch preflight is absent")
    _require(ranks == list(range(len(bindings))), "trusted GPU launch ranks are not contiguous")
    _require(devices == list(range(len(bindings))), "trusted GPU launch devices are not contiguous")
    _require(
        len({str(binding["host"]) for binding in bindings}) == 1,
        "trusted GPU launch spans more than one host",
    )
    return sorted(bindings, key=lambda binding: int(binding["rank"]))


def _marker(payload: bytes, artifact: bytes, label: str) -> None:
    try:
        lines = payload.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise PressurePilotCaseError(f"{label}: completion marker is not ASCII") from error
    _require(
        len(lines) == 3 and lines[0] == "ATHENAK_RESTART_COMPLETE_V1",
        f"{label}: completion marker is malformed",
    )
    _require(lines[1] == f"size={len(artifact)}", f"{label}: completion size drifted")
    _require(lines[2] == f"fnv1a64={_fnv1a64(artifact)}", f"{label}: completion digest drifted")


def _restart_members(
    tree: object,
    inventory: dict[str, dict[str, object]],
    case_id: str,
    index: int,
) -> tuple[list[str], dict[str, object]]:
    suffix = f"{index:05d}"
    raw_manifest = f"output/rst/{case_id}.{suffix}.rst.manifest"
    raw_manifest_marker = raw_manifest + ".complete"
    manifest_payload = _read(tree, inventory, raw_manifest)
    marker_payload = _read(tree, inventory, raw_manifest_marker)
    _marker(marker_payload, manifest_payload, raw_manifest_marker)
    decoded = _object(
        _decode_json(manifest_payload, raw_manifest),
        {"schema", "members"},
        raw_manifest,
    )
    _strict_equal(decoded["schema"], "ATHENAK_RESTART_MANIFEST_V1", f"{raw_manifest}/schema")
    raw_members = _list(decoded["members"], f"{raw_manifest}/members")
    _require(bool(raw_members), f"{raw_manifest}: restart manifest has no members")
    source_paths = [raw_manifest, raw_manifest_marker]
    bundle_members = []
    ranked = []
    expected_name = f"{case_id}.{suffix}.rst"
    for member_index, raw_member in enumerate(raw_members):
        label = f"{raw_manifest}/members[{member_index}]"
        item = _object(raw_member, {"path", "size", "fnv1a64"}, label)
        relative = PurePosixPath(_relative_path(item["path"], f"{label}/path"))
        size = _integer(item["size"], f"{label}/size")
        digest = _text(item["fnv1a64"], f"{label}/fnv1a64")
        _require(size >= 0 and _FNV1A64_PATTERN.fullmatch(digest) is not None, f"{label}: digest is malformed")
        _require(relative.parts[0] == "rst" and relative.name == expected_name, f"{label}: member path drifted")
        if relative.parent == PurePosixPath("rst"):
            rank = None
        else:
            _require(len(relative.parts) == 3, f"{label}: member layout drifted")
            match = _RANK_DIRECTORY_PATTERN.fullmatch(relative.parent.name)
            _require(match is not None, f"{label}: rank directory drifted")
            rank = int(match.group(1))
        ranked.append(rank)
        raw_payload = f"output/{relative.as_posix()}"
        raw_marker = raw_payload + ".complete"
        payload = _read(tree, inventory, raw_payload)
        completion = _read(tree, inventory, raw_marker)
        _marker(completion, payload, raw_marker)
        _require(len(payload) == size and _fnv1a64(payload) == digest, f"{label}: member digest drifted")
        source_paths.extend((raw_payload, raw_marker))
        if index == len(_TIMES) - 1:
            target = f"cases/{case_id}/{relative.as_posix()}"
            bundle_members.append(
                {
                    "artifact": _binding(target, payload),
                    "complete": _binding(target + ".complete", completion),
                }
            )
    shared = all(rank is None for rank in ranked)
    _require(shared or all(rank is not None for rank in ranked), f"{raw_manifest}: restart layout is mixed")
    if shared:
        _require(len(ranked) == 1, f"{raw_manifest}: shared restart must contain one member")
    else:
        _require(ranked == list(range(len(ranked))), f"{raw_manifest}: rank members are not contiguous")
    terminal = {}
    if index == len(_TIMES) - 1:
        target_manifest = f"cases/{case_id}/rst/{case_id}.{suffix}.rst.manifest"
        terminal = {
            "time": _TIMES[index],
            "manifest": _binding(target_manifest, manifest_payload),
            "manifest_complete": _binding(target_manifest + ".complete", marker_payload),
            "members": bundle_members,
        }
    return source_paths, terminal


def _analyze_tree(tree: object, case_id: str) -> dict[str, object]:
    _require(case_id in _CASE_BY_ID, f"unregistered pressure-pilot case: {case_id}")
    inventory = load_inventory(tree)
    inventory_payload = tree.read("artifact_inventory.json")
    inventory_sha256 = _sha256(inventory_payload)
    ps_p0, argv_value = _CASE_BY_ID[case_id]
    allowlist_path, allowlist_payload, allowlist = _allowlist(tree, inventory, case_id)
    stdout_payload = _read(tree, inventory, "athena_stdout.txt")
    stdout_checksum_payload = _read(tree, inventory, "athena_stdout.sha256")
    _require(
        stdout_checksum_payload == (_sha256(stdout_payload) + "\n").encode("ascii"),
        "athena_stdout.sha256 does not bind Athena stdout",
    )
    stderr_payload = _read(tree, inventory, "athena_stderr.txt")
    try:
        stdout = stdout_payload.decode("utf-8")
        stderr = stderr_payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PressurePilotCaseError("Athena stdout or stderr is not UTF-8") from error
    validate_frontier_mpich_diagnostic_stderr(stderr)
    rank_gpu_bindings = _trusted_gpu_preflight(stdout)
    error_history_path = f"output/{case_id}-errs.dat"
    error_history_payload = _read(tree, inventory, error_history_path)

    expected_paths = {
        allowlist_path,
        "athena_stdout.txt",
        "athena_stdout.sha256",
        "athena_stderr.txt",
        error_history_path,
    }
    bundle_members = [
        _member("athena_stdout.txt", f"cases/{case_id}/stdout.txt", stdout_payload)
    ]
    snapshots = []
    for index, time in enumerate(_TIMES):
        snapshot = {"time": time}
        for name, source_path in _snapshot_source_paths(case_id, index).items():
            payload = _read(tree, inventory, source_path)
            target_path = _snapshot_target_paths(case_id, index)[name]
            expected_paths.add(source_path)
            snapshot[name] = _binding(target_path, payload)
            bundle_members.append(_member(source_path, target_path, payload))
        snapshots.append(snapshot)

    terminal_restart = {}
    for index in range(len(_TIMES)):
        restart_paths, terminal = _restart_members(tree, inventory, case_id, index)
        expected_paths.update(restart_paths)
        if terminal:
            terminal_restart = terminal
            terminal_sources = restart_paths
    _require(bool(terminal_restart), f"{case_id}: terminal restart descriptor is absent")
    for source_path in terminal_sources:
        payload = _read(tree, inventory, source_path)
        bundle_members.append(
            _member(source_path, f"cases/{case_id}/{source_path.removeprefix('output/')}", payload)
        )
    _require(set(inventory) == expected_paths, f"{case_id}: raw trampoline tree closure drifted")
    tree.require_tree_closure()
    return {
        "schema_version": 1,
        "record_type": "q011_section54_pressure_pilot_verified_raw_case",
        "evidence_class": EVIDENCE_CLASS,
        "qualification_effect": QUALIFICATION_EFFECT,
        "launch_contract": "trusted_trampoline_athena_argv_v1",
        "case_id": case_id,
        "ps_p0": ps_p0,
        "argv_value": argv_value,
        "artifact_inventory_sha256": inventory_sha256,
        "runtime_artifacts": {
            allowlist_path: _sha256(allowlist_payload),
            "athena_stdout.txt": _sha256(stdout_payload),
            "athena_stdout.sha256": _sha256(stdout_checksum_payload),
            "athena_stderr.txt": _sha256(stderr_payload),
            error_history_path: _sha256(error_history_payload),
        },
        "runtime_profile": allowlist["PIC_FRONTIER_PROFILE"],
        "parallel_ranks": len(rank_gpu_bindings),
        "rank_gpu_bindings": rank_gpu_bindings,
        "manifest_case": {
            "case_id": case_id,
            "ps_p0": ps_p0,
            "overrides": [*_FIXED_OVERRIDES, f"problem/ps_p0={argv_value}"],
            "snapshots": snapshots,
            "stdout": _binding(f"cases/{case_id}/stdout.txt", stdout_payload),
            "terminal_restart": terminal_restart,
        },
        "bundle_members": sorted(bundle_members, key=lambda member: str(member["path"])),
    }


def analyze_case_tree(artifact_dir: str | Path, case_id: str) -> dict[str, object]:
    """Return the strict descriptor for one frozen raw case without publishing it."""
    with StructuredArtifactTree(Path(artifact_dir)) as tree:
        return _analyze_tree(tree, case_id)


def publish_case_descriptor(artifact_dir: str | Path, case_id: str) -> dict[str, object]:
    """Publish one immutable case descriptor into the trampoline analysis directory."""
    with StructuredArtifactTree(Path(artifact_dir)) as tree:
        descriptor = _analyze_tree(tree, case_id)
        write_result_exclusive(tree, CASE_DESCRIPTOR_PATH, descriptor)
        return descriptor


def _descriptor_bytes(tree: object) -> bytes:
    tree.require_tree_closure()
    tree.require_analysis_identity()
    descriptor = os.open(
        "analysis.json",
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=tree.analysis_fd,
    )
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and not before.st_mode & 0o222,
            "published raw-case descriptor is not a read-only regular file",
        )
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            payload = stream.read()
        after = os.fstat(descriptor)
        entry = os.stat("analysis.json", dir_fd=tree.analysis_fd, follow_symlinks=False)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        _require(
            all(getattr(before, name) == getattr(after, name) for name in stable)
            and (entry.st_dev, entry.st_ino) == (after.st_dev, after.st_ino)
            and len(payload) == after.st_size,
            "published raw-case descriptor changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


def verify_published_case_descriptor(
    tree: object,
    case_id: str,
    expected_descriptor_sha256: str,
) -> dict[str, object]:
    """Recompute and verify one descriptor while retaining the caller's root pin."""
    _require(
        _SHA256_PATTERN.fullmatch(expected_descriptor_sha256) is not None,
        "expected raw-case descriptor SHA-256 is malformed",
    )
    load_inventory(tree)
    payload = _descriptor_bytes(tree)
    _require(_sha256(payload) == expected_descriptor_sha256, "raw-case descriptor SHA-256 drifted")
    decoded = _decode_json(payload, "raw-case descriptor")
    _require(type(decoded) is dict, "raw-case descriptor must be an object")
    _require(payload == canonical_json_bytes(decoded), "raw-case descriptor is not canonical JSON")
    recomputed = _analyze_tree(tree, case_id)
    _strict_equal(decoded, recomputed, "raw-case descriptor")
    tree.require_tree_closure()
    tree.require_analysis_identity()
    _require(_descriptor_bytes(tree) == payload, "raw-case descriptor changed after verification")
    return recomputed


def _descriptor_sha256(descriptor: Mapping[str, object]) -> str:
    return _sha256(canonical_json_bytes(dict(descriptor)))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    parser.add_argument("--case-id", required=True, choices=CASE_IDS)
    parser.add_argument("--artifact-dir-fd", type=int)
    parser.add_argument("--verify-artifact-inventory-sha256")
    parser.add_argument("--verify-result-sha256")
    args = parser.parse_args(argv)
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
    with StructuredArtifactTree(
        args.artifact_dir, inherited_root_fd=args.artifact_dir_fd
    ) as tree:
        if args.verify_artifact_inventory_sha256 is not None:
            require_inventory_sha256(tree, args.verify_artifact_inventory_sha256)
        descriptor = _analyze_tree(tree, args.case_id)
        descriptor_sha256 = _descriptor_sha256(descriptor)
        if args.verify_result_sha256 is not None:
            require_inventory_sha256(tree, args.verify_artifact_inventory_sha256)
            if descriptor_sha256 != args.verify_result_sha256:
                raise ValueError("Offline analysis recomputation differs from bound result")
            return 0
        write_result_exclusive(tree, CASE_DESCRIPTOR_PATH, descriptor)
    print(descriptor_sha256)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
