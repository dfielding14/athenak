#!/usr/bin/env python3
"""Registered-execution admission for the Q-043 deposited-current oracle."""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tarfile
from typing import Any, Mapping, Sequence

from tst.publication import analyze_q011_section54_outputs as binary
from tst.publication import (
    q043_bell_current_volume_aware_deposited_current_oracle as oracle,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
AUTHORIZED_ORION_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
RUN_NAMESPACE = "runs/q043_registered_execution_raw_oracle_successor_v1"
SCHEMA_VERSION = 1
SUCCESSOR_ID = "q043_registered_execution_raw_oracle_qualification_successor_v1"
CASE_RECORD_TYPE = "q043_registered_execution_raw_oracle_case_admission"
MATRIX_RECORD_TYPE = "q043_registered_execution_raw_oracle_matrix_qualification"
LAUNCH_CONTRACT_RECORD_TYPE = "q043_registered_execution_launch_contract"
EXECUTION_RECEIPT_RECORD_TYPE = "q043_reconciled_registered_execution_receipt"
TERMINAL_RECEIPT_RECORD_TYPE = "q043_registered_execution_terminal_receipt"
SOURCE_LOCAL_CASE_RECORD_TYPE = (
    "q043_bell_current_volume_aware_deposited_current_raw_case_oracle"
)
SOURCE_LOCAL_MATRIX_RECORD_TYPE = (
    "q043_bell_current_volume_aware_deposited_current_raw_matrix_oracle"
)
REQUIRED_CYCLES = (0, 1)
REQUIRED_FIELDS = tuple(oracle.FIELDS)
EXPECTED_CASE_COUNT = 132
REQUIRED_DIMENSIONS = (1, 2, 3)
REQUIRED_RESOLUTIONS = ("coarse", "fine")
REQUIRED_PPC_VALUES = (1, 4)
REQUIRED_ARTIFICIAL_C_OVER_V_CR_VALUES = (100, 1000, 10000)
REQUIRED_DECOMPOSITIONS_BY_DIMENSION = {
    1: ("single", "split_x1"),
    2: ("single", "split_x1", "split_x2", "split_x1x2"),
    3: ("single", "split_x1", "split_x2", "split_x3", "split_xyz"),
}
REQUIRED_DECOMPOSITION_GRIDS = {
    "single": (1, 1, 1),
    "split_x1": (2, 1, 1),
    "split_x2": (1, 2, 1),
    "split_x1x2": (2, 2, 1),
    "split_x3": (1, 1, 2),
    "split_xyz": (2, 2, 2),
}
SOURCE_LOCAL_RAW_BINDING_KEYS = {"field", "shard_index", "path", "size", "sha256"}
SOURCE_LOCAL_CASE_RESULT_KEYS = {
    "schema_version",
    "record_type",
    "campaign_id",
    "case_id",
    "qualification_effect",
    "launch_authorized",
    "scientific_claim_authorized",
    "publication_authorized",
    "passed",
    "source_local_oracle_check_pass",
    "case_contract",
    "raw_provenance",
    "cross_field_runtime_metadata",
    "representation_derived_absolute_tolerance",
    "configured_volume_mean_j_over_c",
    "measured_volume_mean_j_over_c_vector",
    "measured_guide_projected_volume_mean_j_over_c",
    "measured_volume_mean_prtcl_rho",
    "maximum_transverse_current",
    "maximum_current_spatial_nonuniformity",
    "maximum_rho_spatial_nonuniformity",
    "tolerance_scope",
}
SOURCE_LOCAL_MATRIX_RESULT_KEYS = {
    "schema_version",
    "record_type",
    "campaign_id",
    "qualification_effect",
    "launch_authorized",
    "scientific_claim_authorized",
    "publication_authorized",
    "passed",
    "source_local_oracle_check_pass",
    "case_count",
    "dimensions_verified",
    "resolutions_verified",
    "ppc_verified",
    "decompositions_verified",
    "artificial_c_over_v_cr_verified",
    "required_output_cycle",
    "output_dcycle",
    "output_timing",
    "cross_field_runtime_metadata",
    "initial_state_verified",
    "artificial_c_in_formula",
    "configured_current_formula",
    "maximum_cross_matrix_projected_current_spread",
    "cases",
}
REQUIRED_CANDIDATE_SOURCE_PATHS = frozenset(
    (
        "src/CMakeLists.txt",
        "src/pgen/pgen.cpp",
        "src/pgen/pgen.hpp",
        "src/pgen/tests/q043_bell_current_volume_aware.cpp",
        "tst/publication/q043_bell_current_volume_aware_deposited_current_oracle.py",
        "tst/publication/q043_registered_execution_raw_oracle_qualification_successor_v1.py",
        "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle/"
        "deck_manifest.json",
    )
)
OUTPUT_BLOCKS = tuple(f"output{index}" for index in range(1, len(REQUIRED_FIELDS) + 1))
OUTPUT_BOOKKEEPING_KEYS = frozenset(("file_number", "last_time"))
QUALIFICATION_EFFECT = (
    "registered_execution_raw_oracle_prerequisite_only_no_launch_no_policy_"
    "no_q023_no_q019_no_science_no_publication_authority"
)

_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
_JOB_ID = re.compile(r"[1-9][0-9]*")
_TIME_CYCLE = re.compile(r"^time=([^\s]+) cycle=(0|[1-9][0-9]*)$", re.MULTILINE)
_LIMITS = re.compile(r"^tlim=([^\s]+) nlim=([^\s]+)$", re.MULTILINE)
_RANK_EVIDENCE = re.compile(
    r"^Q043_REGISTERED_EXECUTION case_id=([^\s]+) mpi_world_size=([1-9][0-9]*) "
    r"rank_ids=([0-9]+(?:,[0-9]+)*)$",
    re.MULTILINE,
)
_EXIT_EVIDENCE = re.compile(
    r"^Q043_REGISTERED_EXECUTION_EXIT exit_code=([0-9]+) signal=([0-9]+)$",
    re.MULTILINE,
)


class AdmissionError(ValueError):
    """Reject incomplete, mixed-lineage, synthetic, or authority-bearing evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AdmissionError(message)


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(set(value) == keys, f"{label}: keys drifted")
    return value


def _array(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected array")
    return value


def _text(value: object, *, label: str) -> str:
    _require(
        type(value) is str and bool(value) and not any(char.isspace() for char in value),
        f"{label}: expected nonempty whitespace-free text",
    )
    return value


def _absolute_path(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    path = PurePosixPath(text)
    _require(
        path.is_absolute()
        and path.as_posix() == text
        and all(part not in ("", ".", "..") for part in path.parts),
        f"{label}: expected canonical absolute path",
    )
    return text


def _relative_path(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    path = PurePosixPath(text)
    _require(
        not path.is_absolute()
        and path.as_posix() == text
        and text != "."
        and all(part not in ("", ".", "..") for part in path.parts),
        f"{label}: unsafe root-relative path",
    )
    return text


def _sha256(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _git_commit(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _GIT_COMMIT.fullmatch(value) is not None,
        f"{label}: malformed Git commit",
    )
    return value


def _positive_integer(value: object, *, label: str) -> int:
    _require(type(value) is int and value > 0, f"{label}: expected positive integer")
    return value


def _nonnegative_integer(value: object, *, label: str) -> int:
    _require(
        type(value) is int and value >= 0,
        f"{label}: expected nonnegative integer",
    )
    return value


def _canonical_json_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise AdmissionError("record is not canonical finite JSON") from error


def canonical_sha256(value: object) -> str:
    """Return the deterministic SHA-256 of one finite JSON value."""
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _strict_equal(lvalue, rvalue) for lvalue, rvalue in zip(left, right)
        )
    return left == right


def _pairs_to_unique_object(
    pairs: Sequence[tuple[str, object]], *, label: str
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        _require(key not in result, f"{label}: duplicate JSON key {key!r}")
        result[key] = value
    return result


def _json_payload(payload: bytes, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=lambda pairs: _pairs_to_unique_object(pairs, label=label),
            parse_constant=lambda token: (_ for _ in ()).throw(
                AdmissionError(f"{label}: non-finite JSON constant {token}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AdmissionError(f"{label}: invalid UTF-8 JSON") from error
    _require(type(value) is dict, f"{label}: expected JSON object")
    return value


def _canonical_root(value: str | Path, *, label: str) -> Path:
    path = Path(os.path.abspath(value))
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise AdmissionError(f"{label}: root is unavailable") from error
    _require(resolved == path and path.is_dir(), f"{label}: root uses a symlink alias")
    return path


def _canonical_below(
    value: object, root: Path, *, label: str, directory: bool = False
) -> Path:
    path = Path(_absolute_path(value, label=label))
    try:
        resolved = path.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise AdmissionError(f"{label}: path escaped or is unavailable") from error
    _require(resolved == path, f"{label}: path uses a symlink alias")
    _require(
        path.is_dir() if directory else path.is_file(),
        f"{label}: path has the wrong file type",
    )
    return path


def _source_member(source_root: Path, relative: str, *, label: str) -> Path:
    relative = _relative_path(relative, label=label)
    try:
        path = (source_root / relative).resolve(strict=True)
        path.relative_to(source_root)
    except (OSError, ValueError) as error:
        raise AdmissionError(f"{label}: source member escaped or is unavailable") from error
    _require(path == source_root / relative, f"{label}: source member uses a symlink alias")
    return path


def _read_stable_regular_file(
    path: Path,
    *,
    label: str,
    expected_sha256: str | None = None,
    expected_byte_count: int | None = None,
    executable: bool = False,
) -> tuple[bytes, dict[str, int]]:
    descriptor = None
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        before = os.fstat(descriptor)
        chunks = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
    except OSError as error:
        raise AdmissionError(f"{label}: bound file is unavailable") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
    identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_nlink,
    )
    _require(
        stat.S_ISREG(before.st_mode)
        and not path.is_symlink()
        and before.st_nlink == 1
        and identity
        == (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_nlink,
        )
        == (
            current.st_dev,
            current.st_ino,
            current.st_size,
            current.st_mtime_ns,
            current.st_nlink,
        ),
        f"{label}: bound file is not one stable unaliased regular file",
    )
    _require(
        not executable or bool(before.st_mode & 0o111),
        f"{label}: bound executable lacks an execute bit",
    )
    digest = hashlib.sha256(payload).hexdigest()
    if expected_sha256 is not None:
        _require(digest == expected_sha256, f"{label}: SHA-256 drifted")
    if expected_byte_count is not None:
        _require(len(payload) == expected_byte_count, f"{label}: byte count drifted")
    return payload, {"device": before.st_dev, "inode": before.st_ino}


def _binding(value: object, *, label: str) -> dict[str, object]:
    record = _object(value, {"path", "sha256", "byte_count"}, label=label)
    return {
        "path": _absolute_path(record["path"], label=f"{label}/path"),
        "sha256": _sha256(record["sha256"], label=f"{label}/sha256"),
        "byte_count": _positive_integer(
            record["byte_count"], label=f"{label}/byte_count"
        ),
    }


def _verified_binding(
    value: object,
    *,
    root: Path,
    label: str,
    executable: bool = False,
) -> tuple[dict[str, object], bytes]:
    binding = _binding(value, label=label)
    path = _canonical_below(binding["path"], root, label=f"{label}/path")
    payload, identity = _read_stable_regular_file(
        path,
        label=label,
        expected_sha256=str(binding["sha256"]),
        expected_byte_count=int(binding["byte_count"]),
        executable=executable,
    )
    return {**binding, "filesystem_identity": identity}, payload


def _public_binding(value: Mapping[str, object]) -> dict[str, object]:
    return {key: value[key] for key in ("path", "sha256", "byte_count")}


def _expected_case_axis_rows() -> list[tuple[object, ...]]:
    return [
        (dimension, resolution, ppc, decomposition, artificial_c_over_v_cr)
        for dimension in REQUIRED_DIMENSIONS
        for resolution in REQUIRED_RESOLUTIONS
        for ppc in REQUIRED_PPC_VALUES
        for decomposition in REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension]
        for artificial_c_over_v_cr in REQUIRED_ARTIFICIAL_C_OVER_V_CR_VALUES
    ]


def _validate_case_decomposition_contract(case: Mapping[str, object]) -> None:
    case_id = case.get("case_id")
    dimension = case.get("dimension")
    decomposition = case.get("decomposition")
    _require(
        type(case_id) is str
        and type(dimension) is int
        and dimension in REQUIRED_DECOMPOSITIONS_BY_DIMENSION
        and type(decomposition) is str
        and decomposition in REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension],
        "Q043 case has invalid dimension or decomposition identity",
    )
    _require(
        type(case.get("ppc")) is int and int(case["ppc"]) > 0,
        f"{case_id}: PPC must be a positive integer",
    )
    grid = case.get("meshblock_grid")
    global_nx = case.get("global_nx")
    meshblock_nx = case.get("meshblock_nx")
    expected_grid = REQUIRED_DECOMPOSITION_GRIDS[decomposition]
    _require(
        type(grid) is list
        and tuple(grid) == expected_grid
        and type(global_nx) is list
        and type(meshblock_nx) is list
        and len(global_nx) == len(meshblock_nx) == 3
        and all(type(value) is int and value > 0 for value in global_nx)
        and all(type(value) is int and value > 0 for value in meshblock_nx)
        and all(
            global_count == block_count * grid_count
            for global_count, block_count, grid_count in zip(
                global_nx, meshblock_nx, expected_grid
            )
        ),
        f"{case_id}: exact MeshBlock decomposition drifted",
    )
    expected_partitioned_axes = [
        index for index, count in enumerate(expected_grid, 1) if count > 1
    ]
    _require(
        case.get("partitioned_axes") == expected_partitioned_axes
        and type(case.get("mpi_ranks")) is int
        and int(case["mpi_ranks"]) == math.prod(expected_grid),
        f"{case_id}: dimension-valid one-MeshBlock-per-rank contract drifted",
    )


def _validate_exact_case_matrix(cases: Sequence[Mapping[str, object]]) -> None:
    _require(
        len(cases) == EXPECTED_CASE_COUNT
        and getattr(oracle, "EXPECTED_CASE_COUNT", None) == EXPECTED_CASE_COUNT,
        "Q043 registered admission requires the exact 132-case foundation",
    )
    observed_rows = []
    observed_ids = []
    for index, case in enumerate(cases):
        _require(type(case) is dict, f"Q043 case matrix row {index} is malformed")
        _validate_case_decomposition_contract(case)
        row = (
            case.get("dimension"),
            case.get("resolution"),
            case.get("ppc"),
            case.get("decomposition"),
            case.get("artificial_c_over_v_cr"),
        )
        observed_rows.append(row)
        expected_id = (
            f"q043-current-oracle-d{row[0]}-{row[1]}-ppc{row[2]}-"
            f"{row[3]}-cvr{row[4]}"
        )
        _require(case.get("case_id") == expected_id, "Q043 case ID or axis label drifted")
        observed_ids.append(expected_id)
    _require(
        observed_rows == _expected_case_axis_rows()
        and len(observed_ids) == len(set(observed_ids)),
        "Q043 registered admission matrix axes, order, or uniqueness drifted",
    )


def _case_map() -> dict[str, dict[str, object]]:
    cases = tuple(dict(case) for case in oracle.expected_cases())
    _validate_exact_case_matrix(cases)
    return {str(case["case_id"]): case for case in cases}


def _validate_deck_manifest_matrix(manifest: Mapping[str, object]) -> None:
    expected_cases = list(_case_map())
    expected_decompositions = list(
        dict.fromkeys(
            decomposition
            for dimension in REQUIRED_DIMENSIONS
            for decomposition in REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension]
        )
    )
    expected_axes = {
        "dimensions": list(REQUIRED_DIMENSIONS),
        "resolutions": list(REQUIRED_RESOLUTIONS),
        "ppc": list(REQUIRED_PPC_VALUES),
        "decompositions": expected_decompositions,
        "decompositions_by_dimension": {
            str(dimension): list(REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension])
            for dimension in REQUIRED_DIMENSIONS
        },
        "artificial_c_over_v_cr": list(REQUIRED_ARTIFICIAL_C_OVER_V_CR_VALUES),
    }
    records = manifest.get("cases")
    _require(
        manifest.get("case_count") == EXPECTED_CASE_COUNT
        and type(records) is list
        and len(records) == EXPECTED_CASE_COUNT
        and all(type(record) is dict for record in records)
        and [record.get("case_id") for record in records] == expected_cases
        and _strict_equal(manifest.get("matrix_axes"), expected_axes),
        "Q043 checked-in deck manifest is not the exact 132-case matrix",
    )


def required_candidate_source_paths() -> frozenset[str]:
    """Return the exact Q043 source/archive closure shared by every case."""
    return frozenset(
        {
            *REQUIRED_CANDIDATE_SOURCE_PATHS,
            *(
                "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle/"
                f"{case_id}.athinput"
                for case_id in _case_map()
            ),
        }
    )


def _case_contract(case_id: object) -> dict[str, object]:
    case_id = _text(case_id, label="case_id")
    cases = _case_map()
    _require(case_id in cases, "case_id: unknown Q043 matrix case")
    return cases[case_id]


def _deck_binding(case: Mapping[str, object], source_root: Path) -> dict[str, object]:
    relative = (
        "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle/"
        f"{case['case_id']}.athinput"
    )
    path = _source_member(source_root, relative, label=f"{case['case_id']}/deck")
    payload, _ = _read_stable_regular_file(path, label=f"{case['case_id']}/deck")
    expected = oracle.render_oracle_deck(case).encode("utf-8")
    _require(payload == expected, f"{case['case_id']}: checked-in deck bytes drifted")
    manifest_path = _source_member(
        source_root,
        "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle/"
        "deck_manifest.json",
        label="Q043 deck manifest",
    )
    manifest_payload, _ = _read_stable_regular_file(
        manifest_path, label="Q043 deck manifest"
    )
    manifest = _json_payload(manifest_payload, label="Q043 deck manifest")
    _validate_deck_manifest_matrix(manifest)
    records = [
        item
        for item in manifest.get("cases", [])
        if type(item) is dict and item.get("case_id") == case["case_id"]
    ]
    digest = hashlib.sha256(payload).hexdigest()
    _require(
        len(records) == 1
        and records[0].get("deck_path") == f"{case['case_id']}.athinput"
        and records[0].get("deck_sha256") == digest,
        f"{case['case_id']}: checked-in deck manifest binding drifted",
    )
    return {
        "path": relative,
        "absolute_path": str(path),
        "sha256": digest,
        "byte_count": len(payload),
        "deck_manifest": {
            "path": str(manifest_path.relative_to(source_root)),
            "sha256": hashlib.sha256(manifest_payload).hexdigest(),
            "byte_count": len(manifest_payload),
        },
    }


def _source_archive_closure(
    payload: bytes,
    *,
    source_root: Path,
    required_paths: frozenset[str],
) -> dict[str, object]:
    try:
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:*") as archive:
            members = archive.getmembers()
            observed: dict[str, tarfile.TarInfo] = {}
            for member in members:
                name = member.name.rstrip("/")
                _require(bool(name), "clean-candidate source archive has empty member")
                name = _relative_path(name, label="clean-candidate source archive member")
                _require(
                    name not in observed,
                    "clean-candidate source archive has duplicate member",
                )
                _require(
                    member.isdir() or member.isreg(),
                    "clean-candidate source archive has non-regular member",
                )
                observed[name] = member
            missing = sorted(required_paths - set(observed))
            _require(
                not missing,
                f"clean-candidate source archive closure is incomplete: {missing}",
            )
            bindings = []
            for relative in sorted(required_paths):
                member = observed[relative]
                _require(
                    member.isreg(),
                    f"clean-candidate source archive member is not regular: {relative}",
                )
                stream = archive.extractfile(member)
                _require(
                    stream is not None,
                    f"clean-candidate source archive member is unreadable: {relative}",
                )
                archived = stream.read()
                source_path = _source_member(
                    source_root, relative, label=f"candidate archive closure/{relative}"
                )
                reviewed, _ = _read_stable_regular_file(
                    source_path, label=f"reviewed source/{relative}"
                )
                _require(
                    archived == reviewed,
                    f"clean-candidate source archive member drifted: {relative}",
                )
                bindings.append(
                    {
                        "path": relative,
                        "sha256": hashlib.sha256(archived).hexdigest(),
                        "byte_count": len(archived),
                    }
                )
    except (tarfile.TarError, OSError) as error:
        raise AdmissionError("clean-candidate source archive is malformed") from error
    return {"members": bindings, "sha256": canonical_sha256(bindings)}


def _candidate_binding(
    value: object,
    authorized_root: Path,
    *,
    source_root: Path,
    required_paths: frozenset[str],
) -> dict[str, object]:
    record = _object(
        value,
        {
            "git_commit",
            "source_bundle_sha256",
            "clean_candidate_manifest",
            "source_archive",
            "executable",
            "environment_profile",
        },
        label="candidate_binding",
    )
    git_commit = _git_commit(record["git_commit"], label="candidate_binding/git_commit")
    bundle = _sha256(
        record["source_bundle_sha256"],
        label="candidate_binding/source_bundle_sha256",
    )
    manifest_binding, manifest_payload = _verified_binding(
        record["clean_candidate_manifest"],
        root=authorized_root,
        label="clean-candidate manifest",
    )
    archive_binding, archive_payload = _verified_binding(
        record["source_archive"],
        root=authorized_root,
        label="clean-candidate source archive",
    )
    executable_binding, _ = _verified_binding(
        record["executable"],
        root=authorized_root,
        label="clean-candidate executable",
        executable=True,
    )
    environment_binding, _ = _verified_binding(
        record["environment_profile"],
        root=authorized_root,
        label="registered environment profile",
    )
    manifest = _json_payload(manifest_payload, label="clean-candidate manifest")
    source = manifest.get("source")
    build = manifest.get("build")
    _require(
        manifest.get("schema_version") == 4
        and type(source) is dict
        and type(build) is dict,
        "clean-candidate manifest identity drifted",
    )
    _require(
        source.get("worktree_status") == "clean"
        and source.get("git_commit") == git_commit
        and source.get("source_bundle_sha256") == bundle
        and source.get("archive_path") == archive_binding["path"]
        and source.get("archive_sha256") == archive_binding["sha256"]
        and build.get("source_archive_sha256") == archive_binding["sha256"]
        and build.get("source_bundle_sha256") == bundle
        and build.get("executable_path") == executable_binding["path"]
        and build.get("executable_sha256") == executable_binding["sha256"],
        "clean-candidate manifest candidate cross-link drifted",
    )
    archive_closure = _source_archive_closure(
        archive_payload, source_root=source_root, required_paths=required_paths
    )
    return {
        "git_commit": git_commit,
        "source_bundle_sha256": bundle,
        "clean_candidate_manifest": manifest_binding,
        "source_archive": archive_binding,
        "source_archive_closure": archive_closure,
        "executable": executable_binding,
        "environment_profile": environment_binding,
    }


def _meshblock_grid(case: Mapping[str, object]) -> tuple[int, int, int]:
    _validate_case_decomposition_contract(case)
    values = []
    for global_count, block_count in zip(case["global_nx"], case["meshblock_nx"]):
        _require(
            type(global_count) is int
            and type(block_count) is int
            and global_count > 0
            and block_count > 0
            and global_count % block_count == 0,
            f"{case['case_id']}: invalid MeshBlock decomposition",
        )
        values.append(global_count // block_count)
    grid = tuple(values)
    _require(
        grid == tuple(case["meshblock_grid"])
        and grid == REQUIRED_DECOMPOSITION_GRIDS[str(case["decomposition"])]
        and math.prod(grid) == int(case["mpi_ranks"]),
        f"{case['case_id']}: exact one-MeshBlock-per-rank contract drifted",
    )
    return grid


def _partition_axes(grid: Sequence[int]) -> list[str]:
    return [f"x{index}" for index, count in enumerate(grid, 1) if count > 1]


def expected_mpi_evidence(case: Mapping[str, object]) -> dict[str, object]:
    """Return exact registered MPI and decomposition evidence for one case."""
    ranks = int(case["mpi_ranks"])
    grid = _meshblock_grid(case)
    return {
        "launcher": "srun",
        "mpi_enabled": True,
        "requested_nodes": 1,
        "requested_tasks": ranks,
        "observed_world_size": ranks,
        "observed_rank_ids": list(range(ranks)),
        "meshblock_grid": list(grid),
        "partition_axes": _partition_axes(grid),
        "rank_meshblock_ids": [
            {"rank": rank, "meshblock_ids": [rank]} for rank in range(ranks)
        ],
    }


def _expected_command(
    case: Mapping[str, object],
    *,
    candidate: Mapping[str, object],
    launch_deck: Mapping[str, object],
    raw_root: Path,
) -> list[str]:
    ranks = int(case["mpi_ranks"])
    return [
        "srun",
        "--nodes=1",
        f"--ntasks={ranks}",
        f"--ntasks-per-node={ranks}",
        "--cpus-per-task=1",
        "--gpus-per-task=1",
        "--gpu-bind=closest",
        str(candidate["executable"]["path"]),
        "-i",
        str(launch_deck["path"]),
        "-d",
        str(raw_root),
    ]


def _launch_contract(
    *,
    case: Mapping[str, object],
    candidate: Mapping[str, object],
    deck: Mapping[str, object],
    launch_deck: Mapping[str, object],
    case_root: Path,
    raw_root: Path,
    artifact_dir: Path,
) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": LAUNCH_CONTRACT_RECORD_TYPE,
        "registration_scope": "registered_science",
        "campaign_id": oracle.CAMPAIGN_ID,
        "case_id": case["case_id"],
        "candidate_binding_sha256": canonical_sha256(candidate),
        "clean_candidate_manifest": _public_binding(
            candidate["clean_candidate_manifest"]
        ),
        "source_archive": _public_binding(candidate["source_archive"]),
        "executable": _public_binding(candidate["executable"]),
        "environment_profile": _public_binding(candidate["environment_profile"]),
        "deck": {
            "reviewed_checked_in_deck": {
                key: deck[key]
                for key in ("path", "absolute_path", "sha256", "byte_count")
            },
            "immutable_launch_deck": _public_binding(launch_deck),
        },
        "authorized_orion_case_root": str(case_root),
        "raw_output_root": str(raw_root),
        "artifact_dir": str(artifact_dir),
        "command": _expected_command(
            case, candidate=candidate, launch_deck=launch_deck, raw_root=raw_root
        ),
        "mpi_evidence": expected_mpi_evidence(case),
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "policy_mutation_authorized": False,
    }


def _validate_stdout(
    payload: bytes, *, case: Mapping[str, object]
) -> dict[str, object]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise AdmissionError("stdout: invalid UTF-8") from error
    _require(
        "FATAL ERROR" not in text and "Traceback" not in text,
        "stdout: fatal execution evidence present",
    )
    ranks = int(case["mpi_ranks"])
    rank_matches = _RANK_EVIDENCE.findall(text)
    _require(
        rank_matches
        == [
            (
                case["case_id"],
                str(ranks),
                ",".join(str(rank) for rank in range(ranks)),
            )
        ],
        "stdout: registered rank evidence missing or ambiguous",
    )
    times = _TIME_CYCLE.findall(text)
    limits = _LIMITS.findall(text)
    exits = _EXIT_EVIDENCE.findall(text)
    terminations = [
        line for line in text.splitlines() if line.startswith("Terminating on ")
    ]
    _require(
        len(times) == 1
        and len(limits) == 1
        and exits == [("0", "0")]
        and terminations == ["Terminating on cycle limit"],
        "stdout: terminal success evidence missing or ambiguous",
    )
    try:
        observed_time = float(times[0][0])
        tlim = float(limits[0][0])
        cycle = int(times[0][1])
        nlim = int(limits[0][1])
    except ValueError as error:
        raise AdmissionError("stdout: terminal evidence is non-numeric") from error
    _require(
        math.isfinite(observed_time)
        and observed_time > 0.0
        and math.isfinite(tlim)
        and cycle == 1
        and nlim == 1,
        "stdout: terminal cycle-limit evidence drifted",
    )
    return {
        "termination_reason": "Terminating on cycle limit",
        "terminal_cycle": cycle,
        "observed_time": observed_time,
        "tlim": tlim,
        "nlim": nlim,
        "observed_world_size": ranks,
        "observed_rank_ids": list(range(ranks)),
        "exit_code": 0,
        "signal": 0,
    }


def _execution_binding(
    value: object,
    *,
    case: Mapping[str, object],
    candidate: Mapping[str, object],
    deck: Mapping[str, object],
    authorized_root: Path,
) -> dict[str, object]:
    record = _object(
        value,
        {
            "case_root",
            "raw_output_root",
            "artifact_dir",
            "launch_contract",
            "launch_deck",
            "registered_execution_receipt",
            "stdout_artifact",
            "terminal_receipt",
        },
        label="execution_binding",
    )
    case_root = _canonical_below(
        record["case_root"], authorized_root, label="execution_binding/case_root", directory=True
    )
    raw_root = _canonical_below(
        record["raw_output_root"],
        authorized_root,
        label="execution_binding/raw_output_root",
        directory=True,
    )
    artifact_dir = _canonical_below(
        record["artifact_dir"],
        authorized_root,
        label="execution_binding/artifact_dir",
        directory=True,
    )
    launch_deck_binding, launch_deck_payload = _verified_binding(
        record["launch_deck"], root=artifact_dir, label="immutable launch deck"
    )
    contract_binding, contract_payload = _verified_binding(
        record["launch_contract"], root=artifact_dir, label="launch contract"
    )
    receipt_binding, receipt_payload = _verified_binding(
        record["registered_execution_receipt"],
        root=artifact_dir,
        label="registered execution receipt",
    )
    stdout_binding, stdout_payload = _verified_binding(
        record["stdout_artifact"], root=artifact_dir, label="registered stdout"
    )
    terminal_binding, terminal_payload = _verified_binding(
        record["terminal_receipt"], root=artifact_dir, label="terminal receipt"
    )
    _require(
        launch_deck_binding["path"]
        == str(artifact_dir / "registration/input_deck.athinput")
        and contract_binding["path"]
        == str(artifact_dir / "registration/launch_contract.json")
        and receipt_binding["path"]
        == str(artifact_dir / "registration/registered_execution_receipt.json")
        and stdout_binding["path"] == str(artifact_dir / "stdout/athena_stdout.txt")
        and terminal_binding["path"] == str(artifact_dir / "terminal/terminal_receipt.json"),
        "execution_binding: evidence path layout drifted",
    )
    _require(
        launch_deck_binding["sha256"] == deck["sha256"]
        and launch_deck_binding["byte_count"] == deck["byte_count"]
        and hashlib.sha256(launch_deck_payload).hexdigest() == deck["sha256"],
        "execution_binding: immutable launch deck differs from checked-in deck",
    )
    contract = _json_payload(contract_payload, label="launch contract")
    receipt = _json_payload(receipt_payload, label="registered execution receipt")
    terminal = _json_payload(terminal_payload, label="terminal receipt")
    submission_id = _text(
        receipt.get("submission_id"), label="registered execution receipt/submission_id"
    )
    _require(
        _UUID.fullmatch(submission_id) is not None,
        "registered execution receipt/submission_id: malformed UUID",
    )
    expected_case_root = (
        authorized_root / RUN_NAMESPACE / str(case["case_id"]) / submission_id
    )
    _require(
        case_root == expected_case_root
        and raw_root == case_root / "raw"
        and artifact_dir == case_root / "artifacts",
        "execution_binding: authorized Orion case-root layout drifted",
    )
    expected_contract = _launch_contract(
        case=case,
        candidate=candidate,
        deck=deck,
        launch_deck=launch_deck_binding,
        case_root=case_root,
        raw_root=raw_root,
        artifact_dir=artifact_dir,
    )
    _require(
        _strict_equal(contract, expected_contract),
        "launch contract: candidate, deck, command, MPI, or authority drifted",
    )
    receipt = _object(
        receipt,
        {
            "schema_version",
            "record_type",
            "receipt_role",
            "registration_scope",
            "reconciled",
            "campaign_id",
            "case_id",
            "reservation_id",
            "submission_id",
            "reconciliation_event_sha256",
            "source_commit",
            "source_bundle_sha256",
            "source_archive_sha256",
            "clean_candidate_manifest_sha256",
            "executable_sha256",
            "environment_sha256",
            "deck_sha256",
            "launch_contract_sha256",
            "command",
            "mpi_evidence",
            "slurm_job_id",
            "slurm_terminal_state",
            "slurm_exit_code",
            "raw_output_root",
            "artifact_dir",
            "stdout_sha256",
            "terminal_receipt_sha256",
            "pre_submit_manifest_sha256",
        },
        label="registered execution receipt",
    )
    _require(
        receipt["schema_version"] == SCHEMA_VERSION
        and receipt["record_type"] == EXECUTION_RECEIPT_RECORD_TYPE
        and receipt["receipt_role"] == "immutable_reconciled_registered_execution"
        and receipt["registration_scope"] == "registered_science"
        and receipt["reconciled"] is True
        and receipt["campaign_id"] == oracle.CAMPAIGN_ID
        and receipt["case_id"] == case["case_id"]
        and type(receipt["reservation_id"]) is str
        and _UUID.fullmatch(receipt["reservation_id"]) is not None
        and receipt["submission_id"] == submission_id
        and type(receipt["reconciliation_event_sha256"]) is str
        and _SHA256.fullmatch(receipt["reconciliation_event_sha256"]) is not None
        and receipt["source_commit"] == candidate["git_commit"]
        and receipt["source_bundle_sha256"] == candidate["source_bundle_sha256"]
        and receipt["source_archive_sha256"] == candidate["source_archive"]["sha256"]
        and receipt["clean_candidate_manifest_sha256"]
        == candidate["clean_candidate_manifest"]["sha256"]
        and receipt["executable_sha256"] == candidate["executable"]["sha256"]
        and receipt["environment_sha256"] == candidate["environment_profile"]["sha256"]
        and receipt["deck_sha256"] == deck["sha256"]
        and receipt["launch_contract_sha256"]
        == hashlib.sha256(contract_payload).hexdigest()
        and _strict_equal(receipt["command"], expected_contract["command"])
        and _strict_equal(receipt["mpi_evidence"], expected_contract["mpi_evidence"])
        and type(receipt["slurm_job_id"]) is str
        and _JOB_ID.fullmatch(receipt["slurm_job_id"]) is not None
        and receipt["slurm_terminal_state"] == "COMPLETED"
        and receipt["slurm_exit_code"] == "0:0"
        and receipt["raw_output_root"] == str(raw_root)
        and receipt["artifact_dir"] == str(artifact_dir)
        and receipt["stdout_sha256"] == stdout_binding["sha256"]
        and receipt["terminal_receipt_sha256"] == terminal_binding["sha256"]
        and type(receipt["pre_submit_manifest_sha256"]) is str
        and _SHA256.fullmatch(receipt["pre_submit_manifest_sha256"]) is not None,
        "registered execution receipt immutable cross-link drifted",
    )
    stdout_terminal = _validate_stdout(stdout_payload, case=case)
    expected_terminal = {
        "schema_version": SCHEMA_VERSION,
        "record_type": TERMINAL_RECEIPT_RECORD_TYPE,
        "campaign_id": oracle.CAMPAIGN_ID,
        "case_id": case["case_id"],
        "submission_id": submission_id,
        "slurm_job_id": receipt["slurm_job_id"],
        "slurm_terminal_state": "COMPLETED",
        "slurm_exit_code": "0:0",
        "termination_reason": "cycle_limit",
        "terminal_cycle": 1,
        "observed_world_size": int(case["mpi_ranks"]),
        "stdout_sha256": stdout_binding["sha256"],
    }
    _require(
        _strict_equal(terminal, expected_terminal),
        "terminal receipt: completion evidence drifted",
    )
    return {
        "case_root": str(case_root),
        "raw_output_root": str(raw_root),
        "artifact_dir": str(artifact_dir),
        "launch_deck": launch_deck_binding,
        "launch_contract": contract_binding,
        "registered_execution_receipt": receipt_binding,
        "stdout_artifact": stdout_binding,
        "terminal_receipt": terminal_binding,
        "command": expected_contract["command"],
        "mpi_evidence": expected_contract["mpi_evidence"],
        "registered_execution_identity": {
            key: receipt[key]
            for key in (
                "reservation_id",
                "submission_id",
                "reconciliation_event_sha256",
                "slurm_job_id",
                "pre_submit_manifest_sha256",
            )
        },
        "stdout_terminal_success": stdout_terminal,
    }


def _expected_runtime_parameters(
    case: Mapping[str, object],
) -> dict[str, dict[str, str]]:
    helper = getattr(oracle, "expected_runtime_parameters", None)
    _require(
        callable(helper),
        "hardened Q043 exact frozen-runtime helper is unavailable",
    )
    try:
        expected = helper(case)
    except oracle.ContractError as error:
        raise AdmissionError("hardened Q043 frozen-runtime contract failed") from error
    _require(type(expected) is dict, "hardened Q043 frozen-runtime contract is malformed")
    return expected


def _normalized_runtime_parameters(
    parameters: Mapping[str, Mapping[str, str]],
    *,
    case: Mapping[str, object],
    field: str,
    cycle: int,
) -> dict[str, dict[str, str]]:
    _require(field in REQUIRED_FIELDS, "runtime output field is unknown")
    _require(cycle in REQUIRED_CYCLES, "runtime output cycle is unknown")
    output_blocks = {block for block in parameters if block.startswith("output")}
    _require(
        output_blocks == set(OUTPUT_BLOCKS),
        f"{field}/cycle{cycle}: runtime output block inventory drifted",
    )
    field_index = REQUIRED_FIELDS.index(field) + 1
    normalized = {}
    for block, values in parameters.items():
        _require(
            type(block) is str and type(values) is dict,
            f"{field}/cycle{cycle}: runtime parameter block is malformed",
        )
        numbered_output = block in OUTPUT_BLOCKS
        if numbered_output:
            output_index = int(block[6:])
            _require(
                values.get("file_type") == "bin"
                and values.get("variable") == REQUIRED_FIELDS[output_index - 1]
                and values.get("id") == REQUIRED_FIELDS[output_index - 1]
                and values.get("dcycle") == str(oracle.RAW_ORACLE_DCYCLE)
                and values.get("ghost_zones") == "false"
                and values.get("single_file_per_rank")
                == ("true" if int(case["mpi_ranks"]) > 1 else "false"),
                f"{field}/cycle{cycle}: runtime {block} immutable contract drifted",
            )
            _require(
                OUTPUT_BOOKKEEPING_KEYS <= set(values),
                f"{field}/cycle{cycle}: runtime {block} bookkeeping is incomplete",
            )
            file_number = values["file_number"]
            try:
                parsed_file_number = int(file_number)
            except (TypeError, ValueError) as error:
                raise AdmissionError(
                    f"{field}/cycle{cycle}: runtime {block} file_number is invalid"
                ) from error
            _require(
                type(file_number) is str and file_number == str(parsed_file_number),
                f"{field}/cycle{cycle}: runtime {block} file_number is noncanonical",
            )
            expected_file_number = cycle if output_index <= field_index else cycle + 1
            _require(
                parsed_file_number == expected_file_number,
                f"{field}/cycle{cycle}: runtime {block} sequential file_number drifted",
            )
            try:
                last_time = float(values["last_time"])
            except (TypeError, ValueError) as error:
                raise AdmissionError(
                    f"{field}/cycle{cycle}: runtime {block} last_time is invalid"
                ) from error
            expected_last_time = -1.0 if cycle == 0 and output_index <= field_index else 0.0
            _require(
                math.isfinite(last_time) and last_time == expected_last_time,
                f"{field}/cycle{cycle}: runtime {block} sequential last_time drifted",
            )
        normalized[block] = {
            key: value
            for key, value in values.items()
            if not (numbered_output and key in OUTPUT_BOOKKEEPING_KEYS)
        }
    _require(
        _strict_equal(normalized, _expected_runtime_parameters(case)),
        f"{field}/cycle{cycle}: runtime parameters drifted from exact frozen deck",
    )
    return normalized


def _validate_raw_dataset(
    dataset: binary.AthenaBinaryDataset,
    *,
    case: Mapping[str, object],
    field: str,
    cycle: int,
) -> None:
    _require(dataset.cycle == cycle, f"{field}/cycle{cycle}: cycle drifted")
    _require(
        (cycle == 0 and dataset.time == 0.0) or (cycle == 1 and dataset.time > 0.0),
        f"{field}/cycle{cycle}: observed time drifted",
    )
    _require(
        dataset.variable_names == (field,),
        f"{field}/cycle{cycle}: raw variable schema drifted",
    )
    _require(
        dataset.root_grid_shape == tuple(case["global_nx"])
        and dataset.meshblock_shape == tuple(case["meshblock_nx"]),
        f"{field}/cycle{cycle}: raw geometry drifted",
    )
    _normalized_runtime_parameters(
        dataset.input_parameters, case=case, field=field, cycle=cycle
    )


def _raw_relative_path(
    case: Mapping[str, object], *, field: str, cycle: int, rank: int
) -> str:
    basename = str(case["case_id"]).replace("-", "_")
    filename = f"{basename}.{field}.{cycle:05d}.bin"
    if int(case["mpi_ranks"]) > 1:
        return f"bin/rank_{rank:08d}/{filename}"
    return f"bin/{filename}"


def _raw_artifact(
    value: object,
    *,
    case: Mapping[str, object],
    raw_root: Path,
    index: int,
) -> tuple[dict[str, object], Path]:
    label = f"raw_artifacts[{index}]"
    record = _object(
        value,
        {"path", "sha256", "byte_count", "case_id", "field", "cycle", "rank"},
        label=label,
    )
    field = _text(record["field"], label=f"{label}/field")
    _require(field in REQUIRED_FIELDS, f"{label}: unknown raw field")
    cycle = _nonnegative_integer(record["cycle"], label=f"{label}/cycle")
    rank = _nonnegative_integer(record["rank"], label=f"{label}/rank")
    _require(
        record["case_id"] == case["case_id"]
        and cycle in REQUIRED_CYCLES
        and rank < int(case["mpi_ranks"]),
        f"{label}: mixed or invalid case/cycle/rank lineage",
    )
    relative = _relative_path(record["path"], label=f"{label}/path")
    _require(
        relative == _raw_relative_path(case, field=field, cycle=cycle, rank=rank),
        f"{label}: raw path does not match exact output inventory",
    )
    absolute = _source_member(raw_root, relative, label=label)
    digest = _sha256(record["sha256"], label=f"{label}/sha256")
    byte_count = _positive_integer(record["byte_count"], label=f"{label}/byte_count")
    payload, identity = _read_stable_regular_file(
        absolute,
        label=label,
        expected_sha256=digest,
        expected_byte_count=byte_count,
    )
    try:
        dataset = binary.read_athenak_binary(absolute)
    except binary.AnalysisError as error:
        raise AdmissionError(f"{label}: malformed AthenaK binary output") from error
    _validate_raw_dataset(dataset, case=case, field=field, cycle=cycle)
    return (
        {
            "path": relative,
            "absolute_path": str(absolute),
            "sha256": digest,
            "byte_count": len(payload),
            "case_id": case["case_id"],
            "field": field,
            "cycle": cycle,
            "rank": rank,
            "filesystem_identity": identity,
        },
        absolute,
    )


def _raw_filesystem_inventory(raw_root: Path) -> list[str]:
    files = []
    for path in raw_root.rglob("*"):
        relative = path.relative_to(raw_root).as_posix()
        try:
            metadata = path.lstat()
        except OSError as error:
            raise AdmissionError("raw output tree changed during inventory") from error
        _require(
            not stat.S_ISLNK(metadata.st_mode),
            f"raw output tree contains symlink: {relative}",
        )
        if stat.S_ISDIR(metadata.st_mode):
            continue
        _require(
            stat.S_ISREG(metadata.st_mode),
            f"raw output tree contains non-regular artifact: {relative}",
        )
        files.append(_relative_path(relative, label="raw output tree member"))
    return sorted(files)


def _authorization_boundary() -> dict[str, object]:
    return {
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "policy_mutation_authorized": False,
        "q043_qualification_authorized": False,
        "q023_qualification_authorized": False,
        "q019_qualification_authorized": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "required_next_boundary": (
            "separate_review_of_complete_registered_multidirectional_Q043_matrix"
        ),
    }


def _source_local_insufficiency() -> dict[str, object]:
    return {
        "source_local_or_synthetic_raw_analysis_sufficient": False,
        "source_local_oracle_result_accepted_as_registered_evidence": False,
        "source_local_foundation_role": (
            "analysis_dependency_only_not_registered_execution_evidence"
        ),
        "registered_case_admission_required": True,
        "complete_registered_matrix_required": True,
        "exact_132_case_dimension_valid_matrix_required": True,
        "registered_raw_execution_required_for_every_exact_matrix_case": True,
        "multidirectional_mpi_coverage_required": True,
        "downstream_q023_or_q019_qualification_effect": "none",
    }


def _validate_source_local_case_result(
    value: object, *, case: Mapping[str, object]
) -> dict[str, object]:
    result = _object(
        value, SOURCE_LOCAL_CASE_RESULT_KEYS, label="source-local case oracle result"
    )
    provenance = _array(
        result["raw_provenance"], label="source-local case oracle result/raw_provenance"
    )
    observed_shards = []
    for index, item in enumerate(provenance):
        binding = _object(
            item,
            SOURCE_LOCAL_RAW_BINDING_KEYS,
            label=f"source-local case oracle result/raw_provenance[{index}]",
        )
        field = _text(
            binding["field"],
            label=f"source-local case oracle result/raw_provenance[{index}]/field",
        )
        shard = _nonnegative_integer(
            binding["shard_index"],
            label=f"source-local case oracle result/raw_provenance[{index}]/shard_index",
        )
        _absolute_path(
            binding["path"],
            label=f"source-local case oracle result/raw_provenance[{index}]/path",
        )
        _positive_integer(
            binding["size"],
            label=f"source-local case oracle result/raw_provenance[{index}]/size",
        )
        _sha256(
            binding["sha256"],
            label=f"source-local case oracle result/raw_provenance[{index}]/sha256",
        )
        observed_shards.append((field, shard))
    expected_shards = [
        (field, shard)
        for field in REQUIRED_FIELDS
        for shard in range(int(case["mpi_ranks"]))
    ]
    _require(
        result["schema_version"] == oracle.SCHEMA_VERSION
        and result["record_type"] == SOURCE_LOCAL_CASE_RECORD_TYPE
        and result["campaign_id"] == oracle.CAMPAIGN_ID
        and result["case_id"] == case["case_id"]
        and result["qualification_effect"] == oracle.QUALIFICATION_EFFECT
        and result["launch_authorized"] is False
        and result["scientific_claim_authorized"] is False
        and result["publication_authorized"] is False
        and result["passed"] is False
        and result["source_local_oracle_check_pass"] is True
        and _strict_equal(result["case_contract"], case)
        and result["cross_field_runtime_metadata"] == oracle.RUNTIME_METADATA_CONTRACT
        and result["tolerance_scope"]
        == (
            "source_local_exact_deposition_oracle_representation_roundoff_only_"
            "not_science_acceptance"
        )
        and observed_shards == expected_shards,
        f"{case['case_id']}: source-local case oracle schema or authority drifted",
    )
    return result


def _validate_source_local_matrix_result(value: object) -> dict[str, object]:
    result = _object(
        value, SOURCE_LOCAL_MATRIX_RESULT_KEYS, label="source-local matrix oracle result"
    )
    expected_cases = list(_case_map().values())
    nested = _array(
        result["cases"], label="source-local matrix oracle result/cases"
    )
    _require(
        len(nested) == EXPECTED_CASE_COUNT,
        "source-local matrix oracle result/cases: exact case count drifted",
    )
    for case, item in zip(expected_cases, nested):
        _validate_source_local_case_result(item, case=case)
    expected_decompositions = list(
        dict.fromkeys(
            decomposition
            for dimension in REQUIRED_DIMENSIONS
            for decomposition in REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension]
        )
    )
    _require(
        result["schema_version"] == oracle.SCHEMA_VERSION
        and result["record_type"] == SOURCE_LOCAL_MATRIX_RECORD_TYPE
        and result["campaign_id"] == oracle.CAMPAIGN_ID
        and result["qualification_effect"] == oracle.QUALIFICATION_EFFECT
        and result["launch_authorized"] is False
        and result["scientific_claim_authorized"] is False
        and result["publication_authorized"] is False
        and result["passed"] is False
        and result["source_local_oracle_check_pass"] is True
        and result["case_count"] == EXPECTED_CASE_COUNT
        and _strict_equal(result["dimensions_verified"], list(REQUIRED_DIMENSIONS))
        and _strict_equal(result["resolutions_verified"], list(REQUIRED_RESOLUTIONS))
        and _strict_equal(result["ppc_verified"], list(REQUIRED_PPC_VALUES))
        and _strict_equal(result["decompositions_verified"], expected_decompositions)
        and _strict_equal(
            result["artificial_c_over_v_cr_verified"],
            list(REQUIRED_ARTIFICIAL_C_OVER_V_CR_VALUES),
        )
        and result["required_output_cycle"] == 1
        and result["output_dcycle"] == oracle.RAW_ORACLE_DCYCLE
        and result["output_timing"]
        == "cycle_zero_initialization_and_cycle_one_finalize_only"
        and result["cross_field_runtime_metadata"] == oracle.RUNTIME_METADATA_CONTRACT
        and result["initial_state_verified"]
        == "uniform_zero_perturbation_parallel_stream"
        and result["artificial_c_in_formula"] is False
        and result["configured_current_formula"]
        == "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=2*B_g*k0",
        "source-local matrix oracle schema or authority drifted",
    )
    return result


def build_case_admission(
    *,
    case_id: str,
    candidate_binding: Mapping[str, object],
    execution_binding: Mapping[str, object],
    raw_artifacts: Sequence[Mapping[str, object]],
    source_root: str | Path = REPO_ROOT,
    authorized_orion_root: str | Path = AUTHORIZED_ORION_ROOT,
) -> dict[str, object]:
    """Build one filesystem-backed non-authorizing registered Q043 admission."""
    source = _canonical_root(source_root, label="source_root")
    authorized = _canonical_root(authorized_orion_root, label="authorized_orion_root")
    case = _case_contract(case_id)
    deck = _deck_binding(case, source)
    candidate = _candidate_binding(
        candidate_binding,
        authorized,
        source_root=source,
        required_paths=required_candidate_source_paths(),
    )
    execution = _execution_binding(
        execution_binding,
        case=case,
        candidate=candidate,
        deck=deck,
        authorized_root=authorized,
    )
    _require(type(raw_artifacts) is list, "raw_artifacts: expected array")
    verified = [
        _raw_artifact(
            item,
            case=case,
            raw_root=Path(execution["raw_output_root"]),
            index=index,
        )
        for index, item in enumerate(raw_artifacts)
    ]
    artifacts = [item[0] for item in verified]
    paths = [str(item["absolute_path"]) for item in artifacts]
    identities = [
        (item["filesystem_identity"]["device"], item["filesystem_identity"]["inode"])
        for item in artifacts
    ]
    _require(len(paths) == len(set(paths)), "raw_artifacts: path reused")
    _require(
        len(identities) == len(set(identities)),
        "raw_artifacts: filesystem object reused",
    )
    expected_inventory = [
        (
            cycle,
            field,
            rank,
            _raw_relative_path(case, field=field, cycle=cycle, rank=rank),
        )
        for cycle in REQUIRED_CYCLES
        for field in REQUIRED_FIELDS
        for rank in range(int(case["mpi_ranks"]))
    ]
    observed_inventory = [
        (item["cycle"], item["field"], item["rank"], item["path"]) for item in artifacts
    ]
    _require(
        observed_inventory == expected_inventory,
        "raw_artifacts: expected exactly one cycle-zero and cycle-one field/rank artifact",
    )
    _require(
        _raw_filesystem_inventory(Path(execution["raw_output_root"]))
        == sorted(item["path"] for item in artifacts),
        "raw output filesystem inventory differs from the exact declared inventory",
    )
    cycle_one_paths = {
        field: tuple(
            Path(item["absolute_path"])
            for item in artifacts
            if item["cycle"] == 1 and item["field"] == field
        )
        for field in REQUIRED_FIELDS
    }
    try:
        source_local_result = oracle.analyze_raw_case(str(case["case_id"]), cycle_one_paths)
    except oracle.ContractError as error:
        raise AdmissionError(f"{case['case_id']}: source-local raw oracle failed") from error
    source_local_result = _validate_source_local_case_result(
        source_local_result, case=case
    )
    expected_provenance = {
        (item["absolute_path"], item["byte_count"], item["sha256"])
        for item in artifacts
        if item["cycle"] == 1
    }
    measured_provenance = {
        (item["path"], item["size"], item["sha256"])
        for item in source_local_result["raw_provenance"]
    }
    _require(
        expected_provenance == measured_provenance
        and source_local_result.get("source_local_oracle_check_pass") is True
        and source_local_result.get("passed") is False
        and source_local_result.get("launch_authorized") is False
        and source_local_result.get("scientific_claim_authorized") is False
        and source_local_result.get("publication_authorized") is False,
        f"{case['case_id']}: source-local oracle result or provenance drifted",
    )
    hardened_result = {
        "registered_execution_raw_oracle_check_pass": True,
        "strict_cycle_zero_and_cycle_one_runtime_metadata_check_pass": True,
        "exact_output_inventory_check_pass": True,
        "source_local_oracle_result": source_local_result,
        "source_local_result_sufficient_for_downstream_qualification": False,
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": CASE_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": oracle.CAMPAIGN_ID,
        "status": "registered_case_evidence_admitted_non_authorizing",
        "qualification_effect": QUALIFICATION_EFFECT,
        "case_id": case["case_id"],
        "case_contract": case,
        "deck_binding": deck,
        "candidate_binding": candidate,
        "candidate_binding_sha256": canonical_sha256(candidate),
        "execution_binding": execution,
        "execution_binding_sha256": canonical_sha256(execution),
        "raw_artifacts": artifacts,
        "raw_inventory_sha256": canonical_sha256(artifacts),
        "hardened_raw_oracle_result": hardened_result,
        "hardened_raw_oracle_result_sha256": canonical_sha256(hardened_result),
        "source_local_insufficiency": _source_local_insufficiency(),
        "authorization": _authorization_boundary(),
    }


def _candidate_input(value: Mapping[str, object]) -> dict[str, object]:
    return {
        "git_commit": value["git_commit"],
        "source_bundle_sha256": value["source_bundle_sha256"],
        **{
            key: _public_binding(value[key])
            for key in (
                "clean_candidate_manifest",
                "source_archive",
                "executable",
                "environment_profile",
            )
        },
    }


def _execution_input(value: Mapping[str, object]) -> dict[str, object]:
    return {
        key: value[key]
        for key in ("case_root", "raw_output_root", "artifact_dir")
    } | {
        key: _public_binding(value[key])
        for key in (
            "launch_deck",
            "launch_contract",
            "registered_execution_receipt",
            "stdout_artifact",
            "terminal_receipt",
        )
    }


def _raw_input(value: Mapping[str, object]) -> dict[str, object]:
    return {
        key: value[key]
        for key in ("path", "sha256", "byte_count", "case_id", "field", "cycle", "rank")
    }


def validate_case_admission(
    value: object,
    *,
    source_root: str | Path = REPO_ROOT,
    authorized_orion_root: str | Path = AUTHORIZED_ORION_ROOT,
) -> dict[str, object]:
    """Rebuild and exactly validate one registered Q043 case admission."""
    record = _object(
        value,
        {
            "schema_version",
            "record_type",
            "successor_id",
            "campaign_id",
            "status",
            "qualification_effect",
            "case_id",
            "case_contract",
            "deck_binding",
            "candidate_binding",
            "candidate_binding_sha256",
            "execution_binding",
            "execution_binding_sha256",
            "raw_artifacts",
            "raw_inventory_sha256",
            "hardened_raw_oracle_result",
            "hardened_raw_oracle_result_sha256",
            "source_local_insufficiency",
            "authorization",
        },
        label="case_admission",
    )
    rebuilt = build_case_admission(
        case_id=record["case_id"],
        candidate_binding=_candidate_input(record["candidate_binding"]),
        execution_binding=_execution_input(record["execution_binding"]),
        raw_artifacts=[_raw_input(item) for item in record["raw_artifacts"]],
        source_root=source_root,
        authorized_orion_root=authorized_orion_root,
    )
    _require(
        _strict_equal(record, rebuilt),
        "case_admission: derived fields or authorization boundary drifted",
    )
    return rebuilt


def multidirectional_mpi_coverage_report(
    case_admissions: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Require the exact registered 132-case dimension-valid MPI matrix."""
    _require(type(case_admissions) is list, "case_admissions: expected array")
    expected_ids = list(_case_map())
    _require(
        len(case_admissions) == EXPECTED_CASE_COUNT,
        "registered Q043 MPI coverage requires the exact 132-case matrix",
    )
    records = []
    partition_records = []
    observed_ids = []
    observed_by_dimension: dict[int, list[str]] = {
        dimension: [] for dimension in REQUIRED_DIMENSIONS
    }
    decomposition_counts: dict[str, int] = {}
    for index, admission in enumerate(case_admissions):
        _require(type(admission) is dict, f"case_admissions[{index}]: expected object")
        case = admission.get("case_contract")
        execution = admission.get("execution_binding")
        _require(
            type(case) is dict and type(execution) is dict,
            f"case_admissions[{index}]: missing case or execution binding",
        )
        _validate_case_decomposition_contract(case)
        evidence = execution.get("mpi_evidence")
        _require(type(evidence) is dict, f"case_admissions[{index}]: missing MPI evidence")
        _require(
            _strict_equal(evidence, expected_mpi_evidence(case)),
            f"case_admissions[{index}]: exact registered MPI evidence drifted",
        )
        case_id = str(case["case_id"])
        dimension = int(case["dimension"])
        decomposition = str(case["decomposition"])
        axes = evidence["partition_axes"]
        observed_ids.append(case_id)
        if decomposition not in observed_by_dimension[dimension]:
            observed_by_dimension[dimension].append(decomposition)
        count_key = f"d{dimension}/{decomposition}"
        decomposition_counts[count_key] = decomposition_counts.get(count_key, 0) + 1
        record = {
            "case_id": case_id,
            "dimension": dimension,
            "decomposition": decomposition,
            "meshblock_grid": evidence["meshblock_grid"],
            "partition_axes": axes,
            "observed_world_size": evidence["observed_world_size"],
        }
        records.append(record)
        if axes:
            partition_records.append(record)
    expected_repetitions = (
        len(REQUIRED_RESOLUTIONS)
        * len(REQUIRED_PPC_VALUES)
        * len(REQUIRED_ARTIFICIAL_C_OVER_V_CR_VALUES)
    )
    expected_counts = {
        f"d{dimension}/{decomposition}": expected_repetitions
        for dimension in REQUIRED_DIMENSIONS
        for decomposition in REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension]
    }
    expected_by_dimension = {
        dimension: list(REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension])
        for dimension in REQUIRED_DIMENSIONS
    }
    _require(
        observed_ids == expected_ids
        and observed_by_dimension == expected_by_dimension
        and decomposition_counts == expected_counts,
        "registered Q043 matrix lacks exact dimension-valid decomposition coverage",
    )
    covered = {
        axis for record in partition_records for axis in record["partition_axes"]
    }
    multi_axis = [
        record for record in partition_records if len(record["partition_axes"]) >= 2
    ]
    return {
        "expected_case_count": EXPECTED_CASE_COUNT,
        "registered_case_count": len(records),
        "required_decompositions_by_dimension": {
            str(dimension): list(REQUIRED_DECOMPOSITIONS_BY_DIMENSION[dimension])
            for dimension in REQUIRED_DIMENSIONS
        },
        "registered_decompositions_by_dimension": {
            str(dimension): observed_by_dimension[dimension]
            for dimension in REQUIRED_DIMENSIONS
        },
        "decomposition_case_counts": decomposition_counts,
        "required_partition_axes": ["x1", "x2", "x3"],
        "covered_partition_axes": sorted(covered),
        "multi_axis_partition_required": True,
        "multi_axis_partition_count": len(multi_axis),
        "registered_mpi_partition_records": partition_records,
    }


def build_matrix_qualification(
    *,
    case_admissions: Sequence[Mapping[str, object]],
    source_root: str | Path = REPO_ROOT,
    authorized_orion_root: str | Path = AUTHORIZED_ORION_ROOT,
) -> dict[str, object]:
    """Build the complete non-authorizing registered Q043 matrix prerequisite."""
    _require(type(case_admissions) is list, "case_admissions: expected array")
    admissions = [
        validate_case_admission(
            item,
            source_root=source_root,
            authorized_orion_root=authorized_orion_root,
        )
        for item in case_admissions
    ]
    expected_ids = list(_case_map())
    observed_ids = [str(item["case_id"]) for item in admissions]
    _require(
        observed_ids == expected_ids,
        "registered Q043 matrix is incomplete, unknown, or noncanonical",
    )
    _require(
        len({item["candidate_binding_sha256"] for item in admissions}) == 1,
        "registered Q043 matrix mixes clean candidates",
    )
    raw_paths = [
        item["absolute_path"]
        for admission in admissions
        for item in admission["raw_artifacts"]
    ]
    raw_identities = [
        (item["filesystem_identity"]["device"], item["filesystem_identity"]["inode"])
        for admission in admissions
        for item in admission["raw_artifacts"]
    ]
    _require(len(raw_paths) == len(set(raw_paths)), "registered Q043 raw path reused")
    _require(
        len(raw_identities) == len(set(raw_identities)),
        "registered Q043 raw filesystem object reused",
    )
    execution_artifacts = [
        admission["execution_binding"][key]
        for admission in admissions
        for key in (
            "launch_deck",
            "launch_contract",
            "registered_execution_receipt",
            "stdout_artifact",
            "terminal_receipt",
        )
    ]
    execution_paths = [item["path"] for item in execution_artifacts]
    execution_filesystem_identities = [
        (item["filesystem_identity"]["device"], item["filesystem_identity"]["inode"])
        for item in execution_artifacts
    ]
    _require(
        len(execution_paths) == len(set(execution_paths)),
        "registered Q043 execution evidence path reused",
    )
    _require(
        len(execution_filesystem_identities)
        == len(set(execution_filesystem_identities)),
        "registered Q043 execution evidence filesystem object reused",
    )
    identities = [
        admission["execution_binding"]["registered_execution_identity"]
        for admission in admissions
    ]
    for key in (
        "reservation_id",
        "submission_id",
        "reconciliation_event_sha256",
        "slurm_job_id",
        "pre_submit_manifest_sha256",
    ):
        values = [item[key] for item in identities]
        _require(
            len(values) == len(set(values)),
            f"registered Q043 execution {key} reused",
        )
    mpi_coverage = multidirectional_mpi_coverage_report(admissions)
    raw_cases = {
        str(admission["case_id"]): {
            field: tuple(
                Path(item["absolute_path"])
                for item in admission["raw_artifacts"]
                if item["cycle"] == 1 and item["field"] == field
            )
            for field in REQUIRED_FIELDS
        }
        for admission in admissions
    }
    try:
        source_local_matrix = oracle.analyze_raw_matrix(raw_cases)
    except oracle.ContractError as error:
        raise AdmissionError("complete source-local Q043 matrix recompute failed") from error
    source_local_matrix = _validate_source_local_matrix_result(source_local_matrix)
    _require(
        source_local_matrix.get("source_local_oracle_check_pass") is True
        and source_local_matrix.get("passed") is False
        and source_local_matrix.get("launch_authorized") is False
        and source_local_matrix.get("scientific_claim_authorized") is False
        and source_local_matrix.get("publication_authorized") is False,
        "source-local matrix result authority or status drifted",
    )
    case_bindings = [
        {
            "case_id": item["case_id"],
            "case_admission_sha256": canonical_sha256(item),
            "candidate_binding_sha256": item["candidate_binding_sha256"],
            "execution_binding_sha256": item["execution_binding_sha256"],
            "raw_inventory_sha256": item["raw_inventory_sha256"],
            "hardened_raw_oracle_result_sha256": item[
                "hardened_raw_oracle_result_sha256"
            ],
        }
        for item in admissions
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": MATRIX_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": oracle.CAMPAIGN_ID,
        "status": "complete_registered_Q043_matrix_prerequisite_non_authorizing",
        "qualification_effect": QUALIFICATION_EFFECT,
        "registered_execution_qualification_check_pass": True,
        "case_count": len(admissions),
        "case_admissions": admissions,
        "case_bindings": case_bindings,
        "case_bindings_sha256": canonical_sha256(case_bindings),
        "multidirectional_mpi_coverage": mpi_coverage,
        "source_local_matrix_result": source_local_matrix,
        "source_local_matrix_result_sufficient_for_downstream_qualification": False,
        "source_local_insufficiency": _source_local_insufficiency(),
        "authorization": _authorization_boundary(),
    }


def validate_matrix_qualification(
    value: object,
    *,
    source_root: str | Path = REPO_ROOT,
    authorized_orion_root: str | Path = AUTHORIZED_ORION_ROOT,
) -> dict[str, object]:
    """Rebuild and exactly validate a complete registered Q043 matrix record."""
    _require(
        type(value) is dict and value.get("record_type") == MATRIX_RECORD_TYPE,
        "downstream Q023/Q019 prerequisite requires registered Q043 matrix record",
    )
    rebuilt = build_matrix_qualification(
        case_admissions=value.get("case_admissions"),
        source_root=source_root,
        authorized_orion_root=authorized_orion_root,
    )
    _require(
        _strict_equal(value, rebuilt),
        "matrix_qualification: derived fields or authorization boundary drifted",
    )
    return rebuilt


def validate_downstream_q023_q019_prerequisite(
    value: object,
    *,
    source_root: str | Path = REPO_ROOT,
    authorized_orion_root: str | Path = AUTHORIZED_ORION_ROOT,
) -> dict[str, object]:
    """Reject source-local evidence and accept only an exact registered matrix."""
    record = validate_matrix_qualification(
        value,
        source_root=source_root,
        authorized_orion_root=authorized_orion_root,
    )
    _require(
        record["registered_execution_qualification_check_pass"] is True
        and record["source_local_matrix_result_sufficient_for_downstream_qualification"]
        is False
        and all(
            value is False
            for key, value in record["authorization"].items()
            if key.endswith("_authorized")
        ),
        "registered Q043 downstream prerequisite boundary drifted",
    )
    return record
