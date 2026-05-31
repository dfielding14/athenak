#!/usr/bin/env python3
"""Extract nonqualifying Q-029 Hall-Bell source-local raw projected traces."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import platform
from pathlib import Path, PurePosixPath
import stat
from typing import Any, Callable, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q029-HALL-BELL-LINEAR"
SCHEMA_VERSION = 1
ARTIFACT_ROLE = "q029_source_local_raw_artifact_extraction_scaffold_only"
QUALIFICATION_EFFECT = "none_source_local_raw_projection_scaffold_only"
PROJECTION_CONTRACT = (
    "q023_seed_carrier_geometric_fourier_projection_only_"
    "not_hall_dispersion_or_reference_oracle"
)
PREPARATION_RECORD = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q029_hall_bell_linear_source_local_preparation_2026-05-30.json"
)
PREPARATION_RECORD_SHA256 = (
    "ff6da49895b60f0977e16231c78881e13261529e386ac795a0dee0503c92de88"
)
BINARY_READER = REPO_ROOT / "vis/python/bin_convert_new.py"
BINARY_READER_SHA256 = (
    "5835cad3233c351689c9a961d57e352613a9a731608b80a733f8f00977daac65"
)
PYTHON_VERSION = "3.11.7"
NUMPY_VERSION = "1.26.4"
APPROVED_ARTIFACT_ROOT_IDENTITY = (135357496, 720587399016549689)

_PREPARATION_KEYS = {
    "record_type",
    "schema_version",
    "date",
    "gate",
    "campaign_id",
    "status",
    "qualification_effect",
    "claim_closure",
    "frontier_authorization",
    "launch_status",
    "scope",
    "generator",
    "launch_preparation_grid",
    "artifact_bindings",
    "shared_registration_contract",
    "analyzer_contract",
    "compiled_host_contract",
    "source_local_runtime_smoke",
    "nonqualification_boundary",
    "verification",
    "explicitly_not_claimed",
    "residual_work",
    "q029_disposition",
}
_SMOKE_KEYS = {
    "qualification_effect",
    "source_commit",
    "artifact_root",
    "predecessor_magnetic_only_artifact_root",
    "retention",
    "executable_path",
    "executable_sha256",
    "cycle_zero_initializations",
    "restart_continuation",
}
_RETENTION_KEYS = {
    "file_count",
    "inventory_algorithm",
    "inventory_sha256",
    "writable_entries",
    "status",
}
_INITIALIZATION_KEYS = {
    "dimension",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
    "selected_raw_mhd_w_bcc_path",
    "selected_raw_mhd_w_bcc_sha256",
    "result",
}
_RESTART_KEYS = {
    "dimension",
    "loaded_restart_path",
    "loaded_restart_sha256",
    "uninterrupted_raw_mhd_w_bcc_path",
    "uninterrupted_raw_mhd_w_bcc_sha256",
    "continued_raw_mhd_w_bcc_path",
    "continued_raw_mhd_w_bcc_sha256",
    "exact_parity_fields",
    "max_absolute_field_difference",
    "final_cycle",
    "final_time",
    "result",
}
_BUNDLE_KEYS = {
    "schema_version",
    "campaign_id",
    "artifact_role",
    "qualification_effect",
    "qualifying_evidence",
    "hall_bell_qualification",
    "linear_hall_bell_qualification",
    "nonlinear_hall_bell_qualification",
    "projection_contract",
    "provenance",
    "traces",
}
_PROVENANCE_KEYS = {
    "preparation_record_path",
    "preparation_record_sha256",
    "artifact_root",
    "artifact_root_device",
    "artifact_root_inode",
    "artifact_inventory_sha256",
    "artifact_file_count",
    "artifact_writable_entries",
    "binary_reader_path",
    "binary_reader_sha256",
    "python_version",
    "numpy_version",
}
_TRACE_KEYS = {
    "trace_id",
    "trace_role",
    "dimension",
    "raw_geometry",
    "raw_artifacts",
    "time",
    "num_cycles",
    "magnetic_right_mode_real",
    "magnetic_right_mode_imag",
    "magnetic_left_mode_real",
    "magnetic_left_mode_imag",
    "velocity_right_mode_real",
    "velocity_right_mode_imag",
    "velocity_left_mode_real",
    "velocity_left_mode_imag",
}
_RAW_ARTIFACT_KEYS = {"sample_role", "path", "sha256"}
_RAW_GEOMETRY_KEYS = {"nx", "xmin", "extent"}
_RAW_DATASET_KEYS = {
    "MaxLevel",
    "NumCycles",
    "Time",
    "bcc1",
    "bcc2",
    "bcc3",
    "dens",
    "eint",
    "velx",
    "vely",
    "velz",
    "x1f",
    "x1v",
    "x2f",
    "x2v",
    "x3f",
    "x3v",
}
_FIELD_NAMES = ("dens", "eint", "bcc1", "bcc2", "bcc3", "velx", "vely", "velz")
_SHA256 = frozenset("0123456789abcdef")

_EXPECTED_GEOMETRY = {
    1: {
        "nx": (32, 4, 1),
        "xmin": (0.0, 0.0, 0.0),
        "extent": (1.0, 1.0, 1.0),
    },
    2: {
        "nx": (64, 32, 1),
        "xmin": (0.0, 0.0, 0.0),
        "extent": (math.sqrt(5.0), math.sqrt(1.25), 1.0),
    },
    3: {
        "nx": (128, 64, 32),
        "xmin": (0.0, 0.0, 0.0),
        "extent": (math.sqrt(21.0), math.sqrt(5.25), math.sqrt(1.3125)),
    },
}

_APPROVED_TRACE_SOURCES = (
    {
        "trace_id": "Q029-SOURCE-LOCAL-CYCLE-ZERO-1D",
        "trace_role": "cycle_zero_initialization",
        "dimension": 1,
        "raw_artifacts": (
            {
                "sample_role": "cycle_zero_initialization",
                "path": (
                    "1d/bin/"
                    "pic_q029_hall_bell_linear_1d_candidate.mhd_w_bcc.00000.bin"
                ),
                "sha256": (
                    "ad1c39450f7ac65f381b6254edaedfa9954cd32a1723c7aceae7c8f24578279d"
                ),
            },
        ),
    },
    {
        "trace_id": "Q029-SOURCE-LOCAL-CYCLE-ZERO-2D",
        "trace_role": "cycle_zero_initialization",
        "dimension": 2,
        "raw_artifacts": (
            {
                "sample_role": "cycle_zero_initialization",
                "path": (
                    "2d/bin/"
                    "pic_q029_hall_bell_linear_2d_candidate.mhd_w_bcc.00000.bin"
                ),
                "sha256": (
                    "6e165fdacb2e66910f1a522efc60717926ccde3358d35533a83755314328c227"
                ),
            },
        ),
    },
    {
        "trace_id": "Q029-SOURCE-LOCAL-CYCLE-ZERO-3D",
        "trace_role": "cycle_zero_initialization",
        "dimension": 3,
        "raw_artifacts": (
            {
                "sample_role": "cycle_zero_initialization",
                "path": (
                    "3d/bin/"
                    "pic_q029_hall_bell_linear_3d_candidate.mhd_w_bcc.00000.bin"
                ),
                "sha256": (
                    "9a0d3626135d5c924634853ffbfbfe786d53e2e57762e2b307bffefab2317bfa"
                ),
            },
        ),
    },
    {
        "trace_id": "Q029-SOURCE-LOCAL-RESTART-PARITY-2D",
        "trace_role": "restart_projected_pair_observation_from_predecessor_parity",
        "dimension": 2,
        "raw_artifacts": (
            {
                "sample_role": "uninterrupted",
                "path": (
                    "uninterrupted2d/bin/"
                    "pic_q029_hall_bell_linear_uninterrupted2d."
                    "mhd_w_bcc.00002.bin"
                ),
                "sha256": (
                    "046a41759428bdd587ddc50f4f9215e5bde263b4c73a3ea129b31e9c4709d1df"
                ),
            },
            {
                "sample_role": "continued",
                "path": (
                    "restart2d/bin/"
                    "pic_q029_hall_bell_linear_restart2d_cont."
                    "mhd_w_bcc.00003.bin"
                ),
                "sha256": (
                    "7b3c1a9a88f9edc129d550703623349f72aec86a329090061b93e4a6a8c65975"
                ),
            },
        ),
    },
)


class ContractError(ValueError):
    """Raised when the Q-029 raw extraction scaffold fails closed."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity(stat_result: os.stat_result) -> tuple[int, ...]:
    return (
        stat_result.st_dev,
        stat_result.st_ino,
        stat_result.st_mode,
        stat_result.st_size,
        stat_result.st_mtime_ns,
        stat_result.st_ctime_ns,
    )


def _read_fd_bytes(fd: int) -> bytes:
    chunks = []
    while True:
        chunk = os.read(fd, 1024 * 1024)
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)


def _write_fd_bytes(fd: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(fd, content[offset:])
        if written <= 0:
            raise ContractError("Q-029 sealed raw decode buffer write failed")
        offset += written


def _verified_regular_bytes(path: Path, label: str) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(path, flags)
    except OSError as error:
        raise ContractError(f"Q-029 {label} is unavailable") from error
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode):
            raise ContractError(f"Q-029 {label} is not a regular file")
        content = _read_fd_bytes(fd)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    if _identity(before) != _identity(after):
        raise ContractError(f"Q-029 {label} changed while being read")
    return content


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and set(value).issubset(_SHA256)
    )


def _require_exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != expected:
        raise ContractError(f"{label} keys do not match the Q-029 raw contract")
    return value


def _load_preparation_record(path: Path = PREPARATION_RECORD) -> dict[str, Any]:
    """Load the exact approved Q-029 source-local preparation predecessor."""
    try:
        raw = _verified_regular_bytes(path, "preparation record")
        digest = hashlib.sha256(raw).hexdigest()
        record = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError("Q-029 preparation record is unavailable") from error
    if digest != PREPARATION_RECORD_SHA256:
        raise ContractError("Q-029 preparation record digest mismatch")
    _require_exact_keys(record, _PREPARATION_KEYS, "Q-029 preparation record")
    if (
        record["record_type"] != "q029_hall_bell_linear_source_local_preparation"
        or type(record["schema_version"]) is not int
        or record["schema_version"] != 1
        or record["gate"] != "Q-029"
        or record["campaign_id"] != CAMPAIGN_ID
        or record["claim_closure"] is not False
    ):
        raise ContractError("Q-029 preparation record identity mismatch")
    smoke = _require_exact_keys(
        record["source_local_runtime_smoke"], _SMOKE_KEYS, "Q-029 runtime smoke"
    )
    retention = _require_exact_keys(
        smoke["retention"], _RETENTION_KEYS, "Q-029 runtime-smoke retention"
    )
    if (
        type(retention["file_count"]) is not int
        or retention["file_count"] <= 0
        or retention["inventory_algorithm"]
        != "sha256(sorted lines of '<file_sha256>  <root-relative-path>\\n')"
        or not _is_sha256(retention["inventory_sha256"])
        or retention["writable_entries"] != 0
        or retention["status"] != "pass_recursively_read_only"
    ):
        raise ContractError("Q-029 runtime-smoke retention provenance drift")
    initializations = smoke["cycle_zero_initializations"]
    if not isinstance(initializations, list) or len(initializations) != 3:
        raise ContractError("Q-029 cycle-zero initialization provenance drift")
    for initialization in initializations:
        _require_exact_keys(
            initialization, _INITIALIZATION_KEYS, "Q-029 cycle-zero initialization"
        )
    _require_exact_keys(
        smoke["restart_continuation"], _RESTART_KEYS, "Q-029 restart continuation"
    )
    _validate_preparation_artifact_bindings(smoke)
    return record


def _validate_preparation_artifact_bindings(smoke: dict[str, Any]) -> None:
    initializations = {
        item["dimension"]: item for item in smoke["cycle_zero_initializations"]
    }
    if set(initializations) != {1, 2, 3}:
        raise ContractError("Q-029 cycle-zero initialization dimensions drifted")
    for source in _APPROVED_TRACE_SOURCES[:3]:
        artifact = source["raw_artifacts"][0]
        initialization = initializations[source["dimension"]]
        if (
            initialization["selected_raw_mhd_w_bcc_path"] != artifact["path"]
            or initialization["selected_raw_mhd_w_bcc_sha256"] != artifact["sha256"]
            or initialization["result"] != "pass"
        ):
            raise ContractError("Q-029 approved cycle-zero raw artifact drift")
    restart_source = _APPROVED_TRACE_SOURCES[3]
    uninterrupted, continued = restart_source["raw_artifacts"]
    restart = smoke["restart_continuation"]
    if (
        restart["dimension"] != restart_source["dimension"]
        or restart["uninterrupted_raw_mhd_w_bcc_path"] != uninterrupted["path"]
        or restart["uninterrupted_raw_mhd_w_bcc_sha256"] != uninterrupted["sha256"]
        or restart["continued_raw_mhd_w_bcc_path"] != continued["path"]
        or restart["continued_raw_mhd_w_bcc_sha256"] != continued["sha256"]
        or restart["max_absolute_field_difference"] != 0.0
        or restart["result"] != "pass_exact_array_parity"
    ):
        raise ContractError("Q-029 approved restart-parity raw artifact drift")


def _inventory_sha256_fd(root_fd: int) -> tuple[int, int, str, dict[str, bytes]]:
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    rows = []
    file_bytes = {}
    file_count = 0
    writable_entries = 0

    def visit(directory_fd: int, prefix: str) -> None:
        nonlocal file_count, writable_entries
        before = os.fstat(directory_fd)
        if before.st_mode & 0o222:
            writable_entries += 1
        try:
            names = sorted(os.listdir(directory_fd))
        except OSError as error:
            raise ContractError(
                "Q-029 authorized artifact inventory is unavailable"
            ) from error
        for name in names:
            relative = f"{prefix}/{name}" if prefix else name
            try:
                observed = os.stat(
                    name, dir_fd=directory_fd, follow_symlinks=False
                )
            except OSError as error:
                raise ContractError(
                    "Q-029 authorized artifact inventory is unavailable"
                ) from error
            if stat.S_ISLNK(observed.st_mode):
                raise ContractError("Q-029 authorized artifact root contains a symlink")
            if stat.S_ISDIR(observed.st_mode):
                try:
                    child_fd = os.open(name, flags, dir_fd=directory_fd)
                except OSError as error:
                    raise ContractError(
                        "Q-029 authorized artifact directory changed during inventory"
                    ) from error
                try:
                    if _identity(observed) != _identity(os.fstat(child_fd)):
                        raise ContractError(
                            "Q-029 authorized artifact directory identity drift"
                        )
                    visit(child_fd, relative)
                    if _identity(observed) != _identity(os.fstat(child_fd)):
                        raise ContractError(
                            "Q-029 authorized artifact directory changed during inventory"
                        )
                finally:
                    os.close(child_fd)
            elif stat.S_ISREG(observed.st_mode):
                if observed.st_mode & 0o222:
                    writable_entries += 1
                file_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
                try:
                    file_fd = os.open(name, file_flags, dir_fd=directory_fd)
                except OSError as error:
                    raise ContractError(
                        "Q-029 authorized artifact file changed during inventory"
                    ) from error
                try:
                    if _identity(observed) != _identity(os.fstat(file_fd)):
                        raise ContractError(
                            "Q-029 authorized artifact file identity drift"
                        )
                    raw = _read_fd_bytes(file_fd)
                    if _identity(observed) != _identity(os.fstat(file_fd)):
                        raise ContractError(
                            "Q-029 authorized artifact file changed during inventory"
                        )
                finally:
                    os.close(file_fd)
                rows.append(f"{hashlib.sha256(raw).hexdigest()}  {relative}\n")
                file_bytes[relative] = raw
                file_count += 1
            else:
                raise ContractError(
                    "Q-029 authorized artifact root contains a special file"
                )
        if _identity(before) != _identity(os.fstat(directory_fd)):
            raise ContractError(
                "Q-029 authorized artifact directory changed during inventory"
            )

    visit(root_fd, "")
    digest = hashlib.sha256("".join(rows).encode("utf-8")).hexdigest()
    return file_count, writable_entries, digest, file_bytes


def _inventory_sha256(root: Path) -> tuple[int, int, str]:
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    try:
        root_fd = os.open(root, flags)
    except OSError as error:
        raise ContractError("Q-029 authorized artifact inventory is unavailable") from error
    try:
        file_count, writable_entries, digest, _ = _inventory_sha256_fd(root_fd)
    finally:
        os.close(root_fd)
    return file_count, writable_entries, digest


def _revalidate_authorized_artifact_root(root: Path, root_fd: int) -> None:
    try:
        path_state = os.stat(root, follow_symlinks=False)
        descriptor_state = os.fstat(root_fd)
    except OSError as error:
        raise ContractError("Q-029 authorized artifact root is unavailable") from error
    if (
        (descriptor_state.st_dev, descriptor_state.st_ino)
        != APPROVED_ARTIFACT_ROOT_IDENTITY
        or _identity(path_state) != _identity(descriptor_state)
    ):
        raise ContractError("Q-029 authorized artifact-root identity drift")


def _authorized_artifact_root(
    artifact_root: Path, smoke: dict[str, Any]
) -> tuple[Path, int, dict[str, bytes]]:
    root = Path(artifact_root)
    if not root.is_absolute():
        raise ContractError("Q-029 authorized artifact root must be absolute")
    approved = Path(smoke["artifact_root"])
    if not approved.is_absolute() or root != approved:
        raise ContractError("Q-029 artifact root is not the approved source-local root")
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    try:
        root_fd = os.open(root, flags)
    except OSError as error:
        raise ContractError("Q-029 authorized artifact root is unavailable") from error
    try:
        _revalidate_authorized_artifact_root(root, root_fd)
        file_count, writable_entries, inventory_sha256, file_bytes = (
            _inventory_sha256_fd(root_fd)
        )
        retention = smoke["retention"]
        if (
            file_count != retention["file_count"]
            or writable_entries != retention["writable_entries"]
            or inventory_sha256 != retention["inventory_sha256"]
        ):
            raise ContractError("Q-029 authorized artifact-root provenance drift")
        return root, root_fd, file_bytes
    except Exception:
        os.close(root_fd)
        raise


def _normalized_artifact_bytes(
    file_bytes: dict[str, bytes], relative_path: str, expected_sha256: str
) -> bytes:
    if not isinstance(relative_path, str) or not relative_path:
        raise ContractError("Q-029 raw artifact path is malformed")
    pure = PurePosixPath(relative_path)
    if (
        pure.is_absolute()
        or pure.as_posix() != relative_path
        or any(part in {"", ".", ".."} for part in pure.parts)
    ):
        raise ContractError("Q-029 raw artifact path must be normalized root-relative")
    if not _is_sha256(expected_sha256):
        raise ContractError("Q-029 raw artifact digest is malformed")
    try:
        raw = file_bytes[relative_path]
    except KeyError as error:
        raise ContractError("Q-029 raw artifact path is not an approved file") from error
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ContractError("Q-029 raw artifact digest mismatch")
    if not isinstance(raw, bytes):
        raise ContractError("Q-029 raw artifact path is not an approved file")
    return raw


def _raw_geometry(dimension: int) -> dict[str, list[int] | list[float]]:
    try:
        geometry = _EXPECTED_GEOMETRY[dimension]
    except KeyError as error:
        raise ContractError("Q-029 raw trace dimension is not approved") from error
    return {
        "nx": list(geometry["nx"]),
        "xmin": list(geometry["xmin"]),
        "extent": list(geometry["extent"]),
    }


def _mode_basis(dimension: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if dimension not in _EXPECTED_GEOMETRY:
        raise ContractError("Q-029 raw trace dimension is not approved")
    raw = np.array(
        [1.0, 2.0 if dimension >= 2 else 0.0, 4.0 if dimension >= 3 else 0.0]
    )
    parallel = raw / np.linalg.norm(raw)
    transverse_a = (
        np.array([0.0, 1.0, 0.0])
        if dimension == 1
        else np.array([-parallel[1], parallel[0], 0.0])
    )
    transverse_a /= np.linalg.norm(transverse_a)
    transverse_b = np.cross(parallel, transverse_a)
    return parallel, transverse_a, transverse_b


def _finite_scalar(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ContractError(f"Q-029 raw dataset {label} must be a finite scalar")
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"Q-029 raw dataset {label} must be a finite scalar")
    return measured


def _validate_raw_dataset(
    dataset: dict[str, Any], dimension: int
) -> tuple[list[np.ndarray], tuple[int, int, int]]:
    _require_exact_keys(dataset, _RAW_DATASET_KEYS, "Q-029 raw mhd_w_bcc dataset")
    geometry = _EXPECTED_GEOMETRY[dimension]
    nx = geometry["nx"]
    xmin = geometry["xmin"]
    extent = geometry["extent"]
    coordinates = []
    for axis in range(3):
        faces = np.asarray(dataset[f"x{axis + 1}f"])
        centers = np.asarray(dataset[f"x{axis + 1}v"])
        if (
            faces.ndim != 1
            or centers.ndim != 1
            or not np.issubdtype(faces.dtype, np.floating)
            or not np.issubdtype(centers.dtype, np.floating)
            or not np.all(np.isfinite(faces))
            or not np.all(np.isfinite(centers))
            or faces.size != nx[axis] + 1
            or centers.size != nx[axis]
        ):
            raise ContractError("Q-029 raw mhd_w_bcc geometry drift")
        expected_faces = (
            xmin[axis]
            + np.arange(nx[axis] + 1, dtype=float) * (extent[axis] / nx[axis])
        ).astype(faces.dtype)
        expected_centers = (
            (faces[:-1] + faces[1:]) / np.asarray(2, dtype=faces.dtype)
        ).astype(centers.dtype)
        if (
            not np.array_equal(faces, expected_faces)
            or not np.array_equal(centers, expected_centers)
        ):
            raise ContractError("Q-029 raw mhd_w_bcc geometry drift")
        coordinates.append(centers)

    expected_shape = tuple(reversed(nx))
    for name in _FIELD_NAMES:
        values = np.asarray(dataset[name])
        if (
            values.shape != expected_shape
            or not np.issubdtype(values.dtype, np.floating)
            or not np.all(np.isfinite(values))
        ):
            raise ContractError(f"Q-029 raw mhd_w_bcc {name} field drift")
    max_level = dataset["MaxLevel"]
    num_cycles = dataset["NumCycles"]
    if (
        isinstance(max_level, bool)
        or not isinstance(max_level, (int, np.integer))
        or int(max_level) != 0
    ):
        raise ContractError("Q-029 raw dataset MaxLevel drift")
    if isinstance(num_cycles, bool) or not isinstance(num_cycles, (int, np.integer)):
        raise ContractError("Q-029 raw dataset NumCycles drift")
    _finite_scalar(dataset["Time"], "Time")
    return coordinates, expected_shape


def _project_components(
    dataset: dict[str, Any],
    coordinates: Sequence[np.ndarray],
    dimension: int,
    names: Sequence[str],
) -> tuple[complex, complex]:
    parallel, transverse_a, transverse_b = _mode_basis(dimension)
    vector = np.stack([np.asarray(dataset[name], dtype=float) for name in names])
    x3, x2, x1 = np.meshgrid(
        coordinates[2], coordinates[1], coordinates[0], indexing="ij"
    )
    spatial_phase = 2.0 * math.pi * (
        parallel[0] * x1 + parallel[1] * x2 + parallel[2] * x3
    )
    fourier_weight = np.exp(-1.0j * spatial_phase)
    mode_a = np.mean(np.tensordot(transverse_a, vector, axes=1) * fourier_weight)
    mode_b = np.mean(np.tensordot(transverse_b, vector, axes=1) * fourier_weight)
    right = 0.5 * (mode_a - 1.0j * mode_b)
    left = 0.5 * (mode_a + 1.0j * mode_b)
    if not all(math.isfinite(value) for value in (
        right.real, right.imag, left.real, left.imag
    )):
        raise ContractError("Q-029 projected raw trace contains a nonfinite value")
    return right, left


def _project_dataset(dataset: dict[str, Any], dimension: int) -> dict[str, Any]:
    """Project one exact raw snapshot without applying a physical oracle."""
    coordinates, _ = _validate_raw_dataset(dataset, dimension)
    magnetic_right, magnetic_left = _project_components(
        dataset, coordinates, dimension, ("bcc1", "bcc2", "bcc3")
    )
    velocity_right, velocity_left = _project_components(
        dataset, coordinates, dimension, ("velx", "vely", "velz")
    )
    return {
        "time": _finite_scalar(dataset["Time"], "Time"),
        "num_cycles": int(dataset["NumCycles"]),
        "magnetic_right_mode_real": float(magnetic_right.real),
        "magnetic_right_mode_imag": float(magnetic_right.imag),
        "magnetic_left_mode_real": float(magnetic_left.real),
        "magnetic_left_mode_imag": float(magnetic_left.imag),
        "velocity_right_mode_real": float(velocity_right.real),
        "velocity_right_mode_imag": float(velocity_right.imag),
        "velocity_left_mode_real": float(velocity_left.real),
        "velocity_left_mode_imag": float(velocity_left.imag),
    }


def _binary_reader() -> Callable[[str], dict[str, Any]]:
    raw = _verified_regular_bytes(BINARY_READER, "raw Athena binary reader")
    if hashlib.sha256(raw).hexdigest() != BINARY_READER_SHA256:
        raise ContractError("Q-029 raw Athena binary reader digest mismatch")
    namespace = {
        "__file__": str(BINARY_READER),
        "__name__": "q029_bin_convert_new",
    }
    try:
        exec(compile(raw, str(BINARY_READER), "exec"), namespace)
        reader = namespace["read_binary_as_athdf"]
    except (KeyError, SyntaxError) as error:
        raise ContractError("unable to load the Q-029 raw Athena binary reader")
    if not callable(reader):
        raise ContractError("unable to load the Q-029 raw Athena binary reader")
    return reader


def _read_verified_dataset(
    raw: bytes,
    expected_sha256: str,
    reader: Callable[[str], dict[str, Any]],
) -> dict[str, Any]:
    """Decode one sealed byte copy after checking the approved raw digest."""
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ContractError("Q-029 raw artifact digest mismatch")
    seals = (
        fcntl.F_SEAL_SEAL
        | fcntl.F_SEAL_SHRINK
        | fcntl.F_SEAL_GROW
        | fcntl.F_SEAL_WRITE
    )
    fd = os.memfd_create("q029-raw-decode", os.MFD_ALLOW_SEALING)
    try:
        _write_fd_bytes(fd, raw)
        fcntl.fcntl(fd, fcntl.F_ADD_SEALS, seals)
        if fcntl.fcntl(fd, fcntl.F_GET_SEALS) != seals:
            raise ContractError("Q-029 raw decode buffer sealing failed")
        os.lseek(fd, 0, os.SEEK_SET)
        if _read_fd_bytes(fd) != raw:
            raise ContractError("Q-029 sealed raw decode buffer drift")
        return reader(f"/proc/self/fd/{fd}")
    finally:
        os.close(fd)


def _provenance(root: Path, smoke: dict[str, Any]) -> dict[str, Any]:
    if platform.python_version() != PYTHON_VERSION or np.__version__ != NUMPY_VERSION:
        raise ContractError("Q-029 numerical runtime binding mismatch")
    return {
        "preparation_record_path": str(PREPARATION_RECORD.relative_to(REPO_ROOT)),
        "preparation_record_sha256": PREPARATION_RECORD_SHA256,
        "artifact_root": root.as_posix(),
        "artifact_root_device": APPROVED_ARTIFACT_ROOT_IDENTITY[0],
        "artifact_root_inode": APPROVED_ARTIFACT_ROOT_IDENTITY[1],
        "artifact_inventory_sha256": smoke["retention"]["inventory_sha256"],
        "artifact_file_count": smoke["retention"]["file_count"],
        "artifact_writable_entries": smoke["retention"]["writable_entries"],
        "binary_reader_path": str(BINARY_READER.relative_to(REPO_ROOT)),
        "binary_reader_sha256": BINARY_READER_SHA256,
        "python_version": PYTHON_VERSION,
        "numpy_version": NUMPY_VERSION,
    }


def extract_raw_trace_bundle(
    artifact_root: Path,
) -> dict[str, Any]:
    """Extract observed projections from approved bytes with the pinned decoder."""
    preparation = _load_preparation_record()
    smoke = preparation["source_local_runtime_smoke"]
    root, root_fd, file_bytes = _authorized_artifact_root(artifact_root, smoke)
    try:
        reader = _binary_reader()
        traces = []
        for source in _APPROVED_TRACE_SOURCES:
            rows = []
            raw_artifacts = []
            for artifact in source["raw_artifacts"]:
                raw = _normalized_artifact_bytes(
                    file_bytes, artifact["path"], artifact["sha256"]
                )
                dataset = _read_verified_dataset(raw, artifact["sha256"], reader)
                rows.append(_project_dataset(dataset, source["dimension"]))
                raw_artifacts.append(dict(artifact))
            traces.append(
                {
                "trace_id": source["trace_id"],
                "trace_role": source["trace_role"],
                "dimension": source["dimension"],
                "raw_geometry": _raw_geometry(source["dimension"]),
                "raw_artifacts": raw_artifacts,
                "time": [row["time"] for row in rows],
                "num_cycles": [row["num_cycles"] for row in rows],
                "magnetic_right_mode_real": [
                    row["magnetic_right_mode_real"] for row in rows
                ],
                "magnetic_right_mode_imag": [
                    row["magnetic_right_mode_imag"] for row in rows
                ],
                "magnetic_left_mode_real": [
                    row["magnetic_left_mode_real"] for row in rows
                ],
                "magnetic_left_mode_imag": [
                    row["magnetic_left_mode_imag"] for row in rows
                ],
                "velocity_right_mode_real": [
                    row["velocity_right_mode_real"] for row in rows
                ],
                "velocity_right_mode_imag": [
                    row["velocity_right_mode_imag"] for row in rows
                ],
                "velocity_left_mode_real": [
                    row["velocity_left_mode_real"] for row in rows
                ],
                "velocity_left_mode_imag": [
                    row["velocity_left_mode_imag"] for row in rows
                ],
                }
            )
        result = {
            "schema_version": SCHEMA_VERSION,
            "campaign_id": CAMPAIGN_ID,
            "artifact_role": ARTIFACT_ROLE,
            "qualification_effect": QUALIFICATION_EFFECT,
            "qualifying_evidence": False,
            "hall_bell_qualification": False,
            "linear_hall_bell_qualification": False,
            "nonlinear_hall_bell_qualification": False,
            "projection_contract": PROJECTION_CONTRACT,
            "provenance": _provenance(root, smoke),
            "traces": traces,
        }
        _revalidate_authorized_artifact_root(root, root_fd)
        return result
    finally:
        os.close(root_fd)


def _finite_trace_array(trace: dict[str, Any], key: str, size: int) -> None:
    values = trace[key]
    if not isinstance(values, list) or len(values) != size:
        raise ContractError(f"Q-029 projected trace {key} shape drift")
    for value in values:
        _finite_scalar(value, key)


def validate_raw_trace_bundle(
    bundle: dict[str, Any],
    artifact_root: Path,
) -> dict[str, Any]:
    """Replay exact artifacts and reject extracted-bundle schema or value drift."""
    _require_exact_keys(bundle, _BUNDLE_KEYS, "Q-029 projected raw-trace bundle")
    if (
        type(bundle["schema_version"]) is not int
        or bundle["schema_version"] != SCHEMA_VERSION
        or bundle["campaign_id"] != CAMPAIGN_ID
        or bundle["artifact_role"] != ARTIFACT_ROLE
        or bundle["qualification_effect"] != QUALIFICATION_EFFECT
        or bundle["qualifying_evidence"] is not False
        or bundle["hall_bell_qualification"] is not False
        or bundle["linear_hall_bell_qualification"] is not False
        or bundle["nonlinear_hall_bell_qualification"] is not False
        or bundle["projection_contract"] != PROJECTION_CONTRACT
    ):
        raise ContractError("Q-029 projected raw-trace bundle identity drift")
    _require_exact_keys(bundle["provenance"], _PROVENANCE_KEYS, "Q-029 raw provenance")
    traces = bundle["traces"]
    if not isinstance(traces, list) or len(traces) != len(_APPROVED_TRACE_SOURCES):
        raise ContractError("Q-029 projected raw-trace inventory drift")
    expected = extract_raw_trace_bundle(artifact_root)
    if bundle["provenance"] != expected["provenance"]:
        raise ContractError("Q-029 raw provenance drift")
    for trace, expected_trace in zip(traces, expected["traces"]):
        _require_exact_keys(trace, _TRACE_KEYS, "Q-029 projected trace")
        if (
            trace["trace_id"] != expected_trace["trace_id"]
            or trace["trace_role"] != expected_trace["trace_role"]
            or trace["dimension"] != expected_trace["dimension"]
            or trace["raw_geometry"] != expected_trace["raw_geometry"]
            or trace["raw_artifacts"] != expected_trace["raw_artifacts"]
        ):
            raise ContractError("Q-029 projected trace provenance drift")
        if (
            not isinstance(trace["raw_geometry"], dict)
            or set(trace["raw_geometry"]) != _RAW_GEOMETRY_KEYS
        ):
            raise ContractError("Q-029 projected trace geometry schema drift")
        if any(
            not isinstance(artifact, dict) or set(artifact) != _RAW_ARTIFACT_KEYS
            for artifact in trace["raw_artifacts"]
        ):
            raise ContractError("Q-029 raw artifact provenance schema drift")
        size = len(trace["raw_artifacts"])
        if (
            not isinstance(trace["num_cycles"], list)
            or len(trace["num_cycles"]) != size
            or any(type(value) is not int for value in trace["num_cycles"])
        ):
            raise ContractError("Q-029 projected trace NumCycles drift")
        for key in _TRACE_KEYS - {
            "trace_id",
            "trace_role",
            "dimension",
            "raw_geometry",
            "raw_artifacts",
            "num_cycles",
        }:
            _finite_trace_array(trace, key, size)
        if trace != expected_trace:
            raise ContractError("Q-029 projected raw trace value drift")
    return bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--validate-bundle", type=Path)
    args = parser.parse_args()
    if args.validate_bundle is None:
        result = extract_raw_trace_bundle(args.artifact_root)
    else:
        try:
            bundle = json.loads(args.validate_bundle.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ContractError("Q-029 projected raw-trace bundle is unavailable") from error
        result = validate_raw_trace_bundle(bundle, args.artifact_root)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
