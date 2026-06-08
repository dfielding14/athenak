#!/usr/bin/env python3
"""Focused tests for the corrected Q023 J/c Bell linear predecessor."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from tst.publication import analyze_q023_paper_bell_linear as legacy
from tst.publication import analyze_q023_paper_bell_linear_joverc as bell


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q023_paper_bell_linear_joverc_source_local_preparation_2026-06-06.json"
)
HARNESS = REPO_ROOT / "tst/publication/q023_paper_bell_linear_joverc_host_harness.cpp"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(
    bundle: dict[str, object], member_id: str, epsilon: float
) -> dict[str, object]:
    return next(
        record
        for record in bundle["records"]
        if str(record["member_id"]).startswith(member_id + "-eps")
        and record["epsilon"] == epsilon
    )


def _set_physics_trace(
    record: dict[str, object],
    *,
    growth: float | None = None,
    phase: float | None = None,
    left_ratio: float = 100.0,
) -> None:
    epsilon = float(record["epsilon"])
    expected_phase, expected_growth = bell.theoretical_dispersion(epsilon)
    growth = expected_growth if growth is None else growth
    phase = expected_phase if phase is None else phase
    time = np.linspace(0.0, 6.0 / expected_growth, 49)
    right = np.exp(growth * time) * np.exp(1.0j * phase * time)
    left = right / left_ratio
    velocity_right = (-epsilon - 1.0j * expected_growth) * right
    velocity_left = velocity_right / left_ratio
    interval, change = legacy._fixed_interval_phase_trace(time, right)
    velocity_interval, velocity_change = legacy._fixed_interval_phase_trace(
        time, velocity_right
    )
    paper_interval, paper_change = legacy._fixed_interval_phase_trace(time, right)
    record["physics_trace"] = {
        "normalized_time": time.tolist(),
        "right_mode_real": right.real.tolist(),
        "right_mode_imag": right.imag.tolist(),
        "left_mode_real": left.real.tolist(),
        "left_mode_imag": left.imag.tolist(),
        "velocity_right_mode_real": velocity_right.real.tolist(),
        "velocity_right_mode_imag": velocity_right.imag.tolist(),
        "velocity_left_mode_real": velocity_left.real.tolist(),
        "velocity_left_mode_imag": velocity_left.imag.tolist(),
        "phase_interval": interval.tolist(),
        "phase_change": change.tolist(),
        "velocity_phase_interval": velocity_interval.tolist(),
        "velocity_phase_change": velocity_change.tolist(),
        "paper_delta_u_y_sine_fit_real": right.real.tolist(),
        "paper_delta_u_y_sine_fit_imag": right.imag.tolist(),
        "paper_delta_u_y_phase_interval": paper_interval.tolist(),
        "paper_delta_u_y_phase_change": paper_change.tolist(),
        "paper_volume_averaged_abs_delta_u": np.abs(velocity_right).tolist(),
    }


def _runtime_header(
    member: dict[str, object], snapshot_index: int, variable_output_index: int
) -> bytes:
    text = bell.render_deck(member)
    for index in range(1, 7):
        if index <= 5:
            file_number = (
                snapshot_index if index <= variable_output_index else snapshot_index + 1
            )
            last_time = (
                -1.0
                if snapshot_index == 0 and index <= variable_output_index
                else (
                    snapshot_index * bell.LINEAR_OUTPUT_DT
                    if index > variable_output_index
                    else (snapshot_index - 1) * bell.LINEAR_OUTPUT_DT
                )
            )
        else:
            file_number = 0 if snapshot_index == 0 else 1
            last_time = -1.0 if snapshot_index == 0 else 0.0
        text = text.replace(
            f"<output{index}>\n",
            f"<output{index}>\nfile_number = {file_number}\nlast_time = {last_time:.17g}\n",
            1,
        )
    return text.encode("utf-8")


def _mode_fields(member: dict[str, object], physical_time: float) -> dict[str, np.ndarray]:
    nx1, nx2, nx3 = (int(value) for value in member["global_nx"])
    bounds = tuple(tuple(float(value) for value in axis) for axis in member["bounds"])
    coordinates = [
        np.linspace(lower, upper, count, endpoint=False)
        + 0.5 * (upper - lower) / count
        for (lower, upper), count in zip(bounds, (nx1, nx2, nx3))
    ]
    x3, x2, x1 = np.meshgrid(coordinates[2], coordinates[1], coordinates[0], indexing="ij")
    parallel = np.asarray(bell._mode_basis(int(member["dimension"])))
    transverse_a = np.array([0.0, 1.0, 0.0])
    transverse_b = np.cross(parallel, transverse_a)
    if int(member["dimension"]) > 1:
        transverse_a = np.array([-parallel[1], parallel[0], 0.0])
        transverse_a /= np.linalg.norm(transverse_a)
        transverse_b = np.cross(parallel, transverse_a)
    spatial_phase = bell.K0 * (
        parallel[0] * x1 + parallel[1] * x2 + parallel[2] * x3
    )
    epsilon = float(member["epsilon"])
    _, growth = bell.theoretical_dispersion(epsilon)
    normalized_time = bell.K0 * physical_time
    phase = spatial_phase - epsilon * normalized_time
    amplitude = 1.0e-6 * np.exp(growth * normalized_time)
    magnetic_a = amplitude * np.cos(phase)
    magnetic_b = amplitude * np.sin(phase)
    right = 0.5 * amplitude * np.exp(-1.0j * epsilon * normalized_time)
    velocity_right = complex(-epsilon, -growth) * right
    velocity_a = 2.0 * np.real(velocity_right * np.exp(1.0j * spatial_phase))
    velocity_b = 2.0 * np.real(-1.0j * velocity_right * np.exp(1.0j * spatial_phase))
    magnetic = (
        parallel[:, None, None, None]
        + transverse_a[:, None, None, None] * magnetic_a
        + transverse_b[:, None, None, None] * magnetic_b
    )
    velocity = (
        transverse_a[:, None, None, None] * velocity_a
        + transverse_b[:, None, None, None] * velocity_b
    )
    shape = (nx3, nx2, nx1)
    return {
        "dens": np.ones(shape),
        "eint": np.ones(shape),
        "velx": velocity[0],
        "vely": velocity[1],
        "velz": velocity[2],
        "bcc1": magnetic[0],
        "bcc2": magnetic[1],
        "bcc3": magnetic[2],
    }


def _raw_payload_from_fields(
    member: dict[str, object],
    physical_time: float,
    cycle: int,
    *,
    fields: tuple[str, ...],
    values: dict[str, np.ndarray],
    variable_output_index: int,
    snapshot_index: int | None = None,
) -> bytes:
    meshblock = tuple(int(value) for value in member["meshblock_nx"])
    splits = tuple(int(value) for value in member["decomposition_splits"])
    bounds = tuple(tuple(float(value) for value in axis) for axis in member["bounds"])
    blocks = []
    for logical_x3 in range(splits[2]):
        for logical_x2 in range(splits[1]):
            for logical_x1 in range(splits[0]):
                logical = (logical_x1, logical_x2, logical_x3)
                geometry = tuple(
                    coordinate
                    for axis, split_index in enumerate(logical)
                    for coordinate in (
                        bounds[axis][0]
                        + split_index * (bounds[axis][1] - bounds[axis][0]) / splits[axis],
                        bounds[axis][0]
                        + (split_index + 1)
                        * (bounds[axis][1] - bounds[axis][0])
                        / splits[axis],
                    )
                )
                i0 = logical_x1 * meshblock[0]
                j0 = logical_x2 * meshblock[1]
                k0 = logical_x3 * meshblock[2]
                slices = (
                    slice(k0, k0 + meshblock[2]),
                    slice(j0, j0 + meshblock[1]),
                    slice(i0, i0 + meshblock[0]),
                )
                blocks.append(
                    struct.pack(
                        "<10i",
                        0,
                        meshblock[0] - 1,
                        0,
                        meshblock[1] - 1,
                        0,
                        meshblock[2] - 1,
                        *logical,
                        0,
                    )
                    + struct.pack("<6d", *geometry)
                    + np.concatenate(
                        [
                            np.asarray(values[field][slices], dtype="<f4").ravel()
                            for field in fields
                        ]
                    ).tobytes()
                )
    parameter_header = _runtime_header(
        member,
        cycle if snapshot_index is None else snapshot_index,
        variable_output_index,
    )
    return (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={physical_time:.17g}\n".encode()
        + f"  cycle={cycle}\n".encode()
        + b"  size of location=8\n"
        b"  size of variable=4\n"
        + f"  number of variables={len(fields)}\n".encode()
        + b"  variables:  "
        + b"  ".join(field.encode() for field in fields)
        + b"  \n"
        + f"  header offset={len(parameter_header)}\n".encode()
        + parameter_header
        + b"".join(blocks)
    )


def _raw_payload(
    member: dict[str, object],
    physical_time: float,
    cycle: int,
    *,
    snapshot_index: int | None = None,
) -> bytes:
    fields = ("dens", "eint", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3")
    return _raw_payload_from_fields(
        member,
        physical_time,
        cycle,
        fields=fields,
        values=_mode_fields(member, physical_time),
        variable_output_index=1,
        snapshot_index=snapshot_index,
    )


def _particle_raw_payload(
    member: dict[str, object],
    physical_time: float,
    cycle: int,
    variable: str,
    *,
    snapshot_index: int | None = None,
) -> bytes:
    shape = tuple(reversed(tuple(int(value) for value in member["global_nx"])))
    if cycle == 0:
        value = 0.0
    elif variable == "prtcl_rho":
        value = (
            int(member["ppc"])
            * float(member["deposit_qscale"])
            * float(member["species_charge"])
            / float(member["root_cell_volume"])
        )
    else:
        component = {"prtcl_jx": 0, "prtcl_jy": 1, "prtcl_jz": 2}[variable]
        value = bell.EXPECTED_J_OVER_C * bell._mode_basis(int(member["dimension"]))[
            component
        ]
    return _raw_payload_from_fields(
        member,
        physical_time,
        cycle,
        fields=(variable,),
        values={variable: np.full(shape, value)},
        variable_output_index=bell._RAW_OUTPUT_INDEX[variable],
        snapshot_index=snapshot_index,
    )


def _write_materialized_member(
    root: Path,
    member: dict[str, object],
    dependency: dict[str, object],
    *,
    cycles: list[int] | None = None,
    physical_times: list[float] | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    executable = root / "candidate/athena"
    executable.parent.mkdir(parents=True, exist_ok=True)
    executable.write_bytes(b"clean executable fixture\n")
    stdout = root / "execution/stdout.txt"
    stdout.parent.mkdir(parents=True, exist_ok=True)
    stdout.write_text("athenak_driver_completed_successfully\n", encoding="utf-8")
    if physical_times is None:
        physical_times = [
            *(index * bell.LINEAR_OUTPUT_DT for index in range(88)),
            bell.LINEAR_RUNTIME_TLIM,
        ]
    if len(physical_times) != 89:
        raise ValueError("materialized fixture time inventory must contain 89 entries")
    measured_cycles = list(range(len(physical_times))) if cycles is None else cycles
    if len(measured_cycles) != len(physical_times):
        raise ValueError("materialized fixture cycle inventory must contain 89 entries")
    raw_artifacts = []
    datasets = []
    basename = "q023_joverc_" + str(member["member_id"]).replace("-", "_")
    for output_index, (cycle, physical_time) in enumerate(
        zip(measured_cycles, physical_times)
    ):
        path = root / f"raw/bin/{basename}.mhd_w_bcc.{output_index:05d}.bin"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            _raw_payload(
                member,
                physical_time,
                cycle,
                snapshot_index=output_index,
            )
        )
        raw_artifacts.append(
            {
                "path": path.relative_to(root).as_posix(),
                "sha256": _sha256(path),
                "variable": "mhd_w_bcc",
                "cycle": cycle,
                "time": physical_time,
            }
        )
        datasets.append(
            bell.binary.parse_athenak_binary_bytes(path.read_bytes(), source=str(path))
        )
        for variable in bell._PARTICLE_FIELDS:
            particle_path = (
                root / f"raw/bin/{basename}.{variable}.{output_index:05d}.bin"
            )
            particle_path.write_bytes(
                _particle_raw_payload(
                    member,
                    physical_time,
                    cycle,
                    variable,
                    snapshot_index=output_index,
                )
            )
            raw_artifacts.append(
                {
                    "path": particle_path.relative_to(root).as_posix(),
                    "sha256": _sha256(particle_path),
                    "variable": variable,
                    "cycle": cycle,
                    "time": physical_time,
                }
            )
    physics_trace = bell._physics_trace_from_raw_datasets(datasets, member=member)
    rank_count = math.prod(int(value) for value in member["decomposition_splits"])
    topology_path = root / "execution/rank_topology.json"
    splits = tuple(int(value) for value in member["decomposition_splits"])
    locations = [
        [x1, x2, x3]
        for x3 in range(splits[2])
        for x2 in range(splits[1])
        for x1 in range(splits[0])
    ]
    topology_path.write_text(
        json.dumps(
            {
                "schema_version": bell.SCHEMA_VERSION,
                "record_type": "q023_paper_bell_linear_joverc_mpi_rank_topology_evidence",
                "campaign_id": bell.CAMPAIGN_ID,
                "member_id": member["member_id"],
                "rank_count": rank_count,
                "decomposition_splits": member["decomposition_splits"],
                "ranks": [
                    {"rank": rank, "meshblock_logical_locations": [location]}
                    for rank, location in enumerate(locations)
                ],
                "publication_authorized": False,
            }
        ),
        encoding="utf-8",
    )
    deck_path = (
        bell.DECK_ROOT.relative_to(REPO_ROOT) / str(member["deck_path"])
    ).as_posix()
    receipt_path = root / "execution/receipt.json"
    receipt = {
        "schema_version": bell.SCHEMA_VERSION,
        "record_type": "q023_paper_bell_linear_joverc_registered_execution_receipt",
        "campaign_id": bell.CAMPAIGN_ID,
        "member_id": member["member_id"],
        "deck_path": deck_path,
        "deck_sha256": member["deck_sha256"],
        "source_path": bell.SOURCE_PATH.as_posix(),
        "source_sha256": _sha256(REPO_ROOT / bell.SOURCE_PATH),
        "corrected_eigenmode_header_path": bell.CORRECTED_EIGENMODE_HEADER_PATH.as_posix(),
        "corrected_eigenmode_header_sha256": _sha256(
            REPO_ROOT / bell.CORRECTED_EIGENMODE_HEADER_PATH
        ),
        "executable_path": "candidate/athena",
        "executable_sha256": _sha256(executable),
        "candidate_clean": True,
        "source_clean": True,
        "executable_clean": True,
        "command": {
            "argv": [
                "srun",
                "-n",
                str(rank_count),
                "candidate/athena",
                "-i",
                deck_path,
            ],
            "working_directory": str(root),
        },
        "mpi": {
            "launcher": "srun",
            "rank_count": rank_count,
            "decomposition_splits": member["decomposition_splits"],
            "rank_topology_path": "execution/rank_topology.json",
            "rank_topology_sha256": _sha256(topology_path),
        },
        "termination": {
            "status": "completed",
            "return_code": 0,
            "final_cycle": 100,
            "final_time": bell.LINEAR_RUNTIME_TLIM,
        },
        "stdout": {
            "path": "execution/stdout.txt",
            "sha256": _sha256(stdout),
            "completion_marker": "athenak_driver_completed_successfully",
        },
        "raw_artifacts": raw_artifacts,
        "publication_authorized": False,
    }
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    provenance = bell.synthetic_provenance(member, dependency)
    provenance.update(
        {
            "kind": "registered_execution_trace",
            "authorized_artifact_root": str(root),
            "candidate_clean": True,
            "source_clean": True,
            "executable_clean": True,
            "executable_path": "candidate/athena",
            "executable_sha256": _sha256(executable),
            "registered_execution_receipt_path": "execution/receipt.json",
            "registered_execution_receipt_sha256": _sha256(receipt_path),
            "raw_artifacts": raw_artifacts,
        }
    )
    return provenance, physics_trace


def _q043_bound_dependency_fixture_for_q023_provenance_only() -> dict[str, object]:
    dependency = bell.synthetic_q043_registered_raw_oracle_dependency()
    dependency["registered_admission_digest_bound"] = True
    dependency["registered_admission_schema_bound"] = True
    dependency["registered_execution_qualification_check_pass"] = True
    dependency["registered_raw_oracle_pass"] = True
    dependency["complete_foundational_raw_oracle_matrix_pass"] = True
    dependency["measured_case_count"] = bell.Q043_REQUIRED_CASE_COUNT
    return dependency


def _write_registered_manifest_fixture(
    root: Path,
    case_root: Path,
    member: dict[str, object],
) -> tuple[dict[str, object], Path]:
    submission_id = "11111111-1111-4111-8111-111111111111"
    freeze_id = "22222222-2222-4222-8222-222222222222"
    candidate_root = root / "clean_candidates" / freeze_id
    candidate_root.mkdir(parents=True)
    candidate_executable = candidate_root / "athena"
    candidate_executable.write_bytes(b"registered executable fixture\n")
    executable_sha256 = _sha256(candidate_executable)
    source_archive = candidate_root / "source.tar"
    with tarfile.open(source_archive, "w") as archive:
        for relative in (
            bell.SOURCE_PATH.as_posix(),
            bell.CORRECTED_EIGENMODE_HEADER_PATH.as_posix(),
        ):
            payload = (REPO_ROOT / relative).read_bytes()
            member_info = tarfile.TarInfo(relative)
            member_info.size = len(payload)
            member_info.mode = 0o644
            archive.addfile(member_info, io.BytesIO(payload))
    source_archive_sha256 = _sha256(source_archive)
    source_bundle_sha256 = "c" * 64
    candidate_manifest = candidate_root / "clean_candidate_manifest.json"
    candidate_manifest.write_text(
        json.dumps(
            {
                "schema_version": 4,
                "freeze_id": freeze_id,
                "source": {
                    "archive_path": str(source_archive),
                    "archive_sha256": source_archive_sha256,
                    "source_bundle_sha256": source_bundle_sha256,
                    "git_commit": "a" * 40,
                    "worktree_status": "clean",
                },
                "build": {
                    "source_archive_sha256": source_archive_sha256,
                    "source_bundle_sha256": source_bundle_sha256,
                    "executable_path": str(candidate_executable),
                    "executable_sha256": executable_sha256,
                },
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    candidate_manifest.chmod(0o444)
    source_archive.chmod(0o444)
    candidate_executable.chmod(0o555)
    manifest_root = (
        root
        / "manifests"
        / bell.REGISTERED_CAMPAIGN
        / submission_id
    )
    snapshot = manifest_root / "snapshot"
    snapshot.mkdir(parents=True)
    executable = snapshot / "athena"
    executable.write_bytes(candidate_executable.read_bytes())
    executable.chmod(0o555)
    receipt = {
        "submission_id": submission_id,
        "source_commit": "a" * 40,
        "control_plane_version": "b" * 64,
        "registered_science_authorization_id": "q023-linear-001-v1",
        "clean_candidate_manifest_sha256": _sha256(candidate_manifest),
        "executable_sha256": executable_sha256,
        "source_bundle_sha256": source_bundle_sha256,
        "source_archive_sha256": source_archive_sha256,
    }
    manifest = {
        "schema_version": bell.SCHEMA_VERSION,
        "pic_root": str(root),
        "campaign": bell.REGISTERED_CAMPAIGN,
        "test_id": member["member_id"],
        "submission_id": submission_id,
        "submission_scope": "registered_science",
        "artifact_dir": str(case_root),
        "git_commit": receipt["source_commit"],
        "control_plane_version": receipt["control_plane_version"],
        "registered_science_authorization_id": receipt[
            "registered_science_authorization_id"
        ],
        "clean_candidate_manifest_path": str(candidate_manifest),
        "clean_candidate_manifest_sha256": receipt[
            "clean_candidate_manifest_sha256"
        ],
        "snapshot_files": [
            {
                "role": "executable",
                "path": str(executable),
                "sha256": executable_sha256,
                "source_path": str(candidate_executable),
                "source_sha256": executable_sha256,
            }
        ],
    }
    manifest_path = manifest_root / "pre_submit_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_path.chmod(0o444)
    receipt.update(
        {
            "pre_submit_manifest_path": str(manifest_path),
            "pre_submit_manifest_sha256": _sha256(manifest_path),
        }
    )
    return receipt, executable


class Q023PaperBellLinearJOverCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary_directory = tempfile.TemporaryDirectory()
        cls.harness_binary = (
            Path(cls._temporary_directory.name) / "q023-paper-bell-linear-joverc-host"
        )
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(HARNESS),
                "-o",
                str(cls.harness_binary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        cls.harness_lines = subprocess.run(
            [str(cls.harness_binary)],
            cwd=REPO_ROOT,
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        ).stdout.splitlines()

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary_directory.cleanup()

    def test_source_helpers_enforce_integral_ppc_general_charge_and_j_over_c(self) -> None:
        ppc = [line.split() for line in self.harness_lines if line.startswith("ppc ")]
        self.assertEqual(
            [(float(row[1]), row[2]) for row in ppc],
            [(0.5, "0"), (1.0, "1"), (4.0, "1"), (4.5, "0")],
        )
        closures = [
            line.split() for line in self.harness_lines if line.startswith("closure ")
        ]
        self.assertEqual(len(closures), 3)
        for row in closures:
            self.assertAlmostEqual(float(row[3]), 4.0 * math.pi, places=13)
            self.assertEqual(row[4], "1")
        species = next(
            line.split() for line in self.harness_lines if line.startswith("species ")
        )
        self.assertEqual(species[1:], ["1", "0"])
        velocity_match = next(
            line.split()
            for line in self.harness_lines
            if line.startswith("velocity_match ")
        )
        self.assertEqual(velocity_match[1:], ["1", "0", "0"])
        unstable = next(
            line.split()
            for line in self.harness_lines
            if line.startswith("unstable_phase_zero ")
        )
        historical = next(
            line.split()
            for line in self.harness_lines
            if line.startswith("historical_phase_zero ")
        )
        self.assertAlmostEqual(float(unstable[1]), 1.0e-6)
        self.assertEqual(float(unstable[2]), 0.0)
        self.assertAlmostEqual(float(unstable[3]), -0.4e-6)
        self.assertLess(float(unstable[4]), 0.0)
        self.assertGreater(float(historical[4]), 0.0)

    def test_source_is_unique_registered_and_strictly_downstream_of_q043(self) -> None:
        source = (REPO_ROOT / bell.SOURCE_PATH).read_text(encoding="utf-8")
        self.assertIn("HasPositiveIntegralPPC(ppc)", source)
        self.assertIn("PPC must be a positive integer", source)
        self.assertIn('pin->GetReal("species0", "vx0")', source)
        self.assertIn('pin->GetReal("species0", "vy0")', source)
        self.assertIn('pin->GetReal("species0", "vz0")', source)
        self.assertIn("HasMatchingVelocityVectors", source)
        self.assertIn("actual particle initialization", source)
        self.assertIn(f'"{bell.PARTICLE_VELOCITY_SEMANTICS}"', source)
        self.assertIn('"Q043-BELL-CURRENT-VOLUME-AWARE"', source)
        self.assertIn('"Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE"', source)
        self.assertIn(f'"{bell.Q043_FOUNDATIONAL_BINDING_STATUS}"', source)
        self.assertIn(f'"{bell.Q043_FOUNDATIONAL_LINEAGE_COMMIT}"', source)
        self.assertIn(f'"{bell.Q043_CURRENT_SOURCE_SHA256}"', source)
        self.assertIn(f'"{bell.Q043_RAW_ORACLE_ANALYZER_SHA256}"', source)
        self.assertIn(f'"{bell.Q043_RAW_ORACLE_DECK_MANIFEST_SHA256}"', source)
        self.assertIn(f'"{bell.Q043_SUPERSESSION_SHA256}"', source)
        self.assertIn(f'"{bell.Q043_REGISTERED_ADMISSION_BINDING_STATUS}"', source)
        self.assertIn('"false_until_hardened_successor_lands"', source)
        self.assertNotIn("0bab5b1895bded3685be15d3237094b8a613c2007911818651fa070958ee7439", source)
        self.assertNotIn("8d901a65b58577dd3348b0569a43011f092721c0ee23c6270e27f80a31cbe370", source)
        self.assertNotIn("95a27328efb76e5e2d97b1d1c183abfb54c617df", source)
        self.assertIn('"132"', source)
        self.assertIn('"required_before_linear_qualification"', source)
        self.assertIn('"required_before_any_q023_execution_or_qualification"', source)
        self.assertIn('"dimension_appropriate_multidirectional_matrix_required"', source)
        self.assertIn('"corrected_linear_eigenmode"', source)
        for name in (
            "couple_j_to_efield_coeff",
            "couple_moments_momentum_coeff",
            "couple_moments_energy_coeff",
        ):
            self.assertIn(f'pin->GetReal("particles", "{name}")', source)
        self.assertNotIn('"uniform_current_oracle"', source)
        self.assertNotIn("2.0*b_g*light_speed*k0", source)
        pgen = (REPO_ROOT / "src/pgen/pgen.cpp").read_text(encoding="utf-8")
        self.assertEqual(pgen.count('compare("q023_paper_bell_linear_joverc")'), 2)
        self.assertIn(
            "void Q023PaperBellLinearJOverC(ParameterInput *pin, const bool restart);",
            (REPO_ROOT / "src/pgen/pgen.hpp").read_text(encoding="utf-8"),
        )
        self.assertIn(
            "pgen/tests/q023_paper_bell_linear_joverc.cpp",
            (REPO_ROOT / "src/CMakeLists.txt").read_text(encoding="utf-8"),
        )

    def test_every_deck_matches_the_actual_source_literal_contract(self) -> None:
        source = (REPO_ROOT / bell.SOURCE_PATH).read_text(encoding="utf-8")
        required_strings = {
            (bell.PGEN_NAME if block == "block" else block.strip('"'), name): expected
            for block, name, expected in re.findall(
                r'Q023JOverCRequireString\(\s*pin,\s*(block|"[^"]+")\s*,\s*'
                r'"([^"]+)"\s*,\s*"([^"]+)"\s*\);',
                source,
                re.DOTALL,
            )
        }
        required_booleans = {
            (bell.PGEN_NAME if block == "block" else block.strip('"'), name): expected
            for block, name, expected in re.findall(
                r'Q023JOverCRequireBoolean\(\s*pin,\s*(block|"[^"]+")\s*,\s*'
                r'"([^"]+)"\s*,\s*(true|false)\s*\);',
                source,
                re.DOTALL,
            )
        }
        self.assertGreater(len(required_strings), 20)
        self.assertGreater(len(required_booleans), 5)
        manifest = bell.validate_checked_in_decks()
        for record in manifest["cases"]:
            text = (bell.DECK_ROOT / record["deck_path"]).read_text(encoding="utf-8")
            self.assertNotIn("95a27328efb76e5e2d97b1d1c183abfb54c617df", text)
            blocks = bell.parse_athinput_text(text)
            for (block, name), expected in required_strings.items():
                self.assertEqual(blocks[block][name], expected)
            for (block, name), expected in required_booleans.items():
                self.assertEqual(blocks[block][name], expected)

    def test_checked_in_deck_matrix_binds_volume_aware_multidirectional_contract(
        self,
    ) -> None:
        manifest = bell.validate_checked_in_decks()
        self.assertEqual(manifest["case_count"], 55)
        self.assertEqual(
            manifest["foundational_binding_status"],
            bell.Q043_FOUNDATIONAL_BINDING_STATUS,
        )
        self.assertEqual(
            manifest["foundational_current_lineage_commit"],
            bell.Q043_FOUNDATIONAL_LINEAGE_COMMIT,
        )
        self.assertEqual(
            manifest["foundational_registered_admission_binding_status"],
            bell.Q043_REGISTERED_ADMISSION_BINDING_STATUS,
        )
        self.assertFalse(manifest["foundational_registered_admission_digest_bound"])
        self.assertFalse(manifest["foundational_registered_admission_schema_bound"])
        self.assertEqual(
            manifest["foundational_raw_oracle_case_count"],
            bell.Q043_REQUIRED_CASE_COUNT,
        )
        self.assertFalse(manifest["launch_authorized"])
        self.assertFalse(manifest["qualification_eligible"])
        self.assertEqual(
            manifest["decompositions_by_dimension"]["2"]["split_x1x2"], [2, 2, 1]
        )
        self.assertEqual(
            manifest["decompositions_by_dimension"]["3"]["split_xyz"], [2, 2, 2]
        )
        self.assertEqual(
            manifest["decompositions_by_dimension"]["3"]["split_x2"], [1, 2, 1]
        )
        self.assertEqual(
            {record["ppc"] for record in manifest["cases"]},
            {1},
        )
        self.assertEqual(
            {record["species_mass"] for record in manifest["cases"]},
            {2.0},
        )
        for record in manifest["cases"]:
            self.assertAlmostEqual(
                record["configured_deposited_j_over_c"], 4.0 * math.pi
            )
            self.assertEqual(
                _sha256(bell.DECK_ROOT / record["deck_path"]), record["deck_sha256"]
            )
            blocks = bell.parse_athinput_text(
                (bell.DECK_ROOT / record["deck_path"]).read_text(encoding="utf-8")
            )
            self.assertEqual(
                blocks[bell.PGEN_NAME]["foundational_raw_oracle_id"],
                bell.Q043_RAW_ORACLE_ID,
            )
            self.assertEqual(
                blocks[bell.PGEN_NAME]["foundational_binding_status"],
                bell.Q043_FOUNDATIONAL_BINDING_STATUS,
            )
            self.assertEqual(
                blocks[bell.PGEN_NAME]["foundational_raw_oracle_case_count"],
                str(bell.Q043_REQUIRED_CASE_COUNT),
            )
            self.assertAlmostEqual(record["stream_speed"], 1.0 / record["epsilon"])
            self.assertAlmostEqual(
                record["actual_species_stream_speed"], record["stream_speed"]
            )
            for axis in ("x", "y", "z"):
                self.assertEqual(
                    blocks["species0"][f"v{axis}0"],
                    blocks["particles"][f"cr_v{axis}0"],
                )
            self.assertEqual(blocks[bell.PGEN_NAME]["launch_authorized"], "false")
            self.assertEqual(blocks[bell.PGEN_NAME]["qualification_eligible"], "false")
            self.assertEqual(float(blocks["time"]["tlim"]), bell.LINEAR_RUNTIME_TLIM)
            for index in range(1, 6):
                self.assertEqual(float(blocks[f"output{index}"]["dt"]), bell.LINEAR_OUTPUT_DT)

    def test_registered_manifest_rebinding_changes_only_dependency_metadata(
        self,
    ) -> None:
        dependency = bell.synthetic_q043_registered_raw_oracle_dependency()
        dependency.update(
            {
                "binding_kind": "registered_matrix_qualification",
                "registered_admission_binding_status": (
                    bell.Q043_REGISTERED_MATRIX_BINDING_STATUS
                ),
                "registered_admission_digest_bound": True,
                "registered_admission_schema_bound": True,
                "registered_execution_qualification_check_pass": True,
                "registered_raw_oracle_pass": True,
                "complete_foundational_raw_oracle_matrix_pass": True,
                "registered_matrix_path": "analysis/q043/matrix.json",
                "registered_matrix_sha256": "a" * 64,
                "registered_matrix_record_type": (
                    "q043_registered_execution_raw_oracle_matrix_qualification"
                ),
                "registered_matrix_case_bindings_sha256": "b" * 64,
                "measured_case_count": bell.Q043_REQUIRED_CASE_COUNT,
            }
        )
        pending, pending_decks = bell.build_deck_manifest()
        with patch.object(
            bell, "validate_q043_dependency", return_value=dependency
        ):
            registered, registered_decks = bell.build_deck_manifest(
                q043_registered_raw_oracle_dependency=dependency,
                q043_artifact_root=Path("/registered/q043"),
            )
        self.assertEqual(registered_decks, pending_decks)
        self.assertEqual(registered["cases"], pending["cases"])
        self.assertEqual(
            registered["foundational_registered_admission_binding_status"],
            bell.Q043_REGISTERED_MATRIX_BINDING_STATUS,
        )
        self.assertTrue(
            registered["foundational_registered_admission_digest_bound"]
        )
        self.assertEqual(
            registered["foundational_registered_matrix_sha256"], "a" * 64
        )
        self.assertEqual(
            registered["foundational_registered_dependency_sha256"],
            bell._dependency_digest(dependency),
        )

    def test_fractional_ppc_and_source_string_drift_fail_closed(self) -> None:
        member = bell.expected_deck_members()[0]
        text = bell.render_deck(member)
        with self.assertRaisesRegex(bell.ContractError, "positive integer"):
            bell.validate_rendered_deck(
                member, text.replace("ppc = 1.0", "ppc = 1.5", 1)
            )
        with self.assertRaisesRegex(bell.ContractError, "foundational_raw_oracle_id"):
            bell.validate_rendered_deck(
                member,
                text.replace(
                    f"foundational_raw_oracle_id = {bell.Q043_RAW_ORACLE_ID}",
                    "foundational_raw_oracle_id = stale",
                    1,
                ),
            )
        for name in (
            "couple_j_to_efield_coeff",
            "couple_moments_momentum_coeff",
            "couple_moments_energy_coeff",
        ):
            with self.subTest(name=name):
                with self.assertRaisesRegex(bell.ContractError, name):
                    bell.validate_rendered_deck(
                        member,
                        text.replace(f"{name} = 1.0", f"{name} = 0.0", 1),
                    )

    def test_species_specific_zero_and_different_velocity_overrides_fail_closed(
        self,
    ) -> None:
        member = next(
            value
            for value in bell.expected_deck_members()
            if value["member_id"] == "d1-coarse-reference_x1-eps0p4"
        )
        text = bell.render_deck(member)
        species_block = text.split("<species0>", 1)[1].split("<problem>", 1)[0]
        expected_vx = bell._float_token(float(member["stream_speed"]))
        for label, replacement in (("zero", "0.0"), ("different", "2.0")):
            with self.subTest(label=label):
                mutated_species = species_block.replace(
                    f"vx0 = {expected_vx}", f"vx0 = {replacement}", 1
                )
                with self.assertRaisesRegex(
                    bell.ContractError, "species0 actual initialization velocity drifted"
                ):
                    bell.validate_rendered_deck(
                        member, text.replace(species_block, mutated_species, 1)
                    )
        with self.assertRaisesRegex(
            bell.ContractError,
            "global CR velocity and species0 actual initialization velocity differ",
        ):
            bell.validate_rendered_deck(
                member,
                text.replace(f"cr_vx0 = {expected_vx}", "cr_vx0 = 2.0", 1),
            )

    def test_complete_synthetic_matrix_passes_all_gates_without_authority(self) -> None:
        bundle = bell.synthetic_predecessor_bundle()
        report = bell.analyze_predecessor_bundle(bundle)
        self.assertEqual(report["record_count"], 55)
        self.assertTrue(report["q043_dependency_contract_pass"])
        self.assertFalse(report["q043_registered_raw_oracle_dependency_pass"])
        self.assertFalse(report["registered_q043_foundation_bound"])
        self.assertFalse(report["materialized_q023_trace_matrix_bound"])
        self.assertTrue(report["growth_gates_pass"])
        self.assertTrue(report["signed_phase_gates_pass"])
        self.assertTrue(report["polarization_gates_pass"])
        self.assertTrue(report["growth_fit_quality_gates_pass"])
        self.assertTrue(report["phase_fit_quality_gates_pass"])
        self.assertTrue(report["provenance_gates_pass"])
        self.assertTrue(report["resolution_convergence_gate_pass"])
        self.assertTrue(report["multidirectional_decomposition_gate_pass"])
        self.assertTrue(report["predecessor_contract_pass"])
        for record in report["records"]:
            self.assertAlmostEqual(
                record["expected_signed_phase_frequency_over_k0_ua"],
                -record["epsilon"],
            )
        self.assertFalse(report["linear_qualification_prerequisites_pass"])
        self.assertFalse(report["passed"])
        self.assertFalse(report["launch_authorized"])
        self.assertFalse(report["qualification_eligible"])
        self.assertFalse(report["scientific_claim_authorized"])
        self.assertFalse(report["publication_authorized"])
        self.assertFalse(report["complete_authoritative_q043_foundation_claimed"])
        physics = bell.analyze_physics_trace_matrix(
            [
                {
                    key: record[key]
                    for key in bell._PHYSICS_MATRIX_RECORD_KEYS
                }
                for record in bundle["records"]
            ]
        )
        self.assertTrue(physics["predecessor_contract_pass"])
        self.assertFalse(physics["qualification_eligible"])
        self.assertFalse(physics["passed"])
        self.assertEqual(
            physics["records"],
            [
                {
                    key: record[key]
                    for key in record
                    if key not in {"provenance_kind", "provenance_gate_pass"}
                }
                for record in report["records"]
            ],
        )

    def test_pure_physics_matrix_rejects_noncanonical_order(self) -> None:
        records = [
            {
                key: record[key]
                for key in bell._PHYSICS_MATRIX_RECORD_KEYS
            }
            for record in bell.synthetic_predecessor_bundle()["records"]
        ]
        records[0], records[1] = records[1], records[0]
        with self.assertRaisesRegex(
            bell.ContractError, "incomplete, extra, or noncanonical"
        ):
            bell.analyze_physics_trace_matrix(records)

    def test_growth_signed_phase_and_polarization_fail_scientifically(self) -> None:
        mutations = (
            ("growth", {"growth": bell.theoretical_dispersion(0.4)[1] + 0.04}),
            ("signed_phase", {"phase": 0.4}),
            ("polarization", {"left_ratio": 1.0}),
        )
        for label, values in mutations:
            with self.subTest(label=label):
                bundle = bell.synthetic_predecessor_bundle()
                _set_physics_trace(
                    _record(bundle, "d1-fine-reference_x1", 0.4), **values
                )
                report = bell.analyze_predecessor_bundle(bundle)
                self.assertFalse(report["predecessor_contract_pass"])
                if label == "growth":
                    self.assertFalse(report["growth_gates_pass"])
                elif label == "signed_phase":
                    self.assertFalse(report["signed_phase_gates_pass"])
                else:
                    self.assertFalse(report["polarization_gates_pass"])
                self.assertFalse(report["passed"])

    def test_robust_growth_phase_quality_and_fit_window_polarization_gates(self) -> None:
        bundle = bell.synthetic_predecessor_bundle()
        record = _record(bundle, "d1-fine-reference_x1", 0.4)
        trace = record["physics_trace"]
        modulation = np.exp(0.25 * np.sin(np.linspace(0.0, 12.0 * math.pi, 49)))
        for name in (
            "right_mode_real",
            "right_mode_imag",
            "left_mode_real",
            "left_mode_imag",
            "velocity_right_mode_real",
            "velocity_right_mode_imag",
            "velocity_left_mode_real",
            "velocity_left_mode_imag",
            "paper_delta_u_y_sine_fit_real",
            "paper_delta_u_y_sine_fit_imag",
            "paper_volume_averaged_abs_delta_u",
        ):
            trace[name] = (np.asarray(trace[name]) * modulation).tolist()
        report = bell.analyze_predecessor_bundle(bundle)
        target = next(item for item in report["records"] if item["member_id"] == record["member_id"])
        self.assertLess(target["growth_fit_r2"], bell.MIN_GROWTH_FIT_R2)
        self.assertFalse(report["growth_gates_pass"])
        self.assertFalse(report["growth_fit_quality_gates_pass"])

        bundle = bell.synthetic_predecessor_bundle()
        record = _record(bundle, "d1-fine-reference_x1", 0.4)
        trace = record["physics_trace"]
        time = np.asarray(trace["normalized_time"])
        phase_modulation = 0.6 * np.sin(np.linspace(0.0, 10.0 * math.pi, time.size))
        for real_name, imag_name in (
            ("right_mode_real", "right_mode_imag"),
            ("left_mode_real", "left_mode_imag"),
            ("velocity_right_mode_real", "velocity_right_mode_imag"),
            ("velocity_left_mode_real", "velocity_left_mode_imag"),
            ("paper_delta_u_y_sine_fit_real", "paper_delta_u_y_sine_fit_imag"),
        ):
            mode = np.asarray(trace[real_name]) + 1.0j * np.asarray(trace[imag_name])
            mode *= np.exp(1.0j * phase_modulation)
            trace[real_name] = mode.real.tolist()
            trace[imag_name] = mode.imag.tolist()
        for mode_names, interval_name, change_name in (
            (
                ("right_mode_real", "right_mode_imag"),
                "phase_interval",
                "phase_change",
            ),
            (
                ("velocity_right_mode_real", "velocity_right_mode_imag"),
                "velocity_phase_interval",
                "velocity_phase_change",
            ),
            (
                ("paper_delta_u_y_sine_fit_real", "paper_delta_u_y_sine_fit_imag"),
                "paper_delta_u_y_phase_interval",
                "paper_delta_u_y_phase_change",
            ),
        ):
            mode = np.asarray(trace[mode_names[0]]) + 1.0j * np.asarray(trace[mode_names[1]])
            interval, change = legacy._fixed_interval_phase_trace(time, mode)
            trace[interval_name] = interval.tolist()
            trace[change_name] = change.tolist()
        report = bell.analyze_predecessor_bundle(bundle)
        target = next(item for item in report["records"] if item["member_id"] == record["member_id"])
        self.assertLess(target["phase_fit_r2"], bell.MIN_PHASE_FIT_R2)
        self.assertFalse(report["signed_phase_gates_pass"])
        self.assertFalse(report["phase_fit_quality_gates_pass"])

        bundle = bell.synthetic_predecessor_bundle()
        record = _record(bundle, "d1-fine-reference_x1", 0.4)
        trace = record["physics_trace"]
        _, growth = bell.theoretical_dispersion(0.4)
        time = np.asarray(trace["normalized_time"])
        fit_indices = np.flatnonzero((growth * time >= 1.0) & (growth * time <= 5.0))
        index = int(fit_indices[len(fit_indices) // 2])
        trace["left_mode_real"][index] = trace["right_mode_real"][index]
        trace["left_mode_imag"][index] = trace["right_mode_imag"][index]
        report = bell.analyze_predecessor_bundle(bundle)
        self.assertFalse(report["polarization_gates_pass"])

    def test_resolution_convergence_and_multidirectional_decomposition_fail_closed(
        self,
    ) -> None:
        bundle = bell.synthetic_predecessor_bundle()
        expected_growth = bell.theoretical_dispersion(0.4)[1]
        _set_physics_trace(
            _record(bundle, "d1-coarse-reference_x1", 0.4),
            growth=expected_growth + 0.019,
        )
        _set_physics_trace(
            _record(bundle, "d1-fine-reference_x1", 0.4),
            growth=expected_growth - 0.019,
        )
        report = bell.analyze_predecessor_bundle(bundle)
        self.assertTrue(report["growth_gates_pass"])
        self.assertFalse(report["resolution_convergence_gate_pass"])
        self.assertFalse(report["predecessor_contract_pass"])

        bundle = bell.synthetic_predecessor_bundle()
        _set_physics_trace(
            _record(bundle, "d3-fine-split_xyz", 0.4),
            growth=expected_growth + 0.019,
        )
        report = bell.analyze_predecessor_bundle(bundle)
        self.assertTrue(report["growth_gates_pass"])
        self.assertFalse(report["multidirectional_decomposition_gate_pass"])
        self.assertFalse(report["predecessor_contract_pass"])

    def test_matrix_q043_dependency_and_provenance_drift_fail_closed(self) -> None:
        bundle = bell.synthetic_predecessor_bundle()
        bundle["records"].pop()
        with self.assertRaisesRegex(bell.ContractError, "matrix is incomplete"):
            bell.analyze_predecessor_bundle(bundle)

        bundle = bell.synthetic_predecessor_bundle()
        bundle["q043_registered_raw_oracle_dependency"][
            "complete_foundational_raw_oracle_matrix_pass"
        ] = True
        with self.assertRaisesRegex(bell.ContractError, "must remain unclaimed"):
            bell.analyze_predecessor_bundle(bundle)

        bundle = bell.synthetic_predecessor_bundle()
        bundle["q043_registered_raw_oracle_dependency"]["integration_checkpoint_commit"] = (
            "0" * 39
        )
        with self.assertRaisesRegex(bell.ContractError, "checkpoint commit malformed"):
            bell.analyze_predecessor_bundle(bundle)

        bundle = bell.synthetic_predecessor_bundle()
        bundle["q043_registered_raw_oracle_dependency"]["registered_raw_oracle_pass"] = (
            True
        )
        with self.assertRaisesRegex(bell.ContractError, "must remain unclaimed"):
            bell.analyze_predecessor_bundle(bundle)

        bundle = bell.synthetic_predecessor_bundle()
        bundle["q043_registered_raw_oracle_dependency"][
            "registered_admission_digest_bound"
        ] = True
        with self.assertRaisesRegex(bell.ContractError, "provisional.*binding is forbidden"):
            bell.analyze_predecessor_bundle(bundle)

        bundle = bell.synthetic_predecessor_bundle()
        bundle["records"][0]["provenance"]["deck_sha256"] = "0" * 64
        with self.assertRaisesRegex(bell.ContractError, "deck binding drifted"):
            bell.analyze_predecessor_bundle(bundle)

    def test_registered_q043_dependency_binds_exact_hardened_matrix(self) -> None:
        matrix = {
            "record_type": bell.q043_registered.MATRIX_RECORD_TYPE,
            "case_count": bell.Q043_REQUIRED_CASE_COUNT,
            "case_bindings_sha256": "a" * 64,
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            path = root / "q043-matrix.json"
            path.write_text(
                json.dumps(matrix, sort_keys=True) + "\n", encoding="utf-8"
            )
            with patch.object(
                bell.q043_registered,
                "validate_downstream_q023_q019_prerequisite",
                return_value=matrix,
            ):
                dependency = bell.registered_q043_raw_oracle_dependency(
                    path, artifact_root=root
                )
                self.assertEqual(
                    dependency["binding_kind"], "registered_matrix_qualification"
                )
                self.assertTrue(
                    dependency["complete_foundational_raw_oracle_matrix_pass"]
                )
                self.assertEqual(
                    dependency["registered_matrix_case_bindings_sha256"],
                    matrix["case_bindings_sha256"],
                )
                bell.validate_q043_dependency(
                    dependency, artifact_root=root
                )

                path.write_text("{}\n", encoding="utf-8")
                with self.assertRaisesRegex(bell.ContractError, "digest drifted"):
                    bell.validate_q043_dependency(
                        dependency, artifact_root=root
                    )

    def test_materialized_provenance_requires_clean_bound_files_and_registered_q043(
        self,
    ) -> None:
        dependency = _q043_bound_dependency_fixture_for_q023_provenance_only()
        member = next(
            value
            for value in bell._manifest_members().values()
            if value["member_id"] == "d1-coarse-reference_x1-eps0p4"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            provenance, physics_trace = _write_materialized_member(
                root, member, dependency
            )
            self.assertEqual(
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                ),
                "registered_execution_trace",
            )
            arbitrary_trace = copy.deepcopy(physics_trace)
            arbitrary_trace["right_mode_real"][10] += 1.0e-8
            with self.assertRaisesRegex(
                bell.ContractError, "was not derived from the exact retained raw outputs"
            ):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=arbitrary_trace,
                )
            raw_path = root / provenance["raw_artifacts"][10]["path"]
            raw_path.write_bytes(b"unrelated raw file\n")
            with self.assertRaisesRegex(bell.ContractError, "digest drifted"):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                )

    def test_hardened_registered_receipt_binds_controller_and_raw_inventory(
        self,
    ) -> None:
        dependency = _q043_bound_dependency_fixture_for_q023_provenance_only()
        member = next(
            value
            for value in bell._manifest_members().values()
            if value["member_id"] == "d1-coarse-reference_x1-eps0p4"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            provenance, _ = _write_materialized_member(root, member, dependency)
            rank_count = math.prod(
                int(value) for value in member["decomposition_splits"]
            )
            stdout = root / "athena_stdout.txt"
            stdout.write_text(
                (
                    f"Q023_REGISTERED_EXECUTION case_id={member['member_id']} "
                    f"mpi_world_size={rank_count} "
                    f"rank_ids={','.join(str(rank) for rank in range(rank_count))}\n"
                    "Q023_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0\n"
                    "Terminating on time limit\n"
                    "time=1.750000e+00 cycle=100\n"
                    "tlim=1.750000e+00 nlim=2000000\n"
                ),
                encoding="utf-8",
            )
            raw_inventory = [
                {
                    "path": str(artifact["path"]).removeprefix("raw/"),
                    "sha256": artifact["sha256"],
                    "byte_count": (root / artifact["path"]).stat().st_size,
                    "member_id": member["member_id"],
                    "variable": artifact["variable"],
                    "output_index": index // len(bell._RAW_OUTPUT_VARIABLES),
                }
                for index, artifact in enumerate(provenance["raw_artifacts"])
            ]
            producer = {
                "entrypoint": "reconcile_q023_registered_execution.py",
                "entrypoint_sha256": "1" * 64,
                "launch_trampoline_sha256": "2" * 64,
                "control_plane_version": "3" * 64,
            }
            action = {
                "action_id": member["member_id"],
                "kind": "athena",
                "resources": {"tasks": rank_count},
                "arguments": [
                    {"literal": "-i"},
                    {"snapshot_role": "input-deck"},
                    {"literal": "-d"},
                    {"artifact_directory": "raw"},
                ],
            }
            wrapper = {
                "stdout_sha256": _sha256(stdout),
                "observed_world_size": rank_count,
                "observed_rank_ids": list(range(rank_count)),
                "terminal_cycle": 88,
                "terminal_time": bell.LINEAR_RUNTIME_TLIM,
            }
            receipt = {
                "schema_version": bell.SCHEMA_VERSION,
                "record_type": "q023_reconciled_registered_execution_receipt",
                "receipt_role": "immutable_reconciled_registered_execution",
                "registration_scope": "registered_science",
                "reconciled": True,
                "campaign_id": bell.CAMPAIGN_ID,
                "member_id": member["member_id"],
                "reservation_id": "r",
                "submission_id": "s",
                "reconciliation_event_sha256": "4" * 64,
                "reconciliation_mirror_ack_sha256": "5" * 64,
                "control_plane_version": producer["control_plane_version"],
                "project_home_mirrors": {},
                "producer": producer,
                "registered_science_authorization_id": "q023-linear-001",
                "source_commit": "a" * 40,
                "source_bundle_sha256": "6" * 64,
                "source_archive_sha256": "7" * 64,
                "clean_candidate_manifest_sha256": "8" * 64,
                "executable_sha256": provenance["executable_sha256"],
                "environment_sha256": "9" * 64,
                "deck_sha256": member["deck_sha256"],
                "command_evidence": {
                    "source": "trusted_pre_submit_manifest_and_installed_trampoline",
                    "executor": "trusted_trampoline_athena_argv_v1",
                    "action": action,
                    "launch_trampoline_entrypoint": "launch_trampoline.py",
                    "launch_trampoline_sha256": producer[
                        "launch_trampoline_sha256"
                    ],
                    "trusted_wrapper_evidence": wrapper,
                },
                "mpi_evidence": {
                    "tasks": rank_count,
                    "observed_world_size": rank_count,
                    "observed_rank_ids": list(range(rank_count)),
                },
                "slurm_job_id": "123",
                "slurm_terminal_state": "COMPLETED",
                "slurm_exit_code": "0:0",
                "terminal_cycle": 88,
                "terminal_time": bell.LINEAR_RUNTIME_TLIM,
                "raw_output_root": str(root / "raw"),
                "artifact_dir": str(root),
                "artifact_inventory": {},
                "trampoline_completion_receipt": {},
                "terminal_receipt_sha256": "a" * 64,
                "pre_submit_manifest_path": str(root / "manifest.json"),
                "pre_submit_manifest_sha256": "b" * 64,
                "raw_inventory": raw_inventory,
                "raw_inventory_sha256": bell._sha256_bytes(
                    bell._canonical_json_bytes(raw_inventory)
                ),
            }
            bell._validate_registered_execution_receipt(
                receipt,
                root=root,
                member=member,
                provenance=provenance,
                raw_artifacts=provenance["raw_artifacts"],
            )
            receipt["raw_inventory"][0]["sha256"] = "0" * 64
            receipt["raw_inventory_sha256"] = bell._sha256_bytes(
                bell._canonical_json_bytes(receipt["raw_inventory"])
            )
            with self.assertRaisesRegex(
                bell.ContractError, "raw inventory differs from provenance"
            ):
                bell._validate_registered_execution_receipt(
                    receipt,
                    root=root,
                    member=member,
                    provenance=provenance,
                    raw_artifacts=provenance["raw_artifacts"],
                )

    def test_registered_manifest_binds_only_canonical_immutable_executable(
        self,
    ) -> None:
        member = next(iter(bell._manifest_members().values()))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            case_root = root / "runs/case"
            case_root.mkdir(parents=True)
            receipt, executable = _write_registered_manifest_fixture(
                root, case_root, member
            )
            binding = bell.registered_manifest_executable_binding(
                receipt,
                member=member,
                artifact_root=case_root,
                authorized_manifest_root=root,
            )
            self.assertEqual(binding["path"], str(executable))
            self.assertEqual(binding["sha256"], _sha256(executable))
            self.assertEqual(
                binding["source_bindings"][bell.SOURCE_PATH.as_posix()][
                    "sha256"
                ],
                _sha256(REPO_ROOT / bell.SOURCE_PATH),
            )

            traversing_receipt = copy.deepcopy(receipt)
            traversing_receipt["submission_id"] = "../escape"
            with self.assertRaisesRegex(
                bell.ContractError, "submission ID is not a canonical UUID"
            ):
                bell.registered_manifest_executable_binding(
                    traversing_receipt,
                    member=member,
                    artifact_root=case_root,
                    authorized_manifest_root=root,
                )

            manifest_path = Path(str(receipt["pre_submit_manifest_path"]))
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            alternate = manifest_path.parent / "snapshot/alternate-athena"
            alternate.write_bytes(executable.read_bytes())
            alternate.chmod(0o555)
            manifest["snapshot_files"][0]["path"] = str(alternate)
            manifest_path.chmod(0o644)
            manifest_path.write_text(
                json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
            )
            manifest_path.chmod(0o444)
            receipt["pre_submit_manifest_sha256"] = _sha256(manifest_path)
            with self.assertRaisesRegex(
                bell.ContractError, "executable snapshot record"
            ):
                bell.registered_manifest_executable_binding(
                    receipt,
                    member=member,
                    artifact_root=case_root,
                    authorized_manifest_root=root,
                )

            manifest["snapshot_files"][0]["path"] = str(executable)
            manifest_path.chmod(0o644)
            manifest_path.write_text(
                json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
            )
            manifest_path.chmod(0o444)
            receipt["pre_submit_manifest_sha256"] = _sha256(manifest_path)
            os.link(executable, executable.parent / "hardlink-athena")
            with self.assertRaisesRegex(
                bell.ContractError, "immutable retained file"
            ):
                bell.registered_manifest_executable_binding(
                    receipt,
                    member=member,
                    artifact_root=case_root,
                    authorized_manifest_root=root,
                )

    def test_retained_raw_batch_rejects_namespace_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            raw_bin = root / "raw/bin"
            raw_bin.mkdir(parents=True)
            raw_path = raw_bin / "sample.bin"
            raw_path.write_bytes(b"receipt-bound raw bytes\n")
            raw_path.chmod(0o444)
            raw_bin.chmod(0o555)
            raw_bin.parent.chmod(0o555)
            root.chmod(0o555)
            binding = {
                "path": "raw/bin/sample.bin",
                "sha256": _sha256(raw_path),
                "byte_count": raw_path.stat().st_size,
            }
            batch = bell._RetainedRawBatch(root, [binding])
            try:
                raw_bin.chmod(0o755)
                displaced = raw_bin / "displaced.bin"
                raw_path.rename(displaced)
                raw_path.write_bytes(displaced.read_bytes())
                raw_path.chmod(0o444)
                raw_bin.chmod(0o555)
                with self.assertRaisesRegex(
                    bell.ContractError,
                    "retained raw namespace changed during analysis",
                ):
                    batch.revalidate()
            finally:
                batch.close()

    def test_materialized_receipt_and_actual_mpi_topology_fail_closed(self) -> None:
        dependency = _q043_bound_dependency_fixture_for_q023_provenance_only()
        member = next(
            value
            for value in bell._manifest_members().values()
            if value["member_id"] == "d1-coarse-reference_x1-eps0p4"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            provenance, physics_trace = _write_materialized_member(root, member, dependency)
            receipt_path = root / provenance["registered_execution_receipt_path"]
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["mpi"]["rank_count"] = 1
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            provenance["registered_execution_receipt_sha256"] = _sha256(receipt_path)
            with self.assertRaisesRegex(bell.ContractError, "MPI rank count"):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                )

    def test_materialized_raw_output_counter_progression_fails_closed(self) -> None:
        dependency = _q043_bound_dependency_fixture_for_q023_provenance_only()
        member = next(
            value
            for value in bell._manifest_members().values()
            if value["member_id"] == "d1-coarse-reference_x1-eps0p4"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            provenance, physics_trace = _write_materialized_member(root, member, dependency)
            raw_path = root / provenance["raw_artifacts"][0]["path"]
            raw_path.write_bytes(
                raw_path.read_bytes().replace(
                    b"file_number = 0\n", b"file_number = 9\n", 1
                )
            )
            provenance["raw_artifacts"][0]["sha256"] = _sha256(raw_path)
            receipt_path = root / provenance["registered_execution_receipt_path"]
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["raw_artifacts"] = provenance["raw_artifacts"]
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            provenance["registered_execution_receipt_sha256"] = _sha256(receipt_path)
            with self.assertRaisesRegex(bell.ContractError, "counter progression"):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                )

    def test_materialized_particle_current_inventory_and_physics_fail_closed(
        self,
    ) -> None:
        dependency = _q043_bound_dependency_fixture_for_q023_provenance_only()
        member = next(
            value
            for value in bell._manifest_members().values()
            if value["member_id"] == "d1-coarse-reference_x1-eps0p4"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            provenance, physics_trace = _write_materialized_member(
                root, member, dependency
            )
            provenance["raw_artifacts"] = [
                artifact
                for artifact in provenance["raw_artifacts"]
                if not (
                    artifact["variable"] == "prtcl_jz" and artifact["cycle"] == 1
                )
            ]
            receipt_path = root / provenance["registered_execution_receipt_path"]
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["raw_artifacts"] = provenance["raw_artifacts"]
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            provenance["registered_execution_receipt_sha256"] = _sha256(receipt_path)
            with self.assertRaisesRegex(bell.ContractError, "field inventory"):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                )

            provenance, physics_trace = _write_materialized_member(
                root, member, dependency
            )
            artifact = next(
                value
                for value in provenance["raw_artifacts"]
                if value["variable"] == "prtcl_jy" and value["cycle"] == 1
            )
            raw_path = root / artifact["path"]
            shape = tuple(
                reversed(tuple(int(value) for value in member["global_nx"]))
            )
            transverse = np.zeros(shape)
            transverse.flat[0] = 0.1
            raw_path.write_bytes(
                _raw_payload_from_fields(
                    member,
                    float(artifact["time"]),
                    int(artifact["cycle"]),
                    fields=("prtcl_jy",),
                    values={"prtcl_jy": transverse},
                    variable_output_index=bell._RAW_OUTPUT_INDEX["prtcl_jy"],
                )
            )
            artifact["sha256"] = _sha256(raw_path)
            receipt_path = root / provenance["registered_execution_receipt_path"]
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["raw_artifacts"] = provenance["raw_artifacts"]
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            provenance["registered_execution_receipt_sha256"] = _sha256(receipt_path)
            with self.assertRaisesRegex(
                bell.ContractError, "transverse deposited current"
            ):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                )

            provenance, physics_trace = _write_materialized_member(root, member, dependency)
            topology_path = root / "execution/rank_topology.json"
            topology = json.loads(topology_path.read_text(encoding="utf-8"))
            topology["ranks"][1]["meshblock_logical_locations"] = [[0, 0, 0]]
            topology_path.write_text(json.dumps(topology), encoding="utf-8")
            receipt_path = root / provenance["registered_execution_receipt_path"]
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["mpi"]["rank_topology_sha256"] = _sha256(topology_path)
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            provenance["registered_execution_receipt_sha256"] = _sha256(receipt_path)
            with self.assertRaisesRegex(bell.ContractError, "rank topology"):
                bell._validate_provenance(
                    provenance,
                    member=member,
                    dependency=dependency,
                    artifact_root=root,
                    physics_trace=physics_trace,
                )

    def test_provisional_q043_registered_admission_binding_is_rejected(self) -> None:
        dependency = bell.synthetic_q043_registered_raw_oracle_dependency()
        dependency["binding_kind"] = "registered_matrix_qualification"
        with self.assertRaisesRegex(
            bell.ContractError,
            "registered Q043 dependency requires an absolute artifact root",
        ):
            bell.validate_q043_dependency(dependency)

    def test_readiness_record_binds_exact_non_authorizing_predecessor(self) -> None:
        record = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(record["campaign_id"], bell.CAMPAIGN_ID)
        self.assertEqual(
            record["dependency_contract"]["foundational_raw_oracle_id"],
            bell.Q043_RAW_ORACLE_ID,
        )
        self.assertEqual(
            record["dependency_contract"]["foundational_binding_status"],
            bell.Q043_FOUNDATIONAL_BINDING_STATUS,
        )
        self.assertEqual(
            record["dependency_contract"]["foundational_current_lineage_commit"],
            bell.Q043_FOUNDATIONAL_LINEAGE_COMMIT,
        )
        self.assertEqual(
            record["dependency_contract"]["foundational_current_source_sha256"],
            bell.Q043_CURRENT_SOURCE_SHA256,
        )
        self.assertEqual(
            record["dependency_contract"]["foundational_raw_oracle_analyzer_sha256"],
            bell.Q043_RAW_ORACLE_ANALYZER_SHA256,
        )
        self.assertEqual(
            record["dependency_contract"][
                "foundational_raw_oracle_deck_manifest_sha256"
            ],
            bell.Q043_RAW_ORACLE_DECK_MANIFEST_SHA256,
        )
        self.assertEqual(
            record["dependency_contract"]["foundational_registered_admission_binding_status"],
            bell.Q043_REGISTERED_ADMISSION_BINDING_STATUS,
        )
        self.assertEqual(
            record["dependency_contract"]["foundational_supersession_sha256"],
            bell.Q043_SUPERSESSION_SHA256,
        )
        self.assertFalse(
            record["dependency_contract"]["registered_execution_admission_digest_bound"]
        )
        self.assertFalse(
            record["dependency_contract"]["registered_execution_admission_schema_bound"]
        )
        self.assertFalse(
            record["dependency_contract"]["stale_95a27328e_hash_binding_permitted"]
        )
        self.assertEqual(
            record["dependency_contract"]["required_foundational_case_count"],
            bell.Q043_REQUIRED_CASE_COUNT,
        )
        self.assertFalse(record["authority"]["launch_authorized"])
        self.assertFalse(record["authority"]["qualification_eligible"])
        self.assertFalse(record["authority"]["scientific_claim_authorized"])
        self.assertFalse(record["authority"]["publication_authorized"])
        self.assertEqual(record["matrix_contract"]["deck_case_count"], 55)
        self.assertEqual(record["matrix_contract"]["trace_record_count"], 55)
        self.assertTrue(record["gates"]["signed_phase_frequency_required_without_absolute_value"])
        self.assertEqual(record["gates"]["minimum_growth_fit_r2"], bell.MIN_GROWTH_FIT_R2)
        self.assertEqual(record["gates"]["minimum_phase_fit_r2"], bell.MIN_PHASE_FIT_R2)
        self.assertTrue(record["gates"]["multidirectional_decomposition_required"])
        for relative, digest in record["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), digest)


if __name__ == "__main__":
    unittest.main()
