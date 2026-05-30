"""Q-016 regression for typed PVTK provenance and bounded weighted spectra."""

from __future__ import annotations

import json
import os
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "tst"))

if __package__:
    from .pvtk_particles import read_particle_vtk
    from .q016_particle_spectra import write_weighted_species_cohort_spectra
else:
    from pvtk_particles import read_particle_vtk
    from q016_particle_spectra import write_weighted_species_cohort_spectra

from scripts.particles import pic_q016_particle_provenance as q016_harness


_POINTS = [
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    3.0, 0.0, 0.0,
]
_SCALARS = [
    ("gid", "int", "i", [0, 1, 1, 2]),
    ("ptag", "int", "i", [10, 11, 12, 13]),
    ("species", "int", "i", [0, 0, 1, 1]),
    ("cr_source", "int", "i", [0, 1, 1, 0]),
    ("macro_weight", "float", "f", [1.0, 2.0, 0.5, 1.5]),
    ("birth_time", "float", "f", [0.0, 0.25, 0.75, 0.0]),
    ("deltaf_f0", "float", "f", [1.0, 1.0, 0.5, 2.0]),
    ("deltaf_weight", "float", "f", [0.0, 0.25, -0.5, 1.0]),
]
_VELOCITY = [
    0.5, 0.0, 0.0,
    1.5, 0.0, 0.0,
    2.5, 0.0, 0.0,
    3.5, 0.0, 0.0,
]


def _binary(values: list[float] | list[int], code: str) -> bytes:
    return struct.pack(">" + str(len(values)) + code, *values)


def _select(values: list[float] | list[int], indices: list[int], width: int = 1):
    return [
        values[width*index + offset]
        for index in indices
        for offset in range(width)
    ]


def _scalar_section(name: str, vtk_type: str, code: str, values) -> bytes:
    return (
        f"\nSCALARS {name} {vtk_type}\nLOOKUP_TABLE default\n".encode("ascii")
        + _binary(values, code)
    )


def _fixture_bytes(indices: list[int] | None = None) -> bytes:
    if indices is None:
        indices = list(range(4))
    payload = bytearray(
        b"# vtk DataFile Version 2.0\n"
        b"# Q016 typed fixture\n"
        b"BINARY\n"
        b"DATASET UNSTRUCTURED_GRID\n"
        + f"\nPOINTS {len(indices)} float\n".encode("ascii")
    )
    payload.extend(_binary(_select(_POINTS, indices, width=3), "f"))
    payload.extend(f"\n\nPOINT_DATA {len(indices)}\n".encode("ascii"))
    for name, vtk_type, code, values in _SCALARS:
        payload.extend(_scalar_section(name, vtk_type, code, _select(values, indices)))
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(_binary(_select(_VELOCITY, indices, width=3), "f"))
    return bytes(payload)


def _write_fixture(path: Path, indices: list[int] | None = None) -> None:
    path.write_bytes(_fixture_bytes(indices))


def _independent_histogram(
    speed: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
) -> np.ndarray:
    reconstructed = np.zeros(edges.size - 1, dtype=np.float64)
    for value, weight in zip(speed, weights):
        for index in range(edges.size - 1):
            if edges[index] <= value < edges[index + 1]:
                reconstructed[index] += weight
                break
            if index == edges.size - 2 and value == edges[index + 1]:
                reconstructed[index] += weight
                break
    return reconstructed


class TestQ016ParticleProvenanceSpectra(unittest.TestCase):
    def _assert_reader_rejects(self, contents: bytes, message: str) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = Path(tmpdir) / "q016.part.vtk"
            fixture.write_bytes(contents)
            with self.assertRaisesRegex(ValueError, message):
                read_particle_vtk(fixture)

    def test_typed_reader_and_independent_weighted_reconstruction(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = Path(tmpdir) / "q016.part.vtk"
            output = Path(tmpdir) / "spectra.json"
            _write_fixture(fixture)
            data = read_particle_vtk(fixture)

            self.assertEqual(data.scalars["ptag"].dtype.kind, "i")
            self.assertEqual(data.scalars["species"].tolist(), [0, 0, 1, 1])
            self.assertEqual(data.scalars["cr_source"].tolist(), [0, 1, 1, 0])

            edges = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0])
            birth_edges = np.asarray([0.0, 0.5, 1.0])
            speed = np.linalg.norm(data.vectors["vel"], axis=1)
            macro_weight = data.scalars["macro_weight"]
            delta_weight = data.scalars["deltaf_weight"]

            full = write_weighted_species_cohort_spectra(
                data, output, edges, birth_edges, delta_f_semantics="full_f"
            )
            expected_full = _independent_histogram(speed, macro_weight, edges)
            np.testing.assert_allclose(
                full["all_particles"]["weighted_sum_in_bins"], expected_full
            )
            np.testing.assert_allclose(
                full["all_particles"]["weighted_density_per_bin_width"],
                expected_full / np.diff(edges),
            )
            self.assertEqual(full["normalization"], "weighted_sum_in_bins / diff(bin_edges)")
            self.assertEqual(len(full["groups"]), 4)
            self.assertEqual(json.loads(output.read_text(encoding="utf-8")), full)

            perturbation = write_weighted_species_cohort_spectra(
                data,
                output,
                edges,
                birth_edges,
                delta_f_semantics="delta_f_perturbation",
            )
            expected_delta = _independent_histogram(
                speed, macro_weight * delta_weight, edges
            )
            np.testing.assert_allclose(
                perturbation["all_particles"]["weighted_sum_in_bins"],
                expected_delta,
            )
            self.assertIn("signed sampled perturbation", perturbation["delta_f_definition"])

    def test_reader_rejects_missing_or_mismatched_point_data(self) -> None:
        contents = _fixture_bytes()
        self._assert_reader_rejects(
            contents.replace(b"\n\nPOINT_DATA 4\n", b"\n", 1),
            "before POINT_DATA",
        )
        self._assert_reader_rejects(
            contents.replace(b"POINT_DATA 4", b"POINT_DATA 3", 1),
            "does not match",
        )

    def test_reader_rejects_duplicate_sections(self) -> None:
        contents = _fixture_bytes()
        self._assert_reader_rejects(
            contents.replace(
                b"\n\nPOINT_DATA 4\n",
                b"\n\nPOINT_DATA 4\n\nPOINT_DATA 4\n",
                1,
            ),
            "Duplicate POINT_DATA",
        )
        duplicate_scalar = _scalar_section("ptag", "int", "i", [10, 11, 12, 13])
        self._assert_reader_rejects(
            contents.replace(
                b"\nVECTORS vel float\n",
                duplicate_scalar + b"\nVECTORS vel float\n",
            ),
            "Duplicate SCALARS",
        )
        self._assert_reader_rejects(
            contents + b"\nVECTORS vel float\n" + _binary(_VELOCITY, "f"),
            "Duplicate VECTORS",
        )
        scalar_named_vel = _scalar_section("vel", "float", "f", [1.0, 1.0, 1.0, 1.0])
        self._assert_reader_rejects(
            contents.replace(
                b"\nVECTORS vel float\n",
                scalar_named_vel + b"\nVECTORS vel float\n",
            ),
            "Duplicate VECTORS",
        )

    def test_reader_rejects_vector_type_and_trailing_unknown_content(self) -> None:
        contents = _fixture_bytes()
        self._assert_reader_rejects(
            contents.replace(b"VECTORS vel float", b"VECTORS vel int", 1),
            "Unsupported VECTORS type",
        )
        self._assert_reader_rejects(contents + b"\nFORGED\n", "Unexpected VTK content")

    def test_harness_merges_same_cycle_gid_slices_before_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            pvtk = root / "pvtk"
            pvtk.mkdir()
            basename = "pic_q016_partitioned"
            _write_fixture(
                pvtk / f"{basename}.prtcl_all.7.00002.part.vtk",
                indices=[1, 3],
            )
            _write_fixture(
                pvtk / f"{basename}.prtcl_all.3.00002.part.vtk",
                indices=[0, 2],
            )
            with mock.patch.dict(os.environ, {"ATHENA_Q016_EXE_DIR": str(root)}):
                cycles = q016_harness._particle_vtk_cycles(basename)
                snapshot = q016_harness._snapshot(cycles[2])
            self.assertEqual(sorted(cycles), [2])
            self.assertEqual(snapshot["ptag"].tolist(), [10, 11, 12, 13])
            self.assertEqual(snapshot["gid"].tolist(), [0, 1, 1, 2])

    def test_harness_rejects_duplicate_tags_across_gid_slices(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            first = Path(tmpdir) / "first.part.vtk"
            second = Path(tmpdir) / "second.part.vtk"
            _write_fixture(first, indices=[0, 1])
            _write_fixture(second, indices=[1, 2])
            with self.assertRaisesRegex(RuntimeError, "ptag values must be unique"):
                q016_harness._snapshot([str(first), str(second)])

    def test_harness_launcher_configuration_preserves_serial_default(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(q016_harness._launcher_prefix(), [])
        with mock.patch.dict(
            os.environ,
            {"ATHENA_Q016_NPROC": "4", "ATHENA_Q016_LAUNCHER": "srun --cpu-bind=cores"},
            clear=True,
        ):
            self.assertEqual(
                q016_harness._launcher_prefix(),
                ["srun", "--cpu-bind=cores", "-n", "4"],
            )
        with mock.patch.dict(os.environ, {"ATHENA_Q016_NPROC": "0"}, clear=True):
            with self.assertRaisesRegex(RuntimeError, "positive integer"):
                q016_harness._launcher_prefix()

    def test_harness_rejects_multi_rank_launch_with_serial_executable(self) -> None:
        with mock.patch.dict(os.environ, {"ATHENA_Q016_NPROC": "2"}, clear=True):
            with mock.patch.object(q016_harness, "_athena_mpi_enabled", return_value=False):
                with self.assertRaisesRegex(RuntimeError, "MPI-enabled"):
                    q016_harness._require_supported_rank_configuration()


if __name__ == "__main__":
    unittest.main()
