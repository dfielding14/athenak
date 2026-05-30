"""Q-016 regression for typed PVTK provenance and bounded weighted spectra."""

from __future__ import annotations

import json
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

if __package__:
    from .pvtk_particles import read_particle_vtk
    from .q016_particle_spectra import write_weighted_species_cohort_spectra
else:
    from pvtk_particles import read_particle_vtk
    from q016_particle_spectra import write_weighted_species_cohort_spectra


def _binary(values: list[float] | list[int], code: str) -> bytes:
    return struct.pack(">" + str(len(values)) + code, *values)


def _write_fixture(path: Path) -> None:
    points = [
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        3.0, 0.0, 0.0,
    ]
    scalars = [
        ("gid", "int", "i", [0, 1, 1, 2]),
        ("ptag", "int", "i", [10, 11, 12, 13]),
        ("species", "int", "i", [0, 0, 1, 1]),
        ("cr_source", "int", "i", [0, 1, 1, 0]),
        ("macro_weight", "float", "f", [1.0, 2.0, 0.5, 1.5]),
        ("birth_time", "float", "f", [0.0, 0.25, 0.75, 0.0]),
        ("deltaf_f0", "float", "f", [1.0, 1.0, 0.5, 2.0]),
        ("deltaf_weight", "float", "f", [0.0, 0.25, -0.5, 1.0]),
    ]
    velocity = [
        0.5, 0.0, 0.0,
        1.5, 0.0, 0.0,
        2.5, 0.0, 0.0,
        3.5, 0.0, 0.0,
    ]
    payload = bytearray(
        b"# vtk DataFile Version 2.0\n"
        b"# Q016 typed fixture\n"
        b"BINARY\n"
        b"DATASET UNSTRUCTURED_GRID\n"
        b"\nPOINTS 4 float\n"
    )
    payload.extend(_binary(points, "f"))
    payload.extend(b"\n\nPOINT_DATA 4\n")
    for name, vtk_type, code, values in scalars:
        payload.extend(
            f"\nSCALARS {name} {vtk_type}\nLOOKUP_TABLE default\n".encode("ascii")
        )
        payload.extend(_binary(values, code))
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(_binary(velocity, "f"))
    path.write_bytes(payload)


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


if __name__ == "__main__":
    unittest.main()
