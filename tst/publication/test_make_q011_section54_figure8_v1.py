#!/usr/bin/env python3
"""Focused tests for the strict Q011 Figure 8 analyzer and plotter."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import make_q011_section54_figure8_v1 as figure8
    from tst.publication.pvtk_particles import ParticleVTKData
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import make_q011_section54_figure8_v1 as figure8
    from pvtk_particles import ParticleVTKData


def _dataset(
    source: str = "mhd.bin", *, cycle: int = 17
) -> output_primitives.AthenaBinaryDataset:
    nx = 12
    ny = 4
    x1 = np.linspace(3600.0, 6400.0, nx, endpoint=False) + 2800.0 / (2.0 * nx)
    density_x = np.where(x1 < 5000.0, 4.0, 1.0)
    density = np.broadcast_to(density_x[None, :], (ny, nx)).copy()
    bcc1 = np.broadcast_to((1.0 + 0.001 * x1)[None, :], (ny, nx)).copy()
    zeros = np.zeros((ny, nx))
    fields = {
        "dens": density[None, :, :],
        "velx": zeros[None, :, :],
        "vely": zeros[None, :, :],
        "velz": zeros[None, :, :],
        "eint": np.ones((1, ny, nx)),
        "bcc1": bcc1[None, :, :],
        "bcc2": zeros[None, :, :],
        "bcc3": zeros[None, :, :],
    }
    block = output_primitives.AthenaBinaryBlock(
        index_bounds=(0, nx - 1, 0, ny - 1, 0, 0),
        logical_location=(0, 0, 0),
        level=0,
        geometry=(3600.0, 6400.0, 0.0, 400.0, 0.0, 1.0),
        fields=fields,
    )
    return output_primitives.AthenaBinaryDataset(
        source=source,
        time=500.0,
        cycle=cycle,
        location_size=8,
        variable_size=8,
        variable_names=figure8.MHD_FIELDS,
        input_parameters={
            "problem": {"ps_rho0": "1.0", "ps_b0": "1.0", "ps_u0": "30.0"},
            "particles": {"pic_cr_light_speed": "10000.0", "deposit_qscale": "0.0009"},
        },
        root_grid_shape=(nx, ny, 1),
        meshblock_shape=(nx, ny, 1),
        nghost=0,
        domain_bounds=(3600.0, 6400.0, 0.0, 400.0, 0.0, 1.0),
        blocks=(block,),
    )


def _velocity_from_chi(chi: np.ndarray) -> np.ndarray:
    momentum = 30.0 * np.sqrt(chi)
    velocity = momentum / np.sqrt(1.0 + (momentum / 10000.0) ** 2)
    return np.column_stack(
        (velocity, np.zeros_like(velocity), np.zeros_like(velocity))
    ).astype(np.float32).astype(np.float64)


def _particles() -> ParticleVTKData:
    count = 8
    points = np.column_stack(
        (
            np.linspace(4100.0, 5900.0, count),
            np.linspace(25.0, 375.0, count),
            np.full(count, 0.5),
        )
    ).astype(np.float32).astype(np.float64)
    scalars = {
        "gid": np.arange(count, dtype=np.int64),
        "ptag": np.arange(100, 100 + count, dtype=np.int64),
        "species": np.zeros(count, dtype=np.int64),
        "cr_source": np.asarray([1, 1, 1, 1, 1, 0, 1, 1], dtype=np.int64),
        "macro_weight": np.asarray([1, 1, 1, 1, 1, 1, 0, 1], dtype=np.float32).astype(
            np.float64
        ),
        "birth_time": np.asarray(
            [45, 46, 47, 48, 49, 50, 45, 44], dtype=np.float32
        ).astype(np.float64),
        "deltaf_f0": np.zeros(count, dtype=np.float64),
        "deltaf_weight": np.zeros(count, dtype=np.float64),
    }
    return ParticleVTKData(
        points=points,
        scalars=scalars,
        vectors={"vel": _velocity_from_chi(np.linspace(2.0, 32.0, count))},
    )


def _write_particle_header(path: Path, *, cycle: int = 17, time: float = 500.0) -> None:
    path.write_bytes(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {time}  nranks= 2  cycle={cycle}  "
            "variables=prtcl_all\n"
            "BINARY\n"
        ).encode("ascii")
    )


class Q011Section54Figure8V1Tests(unittest.TestCase):
    def test_analysis_filters_particles_and_closes_histogram(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            particle_path = root / "particles.vtk"
            _write_particle_header(particle_path)
            with mock.patch.object(
                figure8.output_primitives, "read_athenak_binary", return_value=_dataset()
            ), mock.patch.object(figure8, "read_particle_vtk", return_value=_particles()):
                result = figure8.analyze_figure8(root / "mhd.bin", particle_path)
        metrics = result["metrics"]
        self.assertEqual(metrics["snapshot"]["cycle"], 17)
        self.assertEqual(metrics["particle_filter"]["selected_particle_count"], 5)
        self.assertEqual(metrics["particle_filter"]["rejected_wrong_source_count"], 1)
        self.assertEqual(metrics["particle_filter"]["rejected_early_birth_count"], 1)
        self.assertEqual(
            metrics["particle_filter"]["rejected_nonpositive_weight_count"], 1
        )
        self.assertGreater(metrics["phase_space"]["histogrammed_energy_fraction"], 0.0)
        self.assertLessEqual(metrics["phase_space"]["histogrammed_energy_fraction"], 1.0)
        self.assertEqual(
            int(np.sum(metrics["phase_space"]["counts"])),
            metrics["phase_space"]["histogrammed_particle_count"],
        )

    def test_cycle_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            particle_path = root / "particles.vtk"
            _write_particle_header(particle_path, cycle=18)
            with mock.patch.object(
                figure8.output_primitives, "read_athenak_binary", return_value=_dataset()
            ):
                with self.assertRaisesRegex(figure8.Figure8Error, "cycles disagree"):
                    figure8.analyze_figure8(root / "mhd.bin", particle_path)

    def test_make_writes_png_pdf_metrics_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mhd_path = root / "mhd.bin"
            particle_path = root / "particles.vtk"
            mhd_path.write_bytes(b"synthetic mhd binding")
            _write_particle_header(particle_path)
            with mock.patch.object(
                figure8.output_primitives, "read_athenak_binary", return_value=_dataset()
            ), mock.patch.object(figure8, "read_particle_vtk", return_value=_particles()):
                outputs = figure8.make_figure8(
                    mhd_path, particle_path, root / "figures", dpi=72
                )
            self.assertEqual(len(outputs), 4)
            self.assertTrue(
                all(path.is_file() and path.stat().st_size > 0 for path in outputs)
            )
            metrics = json.loads(outputs[2].read_text(encoding="utf-8"))
            manifest = json.loads(outputs[3].read_text(encoding="utf-8"))
            self.assertEqual(metrics["record_type"], figure8.RECORD_TYPE)
            self.assertEqual(manifest["record_type"], figure8.MANIFEST_RECORD_TYPE)
            self.assertFalse(manifest["analysis_parameters"]["cross_cycle_substitution"])
            self.assertEqual({item["path"] for item in manifest["inputs"]}, {
                str(mhd_path.resolve()), str(particle_path.resolve())
            })


if __name__ == "__main__":
    unittest.main()
