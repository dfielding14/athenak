#!/usr/bin/env python3
"""Focused tests for the strict Q011 Figure 8 evolution renderer."""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import matplotlib as mpl
import numpy as np

try:
    from tst.publication import make_q011_section54_evolution_v1 as evolution
    from tst.publication.pvtk_particles import ParticleVTKData
except ModuleNotFoundError:
    import make_q011_section54_evolution_v1 as evolution
    from pvtk_particles import ParticleVTKData


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_products(
    root: Path, indices: tuple[int, ...], *, basename: str = "q011_run"
) -> None:
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "pvtk").mkdir(parents=True, exist_ok=True)
    for index in indices:
        token = f"{index:05d}"
        (root / "bin" / f"{basename}.mhd_w_bcc.{token}.bin").write_bytes(
            f"mhd-{index}".encode("ascii")
        )
        (root / "pvtk" / f"{basename}.prtcl_all.{token}.part.vtk").write_bytes(
            f"particles-{index}".encode("ascii")
        )


def _empty_particles() -> ParticleVTKData:
    integer_names = {"gid", "ptag", "species", "cr_source"}
    scalars = {
        name: np.empty(0, dtype=np.int64 if name in integer_names else np.float64)
        for name in evolution.figure8.PVTK_SCALARS
    }
    return ParticleVTKData(
        points=np.empty((0, 3), dtype=np.float64),
        scalars=scalars,
        vectors={"vel": np.empty((0, 3), dtype=np.float64)},
    )


def _particle_header(path: Path, *, cycle: int = 0, time: float = 0.0) -> None:
    path.write_bytes(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {time}  nranks= 2  cycle={cycle}  "
            "variables=prtcl_all\n"
            "BINARY\n"
        ).encode("ascii")
    )


def _reduction(
    *,
    cycle: int,
    time: float,
    front: float,
    density: np.ndarray,
    bmag: np.ndarray,
    phase_density: np.ndarray,
) -> dict[str, object]:
    x1_faces = np.asarray([front - 800.0, front + 200.0, front + 1200.0])
    return {
        "metrics": {
            "snapshot": {
                "cycle": cycle,
                "mhd_time": time,
                "particle_time": time,
                "particle_nranks": 8,
                "mesh_time_projection": time,
            },
            "shock": {
                "ideal_surface_x1": front - 5.0,
                "detected_front_x1": front,
                "density_gradient": -0.25,
                "offset_from_ideal_surface": 5.0,
            },
            "particle_filter": {"selected_particle_count": 12},
            "phase_space": {"histogrammed_energy_fraction": 0.75},
        },
        "density": np.asarray(density, dtype=np.float64),
        "bmag": np.asarray(bmag, dtype=np.float64),
        "x1_faces": x1_faces,
        "x2_faces": np.asarray([0.0, 1.0, 2.0]),
        "log10_chi_edges": np.asarray([-1.0, 0.0, 1.0]),
        "phase_density": np.asarray(phase_density, dtype=np.float64),
    }


class Q011Section54EvolutionV1Tests(unittest.TestCase):
    def test_discovery_pairs_exact_indices_and_rejects_bad_inventories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "valid"
            _write_products(root, (0, 1, 2))
            pairs = evolution.discover_product_pairs(root)
            self.assertEqual([pair.output_index for pair in pairs], [0, 1, 2])
            self.assertEqual(
                [pair.index_token for pair in pairs], ["00000", "00001", "00002"]
            )

        cases = ("missing", "gap", "duplicate", "mismatched_basename", "malformed")
        for case in cases:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                if case == "gap":
                    _write_products(root, (0, 2))
                else:
                    _write_products(root, (0, 1))
                if case == "missing":
                    (root / "pvtk" / "q011_run.prtcl_all.00001.part.vtk").unlink()
                elif case == "duplicate":
                    (root / "bin" / "other.mhd_w_bcc.00001.bin").write_bytes(b"x")
                elif case == "mismatched_basename":
                    particle = root / "pvtk" / "q011_run.prtcl_all.00001.part.vtk"
                    particle.rename(root / "pvtk" / "other.prtcl_all.00001.part.vtk")
                elif case == "malformed":
                    (root / "bin" / "q011_run.mhd_w_bcc.0000.bin").write_bytes(b"x")
                with self.assertRaises(evolution.EvolutionError):
                    evolution.discover_product_pairs(root)

    def test_cycle_zero_pair_is_explicitly_validated_and_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mhd_path = root / "q011.mhd_w_bcc.00000.bin"
            particle_path = root / "q011.prtcl_all.00000.part.vtk"
            mhd_path.write_bytes(b"mhd")
            _particle_header(particle_path)
            pair = evolution.ProductPair(
                output_index=0,
                index_token="00000",
                basename="q011",
                mhd_path=mhd_path,
                particle_path=particle_path,
            )
            dataset = SimpleNamespace(cycle=0, time=0.0)
            with mock.patch.object(
                evolution.figure8.output_primitives,
                "read_athenak_binary",
                return_value=dataset,
            ), mock.patch.object(
                evolution.figure8,
                "read_particle_vtk",
                return_value=_empty_particles(),
            ):
                excluded = evolution._validate_initial_pair(pair)
            self.assertEqual(
                excluded["reason"], "cycle_zero_t0_no_particle_initial_frame"
            )
            self.assertEqual(excluded["particle_count"], 0)

            with mock.patch.object(
                evolution.figure8.output_primitives,
                "read_athenak_binary",
                return_value=dataset,
            ), mock.patch.object(
                evolution.figure8,
                "read_particle_vtk",
                return_value=_empty_particles(),
            ), mock.patch.object(
                evolution.figure8,
                "_validate_particles",
                return_value={"points": np.zeros((1, 3))},
            ):
                with self.assertRaisesRegex(evolution.EvolutionError, "no-particle"):
                    evolution._validate_initial_pair(pair)

    def test_sequence_orders_by_time_and_writes_hashed_dependency_closed_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run_root = root / "run"
            output_dir = root / "frames"
            _write_products(run_root, (0, 1, 2))
            reductions = {
                1: _reduction(
                    cycle=20,
                    time=200.0,
                    front=2000.0,
                    density=np.asarray([[4.0, 5.0], [6.0, 8.0]]),
                    bmag=np.asarray([[10.0, 20.0], [30.0, 40.0]]),
                    phase_density=np.asarray([[0.01, 1.0], [10.0, 0.0]]),
                ),
                2: _reduction(
                    cycle=10,
                    time=100.0,
                    front=1000.0,
                    density=np.asarray([[0.5, 1.0], [2.0, 3.0]]),
                    bmag=np.asarray([[1.0, 2.0], [3.0, 4.0]]),
                    phase_density=np.asarray([[0.0001, 0.1], [1.0, 2.0]]),
                ),
            }
            rendered: list[tuple[float, str, object]] = []

            def analyze(
                mhd_path: Path, _particle_path: Path, **_kwargs: object
            ) -> object:
                index = int(mhd_path.name.rsplit(".", 2)[1])
                return reductions[index]

            def render(
                frame: evolution.FrameReduction,
                path: Path,
                normalization: object,
                **_kwargs: object,
            ) -> None:
                rendered.append((frame.particle_time, path.name, normalization))
                path.write_bytes(b"\x89PNG\r\n\x1a\n" + str(frame.particle_time).encode())

            excluded = {
                "output_index": 0,
                "index_token": "00000",
                "reason": "cycle_zero_t0_no_particle_initial_frame",
                "cycle": 0,
                "mhd_time": 0.0,
                "particle_time": 0.0,
                "particle_count": 0,
            }
            with mock.patch.object(
                evolution, "_validate_initial_pair", return_value=excluded
            ), mock.patch.object(
                evolution.figure8, "analyze_figure8", side_effect=analyze
            ) as strict_analyzer, mock.patch.object(
                evolution, "render_evolution_frame", side_effect=render
            ):
                outputs = evolution.make_evolution_sequence(
                    run_root, output_dir, dpi=72
                )

            self.assertEqual(strict_analyzer.call_count, 2)
            self.assertEqual([item[0] for item in rendered], [100.0, 200.0])
            self.assertEqual(
                [item[1] for item in rendered],
                [
                    "q011_section54_evolution_v1_frame_0000.png",
                    "q011_section54_evolution_v1_frame_0001.png",
                ],
            )
            self.assertIs(rendered[0][2], rendered[1][2])
            self.assertEqual(len(outputs), 3)
            manifest = json.loads(outputs[-1].read_text(encoding="utf-8"))
            self.assertEqual(manifest["record_type"], evolution.MANIFEST_RECORD_TYPE)
            self.assertEqual(
                [frame["time"]["particle"] for frame in manifest["frames"]],
                [100.0, 200.0],
            )
            self.assertEqual(
                [frame["output_index"] for frame in manifest["frames"]], [2, 1]
            )
            self.assertEqual(
                [frame["cycle"] for frame in manifest["frames"]], [10, 20]
            )
            self.assertEqual(
                [frame["front"]["detected_front_x1"] for frame in manifest["frames"]],
                [1000.0, 2000.0],
            )
            self.assertEqual(
                manifest["analysis_parameters"]["near_shock_x_offsets"],
                [-800.0, 1200.0],
            )
            density_norm = manifest["analysis_parameters"][
                "common_color_normalization"
            ]["density_over_rho0"]
            self.assertEqual((density_norm["vmin"], density_norm["vmax"]), (0.5, 8.0))
            self.assertEqual(
                manifest["excluded_frames"][0]["reason"],
                "cycle_zero_t0_no_particle_initial_frame",
            )

            dependencies = {
                Path(record["path"]).name: record
                for record in manifest["source_dependencies"]
            }
            self.assertEqual(
                set(dependencies),
                {
                    "make_q011_section54_figure8_v1.py",
                    "analyze_q011_section54_outputs.py",
                    "q011_section54_model.py",
                    "q011_section54_particles.py",
                    "pvtk_particles.py",
                },
            )
            for record in dependencies.values():
                self.assertEqual(record["sha256"], _sha256(Path(record["path"])))
            self.assertEqual(
                manifest["generator"]["sha256"],
                _sha256(Path(manifest["generator"]["path"])),
            )
            for record, output in zip(manifest["outputs"], outputs[:-1]):
                self.assertEqual(record["sha256"], _sha256(output))
            for record in manifest["inputs"]:
                for kind in ("mhd_w_bcc", "prtcl_all"):
                    artifact = record[kind]
                    self.assertEqual(artifact["sha256"], _sha256(Path(artifact["path"])))

    def test_real_renderer_uses_agg_and_emits_png(self) -> None:
        reduction = _reduction(
            cycle=10,
            time=100.0,
            front=1000.0,
            density=np.asarray([[1.0, 2.0], [3.0, 4.0]]),
            bmag=np.asarray([[1.0, 2.0], [3.0, 4.0]]),
            phase_density=np.asarray([[0.01, 0.1], [1.0, 2.0]]),
        )
        pair = evolution.ProductPair(
            output_index=1,
            index_token="00001",
            basename="q011",
            mhd_path=Path("mhd.bin"),
            particle_path=Path("particles.vtk"),
        )
        frame = evolution._compact_reduction(pair, reduction)
        self.assertFalse(np.shares_memory(frame.density, reduction["density"]))
        normalization = evolution.common_color_normalization((frame,))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "frame.png"
            evolution.render_evolution_frame(
                frame, output, normalization, dpi=40
            )
            self.assertTrue(output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertEqual(mpl.get_backend().lower(), "agg")

    def test_cli_accepts_run_root_output_dir_and_offset_overrides(self) -> None:
        with mock.patch.object(
            evolution,
            "make_evolution_sequence",
            return_value=[Path("frame.png"), Path("manifest.json")],
        ) as make:
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                result = evolution.main(
                    [
                        "--run-root",
                        "run",
                        "--output-dir",
                        "out",
                        "--x-offset-min",
                        "-700",
                        "--x-offset-max",
                        "900",
                    ]
                )
        self.assertEqual(result, 0)
        self.assertEqual(make.call_args.args, (Path("run"), Path("out")))
        self.assertEqual(make.call_args.kwargs["x_offsets"], (-700.0, 900.0))
        self.assertEqual(
            json.loads(stdout.getvalue())["outputs"],
            ["frame.png", "manifest.json"],
        )


if __name__ == "__main__":
    unittest.main()
