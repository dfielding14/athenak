#!/usr/bin/env python3
"""Focused tests for the Q011 Section 5.4 DSA spectrum tool."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import make_q011_section54_dsa_spectrum_v1 as spectrum_tool
    from tst.publication.pvtk_particles import ParticleVTKData
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import make_q011_section54_dsa_spectrum_v1 as spectrum_tool
    from pvtk_particles import ParticleVTKData


def _input_parameters(
    *, tlim: str = "1200.0", deposit_qscale: str = "0.0009"
) -> dict[str, dict[str, str]]:
    return {
        "job": {"basename": "q011_exact_same_model"},
        "mesh": {
            "nghost": "2",
            "nx1": "4000",
            "nx2": "260",
            "nx3": "1",
            "x1min": "0.0",
            "x1max": "48000.0",
            "x2min": "0.0",
            "x2max": "3120.0",
            "x3min": "0.0",
            "x3max": "1.0",
        },
        "meshblock": {"nx1": "20", "nx2": "20", "nx3": "1"},
        "time": {"integrator": "rk2", "tlim": tlim},
        "mhd": {"eos": "ideal", "gamma": "1.66666666667"},
        "particles": {
            "particle_type": "cosmic_ray",
            "pusher": "boris_tsc",
            "nspecies": "1",
            "deposit_moments": "true",
            "deposit_order": "2",
            "deposit_qscale": deposit_qscale,
            "couple_moments_to_mhd": "true",
            "couple_moments_momentum_to_mhd": "true",
            "couple_moments_energy_to_mhd": "true",
            "couple_moments_momentum_coeff": "1.0",
            "couple_moments_energy_coeff": "1.0",
            "pic_physical_mode": "paper_mhd_pic_vl2_tsc",
            "pic_background_mode": "coupled",
            "pic_feedback_mode": "coupled",
            "pic_interp_scheme": "tsc",
            "pic_cr_light_speed": "10000.0",
            "pic_cr_initial_state": "momentum",
            "pic_deltaf_mode": "off",
            "pic_random_seed": "23050101",
        },
        "species0": {"mass": "1.0", "charge": "1.0"},
        "problem": {
            "pgen_name": "pic_parallel_shock",
            "ps_rho0": "1.0",
            "ps_p0": "1.0",
            "ps_u0": "30.0",
            "ps_b0": "1.0",
            "ps_eta": "1.0e-3",
            "ps_vinj_over_u0": "3.16227766017",
            "ps_shock_speed_model": "ideal_surface",
            "ps_remove_birth_time_before": "45.0",
            "ps_enable_injection": "true",
            "ps_enable_gas_subtraction": "true",
            "ps_inject_species": "0",
        },
        "output1": {
            "file_type": "bin",
            "variable": "mhd_w_bcc",
            "file_number": "5",
            "last_time": "500.0",
        },
    }


def _dataset(
    time: float,
    cycle: int,
    *,
    parameters: dict[str, dict[str, str]] | None = None,
) -> output_primitives.AthenaBinaryDataset:
    return output_primitives.AthenaBinaryDataset(
        source=f"mhd.{cycle}.bin",
        time=time,
        cycle=cycle,
        location_size=8,
        variable_size=8,
        variable_names=spectrum_tool.MHD_FIELDS,
        input_parameters=(
            _input_parameters(tlim="500.0" if time < 1200.0 else "1200.0")
            if parameters is None
            else parameters
        ),
        root_grid_shape=(4000, 260, 1),
        meshblock_shape=(20, 20, 1),
        nghost=2,
        domain_bounds=(0.0, 48000.0, 0.0, 3120.0, 0.0, 1.0),
        blocks=(),
    )


def _velocity_from_chi(chi: np.ndarray) -> np.ndarray:
    momentum = 30.0 * np.sqrt(chi)
    velocity = momentum / np.sqrt(1.0 + (momentum / 10000.0) ** 2)
    return np.column_stack(
        (velocity, np.zeros_like(velocity), np.zeros_like(velocity))
    ).astype(np.float32).astype(np.float64)


def _particles(
    time: float,
    *,
    scale: float = 1.0,
    selected_chi: np.ndarray | None = None,
) -> ParticleVTKData:
    count = 18
    surface = 10.0 * time
    x1 = np.full(count, surface - 100.0)
    x1[3] = surface + 10.0
    x1[4] = surface
    points = np.column_stack((x1, np.linspace(0.0, 170.0, count), np.full(count, 0.5)))
    points = points.astype(np.float32).astype(np.float64)
    sources = np.ones(count, dtype=np.int64)
    sources[0] = 0
    births = np.full(count, 50.0, dtype=np.float32).astype(np.float64)
    births[1] = 44.0
    weights = np.linspace(1.0, 2.7, count).astype(np.float32).astype(np.float64)
    weights[2] = 0.0
    chi = np.geomspace(2.0, 256.0 * scale, count)
    if selected_chi is not None:
        selected_values = np.asarray(selected_chi, dtype=np.float64)
        if selected_values.shape != (count - 5,):
            raise ValueError("selected_chi must contain 13 values")
        chi[5:] = selected_values
    return ParticleVTKData(
        points=points,
        scalars={
            "gid": np.arange(count, dtype=np.int64),
            "ptag": np.arange(1000, 1000 + count, dtype=np.int64),
            "species": np.zeros(count, dtype=np.int64),
            "cr_source": sources,
            "macro_weight": weights,
            "birth_time": births,
            "deltaf_f0": np.zeros(count, dtype=np.float64),
            "deltaf_weight": np.zeros(count, dtype=np.float64),
        },
        vectors={"vel": _velocity_from_chi(chi)},
    )


def _write_particle_header(path: Path, *, time: float, cycle: int) -> None:
    path.write_bytes(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {time}  nranks= 8  cycle={cycle}  "
            "variables=prtcl_all\n"
            "BINARY\n"
        ).encode("ascii")
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q011Section54DSASpectrumV1Tests(unittest.TestCase):
    def test_analysis_uses_fixed_chi_spectrum_and_preserves_energy_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            early_particle = root / "early.vtk"
            late_particle = root / "late.vtk"
            _write_particle_header(early_particle, time=500.05, cycle=50)
            _write_particle_header(late_particle, time=1200.0, cycle=120)
            datasets = (_dataset(500.05, 50), _dataset(1200.0, 120))
            particles = (_particles(500.05), _particles(1200.0, scale=2.0))
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                side_effect=datasets,
            ), mock.patch.object(
                spectrum_tool, "read_particle_vtk", side_effect=particles
            ):
                reduction = spectrum_tool.analyze_dsa_spectrum(
                    root / "early.bin",
                    early_particle,
                    root / "late.bin",
                    late_particle,
                    energy_bin_count=24,
                    minimum_particles=1,
                )

        metrics = reduction["metrics"]
        self.assertEqual(
            metrics["qualification_effect"],
            spectrum_tool.CANONICAL_QUALIFICATION_EFFECT,
        )
        self.assertEqual(metrics["deposit_qscale"]["requested"], 0.0009)
        self.assertEqual(metrics["deposit_qscale"]["observed"], 0.0009)
        self.assertTrue(metrics["deposit_qscale"]["canonical_run"])
        self.assertTrue(metrics["same_model_validation"]["passed"])
        self.assertEqual(
            metrics["same_model_validation"]["mutable_execution_parameters_ignored"],
            ["time/tlim"],
        )
        self.assertTrue(metrics["runtime_model"]["backreaction"]["energy_feedback"])
        self.assertEqual(metrics["reference_power_law"]["f_chi_exponent"], -1.5)
        self.assertFalse(
            metrics["spectrum_definition"]["common_data_bound_energy_bins"]
        )
        self.assertEqual(
            metrics["spectrum_definition"]["fixed_bin_edges"],
            list(spectrum_tool.particle_primitives.CHI_BIN_EDGES),
        )
        self.assertEqual(
            metrics["spectrum_definition"]["legacy_energy_bin_count_argument"][
                "requested"
            ],
            24,
        )
        for record in metrics["spectra"]:
            census = record["particle_filter"]["disjoint_census"]
            self.assertEqual(census["rejected_wrong_source"]["particle_count"], 1)
            self.assertEqual(census["rejected_early_birth_time"]["particle_count"], 1)
            self.assertEqual(
                census["rejected_nonpositive_macro_weight"]["particle_count"], 1
            )
            self.assertEqual(census["rejected_upstream"]["particle_count"], 1)
            self.assertEqual(census["rejected_on_ideal_surface"]["particle_count"], 1)
            self.assertEqual(record["selected_particle_count"], 13)
            self.assertEqual(record["histogrammed_particle_count"], 13)
            self.assertEqual(record["accounted_particle_count"], 13)
            self.assertAlmostEqual(record["histogram_closure_fraction"], 1.0)
            self.assertEqual(
                record["chi_bin_edges"],
                list(spectrum_tool.particle_primitives.CHI_BIN_EDGES),
            )
            self.assertEqual(len(record["normalized_chi_f_chi"]), 40)
            self.assertEqual(
                record["normalized_chi_f_chi"],
                record["weighted_spectrum"]["normalized_chi_f_chi"],
            )
            self.assertGreater(record["selected_physical_kinetic_energy"], 0.0)
        self.assertEqual(
            metrics["spectra"][0]["chi_bin_edges"],
            metrics["spectra"][1]["chi_bin_edges"],
        )

        selected = np.arange(18) >= 5
        particle_data = particles[0]
        reconstructed_velocity = particle_data.vectors["vel"][selected]
        speed_squared = np.sum(reconstructed_velocity**2, axis=1)
        gamma_squared = 1.0 / (1.0 - speed_squared / 10000.0**2)
        momentum_squared = gamma_squared * speed_squared
        energy = momentum_squared / (np.sqrt(1.0 + momentum_squared / 10000.0**2) + 1.0)
        weights = particle_data.scalars["macro_weight"][selected]
        expected = 0.0009 * float(np.sum(weights * energy))
        self.assertTrue(
            math.isclose(
                metrics["spectra"][0]["selected_physical_kinetic_energy"],
                expected,
                rel_tol=1.0e-12,
            )
        )

    def test_fixed_reducer_accounts_for_underflow_and_overflow(self) -> None:
        selected_chi = np.asarray(
            [
                0.0,
                2.0,
                3.0,
                5.0,
                9.0,
                17.0,
                33.0,
                65.0,
                129.0,
                257.0,
                513.0,
                900.0,
                4096.0,
            ]
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            early_particle = root / "early.vtk"
            late_particle = root / "late.vtk"
            _write_particle_header(early_particle, time=500.0, cycle=50)
            _write_particle_header(late_particle, time=1200.0, cycle=120)
            particle_data = (
                _particles(500.0, selected_chi=selected_chi),
                _particles(1200.0, selected_chi=selected_chi),
            )
            reducer = spectrum_tool.particle_primitives.weighted_spectrum_record
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                side_effect=(_dataset(500.0, 50), _dataset(1200.0, 120)),
            ), mock.patch.object(
                spectrum_tool, "read_particle_vtk", side_effect=particle_data
            ), mock.patch.object(
                spectrum_tool.particle_primitives,
                "weighted_spectrum_record",
                wraps=reducer,
            ) as reduce_spectrum:
                reduction = spectrum_tool.analyze_dsa_spectrum(
                    root / "early.bin",
                    early_particle,
                    root / "late.bin",
                    late_particle,
                    energy_bin_count=24,
                    minimum_particles=1,
                )

        self.assertEqual(reduce_spectrum.call_count, 2)
        selected_weights = particle_data[0].scalars["macro_weight"][5:]
        for record in reduction["metrics"]["spectra"]:
            fixed = record["weighted_spectrum"]
            self.assertEqual(
                fixed["bin_edges"],
                list(spectrum_tool.particle_primitives.CHI_BIN_EDGES),
            )
            self.assertEqual(record["underflow_count"], 1)
            self.assertEqual(record["overflow_count"], 1)
            self.assertAlmostEqual(
                record["underflow_macro_weight"], selected_weights[0]
            )
            self.assertAlmostEqual(
                record["overflow_macro_weight"], selected_weights[-1]
            )
            self.assertEqual(record["histogrammed_particle_count"], 11)
            self.assertEqual(record["accounted_particle_count"], 13)
            self.assertAlmostEqual(
                record["accounted_macro_weight"], np.sum(selected_weights)
            )
            self.assertAlmostEqual(record["histogram_closure_fraction"], 1.0)
            self.assertFalse(record["overflow_gate"]["passed"])
            self.assertEqual(len(record["normalized_chi_f_chi"]), 40)

    def test_render_plots_dimensionless_chi_f_chi(self) -> None:
        centers = [math.sqrt(2.0), math.sqrt(8.0)]
        spectra = []
        for time, displayed in ((500.0, [0.1, 0.2]), (1200.0, [0.3, 0.4])):
            spectra.append(
                {
                    "snapshot": {"particle_time_omega0_inverse": time},
                    "chi_bin_edges": [1.0, 2.0, 4.0],
                    "chi_bin_centers": centers,
                    "normalized_chi_f_chi": displayed,
                    "epsilon_squared_f_epsilon": [100.0, 200.0],
                    "underflow_count": 0,
                    "underflow_macro_weight_fraction": 0.0,
                    "overflow_count": 0,
                    "overflow_macro_weight_fraction": 0.0,
                }
            )
        reduction = {
            "metrics": {
                "spectra": spectra,
                "reference_power_law": {
                    "line_chi_range": [1.0, 4.0],
                    "anchor_chi": 2.0,
                    "anchor_normalized_chi_f_chi": 0.2,
                    "displayed_chi_f_chi_exponent": -0.5,
                },
                "deposit_qscale": {"canonical_run": True},
            }
        }
        with mock.patch.object(spectrum_tool.plt, "close") as close:
            spectrum_tool.render_dsa_spectrum(reduction, (), 72)
        figure = close.call_args.args[0]
        try:
            axis = figure.axes[0]
            np.testing.assert_allclose(axis.lines[0].get_xdata(), centers)
            np.testing.assert_allclose(axis.lines[0].get_ydata(), [0.1, 0.2])
            self.assertIn("chi", axis.get_xlabel())
            self.assertIn("chi", axis.get_ylabel())
        finally:
            spectrum_tool.plt.close(figure)

    def test_noncanonical_qscale_is_rejected_without_explicit_override(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parameters = _input_parameters(
                tlim="500.0", deposit_qscale="0.0144"
            )
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                return_value=_dataset(500.0, 50, parameters=parameters),
            ):
                with self.assertRaisesRegex(
                    spectrum_tool.DSASpectrumError, "does not exactly match"
                ):
                    spectrum_tool.analyze_dsa_spectrum(
                        root / "early.bin",
                        root / "early.vtk",
                        root / "late.bin",
                        root / "late.vtk",
                        minimum_particles=1,
                    )

    def test_cycle_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            particle = root / "early.vtk"
            _write_particle_header(particle, time=500.0, cycle=51)
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                return_value=_dataset(500.0, 50),
            ):
                with self.assertRaisesRegex(spectrum_tool.DSASpectrumError, "cycles disagree"):
                    spectrum_tool.analyze_dsa_spectrum(
                        root / "early.bin",
                        particle,
                        root / "late.bin",
                        root / "late.vtk",
                        minimum_particles=1,
                    )

    def test_same_model_and_full_backreaction_are_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            early_particle = root / "early.vtk"
            late_particle = root / "late.vtk"
            _write_particle_header(early_particle, time=500.0, cycle=50)
            _write_particle_header(late_particle, time=1200.0, cycle=120)
            drifted = _input_parameters()
            drifted["output1"]["file_number"] = "12"
            drifted["output1"]["last_time"] = "1200.0"
            drifted["problem"]["ps_eta"] = "2.0e-3"
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                side_effect=(_dataset(500.0, 50), _dataset(1200.0, 120, parameters=drifted)),
            ), mock.patch.object(
                spectrum_tool,
                "read_particle_vtk",
                side_effect=(_particles(500.0), _particles(1200.0)),
            ):
                with self.assertRaisesRegex(
                    spectrum_tool.DSASpectrumError, "problem/ps_eta"
                ):
                    spectrum_tool.analyze_dsa_spectrum(
                        root / "early.bin",
                        early_particle,
                        root / "late.bin",
                        late_particle,
                        minimum_particles=1,
                    )

            identity_drift = _input_parameters()
            identity_drift["particles"]["pic_random_seed"] = "23050102"
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                side_effect=(
                    _dataset(500.0, 50),
                    _dataset(1200.0, 120, parameters=identity_drift),
                ),
            ), mock.patch.object(
                spectrum_tool,
                "read_particle_vtk",
                side_effect=(_particles(500.0), _particles(1200.0)),
            ):
                with self.assertRaisesRegex(
                    spectrum_tool.DSASpectrumError, "exact same-model identity"
                ):
                    spectrum_tool.analyze_dsa_spectrum(
                        root / "early.bin",
                        early_particle,
                        root / "late.bin",
                        late_particle,
                        minimum_particles=1,
                    )

            feedback_off = _input_parameters()
            feedback_off["particles"]["couple_moments_energy_to_mhd"] = "false"
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                return_value=_dataset(500.0, 50, parameters=feedback_off),
            ):
                with self.assertRaisesRegex(
                    spectrum_tool.DSASpectrumError, "required Q011 coupling"
                ):
                    spectrum_tool.analyze_dsa_spectrum(
                        root / "early.bin",
                        early_particle,
                        root / "late.bin",
                        late_particle,
                        minimum_particles=1,
                    )

    def test_make_writes_hashed_png_pdf_json_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            early_mhd = root / "early.bin"
            late_mhd = root / "late.bin"
            early_particle = root / "early.vtk"
            late_particle = root / "late.vtk"
            early_mhd.write_bytes(b"early mhd binding")
            late_mhd.write_bytes(b"late mhd binding")
            _write_particle_header(early_particle, time=500.0, cycle=50)
            _write_particle_header(late_particle, time=1200.0, cycle=120)
            early_parameters = _input_parameters(
                tlim="500.0", deposit_qscale="0.0144"
            )
            late_parameters = _input_parameters(
                tlim="1200.0", deposit_qscale="0.0144"
            )
            with mock.patch.object(
                spectrum_tool.output_primitives,
                "read_athenak_binary",
                side_effect=(
                    _dataset(500.0, 50, parameters=early_parameters),
                    _dataset(1200.0, 120, parameters=late_parameters),
                ),
            ), mock.patch.object(
                spectrum_tool,
                "read_particle_vtk",
                side_effect=(_particles(500.0), _particles(1200.0, scale=2.0)),
            ):
                outputs = spectrum_tool.make_dsa_spectrum(
                    early_mhd,
                    early_particle,
                    late_mhd,
                    late_particle,
                    root / "figures",
                    energy_bin_count=24,
                    minimum_particles=1,
                    expected_deposit_qscale=0.0144,
                    dpi=72,
                )
            self.assertEqual([path.suffix for path in outputs], [".png", ".pdf", ".json", ".json"])
            self.assertTrue(all(path.is_file() and path.stat().st_size > 0 for path in outputs))
            metrics = json.loads(outputs[2].read_text(encoding="utf-8"))
            manifest = json.loads(outputs[3].read_text(encoding="utf-8"))
            self.assertEqual(metrics["record_type"], spectrum_tool.RECORD_TYPE)
            self.assertEqual(manifest["record_type"], spectrum_tool.MANIFEST_RECORD_TYPE)
            self.assertEqual(
                metrics["qualification_effect"],
                spectrum_tool.EXPLORATORY_QUALIFICATION_EFFECT,
            )
            self.assertEqual(
                manifest["qualification_effect"],
                spectrum_tool.EXPLORATORY_QUALIFICATION_EFFECT,
            )
            self.assertEqual(metrics["deposit_qscale"]["requested"], 0.0144)
            self.assertEqual(metrics["deposit_qscale"]["observed"], 0.0144)
            self.assertTrue(
                metrics["deposit_qscale"]["explicit_noncanonical_override"]
            )
            self.assertEqual(manifest["deposit_qscale"], metrics["deposit_qscale"])
            self.assertEqual(
                manifest["analysis_parameters"]["expected_deposit_qscale"],
                0.0144,
            )
            self.assertEqual(len(manifest["inputs"]), 4)
            self.assertEqual(len(manifest["outputs"]), 3)
            for artifact in manifest["inputs"] + manifest["outputs"]:
                path = Path(artifact["path"])
                self.assertEqual(artifact["sha256"], _sha256(path))
                self.assertEqual(artifact["byte_count"], path.stat().st_size)
            self.assertEqual(
                manifest["analysis_parameters"]["displayed_quantity"],
                "chi * f_chi / total_post_filter_macro_weight",
            )
            self.assertEqual(manifest["analysis_parameters"]["energy_bin_count"], 24)
            self.assertEqual(manifest["analysis_parameters"]["fixed_chi_bin_count"], 40)
            self.assertTrue(
                manifest["analysis_parameters"]["underflow_overflow_accounting"]
            )

    def test_cli_passes_explicit_expected_qscale(self) -> None:
        with mock.patch.object(
            spectrum_tool, "make_dsa_spectrum", return_value=[]
        ) as make, mock.patch("builtins.print"):
            result = spectrum_tool.main(
                [
                    "--early-mhd",
                    "early.bin",
                    "--early-particles",
                    "early.vtk",
                    "--late-mhd",
                    "late.bin",
                    "--late-particles",
                    "late.vtk",
                    "--output-dir",
                    "figures",
                    "--expected-deposit-qscale",
                    "0.0144",
                    "--energy-bin-count",
                    "24",
                ]
            )
        self.assertEqual(result, 0)
        self.assertEqual(make.call_args.kwargs["expected_deposit_qscale"], 0.0144)
        self.assertEqual(make.call_args.kwargs["energy_bin_count"], 24)


if __name__ == "__main__":
    unittest.main()
