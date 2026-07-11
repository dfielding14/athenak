#!/usr/bin/env python3
"""Host/source contract tests for the corrected Q-043 Bell generator."""

from __future__ import annotations

import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
HARNESS = (
    REPO_ROOT / "tst/publication/q043_bell_current_volume_aware_host_harness.cpp"
)
SOURCE = REPO_ROOT / "src/pgen/tests/q043_bell_current_volume_aware.cpp"


class Q043BellCurrentVolumeAwareHostHarnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary_directory = tempfile.TemporaryDirectory()
        binary = Path(cls._temporary_directory.name) / "q043-bell-current-volume-aware-host"
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(HARNESS),
                "-o",
                str(binary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        cls.lines = subprocess.run(
            [str(binary)],
            cwd=REPO_ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()
        cls.binary = binary

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary_directory.cleanup()

    def test_root_cell_volume_aware_closure_across_dimensions_and_resolutions(
        self,
    ) -> None:
        lines = [line.split() for line in self.lines if line.startswith("root_grid ")]
        self.assertEqual(len(lines), 6)
        target = 4.0*math.pi
        for fields in lines:
            root_nx = [int(value) for value in fields[2:5]]
            root_cell_volume = float(fields[5])
            qscale = float(fields[6])
            self.assertGreater(root_cell_volume, 0.0)
            self.assertGreater(qscale, 0.0)
            self.assertAlmostEqual(float(fields[7]), target, places=13)
            self.assertEqual(fields[8], "1")
            if all(value % 2 == 0 for value in root_nx[: int(fields[1])]):
                coarser = next(
                    (
                        other.split()
                        for other in self.lines
                        if other.startswith("root_grid ")
                        and int(other.split()[1]) == int(fields[1])
                        and [int(value) for value in other.split()[2:5]]
                        == [value // 2 if axis < int(fields[1]) else value
                            for axis, value in enumerate(root_nx)]
                    ),
                    None,
                )
                if coarser is not None:
                    refinement = 2 ** int(fields[1])
                    self.assertAlmostEqual(
                        root_cell_volume, float(coarser[5]) / refinement, places=15
                    )
                    self.assertAlmostEqual(qscale, float(coarser[6]) / refinement)

    def test_artificial_light_speed_multiplied_current_is_rejected(self) -> None:
        lines = [
            line.split()
            for line in self.lines
            if line.startswith("light_speed_multiplied ")
        ]
        self.assertEqual(len(lines), 6)
        target = 4.0*math.pi
        for fields in lines:
            self.assertAlmostEqual(float(fields[3]), 2500.0*target, places=9)
            self.assertEqual(fields[4], "0")

    def test_root_cell_volume_is_meshblock_decomposition_invariant(self) -> None:
        lines = [
            line.split() for line in self.lines if line.startswith("decomposition ")
        ]
        self.assertEqual(len(lines), 9)
        volumes = {float(fields[3]) for fields in lines}
        self.assertEqual(len(volumes), 1)

    def test_required_current_is_invariant_to_artificial_light_speed(self) -> None:
        lines = [line.split() for line in self.lines if line.startswith("target ")]
        self.assertEqual([float(fields[1]) for fields in lines], [25.0, 2500.0, 250000.0])
        for fields in lines:
            self.assertAlmostEqual(float(fields[2]), 4.0*math.pi, places=13)

    def test_corrected_q023_unstable_eigenmode_carrier_is_reused(self) -> None:
        lines = [
            line.split() for line in self.lines if line.startswith("shared_carrier ")
        ]
        self.assertEqual(len(lines), 1)
        values = [float(value) for value in lines[0][1:]]
        self.assertTrue(all(math.isfinite(value) for value in values))
        self.assertAlmostEqual(sum(value*value for value in values[:3]), 1.0, places=11)
        self.assertAlmostEqual(sum(value*value for value in values[3:]), 1.0e-12, places=24)
        self.assertGreater(values[2], 0.0)
        self.assertLess(values[5], 0.0)

    def test_uniform_oracle_and_linear_modes_are_disjoint(self) -> None:
        uniform = next(
            line.split() for line in self.lines if line.startswith("uniform_carrier ")
        )
        values = [float(value) for value in uniform[1:]]
        self.assertAlmostEqual(sum(value*value for value in values[:3]), 1.0, places=13)
        self.assertEqual(values[3:], [0.0, 0.0, 0.0])
        modes = next(line.split() for line in self.lines if line.startswith("source_modes "))
        self.assertEqual(modes[1:], ["1", "0", "0", "1"])
        masses = next(
            line.split() for line in self.lines if line.startswith("source_mode_masses ")
        )
        self.assertEqual(masses[1:], ["1", "1", "0", "1"])

    def test_source_rejects_fractional_ppc_before_exact_current_claim(self) -> None:
        ppc = next(line.split() for line in self.lines if line.startswith("ppc_contract "))
        self.assertEqual(ppc[1:], ["1", "1", "0", "0", "0"])
        source = SOURCE.read_text(encoding="utf-8")
        self.assertIn("PositiveIntegralPPCIsValid(ppc)", source)
        self.assertIn("PPC must be a positive integer", source)

    def test_host_contract_rejects_unknown_source_mode(self) -> None:
        result = subprocess.run(
            [str(self.binary), "unknown", "0", "2", "1", "0.5"],
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout.strip(), "invalid")

    def test_source_enforces_new_campaign_and_j_over_c_contract(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        self.assertIn('#include "q023_paper_bell_linear_joverc.hpp"', source)
        self.assertNotIn('#include "q023_paper_bell_linear.cpp"', source)
        self.assertIn('const std::string block = "q043_bell_current_volume_aware";', source)
        self.assertIn('"Q043-BELL-CURRENT-VOLUME-AWARE"', source)
        self.assertIn('"Q023-PAPER-BELL-LINEAR"', source)
        self.assertIn('"deposited_j_over_c_equals_2_b_g_k0"', source)
        self.assertIn('"root_cell_macro_charge_volume_aware"', source)
        self.assertIn('"forbidden_across_dimensions_and_resolutions"', source)
        self.assertIn(
            "return ppc*qscale*species_charge*stream_speed/root_cell_volume;", source
        )
        self.assertIn("return 2.0*b_g*k0;", source)
        self.assertIn("pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min", source)
        self.assertIn("pmy_mesh_->mesh_indcs.nx1", source)
        self.assertIn("pmy_mesh_->mesh_size.dx1*pmy_mesh_->mesh_size.dx2", source)
        self.assertNotIn("2.0*b_g*light_speed*k0", source)
        self.assertIn("PPC*deposit_qscale*species_charge*v_CR/V_root_cell", source)
        self.assertIn("MeshBlock decomposition", source)
        self.assertNotIn("ppc*qscale*(species_charge/species_mass)", source)
        self.assertNotIn("2.0*b_g*light_speed*k0", source)

    def test_source_scopes_mass_one_to_corrected_linear_mode(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        self.assertIn(
            '"species_charge/species_mass must equal omega/b_g; "',
            source,
        )
        self.assertIn(
            '"corrected_linear_eigenmode requires species_mass=1; "',
            source,
        )
        self.assertIn("species_charge/species_mass", source)
        self.assertIn(
            "HasRequiredSpeciesChargeOverMass(species_mass, species_charge, omega/b_g)",
            source,
        )
        self.assertIn("SourceModeSpeciesMassIsValid(source_mode, species_mass)", source)
        self.assertIn(
            '"species_charge_over_species_mass_equals_omega_over_b_g"',
            source,
        )
        self.assertNotIn(
            '"species_charge_over_species_mass_equals_q_over_mc"',
            source,
        )
        self.assertNotIn("ppc*qscale*(species_charge/species_mass)", source)

    def test_source_requires_section52_runtime_constraints(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        required = (
            "requires MHD and particles",
            "requires ideal-MHD EOS",
            "rejects AMR/SMR",
            "requires periodic boundaries",
            '"pic_physical_mode"',
            '"paper_mhd_pic"',
            '"paper_mhd_pic_vl2_tsc"',
            '"pic_enable_2d3v", true',
            '"pic_cr_initial_state", "velocity"',
            'pin->GetString("particles", "pic_cr_hall_mode")',
            'hall_mode.compare("off")',
            'hall_mode.compare("full")',
            "physical CR-Hall off and full modes",
            '"pic_wave_damping_mode", "off"',
            '"pic_deltaf_mode", "off"',
            '"pic_expanding_box_mode", "off"',
            'Q043VolumeAwareRequireClose("omega", omega, 1.0e-6*k0*u_a);',
            '"uniform_current_oracle"',
            '"corrected_linear_eigenmode"',
            '"uniform_zero_perturbation_parallel_stream"',
            '"section52_positive_current_unstable_eigenmode"',
            "SourceModeAmplitudeIsValid(source_mode, amplitude)",
        )
        for contract in required:
            with self.subTest(contract=contract):
                self.assertIn(contract, source)

    def test_builtin_registration_covers_fresh_restart_and_cmake(self) -> None:
        pgen = (REPO_ROOT / "src/pgen/pgen.cpp").read_text(encoding="utf-8")
        self.assertEqual(pgen.count('compare("q043_bell_current_volume_aware")'), 2)
        self.assertEqual(pgen.count("Q043BellCurrentVolumeAware(pin, false);"), 1)
        self.assertEqual(pgen.count("Q043BellCurrentVolumeAware(pin, true);"), 1)
        self.assertIn(
            "void Q043BellCurrentVolumeAware(ParameterInput *pin, const bool restart);",
            (REPO_ROOT / "src/pgen/pgen.hpp").read_text(encoding="utf-8"),
        )
        self.assertIn(
            "pgen/tests/q043_bell_current_volume_aware.cpp",
            (REPO_ROOT / "src/CMakeLists.txt").read_text(encoding="utf-8"),
        )


if __name__ == "__main__":
    unittest.main()
