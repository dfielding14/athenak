#!/usr/bin/env python3
"""Focused tests for the bounded Q-011 Section 5.4 model helpers."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

PUBLICATION_DIR = Path(__file__).resolve().parent
REPO_ROOT = PUBLICATION_DIR.parents[1]
sys.path.insert(0, str(PUBLICATION_DIR))

import q011_section54_model as model  # noqa: E402


def _amr_deck() -> str:
    return """\
<mesh>
nx1 = 4000
x1min = 0.0
x1max = 48000.0
nx2 = 260
x2min = 0.0
x2max = 3120.0
nx3 = 1
x3min = 0.0
x3max = 1.0

<meshblock>
nx1 = 20
nx2 = 20
nx3 = 1

<mesh_refinement>
refinement = adaptive
num_levels = 3

<mhd>
gamma = 1.66666666667

<particles>
pic_physical_mode = paper_mhd_pic_vl2_tsc
pic_cr_light_speed = 10000.0
pic_cr_initial_state = momentum

<problem>
pgen_name = pic_parallel_shock
ps_u0 = 30.0
ps_shock_speed_model = ideal_surface
ps_enable_curvature_amr = true
"""


def _runtime_line(**replacements: str) -> str:
    values = {
        "physical_mode": model.PAPER_PHYSICAL_MODE,
        "state": model.RUNTIME_MOMENTUM_STATE,
        "C": "10000",
        "background": "coupled",
        "feedback": "coupled",
        "restart_schema": "7",
    }
    values.update(replacements)
    return "PIC runtime model: " + " ".join(
        f"{name}={value}" for name, value in values.items()
    )


class Q011Section54VariantBindingTests(unittest.TestCase):
    def test_canonical_variant_identity_rejects_aliases(self) -> None:
        canonical = (
            "coarse_uniform_dx12",
            "three_level_amr_root_dx12_finest_dx3",
            "fine_uniform_dx3",
        )
        for candidate in canonical:
            with self.subTest(candidate=candidate):
                self.assertEqual(model.canonical_variant_identity(candidate), candidate)
        for candidate in (
            "amr",
            "fine_uniform",
            "amr_fiducial",
            "",
            1,
        ):
            with self.subTest(candidate=candidate):
                with self.assertRaises(model.ModelContractError):
                    model.canonical_variant_identity(candidate)

    def test_fine_uniform_cannot_reuse_bare_amr_binding(self) -> None:
        amr = model.variant_binding("three_level_amr_root_dx12_finest_dx3")
        fine = model.variant_binding("fine_uniform_dx3")
        self.assertEqual(amr.deck_path, fine.deck_path)
        self.assertEqual(amr.model_launch_overrides, ())
        self.assertNotEqual(fine.model_launch_overrides, ())
        with self.assertRaisesRegex(
            model.ModelContractError,
            "fine_uniform_dx3: model launch overrides differ",
        ):
            model.require_variant_binding("fine_uniform_dx3", fine.deck_path, ())

    def test_variant_binding_rejects_deck_and_override_drift(self) -> None:
        fine = model.variant_binding("fine_uniform_dx3")
        with self.assertRaisesRegex(model.ModelContractError, "deck path differs"):
            model.require_variant_binding(
                "fine_uniform_dx3",
                "inputs/publication/common_amr_decoy.athinput",
                fine.model_launch_overrides,
            )
        with self.assertRaisesRegex(model.ModelContractError, "duplicate mesh/nx1"):
            model.require_variant_binding(
                "fine_uniform_dx3",
                fine.deck_path,
                ("mesh/nx1=16000", "mesh/nx1=16000"),
            )

    def test_bound_variants_materialize_distinct_typed_layouts(self) -> None:
        coarse = model.parse_bound_variant_deck_contract(
            "coarse_uniform_dx12",
            model.BASE_DECK_PATH,
            _amr_deck(),
            model.variant_binding("coarse_uniform_dx12").model_launch_overrides,
        )
        amr = model.parse_bound_variant_deck_contract(
            "three_level_amr_root_dx12_finest_dx3",
            model.BASE_DECK_PATH,
            _amr_deck(),
            model.variant_binding(
                "three_level_amr_root_dx12_finest_dx3"
            ).model_launch_overrides,
        )
        fine = model.parse_bound_variant_deck_contract(
            "fine_uniform_dx3",
            model.BASE_DECK_PATH,
            _amr_deck(),
            model.variant_binding("fine_uniform_dx3").model_launch_overrides,
        )
        self.assertEqual(coarse.variant, "coarse_uniform_dx12")
        self.assertEqual(coarse.refinement, "none")
        self.assertEqual(coarse.num_levels, 1)
        self.assertAlmostEqual(coarse.root_dx, 12.0)
        self.assertAlmostEqual(coarse.finest_dx, 12.0)
        self.assertEqual(amr.variant, "three_level_amr_root_dx12_finest_dx3")
        self.assertEqual(amr.refinement, "adaptive")
        self.assertEqual(amr.num_levels, 3)
        self.assertAlmostEqual(amr.root_dx, 12.0)
        self.assertAlmostEqual(amr.finest_dx, 3.0)
        self.assertEqual(fine.variant, "fine_uniform_dx3")
        self.assertEqual(fine.refinement, "none")
        self.assertEqual(fine.num_levels, 1)
        self.assertAlmostEqual(fine.root_dx, 3.0)
        self.assertAlmostEqual(fine.finest_dx, 3.0)

    def test_checked_in_active_deck_satisfies_amr_binding(self) -> None:
        deck = REPO_ROOT / model.BASE_DECK_PATH
        contract = model.parse_bound_variant_deck_contract(
            "three_level_amr_root_dx12_finest_dx3",
            model.BASE_DECK_PATH,
            deck.read_text(encoding="utf-8"),
            (),
        )
        self.assertEqual(contract.variant, "three_level_amr_root_dx12_finest_dx3")


class Q011Section54RuntimeIdentityTests(unittest.TestCase):
    def test_runtime_identity_requires_frozen_paper_projection(self) -> None:
        identity = model.parse_runtime_identity("startup\n" + _runtime_line() + "\n")
        self.assertEqual(identity.physical_mode, model.PAPER_PHYSICAL_MODE)
        self.assertEqual(identity.state, "momentum_p_over_m")
        self.assertEqual(identity.light_speed, 10000.0)
        self.assertEqual(identity.restart_schema, 7)

    def test_runtime_identity_rejects_missing_duplicate_and_malformed_fields(self) -> None:
        cases = {
            "missing": _runtime_line().replace(" restart_schema=7", ""),
            "duplicate": _runtime_line() + " C=10000",
            "malformed": _runtime_line() + " malformed",
            "nonfinite_C": _runtime_line(C="nan"),
            "decimal_schema": _runtime_line(restart_schema="7.0"),
        }
        for label, line in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(model.ModelContractError):
                    model.parse_runtime_identity_line(line)

    def test_runtime_identity_rejects_required_value_drift(self) -> None:
        cases = {
            "historical_mode": {"physical_mode": "paper_mhd_pic"},
            "velocity_state": {"state": "velocity"},
            "light_speed": {"C": "9999"},
            "restart_schema": {"restart_schema": "6"},
        }
        for label, replacements in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(model.ModelContractError):
                    model.parse_runtime_identity_line(_runtime_line(**replacements))

    def test_runtime_identity_rejects_duplicate_model_lines(self) -> None:
        line = _runtime_line()
        with self.assertRaisesRegex(model.ModelContractError, "exactly one"):
            model.parse_runtime_identity(line + "\n" + line + "\n")


class Q011Section54IdealSurfaceTests(unittest.TestCase):
    def test_x_ideal_is_ten_times_nonnegative_finite_time(self) -> None:
        self.assertEqual(model.x_ideal(0.0), 0.0)
        self.assertEqual(model.x_ideal(12.5), 125.0)

    def test_x_ideal_rejects_invalid_time(self) -> None:
        for value in (-1.0, math.inf, -math.inf, math.nan, True, 1.0e308):
            with self.subTest(value=value):
                with self.assertRaises(model.ModelContractError):
                    model.x_ideal(value)


class Q011Section54ChiReconstructionTests(unittest.TestCase):
    def test_chi_reconstructs_from_float_physical_velocity(self) -> None:
        expected = (25.0 / (1.0 - 25.0 / 10000.0**2)) / 30.0**2
        self.assertAlmostEqual(
            model.reconstruct_chi_from_physical_speed(
                5.0,
                physical_mode=model.PAPER_PHYSICAL_MODE,
            ),
            expected,
        )
        self.assertAlmostEqual(
            model.reconstruct_chi_from_physical_velocity(
                (3.0, 4.0, 0.0),
                physical_mode=model.PAPER_PHYSICAL_MODE,
            ),
            expected,
        )

    def test_chi_rejects_nonfinite_and_superluminal_velocity(self) -> None:
        vectors = (
            (math.nan, 0.0, 0.0),
            (math.inf, 0.0, 0.0),
            (10000.0, 0.0, 0.0),
            (10001.0, 0.0, 0.0),
            (1.0, 2.0),
            "1,2,3",
        )
        for velocity in vectors:
            with self.subTest(velocity=velocity):
                with self.assertRaises(model.ModelContractError):
                    model.reconstruct_chi_from_physical_velocity(
                        velocity,
                        physical_mode=model.PAPER_PHYSICAL_MODE,
                    )

    def test_chi_rejects_nonpaper_mode_and_negative_speed(self) -> None:
        with self.assertRaisesRegex(model.ModelContractError, "paper physical mode"):
            model.reconstruct_chi_from_physical_velocity(
                (3.0, 4.0, 0.0),
                physical_mode="paper_mhd_pic",
            )
        with self.assertRaisesRegex(model.ModelContractError, "nonnegative magnitude"):
            model.reconstruct_chi_from_physical_speed(
                -1.0,
                physical_mode=model.PAPER_PHYSICAL_MODE,
            )

    def test_float32_projection_uncertainty_is_documented(self) -> None:
        note = model.FLOAT32_PROJECTION_UNCERTAINTY
        self.assertIn("float32", note)
        self.assertIn("not an exact recovery", note)
        self.assertIn("approaches C", note)


class Q011Section54DeckContractTests(unittest.TestCase):
    def test_typed_deck_contract_freezes_model_assumptions(self) -> None:
        contract = model.parse_deck_contract(_amr_deck())
        self.assertEqual(contract.variant, "three_level_amr_root_dx12_finest_dx3")
        self.assertEqual(contract.physical_mode, model.PAPER_PHYSICAL_MODE)
        self.assertEqual(contract.initial_state, "momentum")
        self.assertEqual(contract.light_speed, 10000.0)
        self.assertAlmostEqual(contract.gamma, 5.0 / 3.0)
        self.assertEqual(contract.upstream_speed_u0, 30.0)
        self.assertEqual(contract.shock_speed_model, "ideal_surface")

    def test_typed_deck_contract_rejects_duplicate_and_malformed_values(self) -> None:
        cases = {
            "duplicate": _amr_deck().replace(
                "pic_cr_light_speed = 10000.0",
                "pic_cr_light_speed = 10000.0\npic_cr_light_speed = 10000.0",
            ),
            "malformed_integer": _amr_deck().replace("num_levels = 3", "num_levels = 3.0"),
            "malformed_boolean": _amr_deck().replace(
                "ps_enable_curvature_amr = true",
                "ps_enable_curvature_amr = enabled",
            ),
        }
        for label, deck in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(model.ModelContractError):
                    model.parse_deck_contract(deck)

    def test_typed_deck_contract_rejects_assumption_drift(self) -> None:
        replacements = {
            "physical_mode": (
                "pic_physical_mode = paper_mhd_pic_vl2_tsc",
                "pic_physical_mode = paper_mhd_pic",
            ),
            "light_speed": (
                "pic_cr_light_speed = 10000.0",
                "pic_cr_light_speed = 3.0",
            ),
            "surface_model": (
                "ps_shock_speed_model = ideal_surface",
                "ps_shock_speed_model = finite_mach",
            ),
            "layout": ("nx1 = 4000", "nx1 = 8000"),
        }
        for label, (expected, replacement) in replacements.items():
            with self.subTest(label=label):
                deck = _amr_deck().replace(expected, replacement, 1)
                with self.assertRaises(model.ModelContractError):
                    model.parse_deck_contract(deck)


if __name__ == "__main__":
    unittest.main()
