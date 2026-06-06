#!/usr/bin/env python3
"""Adversarial tests for Q019 schema-7 particle-state diagnostics."""

from __future__ import annotations

import copy
import math
import unittest

import numpy as np

from tst.publication import q019_nonlinear_bell_particle_state as reducer


def _arguments() -> dict[str, object]:
    c = 10.0
    momentum = np.asarray(
        [
            [3.0, 4.0, 0.0],
            [0.0, 0.0, 2.0],
        ],
        dtype=np.float64,
    )
    gamma = np.sqrt(1.0 + np.sum(momentum * momentum, axis=1) / c**2)
    velocity = momentum / gamma[:, None]
    qscale = 2.0
    weights = np.asarray([1.0, 2.0], dtype=np.float64)
    masses = np.asarray([1.0, 3.0], dtype=np.float64)
    species = np.asarray([0, 1], dtype=np.int64)
    macro_mass = qscale * weights * masses[species]
    qom = np.asarray([0.5, -0.25], dtype=np.float64)
    particle_momentum = np.sum(macro_mass[:, None] * momentum, axis=0)
    particle_ke = float(np.sum(macro_mass * (gamma - 1.0) * c**2))
    gas_momentum = np.asarray([1.0, -2.0, 3.0])
    gas_energy = 2.0 + 3.0 + 5.0
    return {
        "restart_schema": 7,
        "state_kind": "mass_normalized_momentum",
        "deltaf_mode": "off",
        "momentum_per_mass": momentum,
        "sampled_magnetic_field": np.asarray(
            [[0.0, 0.0, 2.0], [1.0, 0.0, 0.0]], dtype=np.float64
        ),
        "macro_weight": weights,
        "species": species,
        "particle_q_over_mc": qom[species],
        "species_mass": masses,
        "species_q_over_mc": qom,
        "deposit_qscale": qscale,
        "artificial_light_speed": c,
        "domain_volume": 4.0,
        "gas_bulk_velocity": np.asarray([0.25, -0.5, 0.0], dtype=np.float64),
        "guide_field_direction": np.asarray([1.0, 0.0, 0.0], dtype=np.float64),
        "mhd_budget": {
            "gas_momentum": gas_momentum,
            "gas_kinetic_energy": 2.0,
            "gas_thermal_energy": 3.0,
            "magnetic_energy": 5.0,
            "mhd_total_energy": gas_energy,
        },
        "reference_budget": {
            "total_momentum": (gas_momentum + particle_momentum).tolist(),
            "total_energy": gas_energy + particle_ke,
        },
        "_expected_velocity": velocity,
        "_expected_macro_mass": macro_mass,
    }


def _reduce(arguments: dict[str, object] | None = None) -> dict[str, object]:
    values = arguments or _arguments()
    return reducer.reduce_schema7_particle_state(
        **{key: value for key, value in values.items() if not key.startswith("_")}
    )


class Q019NonlinearBellParticleStateTests(unittest.TestCase):
    def test_reconstructs_current_tensor_energy_gyroradius_and_conservation(self) -> None:
        arguments = _arguments()
        result = _reduce(arguments)
        velocity = arguments["_expected_velocity"]
        macro_mass = arguments["_expected_macro_mass"]
        momentum = arguments["momentum_per_mass"]
        volume = arguments["domain_volume"]
        qom = arguments["particle_q_over_mc"]
        c = arguments["artificial_light_speed"]
        gas_velocity = arguments["gas_bulk_velocity"]

        expected_j_over_c = np.sum(
            (macro_mass * qom)[:, None] * velocity, axis=0
        ) / volume
        expected_rho_over_c = float(np.sum(macro_mass * qom) / volume)
        expected_flux = np.einsum("n,ni,nj->ij", macro_mass, momentum, velocity) / volume
        gamma = np.sqrt(1.0 + np.sum(momentum * momentum, axis=1) / c**2)
        expected_ke = float(np.sum(macro_mass * (gamma - 1.0) * c**2))

        np.testing.assert_allclose(
            result["current"]["deposited_current_j_over_c"], expected_j_over_c
        )
        np.testing.assert_allclose(
            result["current"]["reconstructed_lab_current"], c * expected_j_over_c
        )
        np.testing.assert_allclose(
            result["current"]["reconstructed_gas_frame_current"],
            c * (expected_j_over_c - expected_rho_over_c * gas_velocity),
        )
        np.testing.assert_allclose(
            result["particle_momentum"]["momentum_flux_tensor"], expected_flux
        )
        self.assertAlmostEqual(
            result["particle_kinetic_energy"]["volume_integrated"], expected_ke
        )
        self.assertAlmostEqual(result["gyroradius"]["minimum"], 5.0)
        self.assertAlmostEqual(result["gyroradius"]["maximum"], 8.0)
        self.assertAlmostEqual(result["conservation"]["momentum_residual_l2"], 0.0)
        self.assertAlmostEqual(result["conservation"]["energy_residual"], 0.0)
        self.assertFalse(result["launch_authorized"])
        self.assertFalse(result["scientific_claim_authorized"])

    def test_reference_is_optional_but_absence_is_explicit(self) -> None:
        arguments = _arguments()
        arguments["reference_budget"] = None
        result = _reduce(arguments)
        self.assertFalse(result["conservation"]["reference_bound"])
        self.assertNotIn("energy_residual", result["conservation"])

    def test_schema_state_deltaf_and_species_binding_fail_closed(self) -> None:
        for key, value, message in (
            ("restart_schema", 6, "schema 7"),
            ("state_kind", "velocity", "mass-normalized momentum"),
            ("deltaf_mode", "on", "rejects delta-f"),
        ):
            with self.subTest(key=key):
                arguments = _arguments()
                arguments[key] = value
                with self.assertRaisesRegex(reducer.ParticleStateError, message):
                    _reduce(arguments)

        arguments = _arguments()
        arguments["particle_q_over_mc"] = np.asarray([0.5, 0.25], dtype=np.float64)
        with self.assertRaisesRegex(reducer.ParticleStateError, "species table"):
            _reduce(arguments)

    def test_nonfinite_zero_field_and_bad_budget_fail_closed(self) -> None:
        arguments = _arguments()
        arguments["momentum_per_mass"][0, 0] = math.nan
        with self.assertRaisesRegex(reducer.ParticleStateError, "must be finite"):
            _reduce(arguments)

        arguments = _arguments()
        arguments["sampled_magnetic_field"][0] = 0.0
        with self.assertRaisesRegex(reducer.ParticleStateError, "must be nonzero"):
            _reduce(arguments)

        arguments = _arguments()
        arguments["mhd_budget"]["mhd_total_energy"] = 11.0
        with self.assertRaisesRegex(reducer.ParticleStateError, "inconsistent"):
            _reduce(arguments)

    def test_numeric_aliases_and_shape_drift_fail_closed(self) -> None:
        arguments = _arguments()
        arguments["restart_schema"] = 7.0
        with self.assertRaisesRegex(reducer.ParticleStateError, "schema 7"):
            _reduce(arguments)

        arguments = _arguments()
        arguments["macro_weight"] = arguments["macro_weight"].astype(np.int64)
        with self.assertRaisesRegex(reducer.ParticleStateError, "floating dtype"):
            _reduce(arguments)

        arguments = _arguments()
        arguments["species"] = arguments["species"].astype(np.float64)
        with self.assertRaisesRegex(reducer.ParticleStateError, "integer dtype"):
            _reduce(arguments)

        arguments = copy.deepcopy(_arguments())
        arguments["momentum_per_mass"] = arguments["momentum_per_mass"][:, :2]
        with self.assertRaisesRegex(reducer.ParticleStateError, r"shape \(n, 3\)"):
            _reduce(arguments)


if __name__ == "__main__":
    unittest.main()
