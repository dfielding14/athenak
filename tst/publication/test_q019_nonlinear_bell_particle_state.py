#!/usr/bin/env python3
"""Adversarial tests for Q019 schema-7 particle-state diagnostics."""

from __future__ import annotations

import copy
import math
import struct
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


def _restart_payload(arguments: dict[str, object] | None = None, **overrides: int) -> bytes:
    values = arguments or _arguments()
    count = len(values["species"])
    real_fields = overrides.get("real_fields", 26)
    integer_fields = overrides.get("integer_fields", 4)
    restart_schema = overrides.get("restart_schema", 7)
    state_kind = overrides.get("state_kind", 1)
    deltaf_mode = overrides.get("deltaf_mode", 0)
    metadata = struct.pack(
        "<15i",
        restart_schema,
        1,
        real_fields,
        integer_fields,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        state_kind,
        0,
    )
    model_ints = [0] * 31
    model_ints[0] = deltaf_mode
    model_reals = [0.0] * 37
    model_reals[20] = values["deposit_qscale"]
    real_rows = np.zeros((count, real_fields), dtype="<f8")
    if real_fields >= 26:
        real_rows[:, [1, 3, 5]] = values["momentum_per_mass"]
        real_rows[:, 6] = values["particle_q_over_mc"]
        real_rows[:, [7, 8, 9]] = values["sampled_magnetic_field"]
        real_rows[:, 22] = values["macro_weight"]
    integer_rows = np.zeros((count, integer_fields), dtype="<i4")
    if integer_fields >= 4:
        integer_rows[:, 0] = 0
        integer_rows[:, 1] = np.arange(count)
        integer_rows[:, 2] = values["species"]
        integer_rows[:, 3] = 0
    return (
        b"<job>\nbasename=q019-schema7-fixture\n<par_end>\n"
        + struct.pack("<Q", reducer.restart_layout.PIC_RESTART_MAGIC)
        + metadata
        + struct.pack("<d", values["artificial_light_speed"])
        + struct.pack("<31i", *model_ints)
        + struct.pack("<37d", *model_reals)
        + struct.pack("<Q", count)
        + struct.pack("<i", count)
        + real_rows.tobytes()
        + integer_rows.tobytes()
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

    def test_raw_schema7_extractor_binds_required_particle_and_model_state(self) -> None:
        arguments = _arguments()
        payload = _restart_payload(arguments)
        extracted = reducer.extract_schema7_particle_state(payload, source="fixture.rst")
        np.testing.assert_array_equal(
            extracted["momentum_per_mass"], arguments["momentum_per_mass"]
        )
        np.testing.assert_array_equal(
            extracted["sampled_magnetic_field"], arguments["sampled_magnetic_field"]
        )
        np.testing.assert_array_equal(extracted["macro_weight"], arguments["macro_weight"])
        np.testing.assert_array_equal(extracted["species"], arguments["species"])
        np.testing.assert_array_equal(
            extracted["particle_q_over_mc"], arguments["particle_q_over_mc"]
        )
        self.assertEqual(extracted["deposit_qscale"], arguments["deposit_qscale"])
        self.assertEqual(
            extracted["artificial_light_speed"], arguments["artificial_light_speed"]
        )
        from_payload = reducer.reduce_schema7_restart_payload(
            payload,
            source="fixture.rst",
            **{
                key: arguments[key]
                for key in (
                    "species_mass",
                    "species_q_over_mc",
                    "domain_volume",
                    "gas_bulk_velocity",
                    "guide_field_direction",
                    "mhd_budget",
                    "reference_budget",
                )
            },
        )
        self.assertEqual(from_payload, _reduce(arguments))

    def test_raw_schema7_extractor_rejects_layout_and_model_drift(self) -> None:
        for overrides, message in (
            ({"restart_schema": 6}, "schema-7 restart probe failed"),
            ({"real_fields": 25}, "real-field count drifted"),
            ({"integer_fields": 5}, "integer-field count drifted"),
            ({"state_kind": 0}, "state is not p/m"),
            ({"deltaf_mode": 1}, "delta-f mode is not off"),
        ):
            with self.subTest(overrides=overrides):
                with self.assertRaisesRegex(reducer.ParticleStateError, message):
                    reducer.extract_schema7_particle_state(
                        _restart_payload(**overrides), source="fixture.rst"
                    )
        with self.assertRaisesRegex(reducer.ParticleStateError, "schema-7 restart probe failed"):
            reducer.extract_schema7_particle_state(_restart_payload()[:-8], source="fixture.rst")

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
