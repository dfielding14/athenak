#!/usr/bin/env python3
"""Synthetic tests for the bounded Q-011 restart-continuation tranche."""

from __future__ import annotations

import copy
import math
from pathlib import Path
import struct
import unittest

try:
    from tst.publication import q011_section54_restart as restart
except ModuleNotFoundError:
    import q011_section54_restart as restart


_VALID_LEDGER = {
    "ps_cr_ledger_schema": "3",
    "ps_cr_ledger_complete": "1",
    "ps_mass_reservoir_global": "0.25",
    "ps_injected_cr_count_global": "4.0",
    "ps_injected_cr_mass_global": "4.0",
    "ps_injected_cr_momentum_x1_global": "2.0",
    "ps_injected_cr_momentum_x2_global": "-1.0",
    "ps_injected_cr_momentum_x3_global": "0.5",
    "ps_injected_cr_energy_global": "8.0",
    "ps_removed_excluded_early_cohort": "1",
    "ps_removed_cr_count_global": "2.0",
    "ps_removed_cr_mass_global": "2.0",
    "ps_removed_cr_momentum_x1_global": "1.0",
    "ps_removed_cr_momentum_x2_global": "-0.5",
    "ps_removed_cr_momentum_x3_global": "0.25",
    "ps_removed_cr_energy_global": "3.0",
    "ps_tag_seeded": "1",
    "ps_injection_tag_floor": "10",
    "ps_next_tag": "14",
    "ps_escape_ledger_schema": "1",
    "ps_escape_ledger_complete": "1",
    "ps_escape_audit_calls": "10002",
    "ps_escape_last_audit_time": "500.05",
    "ps_escaped_injected_cr_count_global": "1.0",
    "ps_escaped_injected_cr_mass_global": "1.0",
    "ps_escaped_injected_cr_momentum_x1_global": "0.5",
    "ps_escaped_injected_cr_momentum_x2_global": "-0.25",
    "ps_escaped_injected_cr_momentum_x3_global": "0.125",
    "ps_escaped_injected_cr_energy_global": "2.0",
    "ps_escaped_initial_cr_count_global": "0.0",
}
_CHECKPOINT_OBSERVED_CYCLE = 5001
_CHECKPOINT_OBSERVED_TIME = 500.05


def _restart_payload(
    *,
    restart_schema: int = 7,
    ledger_overrides: dict[str, str | None] | None = None,
    meshblock_particle_counts: list[int] | None = None,
    duplicate_field: tuple[str, str] | None = None,
) -> bytes:
    ledger = dict(_VALID_LEDGER)
    for field, value in (ledger_overrides or {}).items():
        if value is None:
            del ledger[field]
        else:
            ledger[field] = value
    ledger_text = "".join(f"{field}={value}\n" for field, value in ledger.items())
    if duplicate_field is not None:
        ledger_text += f"{duplicate_field[0]}={duplicate_field[1]}\n"
    header = (
        "<job>\n"
        "basename=q011_restart_fixture\n"
        "<problem>\n"
        + ledger_text
        + "<par_end>\n"
    ).encode("ascii")
    meshblock_particle_counts = meshblock_particle_counts or [1]
    meshblock_count = len(meshblock_particle_counts)
    real_fields = 2
    integer_fields = 4
    particle_count = 1
    metadata = struct.pack(
        "<15i",
        restart_schema,
        meshblock_count,
        real_fields,
        integer_fields,
        *([0] * 11),
    )
    model_payload = (
        struct.pack("<d", 500.0)
        + bytes(31 * struct.calcsize("<i"))
        + bytes(37 * struct.calcsize("<d"))
    )
    particle_payload = (
        struct.pack("<Q", particle_count)
        + struct.pack(f"<{meshblock_count}i", *meshblock_particle_counts)
        + struct.pack("<2d", 1.0, 2.0)
        + struct.pack("<4i", 0, 10, 1, 1)
    )
    return (
        header
        + struct.pack("<Q", restart.PIC_RESTART_MAGIC)
        + metadata
        + model_payload
        + particle_payload
    )


def _binding(payload: bytes | None = None) -> dict[str, object]:
    policy = restart.load_preregistration()
    contract = policy["continuation_contract"]
    return restart.bind_checkpoint_for_continuation(
        payload or _restart_payload(),
        checkpoint_nominal_slot_omega0_inverse=contract[
            "checkpoint_nominal_slot_omega0_inverse"
        ],
        checkpoint_observed_committed_cycle=_CHECKPOINT_OBSERVED_CYCLE,
        checkpoint_observed_committed_time_omega0_inverse=_CHECKPOINT_OBSERVED_TIME,
        retained_output_nominal_slots_after_checkpoint_omega0_inverse=contract[
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
        ],
        comparison_tolerances_max_absolute_difference=contract[
            "comparison_tolerances_max_absolute_difference"
        ],
        preregistration=policy,
    )


def _observation(binding: dict[str, object]) -> dict[str, object]:
    outputs = []
    for index, nominal_slot in enumerate(
        binding["retained_output_nominal_slots_after_checkpoint_omega0_inverse"]
    ):
        outputs.append(
            {
                "nominal_slot_omega0_inverse": nominal_slot,
                "observed_committed_cycle": _CHECKPOINT_OBSERVED_CYCLE
                + 100 * (index + 1),
                "observed_committed_time_omega0_inverse": nominal_slot + 0.05,
                "fields": {
                    "rho_bin": [1.0, 2.0 + index],
                    "bmag_bin": [3.0, 4.0 + index],
                    "prtcl_jx_bin": [5.0, 6.0 + index],
                    "j2_bin": [7.0, 8.0 + index],
                    "prtcl_all_pvtk_integer_payload": [0, 10, 1, 1],
                    "prtcl_all_pvtk_float_payload": [1.0, 2.0 + index],
                },
            }
        )
    return {"binding": copy.deepcopy(binding), "outputs_after_checkpoint": outputs}


class Q011Section54RestartPolicyTests(unittest.TestCase):
    def test_checked_in_preregistration_is_exact_and_non_executing(self) -> None:
        policy = restart.load_preregistration()
        self.assertEqual(
            policy["continuation_contract"][
                "checkpoint_nominal_slot_omega0_inverse"
            ],
            500.0,
        )
        self.assertEqual(
            policy["continuation_contract"][
                "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
            ],
            [600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0],
        )
        self.assertEqual(policy["schema_version"], 2)
        self.assertEqual(policy["date"], "2026-06-06")
        self.assertEqual(
            policy["continuation_contract"][
                "comparison_tolerances_max_absolute_difference"
            ]["prtcl_all_pvtk_float_payload"],
            1.0e-6,
        )
        self.assertIn(
            "checkpoint_observed_committed_time_omega0_inverse",
            policy["execution_policy"]["required_before_continuation_execution"],
        )
        self.assertNotIn(
            "paired_output_observed_committed_cycle_and_time",
            policy["execution_policy"]["required_before_continuation_execution"],
        )
        self.assertEqual(
            policy["execution_policy"][
                "required_after_execution_before_parity_result"
            ],
            ["paired_output_observed_committed_cycle_and_time"],
        )
        self.assertFalse(
            policy["execution_policy"]["scheduler_calls_authorized_by_this_record"]
        )
        source = Path(restart.__file__).read_text(encoding="utf-8")
        for forbidden in ("subprocess", "os.system", "sbatch", "srun"):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)

    def test_policy_schema_and_numeric_alias_drift_fail_closed(self) -> None:
        extra = restart.frozen_preregistration()
        extra["unexpected"] = True
        with self.assertRaisesRegex(restart.RestartPolicyError, "schema drift"):
            restart.validate_preregistration(extra)

        alias = restart.frozen_preregistration()
        alias["restart_payload_probe"]["restart_schema"] = 7.0
        with self.assertRaisesRegex(restart.RestartPolicyError, "scalar type drift"):
            restart.validate_preregistration(alias)


class Q011Section54RestartPayloadTests(unittest.TestCase):
    def test_schema7_payload_probe_and_complete_ledgers_extract(self) -> None:
        payload = _restart_payload()
        probe = restart.probe_schema7_restart_payload(payload)
        self.assertEqual(probe.restart_schema, 7)
        self.assertEqual(probe.meshblock_count, 1)
        self.assertEqual(probe.particle_count, 1)
        ledger = restart.extract_startup_shock_ledger(payload)
        self.assertEqual(set(ledger), set(restart.STARTUP_SHOCK_LEDGER_FIELDS))
        self.assertEqual(ledger["ps_cr_ledger_schema"], 3)
        self.assertEqual(ledger["ps_injected_cr_count_global"], 4.0)
        self.assertTrue(ledger["ps_removed_excluded_early_cohort"])
        escape = restart.extract_particle_escape_ledger(payload)
        self.assertEqual(set(escape), set(restart.ESCAPE_LEDGER_FIELDS))
        self.assertEqual(escape["ps_escape_ledger_schema"], 1)
        self.assertEqual(escape["ps_escape_audit_calls"], 10002)
        self.assertEqual(escape["ps_escaped_initial_cr_count_global"], 0.0)

    def test_malformed_restart_schemas_fail_closed(self) -> None:
        for schema in (6, 8):
            with self.subTest(schema=schema):
                with self.assertRaisesRegex(
                    restart.RestartPolicyError, "restart schema is not 7"
                ):
                    restart.probe_schema7_restart_payload(
                        _restart_payload(restart_schema=schema)
                    )

        with self.assertRaisesRegex(
            restart.RestartPolicyError, "truncated particle restart payload"
        ):
            restart.probe_schema7_restart_payload(_restart_payload()[:-1])

        with self.assertRaisesRegex(
            restart.RestartPolicyError, "MeshBlock count table is inconsistent"
        ):
            restart.probe_schema7_restart_payload(
                _restart_payload(meshblock_particle_counts=[0])
            )

    def test_missing_startup_ledger_entries_fail_closed(self) -> None:
        payload = _restart_payload(
            ledger_overrides={"ps_injected_cr_energy_global": None}
        )
        with self.assertRaisesRegex(
            restart.RestartPolicyError,
            "startup shock ledger is missing entries.*ps_injected_cr_energy_global",
        ):
            restart.extract_startup_shock_ledger(payload)

    def test_unrecognized_startup_ledger_boolean_fails_closed(self) -> None:
        payload = _restart_payload(ledger_overrides={"ps_tag_seeded": "yes"})
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "expected 0, 1, false or true"
        ):
            restart.extract_startup_shock_ledger(payload)

    def test_invalid_startup_ledger_tag_progression_fails_closed(self) -> None:
        payload = _restart_payload(ledger_overrides={"ps_next_tag": "15"})
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "invalid next-tag progression"
        ):
            restart.extract_startup_shock_ledger(payload)

    def test_physically_impossible_startup_momentum_energy_fails_closed(self) -> None:
        cases = {
            "ps_injected_cr": {"ps_injected_cr_momentum_x1_global": "1e9"},
            "ps_removed_cr": {"ps_removed_cr_momentum_x1_global": "1e9"},
        }
        for prefix, overrides in cases.items():
            with self.subTest(prefix=prefix):
                with self.assertRaisesRegex(
                    restart.RestartPolicyError,
                    prefix + " energy is below its aggregate momentum lower bound",
                ):
                    restart.extract_startup_shock_ledger(
                        _restart_payload(ledger_overrides=overrides)
                    )

    def test_extreme_finite_startup_momentum_fails_closed(self) -> None:
        cases = {
            "extreme momentum": {"ps_injected_cr_momentum_x1_global": "1e308"},
            "unreliable accumulation count": {
                "ps_injected_cr_count_global": "3000000000000000.0",
                "ps_injected_cr_mass_global": "3000000000000000.0",
            },
        }
        for label, overrides in cases.items():
            with self.subTest(label=label):
                with self.assertRaisesRegex(
                    restart.RestartPolicyError,
                    "ps_injected_cr energy-momentum admissibility bound is invalid",
                ):
                    restart.extract_startup_shock_ledger(
                        _restart_payload(ledger_overrides=overrides)
                    )

    def test_missing_duplicate_and_invalid_escape_ledgers_fail_closed(self) -> None:
        with self.assertRaisesRegex(
            restart.RestartPolicyError,
            "particle escape ledger is missing entries.*ps_escape_audit_calls",
        ):
            restart.extract_particle_escape_ledger(
                _restart_payload(ledger_overrides={"ps_escape_audit_calls": None})
            )

        with self.assertRaisesRegex(
            restart.RestartPolicyError, "duplicate problem parameter"
        ):
            restart.extract_particle_escape_ledger(
                _restart_payload(
                    duplicate_field=("ps_escape_ledger_complete", "0")
                )
            )

        cases = {
            "initial-particle physical escape is unaccounted": {
                "ps_escaped_initial_cr_count_global": "1.0"
            },
            "escaped mass is inconsistent": {
                "ps_escaped_injected_cr_mass_global": "0.5"
            },
            "removed plus escaped count exceeds": {
                "ps_escaped_injected_cr_count_global": "3.0",
                "ps_escaped_injected_cr_mass_global": "3.0",
            },
            "empty escape ledger contains accumulated state": {
                "ps_escaped_injected_cr_count_global": "0.0",
                "ps_escaped_injected_cr_mass_global": "0.0",
            },
            "escaped_injected_cr energy is below": {
                "ps_escaped_injected_cr_momentum_x1_global": "1e9"
            },
        }
        for expected, overrides in cases.items():
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(restart.RestartPolicyError, expected):
                    _binding(_restart_payload(ledger_overrides=overrides))

    def test_escape_chronology_is_bound_to_checkpoint_commit(self) -> None:
        cases = {
            "wrong_calls": {"ps_escape_audit_calls": "10001"},
            "zero_calls": {"ps_escape_audit_calls": "0"},
            "stale_time": {"ps_escape_last_audit_time": "500.0"},
            "one_ulp_time": {
                "ps_escape_last_audit_time": repr(
                    math.nextafter(_CHECKPOINT_OBSERVED_TIME, math.inf)
                )
            },
        }
        for label, overrides in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(restart.RestartPolicyError):
                    _binding(_restart_payload(ledger_overrides=overrides))


class Q011Section54ContinuationParityTests(unittest.TestCase):
    def test_checkpoint_schedule_and_tolerances_bind_before_execution(self) -> None:
        policy = restart.load_preregistration()
        contract = policy["continuation_contract"]
        cases = [
            (
                "checkpoint",
                {
                    "checkpoint_nominal_slot_omega0_inverse": 600.0,
                    "checkpoint_observed_committed_cycle": _CHECKPOINT_OBSERVED_CYCLE,
                    "checkpoint_observed_committed_time_omega0_inverse": (
                        _CHECKPOINT_OBSERVED_TIME
                    ),
                    "retained_output_nominal_slots_after_checkpoint_omega0_inverse": contract[
                        "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
                    ],
                    "comparison_tolerances_max_absolute_difference": contract[
                        "comparison_tolerances_max_absolute_difference"
                    ],
                },
            ),
            (
                "schedule",
                {
                    "checkpoint_nominal_slot_omega0_inverse": contract[
                        "checkpoint_nominal_slot_omega0_inverse"
                    ],
                    "checkpoint_observed_committed_cycle": _CHECKPOINT_OBSERVED_CYCLE,
                    "checkpoint_observed_committed_time_omega0_inverse": (
                        _CHECKPOINT_OBSERVED_TIME
                    ),
                    "retained_output_nominal_slots_after_checkpoint_omega0_inverse": [
                        700.0,
                        800.0,
                        900.0,
                        1000.0,
                        1100.0,
                        1200.0,
                    ],
                    "comparison_tolerances_max_absolute_difference": contract[
                        "comparison_tolerances_max_absolute_difference"
                    ],
                },
            ),
            (
                "tolerances",
                {
                    "checkpoint_nominal_slot_omega0_inverse": contract[
                        "checkpoint_nominal_slot_omega0_inverse"
                    ],
                    "checkpoint_observed_committed_cycle": _CHECKPOINT_OBSERVED_CYCLE,
                    "checkpoint_observed_committed_time_omega0_inverse": (
                        _CHECKPOINT_OBSERVED_TIME
                    ),
                    "retained_output_nominal_slots_after_checkpoint_omega0_inverse": contract[
                        "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
                    ],
                    "comparison_tolerances_max_absolute_difference": {
                        **contract["comparison_tolerances_max_absolute_difference"],
                        "rho_bin": 1.0e-6,
                    },
                },
            ),
        ]
        for label, values in cases:
            with self.subTest(label=label):
                with self.assertRaises(restart.RestartPolicyError):
                    restart.bind_checkpoint_for_continuation(
                        _restart_payload(),
                        preregistration=policy,
                        **values,
                    )

    def test_observed_commit_metadata_must_be_canonical_and_increasing(self) -> None:
        policy = restart.load_preregistration()
        contract = policy["continuation_contract"]
        with self.assertRaisesRegex(restart.RestartPolicyError, "canonical integer"):
            restart.bind_checkpoint_for_continuation(
                _restart_payload(),
                checkpoint_nominal_slot_omega0_inverse=contract[
                    "checkpoint_nominal_slot_omega0_inverse"
                ],
                checkpoint_observed_committed_cycle=True,
                checkpoint_observed_committed_time_omega0_inverse=(
                    _CHECKPOINT_OBSERVED_TIME
                ),
                retained_output_nominal_slots_after_checkpoint_omega0_inverse=contract[
                    "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
                ],
                comparison_tolerances_max_absolute_difference=contract[
                    "comparison_tolerances_max_absolute_difference"
                ],
                preregistration=policy,
            )

        observation = _observation(_binding())
        observation["outputs_after_checkpoint"][2]["observed_committed_time_omega0_inverse"] = (
            observation["outputs_after_checkpoint"][1][
                "observed_committed_time_omega0_inverse"
            ]
        )
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "sequence is not strictly increasing"
        ):
            restart.compare_deterministic_continuation_parity(observation, observation)

    def test_observed_cycle_and_time_parity_mismatches_fail_closed(self) -> None:
        binding = _binding()
        uninterrupted = _observation(binding)

        continued = _observation(binding)
        continued["outputs_after_checkpoint"][1]["observed_committed_cycle"] += 1
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "observed committed cycle parity"
        ):
            restart.compare_deterministic_continuation_parity(
                uninterrupted, continued
            )

        continued = _observation(binding)
        continued["outputs_after_checkpoint"][4][
            "observed_committed_time_omega0_inverse"
        ] += 5.0e-7
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "observed committed time parity"
        ):
            restart.compare_deterministic_continuation_parity(
                uninterrupted, continued
            )

    def test_binding_identity_and_tolerance_failures_fail_closed(self) -> None:
        binding = _binding()
        uninterrupted = _observation(binding)
        continued = _observation(binding)
        continued["binding"]["startup_shock_ledger"][
            "ps_injected_cr_momentum_x1_global"
        ] = 2.5
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "comparison binding identity"
        ):
            restart.compare_deterministic_continuation_parity(
                uninterrupted, continued
            )

        continued = _observation(binding)
        continued["binding"]["particle_escape_ledger"][
            "ps_escaped_injected_cr_momentum_x1_global"
        ] = 0.75
        with self.assertRaisesRegex(
            restart.RestartPolicyError, "comparison binding identity"
        ):
            restart.compare_deterministic_continuation_parity(
                uninterrupted, continued
            )

        continued = _observation(binding)
        continued["outputs_after_checkpoint"][0]["fields"]["rho_bin"][0] += 2.0e-12
        with self.assertRaisesRegex(restart.RestartPolicyError, "tolerance exceeded"):
            restart.compare_deterministic_continuation_parity(
                uninterrupted, continued
            )

    def test_valid_deterministic_continuation_parity(self) -> None:
        binding = _binding()
        uninterrupted = _observation(binding)
        continued = _observation(binding)
        continued["outputs_after_checkpoint"][3]["fields"][
            "prtcl_all_pvtk_float_payload"
        ][1] += 5.0e-7
        result = restart.compare_deterministic_continuation_parity(
            uninterrupted, continued
        )
        self.assertEqual(result["result"], "pass_deterministic_continuation_parity")
        self.assertEqual(result["checkpoint_nominal_slot_omega0_inverse"], 500.0)
        self.assertEqual(
            result["checkpoint_observed_committed_cycle"],
            _CHECKPOINT_OBSERVED_CYCLE,
        )
        self.assertEqual(
            result["checkpoint_observed_committed_time_omega0_inverse"],
            _CHECKPOINT_OBSERVED_TIME,
        )
        self.assertEqual(
            result["paired_output_observed_commits"][0],
            {
                "nominal_slot_omega0_inverse": 600.0,
                "observed_committed_cycle": 5101,
                "observed_committed_time_omega0_inverse": 600.05,
            },
        )
        self.assertAlmostEqual(
            result["maximum_absolute_difference_by_field"][
                "prtcl_all_pvtk_float_payload"
            ],
            5.0e-7,
        )


if __name__ == "__main__":
    unittest.main()
