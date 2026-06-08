#!/usr/bin/env python3
"""Adversarial tests for the Q011 exact conservation-closure successor."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest import mock

import numpy as np

from tst.publication import q011_section54_exact_conservation_closure_successor_v1 as reducer
from tst.publication import q011_section54_restart as restart_layout
from tst.scripts.particles import (
    pic_parallel_shock_exact_conservation_closure_vl2_tsc as runtime,
)


ROOT = Path(__file__).resolve().parents[2]
DECK = (
    ROOT
    / "inputs/publication/"
    "pic_parallel_shock_section54_exact_conservation_successor_v1_vl2_tsc.athinput"
)
READINESS = (
    ROOT
    / "tst/publication/readiness/"
    "q011_section54_exact_conservation_numerical_repair_successor_v3_2026-06-06.json"
)
QSCALE = 9.0e-4
LIGHT_SPEED = 10000.0
INITIAL_MHD = [10.0, 0.0, 0.0, 0.0, 100.0]
CR_STATE = [QSCALE, 0.0, 0.0, 0.0, 0.0]


def _fixture_deck_payload() -> bytes:
    text = DECK.read_text(encoding="ascii")
    replacements = (
        ("nx1       = 4000", "nx1       = 2"),
        ("x1max     = 48000.0", "x1max     = 2.0"),
        ("nx2       = 260", "nx2       = 2"),
        ("x2max     = 3120.0", "x2max     = 2.0"),
        ("nx1       = 20", "nx1       = 2"),
        ("nx2       = 20", "nx2       = 2"),
        ("refinement           = adaptive", "refinement           = none"),
        ("num_levels           = 3", "num_levels           = 1"),
        ("tlim       = 1200.0", "tlim       = 200.0"),
        ("ps_enable_curvature_amr       = true",
         "ps_enable_curvature_amr       = false"),
    )
    for old, new in replacements:
        if text.count(old) != 1:
            raise AssertionError("fixture deck replacement drifted: " + old)
        text = text.replace(old, new, 1)
    return text.encode("ascii")


FIXTURE_DECK = _fixture_deck_payload()


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _artifact(name: str, payload: bytes, lineage: dict[str, str]) -> dict[str, object]:
    return {
        "name": name,
        "payload": payload,
        "sha256": _sha(payload),
        "byte_count": len(payload),
        "lineage": copy.deepcopy(lineage),
    }


def _history(labels: tuple[str, ...], rows: list[list[float]]) -> bytes:
    header = "#  " + " ".join(
        f"[{index}]={label}" for index, label in enumerate(labels, start=1)
    )
    body = [" ".join(f"{value:.16e}" for value in row) for row in rows]
    return ("# Athena++ history data\n" + header + "\n" + "\n".join(body) + "\n").encode(
        "ascii"
    )


def _mhd_rows(*, mass_offset_t200: float = 0.0) -> list[list[float]]:
    initial = [0.0, 1.0, *INITIAL_MHD, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    checkpoint = [
        100.0,
        1.0,
        INITIAL_MHD[0] - CR_STATE[0],
        0.0,
        0.0,
        0.0,
        INITIAL_MHD[4],
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    terminal = list(checkpoint)
    terminal[0] = 200.0
    terminal[2] += mass_offset_t200
    return [initial, checkpoint, terminal]


def _user_rows(
    *, ext_mass_t100: float = 0.0, ext_mass_t200: float = 0.0
) -> list[list[float]]:
    rows = []
    for time in (0.0, 100.0, 200.0):
        values = [time, 1.0] + [0.0] * 12
        if time == 100.0:
            values[2] = ext_mass_t100
        elif time == 200.0:
            values[2] = ext_mass_t200
        rows.append(values)
    return rows


def _ledger(
    *,
    time: float,
    cycle: int,
    overrides: dict[str, str | None] | None = None,
) -> dict[str, str]:
    values = {
        "ps_cr_ledger_schema": "3",
        "ps_cr_ledger_complete": "true",
        "ps_mass_reservoir_global": "0.0",
        "ps_injected_cr_count_global": "1.0",
        "ps_injected_cr_mass_global": repr(QSCALE),
        "ps_injected_cr_momentum_x1_global": "0.0",
        "ps_injected_cr_momentum_x2_global": "0.0",
        "ps_injected_cr_momentum_x3_global": "0.0",
        "ps_injected_cr_energy_global": "0.0",
        "ps_removed_excluded_early_cohort": "true",
        "ps_removed_cr_count_global": "0.0",
        "ps_removed_cr_mass_global": "0.0",
        "ps_removed_cr_momentum_x1_global": "0.0",
        "ps_removed_cr_momentum_x2_global": "0.0",
        "ps_removed_cr_momentum_x3_global": "0.0",
        "ps_removed_cr_energy_global": "0.0",
        "ps_tag_seeded": "true",
        "ps_injection_tag_floor": "10",
        "ps_next_tag": "11",
        "ps_escape_ledger_schema": "1",
        "ps_escape_ledger_complete": "true",
        "ps_escape_audit_calls": str(2 * cycle),
        "ps_escape_last_audit_time": repr(time),
        "ps_escaped_injected_cr_count_global": "0.0",
        "ps_escaped_injected_cr_mass_global": "0.0",
        "ps_escaped_injected_cr_momentum_x1_global": "0.0",
        "ps_escaped_injected_cr_momentum_x2_global": "0.0",
        "ps_escaped_injected_cr_momentum_x3_global": "0.0",
        "ps_escaped_injected_cr_energy_global": "0.0",
        "ps_escaped_initial_cr_count_global": "0.0",
        "ps_conservation_ledger_schema": "1",
        "ps_conservation_ledger_complete": "true",
        "ps_conservation_committed_cycles": str(cycle),
        "ps_conservation_committed_time": repr(time),
    }
    for prefix in reducer._CONS_PREFIXES:
        for component in reducer._COMPONENTS:
            values[f"{prefix}_{component}_global"] = "0.0"
    values["ps_cons_gas_subtracted_mass_global"] = repr(QSCALE)
    for key, value in (overrides or {}).items():
        if value is None:
            values.pop(key, None)
        else:
            values[key] = value
    return values


def _restart_payload(
    *,
    time: float,
    cycle: int,
    deck_payload: bytes | None = None,
    mhd_state: list[float] | None = None,
    ledger_overrides: dict[str, str | None] | None = None,
    duplicate_problem_line: str | None = None,
) -> bytes:
    ledger = _ledger(time=time, cycle=cycle, overrides=ledger_overrides)
    parameter_lines = "".join(f"{key}={value}\n" for key, value in ledger.items())
    if duplicate_problem_line is not None:
        parameter_lines += duplicate_problem_line + "\n"
    deck_text = (deck_payload or FIXTURE_DECK).decode("ascii")
    output_marker = "\n<output1>\n"
    if output_marker not in deck_text:
        raise AssertionError("fixture deck lacks <output1> marker")
    header = (
        deck_text.replace(
            output_marker, "\n" + parameter_lines + output_marker, 1
        )
        + "\n<par_end>\n"
    ).encode("ascii")

    mesh_header = (
        struct.pack("<ii", 1, 0)
        + struct.pack("<9d", 0.0, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0)
        + struct.pack("<19i", 2, 2, 2, 1, 2, 3, 2, 3, 0, 0, *([0] * 9))
        + struct.pack("<19i", 2, 2, 2, 1, 2, 3, 2, 3, 0, 0, *([0] * 9))
        + struct.pack("<ddii", time, 1.0, cycle, 1)
    )
    state = np.asarray(
        mhd_state or [INITIAL_MHD[0] - CR_STATE[0], 0.0, 0.0, 0.0, INITIAL_MHD[4]],
        dtype="<f8",
    )
    mhd = np.zeros((5, 1, 6, 6), dtype="<f8")
    mhd[:, 0, 2:4, 2:4] = state[:, None, None] / 4.0
    face_fields = np.zeros(42 + 42 + 72, dtype="<f8")
    mhd_payload = mhd.tobytes() + face_fields.tobytes()
    mesh_layout = (
        struct.pack("<4i", 0, 0, 0, 0)
        + struct.pack("<f", 1.0)
        + struct.pack("<i", 0)
        + struct.pack("<i", 0)
        + struct.pack("<i", 1)
        + struct.pack("<QiiQ", reducer._MESH_METADATA_MAGIC, 2, 0, 1)
        + struct.pack("<Q", len(mhd_payload))
    )
    metadata = struct.pack(
        "<15i",
        7,
        1,
        26,
        4,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
    )
    model_ints = [0] * 31
    model_reals = [0.0] * 37
    model_reals[20] = QSCALE
    real_row = np.zeros((1, 26), dtype="<f8")
    real_row[0, 6] = 1.0
    real_row[0, 7] = 1.0
    real_row[0, 22] = 1.0
    integer_row = np.asarray([[0, 10, 0, 1]], dtype="<i4")
    particle_payload = (
        struct.pack("<Q", restart_layout.PIC_RESTART_MAGIC)
        + metadata
        + struct.pack("<d", LIGHT_SPEED)
        + struct.pack("<31i", *model_ints)
        + struct.pack("<37d", *model_reals)
        + struct.pack("<Q", 1)
        + struct.pack("<i", 1)
        + real_row.tobytes()
        + integer_row.tobytes()
    )
    return header + mesh_header + mesh_layout + mhd_payload + particle_payload


def _fixture(
    *,
    deck_payload: bytes | None = None,
    mhd_payload: bytes | None = None,
    user_payload: bytes | None = None,
    restart_payloads: list[bytes] | None = None,
) -> dict[str, object]:
    candidate = b"q011 exact conservation candidate fixture\n"
    deck = deck_payload if deck_payload is not None else FIXTURE_DECK
    attempt = b'{"attempt_id":"q011-exact-conservation-fixture"}\n'
    mhd = mhd_payload or _history(reducer._MHD_LABELS, _mhd_rows())
    user = user_payload or _history(reducer._USER_LABELS, _user_rows())
    restarts = restart_payloads or [
        _restart_payload(time=100.0, cycle=10, deck_payload=deck),
        _restart_payload(time=200.0, cycle=20, deck_payload=deck),
    ]
    lineage = {
        "attempt_id": "q011-exact-conservation-fixture",
        "attempt_sha256": _sha(attempt),
        "candidate_commit": "a" * 40,
        "candidate_sha256": _sha(candidate),
        "deck_sha256": _sha(deck),
        "executable_sha256": "e" * 64,
    }
    return {
        "lineage": lineage,
        "candidate_artifact": _artifact("candidate.tar", candidate, lineage),
        "deck_artifact": _artifact(DECK.name, deck, lineage),
        "attempt_artifact": _artifact("attempt.json", attempt, lineage),
        "mhd_history_artifact": _artifact("fixture.mhd.hst", mhd, lineage),
        "user_history_artifact": _artifact("fixture.user.hst", user, lineage),
        "restart_artifacts": [
            _artifact(f"fixture.{index:05d}.rst", payload, lineage)
            for index, payload in enumerate(restarts, start=1)
        ],
    }


def _reduce(fixture: dict[str, object] | None = None) -> dict[str, object]:
    return reducer.reduce_exact_conservation_closure(**(fixture or _fixture()))


class Q011ExactConservationClosureTests(unittest.TestCase):
    def test_valid_fixture_binds_exact_closure_and_grants_no_authority(self) -> None:
        result = _reduce()
        self.assertEqual(result["record_type"], reducer.RECORD_TYPE)
        self.assertEqual(len(result["checkpoints"]), 2)
        closure = result["checkpoints"][-1]["conservation_closure"]
        self.assertLess(closure["absolute_residual"]["mass"], 1.0e-15)
        self.assertEqual(closure["absolute_residual"]["energy"], 0.0)
        self.assertFalse(result["checkpoints"][-1]["instantaneous_partitions"][
            "is_conservation_closure"
        ])
        self.assertFalse(closure["is_instantaneous_partition"])
        self.assertEqual(
            closure["signed_external_contributions"]["ps_cons_gas_subtracted"]["mass"],
            -QSCALE,
        )
        self.assertEqual(
            closure["cumulative_ledger_terms"]["ps_cons_gas_subtracted"]["mass"],
            QSCALE,
        )
        self.assertIn("cumulative_transport_l1", closure["denominators"])
        self.assertIn("max_state_abs", closure["denominators"])
        self.assertEqual(
            result["checkpoint_timing_contract"]["interpretation"],
            "closure_uses_restart_state_and_separately_reports_mhd_state_change_after_history",
        )
        self.assertEqual(
            result["status"],
            "supplied_artifact_self_consistency_reduction_only_not_execution_provenance",
        )
        self.assertTrue(
            any(
                "does not establish build, scheduler, command" in limitation
                for limitation in result["limitations"]
            )
        )
        self.assertEqual(
            result["history_publication_audit"]["mhd"][
                "identical_duplicate_times_collapsed"
            ],
            0,
        )
        self.assertTrue(all(value is False for value in result["authority"].values()))
        self.assertTrue(reducer.canonical_json_bytes(result).endswith(b"\n"))

    def test_nonzero_residual_is_reported_not_reclassified_as_partition(self) -> None:
        post_amr = list(INITIAL_MHD)
        post_amr[0] = INITIAL_MHD[0] - CR_STATE[0] + 1.0
        result = _reduce(
            _fixture(
                restart_payloads=[
                    _restart_payload(time=100.0, cycle=10),
                    _restart_payload(time=200.0, cycle=20, mhd_state=post_amr),
                ]
            )
        )
        closure = result["checkpoints"][-1]["conservation_closure"]
        self.assertAlmostEqual(closure["signed_residual"]["mass"], 1.0)
        self.assertAlmostEqual(closure["absolute_residual"]["mass"], 1.0)
        self.assertAlmostEqual(
            closure["mhd_state_change_after_history_before_restart"]["mass"], 1.0
        )
        self.assertAlmostEqual(
            closure["history_signed_residual"]["mass"], 0.0
        )
        self.assertGreater(
            closure["absolute_normalized_residual"]["max_state_abs"]["mass"], 0.0
        )
        self.assertFalse(closure["is_instantaneous_partition"])

    def test_identical_restart_history_publication_duplicates_are_audited(self) -> None:
        mhd = _history(reducer._MHD_LABELS, _mhd_rows())
        user = _history(reducer._USER_LABELS, _user_rows())
        mhd_restarted = mhd + _history(reducer._MHD_LABELS, [_mhd_rows()[-1]])
        user_restarted = user + _history(reducer._USER_LABELS, [_user_rows()[-1]])
        result = _reduce(_fixture(mhd_payload=mhd_restarted, user_payload=user_restarted))
        for physics in ("mhd", "user"):
            audit = result["history_publication_audit"][physics]
            self.assertEqual(audit["segment_headers"], 2)
            self.assertEqual(audit["identical_duplicate_times_collapsed"], 1)

    def test_history_missing_drifted_duplicate_nonmonotonic_and_grid_drift_fail_closed(
        self,
    ) -> None:
        valid = _history(reducer._MHD_LABELS, _mhd_rows())
        duplicate_header = valid + valid.splitlines(keepends=True)[1]
        duplicate_drift_rows = _mhd_rows()
        duplicate_drift = list(duplicate_drift_rows[-1])
        duplicate_drift[2] += 1.0
        duplicate_drift_rows.append(duplicate_drift)
        nonmonotonic_rows = _mhd_rows()
        nonmonotonic = list(nonmonotonic_rows[-1])
        nonmonotonic[0] = 150.0
        nonmonotonic_rows.append(nonmonotonic)
        missing_checkpoint = _history(
            reducer._MHD_LABELS, [_mhd_rows()[0], _mhd_rows()[2]]
        )
        user_missing_checkpoint = _history(
            reducer._USER_LABELS, [_user_rows()[0], _user_rows()[2]]
        )
        cases = (
            (_fixture(mhd_payload=duplicate_header), "history segment header lacks banner"),
            (
                _fixture(
                    mhd_payload=_history(reducer._MHD_LABELS, duplicate_drift_rows)
                ),
                "duplicate time has state drift",
            ),
            (
                _fixture(mhd_payload=_history(reducer._MHD_LABELS, nonmonotonic_rows)),
                "nonmonotonic time",
            ),
            (
                _fixture(
                    mhd_payload=missing_checkpoint,
                    user_payload=user_missing_checkpoint,
                ),
                "history: does not contain the complete configured checkpoint sequence",
            ),
            (
                _fixture(user_payload=user_missing_checkpoint),
                "MHD and user history time grids differ",
            ),
        )
        for fixture, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    _reduce(fixture)

    def test_restart_commit_sequence_schema_and_duplicate_ledgers_fail_closed(self) -> None:
        duplicate = _restart_payload(
            time=100.0,
            cycle=10,
            duplicate_problem_line="ps_conservation_ledger_schema=1",
        )
        cases = (
            (
                [
                    _restart_payload(time=100.0, cycle=10),
                    _restart_payload(time=100.0, cycle=10),
                ],
                "restart sequence is duplicate or nonmonotonic",
            ),
            (
                [
                    _restart_payload(
                        time=100.0,
                        cycle=10,
                        ledger_overrides={"ps_conservation_committed_cycles": "11"},
                    ),
                    _restart_payload(time=200.0, cycle=20),
                ],
                "restart header/ledger commit discontinuity",
            ),
            (
                [
                    _restart_payload(
                        time=100.0,
                        cycle=10,
                        ledger_overrides={"ps_conservation_ledger_complete": "false"},
                    ),
                    _restart_payload(time=200.0, cycle=20),
                ],
                "conservation ledger is incomplete",
            ),
            (
                [
                    _restart_payload(
                        time=100.0,
                        cycle=10,
                        ledger_overrides={"ps_cons_mhd_boundary_energy_global": None},
                    ),
                    _restart_payload(time=200.0, cycle=20),
                ],
                "conservation ledger is missing",
            ),
            (
                [
                    _restart_payload(
                        time=100.0,
                        cycle=10,
                        ledger_overrides={
                            "ps_cons_particle_escape_energy_global": "1.0"
                        },
                    ),
                    _restart_payload(time=200.0, cycle=20),
                ],
                "particle-escape mass or energy delta is positive",
            ),
            (
                [
                    _restart_payload(
                        time=100.0,
                        cycle=10,
                        ledger_overrides={
                            "ps_cons_particle_reflect_energy_global": "1.0"
                        },
                    ),
                    _restart_payload(time=200.0, cycle=20),
                ],
                "reflecting-wall mass or energy delta is nonzero",
            ),
            ([duplicate, _restart_payload(time=200.0, cycle=20)], "duplicate or empty parameter"),
        )
        for restarts, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    _reduce(_fixture(restart_payloads=restarts))

    def test_sparse_restart_sequence_and_escape_ledger_forgery_fail_closed(self) -> None:
        deck300 = FIXTURE_DECK.replace(b"tlim       = 200.0", b"tlim       = 300.0")
        self.assertNotEqual(deck300, FIXTURE_DECK)
        mhd_rows = _mhd_rows()
        mhd_terminal = list(mhd_rows[-1])
        mhd_terminal[0] = 300.0
        mhd_rows.append(mhd_terminal)
        user_rows = _user_rows()
        user_terminal = list(user_rows[-1])
        user_terminal[0] = 300.0
        user_rows.append(user_terminal)
        with self.assertRaisesRegex(
            reducer.ConservationClosureError,
            "restart artifacts: does not contain the complete configured checkpoint sequence",
        ):
            _reduce(
                _fixture(
                    deck_payload=deck300,
                    mhd_payload=_history(reducer._MHD_LABELS, mhd_rows),
                    user_payload=_history(reducer._USER_LABELS, user_rows),
                    restart_payloads=[
                        _restart_payload(time=100.0, cycle=10, deck_payload=deck300),
                        _restart_payload(time=300.0, cycle=30, deck_payload=deck300),
                    ],
                )
            )

        cases = (
            (
                {"ps_escape_ledger_complete": "false"},
                "physical escape ledger is incomplete",
            ),
            (
                {"ps_escape_audit_calls": "1"},
                "physical escape audit-call chronology is incomplete",
            ),
            (
                {"ps_escape_last_audit_time": "0.0"},
                "physical escape audit time differs from committed time",
            ),
            (
                {"ps_escaped_initial_cr_count_global": "1.0"},
                "initially seeded CR escape is not fully accounted",
            ),
            (
                {
                    "ps_escaped_injected_cr_count_global": "1.0",
                    "ps_escaped_injected_cr_mass_global": repr(QSCALE),
                },
                "generic and reason-coded physical escape differ at mass",
            ),
            (
                {"ps_escape_ledger_schema": None},
                "physical escape ledger is missing",
            ),
        )
        for overrides, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    _reduce(
                        _fixture(
                            restart_payloads=[
                                _restart_payload(
                                    time=100.0,
                                    cycle=10,
                                    ledger_overrides=overrides,
                                ),
                                _restart_payload(time=200.0, cycle=20),
                            ]
                        )
                    )

        duplicate_escape = _restart_payload(
            time=100.0,
            cycle=10,
            duplicate_problem_line="ps_escape_ledger_complete=false",
        )
        with self.assertRaisesRegex(
            reducer.ConservationClosureError, "duplicate or empty parameter"
        ):
            _reduce(
                _fixture(
                    restart_payloads=[
                        duplicate_escape,
                        _restart_payload(time=200.0, cycle=20),
                    ]
                )
            )

    def test_restart_ledger_regression_and_history_mismatch_fail_closed(self) -> None:
        regressed = [
            _restart_payload(
                time=100.0,
                cycle=10,
                ledger_overrides={"ps_cons_gas_subtracted_mass_global": "0.0018"},
            ),
            _restart_payload(time=200.0, cycle=20),
        ]
        with self.assertRaisesRegex(
            reducer.ConservationClosureError, "gas-subtraction ledger regressed"
        ):
            _reduce(
                _fixture(
                    restart_payloads=regressed,
                    user_payload=_history(
                        reducer._USER_LABELS,
                        _user_rows(ext_mass_t100=-QSCALE),
                    ),
                )
            )

        user = _history(reducer._USER_LABELS, _user_rows(ext_mass_t200=1.0))
        with self.assertRaisesRegex(
            reducer.ConservationClosureError,
            "user-history external delta differs from restart ledger",
        ):
            _reduce(_fixture(user_payload=user))

    def test_restart_effective_input_override_drift_fails_closed(self) -> None:
        restart = _restart_payload(time=100.0, cycle=10).replace(
            b"ps_p0                         = 1.0",
            b"ps_p0                         = 0.5",
            1,
        )
        with self.assertRaisesRegex(
            reducer.ConservationClosureError,
            r"restart/deck control drift at <problem>/ps_p0",
        ):
            _reduce(
                _fixture(
                    restart_payloads=[restart, _restart_payload(time=200.0, cycle=20)]
                )
            )

    def test_lineage_hash_byte_count_and_artifact_name_drift_fail_closed(self) -> None:
        cases: list[tuple[dict[str, object], str]] = []
        bad_hash = _fixture()
        bad_hash["mhd_history_artifact"]["sha256"] = "0" * 64
        cases.append((bad_hash, "sha256 mismatch"))
        bad_bytes = _fixture()
        bad_bytes["user_history_artifact"]["byte_count"] += 1
        cases.append((bad_bytes, "byte_count mismatch"))
        bad_lineage = _fixture()
        bad_lineage["restart_artifacts"][0]["lineage"]["attempt_id"] = "mixed"
        cases.append((bad_lineage, "lineage differs"))
        bad_name = _fixture()
        bad_name["mhd_history_artifact"]["name"] = "fixture.hst"
        cases.append((bad_name, "MHD history name drifted"))
        for fixture, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    _reduce(fixture)

    def test_deck_model_and_boundary_drift_fail_closed(self) -> None:
        deck = FIXTURE_DECK.decode("ascii")
        cases = (
            (deck.replace("ps_enable_frame_tracking      = false",
                          "ps_enable_frame_tracking      = true"),
             "ps_enable_frame_tracking does not equal false"),
            (deck.replace("ox1_bc    = inflow", "ox1_bc    = periodic"),
             "ox1_bc is not an exact-ledger physical boundary"),
            (deck.replace("ps_p0                         = 1.0",
                          "ps_p0                         = 0.1"),
             "ps_p0 does not equal 1.0"),
            (deck.replace("gamma       = 1.66666666667",
                          "gamma       = 1.66666666667\nconst_accel = true"),
             "untracked MHD source"),
            (
                deck.replace(
                    "refinement_interval  = 1",
                    "refinement_interval  = 1\nprolong_primitives  = true",
                ),
                "non-conservative primitive-variable AMR prolongation",
            ),
            (
                deck.replace("refinement           = none",
                             "refinement           = adaptive"),
                "fixed uniform mesh",
            ),
            (
                deck.replace("num_levels           = 1", "num_levels           = 2"),
                "one mesh level",
            ),
            (
                deck.replace(
                    "pic_load_balance_cost_per_particle = 0.0",
                    "pic_load_balance_cost_per_particle = 0.25",
                ),
                "pic_load_balance_cost_per_particle does not equal 0.0",
            ),
            (
                deck.replace("ps_enable_curvature_amr       = false",
                             "ps_enable_curvature_amr       = true"),
                "ps_enable_curvature_amr does not equal false",
            ),
        )
        for text, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    _reduce(_fixture(deck_payload=text.encode("ascii")))

        outflow = FIXTURE_DECK.replace(b"ox1_bc    = inflow", b"ox1_bc    = outflow")
        result = _reduce(_fixture(deck_payload=outflow))
        self.assertEqual(len(result["checkpoints"]), 2)

    def test_reducer_rejects_one_dimensional_deck_labeled_2d3v(self) -> None:
        for axis in (b"nx1       = 2", b"nx2       = 2"):
            with self.subTest(axis=axis.decode("ascii").strip()):
                one_dimensional = FIXTURE_DECK.replace(axis, axis[:-1] + b"1")
                with self.assertRaisesRegex(
                    reducer.ConservationClosureError,
                    "requires true 2D x1-x2 geometry",
                ):
                    _reduce(_fixture(deck_payload=one_dimensional))

    def test_restart_mesh_partition_metadata_payload_and_state_drift_fail_closed(
        self,
    ) -> None:
        valid = _restart_payload(time=100.0, cycle=10)
        header_end = valid.index(b"<par_end>\n") + len(b"<par_end>\n")
        fixed_header_bytes = 2 * 4 + 9 * 8 + 2 * 19 * 4 + 2 * 8 + 2 * 4
        logical_location = header_end + fixed_header_bytes
        rank_eachmb = logical_location + 4 * 4 + 4
        metadata = rank_eachmb + 4 + 4 + 4
        data_size = metadata + struct.calcsize("<QiiQ")
        data = data_size + struct.calcsize("<Q")
        pic_magic = struct.pack("<Q", restart_layout.PIC_RESTART_MAGIC)

        bad_rank = bytearray(valid)
        struct.pack_into("<i", bad_rank, rank_eachmb, 1)
        bad_level = bytearray(valid)
        struct.pack_into("<i", bad_level, logical_location + 3 * 4, -1)
        refined_level = bytearray(valid)
        struct.pack_into("<i", refined_level, logical_location + 3 * 4, 1)
        out_of_range_location = bytearray(valid)
        struct.pack_into("<i", out_of_range_location, logical_location, 1)
        nonunit_cost = bytearray(valid)
        struct.pack_into("<f", nonunit_cost, logical_location + 4 * 4, 2.0)
        bad_metadata = bytearray(valid)
        struct.pack_into("<Q", bad_metadata, metadata, 0)
        adaptive_metadata = bytearray(valid)
        struct.pack_into("<i", adaptive_metadata, metadata + 12, 1)
        adaptive_metadata[metadata + 24:metadata + 24] = struct.pack("<i", 0)
        bad_size = bytearray(valid)
        struct.pack_into("<Q", bad_size, data_size, struct.unpack_from("<Q", valid, data_size)[0] + 8)
        bad_boundary = valid.replace(pic_magic, b"\x00" * 8 + pic_magic, 1)
        nonfinite = bytearray(valid)
        struct.pack_into("<d", nonfinite, data + (2 * 6 + 2) * 8, math.nan)

        cases = (
            (bytes(bad_rank), "MeshBlock rank assignment drifted"),
            (bytes(bad_level), "invalid MeshBlock level"),
            (bytes(refined_level), "rejects refined restart topology"),
            (bytes(out_of_range_location), "logical coordinate is out of range"),
            (bytes(nonunit_cost), "rejects non-unit load-balance costs"),
            (bytes(bad_metadata), "mesh metadata schema drifted"),
            (bytes(adaptive_metadata), "rejects adaptive restart metadata"),
            (bytes(bad_size), "untracked or malformed physics payload"),
            (bad_boundary, "MHD/PIC payload boundary drifted"),
            (bytes(nonfinite), "nonfinite restart MHD state"),
        )
        terminal = _restart_payload(time=200.0, cycle=20)
        for payload, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    _reduce(_fixture(restart_payloads=[payload, terminal]))

    def test_fixed_uniform_topology_requires_in_range_unique_complete_tiling(
        self,
    ) -> None:
        valid = ((0, 0, 0, 1), (1, 0, 0, 1), (0, 1, 0, 1), (1, 1, 0, 1))
        reducer._validate_fixed_uniform_topology(
            logical_locations=valid,
            deck_mesh_cells=(4, 2, 1),
            deck_block_cells=(2, 1, 1),
            root_level=1,
            nmb_total=4,
            label="fixture",
        )
        cases = (
            (
                ((0, 0, 0, 1), (1, 0, 0, 1), (0, 1, 0, 1), (2, 1, 0, 1)),
                (4, 2, 1),
                (2, 1, 1),
                1,
                4,
                "logical coordinate is out of range",
            ),
            (
                ((0, 0, 0, 1), (1, 0, 0, 1), (0, 1, 0, 1), (0, 1, 0, 1)),
                (4, 2, 1),
                (2, 1, 1),
                1,
                4,
                "logical coordinate is duplicated",
            ),
            (
                valid[:-1],
                (4, 2, 1),
                (2, 1, 1),
                1,
                4,
                "MeshBlock count is incomplete",
            ),
            (
                valid,
                (3, 2, 1),
                (2, 1, 1),
                1,
                4,
                "MeshBlocks do not tile",
            ),
            (
                valid,
                (4, 2, 1),
                (2, 1, 1),
                0,
                4,
                "root level is inconsistent",
            ),
        )
        for locations, mesh, block, level, count, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(reducer.ConservationClosureError, message):
                    reducer._validate_fixed_uniform_topology(
                        logical_locations=locations,
                        deck_mesh_cells=mesh,
                        deck_block_cells=block,
                        root_level=level,
                        nmb_total=count,
                        label="fixture",
                    )

    def test_nonfinite_history_and_particle_schema_drift_fail_closed(self) -> None:
        mhd = _history(reducer._MHD_LABELS, _mhd_rows()).replace(
            b"1.0000000000000000e+02", b"nan", 1
        )
        with self.assertRaisesRegex(reducer.ConservationClosureError, "noncanonical real"):
            _reduce(_fixture(mhd_payload=mhd))

        restart = bytearray(_restart_payload(time=100.0, cycle=10))
        marker = struct.pack("<Q", restart_layout.PIC_RESTART_MAGIC)
        offset = restart.find(marker) + len(marker)
        struct.pack_into("<i", restart, offset, 6)
        with self.assertRaisesRegex(reducer.ConservationClosureError, "schema-7 probe failed"):
            _reduce(
                _fixture(
                    restart_payloads=[bytes(restart), _restart_payload(time=200.0, cycle=20)]
                )
            )

    def test_runtime_restart_particle_state_comparison_fails_closed(self) -> None:
        valid = _restart_payload(time=100.0, cycle=10)
        comparison = runtime._compare_particle_state(valid, valid, "fixture")
        self.assertEqual(comparison["particle_count"], 1)
        self.assertTrue(comparison["tag_sorted_integer_state_exactly_equal"])
        self.assertEqual(comparison["tag_sorted_real_state_max_absolute_error"], 0.0)

        probe = restart_layout.probe_schema7_restart_payload(valid, source="fixture")
        real_drift = bytearray(valid)
        struct.pack_into("<d", real_drift, probe.particle_real_offset, 1.0)
        integer_drift = bytearray(valid)
        struct.pack_into("<i", integer_drift, probe.particle_integer_offset + 4, 11)
        nonfinite = bytearray(valid)
        struct.pack_into("<d", nonfinite, probe.particle_real_offset, math.nan)
        for payload, message in (
            (bytes(real_drift), "restart particle real state drifted"),
            (bytes(integer_drift), "restart particle integer state drifted"),
            (bytes(nonfinite), "restart particle real state is nonfinite"),
        ):
            with self.subTest(message=message):
                with self.assertRaisesRegex(RuntimeError, message):
                    runtime._compare_particle_state(valid, payload, "fixture")

    def test_spatial_restart_digest_covers_face_fields(self) -> None:
        valid = _restart_payload(time=100.0, cycle=10)
        deck = reducer._deck_physics(FIXTURE_DECK)
        header, aggregate = reducer._restart_mhd_state(valid, deck, "valid")
        modified = bytearray(valid)
        pic_offset = modified.find(struct.pack("<Q", restart_layout.PIC_RESTART_MAGIC))
        self.assertGreater(pic_offset, 0)
        modified[pic_offset - 1] ^= 1
        changed_header, changed_aggregate = reducer._restart_mhd_state(
            bytes(modified), deck, "modified"
        )
        self.assertEqual(aggregate, changed_aggregate)
        self.assertNotEqual(
            header["spatial_mhd_and_face_field_sha256"],
            changed_header["spatial_mhd_and_face_field_sha256"],
        )

    def test_source_and_deck_contract_is_additive_and_pure(self) -> None:
        source = (
            ROOT / "src/pgen/tests/pic_parallel_shock.cpp"
        ).read_text(encoding="utf-8")
        boundary = (ROOT / "src/bvals/bvals_part.cpp").read_text(encoding="utf-8")
        restart_utils = (ROOT / "src/outputs/restart_utils.cpp").read_text(
            encoding="utf-8"
        )
        particles = (ROOT / "src/particles/particles_pushers.cpp").read_text(
            encoding="utf-8"
        )
        reducer_source = Path(reducer.__file__).read_text(encoding="utf-8")
        for required in (
            "ps_enable_conservation_ledger",
            "AccumulateParallelShockMHDBoundaryTransport",
            "CommitParallelShockConservationCycle",
            "ps_conservation_ledger_schema",
            "unledgered EOS floor",
            "prolong_primitives",
            "fixed uniform mesh",
            "runtime load balancing are",
            "ParallelShockExactMeshStateIsFixedUniform",
            "pmesh->max_level != pmesh->root_level",
            "loc.level != pmesh->root_level",
            "pmesh->cost_eachmb[gid] != 1.0F",
            "!pmesh->restart_meta.ncyc_since_ref.empty()",
            "occupied_root_blocks",
            "loc.lx1 < 0 || loc.lx1 >= pmesh->nmb_rootx1",
            "pmy_mesh_->two_d && ppart->pic_enable_2d3v",
            "duplicate or cycle/time-discontinuous",
            "exact conservation history is ",
        ):
            self.assertIn(required, source)
        ledger_parse = source.index(
            'ps_enable_conservation_ledger = pin->GetOrAddBoolean('
        )
        exact_z_guard = source.index(
            "pic_parallel_shock exact conservation ledger expects periodic "
        )
        self.assertGreater(exact_z_guard, ledger_parse)
        self.assertIn(
            "if (ps_enable_conservation_ledger) {",
            source[ledger_parse:exact_z_guard],
        )
        self.assertNotIn("pic_parallel_shock expects periodic z boundaries in 3D.", source)
        self.assertIn("pic_boundary_escape_ledger", boundary)
        self.assertIn('#include "outputs/restart_utils.hpp"', boundary)
        for marker in (
            "or non-physical particle destruction request.",
            "particle send targets as physical escape.",
        ):
            offset = boundary.index(marker)
            self.assertIn(
                "restart_utils::AbortOnFatalError();",
                boundary[offset:offset + 256],
            )
        self.assertIn("MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);", restart_utils)
        self.assertIn("pic_reflecting_boundary_delta", particles)
        self.assertTrue(DECK.is_file())
        for forbidden in ("subprocess", "os.system", "sbatch", "srun"):
            self.assertNotIn(forbidden, reducer_source)
        runtime_source = (
            ROOT
            / "tst/scripts/particles/"
            "pic_parallel_shock_exact_conservation_closure_vl2_tsc.py"
        ).read_text(encoding="utf-8")
        for required in (
            "_run_restart_execution_negative_guards",
            "nonunit_restored_cost",
            "restored_refined_or_adaptive_topology",
            "restored_adaptive_metadata_without_refinement",
            "exact_ledger_rejects_non_2d3v",
            "ledger_disabled_nonperiodic_3d",
            "ledger-disabled nonperiodic 3D new/base physics state differs",
            "ATHENA_PIC_EXACT_CONSERVATION_BASE_CANDIDATE_MANIFEST",
            "ATHENA_PIC_EXACT_CONSERVATION_CANDIDATE_MANIFEST",
            "ATHENA_PIC_EXACT_CONSERVATION_EXPECTED_CANDIDATE_COMMIT",
            "_LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER",
            "_LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER_SHA256",
            "linux_execveat_descriptor_bound_exact_authenticated_bytes",
            "descriptor-bound Athena execution",
            "base and candidate executable paths are identical",
            "authentication failed",
            "_run_athena_expect_fail",
        ):
            self.assertIn(required, runtime_source)
        self.assertIn(
            '"evidence_class": "bounded_source_local_engineering_self_consistency_only"',
            runtime_source,
        )
        self.assertIn('"trusted_execution_provenance": "absent"', runtime_source)
        self.assertFalse(
            (
                ROOT
                / "tst/publication/"
                "q011_section54_exact_conservation_runtime_evidence_successor_v2.py"
            ).exists()
        )
        self.assertFalse(
            (
                ROOT
                / "tst/publication/"
                "test_q011_section54_exact_conservation_runtime_evidence_successor_v2.py"
            ).exists()
        )

    def test_candidate_source_closure_binds_transitive_reducer_dependencies(
        self,
    ) -> None:
        dependencies = (
            "tst/publication/q011_section54_restart.py",
            "tst/publication/q019_nonlinear_bell_particle_state.py",
            "src/particles/particles_data_structs.hpp",
            "tst/scripts/particles/pic_parallel_shock_outer_x1_escape_restart.py",
        )
        self.assertTrue(set(dependencies).issubset(runtime._SOURCE_FILES))
        baseline = runtime._candidate_payload()
        original_read_bytes = Path.read_bytes
        for dependency in dependencies:
            target = (ROOT / dependency).resolve()

            def mutated_read_bytes(path: Path, *, bound_target: Path = target) -> bytes:
                payload = original_read_bytes(path)
                if path.resolve() == bound_target:
                    return payload + b"\n"
                return payload

            with self.subTest(dependency=dependency):
                with mock.patch.object(Path, "read_bytes", mutated_read_bytes):
                    self.assertNotEqual(runtime._candidate_payload(), baseline)

    def test_runtime_parity_case_is_ledger_disabled_nonperiodic_3d(self) -> None:
        self.assertEqual(
            runtime._LEDGER_DISABLED_NONPERIODIC_3D_OVERRIDES,
            (
                "mesh/nx3=4",
                "mesh/x3max=4.0",
                "meshblock/nx3=4",
                "mesh/ix3_bc=outflow",
                "mesh/ox3_bc=outflow",
                "particles/pic_enable_2d3v=false",
                "particles/pic_boundary_conservation_ledger=false",
                "problem/ps_enable_conservation_ledger=false",
                "problem/user_hist=false",
            ),
        )

    def test_runtime_parity_rejects_same_executable_path_and_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = root / "candidate"
            baseline = root / "baseline"
            manifest = root / "clean_candidate_manifest.json"
            candidate.write_bytes(b"candidate-bytes")
            baseline.write_bytes(b"candidate-bytes")
            manifest.write_text("{}\n", encoding="ascii")
            with runtime._bind_executable(
                candidate, label="candidate"
            ) as candidate_bound, runtime._bind_executable(
                candidate, label="same-path"
            ) as same_path:
                with self.assertRaisesRegex(RuntimeError, "paths are identical"):
                    runtime._authenticate_ledger_disabled_parity_base(
                        candidate_executable=candidate_bound,
                        base_executable=same_path,
                        base_manifest=manifest,
                    )
            with runtime._bind_executable(
                candidate, label="candidate"
            ) as candidate_bound, runtime._bind_executable(
                baseline, label="same-bytes"
            ) as same_bytes:
                with self.assertRaisesRegex(RuntimeError, "bytes are identical"):
                    runtime._authenticate_ledger_disabled_parity_base(
                        candidate_executable=candidate_bound,
                        base_executable=same_bytes,
                        base_manifest=manifest,
                    )

    def test_runtime_parity_revalidation_uses_pinned_installed_control_plane(self) -> None:
        completed = mock.Mock(returncode=0, stderr=b"", stdout=b'{"status":"passed"}')
        with mock.patch.object(runtime.subprocess, "run", return_value=completed) as run:
            self.assertEqual(
                runtime._run_clean_candidate_revalidation(
                    Path("/fixed/manifest.json"),
                    "a" * 64,
                    runtime._LEDGER_DISABLED_PARITY_BASE_COMMIT,
                ),
                {"status": "passed"},
            )
        command = run.call_args.args[0]
        self.assertEqual(
            command[0], runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_PYTHON
        )
        self.assertIn(str(runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER), command)
        self.assertIn(runtime._LEDGER_DISABLED_PARITY_BASE_COMMIT, command)
        self.assertEqual(run.call_args.kwargs["cwd"], "/")
        self.assertFalse(run.call_args.kwargs["check"])
        self.assertEqual(len(run.call_args.kwargs["pass_fds"]), 1)

    def test_runtime_parity_rejects_bogus_baseline_and_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = root / "candidate"
            candidate.write_bytes(b"candidate-bytes")
            manifest = root / "clean_candidate_manifest.json"
            manifest.write_text("{}\n", encoding="ascii")
            other = root / "other"
            other.mkdir()
            baseline = other / "athena"
            baseline.write_bytes(b"bogus-baseline")
            with runtime._bind_executable(
                candidate, label="candidate"
            ) as candidate_bound, runtime._bind_executable(
                baseline, label="base"
            ) as base_bound:
                with self.assertRaisesRegex(RuntimeError, "manifest executable"):
                    runtime._authenticate_ledger_disabled_parity_base(
                        candidate_executable=candidate_bound,
                        base_executable=base_bound,
                        base_manifest=manifest,
                    )

            baseline = root / "athena"
            baseline.write_bytes(b"bogus-baseline")
            forged_report = {
                "status": "passed",
                "source": {
                    "git_commit": runtime._LEDGER_DISABLED_PARITY_BASE_COMMIT,
                    "git_tree": "0" * 40,
                },
                "build": {
                    "executable_sha256": hashlib.sha256(
                        baseline.read_bytes()
                    ).hexdigest(),
                    "profile_id": runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
                },
                "clean_candidate_manifest": {
                    "sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                    "path": str(manifest),
                },
                "current_control_plane_version":
                    runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION,
            }
            with mock.patch.object(
                runtime,
                "_run_clean_candidate_revalidation",
                return_value=forged_report,
            ), runtime._bind_executable(
                candidate, label="candidate"
            ) as candidate_bound, runtime._bind_executable(
                baseline, label="base"
            ) as base_bound:
                with self.assertRaisesRegex(RuntimeError, "does not match the pinned"):
                    runtime._authenticate_ledger_disabled_parity_base(
                        candidate_executable=candidate_bound,
                        base_executable=base_bound,
                        base_manifest=manifest,
                    )

    def test_runtime_parity_exact_base_production_profile_is_constructible(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = root / "candidate"
            baseline = root / "athena"
            manifest = root / "clean_candidate_manifest.json"
            candidate.write_bytes(b"candidate-bytes")
            baseline.write_bytes(b"base-production-profile-bytes")
            manifest.write_text("{}\n", encoding="ascii")
            with runtime._bind_executable(
                candidate, label="candidate"
            ) as candidate_bound, runtime._bind_executable(
                baseline, label="base"
            ) as base_bound:
                report = {
                    "status": "passed",
                    "source": {
                        "git_commit": runtime._LEDGER_DISABLED_PARITY_BASE_COMMIT,
                        "git_tree": runtime._LEDGER_DISABLED_PARITY_BASE_TREE,
                    },
                    "build": {
                        "executable_sha256": base_bound.sha256,
                        "profile_id": runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
                        "receipt_control_plane_version": "a" * 64,
                    },
                    "clean_candidate_manifest": {
                        "sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                        "path": str(manifest),
                    },
                    "current_control_plane_version":
                        runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION,
                }
                with mock.patch.object(
                    runtime, "_run_clean_candidate_revalidation", return_value=report
                ):
                    result = runtime._authenticate_ledger_disabled_parity_base(
                        candidate_executable=candidate_bound,
                        base_executable=base_bound,
                        base_manifest=manifest,
                    )
            self.assertEqual(result["executable_sha256"], report["build"]["executable_sha256"])
            self.assertEqual(
                result["build_profile_id"],
                runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
            )
            self.assertEqual(
                result["execution_binding"],
                "linux_execveat_descriptor_bound_exact_authenticated_bytes",
            )

    def test_runtime_parity_candidate_requires_authenticated_production_profile(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_root = root / "candidate-clean"
            candidate_root.mkdir()
            candidate = candidate_root / "athena"
            manifest = candidate_root / "clean_candidate_manifest.json"
            candidate.write_bytes(b"candidate-production-profile-bytes")
            manifest.write_text("{}\n", encoding="ascii")
            expected_commit = "2" * 40
            with runtime._bind_executable(
                candidate, label="candidate"
            ) as candidate_bound:
                report = {
                    "status": "passed",
                    "source": {
                        "git_commit": expected_commit,
                        "git_tree": "3" * 40,
                    },
                    "build": {
                        "executable_sha256": candidate_bound.sha256,
                        "profile_id": runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
                        "receipt_control_plane_version": "a" * 64,
                    },
                    "clean_candidate_manifest": {
                        "sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                        "path": str(manifest),
                    },
                    "current_control_plane_version":
                        runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION,
                }
                with mock.patch.object(
                    runtime, "_run_clean_candidate_revalidation", return_value=report
                ):
                    result = runtime._authenticate_ledger_disabled_parity_candidate(
                        candidate_executable=candidate_bound,
                        candidate_manifest=manifest,
                        expected_commit=expected_commit,
                    )
            self.assertEqual(result["source_commit"], expected_commit)
            self.assertEqual(result["source_tree"], "3" * 40)
            self.assertEqual(
                result["build_profile_id"],
                runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
            )

    def test_bound_athena_rejects_transient_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            work = root / "work"
            source.mkdir()
            work.mkdir()
            executable = source / "athena"
            original = b"#!/bin/sh\nprintf 'bound-original\\n'\n"
            executable.write_bytes(original)
            executable.chmod(0o555)
            deck = work / "input.athinput"
            deck.write_text("# fixture\n", encoding="ascii")
            with runtime._bind_executable(
                executable, label="athena-replacement-probe"
            ) as bound:
                replacement = source / "replacement"
                replacement.write_bytes(b"#!/bin/sh\nprintf 'replacement\\n'\n")
                os.replace(replacement, executable)
                restored = source / "restored"
                restored.write_bytes(original)
                os.replace(restored, executable)
                with self.assertRaisesRegex(RuntimeError, "changed while retained"):
                    runtime._execute_athena(
                        label="bound_replacement_probe",
                        exe_dir=work,
                        basename="bound_probe",
                        deck=deck,
                        bound_executable=bound,
                    )

    def test_bound_athena_rejects_transient_in_place_rewrite_and_restore(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "athena"
            original = b"#!/bin/sh\nprintf 'bound-original\\n'\n"
            executable.write_bytes(original)
            executable.chmod(0o755)
            with runtime._bind_executable(
                executable, label="athena-in-place-rewrite-probe"
            ) as bound:
                executable.write_bytes(b"#!/bin/sh\nprintf 'replacement\\n'\n")
                executable.write_bytes(original)
                with self.assertRaisesRegex(RuntimeError, "changed while retained"):
                    runtime._require_bound_executable_unchanged(
                        bound, label="in-place-rewrite-probe"
                    )

    def test_bound_athena_runs_with_read_only_source_and_writable_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            work = root / "work"
            source.mkdir()
            work.mkdir()
            executable = source / "athena"
            executable.write_bytes(b"#!/bin/sh\ntouch output.marker\n")
            executable.chmod(0o555)
            deck = work / "input.athinput"
            deck.write_text("# fixture\n", encoding="ascii")
            with runtime._bind_executable(
                executable, label="read-only-source-probe"
            ) as bound:
                source.chmod(0o555)
                try:
                    record = runtime._execute_athena(
                        label="read_only_source_probe",
                        exe_dir=work,
                        basename="sealed_probe",
                        deck=deck,
                        bound_executable=bound,
                    )
                finally:
                    source.chmod(0o755)
            self.assertEqual(record["returncode"], 0)
            self.assertTrue((work / "output.marker").is_file())
            self.assertFalse((source / "output.marker").exists())

    def test_descriptor_bound_execution_drops_hostile_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "athena"
            executable.write_bytes(b"#!/bin/sh\nexit 0\n")
            executable.chmod(0o555)
            completed = runtime.subprocess.CompletedProcess([], 0, "", "")
            with runtime._bind_executable(
                executable, label="sanitized-environment-probe"
            ) as bound, mock.patch.dict(
                os.environ,
                {"LD_PRELOAD": "/tmp/hostile.so", "PYTHONPATH": "/tmp/hostile"},
            ), mock.patch.object(
                runtime.subprocess, "run", return_value=completed
            ) as run:
                runtime._run_descriptor_bound_executable(
                    executable=bound,
                    arguments=[str(executable)],
                    cwd=root,
                )
            environment = run.call_args.kwargs["env"]
            self.assertEqual(environment, runtime._SANITIZED_EXECUTION_ENV)
            self.assertNotIn("LD_PRELOAD", environment)
            self.assertNotIn("PYTHONPATH", environment)

    def test_sealed_python_runner_ignores_transient_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            script = root / "runner.py"
            script.write_text("print('sealed-runner-original')\n", encoding="ascii")
            with runtime._seal_python_script(
                script, label="runner-replacement-probe"
            ) as sealed:
                script.write_text("print('replacement')\n", encoding="ascii")
                completed = runtime._run_sealed_python_script(
                    sealed_script=sealed,
                    original_path=script,
                    arguments=[],
                    cwd=str(root),
                    env={"PATH": "/usr/bin:/bin"},
                    timeout=30,
                )
            self.assertEqual(completed.returncode, 0)
            self.assertEqual(completed.stdout, b"sealed-runner-original\n")
            self.assertEqual(completed.stderr, b"")

    def test_paper_vl2_boundary_stage_and_checkpoint_order_are_explicit(self) -> None:
        driver = (ROOT / "src/driver/driver.cpp").read_text(encoding="utf-8")
        tasks = (ROOT / "src/mhd/mhd_tasks.cpp").read_text(encoding="utf-8")
        update = (ROOT / "src/mhd/mhd_update.cpp").read_text(encoding="utf-8")
        source = (ROOT / "src/pgen/tests/pic_parallel_shock.cpp").read_text(
            encoding="utf-8"
        )
        runtime = (
            ROOT
            / "tst/scripts/particles/"
            "pic_parallel_shock_exact_conservation_closure_vl2_tsc.py"
        ).read_text(encoding="utf-8")

        for snippet in (
            "gam0[0] = 0.0;",
            "gam1[0] = 1.0;",
            "beta[0] = 0.5;",
            "gam0[1] = 0.0;",
            "gam1[1] = 1.0;",
            "beta[1] = 1.0;",
        ):
            self.assertIn(snippet, driver)
        self.assertLess(driver.index("user_work_in_loop_func)(pmesh)"),
                        driver.index("pmesh->time = pmesh->time + pmesh->dt"))
        nonrestart_output = driver.index("if (!IsRestartOutput(out)")
        amr = driver.index("AdaptiveMeshRefinement(this, pin)")
        restart_output = driver.index("if (IsRestartOutput(out)", amr)
        self.assertLess(nonrestart_output, amr)
        self.assertLess(amr, restart_output)
        self.assertLess(tasks.index("id.rkupdt"), tasks.index("id.srctrms"))
        self.assertIn(
            "u0_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) "
            "- beta_dt*divf(i);",
            update,
        )
        self.assertIn(
            "ps_conservation_mhd_boundary_cycle_local[n] =\n"
            "          stage_weight*boundary_delta.the_array[n];",
            source,
        )
        self.assertIn(
            "nonzero_boundary_transport_closes_at_roundoff_with_cycle_local_stage2_overwrite",
            runtime,
        )

    def test_readiness_record_is_non_authorizing_and_hash_bound(self) -> None:
        record = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(record["qualification_effect"], "none")
        self.assertIn("trusted_execution_pending", record["status"])
        self.assertEqual(
            record["isolation"]["base_commit"],
            "1d72534619da277fdd6838cabdab0d465d0a8867",
        )
        self.assertEqual(
            record["source_fail_closed_contract"]["qualified_topology"],
            "fixed uniform mesh only",
        )
        self.assertEqual(
            record["source_fail_closed_contract"]["qualified_boundary_geometry"],
            "reflecting inner-x1, inflow or outflow outer-x1, periodic transverse boundaries",
        )
        self.assertEqual(
            record["evidence_boundary"]["source_local_harness"],
            "engineering_self_consistency_only",
        )
        self.assertFalse(record["evidence_boundary"]["trusted_execution_provenance"])
        self.assertFalse(record["evidence_boundary"]["registered_science_evidence"])
        self.assertEqual(
            record["evidence_boundary"]["self_attested_runtime_packet_verifier"],
            "removed_not_authoritative",
        )
        parity = record["source_local_engineering_runtime"]["ledger_disabled_parity"]
        self.assertEqual(
            parity["status"], "pending_authenticated_distinct_baseline_rerun"
        )
        self.assertEqual(
            parity["prior_unauthenticated_parity_result"],
            "invalidated_not_retained",
        )
        self.assertEqual(
            parity["required_base_commit"],
            runtime._LEDGER_DISABLED_PARITY_BASE_COMMIT,
        )
        self.assertEqual(
            parity["required_base_tree"],
            runtime._LEDGER_DISABLED_PARITY_BASE_TREE,
        )
        self.assertEqual(
            parity["required_base_build_profile"],
            runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
        )
        self.assertEqual(
            parity["required_candidate_build_profile"],
            runtime._LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE,
        )
        self.assertEqual(
            parity["required_base_executable_sha256"],
            "derived_from_authenticated_clean_candidate_and_bound_to_execveat_execution",
        )
        self.assertEqual(
            parity["required_control_plane_version"],
            runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION,
        )
        self.assertEqual(
            parity["required_control_plane_runner_sha256"],
            runtime._LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER_SHA256,
        )
        self.assertEqual(
            parity["execution_binding"],
            "linux_execveat_descriptor_bound_exact_authenticated_bytes_with_mutation_watch_sanitized_environment_separate_writable_outputs_and_full_spatial_restart_state_comparison",
        )
        self.assertIn(
            "trusted installed-control-plane receipts",
            record["remaining_limitations"],
        )
        self.assertFalse(record["isolation"]["protected_historical_files_edited"])
        self.assertTrue(all(value is False for value in record["authority"].values()))
        artifact_paths = {artifact["path"] for artifact in record["artifacts"]}
        self.assertTrue(
            {
                "tst/publication/q011_section54_restart.py",
                "tst/publication/q019_nonlinear_bell_particle_state.py",
            }.issubset(artifact_paths)
        )
        for artifact in record["artifacts"]:
            path = ROOT / artifact["path"]
            self.assertTrue(path.is_file(), artifact["path"])
            self.assertEqual(_sha(path.read_bytes()), artifact["sha256"], artifact["path"])


if __name__ == "__main__":
    unittest.main()
