#!/usr/bin/env python3
"""Tests for the Q019 Q023-carrier nonlinear Bell redesign."""

from __future__ import annotations

import copy
import math
from pathlib import Path
import re

import pytest

from tst.publication import q019_q023_carrier_nonlinear_bell_redesign_v1 as redesign
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as q019


SOURCE = (
    Path(__file__).resolve().parents[2]
    / "src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.cpp"
)


def test_redesign_has_exact_staged_inventory_and_no_authority() -> None:
    cases = redesign.expected_cases()
    redesign.validate_cases(cases)
    assert len(cases) == 14
    assert {case["pilot_stage"] for case in cases} == {1, 2, 3}
    assert all(case["saturation_candidate"] is False for case in cases)
    manifest, _ = redesign.build_deck_manifest()
    assert all(value is False for value in manifest["authorization"].values())
    assert manifest["physics"]["external_current_surrogate_not_physical_cr_density"]
    assert manifest["physics"]["exact_locked_current"] is False


def test_fiducial_matches_q023_carrier_and_volume_aware_current() -> None:
    case = next(
        item
        for item in redesign.expected_cases()
        if item["case_id"] == "q019-q023-carrier-s2-window-fiducial-s0"
    )
    assert math.isclose(case["species_charge"], 2.0 * math.pi * 1.0e-6)
    assert math.isclose(case["guide_parallel_stream_speed"], 2.5)
    assert math.isclose(case["rho_cr_over_rho0"], 800000.0)
    assert math.isclose(case["k0_rg0"], 2.5e6)
    assert math.isclose(
        q019.configured_j_over_c(case), 4.0 * math.pi, rel_tol=1.0e-13
    )
    assert case["species_q_over_mc_matches_background"] is False
    assert case["cr_inertia_parameter"] == case["rho_cr_over_rho0"]


def test_real_unmocked_cycle_estimator_is_structurally_feasible() -> None:
    summary = redesign.feasibility_summary()
    assert summary["maximum_case_cycles"] <= redesign.MAXIMUM_STRUCTURAL_CYCLE_COUNT
    assert summary["maximum_root_cells"] <= redesign.MAXIMUM_PILOT_ROOT_CELLS
    assert summary["maximum_seconds_per_cycle_to_fit_cap"] >= 1.0
    assert summary["nonlinear_B_over_B0_safety_envelope"] == 10.0
    assert summary["measured_resource_authority"] is False
    assert summary["execution_authorized"] is False


def test_compiled_runtime_registry_binds_every_redesign_deck() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    registry = dict(
        re.findall(
            r'\{"(q019-q023-carrier-[^"]+)",\s*"([0-9a-f]{64})"\}',
            source[
                source.index("q019_canonical_matrix_identities[]") :
                source.index("Q019CanonicalMatrixFingerprint")
            ],
        )
    )
    expected = {
        str(case["case_id"]): str(case["matrix_identity_fingerprint"])
        for case in redesign.expected_cases()
    }
    assert registry == expected


def test_superseded_fast_particle_design_is_not_structurally_feasible() -> None:
    old = next(
        case
        for case in q019.expected_cases()
        if case["case_id"] == "q019-hr-current-retention-s0"
    )
    minimum_dx = min(
        old["extents"][index] / old["nx"][index]
        for index in range(old["dimension"])
    )
    timestep = old["cfl"] * old["pic_max_cell_cross"] * minimum_dx / old[
        "guide_parallel_stream_speed"
    ]
    cycles = math.ceil(old["terminal_time"] / timestep)
    assert cycles > 100 * redesign.MAXIMUM_STRUCTURAL_CYCLE_COUNT


def test_tampered_current_or_runtime_identity_fails_closed() -> None:
    case = copy.deepcopy(redesign.expected_cases()[0])
    case["deposit_qscale"] *= 1.01
    case["matrix_identity_fingerprint"] = q019.matrix_identity_fingerprint(case)
    with pytest.raises(redesign.RedesignError, match="current closure"):
        redesign.validate_cases((case,) + redesign.expected_cases()[1:])

    valid = redesign.expected_cases()[0]
    text = q019.render_deck(valid).replace(
        "source_lineage = "
        + redesign.SOURCE_LINEAGE,
        "source_lineage = forged",
        1,
    )
    with pytest.raises(q019.ContractError, match="immutable runtime deck semantics"):
        q019.validate_rendered_deck(valid, text)
