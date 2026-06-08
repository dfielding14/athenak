from __future__ import annotations

import json
from copy import deepcopy

import pytest

from tst.publication import q011_escape_event_evidence as evidence


def _header(rank: int = 0, *, basename: str = "segment", cycle: int = 0):
    time = float(cycle)
    return {
        "record_type": "q011_escape_event_stream_header",
        "schema_version": 1,
        "control_fingerprint": "v1:0123456789abcdef",
        "basename": basename,
        "rank": rank,
        "nranks": 1,
        "segment_start_cycle": cycle,
        "segment_start_time": time,
        "state_kind": "momentum_per_mass",
        "light_speed": 10000.0,
        "injected_species": 0,
        "species_mass": 1.0,
        "species_charge": 1.0,
        "q_over_m": 1.0,
        "macro_mass": 0.001,
        "field_interpolation": "tsc_bcc0_w0",
        "mesh_x1min": 0.0,
        "mesh_x1max": 4.0,
        "mesh_x2min": 0.0,
        "mesh_x2max": 1.0,
        "mesh_x3min": 0.0,
        "mesh_x3max": 0.0,
        "escape_audit_calls_at_start": 2 * cycle,
        "escape_last_audit_time_at_start": 0.0 if cycle == 0 else time,
        "global_escape_count_at_start": 0.0,
        "global_escape_mass_at_start": 0.0,
        "global_escape_momentum_x1_at_start": 0.0,
        "global_escape_momentum_x2_at_start": 0.0,
        "global_escape_momentum_x3_at_start": 0.0,
        "global_escape_energy_at_start": 0.0,
        "global_initial_escape_count_at_start": 0.0,
        "global_escape_term_count_at_start": 0.0,
        "global_escape_abs_mass_at_start": 0.0,
        "global_escape_abs_momentum_x1_at_start": 0.0,
        "global_escape_abs_momentum_x2_at_start": 0.0,
        "global_escape_abs_momentum_x3_at_start": 0.0,
        "global_escape_abs_energy_at_start": 0.0,
    }


def _event(tag: int = 7, *, cycle: int = 0):
    state = (3.0, 4.0, 0.0)
    state_squared = sum(value * value for value in state)
    energy = state_squared / ((1.0 + state_squared / 10000.0**2) ** 0.5 + 1.0)
    return {
        "record_type": "q011_escape_event",
        "schema_version": 1,
        "rank": 0,
        "rank_event_index": 0,
        "cycle": cycle,
        "stage": 2,
        "audit_time": float(cycle) + 1.0,
        "tag": tag,
        "source": 1,
        "species": 0,
        "destruction_reason": 1,
        "physical_boundary_mask": 2,
        "parent_gid": 3,
        "x1": 4.01,
        "x2": 0.5,
        "x3": 0.0,
        "state_x1": state[0],
        "state_x2": state[1],
        "state_x3": state[2],
        "state_kind": "momentum_per_mass",
        "light_speed": 10000.0,
        "species_mass": 1.0,
        "species_charge": 1.0,
        "q_over_m": 1.0,
        "macro_weight": 1.0,
        "macro_mass": 0.001,
        "kinetic_energy_per_mass": energy,
        "b1": 1.0,
        "b2": 0.1,
        "b3": 0.2,
        "fluid_v1": -10.0,
        "fluid_v2": 0.0,
        "fluid_v3": 0.0,
        "field_interpolation": "tsc_bcc0_w0",
        "frame_velocity_x1": 0.0,
        "shock_speed": 3.0,
    }


def _stream(header=None, events=None):
    header = _header() if header is None else header
    events = [_event()] if events is None else events
    prefix_lines = [
        json.dumps(header, separators=(",", ":")),
        *(json.dumps(event, separators=(",", ":")) for event in events),
    ]
    prefix = ("".join(line + "\n" for line in prefix_lines)).encode("ascii")
    trailer = {
        "record_type": "q011_escape_event_stream_trailer",
        "schema_version": 1,
        "rank": header["rank"],
        "rank_event_count": len(events),
        "prefix_fnv1a64": f"{evidence._fnv1a64(prefix):016x}",
    }
    return prefix + json.dumps(trailer, separators=(",", ":")).encode("ascii") + b"\n"


def test_parse_and_recompute_raw_event():
    stream = evidence.parse_stream_bytes(_stream())
    reduction = evidence.reduce_streams([stream])
    final = reduction["derived_final_ledger"]
    assert final["count"] == 1.0
    assert final["mass"] == pytest.approx(0.001)
    assert final["momentum_x1"] == pytest.approx(0.003)
    assert final["momentum_x2"] == pytest.approx(0.004)
    assert final["energy"] == pytest.approx(0.001 * _event()["kinetic_energy_per_mass"])


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda header, event: event.update({"source": 0}), "source/species"),
        (lambda header, event: event.update({"x1": 3.99}), "does not cross"),
        (
            lambda header, event: event.update({"kinetic_energy_per_mass": 0.0}),
            "kinetic energy",
        ),
        (lambda header, event: event.update({"macro_weight": 2.0}), "macro_weight"),
        (lambda header, event: event.update({"unknown": 1}), "schema mismatch"),
    ],
)
def test_event_corruption_rejected(mutation, message):
    header = _header()
    event = _event()
    mutation(header, event)
    with pytest.raises(evidence.EscapeEvidenceError, match=message):
        evidence.parse_stream_bytes(_stream(header, [event]))


def test_prefix_tampering_rejected():
    payload = _stream().replace(b'"b1":1.0', b'"b1":2.0')
    with pytest.raises(evidence.EscapeEvidenceError, match="FNV-1a"):
        evidence.parse_stream_bytes(payload)


def test_duplicate_json_field_rejected():
    payload = _stream().replace(
        b'{"record_type":"q011_escape_event_stream_header"',
        b'{"rank":0,"record_type":"q011_escape_event_stream_header"',
        1,
    )
    with pytest.raises(evidence.EscapeEvidenceError, match="duplicate JSON field"):
        evidence.parse_stream_bytes(payload)


def test_duplicate_particle_tag_across_segments_rejected():
    first = evidence.parse_stream_bytes(_stream())
    next_header = _header(basename="continuation", cycle=1)
    first_final = evidence.reduce_streams([first])["derived_final_ledger"]
    for header_name, ledger_name in evidence._START_TO_LEDGER.items():
        next_header[header_name] = first_final[ledger_name]
    second_event = _event(tag=7, cycle=1)
    second = evidence.parse_stream_bytes(_stream(next_header, [second_event]))
    with pytest.raises(evidence.EscapeEvidenceError, match="duplicate escaped"):
        evidence.reduce_streams([first, second])


def test_continuation_must_join_prior_derived_ledger():
    first = evidence.parse_stream_bytes(_stream())
    next_header = _header(basename="continuation", cycle=1)
    second = evidence.parse_stream_bytes(_stream(next_header, [_event(tag=8, cycle=1)]))
    with pytest.raises(evidence.EscapeEvidenceError, match="does not join"):
        evidence.reduce_streams([first, second])


def test_restart_ledger_comparison_is_fail_closed():
    reduction = evidence.reduce_streams([evidence.parse_stream_bytes(_stream())])
    ledger = {
        "ps_escaped_injected_cr_count_global": 1.0,
        "ps_escaped_injected_cr_mass_global": 0.001,
        "ps_escaped_injected_cr_momentum_x1_global": 0.003,
        "ps_escaped_injected_cr_momentum_x2_global": 0.004,
        "ps_escaped_injected_cr_momentum_x3_global": 0.0,
        "ps_escaped_injected_cr_energy_global": reduction["derived_final_ledger"][
            "energy"
        ],
        "ps_escaped_initial_cr_count_global": 0.0,
        "ps_escaped_injected_cr_term_count_global": 1.0,
        "ps_escaped_injected_cr_abs_mass_global": 0.001,
        "ps_escaped_injected_cr_abs_momentum_x1_global": 0.003,
        "ps_escaped_injected_cr_abs_momentum_x2_global": 0.004,
        "ps_escaped_injected_cr_abs_momentum_x3_global": 0.0,
        "ps_escaped_injected_cr_abs_energy_global": reduction["derived_final_ledger"][
            "energy"
        ],
    }
    evidence.assert_matches_restart_ledger(reduction, ledger)
    corrupted = deepcopy(ledger)
    corrupted["ps_escaped_injected_cr_energy_global"] = 0.0
    with pytest.raises(evidence.EscapeEvidenceError, match="does not match"):
        evidence.assert_matches_restart_ledger(reduction, corrupted)
