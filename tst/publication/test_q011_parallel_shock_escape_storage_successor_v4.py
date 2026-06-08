#!/usr/bin/env python3
"""Tests for the Q011 compact escape-evidence storage successor."""

from __future__ import annotations

import pytest

from tst.publication import q011_parallel_shock_escape_storage_successor_v4 as storage


def _quantity(value):
    if isinstance(value, int):
        return value
    return value["numerator"] / value["denominator"]


def test_default_campaign_adds_all_injected_particles_escape_bound():
    report = storage.estimate_storage_with_escape_evidence()
    escape = report["escape_evidence"]
    assert escape["event_bytes"] == 176
    assert escape["maximum_injected_particle_count_per_run"] == 166_400_001
    assert escape["maximum_event_payload_bytes_per_run"] == 29_286_400_176
    assert escape["maximum_framing_bytes_per_run"] == 4_194_304
    assert escape["maximum_escape_evidence_bytes_per_run"] == 29_290_594_480
    assert escape["maximum_escape_evidence_bytes_campaign"] == 702_974_267_520
    assert report["planning_envelope"]["campaign"][
        "reservation_envelope_tb_decimal"
    ] > 11.5


def test_successor_recomputes_envelope_from_predecessor_plus_escape_allowance():
    report = storage.estimate_storage_with_escape_evidence(
        grid_count=2,
        seed_count=3,
        maximum_ranks_per_run=8,
        maximum_segments_per_run=2,
    )
    campaign = report["planning_envelope"]["campaign"]
    predecessor = _quantity(
        campaign["predecessor_logical_bytes_before_escape_evidence"]
    )
    production_mesh = campaign["production_mesh_bin_increment_bytes"]
    escape = campaign["escape_evidence_bytes_allowance"]
    logical = _quantity(campaign["logical_bytes_before_filesystem_overhead"])
    assert logical == pytest.approx(predecessor + production_mesh + escape)
    assert escape == 6 * report["escape_evidence"][
        "maximum_escape_evidence_bytes_per_run"
    ]
    for variant in report["planning_envelope"]["variants"]:
        assert variant["escape_evidence_bytes_per_run_allowance"] == report[
            "escape_evidence"
        ]["maximum_escape_evidence_bytes_per_run"]


def test_rank_and_segment_bounds_are_explicit_positive_parameters():
    report = storage.estimate_storage_with_escape_evidence(
        maximum_ranks_per_run=7, maximum_segments_per_run=3
    )
    assert report["escape_evidence"]["maximum_framing_bytes_per_run"] == (
        7 * 3 * storage.ESCAPE_STREAM_FRAMING_BYTES_PER_RANK_SEGMENT_ALLOWANCE
    )
    with pytest.raises(storage.EscapeStorageError, match="positive integer"):
        storage.estimate_storage_with_escape_evidence(maximum_ranks_per_run=0)
