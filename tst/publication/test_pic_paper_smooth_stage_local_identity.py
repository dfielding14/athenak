#!/usr/bin/env python3
"""Focused source contract for paper_smooth stage-local receiver identity."""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def _source(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="ascii")


def test_record_layout_and_transport_key_bind_stage_local_source_identity() -> None:
    record = _source("src/particles/particles_data_structs.hpp")
    header = _source("src/bvals/bvals.hpp")
    transport = _source("src/bvals/bvals_mom.cpp")

    assert "std::int32_t source_rank;" in record
    assert "std::int32_t source_index;" in record
    assert "std::uint32_t source_rank;" in header
    assert "std::uint32_t source_index;" in header
    assert "sizeof(Record) == 6*sizeof(std::uint32_t) + 12*sizeof(Real)" in transport

    make_key = transport[transport.index("PaperSmoothMomentRecordTransport::MakeKey") :]
    make_key = make_key[: make_key.index("//--------------------------------")]
    assert "record.source_rank" in make_key
    assert "record.source_index" in make_key
    assert "record.ptag" not in make_key


def test_deposition_populates_stage_local_identity() -> None:
    deposition = _source("src/particles/particles_moments.cpp")

    assert "owner_gid, ptag, global_variable::my_rank, p, deposit_flags, 0U" in deposition


def test_runtime_regression_covers_duplicate_tags_and_remote_mpi_receiver() -> None:
    pgen = _source("src/pgen/tests/pic_paper_smooth_tsc_interface.cpp")
    regression = _source("tst/scripts/particles/pic_paper_smooth_tsc_interface.py")
    decks = (
        _source("inputs/tests/pic_paper_smooth_tsc_interface.athinput"),
        _source("inputs/tests/pic_paper_smooth_tsc_interface_3d.athinput"),
    )

    assert '"problem", "duplicate_ptag_pair", false' in pgen
    assert "h_pi(PTAG, p) = 0;" in pgen
    assert "duplicate_ptag_pair ? static_cast<Real>(p + 1) : 1.0" in pgen
    assert all("duplicate_ptag_pair = false" in deck for deck in decks)
    assert "_DUPLICATE_PTAG_WEIGHT_SUM * reference['actual_cells']" in regression
    assert "'mpi2_g:duplicate_ptag_pair_remote_receiver_block'" in regression
