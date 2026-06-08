#!/usr/bin/env python3
"""Adversarial tests for the Q011 retained escape-evidence successor."""

from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
from unittest import mock

import pytest

from tst.publication import q011_escape_event_evidence as evidence
from tst.publication import (
    q011_section54_escape_evidence_admission_successor_v2 as successor,
)
from tst.publication import test_q011_escape_event_evidence as stream_fixture


ATTEMPT_ID = "q011-production-coarse-seed23050101"


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _seal(root: Path) -> None:
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def _restart_payload(ledger: dict[str, float], *, cycle: int = 1, time: float = 1.0) -> bytes:
    parameters = {
        "ps_escape_ledger_schema": "2",
        "ps_escape_ledger_complete": "true",
        "ps_escape_audit_calls": str(2 * cycle),
        "ps_escape_last_audit_time": repr(time),
        "ps_escaped_injected_cr_count_global": repr(ledger["count"]),
        "ps_escaped_injected_cr_mass_global": repr(ledger["mass"]),
        "ps_escaped_injected_cr_momentum_x1_global": repr(ledger["momentum_x1"]),
        "ps_escaped_injected_cr_momentum_x2_global": repr(ledger["momentum_x2"]),
        "ps_escaped_injected_cr_momentum_x3_global": repr(ledger["momentum_x3"]),
        "ps_escaped_injected_cr_energy_global": repr(ledger["energy"]),
        "ps_escaped_initial_cr_count_global": repr(ledger["initial_count"]),
        "ps_escaped_injected_cr_term_count_global": repr(ledger["term_count"]),
        "ps_escaped_injected_cr_abs_mass_global": repr(ledger["abs_mass"]),
        "ps_escaped_injected_cr_abs_momentum_x1_global": repr(
            ledger["abs_momentum_x1"]
        ),
        "ps_escaped_injected_cr_abs_momentum_x2_global": repr(
            ledger["abs_momentum_x2"]
        ),
        "ps_escaped_injected_cr_abs_momentum_x3_global": repr(
            ledger["abs_momentum_x3"]
        ),
        "ps_escaped_injected_cr_abs_energy_global": repr(ledger["abs_energy"]),
    }
    text = [f"<time>\ncycle={cycle}", "<problem>"]
    text.extend(f"{key}={value}" for key, value in parameters.items())
    return ("\n".join(text) + "\n<par_end>\n").encode("ascii") + b"binary"


def _fixture(root: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    raw = root / "raw"
    raw.mkdir()
    header = stream_fixture._header(basename=ATTEMPT_ID)
    stream_payload = stream_fixture._binary_stream(
        header, [stream_fixture._event()]
    )
    stream = evidence.parse_stream_bytes(stream_payload)
    reduction = evidence.reduce_streams([stream])
    restart_payload = _restart_payload(reduction["derived_final_ledger"])
    stream_path = (
        f"{ATTEMPT_ID}.q011_escape_events.cycle00000000.rank00000000.bin"
    )
    restart_path = "rst/rank_00000000/terminal.rst"
    (raw / stream_path).write_bytes(stream_payload)
    (raw / restart_path).parent.mkdir(parents=True)
    (raw / restart_path).write_bytes(restart_payload)
    inventory = [
        {
            "path": path,
            "sha256": _sha(payload),
            "byte_count": len(payload),
        }
        for path, payload in (
            (restart_path, restart_payload),
            (stream_path, stream_payload),
        )
    ]
    admission = {
        "campaign_id": "Q011-SECTION54-PRODUCTION-SCIENCE-SUCCESSOR-V1",
        "claim_id": "CLAIM-PAPER-SHOCK-001",
        "attempt_identity": {
            "attempt_id": ATTEMPT_ID,
            "variant": "coarse_uniform_dx12",
            "seed": 23050101,
        },
        "execution_binding": {"raw_output_root": str(raw)},
        "raw_artifacts": [
            {
                "path": restart_path,
                "sha256": _sha(restart_payload),
                "byte_count": len(restart_payload),
            }
        ],
        "snapshot_bindings": [
            {
                "cycle": 1,
                "observed_committed_time": 1.0,
                "artifact_bindings": {
                    "rst": [
                        {
                            "path": restart_path,
                            "sha256": _sha(restart_payload),
                            "byte_count": len(restart_payload),
                        }
                    ]
                },
            }
        ],
    }
    _seal(raw)
    return admission, inventory


def _build(admission, inventory):
    with mock.patch.object(
        successor.production_v1,
        "validate_attempt_admission",
        return_value=admission,
    ):
        return successor.build_escape_evidence_admission(
            attempt_admission=admission,
            retained_inventory=inventory,
            source_root="/unused",
        )


def test_successor_binds_sealed_inventory_raw_reduction_and_terminal_restart():
    with tempfile.TemporaryDirectory() as temporary:
        admission, inventory = _fixture(Path(temporary))
        record = _build(admission, inventory)
        assert record["retained_file_count"] == 2
        assert record["independent_reduction"]["unique_particle_tags"] == 1
        assert record["terminal_restart"]["cycle"] == 1
        assert all(value is False for value in record["authorization"].values())


def test_mutable_or_incomplete_retained_tree_is_rejected():
    with tempfile.TemporaryDirectory() as temporary:
        admission, inventory = _fixture(Path(temporary))
        raw = Path(admission["execution_binding"]["raw_output_root"])
        raw.chmod(0o755)
        with pytest.raises(successor.EscapeEvidenceAdmissionError, match="not sealed"):
            _build(admission, inventory)

    with tempfile.TemporaryDirectory() as temporary:
        admission, inventory = _fixture(Path(temporary))
        with pytest.raises(
            successor.EscapeEvidenceAdmissionError, match="does not exactly cover"
        ):
            _build(admission, inventory[:-1])


def test_terminal_restart_ledger_mismatch_is_rejected():
    with tempfile.TemporaryDirectory() as temporary:
        admission, inventory = _fixture(Path(temporary))
        raw = Path(admission["execution_binding"]["raw_output_root"])
        restart = raw / str(admission["raw_artifacts"][0]["path"])
        raw.chmod(0o755)
        restart.parent.chmod(0o755)
        restart.chmod(0o644)
        payload = restart.read_bytes().replace(
            b"ps_escaped_injected_cr_count_global=1.0",
            b"ps_escaped_injected_cr_count_global=2.0",
        )
        restart.write_bytes(payload)
        admission["raw_artifacts"][0].update(
            sha256=_sha(payload), byte_count=len(payload)
        )
        admission["snapshot_bindings"][0]["artifact_bindings"]["rst"][0].update(
            sha256=_sha(payload), byte_count=len(payload)
        )
        inventory[0].update(sha256=_sha(payload), byte_count=len(payload))
        _seal(raw)
        with pytest.raises(evidence.EscapeEvidenceError, match="does not match"):
            _build(admission, inventory)


def test_nonzero_first_segment_and_nonbinary_stream_are_rejected():
    with tempfile.TemporaryDirectory() as temporary:
        admission, inventory = _fixture(Path(temporary))
        raw = Path(admission["execution_binding"]["raw_output_root"])
        stream_path = raw / inventory[1]["path"]
        raw.chmod(0o755)
        stream_path.chmod(0o644)
        header = stream_fixture._header(basename=ATTEMPT_ID, cycle=1)
        payload = stream_fixture._binary_stream(
            header, [stream_fixture._event(cycle=1)]
        )
        stream_path.write_bytes(payload)
        inventory[1].update(sha256=_sha(payload), byte_count=len(payload))
        _seal(raw)
        with pytest.raises(evidence.EscapeEvidenceError, match="cycle/time zero"):
            _build(admission, inventory)

    with tempfile.TemporaryDirectory() as temporary:
        admission, inventory = _fixture(Path(temporary))
        raw = Path(admission["execution_binding"]["raw_output_root"])
        stream_path = raw / inventory[1]["path"]
        raw.chmod(0o755)
        stream_path.chmod(0o644)
        payload = stream_fixture._stream(
            stream_fixture._header(basename=ATTEMPT_ID),
            [stream_fixture._event()],
        )
        stream_path.write_bytes(payload)
        inventory[1].update(sha256=_sha(payload), byte_count=len(payload))
        _seal(raw)
        with pytest.raises(
            successor.EscapeEvidenceAdmissionError, match="compact binary"
        ):
            _build(admission, inventory)
