#!/usr/bin/env python3
"""Adversarial tests for pure Q019 registered raw reduction."""

from __future__ import annotations

import copy
import math
import struct

import numpy as np
import pytest

from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as decks
from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_registered_raw_reduction_v1 as reduction
from tst.publication import q019_nonlinear_bell_particle_state as particle


CASE_ID = "q019-fr-runtime-initializer-ppc24-s0"


def _case() -> dict[str, object]:
    return next(row for row in decks.expected_cases() if row["case_id"] == CASE_ID)


def _parameter_header(
    snapshot_index: int, *, output1_last_time: str | None = None
) -> bytes:
    text = decks.render_deck(_case())
    for index in range(1, 14):
        if index == 13:
            text = text.replace(
                f"<output{index}>\n",
                f"<output{index}>\nlast_time = {max(0.0, 0.5 * snapshot_index):.17g}\n",
                1,
            )
        else:
            last_time = (
                output1_last_time
                if index == 1 and output1_last_time is not None
                else f"{max(0.0, 0.5 * snapshot_index):.17g}"
            )
            text = text.replace(
                f"<output{index}>\n",
                (
                    f"<output{index}>\n"
                    f"file_number = {snapshot_index + 1}\n"
                    f"last_time = {last_time}\n"
                ),
                1,
            )
    return text.encode("ascii")


def _binary_payload(
    *,
    fields: tuple[str, ...],
    values: dict[str, np.ndarray],
    cycle: int,
    time: float,
    snapshot_index: int,
    output1_last_time: str | None = None,
) -> bytes:
    case = _case()
    nx1, nx2, nx3 = (int(value) for value in case["nx"])
    mb1, mb2, mb3 = (int(value) for value in case["meshblock_nx"])
    ext1, ext2, ext3 = (float(value) for value in case["extents"])
    blocks = []
    for logical3 in range(nx3 // mb3):
        for logical2 in range(nx2 // mb2):
            for logical1 in range(nx1 // mb1):
                slices = (
                    slice(logical3 * mb3, (logical3 + 1) * mb3),
                    slice(logical2 * mb2, (logical2 + 1) * mb2),
                    slice(logical1 * mb1, (logical1 + 1) * mb1),
                )
                geometry = (
                    logical1 * ext1 / (nx1 // mb1),
                    (logical1 + 1) * ext1 / (nx1 // mb1),
                    logical2 * ext2 / (nx2 // mb2),
                    (logical2 + 1) * ext2 / (nx2 // mb2),
                    logical3 * ext3 / (nx3 // mb3),
                    (logical3 + 1) * ext3 / (nx3 // mb3),
                )
                arrays = np.concatenate(
                    [
                        np.asarray(values[field][slices], dtype="<f4").ravel()
                        for field in fields
                    ]
                )
                blocks.append(
                    struct.pack(
                        "<10i",
                        0,
                        mb1 - 1,
                        0,
                        mb2 - 1,
                        0,
                        mb3 - 1,
                        logical1,
                        logical2,
                        logical3,
                        0,
                    )
                    + struct.pack("<6d", *geometry)
                    + arrays.tobytes()
                )
    header = _parameter_header(
        snapshot_index, output1_last_time=output1_last_time
    )
    return (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={time:.17g}\n".encode()
        + f"  cycle={cycle}\n".encode()
        + b"  size of location=8\n"
        b"  size of variable=4\n"
        + f"  number of variables={len(fields)}\n".encode()
        + b"  variables:  "
        + b"  ".join(name.encode() for name in fields)
        + b"  \n"
        + f"  header offset={len(header)}\n".encode()
        + header
        + b"".join(blocks)
    )


def _fields(cycle: int) -> dict[str, np.ndarray]:
    case = _case()
    shape = tuple(reversed(tuple(int(value) for value in case["nx"])))
    values = {
        "dens": np.ones(shape),
        "eint": np.ones(shape),
        "velx": np.zeros(shape),
        "vely": np.zeros(shape),
        "velz": np.zeros(shape),
        "bcc1": np.ones(shape),
        "bcc2": np.full(shape, 1.0e-5 * (cycle + 1)),
        "bcc3": np.full(shape, 2.0e-5 * (cycle + 1)),
        "prtcl_rho": np.zeros(shape),
        "prtcl_jx": np.zeros(shape),
        "prtcl_jy": np.zeros(shape),
        "prtcl_jz": np.zeros(shape),
        "prtcl_dedt": np.zeros(shape),
        "prtcl_dpxdt": np.zeros(shape),
        "prtcl_dpydt": np.zeros(shape),
        "prtcl_dpzdt": np.zeros(shape),
        "prtcl_ebdot": np.zeros(shape),
    }
    if cycle > 0:
        values["prtcl_rho"].fill(float(case["rho_cr_over_rho0"]) * 10000.0)
        values["prtcl_jx"].fill(float(case["expected_j_over_c"]))
    return values


def _products(cycle: int, time: float, snapshot_index: int) -> dict[str, bytes]:
    values = _fields(cycle)
    products = {
        "mhd_w_bcc": _binary_payload(
            fields=reduction.MHD_FIELDS,
            values=values,
            cycle=cycle,
            time=time,
            snapshot_index=snapshot_index,
        )
    }
    for field in reduction.MOMENT_FIELDS:
        products[field] = _binary_payload(
            fields=(field,),
            values=values,
            cycle=cycle,
            time=time,
            snapshot_index=snapshot_index,
        )
    return products


def _restart_payload() -> bytes:
    case = _case()
    count = 6
    real_rows = np.zeros((count, 26), dtype="<f8")
    speed = float(case["guide_parallel_stream_speed"])
    real_rows[:, 1] = speed
    real_rows[:, 6] = 10000.0
    real_rows[:, 7] = 1.0
    real_rows[:, 22] = 1.0
    integer_rows = np.zeros((count, 4), dtype="<i4")
    integer_rows[:, 1] = np.arange(count)
    metadata = struct.pack(
        "<15i", 7, 1, 26, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0
    )
    model_ints = [0] * 31
    model_reals = [0.0] * 37
    model_reals[20] = float(case["deposit_qscale"])
    return (
        b"<job>\nbasename=q019-raw-reduction-fixture\n<par_end>\n"
        + struct.pack("<Q", particle.restart_layout.PIC_RESTART_MAGIC)
        + metadata
        + struct.pack("<d", 1.0e8)
        + struct.pack("<31i", *model_ints)
        + struct.pack("<37d", *model_reals)
        + struct.pack("<Q", count)
        + struct.pack("<i", count)
        + real_rows.tobytes()
        + integer_rows.tobytes()
    )


def _checkpoint(cycle: int, time: float, index: int) -> dict[str, object]:
    return {
        "cycle": cycle,
        "time": time,
        "binary_products": _products(cycle, time, index),
        "restart_path": f"raw/rst/{CASE_ID}.{index:05d}.rst",
        "restart_payload": _restart_payload(),
    }


def test_matched_checkpoint_reduction_builds_analyzer_inputs() -> None:
    report = reduction.reduce_matched_checkpoints(
        CASE_ID, [_checkpoint(0, 0.0, 0), _checkpoint(7, 0.5, 1)]
    )
    assert report["record_type"] == reduction.RECORD_TYPE
    assert report["matched_checkpoint_count"] == 2
    assert report["chronology"] == [
        {"cycle": 0, "time": 0.0},
        {"cycle": 7, "time": 0.5},
    ]
    assert report["authority"]["qualification_authorized"] is False
    assert report["particle_states"][0]["raw_restart_binding"]["path"].endswith(
        ".00000.rst"
    )
    assert report["particle_reductions"][0]["conservation"]["reference_bound"]
    assert math.isclose(
        report["particle_reductions"][0]["conservation"][
            "relative_energy_residual"
        ],
        0.0,
        abs_tol=1.0e-15,
    )


def test_cross_product_cycle_drift_fails_closed() -> None:
    products = _products(0, 0.0, 0)
    products["prtcl_jx"] = _binary_payload(
        fields=("prtcl_jx",),
        values=_fields(0),
        cycle=1,
        time=0.0,
        snapshot_index=0,
    )
    with pytest.raises(reduction.RawReductionError, match="metadata differs"):
        reduction.compose_snapshot(CASE_ID, products)


def test_embedded_runtime_semantics_drift_fails_closed() -> None:
    products = _products(0, 0.0, 0)
    products["mhd_w_bcc"] = products["mhd_w_bcc"].replace(
        b"reconstruct = plm", b"reconstruct = ppm", 1
    )
    with pytest.raises(
        reduction.RawReductionError, match="immutable matrix row"
    ):
        reduction.compose_snapshot(CASE_ID, products)


def _controller_parameters(artifact_id: str) -> tuple[dict[str, object], dict[str, dict[str, str]]]:
    overlay = next(
        item
        for item in controller.expected_overlays()
        if item["artifact_id"] == artifact_id
    )
    parameters = decks.parse_athinput_text(controller.render_overlay(overlay))
    parameters[controller.CONTROLLER_BLOCK].update(
        {
            "runtime_resolution_samples": "0",
            "runtime_resolution_last_cycle": "0",
            "runtime_resolution_last_time": "0",
            "runtime_resolution_last_B_over_B0": "-1",
            "runtime_resolution_max_B_over_B0": "-1",
            "runtime_controller_triggered": "false",
            "runtime_controller_trigger_failure": "false",
            "runtime_controller_trigger_reason": "0",
            "runtime_controller_trigger_cycle": "0",
            "runtime_controller_trigger_time": "0",
            "runtime_controller_trigger_metric": "-1",
        }
    )
    return overlay, parameters


def test_exact_runtime_controller_overlay_preserves_base_matrix_identity() -> None:
    overlay, parameters = _controller_parameters(
        "q019-controller-pilot-2d-instrumented"
    )
    case = next(
        row
        for row in decks.expected_cases()
        if row["case_id"] == overlay["source_case_id"]
    )
    base, profile = reduction._execution_profile(parameters, case=case)
    assert decks.deck_semantics_payload(base) == decks.matrix_identity_payload(case)
    assert profile["artifact_id"] == overlay["artifact_id"]
    assert profile["authority"] == "excluded_pilot_only"
    assert profile["runtime_state"]["runtime_resolution_samples"] == 0
    assert profile["saturation_evidence_eligible"] is False


def test_runtime_controller_overlay_and_mutable_state_drift_fail_closed() -> None:
    overlay, parameters = _controller_parameters(
        "q019-controller-pilot-2d-instrumented"
    )
    case = next(
        row
        for row in decks.expected_cases()
        if row["case_id"] == overlay["source_case_id"]
    )
    hostile = copy.deepcopy(parameters)
    hostile[controller.CONTROLLER_BLOCK]["pilot_cycle_limit"] = "21"
    with pytest.raises(
        reduction.RawReductionError, match="exact checked-in contract"
    ):
        reduction._execution_profile(hostile, case=case)

    partial = copy.deepcopy(parameters)
    del partial[controller.CONTROLLER_BLOCK]["runtime_resolution_samples"]
    with pytest.raises(reduction.RawReductionError, match="inventory drifted"):
        reduction._execution_profile(partial, case=case)

    inconsistent = copy.deepcopy(parameters)
    inconsistent[controller.CONTROLLER_BLOCK][
        "runtime_controller_triggered"
    ] = "true"
    with pytest.raises(reduction.RawReductionError, match="trigger state drifted"):
        reduction._execution_profile(inconsistent, case=case)


def test_invalid_mutable_output_state_fails_closed() -> None:
    products = _products(0, 0.0, 0)
    products["mhd_w_bcc"] = _binary_payload(
        fields=reduction.MHD_FIELDS,
        values=_fields(0),
        cycle=0,
        time=0.0,
        snapshot_index=0,
        output1_last_time="-2",
    )
    with pytest.raises(
        reduction.RawReductionError,
        match="neither the unwritten sentinel nor a nonnegative time",
    ):
        reduction.compose_snapshot(CASE_ID, products)


def test_missing_product_and_restart_alias_fail_closed() -> None:
    products = _products(0, 0.0, 0)
    products.pop("prtcl_ebdot")
    with pytest.raises(reduction.RawReductionError, match="inventory drifted"):
        reduction.compose_snapshot(CASE_ID, products)
    checkpoints = [_checkpoint(0, 0.0, 0), _checkpoint(7, 0.5, 1)]
    altered = copy.deepcopy(checkpoints)
    altered[1]["restart_path"] = "../substituted.rst"
    with pytest.raises(reduction.RawReductionError, match="path is malformed"):
        reduction.reduce_matched_checkpoints(CASE_ID, altered)


def test_chronology_and_projection_drift_fail_closed() -> None:
    with pytest.raises(reduction.RawReductionError, match="chronology drifted"):
        reduction.reduce_matched_checkpoints(
            CASE_ID, [_checkpoint(1, 0.1, 0), _checkpoint(2, 0.5, 1)]
        )
    second = _checkpoint(7, 0.5, 1)
    second["time"] = 0.6
    with pytest.raises(reduction.RawReductionError, match="projection differs"):
        reduction.reduce_matched_checkpoints(
            CASE_ID, [_checkpoint(0, 0.0, 0), second]
        )
