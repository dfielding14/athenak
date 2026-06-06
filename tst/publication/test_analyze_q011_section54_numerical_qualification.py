#!/usr/bin/env python3
"""Synthetic tests for source-local Q-011 Section 5.4 numerical aggregation."""

from __future__ import annotations

import copy
from contextlib import contextmanager
import hashlib
import os
from pathlib import Path
import sys
import tempfile
import time
from typing import Iterator
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "frontier_control_plane"))

from . import analyze_q011_section54_campaign as admission
from . import analyze_q011_section54_numerical_qualification as qualifier
from . import q011_section54_artifacts as artifacts
from . import q011_section54_particles as particles
from . import test_analyze_q011_section54_campaign as campaign_fixture


def _sha256(character: str) -> str:
    return character * 64


@contextmanager
def _authoritative_pressure_pilot_consumer_fixture(
    authorized_root: Path,
) -> Iterator[mock.Mock]:
    expected_path = (
        authorized_root
        / "publication/q011_section54_pressure_pilot_publication_receipt.json"
    )

    def consume() -> dict[str, object]:
        return {
            "packet_receipt_sha256": "d" * 64,
            "aggregate_receipt_sha256": hashlib.sha256(
                expected_path.read_bytes()
            ).hexdigest(),
            "manifest_sha256": "b" * 64,
            "analysis_result_sha256": "c" * 64,
            "status": "pass_engineering_calibration_only",
        }

    with mock.patch.object(
        qualifier.admission.pressure_selection.historical_pressure_pilot_consumer,
        "consume_exact_historical_production_pressure_pilot",
        side_effect=consume,
    ) as consumer:
        yield consumer


def _identity(index: int, variant: str, seed: int) -> dict[str, object]:
    return {
        "variant": variant,
        "seed": seed,
        "physical_mode": qualifier.PHYSICAL_MODE,
        "attempt_id": f"baseline-{index:03d}-{variant}-seed-{seed}",
    }


def _pairing(identity: dict[str, object]) -> dict[str, object]:
    counterpart = None
    if identity["variant"] == qualifier.AMR_VARIANT:
        counterpart = {"variant": qualifier.FINE_VARIANT, "seed": identity["seed"]}
    elif identity["variant"] == qualifier.FINE_VARIANT:
        counterpart = {"variant": qualifier.AMR_VARIANT, "seed": identity["seed"]}
    return {
        "paired_seed_rule": qualifier.PAIRED_SEED_RULE,
        "pair_key": {"seed": identity["seed"]},
        "amr_fine_uniform_counterpart": counterpart,
        "comparison_status": "schema_wired_not_evaluated_by_artifact_admission_slice",
    }


def _tail_spectrum() -> dict[str, object]:
    edges = np.asarray(particles.CHI_BIN_EDGES)
    centers = np.sqrt(edges[:-1] * edges[1:])
    selected = (centers >= 20.0) & (centers <= 160.0)
    chi = centers[selected]
    macro_weight = 7.0 * chi**-1.5 * np.diff(edges)[selected]
    return particles.weighted_spectrum_record(chi, macro_weight)


def _particle_reductions(*, observed_t500: float = 500.0) -> dict[str, object]:
    early = particles.weighted_spectrum_record([2.0], [1.0])
    late = _tail_spectrum()
    return {
        "t500": {
            "schema_version": particles.SCHEMA_VERSION,
            "record_type": "q011_section54_particle_snapshot_reduction",
            "nominal_slot_time": 500.0,
            "observed_committed_time": observed_t500,
            "snapshot_time_omega0_inverse": observed_t500,
            "weighted_spectrum": early,
        },
        "t1200": {
            "schema_version": particles.SCHEMA_VERSION,
            "record_type": "q011_section54_particle_snapshot_reduction",
            "nominal_slot_time": 1200.0,
            "observed_committed_time": 1200.0,
            "snapshot_time_omega0_inverse": 1200.0,
            "weighted_spectrum": late,
            "late_slope": particles.late_slope_record(late["f_chi"]),
        },
    }


def _spatial_reduction(
    *, passes: bool = True, observed_t500: float = 500.0
) -> dict[str, object]:
    mean = 2.0 if passes else 1.0
    return {
        "schema_version": 1,
        "record_type": "q011_section54_t500_spatial_reduction",
        "nominal_slot_time": 500.0,
        "observed_committed_time": observed_t500,
        "upstream_b_amplification": {
            "nominal_slot_time": 500.0,
            "observed_committed_time": observed_t500,
            "x_ideal_c_over_omega_pi": 0.0,
            "upstream_window_c_over_omega_pi": [120.0, 1200.0],
            "selected_cell_count": 1,
            "selected_area": 1.0,
            "mean_magnetic_magnitude": mean,
            "reference_b0": 1.0,
            "amplification_over_b0": mean,
            "acceptance_range": [1.2, 3.5],
            "passes_gate": passes,
        },
    }


def _paired_spatial_inputs(
    *, center_offset: float = 0.0, observed_t500: float = 500.0
) -> dict[str, object]:
    centers = [-4.5 + center_offset, -1.5 + center_offset, 1.5 + center_offset, 4.5 + center_offset]
    profile = {
        "nominal_slot_time": 500.0,
        "observed_committed_time": observed_t500,
        "x1_centers_c_over_omega_pi": centers,
        "values_x": [1.0, 2.0, 3.0, 4.0],
        "column_areas": [3.0, 3.0, 3.0, 3.0],
    }
    return {
        "nominal_slot_time": 500.0,
        "observed_committed_time": observed_t500,
        "detected_front": {"x_front_c_over_omega_pi": 12.0},
        "upstream_b_amplification": {"amplification_over_b0": 2.0},
        "y_area_weighted_profiles": {
            "rho": {**profile, "quantity": "rho", "source_field": "dens"},
            "bmag": {**profile, "quantity": "bmag", "source_field": "bmag"},
        },
    }


def _admission(identity: dict[str, object], inventory_sha256: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "record_type": admission.RESULT_RECORD_TYPE,
        "campaign_id": admission.CAMPAIGN_ID,
        "qualification_scope": admission.QUALIFICATION_SCOPE,
        "admitted_for_follow_on_numerical_qualification": True,
        "final_claim_closure": False,
        "status": "admitted_for_follow_on_numerical_qualification",
        "failure_reasons": [],
        "admission": {
            "run_identity": copy.deepcopy(identity),
            "preregistration_binding": {
                "sha256": admission.EXPECTED_PREREGISTRATION_SHA256,
                "expected_sha256": admission.EXPECTED_PREREGISTRATION_SHA256,
            },
            "immutable_tree": {"inventory_sha256": inventory_sha256},
            "amr_pairing": _pairing(identity),
            "numerical_qualification_status": "not_evaluated_by_artifact_admission_slice",
        },
    }


def _source_checkpoint_lineage(identity: dict[str, object]) -> dict[str, object]:
    return {
        "nominal_slot_time": 500.0,
        "observed_committed_time": 500.0,
        "retained_attempt_id": identity["attempt_id"],
        "restart_manifest_path": "rst/q011.00500.rst.manifest",
        "restart_member_path": "rst/q011.00500.rst",
        "retained_restart_member_absolute_path": (
            f"/retained/{identity['attempt_id']}/rst/q011.00500.rst"
        ),
        "restart_member_sha256": hashlib.sha256(b"checkpoint").hexdigest(),
    }


def _attempts() -> list[dict[str, object]]:
    attempts = []
    for index, (variant, seed) in enumerate(qualifier.EXPECTED_MATRIX_CELLS):
        inventory_sha256 = f"{index + 1:064x}"
        identity = _identity(index, variant, seed)
        attempts.append(
            qualifier._bind_unit_only_canonical_attempt(
                {
                    "schema_version": 1,
                    "record_type": qualifier.ATTEMPT_RECORD_TYPE,
                    "run_identity": identity,
                    "raw_inventory_sha256": inventory_sha256,
                    "source_checkpoint_lineage": _source_checkpoint_lineage(identity),
                    "admission_result": _admission(identity, inventory_sha256),
                    "particle_reductions": _particle_reductions(),
                    "spatial_reduction": _spatial_reduction(),
                }
            )
        )
    return attempts


def _retained_attempt_bindings(
    attempts: list[dict[str, object]],
) -> list[qualifier._RetainedAttemptBinding]:
    return [
        qualifier._RetainedAttemptBinding(
            wrapper_payload=qualifier._json_bytes(wrapper),
            _key=qualifier._PRODUCTION_BINDING_KEY,
        )
        for wrapper in attempts
    ]


def _attempt_tree_descriptors() -> list[dict[str, object]]:
    return [
        {
            "campaign_root": f"/retained/attempt-{index:03d}",
            "expected_inventory_sha256": f"{index + 1:064x}",
        }
        for index in range(qualifier.EXPECTED_BASELINE_ATTEMPTS)
    ]


def _recompute_descriptors() -> list[dict[str, object]]:
    return [
        {
            "bundle_root": f"/retained/pair-{index:03d}",
            "expected_inventory_sha256": f"{index + 101:064x}",
        }
        for index in range(len(qualifier.QUALIFYING_SEEDS))
    ]


def _retained_pair_bindings(
    attempts: list[dict[str, object]],
) -> list[qualifier._RetainedPairResultBinding]:
    by_sha256 = {wrapper["attempt_sha256"]: wrapper for wrapper in attempts}
    return [
        qualifier._RetainedPairResultBinding(
            record_payload=qualifier._json_bytes(
                {
                    "schema_version": qualifier.SCHEMA_VERSION,
                    "record_type": qualifier.PAIR_RECOMPUTE_RECORD_TYPE,
                    "seed": pair["seed"],
                    "amr_attempt_sha256": pair["amr_attempt_sha256"],
                    "amr_raw_inventory_sha256": by_sha256[
                        pair["amr_attempt_sha256"]
                    ]["attempt"]["raw_inventory_sha256"],
                    "fine_uniform_attempt_sha256": pair[
                        "fine_uniform_attempt_sha256"
                    ],
                    "fine_uniform_raw_inventory_sha256": by_sha256[
                        pair["fine_uniform_attempt_sha256"]
                    ]["attempt"]["raw_inventory_sha256"],
                    "residuals": pair["residuals"],
                }
            ),
            _key=qualifier._PRODUCTION_BINDING_KEY,
        )
        for pair in _pairs(attempts)
    ]


def _retained_restart_binding(
    attempts: list[dict[str, object]],
) -> qualifier._RetainedRestartParityBinding:
    binding = _restart_binding(attempts)
    by_sha256 = {wrapper["attempt_sha256"]: wrapper for wrapper in attempts}
    return qualifier._RetainedRestartParityBinding(
        record_payload=qualifier._json_bytes(
            {
                "schema_version": qualifier.SCHEMA_VERSION,
                "record_type": qualifier.RESTART_RECOMPUTE_RECORD_TYPE,
                "source_attempt_sha256": binding["source_attempt_sha256"],
                "raw_inventory_sha256": by_sha256[
                    binding["source_attempt_sha256"]
                ]["attempt"]["raw_inventory_sha256"],
                "planner_materialization_receipt_member": (
                    "planner/materialization_receipt.json"
                ),
                "retained_campaign_plan_member": "planner/campaign_plan.json",
                "planned_restart_carrier": {
                    "path": qualifier.PLANNED_RESTART_CARRIER_PATH,
                    "sha256": _sha256("1"),
                },
                "source_checkpoint_member": "source/checkpoint.rst",
                "source_checkpoint_sha256": _sha256("2"),
                "source_checkpoint_lineage": by_sha256[
                    binding["source_attempt_sha256"]
                ]["attempt"]["source_checkpoint_lineage"],
                "screen_scope": qualifier.RESTART_SCREEN_SCOPE,
                "uninterrupted": binding["uninterrupted"],
                "continued": binding["continued"],
            }
        ),
        _key=qualifier._PRODUCTION_BINDING_KEY,
    )


def _make_tree_writable(root: Path) -> None:
    for directory, _, filenames in os.walk(root):
        os.chmod(directory, 0o700)
        for filename in filenames:
            os.chmod(Path(directory) / filename, 0o600)


def _publish_recompute_bundle(
    temporary: tempfile.TemporaryDirectory[str],
    *,
    name: str,
    raw_inventory_sha256: str,
    analyzer_paths: tuple[str, ...],
    member: str,
    record: dict[str, object],
    extra_artifacts: dict[str, bytes] | None = None,
) -> tuple[str, str]:
    root = Path(temporary.name) / name
    manifest = artifacts.build_derived_manifest(
        campaign_id=admission.CAMPAIGN_ID,
        raw_artifact_inventory_sha256=raw_inventory_sha256,
        analyzer_bindings={
            analyzer_path: hashlib.sha256(
                (Path(__file__).resolve().parents[2] / analyzer_path).read_bytes()
            ).hexdigest()
            for analyzer_path in analyzer_paths
        },
        reviewer_disposition=_pending(),
    )
    payloads = {
        member: artifacts.canonical_json_bytes(record),
        **(extra_artifacts or {}),
        artifacts.MANIFEST_NAME: artifacts.canonical_json_bytes(manifest),
    }
    inventory_payload = artifacts.canonical_json_bytes(
        {
            "record_type": "q011_section54_derived_artifact_inventory",
            "schema_version": 1,
            "members": [
                {
                    "path": relative,
                    "size": len(payload),
                    "sha256": hashlib.sha256(payload).hexdigest(),
                }
                for relative, payload in sorted(payloads.items())
            ],
        }
    )
    payloads[artifacts.INVENTORY_NAME] = inventory_payload
    for relative, payload in payloads.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        path.chmod(0o444)
    for directory, _, _ in os.walk(root, topdown=False):
        os.chmod(directory, 0o555)
    return str(root), hashlib.sha256(inventory_payload).hexdigest()


def _restart_source(prefix: str, *, branch_role: str) -> dict[str, object]:
    def member(suffix: str) -> dict[str, str]:
        return {
            "bundle_member": f"{prefix}/{suffix}",
            "structured_artifact_member": f"raw/{prefix}/{suffix}",
            "sha256": hashlib.sha256(b"retained").hexdigest(),
        }

    return {
        "branch_role": branch_role,
        "execution_receipt_authority_path": (
            f"/authority/{prefix}/analysis/"
            f"{admission.REGISTERED_EXECUTION_RECEIPT_NAME}"
        ),
        "execution_receipt_member": f"{prefix}/execution_receipt.json",
        "raw_output_root": f"/retained/{prefix}/raw",
        "structured_artifact_inventory_sha256": _sha256("a"),
        "outputs_after_checkpoint": [
            {
                "nominal_slot_omega0_inverse": time,
                "observed_committed_cycle": 5000 + int(time),
                "observed_committed_time_omega0_inverse": time + 0.05,
                "members": {
                    "rho_bin": member(f"{time:.1f}.rho.bin"),
                    "bmag_bin": member(f"{time:.1f}.bmag.bin"),
                    "prtcl_jx_bin": member(f"{time:.1f}.prtcl_jx.bin"),
                    "j2_bin": member(f"{time:.1f}.j2.bin"),
                    "prtcl_all_pvtk": member(f"{time:.1f}.prtcl_all.vtk"),
                },
            }
            for time in (600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0)
        ],
    }


def _planned_execution_inventory() -> dict[str, str]:
    return {
        "planned_restart_carrier_sha256": _sha256("1"),
        "clean_candidate_manifest_sha256": _sha256("2"),
        "executable_sha256": _sha256("3"),
        "paper_deck_sha256": _sha256("4"),
        "environment_profile_sha256": _sha256("5"),
        "control_plane_version": _sha256("6"),
        "git_commit": "a" * 40,
        "restart_carrier_id": "amr-restart-continuation-seed-23050101",
        "restart_authorized_orion_attempt_root": (
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/restart"
        ),
    }


def _shared_execution_receipt(
    *,
    attempt_id: str,
    suffix: str,
    argv: list[str],
    raw_output_root: str,
    artifact_dir: str,
) -> dict[str, object]:
    return {
        "record_type": "q011_section54_reconciled_registered_execution_receipt",
        "schema_version": 1,
        "receipt_role": "immutable_reconciled_registered_execution",
        "registration_scope": "registered_science",
        "reconciled": True,
        "reservation_id": f"10000000-0000-4000-8000-0000000000{suffix}",
        "submission_id": f"00000000-0000-4000-8000-0000000000{suffix}",
        "reconciliation_event_sha256": suffix[-1] * 64,
        "attempt_id": attempt_id,
        "source_commit": _planned_execution_inventory()["git_commit"],
        "executable_sha256": _planned_execution_inventory()["executable_sha256"],
        "deck_sha256": _planned_execution_inventory()["paper_deck_sha256"],
        "environment_sha256": _planned_execution_inventory()[
            "environment_profile_sha256"
        ],
        "control_plane_version": _planned_execution_inventory()[
            "control_plane_version"
        ],
        "argv": argv,
        "slurm_job_id": str(int(suffix)),
        "slurm_terminal_state": "COMPLETED",
        "raw_output_root": raw_output_root,
        "artifact_dir": artifact_dir,
        "planner_retention": {
            "schema_version": 1,
            "retention_role": "q011_section54_deterministic_retained_attempt",
            "attempt_id": attempt_id,
            "authorized_orion_attempt_root": str(Path(raw_output_root).parent),
            "authorized_orion_raw_root": raw_output_root,
            "argv": argv,
        },
        "pre_submit_manifest_sha256": "e" * 64,
    }


def _authoritative_receipt_projection(suffix: str) -> dict[str, str]:
    return {
        "submission_id": f"00000000-0000-4000-8000-0000000000{suffix}",
        "reservation_id": f"10000000-0000-4000-8000-0000000000{suffix}",
        "slurm_job_id": str(int(suffix)),
        "reconciliation_event_sha256": suffix[-1] * 64,
    }


@contextmanager
def _pinned_receipt_fixture(
    path: Path, payload: bytes
) -> Iterator[tuple[Path, bytes, int]]:
    yield path, payload, -1


@contextmanager
def _ledger_snapshot_fixture(
    binding: dict[str, str] | None = None,
) -> Iterator[dict[str, str]]:
    if binding is None:
        raise admission.QualificationError(
            "fixture mirrored ledger drifted",
            code="registered_execution_receipt_ledger_drift",
        )
    yield binding


def _authoritative_receipt_side_effect(
    _branch: dict[str, object],
    _payloads: dict[str, bytes],
    *,
    branch_role: str,
    **_kwargs: object,
) -> dict[str, str]:
    return _authoritative_receipt_projection(
        "01" if branch_role == "uninterrupted_baseline" else "02"
    )


@contextmanager
def _authoritative_receipt_context_side_effect(
    _branch: dict[str, object],
    _payloads: dict[str, bytes],
    *,
    branch_role: str,
    **_kwargs: object,
) -> Iterator[tuple[dict[str, str], int]]:
    yield _authoritative_receipt_projection(
        "01" if branch_role == "uninterrupted_baseline" else "02"
    ), -1


def _retained_restart_fixture(
    source: dict[str, object],
) -> tuple[dict[str, object], dict[str, bytes], dict[str, str]]:
    checkpoint_payload = b"checkpoint"
    checkpoint_sha256 = hashlib.sha256(checkpoint_payload).hexdigest()
    uninterrupted = _restart_source(
        "uninterrupted", branch_role="uninterrupted_baseline"
    )
    continued = _restart_source(
        "continued", branch_role="checkpoint_restart_continuation"
    )
    record = {
        "schema_version": qualifier.SCHEMA_VERSION,
        "record_type": qualifier.RESTART_RECOMPUTE_RECORD_TYPE,
        "source_attempt_sha256": source["attempt_sha256"],
        "raw_inventory_sha256": source["attempt"]["raw_inventory_sha256"],
        "planner_materialization_receipt_member": "planner/materialization_receipt.json",
        "retained_campaign_plan_member": "planner/campaign_plan.json",
        "planned_restart_carrier": {
            "path": qualifier.PLANNED_RESTART_CARRIER_PATH,
            "sha256": _planned_execution_inventory()[
                "planned_restart_carrier_sha256"
            ],
        },
        "source_checkpoint_member": "source/checkpoint.rst",
        "source_checkpoint_sha256": checkpoint_sha256,
        "source_checkpoint_lineage": source["attempt"]["source_checkpoint_lineage"],
        "screen_scope": qualifier.RESTART_SCREEN_SCOPE,
        "uninterrupted": uninterrupted,
        "continued": continued,
    }
    artifacts_by_member = {
        "planner/materialization_receipt.json": artifacts.canonical_json_bytes({}),
        "planner/campaign_plan.json": artifacts.canonical_json_bytes({}),
        "source/checkpoint.rst": checkpoint_payload,
        uninterrupted["execution_receipt_member"]: artifacts.canonical_json_bytes(
            _shared_execution_receipt(
                attempt_id=source["attempt"]["run_identity"]["attempt_id"],
                suffix="01",
                argv=["baseline"],
                raw_output_root="/retained/baseline/raw",
                artifact_dir="/runs/baseline",
            )
        ),
        continued["execution_receipt_member"]: artifacts.canonical_json_bytes(
            _shared_execution_receipt(
                attempt_id=_planned_execution_inventory()["restart_carrier_id"],
                suffix="02",
                argv=["restart"],
                raw_output_root="/retained/restart/raw",
                artifact_dir="/runs/restart",
            )
        ),
    }
    for branch in (uninterrupted, continued):
        for output in branch["outputs_after_checkpoint"]:
            artifacts_by_member.update(
                {
                    member["bundle_member"]: b"retained"
                    for member in output["members"].values()
                }
            )
    return record, artifacts_by_member, _planned_execution_inventory()


def _planner_restart_fixture(
    source: dict[str, object],
) -> tuple[
    dict[str, object],
    dict[str, bytes],
    dict[str, bytes],
    dict[str, object],
    dict[str, object],
]:
    identity = source["attempt"]["run_identity"]
    candidate = {
        "clean_candidate_manifest": {
            "path": "/frozen/clean_candidate_manifest.json",
            "sha256": _planned_execution_inventory()[
                "clean_candidate_manifest_sha256"
            ],
        },
        "executable": {
            "path": "/frozen/athena",
            "sha256": _planned_execution_inventory()["executable_sha256"],
        },
        "environment_profile": {
            "path": "/frozen/environment.sh",
            "sha256": _planned_execution_inventory()["environment_profile_sha256"],
            "control_plane_version": _planned_execution_inventory()[
                "control_plane_version"
            ],
            "reviewed_source": {
                "path": "reviewed/environment.sh",
                "sha256": _sha256("9"),
            },
        },
        "git_commit": _planned_execution_inventory()["git_commit"],
    }
    source_bindings = {
        "restart_preregistration": {
            "path": "bindings/restart_preregistration.json",
            "sha256": admission.EXPECTED_RESTART_PREREGISTRATION_SHA256,
        },
        "paper_deck": {
            "path": "bindings/paper_deck.athinput",
            "sha256": _planned_execution_inventory()["paper_deck_sha256"],
        },
    }
    policy = qualifier.restart.load_preregistration()
    campaign_root = Path(
        "/lustre/orion/ast207/proj-shared/dfielding/PIC/campaigns/q011-fixture"
    )
    carrier_id = f"amr-restart-continuation-seed-{identity['seed']}"
    artifact_root = campaign_root / "restart_continuation" / carrier_id
    descriptor = {
        "attempt_id": identity["attempt_id"],
        "variant": identity["variant"],
        "qualifying_seed": identity["seed"],
        "selected_problem_ps_p0": 0.1,
    }
    descriptor_payload = admission.campaign_planner._json_bytes(descriptor)
    descriptor_binding = {
        "path": f"attempts/{identity['attempt_id']}.json",
        "sha256": hashlib.sha256(descriptor_payload).hexdigest(),
    }
    launch_contract = admission.campaign_planner._restart_launch_contract(
        carrier_id=carrier_id,
        source_attempt=descriptor,
        restart_preregistration=policy,
        candidate=candidate,
        paper_deck_binding=source_bindings["paper_deck"],
        artifact_root=artifact_root,
    )
    launch_payload = admission.campaign_planner._json_bytes(launch_contract)
    launch_binding = {
        "path": f"launch_contracts/restart_continuation/{carrier_id}.json",
        "sha256": hashlib.sha256(launch_payload).hexdigest(),
    }
    continuation = policy["continuation_contract"]
    carrier = {
        "record_type": admission.campaign_planner.RESTART_CARRIER_RECORD_TYPE,
        "schema_version": 1,
        "carrier_id": carrier_id,
        "status": "planned_not_authorized",
        "source_baseline_attempt_id": identity["attempt_id"],
        "variant": identity["variant"],
        "qualifying_seed": identity["seed"],
        "selected_problem_ps_p0": 0.1,
        "authorized_orion_attempt_root": str(artifact_root),
        "restart_preregistration": source_bindings["restart_preregistration"],
        "checkpoint_nominal_slot_omega0_inverse": continuation[
            "checkpoint_nominal_slot_omega0_inverse"
        ],
        "checkpoint_observed_commit_binding_required": True,
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": continuation[
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
        ],
        "retained_output_pairing_policy": continuation[
            "retained_output_pairing_policy"
        ],
        "comparison_tolerances_max_absolute_difference": continuation[
            "comparison_tolerances_max_absolute_difference"
        ],
        "launch_contract": launch_binding,
    }
    carrier_payload = admission.campaign_planner._json_bytes(carrier)
    carrier_binding = {
        "path": qualifier.PLANNED_RESTART_CARRIER_PATH,
        "sha256": hashlib.sha256(carrier_payload).hexdigest(),
    }
    plan_payload = artifacts.canonical_json_bytes(
        {
            "authorized_orion_campaign_root": str(campaign_root),
            "baseline_attempt_descriptors": [descriptor_binding],
            "campaign_matrix": {"qualifying_seeds": list(qualifier.QUALIFYING_SEEDS)},
            "candidate_binding": candidate,
            "selected_pressure": {"selected_case": {"problem_ps_p0": 0.1}},
            "source_bindings": source_bindings,
            "restart_continuation_carrier": carrier_binding,
        }
    )
    record = {
        "planner_materialization_receipt_member": "planner/materialization_receipt.json",
        "retained_campaign_plan_member": "planner/campaign_plan.json",
        "planned_restart_carrier": carrier_binding,
    }
    bundle_payloads = {
        "planner/materialization_receipt.json": artifacts.canonical_json_bytes({}),
        "planner/campaign_plan.json": plan_payload,
    }
    tree_payloads = {
        descriptor_binding["path"]: descriptor_payload,
        carrier_binding["path"]: carrier_payload,
        launch_binding["path"]: launch_payload,
    }
    parsed_source = {"identity": copy.deepcopy(identity)}
    source_attempt_payload = copy.deepcopy(source["attempt"])
    source_attempt_payload["admission_result"]["admission"]["artifact_bindings"] = {
        "campaign_plan": {"sha256": hashlib.sha256(plan_payload).hexdigest()},
        "planner_materialization_receipt": {
            "sha256": hashlib.sha256(
                bundle_payloads["planner/materialization_receipt.json"]
            ).hexdigest()
        },
    }
    return record, bundle_payloads, tree_payloads, parsed_source, source_attempt_payload


def _residuals() -> list[dict[str, object]]:
    return [
        {
            "observable": observable,
            "maximum_absolute_difference": 0.0,
            "relative_mean_absolute": None if mean_bound is None else 0.0,
            "relative_root_mean_square": None if rms_bound is None else 0.0,
        }
        for observable, _, mean_bound, rms_bound in qualifier.PAIRED_RESIDUAL_THRESHOLDS
    ]


def _pair_record(
    source: dict[str, object], counterpart: dict[str, object]
) -> dict[str, object]:
    return {
        "schema_version": qualifier.SCHEMA_VERSION,
        "record_type": qualifier.PAIR_RECOMPUTE_RECORD_TYPE,
        "seed": qualifier.QUALIFYING_SEEDS[0],
        "amr_attempt_sha256": source["attempt_sha256"],
        "amr_raw_inventory_sha256": source["attempt"]["raw_inventory_sha256"],
        "fine_uniform_attempt_sha256": counterpart["attempt_sha256"],
        "fine_uniform_raw_inventory_sha256": counterpart["attempt"][
            "raw_inventory_sha256"
        ],
        "residuals": _residuals(),
    }


def _pairs(attempts: list[dict[str, object]]) -> list[dict[str, object]]:
    by_cell = {
        (
            wrapper["attempt"]["run_identity"]["variant"],
            wrapper["attempt"]["run_identity"]["seed"],
        ): wrapper["attempt_sha256"]
        for wrapper in attempts
    }
    return [
        {
            "seed": seed,
            "amr_attempt_sha256": by_cell[(qualifier.AMR_VARIANT, seed)],
            "fine_uniform_attempt_sha256": by_cell[(qualifier.FINE_VARIANT, seed)],
            "residuals": _residuals(),
        }
        for seed in qualifier.QUALIFYING_SEEDS
    ]


def _restart_binding(attempts: list[dict[str, object]]) -> dict[str, object]:
    source = next(
        wrapper["attempt_sha256"]
        for wrapper in attempts
        if wrapper["attempt"]["run_identity"]["variant"] == qualifier.AMR_VARIANT
    )
    return {
        "source_attempt_sha256": source,
        "uninterrupted": {"synthetic": "uninterrupted"},
        "continued": {"synthetic": "continued"},
    }


def _pending() -> dict[str, object]:
    return {
        "status": "pending_external_review",
        "reviewer": None,
        "reviewed_at_utc": None,
        "notes": "Synthetic aggregate remains pending named external review.",
    }


def _parity_result() -> dict[str, object]:
    return {
        "result": "pass_deterministic_continuation_parity",
        "checkpoint_nominal_slot_omega0_inverse": 500.0,
        "checkpoint_observed_committed_cycle": 5500,
        "checkpoint_observed_committed_time_omega0_inverse": 500.0,
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": [
            600.0,
            700.0,
            800.0,
            900.0,
            1000.0,
            1100.0,
            1200.0,
        ],
        "paired_output_observed_commits": [],
        "maximum_absolute_difference_by_field": {"rho_bin": 0.0},
    }


def _qualify(
    attempts: list[dict[str, object]],
    *,
    pairs: list[dict[str, object]] | None = None,
    inventories: list[str] | None = None,
) -> dict[str, object]:
    with mock.patch.object(
        qualifier.restart,
        "compare_deterministic_continuation_parity",
        return_value=_parity_result(),
    ):
        return qualifier._qualify_unit_only_numerical_aggregate(
            attempts=attempts,
            ordered_raw_inventory_sha256_values=inventories
            or [wrapper["attempt"]["raw_inventory_sha256"] for wrapper in attempts],
            paired_amr_fine_results=pairs or _pairs(attempts),
            restart_parity_binding=_restart_binding(attempts),
            reviewer_disposition=_pending(),
        )


class Q011Section54NumericalQualificationTests(unittest.TestCase):
    def test_raw_reducers_use_observed_committed_time_not_nominal_slot(self) -> None:
        observed_t500 = 500.05

        def product(path: str, nominal: float, observed: float) -> dict[str, object]:
            return {
                "path": path,
                "nominal_slot_time": nominal,
                "observed_committed_time": observed,
            }

        retained = {
            "500.0": {
                "prtcl_all": product("t500.prtcl.vtk", 500.0, observed_t500),
                **{
                    quantity: product(
                        f"t500.{quantity}.bin", 500.0, observed_t500
                    )
                    for quantity in qualifier.spatial.REQUIRED_MESH_QUANTITIES
                },
            },
            "1200.0": {
                "prtcl_all": product("t1200.prtcl.vtk", 1200.0, 1200.0),
            },
        }
        admitted = {
            "run_identity": _identity(
                0, qualifier.GRID_VARIANTS[0], qualifier.QUALIFYING_SEEDS[0]
            ),
            "immutable_tree": {"inventory_sha256": _sha256("1")},
            "retained_snapshot_products": retained,
            "snapshot_payloads": {
                "500.0": {
                    "nominal_slot_time": 500.0,
                    "observed_committed_time": observed_t500,
                },
                "1200.0": {
                    "nominal_slot_time": 1200.0,
                    "observed_committed_time": 1200.0,
                },
            },
        }
        snapshot = mock.Mock()
        member = mock.Mock()
        member.read_bytes.return_value = b"retained"
        snapshot.member_path.return_value = member
        decoded = mock.Mock(
            points=np.zeros((1, 3)),
            scalars={
                "cr_source": np.ones(1, dtype=np.int64),
                "birth_time": np.full(1, 45.0),
                "macro_weight": np.ones(1),
            },
            vectors={"vel": np.zeros((1, 3))},
        )

        def particle_record(**kwargs: object) -> dict[str, object]:
            return {"snapshot_time_omega0_inverse": kwargs["snapshot_time"]}

        with mock.patch.object(
            qualifier.admission, "read_particle_vtk", return_value=decoded
        ), mock.patch.object(
            qualifier.admission,
            "_parse_pvtk_execution_header",
            side_effect=({"time": observed_t500}, {"time": 1200.0}),
        ), mock.patch.object(
            qualifier.admission.output_primitives,
            "parse_athenak_binary_bytes",
            return_value=mock.sentinel.dataset,
        ), mock.patch.object(
            qualifier.particles,
            "reduce_particle_snapshot",
            side_effect=particle_record,
        ) as particle_reducer, mock.patch.object(
            qualifier.particles, "ideal_surface_x1", return_value=123.0
        ) as ideal_surface, mock.patch.object(
            qualifier.spatial,
            "reduce_t500_spatial_snapshot",
            return_value={"observed_committed_time": observed_t500},
        ) as spatial_reducer, mock.patch.object(
            qualifier,
            "_retained_source_checkpoint_lineage",
            return_value={"trusted": "lineage"},
        ):
            result = qualifier._reduce_retained_attempt(
                {"admission": admitted}, snapshot
            )

        self.assertEqual(
            [call.kwargs["snapshot_time"] for call in particle_reducer.call_args_list],
            [observed_t500, 1200.0],
        )
        self.assertEqual(
            result["particle_reductions"]["t500"]["nominal_slot_time"], 500.0
        )
        self.assertEqual(
            result["particle_reductions"]["t500"]["observed_committed_time"],
            observed_t500,
        )
        ideal_surface.assert_called_once_with(observed_t500)
        self.assertEqual(
            spatial_reducer.call_args.kwargs,
            {
                "nominal_slot_time": 500.0,
                "observed_committed_time": observed_t500,
                "x_ideal_c_over_omega_pi": 123.0,
            },
        )

    def test_trusted_pair_recompute_uses_bound_reducer_products(self) -> None:
        amr = {
            "spatial_reduction": _paired_spatial_inputs(),
            "particle_reductions": _particle_reductions(),
        }
        fine = copy.deepcopy(amr)
        residuals = qualifier._recompute_pair_residuals(amr, fine)
        self.assertEqual(
            [record["observable"] for record in residuals],
            [row[0] for row in qualifier.PAIRED_RESIDUAL_THRESHOLDS],
        )
        self.assertTrue(
            all(record["maximum_absolute_difference"] == 0.0 for record in residuals)
        )

        fine = {
            "spatial_reduction": _paired_spatial_inputs(observed_t500=500.1),
            "particle_reductions": _particle_reductions(observed_t500=500.1),
        }
        qualifier._recompute_pair_residuals(amr, fine)

        fine = {
            "spatial_reduction": _paired_spatial_inputs(observed_t500=500.125),
            "particle_reductions": _particle_reductions(observed_t500=500.125),
        }
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "observed-time separation exceeds 0.1",
        ):
            qualifier._recompute_pair_residuals(amr, fine)

        fine = copy.deepcopy(amr)
        fine["spatial_reduction"] = _paired_spatial_inputs(center_offset=0.5)
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "restricted profile areas drifted",
        ):
            qualifier._recompute_pair_residuals(amr, fine)

    def test_trusted_restart_extractor_builds_observation_from_bundle_members(self) -> None:
        source = _restart_source(
            "uninterrupted", branch_role="uninterrupted_baseline"
        )
        payloads = {}
        authority_payloads = {}
        for output in source["outputs_after_checkpoint"]:
            for member in output["members"].values():
                payloads[member["bundle_member"]] = b"retained"
                authority_payloads[member["structured_artifact_member"]] = b"retained"
        policy = qualifier.restart.load_preregistration()
        contract = policy["continuation_contract"]
        binding = {
            "checkpoint_nominal_slot_omega0_inverse": contract[
                "checkpoint_nominal_slot_omega0_inverse"
            ],
            "checkpoint_observed_committed_cycle": 5500,
            "checkpoint_observed_committed_time_omega0_inverse": 500.0,
            "restart_schema": 7,
            "startup_shock_ledger": {},
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse": contract[
                "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
            ],
            "comparison_tolerances_max_absolute_difference": contract[
                "comparison_tolerances_max_absolute_difference"
            ],
        }
        tree = mock.MagicMock()
        tree.__enter__.return_value = tree
        tree.__exit__.return_value = None
        with mock.patch.object(
            qualifier.restart, "bind_checkpoint_for_continuation", return_value=binding
        ) as checkpoint_binder, mock.patch.object(
            qualifier, "_restart_mesh_values", return_value=[1.0]
        ) as mesh_values, mock.patch.object(
            qualifier, "_restart_particle_values", return_value=([1], [2.0])
        ) as particle_values, mock.patch.object(
            qualifier, "_RegisteredRunStructuredArtifactTree", return_value=tree
        ), mock.patch.object(
            qualifier.structured_artifacts, "require_inventory_sha256"
        ), mock.patch.object(
            qualifier.structured_artifacts, "load_inventory", return_value={}
        ), mock.patch.object(
            qualifier.structured_artifacts,
            "read_inventory_bytes",
            side_effect=lambda _tree, _inventory, relative: authority_payloads[relative],
        ):
            observation = qualifier._extract_retained_restart_observation(
                source,
                payloads,
                checkpoint_payload=b"checkpoint",
                checkpoint_commit={
                    "cycle": 5500,
                    "nominal_slot_time": 500.0,
                    "observed_committed_time": 500.0,
                },
                execution_receipt={
                    "raw_output_root": source["raw_output_root"],
                    "artifact_dir": "/authority/uninterrupted",
                },
                artifact_root_descriptor=-1,
                label="uninterrupted",
            )
        self.assertEqual(observation["binding"], binding)
        self.assertEqual(len(observation["outputs_after_checkpoint"]), 7)
        first_source = source["outputs_after_checkpoint"][0]
        first_output = observation["outputs_after_checkpoint"][0]
        self.assertEqual(
            first_output["nominal_slot_omega0_inverse"],
            first_source["nominal_slot_omega0_inverse"],
        )
        self.assertEqual(
            first_output["observed_committed_time_omega0_inverse"],
            first_source["observed_committed_time_omega0_inverse"],
        )
        checkpoint_binder.assert_called_once()
        self.assertEqual(
            checkpoint_binder.call_args.kwargs[
                "checkpoint_observed_committed_time_omega0_inverse"
            ],
            500.0,
        )
        self.assertEqual(
            mesh_values.call_args_list[0].kwargs["observed_committed_time"],
            first_source["observed_committed_time_omega0_inverse"],
        )
        self.assertEqual(
            particle_values.call_args_list[0].kwargs["observed_committed_time"],
            first_source["observed_committed_time_omega0_inverse"],
        )
        self.assertEqual(
            set(first_output["fields"]),
            {
                "rho_bin",
                "bmag_bin",
                "prtcl_jx_bin",
                "j2_bin",
                "prtcl_all_pvtk_integer_payload",
                "prtcl_all_pvtk_float_payload",
            },
        )

    def test_restart_extractor_rejects_forged_bundle_copy_with_equal_shape(self) -> None:
        source = _restart_source(
            "continued", branch_role="checkpoint_restart_continuation"
        )
        first = source["outputs_after_checkpoint"][0]["members"]["rho_bin"]
        payloads = {first["bundle_member"]: b"forged!!"}
        tree = mock.MagicMock()
        tree.__enter__.return_value = tree
        tree.__exit__.return_value = None
        with mock.patch.object(
            qualifier.restart,
            "bind_checkpoint_for_continuation",
            return_value={
                "retained_output_nominal_slots_after_checkpoint_omega0_inverse": [
                    output["nominal_slot_omega0_inverse"]
                    for output in source["outputs_after_checkpoint"]
                ]
            },
        ), mock.patch.object(
            qualifier, "_RegisteredRunStructuredArtifactTree", return_value=tree
        ), mock.patch.object(
            qualifier.structured_artifacts, "require_inventory_sha256"
        ), mock.patch.object(
            qualifier.structured_artifacts, "load_inventory", return_value={}
        ), mock.patch.object(
            qualifier.structured_artifacts,
            "read_inventory_bytes",
            return_value=b"trusted!",
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "retained copy differs from structured raw authority",
        ):
            qualifier._extract_retained_restart_observation(
                source,
                payloads,
                checkpoint_payload=b"checkpoint",
                checkpoint_commit={
                    "cycle": 5500,
                    "nominal_slot_time": 500.0,
                    "observed_committed_time": 500.0,
                },
                execution_receipt={
                    "raw_output_root": source["raw_output_root"],
                    "artifact_dir": "/authority/continued",
                },
                artifact_root_descriptor=-1,
                label="continued",
            )

    def test_registered_run_structured_tree_allows_only_reconciler_receipt(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        artifact_dir = Path(temporary.name) / "runs/q011/submission"
        raw = artifact_dir / "raw"
        analysis = artifact_dir / "analysis"
        raw.mkdir(parents=True)
        analysis.mkdir()
        payload = b"trusted raw output"
        member = raw / "output.bin"
        member.write_bytes(payload)
        member.chmod(0o444)
        receipt = analysis / admission.REGISTERED_EXECUTION_RECEIPT_NAME
        receipt.write_bytes(b"registered execution")
        receipt.chmod(0o444)
        inventory_payload = qualifier.structured_artifacts.canonical_json_bytes(
            {
                "schema_version": 1,
                "files": [
                    {
                        "path": "raw/output.bin",
                        "size": len(payload),
                        "sha256": hashlib.sha256(payload).hexdigest(),
                    }
                ],
            }
        )
        inventory_path = artifact_dir / "artifact_inventory.json"
        inventory_path.write_bytes(inventory_payload)
        inventory_path.chmod(0o444)
        analysis.chmod(0o700)
        raw.chmod(0o555)
        artifact_dir.chmod(0o555)
        self.addCleanup(os.chmod, artifact_dir, 0o755)
        self.addCleanup(os.chmod, raw, 0o755)
        with qualifier._RegisteredRunStructuredArtifactTree(artifact_dir) as tree:
            inventory = qualifier.structured_artifacts.load_inventory(tree)
            self.assertEqual(
                qualifier.structured_artifacts.read_inventory_bytes(
                    tree, inventory, "raw/output.bin"
                ),
                payload,
            )
        unexpected = analysis / "unexpected.json"
        unexpected.write_bytes(b"{}")
        unexpected.chmod(0o444)
        with self.assertRaisesRegex(
            ValueError,
            "unexpected entry",
        ):
            with qualifier._RegisteredRunStructuredArtifactTree(artifact_dir):
                pass

    def test_complete_matrix_passes_numerically_but_remains_pending_review(self) -> None:
        attempts = _attempts()
        result = _qualify(attempts)
        self.assertEqual(result["attempt_count"], 24)
        self.assertTrue(result["numerical_gates_passed"])
        self.assertEqual(result["numerical_gate_status"], "passed")
        self.assertEqual(
            result["restart_parity"]["screen_scope"],
            "unit_only_synthetic_restart_comparator",
        )
        self.assertFalse(result["restart_parity"]["full_state_equivalence_claimed"])
        self.assertEqual(result["status"], "pending_external_review")
        self.assertFalse(result["final_claim_closure"])
        self.assertEqual(
            result["ordered_raw_inventory_sha256_values"],
            [wrapper["attempt"]["raw_inventory_sha256"] for wrapper in attempts],
        )

    def test_missing_duplicate_and_matrix_cell_drift_fail_closed(self) -> None:
        attempts = _attempts()
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "exactly 24"):
            _qualify(attempts[:-1], pairs=_pairs(attempts))

        duplicate = copy.deepcopy(attempts)
        duplicate[1] = copy.deepcopy(duplicate[0])
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "duplicate canonical"):
            _qualify(duplicate)

        drift = copy.deepcopy(attempts)
        drift[1]["attempt"]["run_identity"]["seed"] = qualifier.QUALIFYING_SEEDS[0]
        drift[1]["attempt"]["admission_result"]["admission"]["run_identity"]["seed"] = (
            qualifier.QUALIFYING_SEEDS[0]
        )
        drift[1]["attempt"]["admission_result"]["admission"]["amr_pairing"]["pair_key"]["seed"] = (
            qualifier.QUALIFYING_SEEDS[0]
        )
        drift[1] = qualifier._bind_unit_only_canonical_attempt(drift[1]["attempt"])
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "matrix order drifted"):
            _qualify(drift)

    def test_cross_seed_amr_fine_pair_fails_closed(self) -> None:
        attempts = _attempts()
        pairs = _pairs(attempts)
        pairs[0]["fine_uniform_attempt_sha256"] = pairs[1]["fine_uniform_attempt_sha256"]
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "matched seed"):
            _qualify(attempts, pairs=pairs)

    def test_paired_t500_observed_time_separation_fails_closed(self) -> None:
        attempts = _attempts()
        fine_index = 2 * len(qualifier.QUALIFYING_SEEDS)
        fine = attempts[fine_index]["attempt"]
        fine["particle_reductions"] = _particle_reductions(observed_t500=500.125)
        fine["spatial_reduction"] = _spatial_reduction(observed_t500=500.125)
        fine["source_checkpoint_lineage"]["observed_committed_time"] = 500.125
        attempts[fine_index] = qualifier._bind_unit_only_canonical_attempt(fine)
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "observed-time separation exceeds 0.1",
        ):
            _qualify(attempts)

    def test_attempt_and_pair_numerical_failures_aggregate_without_claim_closure(self) -> None:
        attempts = _attempts()
        attempts[0]["attempt"]["spatial_reduction"] = _spatial_reduction(passes=False)
        attempts[0] = qualifier._bind_unit_only_canonical_attempt(attempts[0]["attempt"])
        result = _qualify(attempts)
        self.assertFalse(result["numerical_gates_passed"])
        self.assertEqual(result["numerical_gate_status"], "failed")
        self.assertEqual(result["status"], "pending_external_review")

        attempts = _attempts()
        pairs = _pairs(attempts)
        pairs[0]["residuals"][0]["maximum_absolute_difference"] = 241.0
        result = _qualify(attempts, pairs=pairs)
        self.assertFalse(result["numerical_gates_passed"])
        self.assertFalse(result["paired_amr_fine_results"][0]["gates_passed"])

    def test_canonical_attempt_tamper_and_inventory_order_tamper_fail_closed(self) -> None:
        attempts = _attempts()
        attempts[0]["attempt"]["raw_inventory_sha256"] = _sha256("f")
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "canonical attempt SHA"):
            _qualify(attempts)

        attempts = _attempts()
        inventories = [
            wrapper["attempt"]["raw_inventory_sha256"] for wrapper in attempts
        ]
        inventories[0], inventories[1] = inventories[1], inventories[0]
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "inventory digest binding"):
            _qualify(attempts, inventories=inventories)

    def test_public_aggregate_rejects_fabricated_attempt_dictionaries(self) -> None:
        attempts = _attempts()
        self.assertFalse(hasattr(qualifier, "bind_canonical_attempt"))
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            r"attempts\[0\]: keys drifted",
        ):
            qualifier.qualify_numerical_aggregate(
                attempts=attempts,
                ordered_raw_inventory_sha256_values=[
                    wrapper["attempt"]["raw_inventory_sha256"]
                    for wrapper in attempts
                ],
                paired_amr_fine_results=_pairs(attempts),
                restart_parity_binding=_restart_binding(attempts),
                reviewer_disposition=_pending(),
            )

    def test_bind_retained_attempt_tree_accepts_admitted_frozen_campaign(self) -> None:
        def reduce_fixture(
            result: dict[str, object],
            _snapshot: object,
        ) -> dict[str, object]:
            admitted = result["admission"]
            attempt = copy.deepcopy(_attempts()[0]["attempt"])
            attempt["run_identity"] = copy.deepcopy(admitted["run_identity"])
            attempt["raw_inventory_sha256"] = admitted["immutable_tree"][
                "inventory_sha256"
            ]
            attempt["source_checkpoint_lineage"] = _source_checkpoint_lineage(
                attempt["run_identity"]
            )
            attempt["admission_result"] = copy.deepcopy(result)
            return attempt

        with campaign_fixture._frozen_fixture() as (
            authorized_root,
            tree,
            inventory_sha256,
        ):
            with _authoritative_pressure_pilot_consumer_fixture(
                authorized_root
            ) as consumer, mock.patch.object(
                qualifier,
                "_reduce_retained_attempt",
                side_effect=reduce_fixture,
            ) as reducer:
                binding = qualifier.bind_retained_attempt_tree(
                    tree,
                    inventory_sha256,
                    authorized_orion_root=authorized_root,
                )

        wrapper = qualifier._load_json_bytes(
            binding.wrapper_payload,
            label="direct retained-attempt binding fixture",
        )
        self.assertIsInstance(binding, qualifier._RetainedAttemptBinding)
        self.assertTrue(
            wrapper["attempt"]["admission_result"][
                "admitted_for_follow_on_numerical_qualification"
            ]
        )
        self.assertEqual(
            wrapper["attempt_sha256"],
            qualifier._canonical_sha256(wrapper["attempt"]),
        )
        consumer.assert_called_once()
        reducer.assert_called_once()

    def test_bind_retained_attempt_tree_rejects_review_packet_binding_drift(
        self,
    ) -> None:
        def drift_review_packet_binding(
            root: Path,
            manifest: dict[str, object],
        ) -> None:
            campaign_fixture._rewrite_attempt_json_artifact(
                root,
                manifest,
                "selected_pressure_receipt",
                lambda receipt: receipt[
                    "published_pressure_pilot_review_packet_receipt"
                ].__setitem__("sha256", "3" * 64),
            )

        observed: dict[str, object] = {}
        qualify_campaign = qualifier.admission.qualify_campaign

        def capture_admission(*args: object, **kwargs: object) -> dict[str, object]:
            result = qualify_campaign(*args, **kwargs)
            observed.update(copy.deepcopy(result))
            return result

        with campaign_fixture._frozen_fixture(drift_review_packet_binding) as (
            authorized_root,
            tree,
            inventory_sha256,
        ):
            with _authoritative_pressure_pilot_consumer_fixture(
                authorized_root
            ) as consumer, mock.patch.object(
                qualifier.admission,
                "qualify_campaign",
                side_effect=capture_admission,
            ), mock.patch.object(qualifier, "_reduce_retained_attempt") as reducer:
                with self.assertRaisesRegex(
                    qualifier.NumericalQualificationError,
                    "retained attempt was not admitted",
                ):
                    qualifier.bind_retained_attempt_tree(
                        tree,
                        inventory_sha256,
                        authorized_orion_root=authorized_root,
                    )

        self.assertEqual(
            observed["failure_reasons"][0]["code"],
            "selected_pressure_receipt_drift",
        )
        self.assertIn(
            "review-packet receipt binding drifted",
            observed["failure_reasons"][0]["message"],
        )
        consumer.assert_not_called()
        reducer.assert_not_called()

    def test_public_aggregate_loads_retained_descriptors_before_passing(self) -> None:
        attempts = _attempts()
        retained_attempts = _retained_attempt_bindings(attempts)
        retained_pairs = _retained_pair_bindings(attempts)
        retained_restart = _retained_restart_binding(attempts)
        restart_source = _restart_binding(attempts)["source_attempt_sha256"]
        with mock.patch.object(
            qualifier, "bind_retained_attempt_tree", side_effect=retained_attempts
        ), mock.patch.object(
            qualifier, "bind_retained_pair_result", side_effect=retained_pairs
        ), mock.patch.object(
            qualifier,
            "bind_retained_restart_parity",
            return_value=retained_restart,
        ), mock.patch.object(
            qualifier.restart,
            "compare_deterministic_continuation_parity",
            return_value=_parity_result(),
        ):
            result = qualifier.qualify_numerical_aggregate(
                attempts=_attempt_tree_descriptors(),
                ordered_raw_inventory_sha256_values=[
                    wrapper["attempt"]["raw_inventory_sha256"]
                    for wrapper in attempts
                ],
                paired_amr_fine_results=_recompute_descriptors(),
                restart_parity_binding={
                    "bundle_root": "/retained/restart",
                    "expected_inventory_sha256": _sha256("f"),
                    "source_attempt_sha256": restart_source,
                },
                reviewer_disposition=_pending(),
            )
        self.assertTrue(result["numerical_gates_passed"])
        self.assertEqual(result["numerical_gate_status"], "passed")
        self.assertEqual(
            result["restart_parity"]["screen_scope"],
            qualifier.RESTART_SCREEN_SCOPE,
        )
        self.assertFalse(result["restart_parity"]["full_state_equivalence_claimed"])

    def test_public_aggregate_rejects_fabricated_pair_dictionaries(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        with mock.patch.object(
            qualifier, "bind_retained_attempt_tree", side_effect=retained
        ):
            with self.assertRaisesRegex(
                qualifier.NumericalQualificationError,
                r"paired_amr_fine_results\[0\]: keys drifted",
            ):
                qualifier.qualify_numerical_aggregate(
                    attempts=_attempt_tree_descriptors(),
                    ordered_raw_inventory_sha256_values=[
                        wrapper["attempt"]["raw_inventory_sha256"]
                        for wrapper in attempts
                    ],
                    paired_amr_fine_results=_pairs(attempts),
                    restart_parity_binding=_restart_binding(attempts),
                    reviewer_disposition=_pending(),
                )

    def test_public_aggregate_rejects_fabricated_restart_dictionary(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        with mock.patch.object(
            qualifier, "bind_retained_attempt_tree", side_effect=retained
        ), mock.patch.object(
            qualifier,
            "bind_retained_pair_result",
            return_value=mock.sentinel.retained_pair,
        ):
            with self.assertRaisesRegex(
                qualifier.NumericalQualificationError,
                "restart_parity_binding: keys drifted",
            ):
                qualifier.qualify_numerical_aggregate(
                    attempts=_attempt_tree_descriptors(),
                    ordered_raw_inventory_sha256_values=[
                        wrapper["attempt"]["raw_inventory_sha256"]
                        for wrapper in attempts
                    ],
                    paired_amr_fine_results=_recompute_descriptors(),
                    restart_parity_binding=_restart_binding(attempts),
                    reviewer_disposition=_pending(),
                )

    def test_retained_recompute_receipts_bind_raw_attempt_inventories(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        amr = retained[len(qualifier.QUALIFYING_SEEDS)]
        fine = retained[2 * len(qualifier.QUALIFYING_SEEDS)]
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)

        pair_record = {
            "schema_version": qualifier.SCHEMA_VERSION,
            "record_type": qualifier.PAIR_RECOMPUTE_RECORD_TYPE,
            "seed": qualifier.QUALIFYING_SEEDS[0],
            "amr_attempt_sha256": source["attempt_sha256"],
            "amr_raw_inventory_sha256": source["attempt"]["raw_inventory_sha256"],
            "fine_uniform_attempt_sha256": counterpart["attempt_sha256"],
            "fine_uniform_raw_inventory_sha256": counterpart["attempt"][
                "raw_inventory_sha256"
            ],
            "residuals": _residuals(),
        }
        pair_root, pair_inventory = _publish_recompute_bundle(
            temporary,
            name="pair",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=pair_record,
        )
        self.addCleanup(_make_tree_writable, Path(pair_root))
        with mock.patch.object(
            qualifier, "_recompute_pair_residuals", return_value=_residuals()
        ):
            self.assertIsInstance(
                qualifier.bind_retained_pair_result(
                    pair_root,
                    pair_inventory,
                    amr_attempt=amr,
                    fine_uniform_attempt=fine,
                    authorized_orion_root=Path(temporary.name),
                ),
                qualifier._RetainedPairResultBinding,
            )

        fabricated = copy.deepcopy(pair_record)
        fabricated["fine_uniform_raw_inventory_sha256"] = _sha256("f")
        bad_root, bad_inventory = _publish_recompute_bundle(
            temporary,
            name="fabricated-pair",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=fabricated,
        )
        self.addCleanup(_make_tree_writable, Path(bad_root))
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "retained pair recompute raw attempt binding drifted",
        ):
            qualifier.bind_retained_pair_result(
                bad_root,
                bad_inventory,
                amr_attempt=amr,
                fine_uniform_attempt=fine,
                authorized_orion_root=Path(temporary.name),
            )

        restart_record, restart_artifacts, planned = _retained_restart_fixture(source)
        restart_root, restart_inventory = _publish_recompute_bundle(
            temporary,
            name="restart",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.RESTART_RECOMPUTE_ANALYZERS,
            member=qualifier.RESTART_RECOMPUTE_MEMBER,
            record=restart_record,
            extra_artifacts=restart_artifacts,
        )
        self.addCleanup(_make_tree_writable, Path(restart_root))
        with mock.patch.object(
            qualifier, "_validate_planned_restart_carrier", return_value=planned
        ), mock.patch.object(
            qualifier,
            "_validated_authoritative_restart_execution_receipt",
            side_effect=_authoritative_receipt_context_side_effect,
        ), mock.patch.object(
            qualifier,
            "_extract_retained_restart_observation",
            side_effect=({"trusted": "uninterrupted"}, {"trusted": "continued"}),
        ):
            self.assertIsInstance(
                qualifier.bind_retained_restart_parity(
                    restart_root,
                    restart_inventory,
                    source_attempt=amr,
                    authorized_orion_root=Path(temporary.name),
                ),
                qualifier._RetainedRestartParityBinding,
            )

    def test_retained_recompute_bundles_must_remain_below_authorized_orion_root(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        outside = tempfile.TemporaryDirectory()
        authorized = tempfile.TemporaryDirectory()
        self.addCleanup(outside.cleanup)
        self.addCleanup(authorized.cleanup)
        pair_record = {
            "schema_version": qualifier.SCHEMA_VERSION,
            "record_type": qualifier.PAIR_RECOMPUTE_RECORD_TYPE,
            "seed": qualifier.QUALIFYING_SEEDS[0],
            "amr_attempt_sha256": source["attempt_sha256"],
            "amr_raw_inventory_sha256": source["attempt"]["raw_inventory_sha256"],
            "fine_uniform_attempt_sha256": counterpart["attempt_sha256"],
            "fine_uniform_raw_inventory_sha256": counterpart["attempt"]["raw_inventory_sha256"],
            "residuals": _residuals(),
        }
        pair_root, pair_inventory = _publish_recompute_bundle(
            outside,
            name="off-root-pair",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=pair_record,
        )
        self.addCleanup(_make_tree_writable, Path(pair_root))
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "retained tree must remain below",
        ):
            qualifier.bind_retained_pair_result(
                pair_root,
                pair_inventory,
                amr_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                fine_uniform_attempt=retained[2 * len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(authorized.name),
            )

    def test_fabricated_sealed_pair_residual_values_fail_closed(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        fabricated_residuals = _residuals()
        fabricated_residuals[0]["maximum_absolute_difference"] = 1.0
        pair_root, pair_inventory = _publish_recompute_bundle(
            temporary,
            name="fabricated-residual-pair",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record={
                "schema_version": qualifier.SCHEMA_VERSION,
                "record_type": qualifier.PAIR_RECOMPUTE_RECORD_TYPE,
                "seed": qualifier.QUALIFYING_SEEDS[0],
                "amr_attempt_sha256": source["attempt_sha256"],
                "amr_raw_inventory_sha256": source["attempt"]["raw_inventory_sha256"],
                "fine_uniform_attempt_sha256": counterpart["attempt_sha256"],
                "fine_uniform_raw_inventory_sha256": counterpart["attempt"]["raw_inventory_sha256"],
                "residuals": fabricated_residuals,
            },
        )
        self.addCleanup(_make_tree_writable, Path(pair_root))
        with mock.patch.object(
            qualifier, "_recompute_pair_residuals", return_value=_residuals()
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "differ from trusted raw-attempt recompute",
        ):
            qualifier.bind_retained_pair_result(
                pair_root,
                pair_inventory,
                amr_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                fine_uniform_attempt=retained[2 * len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_recompute_root_descriptor_binding_race_fails_closed(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        pair_root, pair_inventory = _publish_recompute_bundle(
            temporary,
            name="descriptor-root-binding-race",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=_pair_record(source, counterpart),
        )
        self.addCleanup(_make_tree_writable, Path(pair_root))
        with mock.patch.object(
            qualifier.immutable_orion_tree,
            "_require_root_binding",
            side_effect=qualifier.NumericalQualificationError(
                "retained tree root binding changed during anchored operation"
            ),
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "root binding changed",
        ):
            qualifier.bind_retained_pair_result(
                pair_root,
                pair_inventory,
                amr_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                fine_uniform_attempt=retained[2 * len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_recompute_bundle_empty_directory_fails_tree_closure(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        pair_root, pair_inventory = _publish_recompute_bundle(
            temporary,
            name="empty-directory-closure",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=_pair_record(source, counterpart),
        )
        root = Path(pair_root)
        self.addCleanup(_make_tree_writable, root)
        root.chmod(0o755)
        (root / "undeclared-empty-directory").mkdir(mode=0o555)
        root.chmod(0o555)
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "tree closure drifted",
        ):
            qualifier.bind_retained_pair_result(
                pair_root,
                pair_inventory,
                amr_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                fine_uniform_attempt=retained[2 * len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_recompute_bundle_same_owner_hide_restore_fails_stable_closure_scan(
        self,
    ) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        pair_root, pair_inventory = _publish_recompute_bundle(
            temporary,
            name="hide-restore-closure",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.PAIR_RECOMPUTE_ANALYZERS,
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=_pair_record(source, counterpart),
        )
        root = Path(pair_root)
        self.addCleanup(_make_tree_writable, root)
        real_scan = qualifier.immutable_orion_tree._scan_anchored_tree
        calls = 0

        def hide_restore(descriptor: int, **kwargs: object) -> object:
            nonlocal calls
            calls += 1
            if calls != 2:
                return real_scan(descriptor, **kwargs)
            root.chmod(0o755)
            hidden = root / "same-owner-hide-restore"
            hidden.mkdir(mode=0o555)
            root.chmod(0o555)
            try:
                return real_scan(descriptor, **kwargs)
            finally:
                root.chmod(0o755)
                hidden.rmdir()
                root.chmod(0o555)

        with mock.patch.object(
            qualifier.immutable_orion_tree,
            "_scan_anchored_tree",
            side_effect=hide_restore,
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "changed during recompute bundle closure census",
        ):
            qualifier.bind_retained_pair_result(
                pair_root,
                pair_inventory,
                amr_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                fine_uniform_attempt=retained[2 * len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_recompute_manifest_rejects_extra_analyzer_binding(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        counterpart = attempts[2 * len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        pair_root, pair_inventory = _publish_recompute_bundle(
            temporary,
            name="extra-analyzer-binding",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=(
                *qualifier.PAIR_RECOMPUTE_ANALYZERS,
                qualifier.RESTART_RECOMPUTE_ANALYZER,
            ),
            member=qualifier.PAIR_RECOMPUTE_MEMBER,
            record=_pair_record(source, counterpart),
        )
        self.addCleanup(_make_tree_writable, Path(pair_root))
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "analyzer binding set drifted",
        ):
            qualifier.bind_retained_pair_result(
                pair_root,
                pair_inventory,
                amr_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                fine_uniform_attempt=retained[2 * len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_recompute_analyzer_closure_lists_direct_semantic_dependencies(self) -> None:
        self.assertEqual(
            qualifier.PRESSURE_REVIEW_PACKET_VERIFIER_RECOMPUTE_ANALYZER,
            "tst/publication/frontier_control_plane/q011_pressure_review_packet_verifier.py",
        )
        self.assertTrue(
            {
                qualifier.ARTIFACTS_RECOMPUTE_ANALYZER,
                qualifier.IMMUTABLE_TREE_RECOMPUTE_ANALYZER,
                qualifier.PARTICLES_RECOMPUTE_ANALYZER,
                qualifier.SPATIAL_RECOMPUTE_ANALYZER,
                qualifier.OUTPUT_RECOMPUTE_ANALYZER,
                qualifier.MODEL_RECOMPUTE_ANALYZER,
            }
            <= set(qualifier.PAIR_RECOMPUTE_ANALYZERS)
        )
        self.assertTrue(
            {
                qualifier.CAMPAIGN_RECOMPUTE_ANALYZER,
                qualifier.PLANNER_RECOMPUTE_ANALYZER,
                qualifier.MODEL_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_SELECTION_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_HISTORICAL_CONSUMER_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_REVIEW_PACKET_VERIFIER_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_PILOT_EXECUTION_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_PILOT_PUBLISHER_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_PILOT_RECOMPUTE_ANALYZER,
                qualifier.PRESSURE_PILOT_CASE_RECOMPUTE_ANALYZER,
                qualifier.STRUCTURED_ARTIFACTS_RECOMPUTE_ANALYZER,
                qualifier.CONTROL_PLANE_COMMON_RECOMPUTE_ANALYZER,
                qualifier.CONTROL_PLANE_LEDGER_RECOMPUTE_ANALYZER,
                qualifier.CONTROL_PLANE_ATTESTATION_RECOMPUTE_ANALYZER,
            }
            <= set(qualifier.RESTART_RECOMPUTE_ANALYZERS)
        )

    def test_restart_branches_must_not_alias_retained_outputs(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        record, retained_artifacts, planned = _retained_restart_fixture(source)
        record["continued"]["outputs_after_checkpoint"][0]["members"]["rho_bin"] = (
            record["uninterrupted"]["outputs_after_checkpoint"][0]["members"]["rho_bin"]
        )
        root, inventory = _publish_recompute_bundle(
            temporary,
            name="restart-branch-alias",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.RESTART_RECOMPUTE_ANALYZERS,
            member=qualifier.RESTART_RECOMPUTE_MEMBER,
            record=record,
            extra_artifacts=retained_artifacts,
        )
        self.addCleanup(_make_tree_writable, Path(root))
        with mock.patch.object(
            qualifier, "_validate_planned_restart_carrier", return_value=planned
        ), mock.patch.object(
            qualifier,
            "_validated_authoritative_restart_execution_receipt",
            side_effect=_authoritative_receipt_context_side_effect,
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "branches or metadata members alias",
        ):
            qualifier.bind_retained_restart_parity(
                root,
                inventory,
                source_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_restart_receipt_must_equal_fixed_path_authority_and_bind_mirrored_ledger(
        self,
    ) -> None:
        source = _attempts()[len(qualifier.QUALIFYING_SEEDS)]
        planned = _planned_execution_inventory()
        branch = _restart_source(
            "continued", branch_role="checkpoint_restart_continuation"
        )
        authority = (
            Path("/authorized/runs/restart/analysis")
            / admission.REGISTERED_EXECUTION_RECEIPT_NAME
        )
        branch["execution_receipt_authority_path"] = str(authority)
        payload = b"authoritative registered execution receipt"
        payloads = {branch["execution_receipt_member"]: payload}
        receipt = {
            **_authoritative_receipt_projection("02"),
            "artifact_dir": str(authority.parent.parent),
        }
        source_payload = copy.deepcopy(source["attempt"])
        source_payload["admission_result"]["admission"].update(
            {
                "artifact_bindings": {
                    "executable": {"sha256": planned["executable_sha256"]},
                    "deck": {"sha256": planned["paper_deck_sha256"]},
                },
                "retained_attempt_semantics": {
                    "registered_execution_receipt": {},
                    "registered_execution_ledger_binding": {},
                    "attempt_contract": {"argv": ["baseline"]},
                },
            }
        )
        ledger = {"submission_id": receipt["submission_id"]}
        with mock.patch.object(
            qualifier,
            "_pinned_authoritative_execution_receipt",
            return_value=_pinned_receipt_fixture(authority, payload),
        ), mock.patch.object(
            admission,
            "_validate_registered_execution_receipt",
            return_value=receipt,
        ) as validate_receipt, mock.patch.object(
            admission,
            "_validated_registered_execution_receipt_ledger_snapshot",
            return_value=_ledger_snapshot_fixture(ledger),
        ) as validate_ledger:
            observed = qualifier._validate_authoritative_restart_execution_receipt(
                branch,
                payloads,
                branch_role="checkpoint_restart_continuation",
                source={
                    "identity": source["attempt"]["run_identity"],
                },
                source_attempt_payload=source_payload,
                source_checkpoint_lineage=source["attempt"][
                    "source_checkpoint_lineage"
                ],
                planned=planned,
                authorized_orion_root=Path("/authorized"),
                label="continued",
            )
        self.assertEqual(observed, receipt)
        validate_receipt.assert_called_once()
        validate_ledger.assert_called_once()

        with mock.patch.object(
            qualifier,
            "_pinned_authoritative_execution_receipt",
            return_value=_pinned_receipt_fixture(
                authority, b"different authoritative bytes"
            ),
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "differs from immutable fixed-path authority",
        ):
            qualifier._validate_authoritative_restart_execution_receipt(
                branch,
                payloads,
                branch_role="checkpoint_restart_continuation",
                source={"identity": source["attempt"]["run_identity"]},
                source_attempt_payload=source_payload,
                source_checkpoint_lineage=source["attempt"][
                    "source_checkpoint_lineage"
                ],
                planned=planned,
                authorized_orion_root=Path("/authorized"),
                label="continued",
            )

        with mock.patch.object(
            qualifier,
            "_pinned_authoritative_execution_receipt",
            return_value=_pinned_receipt_fixture(authority, payload),
        ), mock.patch.object(
            admission,
            "_validate_registered_execution_receipt",
            return_value=receipt,
        ), mock.patch.object(
            admission,
            "_validated_registered_execution_receipt_ledger_snapshot",
            return_value=_ledger_snapshot_fixture(),
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "authoritative registered execution receipt failed validation",
        ):
            qualifier._validate_authoritative_restart_execution_receipt(
                branch,
                payloads,
                branch_role="checkpoint_restart_continuation",
                source={"identity": source["attempt"]["run_identity"]},
                source_attempt_payload=source_payload,
                source_checkpoint_lineage=source["attempt"][
                    "source_checkpoint_lineage"
                ],
                planned=planned,
                authorized_orion_root=Path("/authorized"),
                label="continued",
            )

    def test_restart_receipt_ancestry_pin_rejects_transient_hide_restore(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        submission = root / "runs/q011/submission"
        analysis = submission / "analysis"
        analysis.mkdir(parents=True)
        receipt = analysis / admission.REGISTERED_EXECUTION_RECEIPT_NAME
        receipt.write_bytes(b"retained receipt")
        receipt.chmod(0o444)
        hidden = submission.with_name("submission-hidden")
        with self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "ancestry changed",
        ):
            with qualifier._pinned_authoritative_execution_receipt(
                str(receipt),
                authorized_orion_root=root,
                label="transient receipt",
            ):
                time.sleep(0.02)
                submission.rename(hidden)
                hidden.rename(submission)

    def test_restart_source_checkpoint_digest_drift_fails_closed(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        record, retained_artifacts, planned = _retained_restart_fixture(source)
        record["source_checkpoint_sha256"] = _sha256("f")
        root, inventory = _publish_recompute_bundle(
            temporary,
            name="restart-checkpoint-digest-drift",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.RESTART_RECOMPUTE_ANALYZERS,
            member=qualifier.RESTART_RECOMPUTE_MEMBER,
            record=record,
            extra_artifacts=retained_artifacts,
        )
        self.addCleanup(_make_tree_writable, Path(root))
        with mock.patch.object(
            qualifier, "_validate_planned_restart_carrier", return_value=planned
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "source checkpoint lineage drifted",
        ):
            qualifier.bind_retained_restart_parity(
                root,
                inventory,
                source_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_restart_baseline_receipt_authority_path_must_match_admitted_source(
        self,
    ) -> None:
        source = _attempts()[len(qualifier.QUALIFYING_SEEDS)]
        planned = _planned_execution_inventory()
        branch = _restart_source(
            "uninterrupted", branch_role="uninterrupted_baseline"
        )
        wrong_authority = (
            Path("/authorized/runs/wrong/analysis")
            / admission.REGISTERED_EXECUTION_RECEIPT_NAME
        )
        branch["execution_receipt_authority_path"] = str(wrong_authority)
        payload = b"authoritative registered execution receipt"
        source_payload = copy.deepcopy(source["attempt"])
        source_payload["admission_result"]["admission"].update(
            {
                "artifact_bindings": {
                    "executable": {"sha256": planned["executable_sha256"]},
                    "deck": {"sha256": planned["paper_deck_sha256"]},
                },
                "retained_attempt_semantics": {
                    "registered_execution_receipt": {
                        "raw_output_root": "/retained/baseline/raw",
                        "artifact_dir": "/authorized/runs/baseline",
                    },
                    "registered_execution_ledger_binding": {},
                    "attempt_contract": {"argv": ["baseline"]},
                },
            }
        )
        with mock.patch.object(
            qualifier,
            "_pinned_authoritative_execution_receipt",
            return_value=_pinned_receipt_fixture(wrong_authority, payload),
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "differs from admitted source attempt",
        ):
            qualifier._validate_authoritative_restart_execution_receipt(
                branch,
                {branch["execution_receipt_member"]: payload},
                branch_role="uninterrupted_baseline",
                source={"identity": source["attempt"]["run_identity"]},
                source_attempt_payload=source_payload,
                source_checkpoint_lineage=source["attempt"][
                    "source_checkpoint_lineage"
                ],
                planned=planned,
                authorized_orion_root=Path("/authorized"),
                label="uninterrupted",
            )

    def test_restart_source_checkpoint_lineage_substitution_fails_closed(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        record, retained_artifacts, planned = _retained_restart_fixture(source)
        record["source_checkpoint_lineage"]["restart_member_path"] = (
            "rst/self-authored.rst"
        )
        root, inventory = _publish_recompute_bundle(
            temporary,
            name="restart-checkpoint-lineage-substitution",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.RESTART_RECOMPUTE_ANALYZERS,
            member=qualifier.RESTART_RECOMPUTE_MEMBER,
            record=record,
            extra_artifacts=retained_artifacts,
        )
        self.addCleanup(_make_tree_writable, Path(root))
        with mock.patch.object(
            qualifier, "_validate_planned_restart_carrier", return_value=planned
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "source checkpoint lineage drifted",
        ):
            qualifier.bind_retained_restart_parity(
                root,
                inventory,
                source_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )

    def test_restart_planned_carrier_is_reconstructed_from_frozen_planner_tree(
        self,
    ) -> None:
        source = _attempts()[len(qualifier.QUALIFYING_SEEDS)]
        record, payloads, tree_payloads, parsed_source, source_payload = (
            _planner_restart_fixture(source)
        )
        with mock.patch.object(
            admission,
            "_validate_planner_materialization_receipt",
            return_value={"_tree_member_payloads": tree_payloads},
        ):
            planned = qualifier._validate_planned_restart_carrier(
                record,
                payloads,
                source=parsed_source,
                source_attempt_payload=source_payload,
                authorized_orion_root=Path(
                    "/lustre/orion/ast207/proj-shared/dfielding/PIC"
                ),
            )
        self.assertEqual(
            planned["planned_restart_carrier_sha256"],
            record["planned_restart_carrier"]["sha256"],
        )
        self.assertEqual(
            planned["control_plane_version"],
            _planned_execution_inventory()["control_plane_version"],
        )

    def test_restart_self_authored_carrier_rebinding_fails_source_plan_derivation(
        self,
    ) -> None:
        source = _attempts()[len(qualifier.QUALIFYING_SEEDS)]
        record, payloads, tree_payloads, parsed_source, source_payload = (
            _planner_restart_fixture(source)
        )
        carrier_path = record["planned_restart_carrier"]["path"]
        carrier = qualifier._load_json_bytes(
            tree_payloads[carrier_path], label=carrier_path
        )
        carrier["authorized_orion_attempt_root"] = "/self-authored/restart"
        carrier_payload = admission.campaign_planner._json_bytes(carrier)
        rebound = {
            "path": carrier_path,
            "sha256": hashlib.sha256(carrier_payload).hexdigest(),
        }
        record["planned_restart_carrier"] = rebound
        tree_payloads[carrier_path] = carrier_payload
        plan = qualifier._load_json_bytes(
            payloads["planner/campaign_plan.json"], label="campaign plan"
        )
        plan["restart_continuation_carrier"] = rebound
        payloads["planner/campaign_plan.json"] = artifacts.canonical_json_bytes(plan)
        with mock.patch.object(
            admission,
            "_validate_planner_materialization_receipt",
            return_value={"_tree_member_payloads": tree_payloads},
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "planner bytes differ from admitted source attempt",
        ):
            qualifier._validate_planned_restart_carrier(
                record,
                payloads,
                source=parsed_source,
                source_attempt_payload=source_payload,
                authorized_orion_root=Path(
                    "/lustre/orion/ast207/proj-shared/dfielding/PIC"
                ),
            )

    def test_fabricated_sealed_restart_observations_fail_closed(self) -> None:
        attempts = _attempts()
        retained = _retained_attempt_bindings(attempts)
        source = attempts[len(qualifier.QUALIFYING_SEEDS)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        record, restart_artifacts, planned = _retained_restart_fixture(source)
        record["uninterrupted"] = {"synthetic": "caller-authored-values"}
        record["continued"] = {"synthetic": "caller-authored-values"}
        restart_root, restart_inventory = _publish_recompute_bundle(
            temporary,
            name="fabricated-restart-observations",
            raw_inventory_sha256=source["attempt"]["raw_inventory_sha256"],
            analyzer_paths=qualifier.RESTART_RECOMPUTE_ANALYZERS,
            member=qualifier.RESTART_RECOMPUTE_MEMBER,
            record=record,
            extra_artifacts=restart_artifacts,
        )
        self.addCleanup(_make_tree_writable, Path(restart_root))
        with mock.patch.object(
            qualifier, "_validate_planned_restart_carrier", return_value=planned
        ), self.assertRaisesRegex(
            qualifier.NumericalQualificationError,
            "uninterrupted: keys drifted",
        ):
            qualifier.bind_retained_restart_parity(
                restart_root,
                restart_inventory,
                source_attempt=retained[len(qualifier.QUALIFYING_SEEDS)],
                authorized_orion_root=Path(temporary.name),
            )


if __name__ == "__main__":
    unittest.main()
