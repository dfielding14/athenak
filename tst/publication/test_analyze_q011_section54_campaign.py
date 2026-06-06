#!/usr/bin/env python3
"""Focused adversarial tests for the Q-011 Section 5.4 campaign admission gate."""

from __future__ import annotations

from contextlib import contextmanager
import copy
import hashlib
import io
import json
from pathlib import Path
import stat
import struct
import tarfile
import tempfile
from typing import Any, Callable, Iterator
import unittest
from unittest import mock

from tst.publication import analyze_q011_section54_campaign as campaign
from tst.publication import immutable_orion_tree
from tst.publication.frontier_control_plane import control_plane_common


_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_TIMES = tuple(float(value) for value in range(0, 1300, 100))
_RECEIPT = {
    "schema_version": 1,
    "artifact_role": campaign.ARTIFACT_ROLE,
    "qualification_effect": "retained_qualifying_campaign_attempt",
    "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_PLANNER_RECEIPT = {
    "schema_version": 1,
    "artifact_role": campaign._CAMPAIGN_PLAN_ARTIFACT_ROLE,
    "qualification_effect": campaign._CAMPAIGN_PLAN_QUALIFICATION_EFFECT,
    "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_FREEZE_ID = "03a7bd9a-7d4c-4e37-a12b-46de3817eff2"


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fnv1a64(payload: bytes) -> int:
    digest = 14695981039346656037
    for value in payload:
        digest ^= value
        digest = (digest * 1099511628211) & ((1 << 64) - 1)
    return digest


def _make_writable_tree(tree: Path) -> None:
    for path in [tree, *tree.rglob("*")]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | stat.S_IWUSR)


def _elf_executable() -> bytes:
    return b"\x7fELF" + bytes((2, 1, 1)) + b"\0" * 57


def _mesh_bin(kind: str, time: float, cycle: int) -> bytes:
    field = "dens" if kind == "rho" else kind
    parameter_header = (
        "<mesh>\n"
        "nx1=1\nnx2=1\nnx3=1\nnghost=0\n"
        "x1min=0.0\nx1max=1.0\n"
        "x2min=0.0\nx2max=1.0\n"
        "x3min=0.0\nx3max=1.0\n"
        "<meshblock>\n"
        "nx1=1\nnx2=1\nnx3=1\n"
    ).encode("ascii")
    header = (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={format(time, '.15g')}\n".encode("ascii")
        + f"  cycle={cycle}\n".encode("ascii")
        + b"  size of location=8\n"
        + b"  size of variable=4\n"
        + b"  number of variables=1\n"
        + f"  variables:  {field}  \n".encode("ascii")
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header
    )
    block = (
        struct.pack("<10i", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        + struct.pack("<6d", 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        + struct.pack("<f", 1.0)
    )
    return header + block


def _particle_vtk(
    time: float,
    cycle: int,
    *,
    ptags: tuple[int, int] = (10, 11),
    nranks: int = 4,
) -> bytes:
    points = ((100.0, 1.0, 0.0), (200.0, 2.0, 0.0))
    integer_scalars = {
        "gid": (0, 1),
        "ptag": ptags,
        "species": (0, 0),
        "cr_source": (1, 1),
    }
    real_scalars = {
        "macro_weight": (1.0, 2.0),
        "birth_time": (45.0, 46.0),
        "deltaf_f0": (1.0, 1.0),
        "deltaf_weight": (0.0, 0.0),
    }
    velocities = ((3.0, 4.0, 0.0), (6.0, 8.0, 0.0))
    payload = bytearray(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {time}  nranks= {nranks}  "
            f"cycle={cycle}  variables=prtcl_all\n"
            "BINARY\n"
            "DATASET UNSTRUCTURED_GRID\n"
            "\n"
            f"POINTS {len(points)} float\n"
        ).encode("ascii")
    )
    payload.extend(struct.pack(">6f", *(value for point in points for value in point)))
    payload.extend(f"\n\nPOINT_DATA {len(points)}\n".encode("ascii"))
    for name, values in integer_scalars.items():
        payload.extend(f"\nSCALARS {name} int\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(struct.pack(">2i", *values))
    for name, values in real_scalars.items():
        payload.extend(f"\nSCALARS {name} float\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(struct.pack(">2f", *values))
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(
        struct.pack(">6f", *(value for velocity in velocities for value in velocity))
    )
    return bytes(payload)


def _restart_payload(cycle: int) -> bytes:
    return (
        b"<job>\nbasename=q011\n<time>\ncycle="
        + str(cycle).encode("ascii")
        + b"\n<par_end>\n"
        + struct.pack("<4I", cycle, 1, 2, 3)
    )


def _restart_marker(payload: bytes) -> bytes:
    return (
        "ATHENAK_RESTART_COMPLETE_V1\n"
        f"size={len(payload)}\n"
        f"fnv1a64={_fnv1a64(payload):016x}\n"
    ).encode("ascii")


def _runtime_identity_line(**replacements: str) -> str:
    values = dict(campaign._EXPECTED_SECTION54_RUNTIME_PROJECTION)
    values.update(replacements)
    return "PIC runtime model: " + " ".join(
        f"{name}={value}" for name, value in values.items()
    )


def _stdout_telemetry(
    *,
    omit: frozenset[str] = frozenset(),
    include_runtime_identity: bool = True,
    runtime_replacements: dict[str, str] | None = None,
) -> bytes:
    lines = ["AthenaK retained stdout fixture"]
    if include_runtime_identity:
        lines.append(_runtime_identity_line(**(runtime_replacements or {})))
    for name in sorted(campaign._Q017_REQUIRED_NAMES - omit):
        value = 2.0 if name == "schema_version" else 1.0
        lines.append(f"q017.telemetry.{name}={value}")
    lines.extend(
        (
            "Terminating on time limit",
            "time=1200 cycle=120",
            "tlim=1200 nlim=-1",
        )
    )
    return ("\n".join(lines) + "\n").encode("ascii")


def _write_file(root: Path, relative: str, payload: bytes) -> dict[str, str]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {"path": relative, "sha256": _sha256(payload)}


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _binding(path: str, sha256: str) -> dict[str, str]:
    return {"path": path, "sha256": sha256}


def _inventory_payload(members: dict[str, bytes]) -> bytes:
    return "".join(
        f"{_sha256(payload)}  {relative}\n"
        for relative, payload in sorted(members.items())
    ).encode("utf-8")


def _published_pressure_receipt(
    root: Path, descriptors: list[dict[str, Any]]
) -> dict[str, str]:
    publication = root.parent / "publication"
    publication.mkdir(exist_ok=True)
    path = publication / "q011_section54_pressure_pilot_publication_receipt.json"
    payload = _json_bytes(
        {
            "aggregate_bundle": {
                "path": str(publication / "q011_section54_pressure_pilot_bundle"),
                "manifest_sha256": "b" * 64,
            },
            "aggregate_analysis": {
                "path": str(publication / "q011_section54_pressure_pilot_analysis.json"),
                "sha256": "c" * 64,
            },
            "raw_cases": [
                {
                    "case_id": descriptor["case_id"],
                    "artifact_dir": f"/fixture/{descriptor['case_id']}",
                    "descriptor_path": f"/fixture/{descriptor['case_id']}.json",
                    "descriptor_sha256": descriptor["descriptor_sha256"],
                    "artifact_inventory_sha256": f"{index + 100:064x}",
                    "runtime_artifacts": [],
                }
                for index, descriptor in enumerate(descriptors)
            ],
        }
    )
    path.write_bytes(payload)
    path.chmod(0o444)
    return {"path": str(path), "sha256": _sha256(payload)}


def _planner_materialization_receipt(
    root: Path,
    campaign_plan: dict[str, Any],
    source_payloads: dict[str, bytes],
    helper_closure_payload: bytes,
) -> tuple[dict[str, str], dict[str, str]]:
    planner = campaign.campaign_planner
    plan_id = campaign_plan["plan_id"]
    parent = root.parent / "plans"
    parent.mkdir(exist_ok=True)
    planner_root = parent / f"q011-section54-qualifying-campaign-plan-{plan_id}"
    planner_root.mkdir()
    materialized_members: dict[str, bytes] = {}

    def write(relative: str, payload: bytes) -> None:
        path = planner_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        materialized_members[relative] = payload

    for name, binding in campaign_plan["source_bindings"].items():
        write(binding["path"], source_payloads[name])

    campaign_root = Path(campaign_plan["authorized_orion_campaign_root"])
    selected_pressure = campaign_plan["selected_pressure"]["selected_case"][
        "problem_ps_p0"
    ]
    candidate = campaign_plan["candidate_binding"]
    source_bindings = campaign_plan["source_bindings"]
    contract_bindings = []
    descriptor_bindings = []
    descriptors = []
    index = 0
    for variant_id in campaign_plan["campaign_matrix"]["grid_variants"]:
        variant = planner._variant_binding(variant_id)
        for seed in campaign_plan["campaign_matrix"]["qualifying_seeds"]:
            index += 1
            attempt_id = planner._attempt_id(index, variant.variant, seed)
            artifact_root = campaign_root / "baseline" / attempt_id
            contract_path = f"launch_contracts/baseline/{attempt_id}.json"
            contract = planner._baseline_launch_contract(
                attempt_id=attempt_id,
                variant=variant,
                seed=seed,
                selected_ps_p0=selected_pressure,
                candidate=candidate,
                paper_deck_binding=source_bindings["paper_deck"],
                artifact_root=artifact_root,
            )
            contract_payload = _json_bytes(contract)
            write(contract_path, contract_payload)
            contract_binding = _binding(contract_path, _sha256(contract_payload))
            contract_bindings.append(contract_binding)
            descriptor = planner._attempt_descriptor(
                index=index,
                attempt_id=attempt_id,
                variant=variant,
                seed=seed,
                selected_ps_p0=selected_pressure,
                candidate=candidate,
                artifact_root=artifact_root,
                contract_path=contract_path,
                contract_payload=contract_payload,
            )
            descriptor_path = f"attempts/baseline/{attempt_id}.json"
            descriptor_payload = _json_bytes(descriptor)
            write(descriptor_path, descriptor_payload)
            descriptor_bindings.append(
                _binding(descriptor_path, _sha256(descriptor_payload))
            )
            descriptors.append(descriptor)
    campaign_plan["baseline_attempt_descriptors"] = descriptor_bindings

    source_attempt = next(
        descriptor
        for descriptor in descriptors
        if descriptor["variant"] == "three_level_amr_root_dx12_finest_dx3"
        and descriptor["qualifying_seed"]
        == campaign_plan["campaign_matrix"]["qualifying_seeds"][0]
    )
    restart_policy = planner.restart.decode_preregistration(
        source_payloads["restart_preregistration"].decode("utf-8")
    )
    planner.restart.validate_preregistration(restart_policy)
    carrier_id = planner._restart_carrier_id(source_attempt["qualifying_seed"])
    restart_artifact_root = campaign_root / "restart_continuation" / carrier_id
    restart_contract_path = f"launch_contracts/restart_continuation/{carrier_id}.json"
    restart_contract = planner._restart_launch_contract(
        carrier_id=carrier_id,
        source_attempt=source_attempt,
        restart_preregistration=restart_policy,
        candidate=candidate,
        paper_deck_binding=source_bindings["paper_deck"],
        artifact_root=restart_artifact_root,
    )
    restart_contract_payload = _json_bytes(restart_contract)
    write(restart_contract_path, restart_contract_payload)
    restart_contract_binding = _binding(
        restart_contract_path, _sha256(restart_contract_payload)
    )
    restart_carrier = {
        "record_type": planner.RESTART_CARRIER_RECORD_TYPE,
        "schema_version": 1,
        "carrier_id": carrier_id,
        "status": "planned_not_authorized",
        "source_baseline_attempt_id": source_attempt["attempt_id"],
        "variant": source_attempt["variant"],
        "qualifying_seed": source_attempt["qualifying_seed"],
        "selected_problem_ps_p0": selected_pressure,
        "authorized_orion_attempt_root": str(restart_artifact_root),
        "restart_preregistration": source_bindings["restart_preregistration"],
        "checkpoint_nominal_slot_omega0_inverse": restart_policy[
            "continuation_contract"
        ][
            "checkpoint_nominal_slot_omega0_inverse"
        ],
        "checkpoint_observed_commit_binding_required": True,
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": (
            restart_policy["continuation_contract"][
                "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
            ]
        ),
        "retained_output_pairing_policy": restart_policy["continuation_contract"][
            "retained_output_pairing_policy"
        ],
        "comparison_tolerances_max_absolute_difference": restart_policy[
            "continuation_contract"
        ]["comparison_tolerances_max_absolute_difference"],
        "launch_contract": restart_contract_binding,
    }
    restart_carrier_payload = _json_bytes(restart_carrier)
    restart_carrier_path = "restart_continuation/amr_restart_continuation_carrier.json"
    write(restart_carrier_path, restart_carrier_payload)
    campaign_plan["restart_continuation_carrier"] = _binding(
        restart_carrier_path, _sha256(restart_carrier_payload)
    )

    write("helper_source_closure.json", helper_closure_payload)
    campaign_plan["helper_source_closure"] = _binding(
        "helper_source_closure.json", _sha256(helper_closure_payload)
    )
    recompute = planner._independent_recompute_plan(
        plan_id=plan_id,
        campaign_root=campaign_root,
        qualifying_preregistration_binding=source_bindings["qualifying_preregistration"],
    )
    recompute_payload = _json_bytes(recompute)
    write("independent_raw_artifact_recompute_plan.json", recompute_payload)
    campaign_plan["independent_raw_artifact_recompute_plan"] = _binding(
        "independent_raw_artifact_recompute_plan.json", _sha256(recompute_payload)
    )
    fragment = control_plane_common._planner_expected_policy_fragment(
        plan_id=plan_id,
        pic_root=campaign.ORION_BULK_ROOT,
        campaign_root=campaign_root,
        candidate=candidate,
        pressure_receipt_binding=source_bindings["pressure_selection_receipt"],
        contract_bindings=contract_bindings,
        restart_contract_binding=restart_contract_binding,
    )
    fragment_payload = _json_bytes(fragment)
    write("nonauthorizing_policy_fragment.json", fragment_payload)
    campaign_plan["nonauthorizing_policy_fragment"] = _binding(
        "nonauthorizing_policy_fragment.json", _sha256(fragment_payload)
    )
    campaign_plan_payload = _json_bytes(campaign_plan)
    write("campaign_plan.json", campaign_plan_payload)
    materialized_inventory = _inventory_payload(materialized_members)
    internal_receipt_payload = _json_bytes(
        {
            "record_type": campaign._PLANNER_MATERIALIZATION_RECEIPT_RECORD_TYPE,
            "schema_version": 1,
            "plan_id": plan_id,
            "campaign_plan": _binding(
                "campaign_plan.json", _sha256(campaign_plan_payload)
            ),
            "helper_source_closure": _binding(
                "helper_source_closure.json", _sha256(helper_closure_payload)
            ),
            "tree_inventory": {
                "algorithm": (
                    campaign._PLANNER_MATERIALIZED_MEMBER_INVENTORY_ALGORITHM
                ),
                "scope": campaign._PLANNER_MATERIALIZED_MEMBER_INVENTORY_SCOPE,
                "excludes": [
                    campaign._PLANNER_MATERIALIZATION_RECEIPT_NAME,
                    immutable_orion_tree.FREEZE_RECEIPT_NAME,
                    immutable_orion_tree.INVENTORY_NAME,
                ],
                "sha256": _sha256(materialized_inventory),
                "inventoried_file_count": len(materialized_members),
            },
        }
    )
    (planner_root / campaign._PLANNER_MATERIALIZATION_RECEIPT_NAME).write_bytes(
        internal_receipt_payload
    )
    frozen = immutable_orion_tree.freeze_tree(
        planner_root,
        _PLANNER_RECEIPT,
        authorized_root=parent,
    )
    payload = _json_bytes(
        {
            "plan_root": str(planner_root),
            "plan_id": plan_id,
            "campaign_plan_sha256": _sha256(campaign_plan_payload),
            "materialization_receipt": _binding(
                campaign._PLANNER_MATERIALIZATION_RECEIPT_NAME,
                _sha256(internal_receipt_payload),
            ),
            "materialized_member_inventory_sha256": _sha256(materialized_inventory),
            "inventory_sha256": frozen["inventory_sha256"],
            "inventoried_file_count": frozen["inventoried_file_count"],
            "baseline_attempt_count": 24,
            "restart_continuation_carrier_count": 1,
            "recursively_read_only": True,
        }
    )
    return (
        _write_file(
            root, "bindings/q011_section54_campaign_plan.json", campaign_plan_payload
        ),
        _write_file(
            root,
            "bindings/q011_section54_planner_materialization_receipt.json",
            payload,
        ),
    )


def _minimal_planner_materialization_receipt(
    root: Path,
    campaign_plan_payload: bytes,
    helper_closure_payload: bytes,
) -> dict[str, str]:
    """Publish the formerly accepted two-member self-authored planner tree."""
    plan_id = json.loads(campaign_plan_payload)["plan_id"]
    parent = root.parent / "minimal-plans"
    parent.mkdir(exist_ok=True)
    planner_root = parent / f"q011-section54-qualifying-campaign-plan-{plan_id}"
    planner_root.mkdir()
    materialized_members = {
        "campaign_plan.json": campaign_plan_payload,
        "helper_source_closure.json": helper_closure_payload,
    }
    for relative, payload in materialized_members.items():
        (planner_root / relative).write_bytes(payload)
    materialized_inventory = _inventory_payload(materialized_members)
    internal_receipt_payload = _json_bytes(
        {
            "record_type": campaign._PLANNER_MATERIALIZATION_RECEIPT_RECORD_TYPE,
            "schema_version": 1,
            "plan_id": plan_id,
            "campaign_plan": _binding(
                "campaign_plan.json", _sha256(campaign_plan_payload)
            ),
            "helper_source_closure": _binding(
                "helper_source_closure.json", _sha256(helper_closure_payload)
            ),
            "tree_inventory": {
                "algorithm": campaign._PLANNER_MATERIALIZED_MEMBER_INVENTORY_ALGORITHM,
                "scope": campaign._PLANNER_MATERIALIZED_MEMBER_INVENTORY_SCOPE,
                "excludes": [
                    campaign._PLANNER_MATERIALIZATION_RECEIPT_NAME,
                    immutable_orion_tree.FREEZE_RECEIPT_NAME,
                    immutable_orion_tree.INVENTORY_NAME,
                ],
                "sha256": _sha256(materialized_inventory),
                "inventoried_file_count": len(materialized_members),
            },
        }
    )
    (planner_root / campaign._PLANNER_MATERIALIZATION_RECEIPT_NAME).write_bytes(
        internal_receipt_payload
    )
    frozen = immutable_orion_tree.freeze_tree(
        planner_root,
        _PLANNER_RECEIPT,
        authorized_root=parent,
    )
    return _write_file(
        root,
        "bindings/q011_section54_planner_materialization_receipt.json",
        _json_bytes(
            {
                "plan_root": str(planner_root),
                "plan_id": plan_id,
                "campaign_plan_sha256": _sha256(campaign_plan_payload),
                "materialization_receipt": _binding(
                    campaign._PLANNER_MATERIALIZATION_RECEIPT_NAME,
                    _sha256(internal_receipt_payload),
                ),
                "materialized_member_inventory_sha256": _sha256(
                    materialized_inventory
                ),
                "inventory_sha256": frozen["inventory_sha256"],
                "inventoried_file_count": frozen["inventoried_file_count"],
                "baseline_attempt_count": 24,
                "restart_continuation_carrier_count": 1,
                "recursively_read_only": True,
            }
        ),
    )


def _model_overrides(variant: str) -> list[str]:
    return list(campaign._EXPECTED_MODEL_LAUNCH_OVERRIDES[variant])


def _attempt_id(variant: str, seed: int) -> str:
    return campaign._expected_attempt_id({"variant": variant, "seed": seed})


def _contract_argv(
    attempt_id: str, variant: str, seed: int, pressure: float, plan_id: str
) -> list[str]:
    campaign_root = (
        campaign.ORION_BULK_ROOT / "campaigns" / f"q011-section54-{plan_id}"
    )
    attempt_root = f"{campaign_root}/baseline/{attempt_id}"
    return [
        "-i",
        "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
        "-d",
        f"{attempt_root}/raw",
        f"job/basename={attempt_id}",
        f"problem/ps_p0={pressure!r}",
        f"particles/pic_random_seed={seed}",
        f"problem/ps_inject_seed={seed}",
        f"problem/ps_seed_noise_seed={seed}",
        *_model_overrides(variant),
    ]


def _registered_execution_receipt(
    *,
    attempt_id: str,
    source_commit: str,
    executable_sha256: str,
    deck_sha256: str,
    environment_sha256: str,
    control_plane_version: str,
    argv: list[str],
    raw_output_root: str,
    planner_result: dict[str, Any],
    artifact_dir: str = "/fixture/runs/q011/submission",
) -> dict[str, Any]:
    attempt_root = str(Path(raw_output_root).parent)
    return {
        "record_type": "q011_section54_reconciled_registered_execution_receipt",
        "schema_version": 1,
        "receipt_role": "immutable_reconciled_registered_execution",
        "registration_scope": "registered_science",
        "reconciled": True,
        "reservation_id": "11111111-1111-4111-8111-111111111111",
        "submission_id": "22222222-2222-4222-8222-222222222222",
        "reconciliation_event_sha256": "d" * 64,
        "attempt_id": attempt_id,
        "source_commit": source_commit,
        "executable_sha256": executable_sha256,
        "deck_sha256": deck_sha256,
        "environment_sha256": environment_sha256,
        "control_plane_version": control_plane_version,
        "argv": argv,
        "slurm_job_id": "123456",
        "slurm_terminal_state": "COMPLETED",
        "raw_output_root": raw_output_root,
        "artifact_dir": artifact_dir,
        "planner_retention": {
            "schema_version": 1,
            "retention_role": "q011_section54_deterministic_retained_attempt",
            "planner_root": planner_result["plan_root"],
            "planner_inventory_sha256": planner_result["inventory_sha256"],
            "planner_plan_id": planner_result["plan_id"],
            "planner_materialization_receipt": planner_result[
                "materialization_receipt"
            ],
            "attempt_id": attempt_id,
            "authorized_orion_attempt_root": attempt_root,
            "authorized_orion_raw_root": raw_output_root,
            "argv": argv,
        },
        "pre_submit_manifest_sha256": "e" * 64,
    }


def _product(
    root: Path,
    kind: str,
    relative: str,
    payload: bytes,
    snapshot_time: float | None,
) -> dict[str, Any]:
    binding = _write_file(root, relative, payload)
    return {
        "kind": kind,
        "path": binding["path"],
        "sha256": binding["sha256"],
        "nominal_slot_time": snapshot_time,
        "observed_committed_time": snapshot_time,
    }


def _frozen_candidate_path(suffix: str) -> str:
    return (
        campaign.ORION_BULK_ROOT / "clean_candidates" / _FREEZE_ID / suffix
    ).as_posix()


def _clean_candidate_manifest(
    executable: dict[str, str],
    deck: dict[str, str],
    analyzer: dict[str, str],
    *,
    source_archive_sha256: str,
) -> dict[str, Any]:
    return {
        "schema_version": 4,
        "freeze_id": _FREEZE_ID,
        "created_utc": "2026-06-01T12:00:00Z",
        "prepared_artifacts": {
            "inventory_path": campaign.PREPARED_ARTIFACT_INVENTORY_SOURCE_PATH,
            "inventory_sha256": "c" * 64,
            "paper_decks": [
                {"path": campaign.ACTIVE_DECK_SOURCE_PATH, "sha256": deck["sha256"]},
                {"path": "inputs/tests/pic_fixture.athinput", "sha256": "d" * 64},
            ],
            "analyzers": [
                {
                    "path": campaign.CAMPAIGN_ANALYZER_SOURCE_PATH,
                    "sha256": analyzer["sha256"],
                },
                {
                    "path": "tst/publication/analyze_q011_section54_outputs.py",
                    "sha256": "e" * 64,
                },
            ],
        },
        "source": {
            "archive_path": _frozen_candidate_path("source.tar"),
            "archive_sha256": source_archive_sha256,
            "commit_path": _frozen_candidate_path("source.commit"),
            "commit_sha256": "4" * 64,
            "source_bundle_sha256": "2" * 64,
            "git_commit": "1" * 40,
            "git_tree": "5" * 40,
            "worktree_status": "clean",
            "submodule_status": "absent",
            "submodules": [],
        },
        "build": {
            "profile_id": "hip-mpi-release-paper-pic",
            "profile_path": _frozen_candidate_path("build_profile.json"),
            "profile_sha256": "6" * 64,
            "profile_receipt_path": _frozen_candidate_path("profile_receipt.json"),
            "profile_receipt_sha256": "7" * 64,
            "source_archive_sha256": source_archive_sha256,
            "source_commit_sha256": "4" * 64,
            "source_bundle_sha256": "2" * 64,
            "toolchain": "Frontier fixture toolchain",
            "build_invocations_sha256": "8" * 64,
            "executable_path": _frozen_candidate_path("athena"),
            "executable_sha256": executable["sha256"],
        },
    }


def _append_restart_publication(
    root: Path, products: list[dict[str, Any]], time: float, cycle: int
) -> None:
    restart = _product(
        root, "restart", f"rst/q011.{cycle:05d}.rst", _restart_payload(cycle), time
    )
    products.append(restart)
    products.append(
        _product(
            root,
            "restart_complete",
            restart["path"] + ".complete",
            _restart_marker((root / restart["path"]).read_bytes()),
            time,
        )
    )
    manifest_payload = (
        json.dumps(
            {
                "schema": "ATHENAK_RESTART_MANIFEST_V1",
                "members": [
                    {
                        "path": restart["path"],
                        "size": len((root / restart["path"]).read_bytes()),
                        "fnv1a64": f"{_fnv1a64((root / restart['path']).read_bytes()):016x}",
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    manifest = _product(
        root, "restart_manifest", restart["path"] + ".manifest", manifest_payload, time
    )
    products.append(manifest)
    products.append(
        _product(
            root,
            "restart_manifest_complete",
            manifest["path"] + ".complete",
            _restart_marker(manifest_payload),
            time,
        )
    )


def _manifest(
    root: Path,
    planner_mutate: Callable[[dict[str, Any], dict[str, bytes]], None] | None = None,
) -> dict[str, Any]:
    variant = "three_level_amr_root_dx12_finest_dx3"
    seed = 23050101
    attempt_id = _attempt_id(variant, seed)
    selected_pressure = 0.1
    executable = _write_file(root, "bindings/athena", _elf_executable())
    (root / executable["path"]).chmod(0o755)
    deck = _write_file(
        root,
        "bindings/section54.athinput",
        (campaign.REPO_ROOT / campaign.ACTIVE_DECK_SOURCE_PATH).read_bytes(),
    )
    analyzer = _write_file(
        root,
        "bindings/analyze_q011_section54_campaign.py",
        campaign.ANALYZER_PATH.read_bytes(),
    )
    candidate_root = campaign.ORION_BULK_ROOT / "clean_candidates" / _FREEZE_ID
    candidate_root.mkdir(parents=True)
    source_archive_path = candidate_root / "source.tar"
    with tarfile.open(source_archive_path, mode="w") as archive:
        for relative in sorted(
            {
                *control_plane_common.Q011_SECTION54_HELPER_SOURCES,
                *control_plane_common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS.values(),
            }
        ):
            payload = (campaign.REPO_ROOT / relative).read_bytes()
            member = tarfile.TarInfo(relative)
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    source_archive_path.chmod(0o444)
    source_archive_sha256 = _sha256(source_archive_path.read_bytes())
    candidate_payload = (
        json.dumps(
            _clean_candidate_manifest(
                executable,
                deck,
                analyzer,
                source_archive_sha256=source_archive_sha256,
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    (candidate_root / "clean_candidate_manifest.json").write_bytes(candidate_payload)
    (candidate_root / "clean_candidate_manifest.json").chmod(0o444)
    candidate_manifest = _write_file(
        root, "bindings/clean_candidate_manifest.json", candidate_payload
    )
    preregistration = _write_file(
        root,
        "bindings/q011_section54_preregistration.json",
        campaign.PREREGISTRATION_PATH.read_bytes(),
    )
    pressure_descriptors = [
        {
            "case_id": case_id,
            "problem_ps_p0": pressure,
            "descriptor_sha256": f"{index:064x}",
        }
        for index, (case_id, pressure) in enumerate(
            campaign.pressure_selection.REGISTERED_CASES, start=1
        )
    ]
    published_pressure_receipt = _published_pressure_receipt(
        root, pressure_descriptors
    )
    published_review_packet_receipt = {
        "path": str(
            root.parent
            / "publication/q011_section54_pressure_pilot_review_packet_receipt.json"
        ),
        "sha256": "d" * 64,
    }
    authoritative_reanalysis_attestation = {
        "path": str(
            root.parent
            / "pressure_gate_attestations"
            / "20260605T091100Z-q011-section54-pressure-reanalysis-fixture.operator"
            / "attestation.json"
        ),
        "sha256": "e" * 64,
    }
    reviewer_attestation = {
        "path": str(
            root.parent
            / "pressure_gate_attestations"
            / "20260605T091300Z-q011-section54-pressure-selection-fixture.reviewer"
            / "attestation.json"
        ),
        "sha256": "f" * 64,
    }
    selected_pressure_payload = _json_bytes(
        {
            "schema_version": 3,
            "record_type": "q011_section54_pressure_selection_receipt",
            "selection_method": "human_review_only",
            "published_pressure_pilot_receipt": published_pressure_receipt,
            "published_pressure_pilot_review_packet_receipt": (
                published_review_packet_receipt
            ),
            "pilot_bundle_manifest_sha256": "b" * 64,
            "aggregate_pilot_analysis_sha256": "c" * 64,
            "case_descriptors": pressure_descriptors,
            "selected_case": {
                "case_id": "ps_p0_0p10",
                "problem_ps_p0": selected_pressure,
            },
            "authoritative_reanalysis_attestation": (
                authoritative_reanalysis_attestation
            ),
            "reviewer_attestation": reviewer_attestation,
        }
    )
    selected_pressure_receipt = _write_file(
        root,
        "bindings/q011_section54_selected_pressure_receipt.json",
        selected_pressure_payload,
    )
    helper_sources = [
        {
            "path": path,
            "sha256": _sha256((campaign.REPO_ROOT / path).read_bytes()),
        }
        for path in campaign._EXPECTED_HELPER_SOURCE_PATHS
    ]
    environment_profile_payload = (
        campaign.REPO_ROOT
        / control_plane_common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS[
            "environment_profile"
        ]
    ).read_bytes()
    environment_profile = _binding(
        "bindings/environment_profile.sh", _sha256(environment_profile_payload)
    )
    restart_preregistration_payload = campaign.RESTART_PREREGISTRATION_PATH.read_bytes()
    restart_preregistration = _binding(
        "bindings/q011_section54_restart_continuation_preregistration.json",
        _sha256(restart_preregistration_payload),
    )
    paper_deck = _binding(
        "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
        deck["sha256"],
    )
    plan_candidate = {
        "clean_candidate_manifest": {
            "path": _frozen_candidate_path("clean_candidate_manifest.json"),
            "sha256": candidate_manifest["sha256"],
        },
        "freeze_id": _FREEZE_ID,
        "git_commit": "1" * 40,
        "git_tree": "5" * 40,
        "source_archive_sha256": source_archive_sha256,
        "source_commit_sha256": "4" * 64,
        "source_bundle_sha256": "2" * 64,
        "prepared_artifact_inventory_sha256": "c" * 64,
        "validated_submodules": [],
        "build_profile": {
            "path": _frozen_candidate_path("build_profile.json"),
            "sha256": "6" * 64,
        },
        "build_profile_receipt": {
            "path": _frozen_candidate_path("profile_receipt.json"),
            "sha256": "7" * 64,
        },
        "build_invocations_sha256": "8" * 64,
        "executable": {
            "path": _frozen_candidate_path("athena"),
            "sha256": executable["sha256"],
        },
        "environment_profile": {
            "path": str(
                campaign.ORION_BULK_ROOT
                / "control_plane"
                / ("9" * 64)
                / "frontier_pic_environment.sh"
            ),
            "sha256": environment_profile["sha256"],
            "control_plane_version": "9" * 64,
            "reviewed_source": {
                "path": "tst/publication/frontier_control_plane/frontier_pic_environment.sh",
                "sha256": _sha256(
                    (
                        campaign.REPO_ROOT
                        / "tst/publication/frontier_control_plane/frontier_pic_environment.sh"
                    ).read_bytes()
                ),
            },
        },
    }
    planner_source_bindings = {
        "pressure_selection_receipt": _binding(
            "bindings/human_pressure_selection_receipt.json",
            selected_pressure_receipt["sha256"],
        ),
        "clean_candidate_manifest": candidate_manifest,
        "environment_profile": environment_profile,
        "qualifying_preregistration": _binding(
            "bindings/q011_section54_qualifying_campaign_preregistration.json",
            preregistration["sha256"],
        ),
        "restart_preregistration": restart_preregistration,
        "paper_deck": paper_deck,
    }
    campaign_matrix = copy.deepcopy(
        campaign._EXPECTED_POLICY_PROJECTION["campaign_matrix"]
    )
    selected_case = {
        "case_id": "ps_p0_0p10",
        "problem_ps_p0": selected_pressure,
    }
    plan_id = campaign.campaign_planner._digest_value(
        {
            "record_type": campaign.campaign_planner.PLAN_RECORD_TYPE,
            "schema_version": 1,
            "pressure_selection_receipt_sha256": planner_source_bindings[
                "pressure_selection_receipt"
            ]["sha256"],
            "selected_case": selected_case,
            "candidate_binding": plan_candidate,
            "source_binding_sha256": {
                name: binding["sha256"]
                for name, binding in planner_source_bindings.items()
            },
            "helper_source_closure": helper_sources,
            "campaign_matrix": campaign_matrix,
            "authorized_orion_root": str(campaign.ORION_BULK_ROOT),
        }
    )
    planned_campaign_root = str(
        campaign.ORION_BULK_ROOT / "campaigns" / f"q011-section54-{plan_id}"
    )
    planned_attempt_root = f"{planned_campaign_root}/baseline/{attempt_id}"
    helper_closure_payload = _json_bytes(
        {
            "record_type": "q011_section54_helper_source_closure",
            "schema_version": 1,
            "plan_id": plan_id,
            "sources": helper_sources,
        }
    )
    analyzer_helper_source_closure_manifest = _write_file(
        root,
        "bindings/q011_section54_analyzer_helper_source_closure_manifest.json",
        helper_closure_payload,
    )
    campaign_plan_value = {
            "record_type": "q011_section54_qualifying_campaign_execution_plan",
            "schema_version": 1,
            "plan_id": plan_id,
            "artifact_role": (
                "q011_section54_source_local_immutable_qualifying_campaign_plan"
            ),
            "qualification_effect": (
                "plan_only_no_execution_authorization_no_claim_closure"
            ),
            "status": "source_local_immutable_review_plan_only",
            "authorized_orion_root": str(campaign.ORION_BULK_ROOT),
            "authorized_orion_campaign_root": planned_campaign_root,
            "selected_pressure": {
                "selection_method": "human_review_only",
                "selected_case": selected_case,
                "receipt": planner_source_bindings["pressure_selection_receipt"],
            },
            "candidate_binding": plan_candidate,
            "source_bindings": planner_source_bindings,
            "helper_source_closure": _binding(
                "helper_source_closure.json",
                analyzer_helper_source_closure_manifest["sha256"],
            ),
            "campaign_matrix": campaign_matrix,
            "baseline_attempt_count": 24,
            "baseline_attempt_descriptors": [
                _binding(f"attempts/baseline/attempt-{index:03d}.json", f"{index:064x}")
                for index in range(1, 25)
            ],
            "restart_continuation_carrier": _binding(
                "attempts/restart/carrier.json", "d" * 64
            ),
            "independent_raw_artifact_recompute_plan": _binding(
                "independent_recompute_plan.json", "e" * 64
            ),
            "nonauthorizing_policy_fragment": _binding(
                "nonauthorizing_policy_fragment.json", "f" * 64
            ),
            "execution_boundary": {
                "mutates_live_policy": False,
                "scheduler_calls": False,
                "submits_jobs": False,
                "infers_pressure_selection": False,
                "launch_authorized": False,
                "frontier_execution_authorized": False,
                "claim_closure_authorized": False,
            },
            "preregistration_execution_boundary": json.loads(
                campaign.PREREGISTRATION_PATH.read_text(encoding="utf-8")
            )["qualifying_execution_bindings"],
    }
    planner_source_payloads = {
        "pressure_selection_receipt": selected_pressure_payload,
        "clean_candidate_manifest": candidate_payload,
        "environment_profile": environment_profile_payload,
        "qualifying_preregistration": campaign.PREREGISTRATION_PATH.read_bytes(),
        "restart_preregistration": restart_preregistration_payload,
        "paper_deck": (root / deck["path"]).read_bytes(),
    }
    if planner_mutate is not None:
        planner_mutate(campaign_plan_value, planner_source_payloads)
    campaign_plan, planner_materialization_receipt = _planner_materialization_receipt(
        root, campaign_plan_value, planner_source_payloads, helper_closure_payload
    )
    planner_result = json.loads(
        (root / planner_materialization_receipt["path"]).read_text(encoding="utf-8")
    )
    attempt_contract_payload = _json_bytes(
        {
            "record_type": "q011_section54_launch_prohibited_handoff_contract",
            "schema_version": 1,
            "contract_role": "source_local_review_handoff_only",
            "launch_authorized": False,
            "scheduler_submission_authorized": False,
            "live_policy_mutation_authorized": False,
            "attempt_id": attempt_id,
            "variant": variant,
            "qualifying_seed": seed,
            "selected_problem_ps_p0": selected_pressure,
            "executable": plan_candidate["executable"],
            "environment_profile": plan_candidate["environment_profile"],
            "paper_deck": paper_deck,
            "authorized_orion_attempt_root": planned_attempt_root,
            "argv": _contract_argv(
                attempt_id, variant, seed, selected_pressure, plan_id
            ),
            "required_separate_boundary": (
                "review_and_promote_a_registered_frontier_submission_policy_then_use_"
                "the_installed_control_plane_wrapper"
            ),
        }
    )
    attempt_contract = _write_file(
        root, "bindings/q011_section54_attempt_contract.json", attempt_contract_payload
    )
    registered_execution_receipt = _write_file(
        root,
        "bindings/q011_section54_registered_execution_receipt.json",
        _json_bytes(
            _registered_execution_receipt(
                attempt_id=attempt_id,
                source_commit=plan_candidate["git_commit"],
                executable_sha256=plan_candidate["executable"]["sha256"],
                deck_sha256=paper_deck["sha256"],
                environment_sha256=plan_candidate["environment_profile"]["sha256"],
                control_plane_version=plan_candidate["environment_profile"][
                    "control_plane_version"
                ],
                argv=_contract_argv(
                    attempt_id, variant, seed, selected_pressure, plan_id
                ),
                raw_output_root=f"{planned_attempt_root}/raw",
                planner_result=planner_result,
            )
        ),
    )
    products = []
    for cycle, time in enumerate(_TIMES):
        suffix = f"{int(time):05d}"
        for kind in ("rho", "bmag", "prtcl_jx", "j2"):
            products.append(
                _product(
                    root,
                    kind,
                    f"bin/q011.{kind}.{suffix}.bin",
                    _mesh_bin(kind, time, cycle),
                    time,
                )
            )
        products.append(
            _product(
                root,
                "prtcl_all",
                f"pvtk/q011.prtcl_all.{suffix}.part.vtk",
                _particle_vtk(time, cycle),
                time,
            )
        )
        _append_restart_publication(root, products, time, cycle)
    products.append(_product(root, "stdout", "stdout.txt", _stdout_telemetry(), None))
    return {
        "schema_version": 1,
        "record_type": "q011_section54_campaign_run_manifest",
        "campaign_id": campaign.CAMPAIGN_ID,
        "qualification_scope": campaign.QUALIFICATION_SCOPE,
        "authorized_orion_campaign_root": str(root),
        "run_identity": {
            "variant": variant,
            "seed": seed,
            "physical_mode": "paper_mhd_pic_vl2_tsc",
            "attempt_id": attempt_id,
        },
        "candidate_binding": {
            "git_commit": "1" * 40,
            "source_bundle_sha256": "2" * 64,
            "clean_candidate_manifest": candidate_manifest,
        },
        "artifact_bindings": {
            "executable": executable,
            "deck": deck,
            "analyzer": analyzer,
            "preregistration": preregistration,
            "campaign_plan": campaign_plan,
            "planner_materialization_receipt": planner_materialization_receipt,
            "attempt_contract": attempt_contract,
            "selected_pressure_receipt": selected_pressure_receipt,
            "analyzer_helper_source_closure_manifest": (
                analyzer_helper_source_closure_manifest
            ),
            "registered_execution_receipt": registered_execution_receipt,
        },
        "attempt_identity": {
            "campaign_plan_sha256": campaign_plan["sha256"],
            "planner_materialization_receipt_sha256": (
                planner_materialization_receipt["sha256"]
            ),
            "attempt_contract_sha256": attempt_contract["sha256"],
            "selected_pressure_receipt_sha256": selected_pressure_receipt["sha256"],
            "model_launch_overrides": _model_overrides(variant),
            "seed_overrides": {
                "particles/pic_random_seed": 23050101,
                "problem/ps_inject_seed": 23050101,
                "problem/ps_seed_noise_seed": 23050101,
            },
            "analyzer_helper_source_closure_manifest_sha256": (
                analyzer_helper_source_closure_manifest["sha256"]
            ),
            "registered_execution_receipt_sha256": (
                registered_execution_receipt["sha256"]
            ),
        },
        "products": products,
    }


def _find_product(manifest: dict[str, Any], kind: str, time: float | None) -> dict[str, Any]:
    return next(
        product
        for product in manifest["products"]
        if product["kind"] == kind and product["nominal_slot_time"] == time
    )


def _rewrite_product(root: Path, product: dict[str, Any], payload: bytes) -> None:
    (root / product["path"]).write_bytes(payload)
    product["sha256"] = _sha256(payload)


def _rewrite_clean_candidate(
    root: Path,
    manifest: dict[str, Any],
    mutate: Callable[[dict[str, Any]], None],
) -> None:
    binding = manifest["candidate_binding"]["clean_candidate_manifest"]
    path = root / binding["path"]
    candidate = json.loads(path.read_text(encoding="utf-8"))
    mutate(candidate)
    payload = (json.dumps(candidate, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path.write_bytes(payload)
    binding["sha256"] = _sha256(payload)


def _rewrite_bound_artifact(
    root: Path, manifest: dict[str, Any], name: str, payload: bytes
) -> None:
    binding = manifest["artifact_bindings"][name]
    (root / binding["path"]).write_bytes(payload)
    binding["sha256"] = _sha256(payload)


def _rewrite_attempt_json_artifact(
    root: Path,
    manifest: dict[str, Any],
    name: str,
    mutate: Callable[[dict[str, Any]], None],
) -> None:
    binding = manifest["artifact_bindings"][name]
    value = json.loads((root / binding["path"]).read_text(encoding="utf-8"))
    mutate(value)
    _rewrite_bound_artifact(root, manifest, name, _json_bytes(value))
    digest_name = campaign._ATTEMPT_BINDING_SHA256_FIELDS[name]
    manifest["attempt_identity"][digest_name] = manifest["artifact_bindings"][name][
        "sha256"
    ]


def _retarget_attempt(
    root: Path, manifest: dict[str, Any], variant: str, seed: int = 23050101
) -> None:
    attempt_id = _attempt_id(variant, seed)
    plan_id = json.loads(
        (
            root / manifest["artifact_bindings"]["campaign_plan"]["path"]
        ).read_text(encoding="utf-8")
    )["plan_id"]
    manifest["run_identity"].update(
        {"variant": variant, "seed": seed, "attempt_id": attempt_id}
    )
    manifest["attempt_identity"]["model_launch_overrides"] = _model_overrides(variant)
    manifest["attempt_identity"]["seed_overrides"] = {
        name: seed for name in campaign._SEED_OVERRIDE_NAMES
    }

    def mutate_contract(contract: dict[str, Any]) -> None:
        pressure = contract["selected_problem_ps_p0"]
        contract.update(
            {
                "attempt_id": attempt_id,
                "variant": variant,
                "qualifying_seed": seed,
                "authorized_orion_attempt_root": (
                    f"{campaign.ORION_BULK_ROOT}/campaigns/"
                    f"q011-section54-{plan_id}/baseline/{attempt_id}"
                ),
                "argv": _contract_argv(attempt_id, variant, seed, pressure, plan_id),
            }
        )

    _rewrite_attempt_json_artifact(root, manifest, "attempt_contract", mutate_contract)
    contract = json.loads(
        (
            root / manifest["artifact_bindings"]["attempt_contract"]["path"]
        ).read_text(encoding="utf-8")
    )

    def mutate_registered_execution_receipt(receipt: dict[str, Any]) -> None:
        raw_output_root = f"{contract['authorized_orion_attempt_root']}/raw"
        receipt.update(
            {
                "attempt_id": attempt_id,
                "argv": contract["argv"],
                "raw_output_root": raw_output_root,
                "planner_retention": {
                    "schema_version": 1,
                    "retention_role": (
                        "q011_section54_deterministic_retained_attempt"
                    ),
                    "planner_root": receipt["planner_retention"]["planner_root"],
                    "planner_inventory_sha256": receipt["planner_retention"][
                        "planner_inventory_sha256"
                    ],
                    "planner_plan_id": receipt["planner_retention"]["planner_plan_id"],
                    "planner_materialization_receipt": receipt["planner_retention"][
                        "planner_materialization_receipt"
                    ],
                    "attempt_id": attempt_id,
                    "authorized_orion_attempt_root": contract[
                        "authorized_orion_attempt_root"
                    ],
                    "authorized_orion_raw_root": raw_output_root,
                    "argv": contract["argv"],
                },
            }
        )

    _rewrite_attempt_json_artifact(
        root,
        manifest,
        "registered_execution_receipt",
        mutate_registered_execution_receipt,
    )


def _rewrite_restart_payload(
    root: Path, manifest: dict[str, Any], time: float, payload: bytes
) -> None:
    restart = _find_product(manifest, "restart", time)
    restart_marker = _find_product(manifest, "restart_complete", time)
    restart_manifest = _find_product(manifest, "restart_manifest", time)
    manifest_marker = _find_product(manifest, "restart_manifest_complete", time)
    _rewrite_product(root, restart, payload)
    _rewrite_product(root, restart_marker, _restart_marker(payload))
    publication = json.loads((root / restart_manifest["path"]).read_text(encoding="utf-8"))
    publication["members"][0]["size"] = len(payload)
    publication["members"][0]["fnv1a64"] = f"{_fnv1a64(payload):016x}"
    publication_payload = (
        json.dumps(publication, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    _rewrite_product(root, restart_manifest, publication_payload)
    _rewrite_product(root, manifest_marker, _restart_marker(publication_payload))


def _rewrite_slot_observed_time(
    root: Path,
    manifest: dict[str, Any],
    nominal_slot_time: float,
    observed_committed_time: float,
) -> None:
    slot_products = [
        product
        for product in manifest["products"]
        if product["nominal_slot_time"] == nominal_slot_time
    ]
    for product in slot_products:
        product["observed_committed_time"] = observed_committed_time
        if product["kind"] in ("rho", "bmag", "prtcl_jx", "j2"):
            _rewrite_product(
                root,
                product,
                _mesh_bin(
                    product["kind"],
                    campaign._mesh_time_projection(observed_committed_time),
                    int(nominal_slot_time // 100),
                ),
            )
        elif product["kind"] == "prtcl_all":
            _rewrite_product(
                root,
                product,
                _particle_vtk(
                    observed_committed_time,
                    int(nominal_slot_time // 100),
                ),
            )


@contextmanager
def _frozen_fixture(
    mutate: Callable[[Path, dict[str, Any]], None] | None = None,
    *,
    planner_mutate: (
        Callable[[dict[str, Any], dict[str, bytes]], None] | None
    ) = None,
    receipt: dict[str, Any] | None = None,
    external_side_effect: Exception | None = None,
    reanalysis_source_snapshot_side_effect: Exception | None = None,
    wrong_retained_root: bool = False,
) -> Iterator[tuple[Path, Path, str]]:
    with tempfile.TemporaryDirectory() as directory:
        authorized_root = Path(directory)
        with mock.patch.object(campaign, "ORION_BULK_ROOT", authorized_root):
            staging_tree = authorized_root / "campaign-staging"
            staging_tree.mkdir()
            manifest = _manifest(staging_tree, planner_mutate)
            plan = json.loads(
                (
                    staging_tree
                    / manifest["artifact_bindings"]["campaign_plan"]["path"]
                ).read_text(encoding="utf-8")
            )
            planner_materialization = json.loads(
                (
                    staging_tree
                    / manifest["artifact_bindings"]["planner_materialization_receipt"][
                        "path"
                    ]
                ).read_text(encoding="utf-8")
            )
            planner_root = Path(planner_materialization["plan_root"])
            retained_pressure_receipt = json.loads(
                (
                    planner_root
                    / plan["source_bindings"]["pressure_selection_receipt"]["path"]
                ).read_text(encoding="utf-8")
            )
            retained_helper_closure = json.loads(
                (
                    planner_root / plan["helper_source_closure"]["path"]
                ).read_text(encoding="utf-8")
            )["sources"]
            if mutate is not None:
                mutate(staging_tree, manifest)
            if wrong_retained_root:
                tree = (
                    authorized_root
                    / "wrong-retained-root"
                    / manifest["run_identity"]["attempt_id"]
                )
            else:
                tree = (
                    Path(plan["authorized_orion_campaign_root"])
                    / "baseline"
                    / manifest["run_identity"]["attempt_id"]
                )
            tree.parent.mkdir(parents=True, exist_ok=True)
            staging_tree.rename(tree)
            manifest["authorized_orion_campaign_root"] = str(tree)
            (tree / campaign.MANIFEST_NAME).write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            patch_kwargs = (
                {"side_effect": external_side_effect}
                if external_side_effect is not None
                else {
                    "return_value": {
                        "freeze_id": _FREEZE_ID,
                        "referenced_manifest_sha256": manifest["candidate_binding"][
                            "clean_candidate_manifest"
                        ]["sha256"],
                        "retained_executable_sha256": manifest["artifact_bindings"][
                            "executable"
                        ]["sha256"],
                        "validated_submodule_count": 0,
                        "validation": "fixture_control_plane_closure",
                    }
                }
            )

            def consume_pressure_pilot_bundle() -> dict[str, str]:
                expected_path = (
                    authorized_root
                    / "publication/q011_section54_pressure_pilot_publication_receipt.json"
                )
                return {
                    "packet_receipt_sha256": "d" * 64,
                    "aggregate_receipt_sha256": _sha256(expected_path.read_bytes()),
                    "manifest_sha256": "b" * 64,
                    "analysis_result_sha256": "c" * 64,
                    "status": "pass_engineering_calibration_only",
                }

            def expected_pressure_review_packet_bindings(
                authorized_pic_root: Path,
            ) -> tuple[dict[str, str], dict[str, str]]:
                expected_packet = {
                    "path": str(
                        authorized_root
                        / "publication/q011_section54_pressure_pilot_review_packet_receipt.json"
                    ),
                    "sha256": "d" * 64,
                }
                expected_aggregate_path = (
                    authorized_root
                    / "publication/q011_section54_pressure_pilot_publication_receipt.json"
                )
                expected_aggregate = {
                    "path": str(expected_aggregate_path),
                    "sha256": _sha256(expected_aggregate_path.read_bytes()),
                }
                if authorized_pic_root != authorized_root:
                    raise campaign.pressure_selection.pressure_review_packet_verifier.PressureReviewPacketVerificationError(
                        "fixture pressure review-packet root drifted"
                    )
                return expected_packet, expected_aggregate

            def verify_pressure_review_packet(
                path: str | Path,
                *,
                aggregate_receipt_binding: object,
                authorized_pic_root: Path,
            ) -> dict[str, object]:
                expected_packet, expected_aggregate = (
                    expected_pressure_review_packet_bindings(authorized_pic_root)
                )
                if (
                    Path(path) != Path(expected_packet["path"])
                    or aggregate_receipt_binding != expected_aggregate
                ):
                    raise campaign.pressure_selection.pressure_review_packet_verifier.PressureReviewPacketVerificationError(
                        "fixture pressure review-packet binding drifted"
                    )
                return {
                    "receipt_binding": expected_packet,
                    "aggregate_receipt_binding": expected_aggregate,
                    "packet_receipt": {},
                    "aggregate_receipt": aggregate_record,
                    "aggregate_bundle": aggregate_record["aggregate_bundle"],
                    "aggregate_analysis": aggregate_record["aggregate_analysis"],
                    "source_bindings": {},
                    "inventory": {},
                }

            aggregate_record = {
                "aggregate_bundle": {
                    "path": str(
                        authorized_root / "publication" / "pressure-pilot-bundle"
                    ),
                    "manifest_sha256": retained_pressure_receipt[
                        "pilot_bundle_manifest_sha256"
                    ],
                },
                "aggregate_analysis": {
                    "path": str(
                        authorized_root
                        / "publication"
                        / "pressure-pilot-analysis.json"
                    ),
                    "sha256": retained_pressure_receipt[
                        "aggregate_pilot_analysis_sha256"
                    ],
                },
                "raw_cases": [
                    {
                        "case_id": descriptor["case_id"],
                        "artifact_dir": str(
                            authorized_root
                            / "runs"
                            / "pressure-pilot"
                            / descriptor["case_id"]
                        ),
                        "descriptor_path": str(
                            authorized_root
                            / "runs"
                            / "pressure-pilot"
                            / descriptor["case_id"]
                            / "descriptor.json"
                        ),
                        "descriptor_sha256": descriptor["descriptor_sha256"],
                        "artifact_inventory_sha256": "f" * 64,
                        "runtime_artifacts": {},
                    }
                    for descriptor in retained_pressure_receipt["case_descriptors"]
                ],
            }
            pressure_verifier = (
                campaign.pressure_selection.pressure_review_packet_verifier
            )
            expected_reanalysis_binding = retained_pressure_receipt[
                "authoritative_reanalysis_attestation"
            ]
            expected_reviewer_binding = retained_pressure_receipt[
                "reviewer_attestation"
            ]
            helper_by_path = {
                record["path"]: record["sha256"] for record in retained_helper_closure
            }
            reanalysis_source_closure = [
                {"path": path, "sha256": helper_by_path[path]}
                for path in pressure_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
            ]
            reanalysis_source_authorization = {
                "execution_mode": pressure_verifier.PRESSURE_REANALYSIS_EXECUTION_MODE,
                "git_commit": plan["candidate_binding"]["git_commit"],
                "source_archive_sha256": plan["candidate_binding"][
                    "source_archive_sha256"
                ],
                "source_closure_sha256": pressure_verifier._source_closure_digest(
                    reanalysis_source_closure
                ),
                "source_closure": reanalysis_source_closure,
                "historical_production_source_authorization": dict(
                    pressure_verifier.AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION
                ),
            }

            def reanalysis_verification(
                attestation_binding: object,
                *,
                aggregate_receipt_binding: object,
                packet_receipt_binding: object,
                pilot_bundle_manifest_sha256: object,
                aggregate_pilot_analysis_sha256: object,
                authorized_pic_root: Path,
                expected_result: object | None = None,
                now: object | None = None,
                reject_source_snapshot: bool,
            ) -> dict[str, object]:
                del now
                expected_packet, expected_aggregate = (
                    expected_pressure_review_packet_bindings(authorized_pic_root)
                )
                result = consume_pressure_pilot_bundle()
                if (
                    attestation_binding != expected_reanalysis_binding
                    or aggregate_receipt_binding != expected_aggregate
                    or packet_receipt_binding != expected_packet
                    or pilot_bundle_manifest_sha256
                    != retained_pressure_receipt["pilot_bundle_manifest_sha256"]
                    or aggregate_pilot_analysis_sha256
                    != retained_pressure_receipt["aggregate_pilot_analysis_sha256"]
                    or (expected_result is not None and expected_result != result)
                ):
                    raise pressure_verifier.PressureReviewPacketVerificationError(
                        "fixture sealed reanalysis attestation binding drifted"
                    )
                if (
                    reject_source_snapshot
                    and expected_result is None
                    and reanalysis_source_snapshot_side_effect is not None
                ):
                    raise reanalysis_source_snapshot_side_effect
                return {
                    "binding": expected_reanalysis_binding,
                    "attestation": {},
                    "operator_id": "fixture.operator",
                    "recomputed_utc": "2026-06-05T09:10:00Z",
                    "sealed_utc": "2026-06-05T09:11:00Z",
                    "evidence": {
                        "published_pressure_pilot_receipt": expected_aggregate,
                        "published_pressure_pilot_review_packet_receipt": (
                            expected_packet
                        ),
                        "pilot_bundle_manifest_sha256": retained_pressure_receipt[
                            "pilot_bundle_manifest_sha256"
                        ],
                        "aggregate_pilot_analysis_sha256": retained_pressure_receipt[
                            "aggregate_pilot_analysis_sha256"
                        ],
                    },
                    "source_authorization": reanalysis_source_authorization,
                    "result": result,
                }

            def consume_reanalysis_attestation(
                attestation_binding: object,
                **kwargs: object,
            ) -> dict[str, object]:
                return reanalysis_verification(
                    attestation_binding,
                    **kwargs,
                    reject_source_snapshot=True,
                )

            def consume_installed_reanalysis_attestation(
                attestation_binding: object,
                **kwargs: object,
            ) -> dict[str, object]:
                return reanalysis_verification(
                    attestation_binding,
                    **kwargs,
                    reject_source_snapshot=False,
                )

            def reviewer_verification(
                attestation_binding: object,
                *,
                aggregate_receipt_binding: object,
                packet_receipt_binding: object,
                reanalysis_verification: object,
                selected_case: object,
                authorized_pic_root: Path,
            ) -> dict[str, object]:
                expected_packet, expected_aggregate = (
                    expected_pressure_review_packet_bindings(authorized_pic_root)
                )
                if (
                    attestation_binding != expected_reviewer_binding
                    or aggregate_receipt_binding != expected_aggregate
                    or packet_receipt_binding != expected_packet
                    or not isinstance(reanalysis_verification, dict)
                    or reanalysis_verification.get("binding")
                    != expected_reanalysis_binding
                    or selected_case != retained_pressure_receipt["selected_case"]
                ):
                    raise pressure_verifier.PressureReviewPacketVerificationError(
                        "fixture sealed reviewer attestation binding drifted"
                    )
                return {
                    "binding": expected_reviewer_binding,
                    "attestation": {},
                    "reviewer_id": "fixture.reviewer",
                    "reviewed_utc": "2026-06-05T09:12:00Z",
                    "sealed_utc": "2026-06-05T09:13:00Z",
                    "rationale": "Fixture-only reviewed selection.",
                    "selected_case": retained_pressure_receipt["selected_case"],
                    "authoritative_reanalysis_attestation": (
                        expected_reanalysis_binding
                    ),
                }

            def consume_installed_pressure_review_packet(
                path: str | Path,
                *,
                aggregate_receipt_binding: object,
                authorized_pic_root: Path,
            ) -> dict[str, object]:
                expected_packet, expected_aggregate = (
                    expected_pressure_review_packet_bindings(authorized_pic_root)
                )
                if (
                    Path(path) != Path(expected_packet["path"])
                    or aggregate_receipt_binding != expected_aggregate
                ):
                    raise ValueError("fixture installed pressure review-packet drifted")
                return {
                    "receipt_binding": expected_packet,
                    "aggregate_receipt_binding": expected_aggregate,
                    "packet_receipt": {},
                    "aggregate_receipt": aggregate_record,
                    "aggregate_bundle": aggregate_record["aggregate_bundle"],
                    "aggregate_analysis": aggregate_record["aggregate_analysis"],
                    "source_bindings": {},
                    "inventory": {},
                }

            def verify_registered_execution_ledger(
                registered_execution_receipt: dict[str, Any],
                *,
                authorized_pic_root: Path,
                **_kwargs: object,
            ) -> dict[str, str]:
                if authorized_pic_root != authorized_root:
                    raise campaign.QualificationError(
                        "fixture registered-execution ledger root drifted",
                        code="registered_execution_receipt_ledger_drift",
                    )
                return {
                    "reservation_id": registered_execution_receipt["reservation_id"],
                    "submission_id": registered_execution_receipt["submission_id"],
                    "reconciliation_event_sha256": registered_execution_receipt[
                        "reconciliation_event_sha256"
                    ],
                }

            @contextmanager
            def verified_registered_execution_ledger_snapshot(
                registered_execution_receipt: dict[str, Any],
                *,
                authorized_pic_root: Path,
                **kwargs: object,
            ) -> Iterator[dict[str, str]]:
                yield verify_registered_execution_ledger(
                    registered_execution_receipt,
                    authorized_pic_root=authorized_pic_root,
                    **kwargs,
                )

            try:
                frozen = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT if receipt is None else receipt,
                    authorized_root=authorized_root,
                )
                with mock.patch.object(
                    campaign, "_validate_external_clean_candidate_closure", **patch_kwargs
                ), mock.patch.object(
                    campaign.pressure_selection.historical_pressure_pilot_consumer,
                    "consume_exact_historical_production_pressure_pilot",
                    side_effect=consume_pressure_pilot_bundle,
                ), mock.patch.object(
                    campaign.pressure_selection.pressure_review_packet_verifier,
                    "consume_published_pressure_pilot_review_packet",
                    side_effect=verify_pressure_review_packet,
                ), mock.patch.object(
                    campaign.pressure_selection.pressure_review_packet_verifier,
                    "consume_sealed_pressure_reanalysis_attestation",
                    side_effect=consume_reanalysis_attestation,
                ), mock.patch.object(
                    campaign.pressure_selection.pressure_review_packet_verifier,
                    "consume_sealed_pressure_reviewer_attestation",
                    side_effect=reviewer_verification,
                ), mock.patch.object(
                    control_plane_common,
                    "consume_published_pressure_pilot_review_packet",
                    side_effect=consume_installed_pressure_review_packet,
                ), mock.patch.object(
                    control_plane_common,
                    "consume_sealed_pressure_reanalysis_attestation",
                    side_effect=consume_installed_reanalysis_attestation,
                ), mock.patch.dict(
                    control_plane_common.validate_pressure_reanalysis_source_snapshot.__globals__,
                    {
                        "consume_sealed_pressure_reanalysis_attestation": (
                            consume_installed_reanalysis_attestation
                        )
                    },
                ), mock.patch.object(
                    control_plane_common,
                    "consume_sealed_pressure_reviewer_attestation",
                    side_effect=reviewer_verification,
                ), mock.patch.object(
                    campaign,
                    "_validated_registered_execution_receipt_ledger_snapshot",
                    side_effect=verified_registered_execution_ledger_snapshot,
                ):
                    yield authorized_root, tree, frozen["inventory_sha256"]
            finally:
                _make_writable_tree(tree)
                _make_writable_tree(authorized_root / "plans")
                minimal_plans = authorized_root / "minimal-plans"
                if minimal_plans.exists():
                    _make_writable_tree(minimal_plans)


def _qualify(
    authorized_root: Path,
    tree: Path,
    inventory_sha256: str,
) -> dict[str, Any]:
    return campaign.qualify_campaign(
        tree,
        inventory_sha256,
        authorized_orion_root=authorized_root,
    )


class Q011Section54CampaignAdmissionTests(unittest.TestCase):
    def assert_rejected(self, result: dict[str, Any], code: str) -> None:
        self.assertFalse(result["admitted_for_follow_on_numerical_qualification"])
        self.assertIsNone(result["admission"])
        self.assertEqual(result["failure_reasons"][0]["code"], code)

    def test_valid_fixture_is_admitted_without_claiming_numerical_closure(self) -> None:
        with _frozen_fixture() as fixture:
            result = _qualify(*fixture)
        self.assertTrue(result["admitted_for_follow_on_numerical_qualification"])
        self.assertFalse(result["final_claim_closure"])
        self.assertEqual(result["failure_reasons"], [])
        admission = result["admission"]
        self.assertEqual(set(admission["endpoint_products"]), {"500.0", "1200.0"})
        self.assertEqual(len(admission["retained_snapshot_products"]), 13)
        self.assertEqual(len(admission["restart_publications"]), 13)
        self.assertEqual(
            admission["particle_endpoints"]["500.0"]["execution_header"]["nranks"], 4
        )
        self.assertEqual(
            admission["preregistration_binding"]["binding_scope"],
            "complete_retained_bytes_equal_invoked_frozen_policy",
        )
        self.assertEqual(
            admission["numerical_qualification_status"],
            "not_evaluated_by_artifact_admission_slice",
        )
        self.assertEqual(
            admission["attempt_identity"]["campaign_plan_sha256"],
            admission["artifact_bindings"]["campaign_plan"]["sha256"],
        )
        self.assertEqual(
            admission["stdout_telemetry"]["pic_runtime_identity"]["state"],
            "momentum_p_over_m",
        )

    def test_observed_committed_time_overshoot_is_admitted_by_nominal_slot(self) -> None:
        observed = 500.053496123

        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_slot_observed_time(root, manifest, 500.0, observed)

        with _frozen_fixture(mutate) as fixture:
            result = _qualify(*fixture)
        self.assertTrue(result["admitted_for_follow_on_numerical_qualification"])
        payload = result["admission"]["snapshot_payloads"]["500.0"]
        self.assertEqual(payload["nominal_slot_time"], 500.0)
        self.assertEqual(payload["observed_committed_time"], observed)
        self.assertEqual(payload["mesh_bins"]["rho"]["time"], 500.053)

    def test_float32_due_cell_early_observed_time_is_admitted(self) -> None:
        observed = 499.99999

        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_slot_observed_time(root, manifest, 500.0, observed)

        with _frozen_fixture(mutate) as fixture:
            result = _qualify(*fixture)
        self.assertTrue(result["admitted_for_follow_on_numerical_qualification"])
        self.assertEqual(
            result["admission"]["snapshot_payloads"]["500.0"][
                "observed_committed_time"
            ],
            observed,
        )

    def test_selected_t500_lateness_cap_is_fail_closed(self) -> None:
        for observed, admitted in ((500.1, True), (500.100001, False)):
            with self.subTest(observed=observed):
                def mutate(
                    root: Path,
                    manifest: dict[str, Any],
                    *,
                    observed: float = observed,
                ) -> None:
                    _rewrite_slot_observed_time(root, manifest, 500.0, observed)

                with _frozen_fixture(mutate) as fixture:
                    result = _qualify(*fixture)
                self.assertEqual(
                    result["admitted_for_follow_on_numerical_qualification"],
                    admitted,
                )
                if not admitted:
                    self.assert_rejected(result, "snapshot_cadence_drift")

    def test_early_finalization_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            stdout = _find_product(manifest, "stdout", None)
            _rewrite_product(
                root,
                stdout,
                _stdout_telemetry().replace(
                    b"Terminating on time limit\ntime=1200 cycle=120\ntlim=1200 nlim=-1\n",
                    b"Terminating on wall clock limit\ntime=1199 cycle=120\ntlim=1200 nlim=-1\n",
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "terminal_completion_drift")

    def test_ledger_snapshot_stays_pinned_during_raw_tree_validation(self) -> None:
        with _frozen_fixture() as fixture:
            authorized_root, _, _ = fixture
            ledger_parent = authorized_root / "ledger-lifetime-probe"
            ledger_parent.mkdir()
            ledger_path = ledger_parent / "node_hours.jsonl"
            hidden_path = ledger_parent / "hidden-node-hours.jsonl"
            ledger_path.write_text("fixture ledger\n", encoding="utf-8")
            mutation_seen = False
            validate_product_hashes = campaign._validate_product_hashes

            @contextmanager
            def verify_registered_execution_ledger(
                registered_execution_receipt: dict[str, Any],
                *,
                authorized_pic_root: Path,
                **_kwargs: object,
            ) -> Iterator[dict[str, str]]:
                yield {
                    "reservation_id": registered_execution_receipt["reservation_id"],
                    "submission_id": registered_execution_receipt["submission_id"],
                    "reconciliation_event_sha256": registered_execution_receipt[
                        "reconciliation_event_sha256"
                    ],
                }
                if mutation_seen:
                    raise campaign.QualificationError(
                        "fixture ledger namespace changed during raw-tree validation",
                        code="registered_execution_receipt_ledger_drift",
                    )

            def hide_and_restore_ledger(*args: object, **kwargs: object) -> None:
                nonlocal mutation_seen
                ledger_path.rename(hidden_path)
                hidden_path.rename(ledger_path)
                mutation_seen = True
                validate_product_hashes(*args, **kwargs)

            with mock.patch.object(
                campaign,
                "_validated_registered_execution_receipt_ledger_snapshot",
                side_effect=verify_registered_execution_ledger,
            ), mock.patch.object(
                campaign,
                "_validate_product_hashes",
                side_effect=hide_and_restore_ledger,
            ):
                result = _qualify(*fixture)
        self.assertTrue(mutation_seen)
        self.assert_rejected(
            result, "registered_execution_receipt_ledger_drift"
        )

    def test_all_preregistered_canonical_variants_are_admitted(self) -> None:
        for variant in (
            "coarse_uniform_dx12",
            "three_level_amr_root_dx12_finest_dx3",
            "fine_uniform_dx3",
        ):
            with self.subTest(variant=variant):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, variant: str = variant
                ) -> None:
                    _retarget_attempt(root, manifest, variant)

                with _frozen_fixture(mutate) as fixture:
                    self.assertTrue(
                        _qualify(*fixture)["admitted_for_follow_on_numerical_qualification"]
                    )

    def test_legacy_variant_aliases_are_rejected(self) -> None:
        for variant in ("amr", "fine_uniform"):
            with self.subTest(variant=variant):
                def mutate(
                    _root: Path, manifest: dict[str, Any], *, variant: str = variant
                ) -> None:
                    manifest["run_identity"]["variant"] = variant

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(_qualify(*fixture), "invalid_run_identity")

    def test_attempt_identity_and_retained_binding_omissions_are_rejected(self) -> None:
        for section, name in (
            ("attempt_identity", "campaign_plan_sha256"),
            ("attempt_identity", "planner_materialization_receipt_sha256"),
            ("artifact_bindings", "campaign_plan"),
            ("artifact_bindings", "planner_materialization_receipt"),
        ):
            with self.subTest(section=section, name=name):
                def mutate(
                    _root: Path,
                    manifest: dict[str, Any],
                    *,
                    section: str = section,
                    name: str = name,
                ) -> None:
                    del manifest[section][name]

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(_qualify(*fixture), "schema_key_error")

    def test_attempt_digest_cross_binding_drift_is_rejected(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["attempt_identity"]["attempt_contract_sha256"] = "f" * 64

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "attempt_binding_drift")

    def test_duplicate_retained_binding_paths_are_rejected(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            plan = manifest["artifact_bindings"]["campaign_plan"]
            contract = manifest["artifact_bindings"]["attempt_contract"]
            contract["path"] = plan["path"]
            contract["sha256"] = plan["sha256"]
            manifest["attempt_identity"]["attempt_contract_sha256"] = plan["sha256"]

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "duplicate_declared_path")

    def test_model_override_order_and_duplicates_are_rejected(self) -> None:
        fine = campaign.frozen_model.variant_binding("fine_uniform_dx3")
        cases = {
            "order": list(reversed(fine.model_launch_overrides)),
            "duplicate": [fine.model_launch_overrides[0], fine.model_launch_overrides[0]],
        }
        for label, overrides in cases.items():
            with self.subTest(label=label):
                def mutate(
                    _root: Path,
                    manifest: dict[str, Any],
                    *,
                    overrides: list[str] = overrides,
                ) -> None:
                    manifest["run_identity"]["variant"] = "fine_uniform_dx3"
                    manifest["attempt_identity"]["model_launch_overrides"] = overrides

                with _frozen_fixture(mutate) as fixture:
                    expected = (
                        "duplicate_attempt_binding"
                        if label == "duplicate"
                        else "invalid_attempt_identity"
                    )
                    self.assert_rejected(_qualify(*fixture), expected)

    def test_seed_override_alias_and_drift_are_rejected(self) -> None:
        for value, expected in ((True, "schema_type_error"), (23050102, "invalid_attempt_identity")):
            with self.subTest(value=value):
                def mutate(
                    _root: Path, manifest: dict[str, Any], *, value: object = value
                ) -> None:
                    manifest["attempt_identity"]["seed_overrides"][
                        "problem/ps_inject_seed"
                    ] = value

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(_qualify(*fixture), expected)

    def test_stdout_runtime_identity_drift_is_rejected(self) -> None:
        complete = _runtime_identity_line()
        for payload in (
            _stdout_telemetry(include_runtime_identity=False),
            _stdout_telemetry(runtime_replacements={"state": "velocity"}),
            _stdout_telemetry().replace(b" induction=ideal_mhd_only", b""),
            _stdout_telemetry().replace(
                complete.encode("ascii"),
                (complete + " unexpected_projection=on").encode("ascii"),
            ),
            _stdout_telemetry().replace(
                b" background=coupled feedback=coupled",
                b" feedback=coupled background=coupled",
            ),
        ):
            with self.subTest(payload=payload):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, payload: bytes = payload
                ) -> None:
                    _rewrite_product(root, _find_product(manifest, "stdout", None), payload)

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(_qualify(*fixture), "invalid_runtime_identity")

    def test_semantic_campaign_plan_crosslink_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_attempt_json_artifact(
                root,
                manifest,
                "campaign_plan",
                lambda plan: plan["source_bindings"]["paper_deck"].__setitem__(
                    "sha256", "0" * 64
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(
                _qualify(*fixture), "planner_materialization_receipt_drift"
            )

    def test_self_authored_campaign_plan_graph_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_attempt_json_artifact(
                root,
                manifest,
                "campaign_plan",
                lambda plan: plan["nonauthorizing_policy_fragment"].__setitem__(
                    "sha256", "0" * 64
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(
                _qualify(*fixture), "planner_materialization_receipt_drift"
            )

    def test_minimal_planner_tree_with_dangling_children_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            campaign_plan = manifest["artifact_bindings"]["campaign_plan"]
            helper_closure = manifest["artifact_bindings"][
                "analyzer_helper_source_closure_manifest"
            ]
            receipt = _minimal_planner_materialization_receipt(
                root,
                (root / campaign_plan["path"]).read_bytes(),
                (root / helper_closure["path"]).read_bytes(),
            )
            manifest["artifact_bindings"]["planner_materialization_receipt"] = receipt
            manifest["attempt_identity"][
                "planner_materialization_receipt_sha256"
            ] = receipt["sha256"]

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(
                _qualify(*fixture), "planner_materialization_receipt_drift"
            )

    def test_planner_environment_profile_rebinding_is_rejected(self) -> None:
        def mutate(
            plan: dict[str, Any], source_payloads: dict[str, bytes]
        ) -> None:
            payload = b"self-authored environment profile\n"
            source_payloads["environment_profile"] = payload
            plan["source_bindings"]["environment_profile"]["sha256"] = _sha256(payload)

        with _frozen_fixture(planner_mutate=mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "campaign_plan_crosslink_drift")

    def test_planner_restart_preregistration_rebinding_is_rejected(self) -> None:
        def mutate(
            plan: dict[str, Any], source_payloads: dict[str, bytes]
        ) -> None:
            payload = source_payloads["restart_preregistration"] + b"\n"
            source_payloads["restart_preregistration"] = payload
            plan["source_bindings"]["restart_preregistration"]["sha256"] = _sha256(
                payload
            )

        with _frozen_fixture(planner_mutate=mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "campaign_plan_crosslink_drift")

    def test_semantic_attempt_contract_crosslink_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_attempt_json_artifact(
                root,
                manifest,
                "attempt_contract",
                lambda contract: contract.__setitem__("qualifying_seed", 23050102),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "attempt_contract_crosslink_drift")

    def test_wrong_retained_root_is_rejected(self) -> None:
        with _frozen_fixture(wrong_retained_root=True) as fixture:
            self.assert_rejected(_qualify(*fixture), "retained_root_binding_drift")

    def test_registered_execution_receipt_raw_root_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_attempt_json_artifact(
                root,
                manifest,
                "registered_execution_receipt",
                lambda receipt: receipt.__setitem__(
                    "raw_output_root", "/fixture/arbitrarily-relabeled-raw"
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(
                _qualify(*fixture), "registered_execution_receipt_drift"
            )

    def test_semantic_selected_pressure_receipt_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_attempt_json_artifact(
                root,
                manifest,
                "selected_pressure_receipt",
                lambda receipt: receipt["selected_case"].__setitem__(
                    "problem_ps_p0", 0.2
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "selected_pressure_receipt_drift")

    def test_pressure_reanalysis_source_snapshot_rejection_is_rejected(self) -> None:
        error_type = (
            campaign.pressure_selection.pressure_review_packet_verifier
            .PressureReviewPacketVerificationError
        )
        with _frozen_fixture(
            reanalysis_source_snapshot_side_effect=error_type(
                "fixture source-snapshot rejection"
            )
        ) as fixture:
            self.assert_rejected(_qualify(*fixture), "selected_pressure_receipt_drift")

    def test_selected_pressure_receipt_resource_limits_fail_closed(self) -> None:
        cases = (
            (
                "oversized",
                b" "
                * (
                    campaign.pressure_selection.MAX_PRESSURE_SELECTION_RECEIPT_BYTES
                    + 1
                ),
            ),
            (
                "deeply_nested",
                b'{"nested":'
                + (b"[" * 2000)
                + b"0"
                + (b"]" * 2000)
                + b"}\n",
            ),
        )
        for label, payload in cases:
            with self.subTest(label=label):
                def mutate(root: Path, manifest: dict[str, Any]) -> None:
                    _rewrite_bound_artifact(
                        root,
                        manifest,
                        "selected_pressure_receipt",
                        payload,
                    )
                    manifest["attempt_identity"][
                        "selected_pressure_receipt_sha256"
                    ] = manifest["artifact_bindings"]["selected_pressure_receipt"][
                        "sha256"
                    ]

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(
                        _qualify(*fixture),
                        "selected_pressure_receipt_drift",
                    )

    def test_fake_pressure_and_attestation_bindings_are_rejected(self) -> None:
        cases = (
            ("published_pressure_pilot_receipt", "path", "/fixture/fake_receipt.json"),
            ("published_pressure_pilot_receipt", "sha256", "0" * 64),
            (
                "published_pressure_pilot_review_packet_receipt",
                "path",
                "/fixture/fake_review_packet_receipt.json",
            ),
            ("published_pressure_pilot_review_packet_receipt", "sha256", "3" * 64),
            (
                "authoritative_reanalysis_attestation",
                "path",
                "/fixture/fake_reanalysis_attestation.json",
            ),
            ("authoritative_reanalysis_attestation", "sha256", "4" * 64),
            (
                "reviewer_attestation",
                "path",
                "/fixture/fake_reviewer_attestation.json",
            ),
            ("reviewer_attestation", "sha256", "5" * 64),
            (None, "pilot_bundle_manifest_sha256", "1" * 64),
            (None, "aggregate_pilot_analysis_sha256", "2" * 64),
        )
        for section, field, value in cases:
            with self.subTest(section=section, field=field):
                def mutate(
                    root: Path,
                    manifest: dict[str, Any],
                    *,
                    section: str | None = section,
                    field: str = field,
                    value: str = value,
                ) -> None:
                    def rewrite(receipt: dict[str, Any]) -> None:
                        target = receipt if section is None else receipt[section]
                        target[field] = value

                    _rewrite_attempt_json_artifact(
                        root,
                        manifest,
                        "selected_pressure_receipt",
                        rewrite,
                    )

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(
                        _qualify(*fixture), "selected_pressure_receipt_drift"
                    )

    def test_helper_source_closure_requires_model_and_snapshot_members(self) -> None:
        for path in (
            "tst/publication/q011_section54_model.py",
            "tst/publication/q011_section54_historical_pressure_pilot_consumer.py",
            "tst/publication/frontier_control_plane/q011_pressure_review_packet_verifier.py",
            "tst/publication/frontier_control_plane/operator_attestation.py",
        ):
            with self.subTest(path=path):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, path: str = path
                ) -> None:
                    def remove_source(closure: dict[str, Any]) -> None:
                        closure["sources"] = [
                            source
                            for source in closure["sources"]
                            if source["path"] != path
                        ]

                    _rewrite_attempt_json_artifact(
                        root,
                        manifest,
                        "analyzer_helper_source_closure_manifest",
                        remove_source,
                    )

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(
                        _qualify(*fixture), "helper_source_closure_drift"
                    )

    def test_path_escape_is_rejected_without_partial_success(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["artifact_bindings"]["deck"]["path"] = "../outside.athinput"

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "unsafe_relative_path")

    def test_fictional_text_executable_is_rejected_even_when_self_consistent(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_bound_artifact(root, manifest, "executable", b"fictional executable\n")
            digest = manifest["artifact_bindings"]["executable"]["sha256"]
            _rewrite_clean_candidate(
                root,
                manifest,
                lambda candidate: candidate["build"].__setitem__(
                    "executable_sha256", digest
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_executable_elf")

    def test_external_clean_candidate_closure_failure_is_rejected(self) -> None:
        error = campaign.QualificationError(
            "fixture external candidate closure failed",
            code="external_clean_candidate_closure_error",
        )
        with _frozen_fixture(external_side_effect=error) as fixture:
            self.assert_rejected(
                _qualify(*fixture), "external_clean_candidate_closure_error"
            )

    def test_fictional_plain_mesh_bin_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_product(root, _find_product(manifest, "rho", 500.0), b"rho\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_mesh_bin")

    def test_valid_mesh_bin_with_wrong_observable_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "rho", 500.0)
            _rewrite_product(root, product, _mesh_bin("bmag", 500.0, 5))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_mesh_bin")

    def test_duplicate_ptag_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "prtcl_all", 500.0)
            _rewrite_product(root, product, _particle_vtk(500.0, 5, ptags=(10, 10)))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_prtcl_all")

    def test_pvtk_embedded_time_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "prtcl_all", 500.0)
            _rewrite_product(root, product, _particle_vtk(999.0, 5))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "snapshot_metadata_drift")

    def test_pvtk_manifest_time_must_directly_match_embedded_time(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "prtcl_all", 500.0)
            product["observed_committed_time"] = 500.00000075
            _rewrite_product(root, product, _particle_vtk(499.99999925, 5))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "snapshot_metadata_drift")

    def test_wrong_freeze_receipt_semantics_are_rejected(self) -> None:
        receipt = dict(_RECEIPT)
        receipt["qualification_effect"] = "fictional_effect"
        with _frozen_fixture(receipt=receipt) as fixture:
            self.assert_rejected(_qualify(*fixture), "freeze_receipt_semantics_drift")

    def test_cadence_omission_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "j2", 100.0)
            manifest["products"].remove(product)
            (root / product["path"]).unlink()

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "missing_snapshot_cadence")

    def test_duplicate_endpoint_is_rejected_as_ambiguous_cadence(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            manifest["products"].append(
                _product(
                    root,
                    "rho",
                    "bin/q011.rho.00500.duplicate.bin",
                    _mesh_bin("rho", 500.0, 5),
                    500.0,
                )
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "ambiguous_snapshot_cadence")

    def test_fictional_plain_restart_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_restart_payload(root, manifest, 500.0, b"restart\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_restart_payload")

    def test_restart_marker_corruption_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            marker = _find_product(manifest, "restart_complete", 500.0)
            _rewrite_product(root, marker, b"restart complete\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "restart_marker_corruption")

    def test_fictional_plain_stdout_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_product(root, _find_product(manifest, "stdout", None), b"completed\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_q017_telemetry")

    def test_stdout_telemetry_omission_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_product(
                root,
                _find_product(manifest, "stdout", None),
                _stdout_telemetry(omit=frozenset({"particles.total"})),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_q017_telemetry")

    def test_complete_preregistration_checksum_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            payload = campaign.PREREGISTRATION_PATH.read_bytes() + b" "
            _rewrite_bound_artifact(root, manifest, "preregistration", payload)

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "hash_drift")

    def test_forged_outer_candidate_identity_is_rejected(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["candidate_binding"]["git_commit"] = "9" * 40

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "clean_candidate_binding_drift")

    def test_forged_self_consistent_deck_or_analyzer_is_rejected(self) -> None:
        for name in ("deck", "analyzer"):
            with self.subTest(name=name):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, name: str = name
                ) -> None:
                    _rewrite_bound_artifact(
                        root, manifest, name, f"forged {name}\n".encode("ascii")
                    )

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(
                        _qualify(*fixture), "clean_candidate_binding_drift"
                    )

    def test_clean_candidate_schema_aliases_are_rejected(self) -> None:
        for value in (True, 4.0, "4"):
            with self.subTest(value=value):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, value: object = value
                ) -> None:
                    _rewrite_clean_candidate(
                        root,
                        manifest,
                        lambda candidate: candidate.__setitem__("schema_version", value),
                    )

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(_qualify(*fixture), "schema_type_error")

    def test_duplicate_prepared_path_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            def duplicate(candidate: dict[str, Any]) -> None:
                decks = candidate["prepared_artifacts"]["paper_decks"]
                decks.append(copy.deepcopy(decks[0]))

            _rewrite_clean_candidate(root, manifest, duplicate)

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "clean_candidate_schema_error")

    def test_boolean_seed_alias_is_rejected(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["run_identity"]["seed"] = True

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "schema_type_error")

    def test_weighted_spectrum_utility_uses_preregistered_fixed_bins(self) -> None:
        policy = json.loads(campaign.PREREGISTRATION_PATH.read_text(encoding="utf-8"))
        report = campaign.weighted_spectrum_from_chi(
            [0.5, 1.0, 2.0, 2048.0],
            [3.0, 4.0, 5.0, 6.0],
            policy,
        )
        self.assertEqual(report["underflow_count"], 1)
        self.assertEqual(report["overflow_count"], 1)
        self.assertEqual(report["underflow_weight"], 3.0)
        self.assertEqual(report["overflow_weight"], 6.0)
        self.assertEqual(report["macro_weight_in_bins"], 9.0)
        self.assertEqual(report["admitted_weight"], 18.0)
        self.assertEqual(report["overflow_macro_weight_fraction"], 1.0 / 3.0)
        self.assertFalse(report["overflow_gate"]["passed"])


if __name__ == "__main__":
    unittest.main()
