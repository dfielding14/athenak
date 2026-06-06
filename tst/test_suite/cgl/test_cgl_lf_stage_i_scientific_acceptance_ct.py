"""Adversarial native-restart CT-divB tests for Stage I scientific acceptance."""

from __future__ import annotations

from array import array
from copy import deepcopy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import struct
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
INVENTORY_BUILDER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_ct_inventory.py"


def load_utility():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_scientific_acceptance_ct", UTILITY
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


acceptance = load_utility()


def load_inventory_builder():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_ct_inventory", INVENTORY_BUILDER
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


inventory_builder = load_inventory_builder()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def doubles(values: list[float]) -> bytes:
    result = array("d", values)
    if sys.byteorder != "little":
        result.byteswap()
    return result.tobytes()


def flat(k: int, j: int, i: int, nj: int, ni: int) -> int:
    return (k * nj + j) * ni + i


def region_indices(
    ng: int,
    nx1: int,
    nx2: int,
    nx3: int,
    *,
    coarse: bool,
    invalid_active: bool = False,
) -> list[int]:
    active = [
        ng,
        ng + nx1 - 1 + int(invalid_active),
        ng,
        ng + nx2 - 1,
        ng,
        ng + nx3 - 1,
    ]
    if not coarse:
        return [ng, nx1, nx2, nx3, *active, *([0] * 9)]
    cnx1, cnx2, cnx3 = nx1 // 2, nx2 // 2, nx3 // 2
    return [
        ng,
        nx1,
        nx2,
        nx3,
        *active,
        cnx1,
        cnx2,
        cnx3,
        ng,
        ng + cnx1 - 1,
        ng,
        ng + cnx2 - 1,
        ng,
        ng + cnx3 - 1,
    ]


def expected_local_blocks(nmb_total: int, rank_count: int, rank: int) -> int:
    base, remainder = divmod(nmb_total, rank_count)
    return base + int(rank >= rank_count - remainder and remainder > 0)


def write_restart(
    path: Path,
    *,
    rank: int = 0,
    rank_count: int = 2,
    nmb_total: int = 5,
    time_value: float = 9.0,
    divergent: bool = False,
    copied_rank_bytes: bool = False,
    reverse_locations: bool = False,
    duplicate_locations: bool = False,
    invalid_active: bool = False,
    one_dimensional: bool = False,
    bfloor: float = 1.0e-10,
    amr: bool = False,
) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    ng = 1
    block_nx1 = 2
    block_nx2 = 1 if one_dimensional else 2
    block_nx3 = 1 if one_dimensional else 2
    mesh_nx1 = block_nx1 * nmb_total
    mesh_nx2 = block_nx2
    mesh_nx3 = block_nx3
    parameters = (
        "<mesh>\n"
        f"nghost = {ng}\n"
        f"nx1 = {mesh_nx1}\n"
        f"nx2 = {mesh_nx2}\n"
        f"nx3 = {mesh_nx3}\n"
        "<meshblock>\n"
        f"nx1 = {block_nx1}\n"
        f"nx2 = {block_nx2}\n"
        f"nx3 = {block_nx3}\n"
        "<mhd>\n"
        "eos = cgl\n"
        "nscalars = 0\n"
        f"bfloor = {bfloor:.17g}\n"
        "<output3>\n"
        "single_file_per_rank = true\n"
        "<time>\n"
        f"restart_time = {time_value:.17g}\n"
        "<par_end>\n"
    ).encode()
    dx = 0.5
    region = [
        0.0,
        0.0,
        0.0,
        mesh_nx1 * dx,
        mesh_nx2 * dx,
        mesh_nx3 * dx,
        dx,
        dx,
        dx,
    ]
    mesh_indices = region_indices(
        ng, mesh_nx1, mesh_nx2, mesh_nx3, coarse=False
    )
    block_indices = region_indices(
        ng,
        block_nx1,
        block_nx2,
        block_nx3,
        coarse=True,
        invalid_active=invalid_active,
    )
    root_level = (nmb_total - 1).bit_length()
    header = struct.pack(
        acceptance.RESTART_HEADER_FORMAT,
        nmb_total,
        root_level,
        *region,
        *mesh_indices,
        *block_indices,
        time_value,
        0.01,
        100,
    )
    locations = [
        (index, 0, 0, root_level + int(amr)) for index in range(nmb_total)
    ]
    if reverse_locations:
        locations.reverse()
    if duplicate_locations:
        locations[-1] = locations[0]
    meshblock_list = (
        b"".join(struct.pack("<4i", *location) for location in locations)
        + b"".join(struct.pack("<f", 1.0) for _ in locations)
    )

    nout1 = block_nx1 + 2 * ng
    nout2 = block_nx2 + 2 * ng
    nout3 = block_nx3 + 2 * ng
    cells = nout1 * nout2 * nout3
    x1_count = (nout1 + 1) * nout2 * nout3
    x2_count = nout1 * (nout2 + 1) * nout3
    x3_count = nout1 * nout2 * (nout3 + 1)
    cell_centered = [0.0] * (6 * cells)
    x1f = [1.0] * x1_count
    x2f = [0.0] * x2_count
    x3f = [0.0] * x3_count
    if divergent:
        x1f[flat(ng, ng, ng + 1, nout2, nout1 + 1)] = 2.0
    block_data = doubles(cell_centered + x1f + x2f + x3f)
    local_blocks = expected_local_blocks(nmb_total, rank_count, rank)
    rank_token = 0 if copied_rank_bytes else rank
    internal_state = f"qualified-synthetic-rank-{rank_token:08d}".encode()
    payload = (
        parameters
        + header
        + meshblock_list
        + internal_state
        + struct.pack("<Q", len(block_data))
        + block_data * local_blocks
    )
    path.write_bytes(payload)
    return sha256(path)


@pytest.fixture
def policy():
    value = acceptance.load_validated_policy(
        acceptance.DEFAULT_CRITERIA, acceptance.DEFAULT_CRITERIA_REVIEW
    )
    value["replay_tools_approved"] = True
    return value


def build_inventory(
    policy,
    tmp_path: Path,
    *,
    nmb_total: int = 5,
    rank_count: int = 2,
    copied_rank_bytes: bool = False,
    divergent_state_rank: tuple[float, int] | None = None,
    reverse_rank: int | None = None,
    bfloor_by_rank: dict[int, float] | None = None,
) -> tuple[Path, dict[float, list[Path]]]:
    executable = "a" * 64
    case = acceptance.case_manifest_record(policy, "R02")
    segment_paths: list[Path] = []
    state_records: list[dict[str, object]] = []
    restart_paths: dict[float, list[Path]] = {}
    for time_value in (9.0, 10.0):
        ranks: list[dict[str, object]] = []
        paths: list[Path] = []
        for rank in range(rank_count):
            path = (
                tmp_path
                / f"state-{time_value:g}"
                / "rst"
                / f"rank_{rank:08d}"
                / "state.rst"
            )
            write_restart(
                path,
                rank=rank,
                rank_count=rank_count,
                nmb_total=nmb_total,
                time_value=time_value,
                copied_rank_bytes=copied_rank_bytes,
                divergent=divergent_state_rank == (time_value, rank),
                reverse_locations=reverse_rank == rank,
                bfloor=(bfloor_by_rank or {}).get(rank, 1.0e-10),
            )
            paths.append(path)
            ranks.append({**binding(path), "rank": rank})
        restart_paths[time_value] = paths
        inspected_ranks = [
            {key: value for key, value in rank_record.items() if key != "rank"}
            for rank_record in ranks
        ]
        terminal = {
            **inspected_ranks[0],
            "storage": "per_rank",
            "rank_files": inspected_ranks,
        }
        segment_path = tmp_path / f"state-{time_value:g}" / "manifest.json"
        segment = {
            "accounting": {
                "result": "accepted",
                "case_id": "R02",
                "case_name": case["name"],
                "segment": f"state-{time_value:g}",
                "executable_sha256": executable,
            },
            "allocation": {"nodes": 1, "ranks_per_node": rank_count},
            "command": {
                "executable_sha256": executable,
                "matrix_sha256": policy["verified_sources"]["stage_i_manifest"]["sha256"],
            },
            "scientific_inspection": {
                "schema_version": 3,
                "accepted": True,
                "case_id": "R02",
                "segment": f"state-{time_value:g}",
                "manifest": str(segment_path.resolve()),
                "final_time": time_value,
                "required_time": time_value,
                "terminal_restart_time": time_value,
                "restart_times": [time_value],
                "checks": {
                    "required_time_reached": True,
                    "restart_retained": True,
                    "terminal_restart_physical_time_matches_final": True,
                },
                "terminal_restart": terminal,
                "restarts": [terminal],
            },
        }
        write_json(segment_path, segment)
        segment_paths.append(segment_path)
        state_records.append({
            "time": time_value,
            "accepted_segment_manifest": binding(segment_path),
            "rank_files": ranks,
        })
    bundle_path = tmp_path / "bundle.json"
    write_json(bundle_path, {
        "workflow": "paper-mks24-stage-i-production",
        "status": "accepted_for_analysis",
        "production_case_id": "R02",
        "required_final_time": 10.0,
        "accepted_final_time": 10.0,
        "cases": [{
            "name": case["name"],
            "input": case["input"],
            "status": "passed",
        }],
        "production_segment_manifests": [
            str(path.resolve()) for path in segment_paths
        ],
    })
    inventory_path = tmp_path / "inventory.json"
    write_json(inventory_path, {
        "schema_version": 2,
        "record_type": "stage-i-restart-ct-inventory",
        "case_id": "R02",
        "coverage_complete": True,
        "required_state_times": [9.0, 10.0],
        "executable_sha256": executable,
        "accepted_bundle_manifest": binding(bundle_path),
        "states": state_records,
    })
    return inventory_path, restart_paths


def test_offline_multistate_inventory_tests_mechanics_but_is_non_authorizing(
    policy, tmp_path
):
    inventory_path, _ = build_inventory(policy, tmp_path, nmb_total=5, rank_count=2)
    evidence = acceptance.audit_ct_divb(policy, "R02", [], inventory_path)
    acceptance.verify_evidence_digest(evidence, "CT evidence")
    assert evidence["result"] == "inconclusive"
    assert evidence["numerical_result"] == "pass"
    assert evidence["campaign_authority_eligible"] is False
    assert evidence["canonical_campaign_context"] is None
    assert evidence["coverage_complete"] is True
    assert evidence["required_state_times"] == [9, 10]
    assert evidence["state_count"] == 2
    assert evidence["restart_count"] == 4
    assert evidence["sampled_meshblocks"] == 10
    assert evidence["state_meshblocks"] == 10
    assert evidence["maximum_normalized_ct_divb"] == 0.0
    assert [
        record["local_meshblocks"]
        for record in evidence["state_audits"][0]["restart_audits"]
    ] == [2, 3]


def test_pending_independent_review_blocks_ct_inventory_replay(tmp_path):
    pending = acceptance.load_validated_policy(
        acceptance.DEFAULT_CRITERIA, acceptance.DEFAULT_CRITERIA_REVIEW
    )
    inventory_path, _ = build_inventory(pending, tmp_path, nmb_total=5, rank_count=2)
    with pytest.raises(acceptance.AcceptanceError, match="pending independent replay review"):
        acceptance.validate_ct_inventory(pending, "R02", inventory_path)


def test_builder_deterministically_reconstructs_inventory_accepted_by_auditor(
    policy, tmp_path
):
    expected_path, _ = build_inventory(policy, tmp_path / "accepted")
    expected = json.loads(expected_path.read_text())
    bundle_path = Path(expected["accepted_bundle_manifest"]["path"])
    observed = inventory_builder.build_ct_inventory(
        policy, "R02", bundle_path, sha256(bundle_path), "offline"
    )
    assert observed == expected

    observed_path = tmp_path / "observed-inventory.json"
    write_json(observed_path, observed)
    _, _, states, _, canonical = acceptance.validate_ct_inventory(
        policy, "R02", observed_path
    )
    assert [state["time"] for state in states] == [9.0, 10.0]
    assert canonical is None
    evidence = acceptance.audit_ct_divb(policy, "R02", [], observed_path)
    assert evidence["coverage_complete"] is True
    assert evidence["numerical_result"] == "pass"


def test_builder_cli_emits_only_accepted_bundle_derived_inventory(policy, tmp_path):
    expected_path, _ = build_inventory(policy, tmp_path / "accepted")
    expected = json.loads(expected_path.read_text())
    bundle_path = Path(expected["accepted_bundle_manifest"]["path"])
    output = tmp_path / "cli-inventory.json"
    assert inventory_builder.main([
        "--case-id", "R02",
        "--bundle-manifest", str(bundle_path),
        "--expected-bundle-sha256", sha256(bundle_path),
        "--authority-mode", "offline",
        "--output", str(output),
    ]) == 0
    assert json.loads(output.read_text()) == expected
    assert inventory_builder.main([
        "--case-id", "R02",
        "--bundle-manifest", str(bundle_path),
        "--expected-bundle-sha256", sha256(bundle_path),
        "--authority-mode", "offline",
        "--output", str(output),
    ]) == 2


def test_builder_candidate_write_rejects_canonical_root(tmp_path, monkeypatch):
    monkeypatch.setattr(inventory_builder.acceptance, "CANONICAL_CAMPAIGN_ROOT", tmp_path)
    with pytest.raises(inventory_builder.CtInventoryError, match="canonical root"):
        inventory_builder.write_candidate(
            tmp_path / "forbidden-inventory.json",
            {"schema_version": 2},
        )


def test_builder_rejects_nonaccepted_or_incomplete_bundle_lineage(policy, tmp_path):
    expected_path, _ = build_inventory(policy, tmp_path / "nonaccepted")
    expected = json.loads(expected_path.read_text())
    bundle_path = Path(expected["accepted_bundle_manifest"]["path"])
    segment_path = Path(expected["states"][0]["accepted_segment_manifest"]["path"])
    segment = json.loads(segment_path.read_text())
    segment["accounting"]["result"] = "clean_partial"
    write_json(segment_path, segment)
    with pytest.raises(
        inventory_builder.CtInventoryError, match="not an accepted CT source"
    ):
        inventory_builder.build_ct_inventory(
            policy, "R02", bundle_path, sha256(bundle_path), "offline"
        )

    expected_path, _ = build_inventory(policy, tmp_path / "incomplete")
    expected = json.loads(expected_path.read_text())
    bundle_path = Path(expected["accepted_bundle_manifest"]["path"])
    bundle = json.loads(bundle_path.read_text())
    bundle["production_segment_manifests"] = bundle["production_segment_manifests"][1:]
    write_json(bundle_path, bundle)
    with pytest.raises(inventory_builder.CtInventoryError, match="exactly one t=9 state"):
        inventory_builder.build_ct_inventory(
            policy, "R02", bundle_path, sha256(bundle_path), "offline"
        )


def test_builder_rejects_copied_rank_bytes_and_false_canonical_authority(policy, tmp_path):
    expected_path, _ = build_inventory(
        policy, tmp_path / "copied", nmb_total=4, copied_rank_bytes=True
    )
    expected = json.loads(expected_path.read_text())
    bundle_path = Path(expected["accepted_bundle_manifest"]["path"])
    with pytest.raises(inventory_builder.CtInventoryError, match="copied rank bytes"):
        inventory_builder.build_ct_inventory(
            policy, "R02", bundle_path, sha256(bundle_path), "offline"
        )

    expected_path, _ = build_inventory(policy, tmp_path / "noncanonical")
    expected = json.loads(expected_path.read_text())
    bundle_path = Path(expected["accepted_bundle_manifest"]["path"])
    with pytest.raises(inventory_builder.CtInventoryError, match="exact canonical bundle path"):
        inventory_builder.build_ct_inventory(
            policy, "R02", bundle_path, sha256(bundle_path), "canonical"
        )


def test_one_terminal_state_never_false_passes(policy, tmp_path):
    path = tmp_path / "state-9" / "rank_00000000" / "state.rst"
    digest = write_restart(path, rank_count=1, nmb_total=1, time_value=9.0)
    evidence = acceptance.audit_ct_divb(policy, "R02", [f"{path}={digest}"], None)
    assert evidence["result"] == "inconclusive"
    assert evidence["coverage_complete"] is False
    assert evidence["state_count"] == 1

    inventory_path, _ = build_inventory(policy, tmp_path / "bad-inventory")
    inventory = json.loads(inventory_path.read_text())
    inventory["required_state_times"] = [9.0]
    inventory["states"] = inventory["states"][:1]
    write_json(inventory_path, inventory)
    with pytest.raises(acceptance.AcceptanceError, match="required state times differ"):
        acceptance.audit_ct_divb(policy, "R02", [], inventory_path)


def test_divergent_required_state_fails(policy, tmp_path):
    inventory_path, _ = build_inventory(
        policy, tmp_path, divergent_state_rank=(10.0, 1)
    )
    evidence = acceptance.audit_ct_divb(policy, "R02", [], inventory_path)
    assert evidence["result"] == "fail"
    assert evidence["maximum_normalized_ct_divb"] > 1.0e-12


def test_inventory_rejects_partial_or_noncanonical_rank_set(policy, tmp_path):
    inventory_path, _ = build_inventory(policy, tmp_path)
    inventory = json.loads(inventory_path.read_text())
    inventory["states"][0]["rank_files"].pop()
    write_json(inventory_path, inventory)
    with pytest.raises(acceptance.AcceptanceError, match="rank count differs"):
        acceptance.audit_ct_divb(policy, "R02", [], inventory_path)

    inventory_path, _ = build_inventory(policy, tmp_path / "wrong-rank")
    inventory = json.loads(inventory_path.read_text())
    inventory["states"][0]["rank_files"][0]["rank"] = 1
    write_json(inventory_path, inventory)
    with pytest.raises(acceptance.AcceptanceError, match="canonical path"):
        acceptance.audit_ct_divb(policy, "R02", [], inventory_path)


def test_inventory_rejects_copied_duplicate_rank_bytes(policy, tmp_path):
    inventory_path, _ = build_inventory(
        policy, tmp_path, nmb_total=4, copied_rank_bytes=True
    )
    with pytest.raises(acceptance.AcceptanceError, match="copied rank bytes"):
        acceptance.audit_ct_divb(policy, "R02", [], inventory_path)


def test_sibling_contract_and_global_location_consistency_are_exact(policy, tmp_path):
    inventory_path, _ = build_inventory(
        policy, tmp_path / "bfloor", bfloor_by_rank={1: 2.0e-10}
    )
    with pytest.raises(acceptance.AcceptanceError, match="ABI/geometry contract"):
        acceptance.audit_ct_divb(policy, "R02", [], inventory_path)

    inventory_path, _ = build_inventory(
        policy, tmp_path / "locations", reverse_rank=1
    )
    with pytest.raises(acceptance.AcceptanceError, match="global logical locations"):
        acceptance.audit_ct_divb(policy, "R02", [], inventory_path)


def test_parser_rejects_duplicate_locations_invalid_ranges_and_data_size(
    policy, tmp_path
):
    duplicate = tmp_path / "duplicate" / "rank_00000000" / "state.rst"
    digest = write_restart(
        duplicate, rank_count=1, nmb_total=2, duplicate_locations=True
    )
    with pytest.raises(acceptance.AcceptanceError, match="logical locations are duplicated"):
        acceptance.audit_ct_divb(policy, "R02", [f"{duplicate}={digest}"], None)

    invalid = tmp_path / "invalid-active" / "rank_00000000" / "state.rst"
    digest = write_restart(invalid, rank_count=1, nmb_total=1, invalid_active=True)
    with pytest.raises(acceptance.AcceptanceError, match="active ranges"):
        acceptance.audit_ct_divb(policy, "R02", [f"{invalid}={digest}"], None)

    truncated = tmp_path / "truncated" / "rank_00000000" / "state.rst"
    write_restart(truncated, rank_count=1, nmb_total=1)
    truncated.write_bytes(truncated.read_bytes()[:-8])
    digest = sha256(truncated)
    with pytest.raises(acceptance.AcceptanceError, match="variable-data boundary"):
        acceptance.audit_ct_divb(policy, "R02", [f"{truncated}={digest}"], None)


def test_parser_explicitly_rejects_non_3d_false_pass_case(policy, tmp_path):
    path = tmp_path / "one-dimensional" / "rank_00000000" / "state.rst"
    digest = write_restart(
        path, rank_count=1, nmb_total=1, one_dimensional=True
    )
    with pytest.raises(acceptance.AcceptanceError, match="three-dimensional state"):
        acceptance.audit_ct_divb(policy, "R02", [f"{path}={digest}"], None)


def test_parser_rejects_amr_state(policy, tmp_path):
    path = tmp_path / "amr" / "rank_00000000" / "state.rst"
    digest = write_restart(path, rank_count=1, nmb_total=1, amr=True)
    with pytest.raises(acceptance.AcceptanceError, match="rejects AMR"):
        acceptance.audit_ct_divb(policy, "R02", [f"{path}={digest}"], None)


def test_hash_read_hash_catches_equal_size_restored_mtime_mutation(
    policy, tmp_path, monkeypatch
):
    path = tmp_path / "mutation" / "rank_00000000" / "state.rst"
    digest = write_restart(path, rank_count=1, nmb_total=1)
    original = acceptance.descriptor_sha256
    profile = path.stat()
    calls = 0

    def mutate_after_initial_hash(descriptor: int) -> str:
        nonlocal calls
        observed = original(descriptor)
        calls += 1
        if calls == 1:
            payload = path.read_bytes()
            marker = b"qualified-synthetic-rank-00000000"
            offset = payload.index(marker)
            with path.open("r+b") as stream:
                stream.seek(offset)
                stream.write(b"Q" + marker[1:])
                stream.flush()
                os.fsync(stream.fileno())
            os.utime(path, ns=(profile.st_atime_ns, profile.st_mtime_ns))
        return observed

    monkeypatch.setattr(acceptance, "descriptor_sha256", mutate_after_initial_hash)
    with pytest.raises(acceptance.AcceptanceError, match="hash/read/hash audit"):
        acceptance.audit_ct_divb(policy, "R02", [f"{path}={digest}"], None)


def test_inventory_and_supplied_restart_sets_are_exact(policy, tmp_path):
    inventory_path, restarts = build_inventory(policy, tmp_path)
    path = restarts[9.0][0]
    with pytest.raises(acceptance.AcceptanceError, match="differ from the bound CT inventory"):
        acceptance.audit_ct_divb(
            policy, "R02", [f"{path}={sha256(path)}"], inventory_path
        )


def test_verify_evidence_independently_reaudits_ct_semantics(policy, tmp_path):
    inventory_path, _ = build_inventory(policy, tmp_path)
    evidence = acceptance.audit_ct_divb(policy, "R02", [], inventory_path)
    candidate = tmp_path / "ct-evidence.json"
    write_json(candidate, evidence)
    verification = acceptance.verify_evidence(policy, candidate, sha256(candidate))
    assert verification["verified"] is True

    forged = deepcopy(evidence)
    forged.pop("evidence_digest")
    forged["maximum_normalized_ct_divb"] = 0.5e-12
    forged = acceptance.seal_evidence(forged)
    forged_path = tmp_path / "forged-ct-evidence.json"
    write_json(forged_path, forged)
    with pytest.raises(acceptance.AcceptanceError, match="semantic contents differ"):
        acceptance.verify_evidence(policy, forged_path, sha256(forged_path))
