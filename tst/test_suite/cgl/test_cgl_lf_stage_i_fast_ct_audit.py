"""Focused regression tests for the direct-fast Stage I CT snapshot audit."""

from __future__ import annotations

from array import array
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
ADAPTER_PATH = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_ct_audit.py"


def load_adapter():
    name = "cgl_lf_stage_i_fast_ct_audit_test"
    spec = importlib.util.spec_from_file_location(name, ADAPTER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def adapter():
    return load_adapter()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_binary_metadata(
    path: Path,
    *,
    time: float = 1.0,
    cycle: int = 10,
    variables: tuple[str, ...] = (
        "dens", "velx", "vely", "velz", "eint", "p_perp", "bcc1", "bcc2", "bcc3"
    ),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "<mesh>\n"
        "nx1=4\nnx2=4\nnx3=4\n"
        "<meshblock>\n"
        "nx1=2\nnx2=2\nnx3=2\n"
    ).encode("utf-8")
    payload = (
        b"Athena binary output version=1.1\n"
        b"pheader_count=5\n"
        + f"time={time:.17g}\n".encode()
        + f"cycle={cycle}\n".encode()
        + b"size of location=8\n"
        + b"size of variable=8\n"
        + f"nvars={len(variables)}\n".encode()
        + ("Variables= " + " ".join(variables) + "\n").encode()
        + f"header size={len(header)}\n".encode()
        + header
        + b"fixture-field-payload\n"
    )
    path.write_bytes(payload)


def little_endian_doubles(values: list[float]) -> bytes:
    result = array("d", values)
    if sys.byteorder != "little":
        result.byteswap()
    return result.tobytes()


def region_indices(ng: int, nx1: int, nx2: int, nx3: int, *, coarse: bool) -> list[int]:
    active = [ng, ng + nx1 - 1, ng, ng + nx2 - 1, ng, ng + nx3 - 1]
    if not coarse:
        return [ng, nx1, nx2, nx3, *active, *([0] * 9)]
    return [
        ng, nx1, nx2, nx3, *active,
        nx1 // 2, nx2 // 2, nx3 // 2,
        ng, ng + nx1 // 2 - 1,
        ng, ng + nx2 // 2 - 1,
        ng, ng + nx3 // 2 - 1,
    ]


def write_native_restart(path: Path, time_value: float, *, divergent: bool = False) -> None:
    """Write one minimal qualified uniform 3-D native restart for the reviewed parser."""

    path.parent.mkdir(parents=True, exist_ok=True)
    ng = 1
    nx1 = nx2 = nx3 = 2
    parameters = (
        "<mesh>\n"
        f"nghost = {ng}\nnx1 = {nx1}\nnx2 = {nx2}\nnx3 = {nx3}\n"
        "<meshblock>\n"
        f"nx1 = {nx1}\nnx2 = {nx2}\nnx3 = {nx3}\n"
        "<mhd>\neos = cgl\nnscalars = 0\nbfloor = 1e-10\n"
        "<output3>\nsingle_file_per_rank = true\n"
        "<time>\n"
        f"restart_time = {time_value:.17g}\n"
        "<par_end>\n"
    ).encode()
    dx = 0.5
    header = struct.pack(
        "<ii9d19i19iddi",
        1,
        0,
        0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx, dx, dx,
        *region_indices(ng, nx1, nx2, nx3, coarse=False),
        *region_indices(ng, nx1, nx2, nx3, coarse=True),
        time_value,
        0.01,
        100,
    )
    meshblock_list = struct.pack("<4i", 0, 0, 0, 0) + struct.pack("<f", 1.0)
    nout1 = nout2 = nout3 = 4
    cell_count = nout1 * nout2 * nout3
    x1_count = (nout1 + 1) * nout2 * nout3
    x2_count = nout1 * (nout2 + 1) * nout3
    x3_count = nout1 * nout2 * (nout3 + 1)
    x1f = [1.0] * x1_count
    if divergent:
        x1f[(ng * nout2 + ng) * (nout1 + 1) + ng + 1] = 2.0
    block_data = little_endian_doubles(
        [0.0] * (6 * cell_count)
        + x1f
        + [0.0] * x2_count
        + [0.0] * x3_count
    )
    path.write_bytes(
        parameters
        + header
        + meshblock_list
        + b"direct-fast-fixture"
        + struct.pack("<Q", len(block_data))
        + block_data
    )


def build_campaign(
    tmp_path: Path,
    *,
    case_present: bool = True,
    complete_snapshot: bool = True,
    rank_count: int = 1,
    rank_time_offset: float = 0.0,
    restart_times: tuple[float, ...] = (),
    divergent_restart_time: float | None = None,
) -> Path:
    assembled = tmp_path / "assembled"
    matrix = tmp_path / "matrix.json"
    reporter = tmp_path / "reporter.py"
    input_path = tmp_path / "case.athinput"
    segment_dir = tmp_path / "fast_s000_t0_to_t10"
    output_dir = segment_dir / "output"
    manifest = segment_dir / "manifest/fast_run.json"
    mhd = assembled / "cases/R03/history/case.mhd.hst"
    user = assembled / "cases/R03/history/case.user.hst"
    matrix.write_text("{}\n", encoding="utf-8")
    reporter.write_text("# reporter\n", encoding="utf-8")
    input_path.write_text("<mhd>\n", encoding="utf-8")
    write_json(manifest, {
        "schema_version": 1,
        "case_id": "R03",
        "case_name": "fixture",
        "matrix_sha256": sha256(matrix),
        "input_sha256": sha256(input_path),
        "executable_sha256": "a" * 64,
        "nodes": 1,
        "ranks_per_node": rank_count,
        "ranks": rank_count,
        "run_dir": str(segment_dir.resolve()),
        "output_dir": str(output_dir.resolve()),
    })
    mhd.parent.mkdir(parents=True, exist_ok=True)
    mhd.write_text("# [1]=time [2]=mass\n0 1\n", encoding="utf-8")
    user.write_text("# [1]=time [2]=mass\n0 1\n", encoding="utf-8")
    for index, time_value in enumerate(restart_times):
        write_native_restart(
            output_dir / "rst/rank_00000000" / f"case.{index:05d}.rst",
            time_value,
            divergent=time_value == divergent_restart_time,
        )

    rank_files: list[dict[str, object]] = []
    for rank in range(rank_count):
        path = tmp_path / "snapshots" / f"rank_{rank:08d}" / "case.00000.bin"
        write_binary_metadata(path, time=1.0 + (rank_time_offset if rank else 0.0))
        rank_files.append({"path": str(path.resolve()), "size_bytes": path.stat().st_size})
    snapshot_index = {
        "schema_version": 1,
        "snapshot_count": 1,
        "complete_snapshot_count": int(complete_snapshot),
        "time_first": 1.0,
        "time_last": 1.0,
        "snapshots": [{
            "segment": "fixture",
            "lineage_order": 0,
            "time": 1.0,
            "representative": rank_files[0]["path"],
            "expected_ranks": rank_count,
            "rank_files": rank_files,
            "complete": complete_snapshot,
            "missing_rank_files": [],
            "empty_rank_files": [],
        }],
        "duplicates": [],
        "warnings": [],
    }
    snapshot_path = assembled / "cases/R03/snapshots.json"
    write_json(snapshot_path, snapshot_index)
    case = {
        "schema_version": 1,
        "case_id": "R03",
        "case_name": "fixture",
        "status": "partial",
        "errors": [],
        "input": binding(input_path),
        "lineage_identities": {
            "matrix_sha256": [sha256(matrix)],
            "input_sha256": [sha256(input_path)],
            "executable_sha256": ["a" * 64],
        },
        "lineage": [{
            "kind": "fast",
            "order": 0,
            "segment": segment_dir.name,
            "segment_dir": str(segment_dir.resolve()),
            "output": str(output_dir.resolve()),
            "ranks": rank_count,
            "case_id": "R03",
            "case_name": "fixture",
            "matrix_sha256": sha256(matrix),
            "input_sha256": sha256(input_path),
            "executable_sha256": "a" * 64,
            "manifest": binding(manifest),
        }],
        "histories": {
            "mhd": {"available": True, "binding": binding(mhd)},
            "user": {"available": True, "binding": binding(user)},
        },
        "snapshots": {
            "path": str(snapshot_path.resolve()),
            "snapshot_count": 1,
            "complete_snapshot_count": int(complete_snapshot),
        },
    }
    write_json(assembled / "cases/R03/lineage.json", case)
    inventory = {
        "schema_version": 1,
        "root": str(tmp_path.resolve()),
        "output": str(assembled.resolve()),
        "matrix": binding(matrix),
        "adapter": binding(reporter),
        "cases": {"R03": case} if case_present else {},
    }
    inventory_path = assembled / "inventory.json"
    write_json(inventory_path, inventory)
    return inventory_path


def build_accepted_r02_campaign(tmp_path: Path) -> tuple[Path, dict[float, Path]]:
    assembled = tmp_path / "assembled"
    epoch = tmp_path / "runs/E03-forcing-policy"
    bundle_root = epoch / "bundles/R02"
    matrix = tmp_path / "matrix.json"
    reporter = tmp_path / "reporter.py"
    input_path = tmp_path / "case.athinput"
    matrix.write_bytes(
        (REPOSITORY / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json").read_bytes()
    )
    reporter.write_text("# reporter\n", encoding="utf-8")
    input_path.write_text("<mhd>\n", encoding="utf-8")
    case_name = "paper_standard_active_alfvenic_beta10"
    executable_sha = "a" * 64
    segment_paths: list[str] = []
    restart_paths: dict[float, Path] = {}
    for index, time_value in enumerate((9.0, 10.0)):
        segment_name = f"s{index:02d}_t{time_value:g}"
        segment_dir = epoch / "R02" / segment_name
        segment_manifest = segment_dir / "manifest/prepared_run.json"
        restart = (
            segment_dir / "output/rst/rank_00000000"
            / f"{case_name}.{index:05d}.rst"
        )
        write_native_restart(restart, time_value)
        restart_record = {
            "path": str(restart.resolve()),
            "sha256": sha256(restart),
            "size_bytes": restart.stat().st_size,
        }
        terminal = {
            **restart_record,
            "storage": "per_rank",
            "rank_files": [restart_record],
        }
        write_json(segment_manifest, {
            "allocation": {"nodes": 1, "ranks_per_node": 1},
            "command": {
                "matrix_sha256": sha256(matrix),
                "executable_sha256": executable_sha,
            },
            "accounting": {
                "result": "accepted",
                "case_id": "R02",
                "case_name": case_name,
                "segment": segment_name,
                "executable_sha256": executable_sha,
            },
            "scientific_inspection": {
                "schema_version": 3,
                "accepted": True,
                "case_id": "R02",
                "segment": segment_name,
                "manifest": str(segment_manifest.resolve()),
                "required_time": time_value,
                "final_time": time_value,
                "terminal_restart_time": time_value,
                "checks": {
                    "required_time_reached": True,
                    "restart_retained": True,
                    "terminal_restart_physical_time_matches_final": True,
                },
                "restart_times": [time_value],
                "restarts": [terminal],
                "terminal_restart": terminal,
            },
        })
        segment_paths.append(str(segment_manifest.resolve()))
        restart_paths[time_value] = restart

    bundle_manifest = bundle_root / "manifest.json"
    write_json(bundle_manifest, {
        "workflow": "paper-mks24-stage-i-production",
        "status": "accepted_for_analysis",
        "production_case_id": "R02",
        "required_final_time": 10.0,
        "accepted_final_time": 10.0,
        "cases": [{
            "name": case_name,
            "input": "inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta10.athinput",
            "status": "passed",
        }],
        "production_segment_manifests": segment_paths,
    })

    mhd = assembled / "cases/R02/history/case.mhd.hst"
    user = assembled / "cases/R02/history/case.user.hst"
    mhd.parent.mkdir(parents=True, exist_ok=True)
    mhd.write_text("# [1]=time [2]=mass\n0 1\n", encoding="utf-8")
    user.write_text("# [1]=time [2]=mass\n0 1\n", encoding="utf-8")
    snapshot = tmp_path / "snapshots/rank_00000000/case.00000.bin"
    write_binary_metadata(snapshot, time=10.0)
    snapshot_path = assembled / "cases/R02/snapshots.json"
    write_json(snapshot_path, {
        "schema_version": 1,
        "snapshot_count": 1,
        "complete_snapshot_count": 1,
        "time_first": 10.0,
        "time_last": 10.0,
        "snapshots": [{
            "segment": "R02 accepted bundle",
            "lineage_order": 0,
            "time": 10.0,
            "representative": str(snapshot.resolve()),
            "expected_ranks": 1,
            "rank_files": [{
                "path": str(snapshot.resolve()),
                "size_bytes": snapshot.stat().st_size,
            }],
            "complete": True,
            "missing_rank_files": [],
            "empty_rank_files": [],
        }],
        "duplicates": [],
        "warnings": [],
    })
    case = {
        "schema_version": 1,
        "case_id": "R02",
        "case_name": case_name,
        "status": "complete",
        "errors": [],
        "input": binding(input_path),
        "lineage_identities": {
            "matrix_sha256": [sha256(matrix)],
            "input_sha256": [sha256(input_path)],
            "executable_sha256": [executable_sha],
        },
        "lineage": [{
            "kind": "accepted_r02_bundle",
            "order": 0,
            "segment": "R02 accepted bundle",
            "state": "accepted_for_analysis",
            "output": str(bundle_root.resolve()),
            "ranks": 1,
            "manifest": binding(bundle_manifest),
        }],
        "histories": {
            "mhd": {"available": True, "binding": binding(mhd)},
            "user": {"available": True, "binding": binding(user)},
        },
        "snapshots": {
            "path": str(snapshot_path.resolve()),
            "snapshot_count": 1,
            "complete_snapshot_count": 1,
        },
    }
    write_json(assembled / "cases/R02/lineage.json", case)
    inventory_path = assembled / "inventory.json"
    write_json(inventory_path, {
        "schema_version": 1,
        "root": str(tmp_path.resolve()),
        "output": str(assembled.resolve()),
        "matrix": binding(matrix),
        "adapter": binding(reporter),
        "cases": {"R02": case},
    })
    return inventory_path, restart_paths


def attach_exact_t9_replay(
    adapter, inventory_path: Path
) -> Path:
    """Attach one fully bound supplemental exact-t=9 replay to the fixture."""

    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    case = inventory["cases"]["R03"]
    segment = case["lineage"][0]
    segment_manifest_path = Path(segment["manifest"]["path"])
    segment_manifest = json.loads(segment_manifest_path.read_text(encoding="utf-8"))
    executable = inventory_path.parent.parent / "athena"
    executable.write_bytes(b"qualified executable fixture\n")
    executable_binding = binding(executable)
    segment_manifest["executable_sha256"] = executable_binding["sha256"]
    write_json(segment_manifest_path, segment_manifest)
    segment["manifest"] = binding(segment_manifest_path)
    segment["executable_sha256"] = executable_binding["sha256"]
    case["lineage_identities"]["executable_sha256"] = [executable_binding["sha256"]]
    case["status"] = "complete"
    lineage_path = inventory_path.parent / "cases/R03/lineage.json"
    write_json(lineage_path, case)
    inventory["cases"]["R03"] = case
    write_json(inventory_path, inventory)
    primary_binding = binding(inventory_path)
    lineage_binding = binding(lineage_path)

    output = Path(segment["output"])
    parent_path = output / "rst/rank_00000000/case.00000.rst"
    parent_rank = {**binding(parent_path), "rank": 0}
    parent = {
        "time": 8.5,
        "name": parent_path.name,
        "rank_count": 1,
        "root": str((output / "rst").resolve()),
        "source_segment_manifest": segment["manifest"],
        "source_lineage_order": 0,
        "source_segment": segment["segment"],
        "rank_files": [parent_rank],
    }

    replay_root = inventory_path.parent.parent / "ct-replay"
    run_dir = replay_root / "cases/R03/t9"
    replay_output = run_dir / "output"
    terminal_path = replay_output / "rst/rank_00000000/case.t9.rst"
    write_native_restart(terminal_path, 9.0)
    terminal = {
        "time": 9.0,
        "name": terminal_path.name,
        "rank_count": 1,
        "rank_files": [{**binding(terminal_path), "rank": 0}],
    }
    mhd = replay_output / "case.mhd.hst"
    user = replay_output / "case.user.hst"
    mhd.parent.mkdir(parents=True, exist_ok=True)
    mhd.write_text("# [1]=time [2]=mass\n8.5 1\n9 1\n", encoding="utf-8")
    user.write_text("# [1]=time [2]=mass\n8.5 1\n9 1\n", encoding="utf-8")
    exit_path = run_dir / "run_exit_code"
    environment_path = run_dir / "run_environment.txt"
    exit_path.write_text("0\n", encoding="utf-8")
    environment_path.write_text("fixture\n", encoding="utf-8")

    execution = {
        "matrix": inventory["matrix"],
        "executable": executable_binding,
        "input": case["input"],
    }
    allocation = {"nodes": 1, "ranks_per_node": 1, "ranks": 1}
    run_manifest_path = run_dir / "replay_run.json"
    write_json(run_manifest_path, {
        "schema_version": 1,
        "record_type": "stage-i-direct-fast-ct-exact-state-replay-run",
        "case_id": "R03",
        "case_name": case["case_name"],
        "purpose": "CT-only exact-state replay; excluded from accepted science histories",
        "primary_inventory": primary_binding,
        "primary_case_lineage": lineage_binding,
        "execution_identity": execution,
        "allocation": allocation,
        "parent_restart": parent,
        "target_time": 9.0,
        "command_line_overrides": ["time/tlim=9.0"],
        "run_basename": "fixture",
        "paths": {
            "run_dir": str(run_dir.resolve()),
            "output_dir": str(replay_output.resolve()),
            "slurm_log": str((run_dir / "fixture.%j.log").resolve()),
        },
        "tools": {
            "replay_tool": binding(adapter.CT_REPLAY_PATH),
            "restart_parser": binding(
                adapter.CT_REPLAY_PATH.with_name("cgl_lf_stage_i_fast.py")
            ),
        },
        "job": {
            "account": "AST207",
            "partition": "extended",
            "walltime": "04:00:00",
            "job_id": "12345",
        },
    })
    completion_path = run_dir / "completion.json"
    write_json(completion_path, {
        "schema_version": 1,
        "record_type": "stage-i-direct-fast-ct-exact-state-replay-completion",
        "case_id": "R03",
        "case_name": case["case_name"],
        "purpose": "CT-only exact-state replay; excluded from accepted science histories",
        "target_time": 9.0,
        "command_line_overrides": ["time/tlim=9.0"],
        "primary_inventory": primary_binding,
        "primary_case_lineage": lineage_binding,
        "run_manifest": binding(run_manifest_path),
        "source_segment_manifest": segment["manifest"],
        "execution_identity": execution,
        "allocation": allocation,
        "parent_restart": parent,
        "scheduler": {
            "job_id": "12345",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "start": "2026-06-09T00:00:00",
            "end": "2026-06-09T00:01:00",
        },
        "run_exit_code": binding(exit_path),
        "run_environment": binding(environment_path),
        "scientific_sanity": {
            "final_time": 9.0,
            "strict_lf_failure_maxima": {"lf_nonfin": 0.0},
            "mass_relative_drift": 0.0,
            "mhd_user_mass_relative_mismatch": 0.0,
            "mhd_history": binding(mhd),
            "user_history": binding(user),
            "terminal_restart": terminal,
        },
    })
    plan_path = replay_root / "plan.json"
    write_json(plan_path, {"fixture": True})
    replay_inventory_path = replay_root / "inventory.json"
    write_json(replay_inventory_path, {
        "schema_version": 1,
        "record_type": "stage-i-direct-fast-ct-exact-state-replay-inventory",
        "purpose": "supplemental CT-only exact-t=9 evidence",
        "target_time": 9.0,
        "primary_inventory": primary_binding,
        "plan": binding(plan_path),
        "replay_tool": binding(adapter.CT_REPLAY_PATH),
        "cases": {"R03": binding(completion_path)},
    })
    return replay_inventory_path


def test_authenticated_bcc_snapshot_is_inconclusive_and_deterministic(adapter, tmp_path):
    inventory = build_campaign(tmp_path)
    observed = adapter.build_audit(inventory, sha256(inventory), [], "all")
    case = observed["cases"]["R03"]

    assert observed["result"] == "inconclusive"
    assert observed["summary"]["ct_pass_case_count"] == 0
    assert case["audit_status"] == "authenticated_inconclusive"
    assert case["provenance_authenticated"] is True
    assert case["snapshot_content_authenticated"] is True
    assert case["ct_evidence_available"] is False
    assert case["ct_claim_supported"] is False
    assert case["snapshot_audits"][0]["metadata"][
        "cell_centered_magnetic_triplet_complete"
    ] is True
    assert case["snapshot_audits"][0]["metadata"]["face_centered_ct_state_available"] is False
    assert case["snapshot_audits"][0]["ct_divb_measurement"] is None
    native = case["native_restart_ct"]
    assert native["result"] == "inconclusive"
    assert native["ct_evidence_available"] is False
    assert native["public_audit_ct_divb_without_inventory"][
        "can_establish_complete_multirank_state_coverage"
    ] is False
    assert "do not declare produced restart groups" in native["restart_inventory_limitation"]
    assert any("exact required native restart state t=9" in item for item in native["blockers"])

    repeated = adapter.build_audit(inventory, sha256(inventory), [], "all")
    assert adapter.canonical_json(observed) == adapter.canonical_json(repeated)
    assert adapter.render_csv(observed) == adapter.render_csv(repeated)
    assert adapter.render_markdown(observed) == adapter.render_markdown(repeated)
    output = tmp_path / "audit-output"
    adapter.write_outputs(output, observed)
    assert json.loads((output / "ct_audit.json").read_text(encoding="utf-8")) == observed
    assert (output / "ct_audit.csv").read_text(encoding="utf-8") == adapter.render_csv(observed)
    assert (output / "ct_audit.md").read_text(encoding="utf-8") == adapter.render_markdown(observed)


def test_exact_native_t9_t10_restart_groups_use_reviewed_ct_parser(adapter, tmp_path):
    inventory = build_campaign(tmp_path, restart_times=(9.0, 10.0))
    observed = adapter.build_audit(inventory, sha256(inventory), [], "all")
    case = observed["cases"]["R03"]
    native = case["native_restart_ct"]

    assert observed["result"] == "pass"
    assert case["audit_status"] == "authenticated_complete"
    assert case["ct_result"] == "pass"
    assert case["ct_evidence_available"] is True
    assert case["ct_claim_supported"] is True
    assert native["coverage_complete"] is True
    assert native["audited_state_times"] == [9.0, 10.0]
    assert native["maximum_normalized_ct_divb"] == 0.0
    assert native["meshblock_coverage_complete"] is True
    assert native["sampled_meshblock_fraction"] == 1.0
    assert native["campaign_authority_eligible"] is False
    assert native["release_authorizing"] is False
    assert native["blockers"] == []
    assert len(native["state_audits"]) == 2
    assert all(state["meshblock_coverage_complete"] for state in native["state_audits"])

    acceptance, policy, _, _, _ = adapter.load_native_ct_authority()
    supplied = [
        f"{restart['path']}={restart['sha256']}"
        for state in native["state_audits"]
        for restart in state["restart_audits"]
    ]
    public = acceptance.audit_ct_divb(policy, "R03", supplied, None)
    assert public["numerical_result"] == "pass"
    assert public["coverage_complete"] is False
    assert public["result"] == "inconclusive"


def test_nonexact_restart_time_documents_exact_state_blocker(adapter, tmp_path):
    inventory = build_campaign(tmp_path, restart_times=(9.0001, 10.0))
    observed = adapter.build_audit(inventory, sha256(inventory), [], "all")
    native = observed["cases"]["R03"]["native_restart_ct"]

    assert observed["result"] == "inconclusive"
    assert native["numerical_result"] == "pass"
    assert native["coverage_complete"] is False
    assert native["audited_state_times"] == [10.0]
    assert native["missing_required_state_times"] == [9.0]
    assert native["discovered_state_times"] == [9.0001, 10.0]
    assert any("exact required native restart state t=9" in item for item in native["blockers"])


def test_exact_t9_replay_supplements_selected_lineage_t10(adapter, tmp_path):
    inventory = build_campaign(tmp_path, restart_times=(8.5, 10.0))
    replay_inventory = attach_exact_t9_replay(adapter, inventory)

    observed = adapter.build_audit(
        inventory,
        sha256(inventory),
        ["R03"],
        "all",
        replay_inventory,
        sha256(replay_inventory),
    )
    native = observed["cases"]["R03"]["native_restart_ct"]

    assert observed["result"] == "pass"
    assert native["audited_state_times"] == [9.0, 10.0]
    assert native["coverage_complete"] is True
    assert len(native["authenticated_exact_state_replays"]) == 1
    assert [state["source"]["source_kind"] for state in native["state_audits"]] == [
        "exact_state_replay",
        "direct_fast",
    ]


def test_native_face_field_divergence_is_reported_as_failure(adapter, tmp_path):
    inventory = build_campaign(
        tmp_path, restart_times=(9.0, 10.0), divergent_restart_time=9.0
    )
    observed = adapter.build_audit(inventory, sha256(inventory), [], "all")
    case = observed["cases"]["R03"]
    native = case["native_restart_ct"]

    assert observed["result"] == "fail"
    assert case["ct_result"] == "fail"
    assert case["ct_claim_supported"] is True
    assert native["coverage_complete"] is True
    assert native["numerical_result"] == "fail"
    assert native["maximum_normalized_ct_divb"] > native["normalized_ct_divb_lt"]


def test_partial_campaign_and_incomplete_snapshot_remain_inconclusive(adapter, tmp_path):
    inventory = build_campaign(tmp_path, complete_snapshot=False)
    observed = adapter.build_audit(inventory, sha256(inventory), ["R03", "R04"], "all")

    assert observed["result"] == "inconclusive"
    assert observed["cases"]["R03"]["audit_status"] == "partial_inconclusive"
    assert observed["cases"]["R03"]["snapshot_summary"]["authenticated_snapshot_count"] == 0
    assert observed["cases"]["R04"]["audit_status"] == "partial_inconclusive"
    assert observed["cases"]["R04"]["assembled_status"] == "not_present"


def test_inventory_and_case_provenance_authentication_fail_closed(adapter, tmp_path):
    inventory = build_campaign(tmp_path)
    with pytest.raises(adapter.FastCtAuditError, match="inventory SHA-256"):
        adapter.build_audit(inventory, "0" * 64, [], "all")

    lineage_path = inventory.parent / "cases/R03/lineage.json"
    lineage = json.loads(lineage_path.read_text(encoding="utf-8"))
    lineage["status"] = "complete"
    write_json(lineage_path, lineage)
    observed = adapter.build_audit(inventory, sha256(inventory), [], "all")
    case = observed["cases"]["R03"]
    assert observed["result"] == "authentication_failed"
    assert case["audit_status"] == "authentication_failed"
    assert case["provenance_authenticated"] is False


def test_rank_metadata_mismatch_is_not_reported_as_ct_evidence(adapter, tmp_path):
    inventory = build_campaign(tmp_path, rank_count=2, rank_time_offset=0.25)
    observed = adapter.build_audit(inventory, sha256(inventory), [], "all")
    case = observed["cases"]["R03"]

    assert observed["result"] == "authentication_failed"
    assert case["ct_evidence_available"] is False
    assert case["ct_claim_supported"] is False
    assert "differs" in case["errors"][0]


def test_accepted_r02_bundle_reports_numerical_ct_under_authority_blocker(
    adapter, tmp_path, monkeypatch
):
    inventory, _ = build_accepted_r02_campaign(tmp_path)
    acceptance, _, _, _, _ = adapter.load_native_ct_authority()
    blocker = "published F118 audit F116 predecessor differs"

    def blocked_authority(_policy):
        raise acceptance.AcceptanceError(blocker)

    monkeypatch.setattr(acceptance, "current_source_archive_catalog", blocked_authority)
    observed = adapter.build_audit(inventory, sha256(inventory), ["R02"], "latest")
    native = observed["cases"]["R02"]["native_restart_ct"]

    assert observed["result"] == "pass"
    assert native["result"] == "pass"
    assert native["audited_state_times"] == [9.0, 10.0]
    assert native["maximum_normalized_ct_divb"] == 0.0
    assert native["coverage_complete"] is True
    assert native["campaign_authority_eligible"] is False
    assert native["release_authorizing"] is False
    assert native["authority_chain"]["status"] == "blocked"
    assert native["authority_blockers"] == [blocker]
    assert len(native["authenticated_accepted_bundles"]) == 1
    bundle = native["authenticated_accepted_bundles"][0]
    assert bundle["accepted_segment_count"] == 2
    assert all(item["manifest"]["sha256"] for item in bundle["accepted_segments"])
    assert {
        state["source"]["source_kind"] for state in native["state_audits"]
    } == {"accepted_r02_bundle"}
    assert blocker in adapter.render_csv(observed)
    assert f"**R02 authority-chain blocker:** {blocker}" in adapter.render_markdown(observed)


def test_f118_audit_predecessor_digest_schema_is_normalized_exactly(
    adapter, tmp_path, monkeypatch
):
    acceptance, _, _, _, _ = adapter.load_native_ct_authority()
    audit_path = tmp_path / "F118.publication_audit.json"
    canonical = {
        "evidence": "a" * 64,
        "provenance_review": "b" * 64,
        "plasma_review": "c" * 64,
        "publication_audit": "d" * 64,
    }
    published = {f"{key}_sha256": value for key, value in canonical.items()}
    write_json(audit_path, {"historical_f116_authority": published})
    policy = {
        "source_catalog_policy": {
            "successor_paths": {"publication_audit": str(audit_path)}
        }
    }
    original_load_json = acceptance.load_json

    def current_source_archive_catalog(_policy):
        audit, audit_binding = acceptance.load_json(
            audit_path, "published F118 source-authority publication_audit"
        )
        assert audit["historical_f116_authority"] == canonical
        return {}, {"checkpoint": "F-118"}, [audit_binding]

    monkeypatch.setattr(
        acceptance, "current_source_archive_catalog", current_source_archive_catalog
    )
    observed = adapter.r02_authority_chain_record(acceptance, policy)

    assert observed["status"] == "authenticated_non_authorizing"
    assert observed["authority"]["checkpoint"] == "F-118"
    assert observed["bindings"] == [binding(audit_path)]
    assert acceptance.load_json is original_load_json
    assert json.loads(audit_path.read_text(encoding="utf-8"))[
        "historical_f116_authority"
    ] == published


@pytest.mark.parametrize(
    "historical",
    [
        {
            "evidence_sha256": "a" * 64,
            "provenance_review_sha256": "b" * 64,
            "plasma_review_sha256": "c" * 64,
        },
        {
            "evidence_sha256": "a" * 64,
            "provenance_review_sha256": "b" * 64,
            "plasma_review_sha256": "c" * 64,
            "publication_audit_sha256": "d" * 64,
            "unexpected_sha256": "e" * 64,
        },
        {
            "evidence_sha256": "not-a-digest",
            "provenance_review_sha256": "b" * 64,
            "plasma_review_sha256": "c" * 64,
            "publication_audit_sha256": "d" * 64,
        },
    ],
)
def test_f118_audit_predecessor_digest_normalization_fails_closed(
    adapter, tmp_path, monkeypatch, historical
):
    acceptance, _, _, _, _ = adapter.load_native_ct_authority()
    audit_path = tmp_path / "F118.publication_audit.json"
    write_json(audit_path, {"historical_f116_authority": historical})
    policy = {
        "source_catalog_policy": {
            "successor_paths": {"publication_audit": str(audit_path)}
        }
    }
    original_load_json = acceptance.load_json

    def current_source_archive_catalog(_policy):
        acceptance.load_json(
            audit_path, "published F118 source-authority publication_audit"
        )
        raise AssertionError("malformed F118 predecessor digest schema was accepted")

    monkeypatch.setattr(
        acceptance, "current_source_archive_catalog", current_source_archive_catalog
    )
    observed = adapter.r02_authority_chain_record(acceptance, policy)

    assert observed["status"] == "blocked"
    assert "published F118" in observed["blocker"]
    assert acceptance.load_json is original_load_json


def test_accepted_r02_declared_restart_bytes_fail_closed(adapter, tmp_path, monkeypatch):
    inventory, restarts = build_accepted_r02_campaign(tmp_path)
    acceptance, _, _, _, _ = adapter.load_native_ct_authority()
    monkeypatch.setattr(
        acceptance,
        "current_source_archive_catalog",
        lambda _policy: (_ for _ in ()).throw(
            acceptance.AcceptanceError("authority chain fixture blocker")
        ),
    )
    restarts[9.0].write_bytes(restarts[9.0].read_bytes() + b"tampered")
    observed = adapter.build_audit(inventory, sha256(inventory), ["R02"], "latest")
    case = observed["cases"]["R02"]
    native = case["native_restart_ct"]

    assert observed["result"] == "authentication_failed"
    assert case["audit_status"] == "authentication_failed"
    assert native["status"] == "authentication_failed"
    assert native["audited_state_times"] == [10.0]
    assert any("differs from its authority" in error for error in native["errors"])
