"""Focused regressions for the hardened Figure 13 instability product."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
ANALYZER_PATH = REPO_ROOT / "scripts" / "analyze_cgl_lf_paper.py"
WORKFLOW_PATH = REPO_ROOT / "scripts" / "cgl_lf_workflow.py"
TEST_REVISION = "0123456789abcdef0123456789abcdef01234567"
TEST_EXECUTABLE_SHA256 = hashlib.sha256(b"test executable").hexdigest()


def athena_snapshot_payload(time, cycle=0, rank=0):
    """Return a scaled production-topology rank-local Athena snapshot."""

    variables = (
        "dens", "velx", "vely", "velz", "eint", "p_perp",
        "bcc1", "bcc2", "bcc3",
    )
    parameters = (
        "<mesh>\n"
        "nx1=6\nnx2=6\nnx3=12\nnghost=0\n"
        "x1min=0\nx1max=1\nx2min=0\nx2max=1\nx3min=0\nx3max=2\n"
        "<meshblock>\n"
        "nx1=1\nnx2=1\nnx3=2\n"
    ).encode()
    preheader = (
        "Athena binary output version=1.1\n"
        "size of preheader=5\n"
        f"time={time:.17g}\n"
        f"cycle={cycle}\n"
        "size of location=8\n"
        "size of variable=8\n"
        f"number of variables={len(variables)}\n"
        f"variables: {' '.join(variables)}\n"
        f"header offset={len(parameters)}\n"
    ).encode()
    values = (1.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0)
    cell_values = tuple(
        value for value in values for _ in range(2)
    )
    meshblocks = []
    for k in range(6):
        for j in range(6):
            for i in range(6):
                linear = i + 6 * (j + 6 * k)
                if linear % 2 != rank:
                    continue
                geometry = (
                    i / 6.0,
                    (i + 1) / 6.0,
                    j / 6.0,
                    (j + 1) / 6.0,
                    k / 3.0,
                    (k + 1) / 3.0,
                )
                meshblocks.append(
                    struct.pack("<6i", 0, 0, 0, 0, 0, 1)
                    + struct.pack("<4i", i, j, k, 0)
                    + struct.pack("<6d", *geometry)
                    + struct.pack(f"<{len(cell_values)}d", *cell_values)
                )
    assert len(meshblocks) == 108
    return preheader + parameters + b"".join(meshblocks)


def athena_restart_payload(time, meshblocks=2, local_meshblocks=1):
    """Return one deliberately tiny restart for the lightweight unit adapter."""

    parameter_dump = (
        "#------------------------- PAR_DUMP -------------------------\n"
        "<job>\n"
        "basename = test\n"
        "<mesh>\n"
        "nx1 = 2\n"
        "<meshblock>\n"
        "nx1 = 1\n"
        "<time>\n"
        f"restart_time = {time:.17g}\n"
        "<par_end>\n"
    ).encode()
    header = bytearray(252)
    struct.pack_into("<ii", header, 0, meshblocks, 0)
    struct.pack_into("<d", header, 232, time)
    struct.pack_into("<d", header, 240, 0.01)
    struct.pack_into("<i", header, 248, 1)
    locations = b"".join(
        struct.pack("<4i", index, 0, 0, 0) for index in range(meshblocks)
    )
    costs = struct.pack(f"<{meshblocks}f", *([1.0] * meshblocks))
    turbulence = bytearray(248)
    struct.pack_into("<i", turbulence, 0, 3)
    struct.pack_into("<i", turbulence, 4, 1)
    struct.pack_into("<i", turbulence, 23 * 4, 1)
    rng_state = bytes(296)
    amplitudes = struct.pack("<6d", *([0.0] * 6))
    injected_work = struct.pack("<d", 0.0)
    lf_diagnostics = struct.pack("<18d", *([0.0] * 18))
    data_size = 16
    meshblock_payload = struct.pack(
        f"<{2 * local_meshblocks}d", *([1.0] * (2 * local_meshblocks))
    )
    return (
        parameter_dump
        + header
        + locations
        + costs
        + turbulence
        + rng_state
        + amplitudes
        + injected_work
        + lf_diagnostics
        + struct.pack("<Q", data_size)
        + meshblock_payload
    )


def load_module(name, path):
    """Import one repository script without depending on the test-suite cwd."""

    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def install_lightweight_unit_restart_adapter(analyzer):
    """Keep broad unit fixtures small while retaining the real parser for contract tests."""

    analyzer.real_restart_loadability_evidence = analyzer.restart_loadability_evidence

    def restart_loadability_evidence(
        path, expected_time, restart_binary_abi, input_contract, runtime_contract,
    ):
        del input_contract, runtime_contract
        text, parameter_size = analyzer.restart_parameter_dump(path)
        marker = analyzer.restart_marker_time(text, path)
        contract = restart_binary_abi["contract"]
        with path.open("rb") as stream:
            stream.seek(parameter_size)
            header = stream.read(int(contract["mesh_header_size"]))
            if len(header) != int(contract["mesh_header_size"]):
                raise ValueError(f"restart mesh header is truncated: {path}")
            meshblocks = struct.unpack_from("<i", header)[0]
            binary_time = struct.unpack_from("<d", header, 232)[0]
            dt = struct.unpack_from("<d", header, 240)[0]
            cycle = struct.unpack_from("<i", header, 248)[0]
            locations = stream.read(meshblocks * int(contract["logical_location_size"]))
            costs = stream.read(meshblocks * int(contract["cost_size"]))
            metadata = stream.read(int(contract["turbulence_metadata_size"]))
            rng = stream.read(int(contract["rng_state_size"]))
            mode_count = struct.unpack_from("<i", metadata, 4)[0]
            amplitudes = stream.read(6 * mode_count * 8)
            injected = stream.read(8)
            diagnostics = stream.read(int(contract["lf_diagnostic_count"]) * 8)
            encoded_size = stream.read(8)
            if min(
                len(locations), len(costs), len(metadata), len(rng), len(amplitudes),
                len(injected), len(diagnostics), len(encoded_size),
            ) <= 0:
                raise ValueError(f"restart binary layout is truncated: {path}")
            data_size = struct.unpack("<Q", encoded_size)[0]
            payload = stream.read()
        if (
            not analyzer.figure_13_time_close(marker, binary_time)
            or not analyzer.figure_13_time_close(binary_time, expected_time)
            or data_size <= 0
            or len(payload) == 0
            or len(payload) % data_size
        ):
            raise ValueError(f"restart meshblock payload is truncated: {path}")
        if not np.isfinite(np.frombuffer(payload, dtype="<f8")).all():
            raise ValueError(f"restart meshblock payload is nonfinite: {path}")
        local_blocks = len(payload) // data_size
        return {
            "binary_layout_validation": "lightweight-unit-adapter",
            "restart_binary_abi_sha256": restart_binary_abi["contract_sha256"],
            "parameter_dump_size": parameter_size,
            "parameter_dump_sha256": hashlib.sha256(
                (text + "<par_end>\n").encode()
            ).hexdigest(),
            "meshblock_count": meshblocks,
            "local_meshblock_count": local_blocks,
            "binary_time": binary_time,
            "dt": dt,
            "cycle": cycle,
            "logical_location_count": meshblocks,
            "mode_count": mode_count,
            "variable_data_size": data_size,
            "logical_location_inventory_sha256": hashlib.sha256(locations).hexdigest(),
            "cost_inventory_sha256": hashlib.sha256(costs).hexdigest(),
            "turbulence_metadata_sha256": hashlib.sha256(metadata).hexdigest(),
            "turbulence_metadata_contract": {"unit_adapter": True},
            "rng_lifecycle": "unit-adapter",
            "rng_n_updates": 0,
            "rng_canonical_continuation_state_sha256": hashlib.sha256(rng).hexdigest(),
            "rng_authenticated_field_count": 0,
            "turbulence_amplitudes_sha256": hashlib.sha256(amplitudes).hexdigest(),
            "injected_work": struct.unpack("<d", injected)[0],
            "lf_diagnostics": list(struct.unpack("<18d", diagnostics)),
            "mesh_region_contract": "unit-adapter",
        }

    analyzer.restart_loadability_evidence = restart_loadability_evidence


@pytest.fixture(scope="module")
def analyzer():
    module = load_module("analyze_cgl_lf_paper_firehose_test", ANALYZER_PATH)
    install_lightweight_unit_restart_adapter(module)
    return module


@pytest.fixture(scope="module")
def workflow():
    return load_module("cgl_lf_workflow_firehose_test", WORKFLOW_PATH)


def occupancy_fields(beta_delta):
    """Construct positive-pressure, unit-B fields with selected beta-Delta."""

    beta_delta = np.asarray(beta_delta, dtype=float)
    ppar = np.full(beta_delta.shape, 2.0)
    return {
        "eint": ppar,
        "p_perp": ppar + 0.5 * beta_delta,
        "bcc1": np.ones(beta_delta.shape),
        "bcc2": np.zeros(beta_delta.shape),
        "bcc3": np.zeros(beta_delta.shape),
    }


def fake_snapshot_provenance(analyzer, case_id, time, salt=""):
    """Return deterministic synthetic provenance accepted by the product."""

    payload = f"{case_id}:{time:.2f}:{salt}".encode()
    files = [{
        "path": f"/archive/{case_id}/rank_00000000/snapshot.{time:.2f}.bin",
        "symlink_target": None,
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }]
    return {
        "representative_path": files[0]["path"],
        "layout": "rank_local_siblings",
        "expected_rank_count": 1,
        "rank_directory_names": ["rank_00000000"],
        "files": files,
        "aggregate_sha256": analyzer.canonical_json_sha256(files),
    }


def sample_record(
    analyzer, time, beta_delta, case_id="R15", salt="",
    lengths=(1.0, 1.0, 2.0), grid_shape=None,
):
    """Build an exact categorical occupancy record, optionally tiled to a grid."""

    fields = occupancy_fields(beta_delta)
    occupancy = analyzer.firehose_threshold_occupancy(fields)
    shape = fields["eint"].shape
    if grid_shape is not None:
        target_count = int(np.prod(grid_shape))
        source_count = int(np.prod(shape))
        if target_count % source_count:
            raise ValueError("categorical occupancy fixture does not tile target grid")
        scale = target_count // source_count
        for semantics in occupancy["occupancy"].values():
            for criterion in semantics.values():
                criterion["cell_count"] *= scale
        occupancy["normalization"]["total_cell_count"] = target_count
        shape = tuple(grid_shape)
    return {
        "time": time,
        "shape_z_y_x": list(shape),
        "lengths_x_y_z": list(lengths),
        "snapshot_provenance": fake_snapshot_provenance(
            analyzer, case_id, time, salt
        ),
        "firehose_threshold_occupancy": occupancy,
    }


def analyzed_case(
    analyzer, case_id, beta_delta, times=None, salt="", include_provenance=True,
    exact_grid=True,
):
    """Return one analyzed-case fixture with an occupancy ensemble."""

    if times is None:
        times = analyzer.FIGURE_13_SNAPSHOT_TIMES
    records = [
        sample_record(
            analyzer,
            time,
            beta_delta,
            case_id,
            salt,
            grid_shape=(
                analyzer.FIGURE_13_GRID_SHAPE_Z_Y_X if exact_grid else None
            ),
        )
        for time in times
    ]
    if not include_provenance:
        for record in records:
            del record["snapshot_provenance"]
    occupancy = analyzer.average_firehose_threshold_occupancy(records)
    occupancy["analysis_window"].update({
        "requested_time_start": analyzer.FIGURE_13_WINDOW_START,
        "requested_time_end": analyzer.FIGURE_13_WINDOW_END,
    })
    return {"snapshot_ensemble": {"firehose_threshold_occupancy": occupancy}}


def figure_13_case(analyzer, workflow, bundle, case_id):
    """Create one exact archived case and execution-input copy."""

    contract = analyzer.FIGURE_13_CASE_CONTRACTS[case_id]
    source = REPO_ROOT / contract["input"]
    execution_name = f"inputs/{case_id}.athinput"
    execution = bundle / execution_name
    execution.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, execution)
    model = workflow.model_choices(source.read_text(encoding="utf-8"), [])
    assert analyzer.canonical_json_sha256(model) == contract["model_choices_sha256"]
    return {
        "case_id": case_id,
        "name": contract["name"],
        "input": contract["input"],
        "execution_input": execution_name,
        "model_choices": model,
        "status": "passed",
        "lf_active": True,
        "amr": False,
        "outputs": {},
    }


def write_manifest(bundle, analyzer, cases, production_case_id=None):
    """Retain one deterministic bundle manifest."""

    bundle.mkdir(parents=True, exist_ok=True)
    manifest = {
        "workflow": analyzer.STAGE_I_PRODUCTION_WORKFLOW,
        "execution_epoch": "E03-forcing-policy",
        "git_revision": "0123456789abcdef",
        "cases": cases,
    }
    if production_case_id is not None:
        manifest["production_case_id"] = production_case_id
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def bind_snapshot_outputs(case, analyzed):
    """Bind an analyzed fixture to the case manifest paths it represents."""

    occupancy = analyzed["snapshot_ensemble"]["firehose_threshold_occupancy"]
    case["outputs"]["snapshot_paths"] = [
        record["representative_path"]
        for record in occupancy.get("snapshot_provenance", [])
    ]


def invocation(analyzer, bundle, output_dir, time_start=None, time_end=None):
    """Return valid analyzer/config provenance for a direct product call."""

    return analyzer.analysis_invocation_provenance({
        "snapshots": [],
        "bundle": str(bundle),
        "history": [],
        "lf_history": [],
        "output_dir": str(output_dir),
        "pdf_bins": 64,
        "alignment_shells": [1, 2, 3],
        "eddy_samples": 0,
        "eddy_bins": 20,
        "eddy_seed": 0,
        "time_start": time_start,
        "time_end": time_end,
        "reference_curves": [],
        "stage_i_manifest": None,
        "allow_partial_reference_cases": False,
        "figure_13_only": False,
        "synthetic_test": False,
    })


def complete_product_fixture(analyzer, workflow, tmp_path, salt=""):
    """Return a complete four-case fixture authenticated by controller artifacts."""

    label = salt or "base"
    canonical_root = tmp_path / "canonical"
    analyzer.STAGE_I_CANONICAL_ROOT = canonical_root
    executable = canonical_root / "build" / "test" / "src" / "athena"
    executable.parent.mkdir(parents=True, exist_ok=True)
    executable.write_bytes(b"test executable")
    executable.chmod(0o755)
    production_abi = next(iter(
        analyzer.FIGURE_13_QUALIFIED_RESTART_BINARY_ABIS.values()
    ))
    analyzer.FIGURE_13_QUALIFIED_RESTART_BINARY_ABIS[
        (TEST_REVISION, TEST_EXECUTABLE_SHA256)
    ] = dict(production_abi)
    restart_load_evidence_path = (
        canonical_root / analyzer.FIGURE_13_QUALIFIED_RESTART_LOAD_EVIDENCE_RELATIVE
    )
    restart_load_evidence_path.parent.mkdir(parents=True, exist_ok=True)
    restart_load_evidence = {
        "schema_version": 2,
        "acceptance": {
            "qualification_role": analyzer.FIGURE_13_QUALIFIED_RESTART_LOAD_ROLE,
            "staged_restart_boundary": {
                "content_set_matches_g026_selected_boundary": True,
                "restart_time": 0.01,
            },
            "terminal_identity": {"within_tolerance": True},
        },
        "checks": {"terminal_restart_metadata": {"restart_time": 0.02}},
        "execution": {"git_revision": TEST_REVISION, "result": "passed"},
        "provenance": {
            "executable": {
                "path": str(executable.resolve()),
                "sha256": TEST_EXECUTABLE_SHA256,
            },
        },
        "review": {
            "derived_by_read_only_inspection": True,
            "local_review": "passed",
        },
    }
    restart_load_evidence_path.write_text(
        json.dumps(restart_load_evidence, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    restart_load_evidence_sha256 = hashlib.sha256(
        restart_load_evidence_path.read_bytes()
    ).hexdigest()
    build_manifest = canonical_root / "runs" / "build-manifests" / "test"
    build_manifest.mkdir(parents=True, exist_ok=True)
    approval_token = {
        "schema_version": 1,
        "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
        "approved_executable": str(executable),
        "approved_executable_revision": TEST_REVISION,
        "approved_executable_sha256": TEST_EXECUTABLE_SHA256,
        "build_manifest": str(build_manifest),
        "review_notes": f"Approved qualification; g027 {restart_load_evidence_sha256};",
    }
    approval_path = (
        canonical_root / "accounting"
        / "mks24_stage_i_E03_forcing_policy_qualification_approval.json"
    )
    approval_path.parent.mkdir(parents=True, exist_ok=True)
    approval_path.write_text(
        json.dumps(approval_token, sort_keys=True) + "\n", encoding="utf-8"
    )
    qualification_approval = {
        "path": str(approval_path),
        "sha256": hashlib.sha256(approval_path.read_bytes()).hexdigest(),
        "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
        "approved_executable_revision": TEST_REVISION,
        "approved_executable_sha256": TEST_EXECUTABLE_SHA256,
        "token": approval_token,
    }
    bundle = (
        canonical_root / analyzer.STAGE_I_RUNS_RELATIVE
        / analyzer.FIGURE_13_EXECUTION_EPOCH / "bundles" / f"bundle_{label}"
    )
    cases = [
        figure_13_case(analyzer, workflow, bundle, case_id)
        for case_id in analyzer.FIGURE_13_FIREHOSE_CASES
    ]
    beta_delta = [[[-3.0, -1.5, 0.0, 2.0]]]
    analyzed = {
        case["name"]: analyzed_case(
            analyzer, case["case_id"], beta_delta, salt=salt
        )
        for case in cases
    }
    source_bundle = canonical_root / "source-archives" / f"{label}.bundle"
    source_bundle.parent.mkdir(parents=True, exist_ok=True)
    source_bundle.write_bytes(f"source bundle {label}".encode())
    source_bundle_sha256 = analyzer.sha256_file(source_bundle)
    controller_artifacts = []
    for case in cases:
        case_id = case["case_id"]
        segment = f"s_{label}"
        segment_root = (
            canonical_root / analyzer.STAGE_I_RUNS_RELATIVE
            / analyzer.FIGURE_13_EXECUTION_EPOCH / case_id / segment
        )
        artifact_path = segment_root / "manifest" / "prepared_run.json"
        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        submitted_input = artifact_path.parent / "submitted_input.athinput"
        shutil.copyfile(
            REPO_ROOT / analyzer.FIGURE_13_CASE_CONTRACTS[case_id]["input"],
            submitted_input,
        )
        provenance = []
        inspected_snapshots = []
        declared_paths = []
        for index, time in enumerate(analyzer.FIGURE_13_SNAPSHOT_TIMES):
            filename = f"{case_id}.out2.{index:05d}.bin"
            source_rank_files = []
            bundle_rank_files = []
            inspected_rank_files = []
            for rank in range(2):
                source_rank = (
                    segment_root / "output" / "bin" / f"rank_{rank:08d}"
                    / filename
                )
                source_rank.parent.mkdir(parents=True, exist_ok=True)
                source_rank.write_bytes(athena_snapshot_payload(time, index, rank))
                bundle_rank = (
                    bundle / "cases" / case["name"] / "bin"
                    / f"rank_{rank:08d}" / filename
                )
                bundle_rank.parent.mkdir(parents=True, exist_ok=True)
                bundle_rank.symlink_to(source_rank)
                source_rank_files.append(source_rank)
                bundle_rank_files.append(bundle_rank)
                inspected_rank_files.append({
                    "path": str(source_rank),
                    "sha256": analyzer.sha256_file(source_rank),
                    "size_bytes": source_rank.stat().st_size,
                })
            snapshot_provenance = analyzer.snapshot_digest_provenance(
                bundle_rank_files[0], expected_ranks=2
            )
            snapshot_provenance["snapshot_time"] = time
            provenance.append(snapshot_provenance)
            declared_paths.append(str(bundle_rank_files[0].relative_to(bundle)))
            inspected_snapshots.append({
                "path": str(source_rank_files[0]),
                "rank_files": inspected_rank_files,
            })
        occupancy = analyzed[case["name"]]["snapshot_ensemble"][
            "firehose_threshold_occupancy"
        ]
        occupancy["snapshot_provenance"] = provenance
        occupancy["snapshot_provenance_sha256"] = analyzer.canonical_json_sha256(
            provenance
        )
        case["outputs"]["snapshot_paths"] = declared_paths
        terminal_restart_rank_files = []
        for rank in range(2):
            restart = (
                segment_root / "output" / "rst" / f"rank_{rank:08d}"
                / f"{case_id}.terminal.rst"
            )
            restart.parent.mkdir(parents=True, exist_ok=True)
            restart.write_bytes(
                athena_restart_payload(analyzer.FIGURE_13_ACCEPTED_FINAL_TIME)
            )
            terminal_restart_rank_files.append({
                "path": str(restart),
                "sha256": analyzer.sha256_file(restart),
                "size_bytes": restart.stat().st_size,
            })
        terminal_restart = {
            **terminal_restart_rank_files[0],
            "storage": "per_rank",
            "rank_files": terminal_restart_rank_files,
        }
        artifact = {
            "project_root": str(canonical_root),
            "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
            "state": "recorded",
            "run": {
                "case_id": case_id,
                "case_name": case["name"],
                "resolution": analyzer.FIGURE_13_CASE_CONTRACTS[case_id]["resolution"],
                "run_basename": f"{case_id}_{segment}",
                "segment": segment,
            },
            "allocation": {"nodes": 1, "ranks_per_node": 2},
            "command": {
                "input_revision": TEST_REVISION,
                "executable": str(executable),
                "executable_revision": TEST_REVISION,
                "executable_sha256": TEST_EXECUTABLE_SHA256,
                "build_manifest": str(build_manifest),
                "qualification_approval": qualification_approval,
                "input_file": str(submitted_input),
                "input_sha256": analyzer.FIGURE_13_CASE_CONTRACTS[case_id][
                    "input_sha256"
                ],
                "matrix_sha256": analyzer.FIGURE_13_MATRIX_SHA256,
                "time_tlim_target": analyzer.FIGURE_13_ACCEPTED_FINAL_TIME,
                "source_bundle": {
                    "path": str(source_bundle),
                    "sha256": source_bundle_sha256,
                    "verified_revisions": [TEST_REVISION],
                },
                "parent_segment": None,
            },
            "accounting": {
                "case_id": case_id,
                "case_name": case["name"],
                "segment": segment,
                "state": "COMPLETED",
                "exit_code": "0:0",
                "result": "accepted",
                "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
                "input_revision": TEST_REVISION,
                "executable_revision": TEST_REVISION,
                "executable_sha256": TEST_EXECUTABLE_SHA256,
            },
            "scientific_inspection": {
                "case_id": case_id,
                "segment": segment,
                "manifest": str(artifact_path),
                "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
                "accepted": True,
                "clean_for_continuation": True,
                "final_time": analyzer.FIGURE_13_ACCEPTED_FINAL_TIME,
                "terminal_restart_time": analyzer.FIGURE_13_ACCEPTED_FINAL_TIME,
                "terminal_restart": terminal_restart,
                "restart_times": [analyzer.FIGURE_13_ACCEPTED_FINAL_TIME],
                "restarts": [terminal_restart],
                "snapshots": inspected_snapshots,
            },
        }
        artifact_path.write_text(
            json.dumps(artifact, sort_keys=True) + "\n", encoding="utf-8"
        )
        controller_artifacts.append(str(artifact_path))
    manifest = {
        "workflow": analyzer.STAGE_I_PRODUCTION_WORKFLOW,
        "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
        "git_revision": TEST_REVISION,
        "status": "accepted_for_analysis",
        "production_segment_manifests": controller_artifacts,
        "cases": cases,
    }
    bundle.mkdir(parents=True, exist_ok=True)
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    provenance = invocation(analyzer, bundle, tmp_path / "analysis")
    return bundle, cases, manifest, analyzed, provenance


def product_call(analyzer, bundle, manifest, cases, analyzed, provenance):
    return analyzer.figure_13_alternate_firehose_occupancy(
        bundle, manifest, cases, analyzed, provenance
    )


def add_parent_segment(analyzer, bundle, manifest, artifact_index=0):
    """Insert one exact clean-partial parent before an accepted fixture segment."""

    child_path = Path(manifest["production_segment_manifests"][artifact_index])
    child = json.loads(child_path.read_text(encoding="utf-8"))
    case_root = child_path.parents[2]
    parent_segment = "s_parent_t9"
    parent_path = case_root / parent_segment / "manifest" / "prepared_run.json"
    parent_path.parent.mkdir(parents=True, exist_ok=True)
    parent_input = parent_path.parent / "submitted_input.athinput"
    shutil.copyfile(child["command"]["input_file"], parent_input)
    parent_restart_rank_files = []
    expected_ranks = (
        int(child["allocation"]["nodes"])
        * int(child["allocation"]["ranks_per_node"])
    )
    for rank in range(expected_ranks):
        restart = (
            case_root / parent_segment / "output" / "rst" / f"rank_{rank:08d}"
            / f"{child['run']['case_id']}.parent.rst"
        )
        restart.parent.mkdir(parents=True, exist_ok=True)
        restart.write_bytes(athena_restart_payload(9.0))
        parent_restart_rank_files.append({
            "path": str(restart),
            "sha256": analyzer.sha256_file(restart),
            "size_bytes": restart.stat().st_size,
        })
    parent_restart = {
        **parent_restart_rank_files[0],
        "storage": "per_rank",
        "rank_files": parent_restart_rank_files,
    }
    parent = json.loads(json.dumps(child))
    parent["run"]["segment"] = parent_segment
    parent["run"]["run_basename"] = f"{child['run']['case_id']}_{parent_segment}"
    parent["command"]["parent_segment"] = None
    parent["command"]["input_file"] = str(parent_input)
    parent["command"]["time_tlim_target"] = 9.0
    parent["accounting"]["segment"] = parent_segment
    parent["accounting"]["result"] = "clean_partial"
    parent["scientific_inspection"].update({
        "segment": parent_segment,
        "manifest": str(parent_path),
        "accepted": False,
        "clean_for_continuation": True,
        "final_time": 9.0,
        "terminal_restart_time": 9.0,
        "terminal_restart": parent_restart,
        "restart_times": [9.0],
        "restarts": [parent_restart],
    })
    child["command"]["parent_segment"] = {
        "manifest": str(parent_path),
        "execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
        "case_id": child["run"]["case_id"],
        "segment": parent_segment,
        "result": "clean_partial",
        "final_time": 9.0,
        "restart_time": 9.0,
        "input_sha256": child["command"]["input_sha256"],
        "executable_sha256": child["command"]["executable_sha256"],
        "restart_sha256": parent_restart["sha256"],
        "restart_files": [
            record["path"] for record in parent_restart_rank_files
        ],
    }
    archived_restart_rank_files = []
    for rank, source in enumerate(parent_restart_rank_files):
        archived = (
            child_path.parent / "submitted_restart" / f"rank_{rank:08d}"
            / Path(source["path"]).name
        )
        archived.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source["path"], archived)
        archived_restart_rank_files.append({
            "path": str(archived),
            "sha256": analyzer.sha256_file(archived),
            "size_bytes": archived.stat().st_size,
        })
    child["command"].update({
        "source_restart_file": parent_restart["path"],
        "restart_file": archived_restart_rank_files[0]["path"],
        "restart_sha256": parent_restart["sha256"],
        "restart_files": archived_restart_rank_files,
    })
    parent_path.write_text(
        json.dumps(parent, sort_keys=True) + "\n", encoding="utf-8"
    )
    child_path.write_text(
        json.dumps(child, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest["production_segment_manifests"].append(str(parent_path))
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    return child_path


def test_instability_masks_measure_total_and_distinguish_semantics(analyzer):
    values = np.asarray([[[-3.0, -2.0, -1.5, -1.4, 0.0, 1.0, 1.1]]])
    product = analyzer.firehose_threshold_occupancy(occupancy_fields(values))
    strict = product["occupancy"]["strict_manuscript"]
    inclusive = product["occupancy"]["inclusive_solver"]

    assert strict["mirror"]["cell_count"] == 1
    assert strict["parallel_firehose"]["cell_count"] == 1
    assert strict["oblique_firehose"]["cell_count"] == 3
    assert strict["parallel_total_unstable"]["cell_count"] == 2
    assert strict["oblique_total_unstable"]["cell_count"] == 4
    assert strict["oblique_only_reclassification"]["cell_count"] == 2

    assert inclusive["mirror"]["cell_count"] == 2
    assert inclusive["parallel_firehose"]["cell_count"] == 2
    assert inclusive["oblique_firehose"]["cell_count"] == 4
    assert inclusive["parallel_total_unstable"]["cell_count"] == 4
    assert inclusive["oblique_total_unstable"]["cell_count"] == 6
    assert inclusive["oblique_only_reclassification"]["cell_count"] == 2
    definitions = product["threshold_definitions"]["semantics"]
    assert definitions["strict_manuscript"]["parallel_firehose"][
        "comparison_operator"
    ] == "<"
    assert definitions["inclusive_solver"]["parallel_firehose"][
        "comparison_operator"
    ] == "<="


def test_instability_occupancy_fails_closed_on_invalid_fields(analyzer):
    fields = occupancy_fields(np.zeros((1, 1, 2)))
    del fields["p_perp"]
    with pytest.raises(ValueError, match="requires retained snapshot fields"):
        analyzer.firehose_threshold_occupancy(fields)

    fields = occupancy_fields(np.zeros((1, 1, 2)))
    fields["bcc2"][0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="field is nonfinite: bcc2"):
        analyzer.firehose_threshold_occupancy(fields)

    fields = occupancy_fields(np.zeros((1, 1, 2)))
    fields["bcc3"] = np.zeros((1, 2, 2))
    with pytest.raises(ValueError, match="inconsistent grids"):
        analyzer.firehose_threshold_occupancy(fields)

    fields = occupancy_fields(np.zeros((1, 1, 2)))
    fields["bcc1"][0, 0, 0] = 0.0
    with pytest.raises(ValueError, match="finite positive magnetic field strength"):
        analyzer.firehose_threshold_occupancy(fields)

    fields = occupancy_fields(np.zeros((1, 1, 2)))
    fields["bcc1"][0, 0, 0] = np.finfo(float).max
    with pytest.raises(ValueError, match="finite positive magnetic field strength"):
        analyzer.firehose_threshold_occupancy(fields)

    for pressure in ("eint", "p_perp"):
        fields = occupancy_fields(np.zeros((1, 1, 2)))
        fields[pressure][0, 0, 0] = 0.0
        with pytest.raises(ValueError, match="positive p_parallel and p_perp"):
            analyzer.firehose_threshold_occupancy(fields)


def test_trapezoidal_average_is_primary_and_snapshot_mean_secondary(analyzer):
    records = [
        sample_record(analyzer, 8.0, [[[0.0, 0.0]]]),
        sample_record(analyzer, 8.5, [[[-3.0, -3.0]]]),
        sample_record(analyzer, 10.0, [[[-3.0, -3.0]]]),
    ]
    product = analyzer.average_firehose_threshold_occupancy(records)
    parallel = product["occupancy"]["strict_manuscript"][
        "parallel_total_unstable"
    ]

    assert parallel["snapshot_volume_fractions"] == [0.0, 1.0, 1.0]
    assert parallel["comparison_volume_fraction"] == pytest.approx(0.875)
    assert parallel["time_average_volume_fraction"] == pytest.approx(0.875)
    assert parallel["equal_snapshot_mean_volume_fraction"] == pytest.approx(2.0 / 3.0)
    temporal = product["normalization"]["temporal"]
    assert temporal["primary_comparison_kind"] == "trapezoidal_time_average"
    assert temporal["secondary_kind"] == "equal_weight_retained_snapshot_mean"

    inconsistent = [dict(record) for record in records]
    inconsistent[-1] = {**inconsistent[-1], "shape_z_y_x": [1, 2, 2]}
    with pytest.raises(ValueError, match="inconsistent grids"):
        analyzer.average_firehose_threshold_occupancy(inconsistent)


def test_occupancy_averaging_rejects_tampered_thresholds_and_count_fractions(
    analyzer,
):
    records = [
        sample_record(analyzer, 8.0, [[[-3.0, 0.0]]]),
        sample_record(analyzer, 10.0, [[[-3.0, 0.0]]]),
    ]
    records[0]["firehose_threshold_occupancy"]["threshold_definitions"][
        "semantics"
    ]["strict_manuscript"]["parallel_firehose"]["comparison_operator"] = "<="
    with pytest.raises(ValueError, match="threshold definitions do not match"):
        analyzer.average_firehose_threshold_occupancy(records)

    records = [
        sample_record(analyzer, 8.0, [[[-3.0, 0.0]]]),
        sample_record(analyzer, 10.0, [[[-3.0, 0.0]]]),
    ]
    records[0]["firehose_threshold_occupancy"]["occupancy"][
        "strict_manuscript"
    ]["parallel_firehose"]["volume_fraction"] = 0.25
    with pytest.raises(ValueError, match="count/fraction differs"):
        analyzer.average_firehose_threshold_occupancy(records)


def test_complete_product_is_total_unstable_deterministic_and_fully_bound(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )
    repeated = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )

    assert product == repeated
    assert product["complete"] and product["comparison_ready"]
    assert product["available_case_ids"] == ["R03", "R07", "R14", "R15"]
    assert product["comparison_window"]["snapshot_times"] == list(
        analyzer.FIGURE_13_SNAPSHOT_TIMES
    )
    r15 = product["cases"]["R15"]
    assert r15["parallel_total_unstable_volume_fraction"] == 0.5
    assert r15["oblique_total_unstable_volume_fraction"] == 0.75
    assert r15["mirror_volume_fraction"] == 0.25
    assert r15["parallel_firehose_volume_fraction"] == 0.25
    assert r15["oblique_firehose_volume_fraction"] == 0.5
    assert r15["oblique_only_reclassified_volume_fraction"] == 0.25
    assert len(r15["snapshot_provenance"]) == 9
    spatial = r15["normalization"]["spatial"]
    assert spatial["total_cell_count"] == 384 * 192 * 192
    assert spatial["physical_volume"] == 2.0
    assert spatial["cell_volume"] == pytest.approx(2.0 / (384 * 192 * 192))
    assert product["matrix_provenance"]["tracked_matrix_sha256"] == (
        analyzer.FIGURE_13_MATRIX_SHA256
    )
    assert product["analysis_provenance"]["analyzer_sha256"] == analyzer.sha256_file(
        ANALYZER_PATH
    )
    assert product["analysis_provenance"]["bin_convert"]["source_path"] == str(
        (REPO_ROOT / "vis" / "python" / "bin_convert.py").resolve()
    )
    assert product["analysis_provenance"]["bin_convert"]["source_sha256"] == (
        analyzer.sha256_file(REPO_ROOT / "vis" / "python" / "bin_convert.py")
    )
    assert product["execution_authentication"]["available"]
    assert product["execution_authentication"]["authenticated_execution_epoch"] == (
        analyzer.FIGURE_13_EXECUTION_EPOCH
    )
    assert product["execution_authentication"]["authenticated_git_revision"] == (
        TEST_REVISION
    )
    assert product["execution_authentication"]["authenticated_bundle_status"] == (
        "accepted_for_analysis"
    )
    controller = product["execution_authentication"]["controller_artifacts"]["R15"][0]
    assert controller["terminal_restart_loadability"][0]["binary_time"] == 10.0
    assert controller["executable_sha256"] == TEST_EXECUTABLE_SHA256
    assert analyzer.snapshot_time(
        Path(r15["snapshot_provenance"][0]["representative_path"])
    ) == pytest.approx(8.0)
    first_snapshot = r15["snapshot_provenance"][0]
    fields, lengths, snapshot_time = analyzer.read_snapshot(
        Path(first_snapshot["representative_path"]),
        [Path(record["path"]) for record in first_snapshot["files"]],
    )
    raw = analyzer.read_exact_rank_set_binary([
        Path(record["path"]) for record in first_snapshot["files"]
    ])
    assert raw["n_mbs"] == 216
    assert (raw["Nx1"], raw["Nx2"], raw["Nx3"]) == (6, 6, 12)
    assert (raw["nx1_mb"], raw["nx2_mb"], raw["nx3_mb"]) == (1, 1, 2)
    assert fields["eint"].shape == (12, 6, 6)
    assert lengths == (1.0, 1.0, 2.0)
    assert snapshot_time == pytest.approx(8.0)
    assert product["analysis_provenance"]["configuration_sha256"] == (
        analyzer.canonical_json_sha256(
            product["analysis_provenance"]["configuration"]
        )
    )
    assert product["product_sha256"] == analyzer.canonical_json_sha256({
        key: value for key, value in product.items() if key != "product_sha256"
    })
    published = product["published_quantitative_comparison"]
    assert published["case_id"] == "R15"
    assert published["parallel_total_unstable_volume_fraction"] == 0.5
    assert published["oblique_total_unstable_volume_fraction"] == 0.75


@pytest.mark.parametrize(
    "mode", ["missing_cadence", "endpoints", "cli_override", "wrong_grid"]
)
def test_comparison_ready_rejects_inexact_windows_and_cli_overrides(
    analyzer, workflow, tmp_path, mode,
):
    bundle = tmp_path / mode
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    if mode == "missing_cadence":
        times = [
            time for time in analyzer.FIGURE_13_SNAPSHOT_TIMES
            if time != 9.0
        ]
        provenance = invocation(analyzer, bundle, tmp_path / "analysis")
    elif mode == "endpoints":
        times = [8.0, 10.0]
        provenance = invocation(analyzer, bundle, tmp_path / "analysis")
    else:
        times = analyzer.FIGURE_13_SNAPSHOT_TIMES
        provenance = invocation(analyzer, bundle, tmp_path / "analysis")
        if mode == "cli_override":
            provenance = invocation(
                analyzer,
                bundle,
                tmp_path / "analysis",
                time_start=8.0,
                time_end=10.0,
            )
    analyzed = {
        case["name"]: analyzed_case(analyzer, "R15", [[[-3.0, 0.0]]], times)
    }
    if mode == "wrong_grid":
        analyzed[case["name"]]["snapshot_ensemble"][
            "firehose_threshold_occupancy"
        ]["grid"]["shape_z_y_x"] = [1, 1, 2]
    bind_snapshot_outputs(case, analyzed[case["name"]])
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    product = product_call(
        analyzer, bundle, manifest, [case], analyzed, provenance
    )

    assert not product["cases"]["R15"]["available"]
    assert not product["comparison_ready"]
    if mode == "cli_override":
        assert "forbids CLI" in product["cases"]["R15"]["reason"]
    elif mode == "wrong_grid":
        assert "exact 192x192x384" in product["cases"]["R15"]["reason"]
    else:
        assert "exactly cover t=8..10" in product["cases"]["R15"]["reason"]


def test_cadence_accepts_negligible_time_noise_but_rejects_duplicates(
    analyzer, workflow, tmp_path,
):
    bundle = tmp_path / "cadence_tolerance"
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    noise = analyzer.FIGURE_13_TIME_TOLERANCE / 10.0
    times = [
        nominal + (noise if index % 2 else -noise)
        for index, nominal in enumerate(analyzer.FIGURE_13_SNAPSHOT_TIMES)
    ]
    analyzed = {case["name"]: analyzed_case(analyzer, "R15", [[[-3.0, 0.0]]], times)}
    bind_snapshot_outputs(case, analyzed[case["name"]])
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    product = product_call(
        analyzer,
        bundle,
        manifest,
        [case],
        analyzed,
        invocation(analyzer, bundle, tmp_path / "analysis"),
    )
    accepted = product["cases"]["R15"]
    assert accepted["available"]
    assert accepted["analysis_window"]["cadence_validation"][
        "maximum_nominal_time_offset"
    ] <= noise * 1.01

    window = analyzed[case["name"]]["snapshot_ensemble"][
        "firehose_threshold_occupancy"
    ]["analysis_window"]
    window["snapshot_times"][4] = window["snapshot_times"][3]
    product = product_call(
        analyzer,
        bundle,
        manifest,
        [case],
        analyzed,
        invocation(analyzer, bundle, tmp_path / "analysis"),
    )
    assert not product["cases"]["R15"]["available"]
    assert "exactly cover t=8..10" in product["cases"]["R15"]["reason"]


def test_comparison_ready_cases_accept_distinct_negligible_timestamp_noise(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    noise = analyzer.FIGURE_13_TIME_TOLERANCE / 20.0
    for case_index, case in enumerate(cases):
        occupancy = analyzed[case["name"]]["snapshot_ensemble"][
            "firehose_threshold_occupancy"
        ]
        times = [
            nominal + noise * ((index + case_index) % 3 - 1)
            for index, nominal in enumerate(analyzer.FIGURE_13_SNAPSHOT_TIMES)
        ]
        occupancy["analysis_window"].update({
            "selected_time_first": times[0],
            "selected_time_last": times[-1],
            "snapshot_times": times,
        })
        for record, time in zip(occupancy["snapshot_provenance"], times):
            record["snapshot_time"] = time
        occupancy["snapshot_provenance_sha256"] = analyzer.canonical_json_sha256(
            occupancy["snapshot_provenance"]
        )
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )
    assert product["comparison_ready"]
    assert product["available_case_ids"] == list(analyzer.FIGURE_13_FIREHOSE_CASES)


@pytest.mark.parametrize(
    "tamper", ["identity", "model", "execution_input", "case_metadata"]
)
def test_exact_case_identity_and_hash_tampering_fails_closed(
    analyzer, workflow, tmp_path, tamper,
):
    bundle = tmp_path / tamper
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    if tamper == "identity":
        case["input"] = analyzer.FIGURE_13_CASE_CONTRACTS["R14"]["input"]
        match = "exact matrix identity"
    elif tamper == "model":
        case["model_choices"]["passive_delta"] = "true"
        match = "model choices checksum"
    elif tamper == "execution_input":
        (bundle / case["execution_input"]).write_text("tampered\n", encoding="utf-8")
        match = "execution input checksum"
    else:
        case["amr"] = True
        match = "accepted LF case metadata"
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    analyzed = {case["name"]: {"snapshot_ensemble": {"snapshot_count": 0}}}

    with pytest.raises(ValueError, match=match):
        product_call(
            analyzer,
            bundle,
            manifest,
            [case],
            analyzed,
            invocation(analyzer, bundle, tmp_path / "analysis"),
        )


def test_exact_case_rejects_checksum_valid_execution_input_outside_bundle(
    analyzer, workflow, tmp_path,
):
    bundle = tmp_path / "execution_input_escape"
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    outside = tmp_path / "outside.athinput"
    shutil.copyfile(
        REPO_ROOT / analyzer.FIGURE_13_CASE_CONTRACTS["R15"]["input"], outside
    )
    case["execution_input"] = str(outside)
    manifest = write_manifest(bundle, analyzer, [case], "R15")

    with pytest.raises(ValueError, match="execution input path escapes its bundle"):
        product_call(
            analyzer,
            bundle,
            manifest,
            [case],
            {case["name"]: {"snapshot_ensemble": {"snapshot_count": 0}}},
            invocation(analyzer, bundle, tmp_path / "analysis"),
        )


def test_snapshot_provenance_binds_all_rank_siblings_and_detects_mutation(
    analyzer, tmp_path, monkeypatch,
):
    rank0 = tmp_path / "bin" / "rank_00000000"
    rank1 = tmp_path / "bin" / "rank_00000001"
    rank0.mkdir(parents=True)
    rank1.mkdir(parents=True)
    representative = rank0 / "snapshot.bin"
    sibling = rank1 / "snapshot.bin"
    representative.write_bytes(
        (
            "Athena binary output version=1.1\n"
            "  size of preheader=5\n"
            "  time=8.0\n"
            "  cycle=0\n"
            "  size of location=8\n"
            "  size of variable=4\n"
        ).encode()
    )
    sibling_target = tmp_path / "rank-one-source.bin"
    sibling_target.write_bytes(b"rank-one")
    sibling.symlink_to(sibling_target)

    provenance = analyzer.snapshot_digest_provenance(representative)
    assert len(provenance["files"]) == 2
    assert provenance["expected_rank_count"] == 2
    assert provenance["rank_directory_names"] == [
        "rank_00000000",
        "rank_00000001",
    ]
    assert provenance["snapshot_time"] == 8.0
    assert provenance["files"][1]["symlink_target"] == str(sibling_target)
    assert analyzer.validate_snapshot_provenance_record(
        provenance, "unit snapshot"
    ) == provenance
    analyzer.revalidate_snapshot_provenance_record(
        provenance, expected_ranks=2, context="unit snapshot"
    )
    injected = tmp_path / "bin" / "rank_injected"
    injected.mkdir()
    with pytest.raises(ValueError, match="injected rank entries"):
        analyzer.snapshot_digest_provenance(representative, expected_ranks=2)
    injected.rmdir()
    rank2 = tmp_path / "bin" / "rank_00000002"
    rank1.rename(rank2)
    with pytest.raises(ValueError, match="exact expected contiguous set"):
        analyzer.snapshot_digest_provenance(representative, expected_ranks=2)
    rank2.rename(rank1)
    with pytest.raises(ValueError, match="exact expected contiguous set"):
        analyzer.snapshot_digest_provenance(representative, expected_ranks=3)
    missing_rank = tmp_path / "bin" / "rank_00000002"
    missing_rank.mkdir()
    with pytest.raises(ValueError, match="lacks rank sibling"):
        analyzer.snapshot_digest_provenance(representative)
    missing_rank.rmdir()
    sibling.write_bytes(b"rank-one-mutated")
    assert analyzer.snapshot_digest_provenance(representative)[
        "aggregate_sha256"
    ] != provenance["aggregate_sha256"]
    with pytest.raises(ValueError, match="changed before publication"):
        analyzer.revalidate_snapshot_provenance_record(
            provenance, expected_ranks=2, context="unit snapshot"
        )

    sibling.write_bytes(b"stable-before-analysis")
    reads = 0

    exact_sets = []

    def fake_read(_path, rank_files):
        nonlocal reads
        exact_sets.append(list(rank_files))
        reads += 1
        if reads == 2:
            sibling.write_bytes(b"changed-during-analysis")
        return {}, (1.0, 1.0, 1.0), 8.0

    monkeypatch.setattr(analyzer, "read_snapshot", fake_read)
    monkeypatch.setattr(analyzer, "pdf_fields", lambda *_: {})
    monkeypatch.setattr(analyzer, "pressure_density_fields", lambda *_: {})
    monkeypatch.setattr(analyzer, "analyze_fields", lambda *_, **__: {})
    monkeypatch.setattr(
        analyzer, "average_snapshot_records", lambda records: {"count": len(records)}
    )
    with pytest.raises(ValueError, match="changed while being analyzed"):
        analyzer.analyze_snapshot_paths(
            [representative],
            bins=4,
            alignment_shells=[],
            time_start=8.0,
            time_end=8.0,
            expected_ranks_by_path={str(representative): 2},
        )
    assert exact_sets == [[representative, sibling], [representative, sibling]]


def test_exact_rank_reader_preserves_bin_convert_combination_contract(
    analyzer, tmp_path, monkeypatch,
):
    rank_files = [tmp_path / f"rank_{rank:08d}.bin" for rank in range(2)]

    def payload(rank):
        metadata = {
            "header": ["<mesh>", "nx1=2"],
            "time": 8.0,
            "cycle": 10,
            "var_names": ["dens"],
            "Nx1": 2,
            "Nx2": 1,
            "Nx3": 1,
            "nvars": 1,
            "x1min": 0.0,
            "x1max": 1.0,
            "x2min": 0.0,
            "x2max": 1.0,
            "x3min": 0.0,
            "x3max": 1.0,
            "nx1_mb": 1,
            "nx2_mb": 1,
            "nx3_mb": 1,
            "nx1_out_mb": 1,
            "nx2_out_mb": 1,
            "nx3_out_mb": 1,
        }
        return {
            **metadata,
            "n_mbs": 1,
            "mb_index": np.asarray([[rank, rank, 0, 0, 0, 0]]),
            "mb_logical": np.asarray([[rank, 0, 0, 0]]),
            "mb_geometry": np.asarray([
                [0.5 * rank, 0.5 * (rank + 1), 0.0, 1.0, 0.0, 1.0]
            ]),
            "mb_data": {"dens": np.asarray([[[[rank + 1.0]]]])},
        }

    payloads = {str(path): payload(rank) for rank, path in enumerate(rank_files)}
    monkeypatch.setattr(
        analyzer.bin_convert, "read_binary", lambda path: payloads[path]
    )
    combined = analyzer.read_exact_rank_set_binary(rank_files)
    assert combined["n_mbs"] == 2
    assert combined["mb_logical"][:, 0].tolist() == [0, 1]
    assert combined["mb_data"]["dens"].reshape(-1).tolist() == [1.0, 2.0]

    payloads[str(rank_files[1])]["time"] = 8.25
    with pytest.raises(ValueError, match="metadata differs.*time"):
        analyzer.read_exact_rank_set_binary(rank_files)

    payloads[str(rank_files[1])]["time"] = 8.0
    payloads[str(rank_files[1])]["mb_logical"] = np.asarray([[0, 0, 0, 0]])
    with pytest.raises(ValueError, match="duplicate logical meshblocks"):
        analyzer.read_exact_rank_set_binary(rank_files)

    payloads[str(rank_files[1])]["mb_logical"] = np.asarray([[1, 0, 0, 0]])
    with pytest.raises(ValueError, match="uniform meshblock coverage is incomplete"):
        analyzer.read_exact_rank_set_binary(rank_files[:1])

    payloads[str(rank_files[0])]["nx1_mb"] = 0
    payloads[str(rank_files[0])]["nx1_out_mb"] = 0
    with pytest.raises(ValueError, match="meshblock dimensions are invalid"):
        analyzer.read_exact_rank_set_binary(rank_files[:1])


def test_exact_rank_conversion_does_not_mutate_global_bin_convert_reader(
    analyzer, monkeypatch,
):
    raw = {"identity": "exact authenticated rank set"}

    def forbidden_wildcard_reader(_filename):
        raise AssertionError("global wildcard reader must not be called")

    def fake_converter(filename, quantities=None, dtype=None):
        del filename, quantities, dtype
        return read_all_ranks_binary("ignored")  # noqa: F821

    monkeypatch.setattr(
        analyzer.bin_convert, "read_all_ranks_binary", forbidden_wildcard_reader
    )
    monkeypatch.setattr(
        analyzer.bin_convert, "read_all_ranks_binary_as_athdf", fake_converter
    )

    assert analyzer.exact_rank_set_as_athdf(Path("rank0.bin"), raw) is raw
    assert analyzer.bin_convert.read_all_ranks_binary is forbidden_wildcard_reader


def test_shared_snapshot_reader_keeps_historical_one_argument_contract(
    analyzer, tmp_path, monkeypatch,
):
    path = tmp_path / "shared.bin"
    path.write_bytes(
        (
            "Athena binary output version=1.1\n"
            "  size of preheader=5\n"
            "  time=8.0\n"
            "  cycle=0\n"
            "  size of location=8\n"
            "  size of variable=4\n"
        ).encode()
    )
    reads = []

    def fake_read(shared_path):
        reads.append(shared_path)
        return {}, (1.0, 1.0, 1.0), 8.0

    monkeypatch.setattr(analyzer, "read_snapshot", fake_read)
    monkeypatch.setattr(analyzer, "pdf_fields", lambda *_: {})
    monkeypatch.setattr(analyzer, "pressure_density_fields", lambda *_: {})
    monkeypatch.setattr(analyzer, "analyze_fields", lambda *_, **__: {})
    monkeypatch.setattr(
        analyzer, "average_snapshot_records", lambda records: {"count": len(records)}
    )
    records, ensemble = analyzer.analyze_snapshot_paths(
        [path], bins=4, alignment_shells=[], time_start=8.0, time_end=8.0
    )
    assert list(records) == [str(path)]
    assert ensemble["count"] == 1
    assert reads == [path, path]
    provenance = records[str(path)]["snapshot_provenance"]
    assert provenance["snapshot_time"] == 8.0
    assert analyzer.validate_snapshot_provenance_record(
        provenance, "direct-fast shared snapshot"
    ) == provenance
    analyzer.revalidate_snapshot_provenance_record(
        provenance, expected_ranks=1, context="direct-fast shared snapshot"
    )


def test_rank_local_snapshot_reader_rejects_inexact_explicit_rank_set(
    analyzer, tmp_path,
):
    rank0 = tmp_path / "bin" / "rank_00000000"
    rank1 = tmp_path / "bin" / "rank_00000001"
    rank0.mkdir(parents=True)
    rank1.mkdir(parents=True)
    representative = rank0 / "snapshot.bin"
    sibling = rank1 / "snapshot.bin"
    representative.write_bytes(b"rank zero")
    sibling.write_bytes(b"rank one")

    with pytest.raises(ValueError, match="exact expected contiguous set"):
        analyzer.read_snapshot(representative, [representative])


def test_authenticated_bundle_rejects_rank_injection_and_source_substitution(
    analyzer, workflow, tmp_path,
):
    fixture = complete_product_fixture(analyzer, workflow, tmp_path)
    bundle, cases, manifest, _, _ = fixture
    representative = bundle / cases[0]["outputs"]["snapshot_paths"][0]
    injected = representative.parent.parent / "rank_00000002"
    injected.mkdir()
    (injected / representative.name).write_bytes(b"injected")
    with pytest.raises(ValueError, match="exact expected contiguous set"):
        analyzer.figure_13_execution_authentication(bundle, manifest, cases)
    shutil.rmtree(injected)

    rank1 = representative.parent.parent / "rank_00000001" / representative.name
    rank1.unlink()
    substitute = tmp_path / "substitute.bin"
    substitute.write_bytes(b"substitute")
    rank1.symlink_to(substitute)
    with pytest.raises(ValueError, match="rank siblings do not match"):
        analyzer.figure_13_execution_authentication(bundle, manifest, cases)


def test_authenticated_bundle_rejects_snapshot_path_escape(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, _, _ = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    cases[0]["outputs"]["snapshot_paths"][0] = str(
        (bundle / cases[0]["outputs"]["snapshot_paths"][0]).resolve()
    )
    with pytest.raises(ValueError, match="snapshot path escapes its bundle"):
        analyzer.figure_13_execution_authentication(bundle, manifest, cases)


def test_missing_or_tampered_snapshot_provenance_fails_closed(
    analyzer, workflow, tmp_path,
):
    bundle = tmp_path / "provenance"
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    analyzed = {
        case["name"]: analyzed_case(
            analyzer, "R15", [[[-3.0, 0.0]]], include_provenance=False
        )
    }
    product = product_call(
        analyzer,
        bundle,
        manifest,
        [case],
        analyzed,
        invocation(analyzer, bundle, tmp_path / "analysis"),
    )
    assert not product["cases"]["R15"]["available"]
    assert "provenance is required" in product["cases"]["R15"]["reason"]

    analyzed = {case["name"]: analyzed_case(analyzer, "R15", [[[-3.0, 0.0]]])}
    bind_snapshot_outputs(case, analyzed[case["name"]])
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    occupancy = analyzed[case["name"]]["snapshot_ensemble"][
        "firehose_threshold_occupancy"
    ]
    occupancy["snapshot_provenance"][0]["files"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checksum does not match"):
        product_call(
            analyzer,
            bundle,
            manifest,
            [case],
            analyzed,
            invocation(analyzer, bundle, tmp_path / "analysis"),
        )


def test_product_rejects_firehose_only_substitution_for_total_unstable(
    analyzer, workflow, tmp_path,
):
    bundle = tmp_path / "firehose_only"
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    analyzed = {
        case["name"]: analyzed_case(analyzer, "R15", [[[-3.0, 0.0, 2.0]]])
    }
    bind_snapshot_outputs(case, analyzed[case["name"]])
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    strict = analyzed[case["name"]]["snapshot_ensemble"][
        "firehose_threshold_occupancy"
    ]["occupancy"]["strict_manuscript"]
    for threshold in ("parallel", "oblique"):
        strict[f"{threshold}_total_unstable"] = json.loads(json.dumps(
            strict[f"{threshold}_firehose"]
        ))

    with pytest.raises(ValueError, match="total is not mirror OR"):
        product_call(
            analyzer,
            bundle,
            manifest,
            [case],
            analyzed,
            invocation(analyzer, bundle, tmp_path / "analysis"),
        )


def test_product_recomputes_temporal_summaries_from_snapshot_series(
    analyzer, workflow, tmp_path,
):
    bundle = tmp_path / "temporal_summary"
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    analyzed = {
        case["name"]: analyzed_case(analyzer, "R15", [[[-3.0, 0.0]]])
    }
    bind_snapshot_outputs(case, analyzed[case["name"]])
    manifest = write_manifest(bundle, analyzer, [case], "R15")
    strict = analyzed[case["name"]]["snapshot_ensemble"][
        "firehose_threshold_occupancy"
    ]["occupancy"]["strict_manuscript"]
    for name in (
        "parallel_firehose",
        "oblique_firehose",
        "parallel_total_unstable",
        "oblique_total_unstable",
    ):
        strict[name]["comparison_volume_fraction"] = 0.4
        strict[name]["time_average_volume_fraction"] = 0.4

    with pytest.raises(ValueError, match="does not match its snapshot series"):
        product_call(
            analyzer,
            bundle,
            manifest,
            [case],
            analyzed,
            invocation(analyzer, bundle, tmp_path / "analysis"),
        )


def test_snapshot_and_configuration_provenance_change_product_digest(
    analyzer, workflow, tmp_path,
):
    first = complete_product_fixture(analyzer, workflow, tmp_path, salt="first")
    second = complete_product_fixture(analyzer, workflow, tmp_path, salt="second")
    bundle_a, cases_a, manifest_a, analyzed_a, provenance_a = first
    bundle_b, cases_b, manifest_b, analyzed_b, provenance_b = second
    product_a = product_call(
        analyzer, bundle_a, manifest_a, cases_a, analyzed_a, provenance_a
    )
    product_b = product_call(
        analyzer, bundle_b, manifest_b, cases_b, analyzed_b, provenance_b
    )

    assert product_a["cases"]["R15"]["snapshot_provenance_sha256"] != (
        product_b["cases"]["R15"]["snapshot_provenance_sha256"]
    )
    assert product_a["analysis_provenance"]["configuration_sha256"] != (
        product_b["analysis_provenance"]["configuration_sha256"]
    )
    assert product_a["product_sha256"] != product_b["product_sha256"]

    tampered = dict(first[-1])
    tampered["configuration_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="configuration checksum"):
        product_call(
            analyzer, bundle_a, manifest_a, cases_a, analyzed_a, tampered
        )

    (bundle_a / "manifest.json").write_text(
        json.dumps({**manifest_a, "changed": True}, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="manifest changed"):
        product_call(
            analyzer, bundle_a, manifest_a, cases_a, analyzed_a, provenance_a
        )


def test_bin_convert_digest_and_imported_function_identity_fail_closed(
    analyzer, tmp_path, monkeypatch,
):
    provenance = invocation(analyzer, tmp_path / "bundle", tmp_path / "analysis")
    tampered = json.loads(json.dumps(provenance))
    tampered["bin_convert"]["source_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="implementation provenance"):
        analyzer.validate_analysis_invocation_provenance(tampered)

    monkeypatch.setattr(analyzer.bin_convert, "read_binary", lambda *_: {})
    with pytest.raises(ValueError, match="function identity differs: read_binary"):
        analyzer.validate_analysis_invocation_provenance(provenance)


def test_complete_data_without_controller_evidence_is_not_comparison_ready(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    manifest = {
        key: value
        for key, value in manifest.items()
        if key != "production_segment_manifests"
    }
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )

    assert product["complete"]
    assert not product["comparison_ready"]
    assert not product["execution_authentication"]["available"]
    assert any(
        "lacks retained production controller artifacts" in blocker
        for blocker in product["comparison_ready_blockers"]
    )
    assert not product["published_quantitative_comparison"]["available"]


@pytest.mark.parametrize(
    ("field", "value", "match"),
    [
        ("status", "rejected", "epoch, revision, or accepted status"),
        ("execution_epoch", "E99-untrusted", "epoch, revision, or accepted status"),
        ("git_revision", "f" * 40, "revision or input identity"),
    ],
)
def test_bundle_execution_identity_tampering_fails_closed(
    analyzer, workflow, tmp_path, field, value, match,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    manifest[field] = value
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match=match):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)


@pytest.mark.parametrize(
    ("section", "field", "value", "match"),
    [
        ("root", "execution_epoch", "E99-untrusted", "qualifying recorded execution"),
        ("command", "input_revision", "f" * 40, "revision or input identity"),
        ("command", "executable_revision", "f" * 40, "revision or input identity"),
        ("command", "executable_sha256", "0" * 64, "executable size or checksum"),
        ("accounting", "execution_epoch", "E99-untrusted", "execution code identity"),
        ("accounting", "result", "failed", "qualifying recorded execution"),
    ],
)
def test_controller_execution_identity_tampering_fails_closed(
    analyzer, workflow, tmp_path, section, field, value, match,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    artifact_path = Path(manifest["production_segment_manifests"][0])
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    target = artifact if section == "root" else artifact[section]
    target[field] = value
    artifact_path.write_text(
        json.dumps(artifact, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match=match):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)


@pytest.mark.parametrize(
    ("tamper", "match"),
    [
        ("actual_bytes", "size or checksum"),
        ("declared_size", "size or checksum"),
        ("loadability", "loadable"),
        ("duplicate_terminal", "incomplete or ambiguous|not unique"),
        ("duplicate_rank_path", "incomplete or ambiguous"),
    ],
)
def test_controller_terminal_restart_bytes_uniqueness_and_loadability(
    analyzer, workflow, tmp_path, tamper, match,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    artifact_path = Path(manifest["production_segment_manifests"][0])
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    inspection = artifact["scientific_inspection"]
    terminal = inspection["terminal_restart"]
    retained_terminal = inspection["restarts"][0]
    rank_path = Path(terminal["rank_files"][0]["path"])
    if tamper == "actual_bytes":
        rank_path.write_bytes(rank_path.read_bytes() + b"changed")
    elif tamper == "declared_size":
        for record in (terminal["rank_files"][0], retained_terminal["rank_files"][0]):
            record["size_bytes"] += 1
        terminal["size_bytes"] += 1
        retained_terminal["size_bytes"] += 1
    elif tamper == "loadability":
        rank_path.write_bytes(b"not a loadable Athena restart")
        digest = analyzer.sha256_file(rank_path)
        for record in (terminal["rank_files"][0], retained_terminal["rank_files"][0]):
            record["sha256"] = digest
            record["size_bytes"] = rank_path.stat().st_size
        for record in (terminal, retained_terminal):
            record["sha256"] = digest
            record["size_bytes"] = rank_path.stat().st_size
    elif tamper == "duplicate_terminal":
        inspection["restarts"].append(json.loads(json.dumps(terminal)))
        inspection["restart_times"].append(inspection["terminal_restart_time"])
    else:
        duplicate = json.loads(json.dumps(terminal["rank_files"][0]))
        terminal["rank_files"][1] = duplicate
        retained_terminal["rank_files"][1] = duplicate
    artifact_path.write_text(
        json.dumps(artifact, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match=match):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)


def test_real_restart_contract_rejects_fabricated_payload_and_parameter_dump(
    analyzer, tmp_path,
):
    contract = analyzer.json_compatible(next(iter(
        analyzer.FIGURE_13_QUALIFIED_RESTART_BINARY_ABIS.values()
    )))
    abi = {
        "executable_revision": TEST_REVISION,
        "executable_sha256": TEST_EXECUTABLE_SHA256,
        "contract": contract,
        "contract_sha256": analyzer.canonical_json_sha256(contract),
    }
    blocks = analyzer.stage_i_validator.parse_athinput(
        (
            REPO_ROOT
            / analyzer.FIGURE_13_CASE_CONTRACTS["R03"]["input"]
        ).read_bytes(),
        "test R03 input",
    )
    input_contract = analyzer.stage_i_validator.require_input_contract(
        blocks,
        "R03",
        {
            "name": analyzer.FIGURE_13_CASE_CONTRACTS["R03"]["name"],
            "input": analyzer.FIGURE_13_CASE_CONTRACTS["R03"]["input"],
            "resolution": analyzer.FIGURE_13_CASE_CONTRACTS["R03"]["resolution"],
        },
    )
    assert input_contract["restart_data_size"] == 6_394_752
    target = tmp_path / "fabricated.rst"
    target.write_bytes(athena_restart_payload(10.0))
    with pytest.raises(
        ValueError,
        match="exact production loadability failed: restart parameter dump "
        "parameter inventory differs from qualified contract",
    ):
        analyzer.real_restart_loadability_evidence(
            target,
            10.0,
            abi,
            input_contract,
            {"run_basename": "test", "target_time": 10.0},
        )


def test_real_restart_contract_accepts_retained_production_rank_set(
    analyzer, tmp_path,
):
    manifest_path = Path(
        "/lustre/orion/ast207/proj-shared/dfielding/CGL/runs/mks24-stage-i/"
        "E03-forcing-policy/R03/s00_rankio_t0_t0p5/manifest/prepared_run.json"
    )
    if not manifest_path.is_file():
        pytest.skip("retained production R03 restart rank set is unavailable")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    command = manifest["command"]
    run = manifest["run"]
    blocks = analyzer.stage_i_validator.parse_athinput(
        Path(command["input_file"]).read_bytes(), "retained production R03 input"
    )
    input_contract = analyzer.stage_i_validator.require_input_contract(
        blocks,
        "R03",
        {
            "name": run["case_name"],
            "input": analyzer.FIGURE_13_CASE_CONTRACTS["R03"]["input"],
            "resolution": run["resolution"],
        },
    )
    abi = analyzer.qualified_restart_binary_abi(
        command["executable_revision"],
        command["executable_sha256"],
        "retained production R03 restart",
    )
    rank_files = sorted(
        manifest_path.parents[1].glob("output/rst/rank_*/*.00001.rst")
    )
    if len(rank_files) != 8:
        pytest.skip("retained production R03 restart rank set is incomplete")
    evidence = [
        analyzer.real_restart_loadability_evidence(
            path,
            0.31282347945569927,
            abi,
            {
                **input_contract,
            },
            {
                "run_basename": run["run_basename"],
                "target_time": command["time_tlim_target"],
            },
        )
        for path in rank_files
    ]
    assert sum(item["local_meshblock_count"] for item in evidence) == 216
    assert {item["meshblock_count"] for item in evidence} == {216}
    assert {item["variable_data_size"] for item in evidence} == {6_394_752}
    assert {item["mode_count"] for item in evidence} == {11}
    assert {
        item["turbulence_metadata_contract"]["validation"] for item in evidence
    } == {"all-runtime-configuration-fields-exact"}

    source = rank_files[0]
    _, parameter_size = analyzer.restart_parameter_dump(source)
    contract = abi["contract"]
    metadata_offset = (
        parameter_size
        + int(contract["mesh_header_size"])
        + 216 * int(contract["logical_location_size"])
        + 216 * int(contract["cost_size"])
    )
    with source.open("rb") as stream:
        stream.seek(metadata_offset + 4)
        mode_count = struct.unpack("<i", stream.read(4))[0]
        stream.seek(0)
        prefix = stream.read(
            metadata_offset
            + int(contract["turbulence_metadata_size"])
            + int(contract["rng_state_size"])
            + 6 * mode_count * 8
            + 8
            + int(contract["lf_diagnostic_count"]) * 8
        )

    fabricated = tmp_path / "fabricated-data-size.rst"
    fabricated.write_bytes(prefix + struct.pack("<Q", 16) + b"\0" * 16)
    with pytest.raises(
        ValueError,
        match="variable-data size differs from archived input",
    ):
        analyzer.real_restart_loadability_evidence(
            fabricated,
            0.31282347945569927,
            abi,
            input_contract,
            {
                "run_basename": run["run_basename"],
                "target_time": command["time_tlim_target"],
            },
        )

    metadata_tamper = bytearray(prefix)
    struct.pack_into("<i", metadata_tamper, metadata_offset + 5 * 4, 0)
    tampered = tmp_path / "tampered-turbulence-metadata.rst"
    tampered.write_bytes(metadata_tamper)
    with pytest.raises(
        ValueError,
        match="turbulence metadata differs from exact runtime configuration",
    ):
        analyzer.validate_restart_turbulence_metadata(
            tampered, abi, input_contract
        )


@pytest.mark.parametrize(
    "tamper", ["bytes", "path", "approval", "restart_load_evidence"]
)
def test_controller_authenticates_live_executable_against_approval(
    analyzer, workflow, tmp_path, tamper,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    artifact_path = Path(manifest["production_segment_manifests"][0])
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    executable = Path(artifact["command"]["executable"])
    if tamper == "bytes":
        executable.write_bytes(b"different executable")
        match = "executable size or checksum"
    elif tamper == "path":
        substitute = executable.parent / "substitute-athena"
        shutil.copyfile(executable, substitute)
        substitute.chmod(0o755)
        artifact["command"]["executable"] = str(substitute)
        artifact_path.write_text(
            json.dumps(artifact, sort_keys=True) + "\n", encoding="utf-8"
        )
        match = "qualification approval does not bind executable"
    elif tamper == "approval":
        approval = Path(artifact["command"]["qualification_approval"]["path"])
        approval.write_text(
            approval.read_text(encoding="utf-8") + "\n", encoding="utf-8"
        )
        match = "qualification approval does not bind executable"
    else:
        evidence = (
            analyzer.STAGE_I_CANONICAL_ROOT
            / analyzer.FIGURE_13_QUALIFIED_RESTART_LOAD_EVIDENCE_RELATIVE
        )
        evidence.write_text(
            evidence.read_text(encoding="utf-8") + "\n", encoding="utf-8"
        )
        match = "qualified restart-load evidence does not bind"
    with pytest.raises(ValueError, match=match):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("parent_segment", "execution_epoch", "E99-untrusted"),
        ("parent_segment", "input_sha256", "0" * 64),
        ("parent_segment", "executable_sha256", "0" * 64),
        ("parent_segment", "restart_sha256", "0" * 64),
        ("parent_segment", "restart_files", []),
        ("parent_segment", "restart_time", 8.5),
        ("command", "source_restart_file", "/tmp/unbound.rst"),
        ("command", "restart_file", "/tmp/unbound.rst"),
        ("command", "restart_sha256", "0" * 64),
        ("command", "restart_files", []),
    ],
)
def test_controller_lineage_binds_exact_parent_restart_and_execution_identity(
    analyzer, workflow, tmp_path, section, field, value,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    child_path = add_parent_segment(analyzer, bundle, manifest)
    assert product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )["comparison_ready"]

    child = json.loads(child_path.read_text(encoding="utf-8"))
    target = child["command"]["parent_segment"] if section == "parent_segment" else (
        child["command"]
    )
    target[field] = value
    child_path.write_text(
        json.dumps(child, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="unauthenticated parent"):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)


def test_controller_lineage_rejects_unbound_root_restart(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    artifact_path = Path(manifest["production_segment_manifests"][0])
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    artifact["command"]["restart_sha256"] = "0" * 64
    artifact_path.write_text(
        json.dumps(artifact, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="unauthenticated root restart"):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)


def test_snapshot_mutation_after_analysis_blocks_ready_product_and_publication(
    analyzer, workflow, tmp_path,
):
    fixture = complete_product_fixture(analyzer, workflow, tmp_path, salt="product")
    bundle, cases, manifest, analyzed, provenance = fixture
    rank_file = Path(
        analyzed[cases[0]["name"]]["snapshot_ensemble"][
            "firehose_threshold_occupancy"
        ]["snapshot_provenance"][0]["files"][1]["path"]
    )
    rank_file.write_bytes(b"changed after analysis")
    with pytest.raises(ValueError, match="snapshot bytes changed before publication"):
        product_call(analyzer, bundle, manifest, cases, analyzed, provenance)

    fixture = complete_product_fixture(analyzer, workflow, tmp_path, salt="publication")
    bundle, cases, manifest, analyzed, provenance = fixture
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )
    rank_file = Path(
        product["cases"]["R15"]["snapshot_provenance"][0]["files"][0]["path"]
    )
    rank_file.write_bytes(b"changed immediately before publication")
    output = tmp_path / "publication-output"
    output.mkdir()
    named = output / analyzer.FIGURE_13_FIREHOSE_PRODUCT_NAME
    named.write_bytes(b"previous qualified product\n")
    with pytest.raises(ValueError, match="snapshot bytes changed before publication"):
        analyzer.write_figure_13_alternate_firehose_occupancy(product, output)
    assert named.read_bytes() == b"previous qualified product\n"


@pytest.mark.parametrize(
    "artifact_kind",
    [
        "controller",
        "source_bundle",
        "executable",
        "approval",
        "restart_qualification",
        "terminal_restart",
    ],
)
def test_final_publication_revalidates_retained_execution_artifacts(
    analyzer, workflow, tmp_path, artifact_kind,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )
    controller = product["execution_authentication"]["controller_artifacts"][
        "R03"
    ][0]
    if artifact_kind == "controller":
        path = Path(controller["path"])
        path.write_text(
            path.read_text(encoding="utf-8") + "\n", encoding="utf-8"
        )
        match = "controller artifact changed"
    else:
        if artifact_kind == "source_bundle":
            path = Path(controller["source_bundle_path"])
            path.write_bytes(b"changed source bundle")
            match = "source bundle changed"
        elif artifact_kind == "executable":
            path = Path(controller["executable_path"])
            path.write_bytes(b"changed executable")
            match = "executable size or checksum"
        elif artifact_kind == "approval":
            path = Path(controller["qualification_approval_path"])
            path.write_text(
                path.read_text(encoding="utf-8") + "\n", encoding="utf-8"
            )
            match = "executable approval changed"
        elif artifact_kind == "restart_qualification":
            path = Path(controller["qualified_restart_load_evidence"]["path"])
            path.write_text(
                path.read_text(encoding="utf-8") + "\n", encoding="utf-8"
            )
            match = "qualified restart-load evidence does not bind"
        else:
            path = Path(controller["terminal_restart_rank_files"][0]["path"])
            path.write_bytes(path.read_bytes() + b"changed")
            match = "size or checksum"
    with pytest.raises(ValueError, match=match):
        analyzer.revalidate_figure_13_product_for_publication(product)


def test_final_publication_revalidates_child_restart_archive(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    add_parent_segment(analyzer, bundle, manifest)
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )
    child = product["execution_authentication"]["controller_artifacts"]["R03"][0]
    archived = Path(child["restart_archive_rank_files"][0]["path"])
    archived.write_bytes(archived.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="size or checksum"):
        analyzer.revalidate_figure_13_product_for_publication(product)


def test_final_publication_rejects_product_mutation(
    analyzer, workflow, tmp_path,
):
    bundle, cases, manifest, analyzed, provenance = complete_product_fixture(
        analyzer, workflow, tmp_path
    )
    product = product_call(
        analyzer, bundle, manifest, cases, analyzed, provenance
    )
    product["cases"]["R15"]["parallel_total_unstable_volume_fraction"] = 0.0
    with pytest.raises(ValueError, match="product checksum changed"):
        analyzer.revalidate_figure_13_product_for_publication(product)


def test_publication_generation_is_crash_atomic_for_intermediate_readers(
    analyzer, tmp_path, monkeypatch,
):
    output = tmp_path / "transaction"
    old = {
        "generation": "old diagnostics",
        "figure_13_alternate_firehose_occupancy": {
            "comparison_ready": False,
            "generation": "old product",
        },
    }
    analyzer.write_figure_13_publication_generation(old, output)
    old_selected = analyzer.read_analysis_publication(output)
    old_generation = old_selected["generation_dir"]
    real_replace = analyzer.os.replace
    observations = []

    def crash_after_pointer_switch(source, destination):
        if Path(destination).name == analyzer.ANALYSIS_PUBLICATION_CURRENT_NAME:
            observations.append(
                analyzer.read_analysis_publication(output)["generation_dir"]
            )
            real_replace(source, destination)
            observations.append(
                analyzer.read_analysis_publication(output)["generation_dir"]
            )
            raise OSError("injected crash after atomic pointer switch")
        return real_replace(source, destination)

    monkeypatch.setattr(analyzer.os, "replace", crash_after_pointer_switch)
    result = {
        "generation": "new diagnostics",
        "figure_13_alternate_firehose_occupancy": {
            "comparison_ready": False,
            "generation": "new product",
        },
    }
    with pytest.raises(OSError, match="injected crash after atomic pointer switch"):
        analyzer.write_figure_13_publication_generation(result, output)

    selected = analyzer.read_analysis_publication(output)
    assert observations[0] == old_generation
    assert observations[1] == selected["generation_dir"]
    assert selected["generation_dir"] != old_generation
    assert selected["manifest"]["mode"] == "figure-13"
    assert selected["diagnostics"].stat().st_mode & 0o777 == 0o444
    assert selected["named_product"].stat().st_mode & 0o777 == 0o444
    assert selected["generation_dir"].stat().st_mode & 0o777 == 0o555
    assert not list(output.glob(f".{analyzer.ANALYSIS_PUBLICATION_CURRENT_NAME}.*"))


def test_publication_crash_before_pointer_switch_preserves_old_generation(
    analyzer, tmp_path, monkeypatch,
):
    output = tmp_path / "pre-switch-crash"
    analyzer.write_analysis_publication_generation({"generation": "old"}, output)
    old = analyzer.read_analysis_publication(output)
    old_bytes = old["diagnostics"].read_bytes()
    real_replace = analyzer.os.replace

    def crash_before_pointer_switch(source, destination):
        if Path(destination).name == analyzer.ANALYSIS_PUBLICATION_CURRENT_NAME:
            raise OSError("injected crash before atomic pointer switch")
        return real_replace(source, destination)

    monkeypatch.setattr(analyzer.os, "replace", crash_before_pointer_switch)
    with pytest.raises(OSError, match="injected crash before atomic pointer switch"):
        analyzer.write_analysis_publication_generation({"generation": "new"}, output)

    selected = analyzer.read_analysis_publication(output)
    assert selected["generation_dir"] == old["generation_dir"]
    assert selected["diagnostics"].read_bytes() == old_bytes


def test_non_figure_publication_atomically_withdraws_named_product(
    analyzer, tmp_path, monkeypatch,
):
    output = tmp_path / "mode-switch"
    analyzer.write_figure_13_publication_generation({
        "figure_13_alternate_firehose_occupancy": {
            "comparison_ready": False,
            "generation": "figure product",
        },
    }, output)
    named = output / analyzer.FIGURE_13_FIREHOSE_PRODUCT_NAME
    assert named.is_file()
    old_generation = analyzer.read_analysis_publication(output)["generation_dir"]

    visibility = []
    real_replace = analyzer.os.replace

    def observe_pointer_switch(source, destination):
        if Path(destination).name == analyzer.ANALYSIS_PUBLICATION_CURRENT_NAME:
            visibility.append(named.is_file())
            real_replace(source, destination)
            visibility.append(named.is_file())
            return None
        return real_replace(source, destination)

    monkeypatch.setattr(analyzer.os, "replace", observe_pointer_switch)
    analyzer.write_analysis_publication_generation(
        {"generation": "diagnostics only"}, output
    )
    selected = analyzer.read_analysis_publication(output)
    assert selected["manifest"]["mode"] == "diagnostics-only"
    assert selected["named_product"] is None
    assert not named.exists()
    assert named.is_symlink()
    assert old_generation.is_dir()
    assert visibility == [True, False]


def test_publication_refuses_non_atomic_legacy_alias_migration(
    analyzer, tmp_path,
):
    output = tmp_path / "legacy-publication"
    output.mkdir()
    diagnostics = output / "diagnostics.json"
    named = output / analyzer.FIGURE_13_FIREHOSE_PRODUCT_NAME
    diagnostics.write_bytes(b'{"legacy":"diagnostics"}\n')
    named.write_bytes(b'{"legacy":"named product"}\n')
    before = (diagnostics.read_bytes(), named.read_bytes())

    with pytest.raises(ValueError, match="refuses non-atomic legacy alias migration"):
        analyzer.write_analysis_publication_generation(
            {"generation": "new diagnostics-only generation"}, output
        )

    assert (diagnostics.read_bytes(), named.read_bytes()) == before
    assert not (output / analyzer.ANALYSIS_PUBLICATION_CURRENT_NAME).exists()


def test_publication_reader_rejects_generation_profile_mutation(
    analyzer, tmp_path,
):
    output = tmp_path / "profile-mutation"
    analyzer.write_analysis_publication_generation({"generation": "stable"}, output)
    selected = analyzer.read_analysis_publication(output)
    selected["diagnostics"].chmod(0o644)
    with pytest.raises(ValueError, match="member profile differs|modes differ"):
        analyzer.read_analysis_publication(output)


def test_cli_atomic_named_product_update_failure_preservation_and_mode_withdrawal(
    analyzer, workflow, tmp_path,
):
    output = tmp_path / "analysis"
    bundle = tmp_path / "figure13_bundle"
    case = figure_13_case(analyzer, workflow, bundle, "R15")
    snapshot_paths = []
    for index, time in enumerate(analyzer.FIGURE_13_SNAPSHOT_TIMES):
        filename = f"R15.out2.{index:05d}.bin"
        for rank in range(2):
            path = (
                bundle / "cases" / case["name"] / "bin"
                / f"rank_{rank:08d}" / filename
            )
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(athena_snapshot_payload(time, index, rank))
        snapshot_paths.append(str(
            (
                Path("cases") / case["name"] / "bin" / "rank_00000000"
                / filename
            )
        ))
    case["outputs"]["snapshot_paths"] = snapshot_paths
    write_manifest(bundle, analyzer, [case], "R15")
    command = [
        sys.executable,
        str(ANALYZER_PATH),
        "--bundle",
        str(bundle),
        "--output-dir",
        str(output),
    ]
    result = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    named = output / analyzer.FIGURE_13_FIREHOSE_PRODUCT_NAME
    assert named.is_file()
    product = json.loads(named.read_text(encoding="utf-8"))
    assert not product["comparison_ready"]
    diagnostics = json.loads(
        (output / "diagnostics.json").read_text(encoding="utf-8")
    )
    assert product == diagnostics["figure_13_alternate_firehose_occupancy"]
    assert product["publication_generation_sha256"] == (
        diagnostics["figure_13_publication_generation_sha256"]
    )
    assert product["publication_generation_sha256"] == (
        analyzer.figure_13_publication_generation_sha256(diagnostics)
    )
    assert product["product_sha256"] == analyzer.canonical_json_sha256({
        key: value for key, value in product.items() if key != "product_sha256"
    })
    reconstructed = diagnostics["cases"][case["name"]]["snapshot_ensemble"][
        "firehose_threshold_occupancy"
    ]
    assert reconstructed["grid"]["shape_z_y_x"] == [12, 6, 6]
    assert reconstructed["normalization"]["spatial"]["total_cell_count"] == 432
    assert len(diagnostics["cases"][case["name"]]["snapshots"]) == 9
    assert "exact 192x192x384" in product["cases"]["R15"]["reason"]
    retained_product = named.read_bytes()
    retained_generation = analyzer.read_analysis_publication(output)["generation_dir"]

    result = subprocess.run(
        [
            sys.executable,
            str(ANALYZER_PATH),
            "--bundle",
            str(tmp_path / "missing_bundle"),
            "--output-dir",
            str(output),
        ],
        cwd=REPO_ROOT,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert named.read_bytes() == retained_product
    assert analyzer.read_analysis_publication(output)["generation_dir"] == (
        retained_generation
    )

    other_bundle = tmp_path / "other_bundle"
    write_manifest(other_bundle, analyzer, [{
        "case_id": "R02",
        "name": "not_a_figure_13_case",
        "outputs": {},
    }])
    result = subprocess.run(
        [
            sys.executable,
            str(ANALYZER_PATH),
            "--bundle",
            str(other_bundle),
            "--output-dir",
            str(output),
        ],
        cwd=REPO_ROOT,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    selected = analyzer.read_analysis_publication(output)
    assert selected["manifest"]["mode"] == "diagnostics-only"
    assert selected["named_product"] is None
    assert not named.exists()
    assert named.is_symlink()
    assert retained_generation.is_dir()


def test_four_case_comparison_ready_cli_reconstructs_production_shaped_snapshots(
    analyzer, workflow, tmp_path, monkeypatch,
):
    retained_root = Path(
        "/lustre/orion/ast207/proj-shared/dfielding/CGL/runs/mks24-stage-i/"
        "E03-forcing-policy/R02"
    )
    if not retained_root.is_dir():
        pytest.skip("retained production-shaped Athena snapshots are unavailable")
    representatives = {}
    for path in sorted(retained_root.glob("*/output/bin/rank_00000000/*.bin")):
        time = analyzer.snapshot_time(path)
        for nominal in analyzer.FIGURE_13_SNAPSHOT_TIMES:
            if analyzer.figure_13_time_close(time, nominal):
                representatives.setdefault(nominal, path)
    if 9.5 not in representatives and 9.25 in representatives:
        source = representatives[9.25]
        source_ranks = analyzer.snapshot_sibling_paths(source, expected_ranks=8)
        gap_root = tmp_path / "production-shaped-gap-fill"
        for rank, source_rank in enumerate(source_ranks):
            destination = gap_root / f"rank_{rank:08d}" / "gap-fill.bin"
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_rank, destination)
            with destination.open("r+b") as stream:
                preheader = stream.read(512)
                old = b"time=9.2500000000000000e+00"
                new = b"time=9.5000000000000000e+00"
                assert preheader.count(old) == 1 and len(old) == len(new)
                stream.seek(0)
                stream.write(preheader.replace(old, new, 1))
        representatives[9.5] = gap_root / "rank_00000000" / "gap-fill.bin"
    assert set(representatives) == set(analyzer.FIGURE_13_SNAPSHOT_TIMES)

    bundle = tmp_path / "production-shaped-four-case-bundle"
    cases = [
        figure_13_case(analyzer, workflow, bundle, case_id)
        for case_id in analyzer.FIGURE_13_FIREHOSE_CASES
    ]
    ordered_paths = [
        representatives[time] for time in analyzer.FIGURE_13_SNAPSHOT_TIMES
    ]
    for case in cases:
        case["outputs"]["snapshot_paths"] = [str(path) for path in ordered_paths]
    manifest = write_manifest(bundle, analyzer, cases)
    manifest.update({
        "status": "accepted_for_analysis",
        "git_revision": TEST_REVISION,
    })
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
    )

    provenance = {
        str(path): analyzer.snapshot_digest_provenance(path, expected_ranks=8)
        for path in ordered_paths
    }
    authentication = {
        "available": True,
        "canonical_root": str(analyzer.STAGE_I_CANONICAL_ROOT),
        "required_execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
        "required_final_time": analyzer.FIGURE_13_ACCEPTED_FINAL_TIME,
        "authenticated_execution_epoch": analyzer.FIGURE_13_EXECUTION_EPOCH,
        "authenticated_git_revision": TEST_REVISION,
        "authenticated_bundle_status": "accepted_for_analysis",
        "snapshot_expected_ranks": {str(path): 8 for path in ordered_paths},
        "snapshot_controller_inventory": {
            str(path): {
                "expected_rank_count": 8,
                "rank_files": [
                    {
                        "path": record["path"],
                        "sha256": record["sha256"],
                        "size_bytes": record["size_bytes"],
                    }
                    for record in provenance[str(path)]["files"]
                ],
            }
            for path in ordered_paths
        },
        "controller_artifacts": {},
    }
    monkeypatch.setattr(
        analyzer,
        "figure_13_execution_authentication",
        lambda bundle_path, bundle_manifest, bundle_cases: authentication,
    )
    output = tmp_path / "production-shaped-cli-output"
    monkeypatch.setattr(sys, "argv", [
        str(ANALYZER_PATH),
        "--bundle", str(bundle),
        "--output-dir", str(output),
        "--figure-13-only",
    ])
    assert analyzer.main() == 0

    selected = analyzer.read_analysis_publication(output)
    diagnostics = json.loads(selected["diagnostics"].read_text(encoding="utf-8"))
    product = json.loads(selected["named_product"].read_text(encoding="utf-8"))
    assert selected["manifest"]["mode"] == "figure-13"
    assert product["complete"] and product["comparison_ready"]
    assert product["available_case_ids"] == ["R03", "R07", "R14", "R15"]
    assert product["comparison_grid"]["shape_z_y_x"] == [384, 192, 192]
    assert product["cases"]["R15"]["normalization"]["spatial"][
        "total_cell_count"
    ] == 384 * 192 * 192
    assert all(
        diagnostics["cases"][case["name"]]["snapshot_ensemble"][
            "firehose_threshold_occupancy"
        ]["grid"]["shape_z_y_x"] == [384, 192, 192]
        for case in cases
    )
