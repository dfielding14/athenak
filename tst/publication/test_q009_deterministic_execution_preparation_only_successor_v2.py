#!/usr/bin/env python3
"""Focused tests for the deterministic, non-authorizing Q009 preparation."""

from __future__ import annotations

import ast
import builtins
import copy
import hashlib
import inspect
import json
from pathlib import Path
import stat
import subprocess
import sys
import types

import pytest

from tst.publication import (
    q009_deterministic_execution_preparation_only_successor_v2 as preparation,
)


READINESS_PATH = (
    Path(__file__).resolve().parent
    / "readiness/q009_deterministic_execution_preparation_only_successor_v2_2026-06-06.json"
)
LIVE_ORION_POLICY_PATH = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/PIC/policy/storage_policy.json"
)


def _json(payload: bytes) -> dict[str, object]:
    value = json.loads(payload)
    assert isinstance(value, dict)
    return value


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("utf-8")


def _deck_blocks(text: str) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    current: str | None = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            current = line[1:-1].strip()
            blocks[current] = {}
            continue
        assert current is not None
        name, separator, value = line.partition("=")
        assert separator == "="
        blocks[current][name.strip()] = value.strip()
    return blocks


def _materialization() -> tuple[dict[str, object], list[dict[str, object]], dict[str, bytes]]:
    manifest, files = preparation.build_materialization()
    launches = [
        _json(files[record["launch_candidate"]["path"]])
        for record in manifest["stage_records"]
    ]
    return manifest, launches, files


def _test_private_module(name: str, path: Path, payload: bytes) -> types.ModuleType:
    module = types.ModuleType(name)
    module.__file__ = str(path)
    exec(compile(payload, module.__file__, "exec", dont_inherit=True), module.__dict__)
    return module


def _fresh_process(script: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", script],
        cwd=preparation.REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def _source_tree() -> ast.Module:
    return ast.parse(Path(preparation.__file__).read_text(encoding="utf-8"))


def _module_callables() -> set[str]:
    return {
        name
        for name, value in vars(preparation).items()
        if inspect.isfunction(value) and value.__module__ == preparation.__name__
    }


def _minimal_callable(public: object, name: str) -> object:
    assert inspect.isfunction(public)
    for cell in public.__closure__ or ():
        value = cell.cell_contents
        if (
            inspect.isfunction(value)
            and value.__module__ == preparation.MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID
            and value.__name__ == name
        ):
            return value
    raise AssertionError(f"missing isolated minimal callable: {name}")


def _closure_callable(public: object, name: str) -> object:
    pending = [public]
    visited: set[int] = set()
    while pending:
        function = pending.pop()
        if id(function) in visited or not inspect.isfunction(function):
            continue
        visited.add(id(function))
        if function.__name__ == name:
            return function
        for cell in function.__closure__ or ():
            value = cell.cell_contents
            if inspect.isfunction(value):
                pending.append(value)
    raise AssertionError(f"missing closure callable: {name}")


def _noop_code_like(function: object) -> object:
    assert inspect.isfunction(function)
    names = [f"value_{index}" for index in range(len(function.__code__.co_freevars))]
    parameters = ", ".join(names)
    source = f"def outer({parameters}):\n    def inner():\n"
    if names:
        source += "        _ = (" + ", ".join(names) + ",)\n"
    source += "        return None\n    return inner\n"
    namespace: dict[str, object] = {}
    exec(source, namespace)
    return namespace["outer"](*([None] * len(names))).__code__


def _filesystem_mutator_calls() -> set[str]:
    forbidden = {
        "chmod",
        "fchmod",
        "fsync",
        "link",
        "mkdir",
        "rename",
        "rmdir",
        "symlink",
        "unlink",
        "write",
    }
    return {
        node.func.attr
        for node in ast.walk(_source_tree())
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in forbidden
    }


def _live_storage_policy_validator() -> tuple[object, dict[str, object], str]:
    policy_payload = LIVE_ORION_POLICY_PATH.read_bytes()
    policy = json.loads(policy_payload)
    project_root = Path(policy["olcf_side_storage"]["project_home_mirror_root"])
    assert policy_payload == (project_root / "policy/storage_policy.json").read_bytes()
    version = policy["olcf_side_storage"]["installed_control_plane_version"]
    assert isinstance(version, str)
    roots = [
        LIVE_ORION_POLICY_PATH.parents[1] / "control_plane" / version,
        project_root / "control_plane" / version,
    ]
    installed_payloads: dict[str, bytes] | None = None
    required = (
        "control_plane_common.py",
        "operator_attestation.py",
        "q011_pressure_review_packet_verifier.py",
    )
    for root in roots:
        assert root.resolve(strict=True) == root
        assert not root.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
        inventory = json.loads((root / "inventory.json").read_text(encoding="utf-8"))
        assert inventory["version"] == version
        bindings = {record["path"]: record["sha256"] for record in inventory["files"]}
        payloads = {name: (root / name).read_bytes() for name in required}
        for name, payload in payloads.items():
            assert hashlib.sha256(payload).hexdigest() == bindings[name]
        if installed_payloads is not None:
            assert payloads == installed_payloads
        installed_payloads = payloads
    assert installed_payloads is not None
    operator = _test_private_module(
        "_q009_test_live_operator_attestation",
        roots[0] / "operator_attestation.py",
        installed_payloads["operator_attestation.py"],
    )
    pressure = _test_private_module(
        "_q009_test_live_pressure_verifier",
        roots[0] / "q011_pressure_review_packet_verifier.py",
        installed_payloads["q011_pressure_review_packet_verifier.py"],
    )
    aliases = {
        "operator_attestation": operator,
        "q011_pressure_review_packet_verifier": pressure,
    }
    absent = object()
    previous = {name: sys.modules.get(name, absent) for name in aliases}
    sys.modules.update(aliases)
    try:
        common = _test_private_module(
            "_q009_test_live_control_plane_common",
            roots[0] / "control_plane_common.py",
            installed_payloads["control_plane_common.py"],
        )
    finally:
        for name, prior in previous.items():
            if prior is absent:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = prior
    return common.validate_storage_policy, policy, version


def _checked_in_control_plane_common() -> types.ModuleType:
    root = Path(preparation.REPO_ROOT) / "tst/publication/frontier_control_plane"
    payloads = {
        name: (root / name).read_bytes()
        for name in (
            "control_plane_common.py",
            "operator_attestation.py",
            "q011_pressure_review_packet_verifier.py",
        )
    }
    operator = _test_private_module(
        "_q009_test_checked_in_operator",
        root / "operator_attestation.py",
        payloads["operator_attestation.py"],
    )
    pressure = _test_private_module(
        "_q009_test_checked_in_pressure",
        root / "q011_pressure_review_packet_verifier.py",
        payloads["q011_pressure_review_packet_verifier.py"],
    )
    aliases = {
        "operator_attestation": operator,
        "q011_pressure_review_packet_verifier": pressure,
    }
    absent = object()
    previous = {name: sys.modules.get(name, absent) for name in aliases}
    sys.modules.update(aliases)
    try:
        return _test_private_module(
            "_q009_test_checked_in_control_plane_common",
            root / "control_plane_common.py",
            payloads["control_plane_common.py"],
        )
    finally:
        for name, prior in previous.items():
            if prior is absent:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = prior


def _assert_authority_false(value: object) -> None:
    if isinstance(value, dict):
        for key, member in value.items():
            if key.endswith("_authorized") or key.endswith("_qualified"):
                assert member is False, key
            _assert_authority_false(member)
    elif isinstance(value, list):
        for member in value:
            _assert_authority_false(member)


def test_exact_five_case_seven_stage_materialization_is_deterministic() -> None:
    manifest, launches, files = _materialization()
    assert manifest["logical_case_count"] == 5
    assert manifest["scheduler_stage_count"] == 7
    assert len(launches) == 7
    assert len(files) == 25
    assert sum(
        item["launch_contract_syntax_validated_by_minimal_projection"]
        for item in launches
    ) == 5
    second_manifest, second_files = preparation.build_materialization()
    assert manifest == second_manifest
    assert files == second_files
    preparation.validate_materialization(manifest, files)


def test_exact_decks_disable_physical_outputs_and_bind_case_stage_inputs() -> None:
    manifest, launches, files = _materialization()
    for launch, record in zip(launches, manifest["stage_records"]):
        blocks = _deck_blocks(files[record["deck"]["path"]].decode())
        stage = launch["identity"]["scheduler_stage"]
        outputs = {key: value for key, value in blocks.items() if key.startswith("output")}
        if stage == "checkpoint":
            assert outputs == {
                "output1": {
                    "file_type": "rst",
                    "id": "q009_checkpoint",
                    "dt": "44.0",
                    "single_file_per_rank": "false",
                }
            }
        else:
            assert outputs == {}
        case = launch["logical_case_contract"]
        assert int(blocks["particles"]["pic_random_seed"]) == case["random_seed"]
        assert float(blocks["particles"]["pic_load_balance_cost_per_particle"]) == (
            case["load_balance"]["pic_load_balance_cost_per_particle"]
        )


def test_review_candidates_are_syntactic_only_and_all_authority_is_false() -> None:
    manifest, launches, files = _materialization()
    _assert_authority_false(manifest)
    _assert_authority_false(
        [json.loads(payload) for path, payload in files.items() if path.endswith(".json")]
    )
    assert all(value is False for value in preparation.AUTHORIZATION_BOUNDARY.values())
    assert all(value is False for value in preparation.ABSENT_CAPABILITIES.values())
    assert all(value is False for value in manifest["execution_boundary"].values())
    assert manifest["authority_semantics"] == preparation.AUTHORITY_SEMANTICS
    assert manifest["threat_model"] == preparation.THREAT_MODEL
    assert manifest["threat_model"]["privileged_same_process_python_trusted"] is True
    assert manifest["threat_model"]["privileged_same_process_python_resisted"] is False
    assert manifest["threat_model"]["all_public_outputs_can_grant_authority"] is False
    for launch in launches:
        assert launch["authorization"] == preparation.AUTHORIZATION_BOUNDARY
        assert launch["authority_semantics"] == preparation.AUTHORITY_SEMANTICS
        if launch["launch_contract_syntax_validated_by_minimal_projection"]:
            assert preparation.validate_launch_contract(
                launch["launch_contract_candidate"]
            ) == launch["launch_contract_candidate"]


def test_minimal_syntax_projection_matches_captured_reference_for_generated_contracts() -> None:
    reference = _checked_in_control_plane_common()
    _, launches, _ = _materialization()
    for launch in launches:
        candidate = launch["launch_contract_candidate"]
        if candidate is None:
            continue
        assert preparation.validate_launch_contract(candidate) == (
            reference.validate_launch_contract(candidate)
        )
        assert preparation.launch_contract_sha256(candidate) == (
            reference.launch_contract_sha256(candidate)
        )


def test_no_receipt_or_admission_implementation_or_generated_member_exists() -> None:
    manifest, _, files = _materialization()
    callables = {
        name
        for name, value in vars(preparation).items()
        if inspect.isfunction(value) and value.__module__ == preparation.__name__
    }
    forbidden_callables = {
        "produce_receipt",
        "validate_receipt",
        "produce_from_run_root",
        "admit_receipt",
        "admit_run",
    }
    assert callables.isdisjoint(forbidden_callables)
    assert not any("receipt" in path or "admission" in path for path in files)
    assert not any(
        "receipt" in key or "admission" in key for key in manifest["source_bindings"]
    )
    assert manifest["minimal_structured_syntax_compatibility"][
        "receipt_producer_validator_present"
    ] is False
    assert manifest["minimal_structured_syntax_compatibility"][
        "admission_bridge_present"
    ] is False
    source = Path(preparation.__file__).read_text(encoding="utf-8")
    assert "q009_dynamic_amr_aggregate_receipt_schema_v1.json" not in source
    assert "q009_dynamic_amr_aggregate_receipt_successor_v1.py" not in source


def test_every_required_execution_blocker_is_explicit_and_unresolved() -> None:
    manifest, launches, _ = _materialization()
    common = set(preparation.COMMON_RUNTIME_BLOCKERS)
    for launch in launches:
        assert common <= set(launch["required_blockers"])
    migration = [
        launch
        for launch in launches
        if preparation.MIGRATION_DESIGN_BLOCKER in launch["required_blockers"]
    ]
    continuations = [
        launch
        for launch in launches
        if preparation.RESTART_CONTINUATION_BLOCKER in launch["required_blockers"]
    ]
    assert len(migration) == 2
    assert len(continuations) == 2
    compatibility = manifest["minimal_structured_syntax_compatibility"]
    assert compatibility["contract_complete_event_per_rank_per_level_telemetry_present"] is False
    assert compatibility["pairwise_lb_on_off_acceptance_present"] is False
    assert compatibility["trusted_unique_rst_restart_launching_supported"] is False
    assert compatibility[
        "hardened_external_completion_root_inventory_receipts_present"
    ] is False


def test_restart_continuations_remain_launch_prohibited_handoffs() -> None:
    manifest, launches, files = _materialization()
    blocked = [
        launch
        for launch in launches
        if not launch["launch_contract_syntax_validated_by_minimal_projection"]
    ]
    assert len(blocked) == 2
    assert all(launch["identity"]["scheduler_stage"] == "continuation" for launch in blocked)
    assert all(launch["launch_contract_candidate"] is None for launch in blocked)
    handoffs = [
        _json(files[record["restart_continuation_handoff"]["path"]])
        for record in manifest["stage_records"]
        if "restart_continuation_handoff" in record
    ]
    assert len(handoffs) == 2
    assert all(handoff["required_blocker"] == preparation.RESTART_CONTINUATION_BLOCKER for handoff in handoffs)
    assert all(handoff["authorization"] == preparation.AUTHORIZATION_BOUNDARY for handoff in handoffs)
    candidate = copy.deepcopy(launches[2]["launch_contract_candidate"])
    candidate["actions"][0]["arguments"] = [{"literal": "-r"}, {"literal": "/tmp/restart"}]
    with pytest.raises(ValueError, match="restart input is not authorized"):
        preparation.validate_launch_contract(candidate)


def test_budget_serial_order_and_policy_review_fragments_remain_fail_closed() -> None:
    manifest, launches, files = _materialization()
    budget = _json(files["batch_budget_accounting_input.json"])
    assert budget["contract_hard_cap_node_hours"] == 100.0
    assert budget["contract_hard_cap_artifact_storage_gib"] == 200.0
    assert budget["maximum_live_q009_submissions"] == 1
    assert budget["maximum_logical_attempts_per_case"] == 1
    assert budget["maximum_retries_per_case"] == 0
    assert budget["empty_user_queue_required_before_every_serial_stage"]
    assert launches[-1]["prior_stage_closures_required"] == [
        launch["identity"]["stage_id"] for launch in launches[:-1]
    ]
    policy = _json(files["registered_policy_slice_review_fragment.json"])
    assert policy["status"] == "review_fragment_only_not_live_policy"
    assert policy["live_policy_mutation_authorized"] is False
    assert policy["authorization"] == preparation.AUTHORIZATION_BOUNDARY


def test_fabricated_or_nonexistent_operational_bindings_have_no_input_surface() -> None:
    fabricated = {
        "source_commit": "b" * 40,
        "executable_path": "/does/not/exist/athena",
        "executable_sha256": "f" * 64,
    }
    assert not hasattr(preparation, "validate_final_bindings")
    assert not hasattr(preparation, "validate_exact_launch_review_candidate")
    with pytest.raises(TypeError):
        preparation.build_materialization(fabricated)
    with pytest.raises(TypeError):
        preparation.validate_materialization({}, {}, fabricated)
    with pytest.raises(TypeError):
        preparation.build_materialization(final_bindings=fabricated)
    with pytest.raises(TypeError):
        preparation.validate_materialization({}, {}, final_bindings=fabricated)
    assert not hasattr(preparation, "materialize_review_bundle")
    manifest, launches, files = _materialization()
    serialized = json.dumps(manifest) + b"".join(files.values()).decode()
    assert "selected_final_bindings" not in serialized
    assert "PENDING_FINAL_" not in serialized
    assert manifest["operational_binding_input_supported"] is False
    assert manifest["exact_launch_candidate_validation_supported"] is False
    assert all(launch["operational_bindings_accepted"] is False for launch in launches)
    for launch in launches:
        assert launch["unresolved_operational_bindings"] == list(
            preparation.UNRESOLVED_OPERATIONAL_BINDINGS
        )


def test_persistent_review_bundle_publication_and_write_helpers_are_absent() -> None:
    forbidden_callables = {
        "_canonical_review_root_descriptor",
        "_exact_frozen_tree_snapshot_at",
        "_freeze_tree_at",
        "_open_relative_directory",
        "_remove_created_tree_at",
        "_require_materialized_root_identity",
        "_scan_exact_frozen_directory_at",
        "_validated_review_output",
        "_verify_exact_frozen_tree_at",
        "_write_new_at",
        "main",
        "materialize_review_bundle",
    }
    assert _module_callables() == {
        "_is_exact_json_graph",
        "_require",
        "_require_trusted_runtime",
        "build_materialization",
        "launch_contract_sha256",
        "validate_launch_contract",
        "validate_materialization",
    }
    assert _module_callables().isdisjoint(forbidden_callables)
    assert not _filesystem_mutator_calls()
    read_calls = [
        node
        for node in ast.walk(_source_tree())
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "reader"
    ]
    assert len(read_calls) == 1
    assert isinstance(read_calls[0].args[1], ast.Constant)
    assert read_calls[0].args[1].value == "rb"
    assert not any(
        name in vars(preparation)
        for name in (
            "_BUILTIN_FUNCTION_TYPE",
            "_FUNCTION_TYPE",
            "_MAPPING_PROXY_TYPE",
            "_PY_TPFLAGS_HEAPTYPE",
            "_TYPE_TYPE",
            "_TRUSTED_COMPILE",
            "_TRUSTED_EXEC",
            "_TRUSTED_OPEN",
            "_JSON_ESCAPES",
            "_canonical_json",
            "_decode_json",
            "_is_exact_file_graph",
            "_json_bytes",
            "_parse_deck",
            "_sha256_bytes",
            "_stable_repository_bytes",
            "_strict_equal",
            "_minimal_launch_syntax_surfaces",
        )
    )
    source = Path(preparation.__file__).read_text(encoding="utf-8")
    assert "_SAFE_REVIEW_BUNDLE_ID" not in source
    assert "_validated_review_output" not in source
    manifest, _ = preparation.build_materialization()
    compatibility = manifest["minimal_structured_syntax_compatibility"]
    assert compatibility["persistent_review_bundle_publication_api_present"] is False
    assert compatibility["filesystem_write_or_destructive_cleanup_path_present"] is False
    assert compatibility["import_time_preparation_graph_cached"] is True
    assert compatibility["import_time_preparation_graph_authoritative"] is False
    assert compatibility["ordinary_module_global_drift_changes_cached_preparation"] is False
    assert compatibility[
        "minimal_validator_builtin_mapping_recoverable_and_mutable_by_privileged_same_process_python"
    ] is True
    assert compatibility["in_process_integrity_guards_are_authority_boundary"] is False
    assert compatibility["validator_process_isolation_present"] is False
    assert compatibility["authoritative_output_or_receipt_produced"] is False
    assert compatibility["all_public_outputs_non_authoritative_review_inputs"] is True
    assert compatibility["authority_requires_separate_external_verifier_and_admission"] is True
    assert compatibility["retained_mutable_trusted_primitive_globals_present"] is False
    assert compatibility["post_freeze_construction_or_repository_read_helpers_exposed"] is False
    assert compatibility["post_freeze_mutable_primitive_verification_roots_exposed"] is False
    assert compatibility["retained_module_local_callable_count"] == 7
    assert compatibility[
        "retained_runtime_guard_uses_identity_pinned_non_write_non_exec_primitives"
    ] is True
    assert compatibility["privileged_same_process_can_synchronize_guard_roots"] is True
    assert compatibility[
        "spoofed_heap_type_static_builtin_wrappers_fail_before_execution"
    ] is True
    assert compatibility["reachable_internal_guard_code_identity_checks_present"] is True
    assert compatibility["cyclic_or_aliased_container_graphs_rejected_as_non_json"] is True
    assert compatibility["hostile_loader_file_rejected_before_protocol_calls"] is True
    assert preparation.ABSENT_CAPABILITIES[
        "persistent_review_bundle_publication_implemented"
    ] is False


def test_control_plane_module_poisoning_cannot_change_minimal_syntax_projection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    poison = types.ModuleType("control_plane_common")
    poison.TRUSTED_LAUNCH_EXECUTOR = "poisoned"
    poison.validate_launch_contract = lambda value: {"poisoned": value}
    poison.launch_contract_sha256 = lambda value: "0" * 64
    operator_poison = types.ModuleType("operator_attestation")
    pressure_poison = types.ModuleType("q011_pressure_review_packet_verifier")
    monkeypatch.setitem(sys.modules, "control_plane_common", poison)
    monkeypatch.setitem(sys.modules, "operator_attestation", operator_poison)
    monkeypatch.setitem(
        sys.modules, "q011_pressure_review_packet_verifier", pressure_poison
    )
    validator = preparation.validate_launch_contract
    digest = preparation.launch_contract_sha256
    bindings = preparation.CAPTURED_CONTROL_PLANE_REFERENCE_BINDINGS
    assert preparation.TRUSTED_LAUNCH_EXECUTOR == "trusted_trampoline_athena_argv_v1"
    assert validator is not poison.validate_launch_contract
    assert digest is not poison.launch_contract_sha256
    assert _minimal_callable(validator, "validate_launch_contract").__module__ == (
        preparation.MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID
    )
    assert _minimal_callable(digest, "launch_contract_sha256").__module__ == (
        preparation.MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID
    )
    assert sys.modules["operator_attestation"] is operator_poison
    assert sys.modules["q011_pressure_review_packet_verifier"] is pressure_poison
    for binding in bindings.values():
        payload = (Path(preparation.REPO_ROOT) / binding["path"]).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == binding["sha256"]
        assert len(payload) == binding["byte_count"]
        assert binding["loading_role"] == "captured_verified_reference_not_executed"


def test_transitive_sys_modules_poisoning_cannot_change_minimal_syntax_semantics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, launches, _ = _materialization()
    selected = next(
        launch
        for launch in launches
        if launch["launch_contract_syntax_validated_by_minimal_projection"]
    )
    valid = selected["launch_contract_candidate"]
    expected_digest = selected["launch_contract_sha256"]
    poison = types.ModuleType("poison")
    poison.fullmatch = lambda *args, **kwargs: object()
    poison.dumps = lambda *args, **kwargs: "{}"
    poison.sha256 = lambda *args, **kwargs: types.SimpleNamespace(
        hexdigest=lambda: "0" * 64
    )
    for name in ("re", "json", "hashlib"):
        monkeypatch.setitem(sys.modules, name, poison)
    validator = preparation.validate_launch_contract
    digest = preparation.launch_contract_sha256
    malformed = copy.deepcopy(valid)
    malformed["actions"][0]["action_id"] = "BAD SPACE"
    with pytest.raises(ValueError, match="action ID is malformed"):
        validator(malformed)
    assert validator(valid) == valid
    assert digest(valid) == expected_digest
    assert digest(valid) != "0" * 64


def test_preimport_hashlib_poisoning_cannot_change_minimal_syntax_semantics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, launches, _ = _materialization()
    expected = next(
        launch
        for launch in launches
        if launch["launch_contract_syntax_validated_by_minimal_projection"]
    )
    poison = types.ModuleType("hashlib")
    poison.sha256 = lambda *args, **kwargs: types.SimpleNamespace(
        hexdigest=lambda: "0" * 64
    )
    monkeypatch.setitem(sys.modules, "hashlib", poison)
    payload = Path(preparation.__file__).read_bytes()
    fresh = _test_private_module(
        "_q009_test_preimport_hashlib_poison",
        Path(preparation.__file__),
        payload,
    )
    validator = fresh.validate_launch_contract
    digest = fresh.launch_contract_sha256
    candidate = expected["launch_contract_candidate"]
    assert "hashlib" not in fresh.__dict__
    assert not hasattr(fresh, "_sha256_bytes")
    assert not hasattr(fresh, "_stable_repository_bytes")
    assert validator(candidate) == candidate
    assert digest(candidate) == expected["launch_contract_sha256"]
    assert digest(candidate) != "0" * 64
    manifest, files = fresh.build_materialization()
    fresh.validate_materialization(manifest, files)


def test_fresh_import_poisoned_standard_modules_have_no_authority_or_side_effects() -> None:
    probe = _fresh_process(
        r'''
import sys
import types

calls = []

class Poison(types.ModuleType):
    def __getattr__(self, name):
        calls.append((self.__name__, name))
        raise AssertionError("poisoned module attribute was accessed")

for name in (
    "argparse",
    "hashlib",
    "importlib",
    "io",
    "json",
    "math",
    "os",
    "pathlib",
    "re",
    "stat",
    "types",
    "typing",
):
    sys.modules[name] = Poison(name)

from tst.publication import q009_deterministic_execution_preparation_only_successor_v2 as q

manifest, files = q.build_materialization()
q.validate_materialization(manifest, files)
assert not calls
assert all(value is False for value in manifest["execution_boundary"].values())
assert not hasattr(q, "_sha256_bytes")
assert not hasattr(q, "_stable_repository_bytes")
print("fresh-standard-poison-probe-passed")
'''
    )
    assert probe.returncode == 0, probe.stderr
    assert probe.stdout.strip() == "fresh-standard-poison-probe-passed"


def test_hostile_loader_file_is_rejected_before_protocol_calls(tmp_path: Path) -> None:
    sentinel = tmp_path / "hostile-loader-file-side-effect"
    calls: list[str] = []

    class HostileFile:
        def rsplit(self, *args: object, **kwargs: object) -> object:
            calls.append("rsplit")
            sentinel.write_text("called", encoding="utf-8")
            return []

        def __str__(self) -> str:
            calls.append("__str__")
            sentinel.write_text("called", encoding="utf-8")
            return "forged"

    path = Path(preparation.__file__)
    module = types.ModuleType("_q009_hostile_loader_file")
    module.__file__ = HostileFile()
    with pytest.raises(BaseException):
        exec(
            compile(path.read_bytes(), str(path), "exec", dont_inherit=True),
            module.__dict__,
        )
    assert not calls
    assert not sentinel.exists()


def test_fresh_import_hostile_builtin_all_fails_closed_without_calling_it() -> None:
    probe = _fresh_process(
        r'''
import builtins

calls = []

def hostile_all(*args, **kwargs):
    calls.append((args, kwargs))
    return True

builtins.all = hostile_all
try:
    from tst.publication import q009_deterministic_execution_preparation_only_successor_v2
except BaseException:
    assert not calls
    print("fresh-builtin-poison-probe-failed-closed")
else:
    raise AssertionError("poisoned builtins.all was accepted")
'''
    )
    assert probe.returncode == 0, probe.stderr
    assert probe.stdout.strip() == "fresh-builtin-poison-probe-failed-closed"


def test_post_import_hostile_builtin_all_fails_closed_without_calling_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, launches, files = _materialization()
    selected = next(
        launch
        for launch in launches
        if launch["launch_contract_syntax_validated_by_minimal_projection"]
    )
    valid = selected["launch_contract_candidate"]
    malformed = copy.deepcopy(valid)
    malformed["actions"][0]["action_id"] = "BAD SPACE"
    calls: list[object] = []

    def hostile_all(value: object) -> bool:
        calls.append(value)
        return True

    with monkeypatch.context() as poisoned:
        poisoned.setattr(builtins, "all", hostile_all)
        for operation in (
            lambda: preparation.validate_launch_contract(valid),
            lambda: preparation.launch_contract_sha256(valid),
            lambda: preparation.validate_launch_contract(malformed),
            preparation.build_materialization,
            lambda: preparation.validate_materialization(manifest, files),
            lambda: preparation._is_exact_json_graph(valid),
            lambda: preparation._require(True, "post-import builtin poison probe"),
        ):
            try:
                operation()
            except preparation.PreparationError:
                pass
            else:
                raise AssertionError("operation accepted poisoned runtime primitives")
    assert not calls
    assert not hasattr(preparation, "_stable_repository_bytes")


def test_ordinary_module_global_drift_does_not_change_cached_preparation() -> None:
    manifest, files = preparation.build_materialization()
    with pytest.raises(TypeError):
        preparation.AUTHORIZATION_BOUNDARY["publication_authorized"] = True
    forged = copy.deepcopy(manifest)
    forged["execution_boundary"]["publication_authorized"] = True
    original_validator = preparation.validate_launch_contract
    original_digest = preparation.launch_contract_sha256
    original_boundary = preparation.AUTHORIZATION_BOUNDARY
    try:
        preparation.validate_launch_contract = lambda value: value
        preparation.launch_contract_sha256 = lambda value: "0" * 64
        preparation.AUTHORIZATION_BOUNDARY = {
            key: key == "publication_authorized" for key in original_boundary
        }
        rebuilt, rebuilt_files = preparation.build_materialization()
        assert rebuilt == manifest
        assert rebuilt_files == files
        assert rebuilt["execution_boundary"]["publication_authorized"] is False
        with pytest.raises(preparation.PreparationError, match="manifest drifted"):
            preparation.validate_materialization(forged, files)
    finally:
        preparation.validate_launch_contract = original_validator
        preparation.launch_contract_sha256 = original_digest
        preparation.AUTHORIZATION_BOUNDARY = original_boundary


def test_privileged_same_process_mutation_is_explicitly_outside_contract() -> None:
    probe = _fresh_process(
        r'''
import copy
import gc
import json

from tst.publication import q009_deterministic_execution_preparation_only_successor_v2 as q

assert q.THREAT_MODEL["privileged_same_process_python_trusted"] is True
assert q.THREAT_MODEL["privileged_same_process_python_resisted"] is False
assert q.THREAT_MODEL["closure_or_code_object_mutation_resisted"] is False
assert q.THREAT_MODEL["mapping_proxy_backing_dict_recovery_resisted"] is False
assert q.THREAT_MODEL["all_public_outputs_can_grant_authority"] is False

manifest, files = q.build_materialization()
launch_path = manifest["stage_records"][0]["launch_candidate"]["path"]
malformed = json.loads(files[launch_path])
malformed["launch_contract_candidate"]["actions"][0]["action_id"] = "bad space"

wrapper_cells = dict(
    zip(
        q.validate_launch_contract.__code__.co_freevars,
        q.validate_launch_contract.__closure__,
    )
)
minimal_validator = wrapper_cells["validator"].cell_contents
builtin_backing = next(
    value
    for value in gc.get_referents(minimal_validator.__builtins__)
    if type(value) is dict
)
builtin_backing["all"] = lambda value: True
assert q.validate_launch_contract(malformed["launch_contract_candidate"])

forged = copy.deepcopy(manifest)
forged["execution_boundary"]["publication_authorized"] = True
forged["authority_semantics"] = "forged_same_process_value"
forged_payload = (
    json.dumps(forged, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    + "\n"
).encode("utf-8")
for public in (q.build_materialization, q.validate_materialization):
    cells = dict(zip(public.__code__.co_freevars, public.__closure__))
    cells["manifest_bytes"].cell_contents = forged_payload
rebuilt, rebuilt_files = q.build_materialization()
assert rebuilt["execution_boundary"]["publication_authorized"] is True
q.validate_materialization(rebuilt, rebuilt_files)

for forbidden in (
    "admit_receipt",
    "admit_run",
    "main",
    "materialize_review_bundle",
    "produce_receipt",
    "promote_active_policy",
    "submit",
    "validate_receipt",
):
    assert not hasattr(q, forbidden)
print("privileged-same-process-mutation-outside-non-authoritative-contract")
'''
    )
    assert probe.returncode == 0, probe.stderr
    assert probe.stdout.strip() == (
        "privileged-same-process-mutation-outside-non-authoritative-contract"
    )


def test_minimal_validator_guards_detect_ordinary_unsynchronized_drift() -> None:
    _, launches, _ = _materialization()
    valid = next(
        launch["launch_contract_candidate"]
        for launch in launches
        if launch["launch_contract_syntax_validated_by_minimal_projection"]
    )
    malformed = copy.deepcopy(valid)
    malformed["actions"][0]["action_id"] = "BAD SPACE"
    internal = _minimal_callable(
        preparation.validate_launch_contract, "validate_launch_contract"
    )
    captured_builtins = internal.__builtins__
    calls: list[object] = []

    def hostile_all(value: object) -> bool:
        calls.append(value)
        return True

    with pytest.raises(TypeError):
        captured_builtins["all"] = hostile_all
    with pytest.raises(ValueError, match="action ID is malformed"):
        preparation.validate_launch_contract(malformed)
    assert not calls

    namespace = internal.__globals__
    original = namespace["__builtins__"]
    try:
        namespace["__builtins__"] = {"all": hostile_all}
        with pytest.raises(
            preparation.PreparationError,
            match="validator namespace drifted",
        ):
            preparation.validate_launch_contract(valid)
    finally:
        namespace["__builtins__"] = original
    assert not calls

    original_helper = namespace["_lower_identifier"]
    try:
        namespace["_lower_identifier"] = lambda value: True
        with pytest.raises(preparation.PreparationError, match="validator code drifted"):
            preparation.validate_launch_contract(valid)
    finally:
        namespace["_lower_identifier"] = original_helper

    original_code = internal.__code__
    try:
        internal.__code__ = (lambda value: value).__code__
        with pytest.raises(preparation.PreparationError, match="validator code drifted"):
            preparation.validate_launch_contract(valid)
    finally:
        internal.__code__ = original_code


def test_guard_code_and_comparator_checks_detect_ordinary_drift() -> None:
    manifest, launches, files = _materialization()
    valid = next(
        launch["launch_contract_candidate"]
        for launch in launches
        if launch["launch_contract_syntax_validated_by_minimal_projection"]
    )
    internal = _minimal_callable(
        preparation.validate_launch_contract, "validate_launch_contract"
    )
    namespace = internal.__globals__
    validator_guard = _closure_callable(
        preparation.validate_launch_contract, "require_exact_validator_state"
    )
    runtime_guard = _closure_callable(validator_guard, "require_exact_runtime")

    root_runtime_guard = preparation._require_trusted_runtime
    original_root_runtime_guard_code = root_runtime_guard.__code__
    try:
        root_runtime_guard.__code__ = _noop_code_like(root_runtime_guard)
        with pytest.raises(
            preparation.PreparationError,
            match="runtime guard drifted",
        ):
            preparation.validate_launch_contract(valid)
    finally:
        root_runtime_guard.__code__ = original_root_runtime_guard_code

    original_runtime_guard_code = runtime_guard.__code__
    try:
        runtime_guard.__code__ = _noop_code_like(runtime_guard)
        with pytest.raises(
            preparation.PreparationError,
            match="runtime guard drifted",
        ):
            preparation.validate_launch_contract(valid)
    finally:
        runtime_guard.__code__ = original_runtime_guard_code

    original_validator_guard_code = validator_guard.__code__
    original_helper = namespace["_lower_identifier"]
    try:
        validator_guard.__code__ = _noop_code_like(validator_guard)
        namespace["_lower_identifier"] = lambda value: True
        malformed = copy.deepcopy(valid)
        malformed["actions"][0]["action_id"] = "BAD SPACE"
        with pytest.raises(
            preparation.PreparationError,
            match="validator guard drifted",
        ):
            preparation.validate_launch_contract(malformed)
    finally:
        namespace["_lower_identifier"] = original_helper
        validator_guard.__code__ = original_validator_guard_code

    surface_guard = _closure_callable(
        preparation.build_materialization, "require_exact_surface_state"
    )
    original_surface_guard_code = surface_guard.__code__
    try:
        surface_guard.__code__ = _noop_code_like(surface_guard)
        with pytest.raises(
            preparation.PreparationError,
            match="frozen preparation guard drifted",
        ):
            preparation.build_materialization()
    finally:
        surface_guard.__code__ = original_surface_guard_code

    calls: list[object] = []
    original_exact_json_graph = preparation._is_exact_json_graph
    try:
        preparation._is_exact_json_graph = lambda value: calls.append(value) or True
        with pytest.raises(
            preparation.PreparationError,
            match="frozen preparation surface drifted",
        ):
            preparation.validate_materialization(manifest, files)
    finally:
        preparation._is_exact_json_graph = original_exact_json_graph
    assert not calls

    original_json_hex = preparation._JSON_HEX
    try:
        preparation._JSON_HEX = object()
        with pytest.raises(
            preparation.PreparationError,
            match="frozen preparation surface drifted",
        ):
            preparation.build_materialization()
    finally:
        preparation._JSON_HEX = original_json_hex


def test_hostile_non_exact_json_is_rejected_without_protocol_calls(tmp_path: Path) -> None:
    sentinel = tmp_path / "hostile-action-id-side-effect"
    calls: list[str] = []

    class Hostile:
        def __str__(self) -> str:
            calls.append("__str__")
            sentinel.write_text("called", encoding="utf-8")
            return "forged"

        def __iter__(self) -> object:
            calls.append("__iter__")
            return iter(())

        def __eq__(self, other: object) -> bool:
            calls.append("__eq__")
            return True

    class HostileDict(dict):
        def __iter__(self) -> object:
            calls.append("dict.__iter__")
            return super().__iter__()

        def items(self) -> object:
            calls.append("dict.items")
            return super().items()

    class HostileStr(str):
        def __str__(self) -> str:
            calls.append("str.__str__")
            return "forged"

        def __eq__(self, other: object) -> bool:
            calls.append("str.__eq__")
            return True

        __hash__ = str.__hash__

    _, launches, files = _materialization()
    candidate = copy.deepcopy(
        next(
            launch["launch_contract_candidate"]
            for launch in launches
            if launch["launch_contract_syntax_validated_by_minimal_projection"]
        )
    )
    candidate["actions"][0]["action_id"] = Hostile()
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_launch_contract(candidate)
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.launch_contract_sha256(candidate)
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_launch_contract(HostileDict(candidate))
    subclassed_id = copy.deepcopy(candidate)
    subclassed_id["actions"][0]["action_id"] = HostileStr("forged")
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_launch_contract(subclassed_id)

    manifest, _ = preparation.build_materialization()
    hostile_manifest = copy.deepcopy(manifest)
    hostile_manifest["status"] = Hostile()
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_materialization(hostile_manifest, files)
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_materialization(Hostile(), files)
    hostile_files = dict(files)
    hostile_files["hostile"] = Hostile()
    with pytest.raises(ValueError, match="exact byte graph"):
        preparation.validate_materialization(manifest, hostile_files)
    with pytest.raises(ValueError, match="exact byte graph"):
        preparation.validate_materialization(manifest, Hostile())
    cyclic: list[object] = []
    cyclic.append(cyclic)
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_launch_contract(cyclic)
    aliased: dict[str, object] = {}
    with pytest.raises(ValueError, match="exact JSON graph"):
        preparation.validate_materialization(
            {"left": aliased, "right": aliased},
            files,
        )
    assert not hasattr(preparation, "_json_bytes")
    assert not hasattr(preparation, "_strict_equal")
    assert not calls
    assert not sentinel.exists()


def test_fresh_import_hostile_build_class_fails_before_wrapper_execution() -> None:
    probe = _fresh_process(
        r'''
import builtins

calls = []
original = builtins.__build_class__

def hostile_build_class(*args, **kwargs):
    calls.append((args, kwargs))
    return original(*args, **kwargs)

builtins.__build_class__ = hostile_build_class
try:
    from tst.publication import q009_deterministic_execution_preparation_only_successor_v2
except BaseException:
    assert not calls
    print("fresh-build-class-poison-probe-failed-closed")
else:
    raise AssertionError("hostile __build_class__ was accepted")
'''
    )
    assert probe.returncode == 0, probe.stderr
    assert probe.stdout.strip() == "fresh-build-class-poison-probe-failed-closed"
    assert not any(isinstance(node, ast.ClassDef) for node in ast.walk(_source_tree()))


def test_fresh_import_spoofed_static_builtins_fail_before_execution() -> None:
    probe = _fresh_process(
        r'''
import builtins
import types

source_path = "tst/publication/q009_deterministic_execution_preparation_only_successor_v2.py"
with open(source_path, "rb") as source:
    code = compile(source.read(), source_path, "exec", dont_inherit=True)

calls = []
accepted = []
originals = {
    name: getattr(builtins, name)
    for name in (
        "OSError",
        "UnicodeDecodeError",
        "ValueError",
        "enumerate",
        "range",
        "zip",
    )
}

class Spoof:
    __module__ = "builtins"

    def __new__(cls, *args, **kwargs):
        calls.append((cls.__name__, args, kwargs))
        return originals[cls.__name__](*args, **kwargs)

for name, original in originals.items():
    Spoof.__name__ = name
    setattr(builtins, name, Spoof)
    module = types.ModuleType("_q009_spoofed_static_builtin_" + name)
    module.__file__ = source_path
    try:
        exec(code, module.__dict__)
    except BaseException:
        pass
    else:
        accepted.append(name)
    finally:
        setattr(builtins, name, original)

assert not calls
assert not accepted
print("fresh-static-builtin-spoof-probe-failed-closed")
'''
    )
    assert probe.returncode == 0, probe.stderr
    assert probe.stdout.strip() == "fresh-static-builtin-spoof-probe-failed-closed"


def test_post_import_unsynchronized_primitive_drift_fails_before_execution() -> None:
    probe = _fresh_process(
        r'''
import builtins

from tst.publication import q009_deterministic_execution_preparation_only_successor_v2 as q

removed_roots_and_helpers = (
    "_BUILTIN_FUNCTION_TYPE",
    "_FUNCTION_TYPE",
    "_MAPPING_PROXY_TYPE",
    "_PY_TPFLAGS_HEAPTYPE",
    "_TYPE_TYPE",
    "_json_bytes",
    "_parse_deck",
    "_sha256_bytes",
    "_stable_repository_bytes",
    "_strict_equal",
)
assert all(not hasattr(q, name) for name in removed_roots_and_helpers)

manifest, files = q.build_materialization()
launch_path = manifest["stage_records"][0]["launch_candidate"]["path"]
import json
candidate = json.loads(files[launch_path])
calls = []

def hostile_all(*args, **kwargs):
    calls.append("all")
    return True

def hostile_open(*args, **kwargs):
    calls.append("open")
    raise AssertionError("hostile open called")

builtins.all = hostile_all
builtins.open = hostile_open
q._TRUSTED_ALL = hostile_all
for operation in (
    q.build_materialization,
    lambda: q.validate_materialization(manifest, files),
    lambda: q.validate_launch_contract(candidate),
    lambda: q.launch_contract_sha256(candidate),
    lambda: q._is_exact_json_graph(candidate),
    lambda: q._require(True, "primitive drift probe"),
):
    try:
        operation()
    except BaseException:
        pass
    else:
        raise AssertionError("post-import primitive drift entered retained surface")
assert not calls
print("post-import-primitive-root-drift-probe-failed-closed")
'''
    )
    assert probe.returncode == 0, probe.stderr
    assert probe.stdout.strip() == "post-import-primitive-root-drift-probe-failed-closed"


def test_no_operational_control_plane_or_subprocess_surface_is_exposed() -> None:
    assert not hasattr(preparation, "_CONTROL_PLANE_MODULE")
    assert not hasattr(preparation, "_private_module")
    assert not hasattr(preparation, "subprocess")
    assert not any(
        hasattr(preparation, name)
        for name in (
            "TRUSTED_SBATCH",
            "admit_receipt",
            "admit_run",
            "produce_receipt",
            "promote_active_policy",
            "submit",
            "validate_receipt",
        )
    )
    tree = _source_tree()
    imports = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
    ]
    assert imports == []
    assert not any(
        isinstance(node, (ast.Import, ast.ImportFrom))
        and (
            getattr(node, "module", None) == "subprocess"
            or any(alias.name == "subprocess" for alias in getattr(node, "names", ()))
        )
        for node in ast.walk(tree)
    )
    assert _minimal_callable(
        preparation.validate_launch_contract, "validate_launch_contract"
    ).__module__ == (
        preparation.MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID
    )


def test_returned_minimal_syntax_callables_expose_no_operational_namespace() -> None:
    validator = _minimal_callable(
        preparation.validate_launch_contract, "validate_launch_contract"
    )
    digest = _minimal_callable(
        preparation.launch_contract_sha256, "launch_contract_sha256"
    )
    forbidden = {
        "atomic_write_bytes",
        "durable_replace_tree",
        "os",
        "promote_active_policy",
        "subprocess",
        "validate_storage_policy",
    }
    for function in (validator, digest):
        assert forbidden.isdisjoint(function.__globals__)
        assert not any(
            isinstance(value, types.ModuleType)
            for value in function.__globals__.values()
        )
        assert all(
            not isinstance(value, types.FunctionType)
            or value.__globals__ is function.__globals__
            for value in function.__globals__.values()
        )
        captured_builtins = function.__builtins__
        assert captured_builtins is function.__globals__["__builtins__"]
        assert "__import__" not in captured_builtins
        assert "open" not in captured_builtins


def test_rename_in_replacement_cleanup_attack_is_impossible_by_construction(
    tmp_path: Path,
) -> None:
    replacement = tmp_path / "rename-in-replacement"
    replacement.mkdir()
    sentinel = replacement / "must-survive.txt"
    sentinel.write_bytes(b"rename-in replacement remains untouched\n")
    identity = (replacement.stat().st_dev, replacement.stat().st_ino)
    manifest, files = preparation.build_materialization()
    preparation.validate_materialization(manifest, files)
    assert not _filesystem_mutator_calls()
    assert not any("cleanup" in name or "remove_created" in name for name in _module_callables())
    assert (replacement.stat().st_dev, replacement.stat().st_ino) == identity
    assert sentinel.read_bytes() == b"rename-in replacement remains untouched\n"


def test_post_final_verification_rewrite_attack_is_impossible_by_construction(
    tmp_path: Path,
) -> None:
    sentinel = tmp_path / "post-final-rewrite-target.json"
    sentinel.write_bytes(b'{"authority":false}\n')
    before = sentinel.read_bytes()
    manifest, files = preparation.build_materialization()
    preparation.validate_materialization(manifest, files)
    assert not hasattr(preparation, "materialize_review_bundle")
    assert not hasattr(preparation, "_verify_exact_frozen_tree_at")
    assert not _filesystem_mutator_calls()
    assert sentinel.read_bytes() == before


def test_live_storage_policy_validator_rejects_every_policy_review_candidate() -> None:
    validator, live_policy, version = _live_storage_policy_validator()
    assert validator(copy.deepcopy(live_policy), control_plane_version=version) == live_policy
    manifest, _, files = _materialization()
    for record in manifest["stage_records"]:
        candidate = _json(files[record["policy_candidate"]["path"]])[
            "policy_slice_candidate"
        ]
        policy = copy.deepcopy(live_policy)
        policy["registered_science_slices"] = [candidate]
        with pytest.raises(
            ValueError, match="registered-science authorization is malformed"
        ):
            validator(policy, control_plane_version=version)


def test_authority_drift_in_materialized_member_fails_exact_validation() -> None:
    manifest, _, files = _materialization()
    drifted = dict(files)
    path = manifest["stage_records"][0]["launch_candidate"]["path"]
    launch = _json(drifted[path])
    launch["authorization"]["publication_authorized"] = True
    drifted[path] = _json_bytes(launch)
    with pytest.raises(preparation.PreparationError, match="member drifted"):
        preparation.validate_materialization(manifest, drifted)


def test_readiness_binds_only_preparation_sources_and_reports_all_blockers() -> None:
    readiness = json.loads(READINESS_PATH.read_text(encoding="utf-8"))
    assert readiness["successor_id"] == preparation.SUCCESSOR_ID
    assert readiness["status"] == "deterministic_execution_preparation_complete_execution_blocked"
    assert readiness["exact_matrix"]["logical_case_count"] == 5
    assert readiness["exact_matrix"]["scheduler_stage_count"] == 7
    assert preparation.RESTART_CONTINUATION_BLOCKER in readiness["explicit_blockers"]
    assert preparation.MIGRATION_DESIGN_BLOCKER in readiness["explicit_blockers"]
    assert all(value is False for value in readiness["authorization"].values())
    assert all(value is False for value in readiness["absent_capabilities"].values())
    assert readiness["threat_model"] == dict(preparation.THREAT_MODEL)
    assert readiness["authority_semantics"] == preparation.AUTHORITY_SEMANTICS
    output_contract = readiness["materialized_output_contract"]
    assert output_contract["deterministic_in_memory_preparation_graph_supported"] is True
    assert output_contract["import_time_preparation_graph_cached"] is True
    assert output_contract["import_time_preparation_graph_authoritative"] is False
    assert output_contract["privileged_same_process_mutation_resisted"] is False
    assert output_contract["privileged_same_process_mutation_outside_contract"] is True
    assert output_contract["in_process_validation_can_produce_authority"] is False
    assert output_contract["external_isolated_verifier_present"] is False
    assert output_contract["persistent_filesystem_publication_supported"] is False
    assert output_contract["persistent_publication_api_present"] is False
    assert output_contract["persistence_helpers_present"] is False
    assert output_contract["filesystem_write_calls_present"] is False
    assert output_contract["retained_write_capable_file_primitive"] is False
    assert output_contract["retained_exec_or_compile_construction_surface"] is False
    assert output_contract["post_freeze_construction_or_repository_read_helpers_exposed"] is False
    assert output_contract["post_freeze_mutable_primitive_verification_roots_exposed"] is False
    assert output_contract["retained_module_local_callable_count"] == 7
    assert output_contract["post_final_stat_or_scan_rewrite_window_reachable"] is False
    assert output_contract["destructive_failure_cleanup_implemented"] is False
    assert output_contract[
        "name_based_unlink_or_rmdir_reachable_from_publication_or_cleanup"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "all_seven_policy_review_candidates_rejected_by_paired_installed_live_storage_policy_validator"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "captured_control_plane_reference_executed"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "transitive_sys_modules_imports_used_by_minimal_validator"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "preimport_hashlib_sys_modules_poisoning_affects_minimal_validator"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "preimport_standard_module_sys_modules_poisoning_affects_preparation"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "minimal_validator_uses_captured_builtin_primitives"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "retained_mutable_trusted_primitive_globals_present"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "post_freeze_construction_or_repository_read_helpers_exposed"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "post_freeze_mutable_primitive_verification_roots_exposed"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "retained_module_local_callable_count"
    ] == 7
    assert readiness["minimal_structured_syntax_compatibility"][
        "retained_runtime_guard_uses_identity_pinned_non_write_non_exec_primitives"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "minimal_validator_builtin_mapping_recoverable_and_mutable_by_privileged_same_process_python"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "minimal_validator_namespace_and_code_identity_checks_present"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "privileged_same_process_can_synchronize_guard_roots"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "spoofed_heap_type_static_builtin_wrappers_fail_before_execution"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "reachable_internal_guard_code_identity_checks_present"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "non_exact_json_rejected_before_coercion_or_protocol_calls"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "cyclic_or_aliased_container_graphs_rejected_as_non_json"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "hostile_loader_file_rejected_before_protocol_calls"
    ] is True
    assert readiness["minimal_structured_syntax_compatibility"][
        "fresh_import_hostile_build_class_wrapper_called"
    ] is False
    assert readiness["minimal_structured_syntax_compatibility"][
        "preverification_imports_or_class_definitions_present"
    ] is False
    assert not any(
        "receipt" in key or "admission" in key for key in readiness["source_bindings"]
    )
    for group in ("source_bindings", "foundation_bindings"):
        for binding in readiness[group].values():
            source = Path(preparation.REPO_ROOT) / binding["path"]
            assert hashlib.sha256(source.read_bytes()).hexdigest() == binding["sha256"]
    external = readiness["non_authorizing_external_rejection_probe_bindings"]
    assert external["paired_installed_control_plane_version"] == (
        json.loads(LIVE_ORION_POLICY_PATH.read_text(encoding="utf-8"))[
            "olcf_side_storage"
        ]["installed_control_plane_version"]
    )
    for name, binding in external.items():
        if name == "paired_installed_control_plane_version":
            continue
        source = Path(binding["path"])
        assert hashlib.sha256(source.read_bytes()).hexdigest() == binding["sha256"]
