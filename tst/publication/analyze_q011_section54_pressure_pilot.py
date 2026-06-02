#!/usr/bin/env python3
"""Strict offline analyzer for nonqualifying Q-011 pressure-sensitivity pilots.

The four short Frontier pilots are engineering calibration only.  This module
does not authorize execution, qualify Section 5.4 physics, or make a Sun-Bai
claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from . import q011_section54_pressure_pilot_execution as execution
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_pressure_pilot_execution as execution
    from pvtk_particles import ParticleVTKData, read_particle_vtk


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_NAME = "pressure_pilot_manifest.json"
PREREGISTRATION_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_pressure_pilot_preregistration_2026-06-01.json"
)
REGISTERED_EXECUTION_PREREGISTRATION_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json"
)
ACTIVE_DECK_PATH = (
    REPO_ROOT
    / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
AUTHORIZED_PUBLICATION_ROOT = (
    execution.AUTHORIZED_PIC_ROOT / "publication"
)
RESULT_RECORD_TYPE = "q011_section54_pressure_pilot_analysis"
EVIDENCE_CLASS = "engineering_calibration_only"
QUALIFICATION_EFFECT = "none_no_sun_bai_claim_no_execution_authorization"

_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_FNV1A64_PATTERN = re.compile(r"[0-9a-f]{16}")
_RANK_DIRECTORY_PATTERN = re.compile(r"rank_([0-9]{8})")
_PVTK_PROVENANCE_PATTERN = re.compile(
    rb"# AthenaK particle data at time= ([^ ]+)  nranks= ([0-9]+)  "
    rb"cycle=([0-9]+)  variables=prtcl_all"
)
_Q017_PATTERN = re.compile(
    r"q017\.telemetry\.([A-Za-z0-9_.]+)="
    r"([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)"
)
_NUMBER_PATTERN = re.compile(
    r"(?<![A-Za-z_])"
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
)
_NONFINITE_PATTERN = re.compile(r"(?i)(?:^|[^A-Za-z_])(?:nan|inf)(?:[^A-Za-z_]|$)")
_MHD_FIELDS = ("dens", "eint", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3")
_SCALAR_BIN_FIELDS = {
    "bmag": ("bmag",),
    "prtcl_jx": ("prtcl_jx",),
    "j2": ("j2",),
}
_PVTK_SCALARS = {
    "gid",
    "ptag",
    "species",
    "cr_source",
    "macro_weight",
    "birth_time",
    "deltaf_f0",
    "deltaf_weight",
}
_PVTK_SCALAR_LIST = [
    "gid",
    "ptag",
    "species",
    "cr_source",
    "macro_weight",
    "birth_time",
    "deltaf_f0",
    "deltaf_weight",
]
_Q017_REQUIRED = {
    "schema_version",
    "mpi.ranks",
    "cycles",
    "meshblocks.total",
    "mesh.active_cells",
    "particles.total",
    "throughput.zone_cycles_per_second",
    "throughput.particle_updates_per_second",
    "load.meshblock_efficiency",
    "load.particle_efficiency",
    "amr.enabled",
    "particle_memory.resident_records.bytes_total",
    "particle_memory.invalid_records",
    (
        "particle_memory.athenak_owned_tracked_kokkos_views."
        "allocated_snapshot_bytes_total"
    ),
    (
        "particle_memory.athenak_owned_tracked_kokkos_views."
        "allocated_high_water_bytes_rank_max"
    ),
}
_Q017_REQUIRED_LIST = [
    "schema_version",
    "mpi.ranks",
    "cycles",
    "meshblocks.total",
    "mesh.active_cells",
    "particles.total",
    "throughput.zone_cycles_per_second",
    "throughput.particle_updates_per_second",
    "load.meshblock_efficiency",
    "load.particle_efficiency",
    "amr.enabled",
    "particle_memory.resident_records.bytes_total",
    "particle_memory.invalid_records",
    (
        "particle_memory.athenak_owned_tracked_kokkos_views."
        "allocated_snapshot_bytes_total"
    ),
    (
        "particle_memory.athenak_owned_tracked_kokkos_views."
        "allocated_high_water_bytes_rank_max"
    ),
]
_SHOCK_PREFIXES = (
    "pic_parallel_shock feedback_diag_cfg:",
    "pic_parallel_shock feedback_diag:",
    "pic_parallel_shock source_transaction_diag:",
    "pic_parallel_shock removed_cr_sink:",
)
_CASES = (
    ("ps_p0_1p00", 1.0, "1.0"),
    ("ps_p0_0p05", 0.05, "0.05"),
    ("ps_p0_0p10", 0.1, "0.10"),
    ("ps_p0_0p20", 0.2, "0.20"),
)
_TIMES = (0.0, 15.0, 30.0, 45.0, 60.0)
_FIXED_OVERRIDES = (
    "mesh/nx1=100",
    "mesh/x1max=1200",
    "mesh/nx2=20",
    "mesh/x2max=240",
    "mesh_refinement/refinement=none",
    "mesh_refinement/num_levels=1",
    "time/tlim=60",
    "time/nlim=4096",
    "time/ndiag=50",
    "problem/ps_enable_curvature_amr=false",
    "problem/ps_feedback_diag_dcycle=50",
    "output1/variable=mhd_w_bcc",
    "output1/id=mhd_w_bcc",
    "output1/dt=15",
    "output2/dt=15",
    "output3/dt=15",
    "output4/dt=15",
    "output5/dt=15",
    "output6/dt=15",
)
_RUNTIME_PARAMETERS = {
    ("mesh", "nx1"): "100",
    ("mesh", "x1min"): "0.0",
    ("mesh", "x1max"): "1200",
    ("mesh", "nx2"): "20",
    ("mesh", "x2min"): "0.0",
    ("mesh", "x2max"): "240",
    ("mesh", "nx3"): "1",
    ("mesh", "x3min"): "0.0",
    ("mesh", "x3max"): "1.0",
    ("mesh_refinement", "refinement"): "none",
    ("mesh_refinement", "num_levels"): "1",
    ("time", "tlim"): "60",
    ("time", "nlim"): "4096",
    ("time", "ndiag"): "50",
    ("mhd", "gamma"): "1.66666666667",
    ("particles", "pic_physical_mode"): "paper_mhd_pic_vl2_tsc",
    ("problem", "ps_enable_curvature_amr"): "false",
    ("problem", "ps_feedback_diag_dcycle"): "50",
    ("output1", "variable"): "mhd_w_bcc",
    ("output1", "id"): "mhd_w_bcc",
    ("output1", "dt"): "15",
    ("output2", "variable"): "mhd_bmag",
    ("output2", "id"): "bmag",
    ("output2", "dt"): "15",
    ("output3", "variable"): "prtcl_jx",
    ("output3", "id"): "prtcl_jx",
    ("output3", "dt"): "15",
    ("output4", "variable"): "mhd_j2",
    ("output4", "id"): "j2",
    ("output4", "dt"): "15",
    ("output5", "variable"): "prtcl_all",
    ("output5", "id"): "prtcl_all",
    ("output5", "dt"): "15",
    ("output6", "file_type"): "rst",
    ("output6", "dt"): "15",
}


class PilotAnalysisError(ValueError):
    """Raised when the pressure-pilot bundle fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PilotAnalysisError(message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fnv1a64(payload: bytes) -> str:
    value = 14695981039346656037
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def _reject_constant(value: str) -> None:
    raise PilotAnalysisError(f"JSON constant is forbidden: {value}")


def _decode_json(payload: bytes, label: str) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for name, value in pairs:
            _require(name not in result, f"{label}: duplicate JSON key {name!r}")
            result[name] = value
        return result

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PilotAnalysisError(f"{label}: JSON is not UTF-8") from error
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=_reject_constant,
        )
    except json.JSONDecodeError as error:
        raise PilotAnalysisError(f"{label}: malformed JSON") from error


def _object(value: object, keys: set[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    mapping = value
    _require(set(mapping) == keys, f"{label}: keys drifted")
    return mapping


def _list(value: object, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected list")
    return value


def _text(value: object, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected nonempty text")
    return value


def _float(value: object, label: str) -> float:
    _require(type(value) is float, f"{label}: expected JSON float")
    _require(math.isfinite(value), f"{label}: expected finite float")
    return value


def _int(value: object, label: str) -> int:
    _require(type(value) is int, f"{label}: expected JSON integer")
    return value


def _sha256(value: object, label: str) -> str:
    text = _text(value, label)
    _require(_SHA256_PATTERN.fullmatch(text) is not None, f"{label}: invalid SHA-256")
    return text


def _relative_path(value: object, label: str) -> str:
    text = _text(value, label)
    candidate = PurePosixPath(text)
    _require(
        not candidate.is_absolute()
        and candidate.as_posix() == text
        and text != "."
        and all(part not in ("", ".", "..") for part in candidate.parts),
        f"{label}: unsafe root-relative path",
    )
    return text


def _binding(value: object, label: str) -> dict[str, str]:
    item = _object(value, {"path", "sha256"}, label)
    return {
        "path": _relative_path(item["path"], f"{label}/path"),
        "sha256": _sha256(item["sha256"], f"{label}/sha256"),
    }


def _regular_bytes(path: Path, label: str) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(path, flags)
    except OSError as error:
        raise PilotAnalysisError(f"{label}: unable to open retained regular file") from error
    try:
        before = os.fstat(fd)
        _require(stat.S_ISREG(before.st_mode), f"{label}: expected regular file")
        chunks = []
        while chunk := os.read(fd, 1024 * 1024):
            chunks.append(chunk)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    identity = lambda item: (
        item.st_dev,
        item.st_ino,
        item.st_mode,
        item.st_size,
        item.st_mtime_ns,
        item.st_ctime_ns,
    )
    _require(identity(before) == identity(after), f"{label}: file changed while reading")
    return b"".join(chunks)


def _source_binding(path: str, digest: str) -> dict[str, str]:
    return {"path": path, "sha256": digest}


def _strict_equal(actual: object, expected: object, label: str) -> None:
    _require(type(actual) is type(expected), f"{label}: primitive type drifted")
    if type(expected) is dict:
        _require(set(actual) == set(expected), f"{label}: object keys drifted")
        for key in sorted(expected):
            _strict_equal(actual[key], expected[key], f"{label}/{key}")
    elif type(expected) is list:
        _require(len(actual) == len(expected), f"{label}: list length drifted")
        for index, (left, right) in enumerate(zip(actual, expected)):
            _strict_equal(left, right, f"{label}[{index}]")
    else:
        _require(actual == expected, f"{label}: value drifted")


def _load_policy() -> dict[str, Any]:
    payload = _regular_bytes(PREREGISTRATION_PATH, "pressure-pilot preregistration")
    policy = _object(
        _decode_json(payload, "pressure-pilot preregistration"),
        {
            "record_type",
            "schema_version",
            "date",
            "gate",
            "classification",
            "qualification_effect",
            "scope",
            "active_deck_binding",
            "analysis_bindings",
            "pilot_contract",
            "execution_policy",
            "limitations",
        },
        "pressure-pilot preregistration",
    )
    _strict_equal(
        {key: policy[key] for key in (
            "record_type",
            "schema_version",
            "date",
            "gate",
            "classification",
            "qualification_effect",
        )},
        {
            "record_type": "q011_section54_pressure_pilot_preregistration",
            "schema_version": 1,
            "date": "2026-06-01",
            "gate": "Q-011",
            "classification": EVIDENCE_CLASS,
            "qualification_effect": QUALIFICATION_EFFECT,
        },
        "pressure-pilot preregistration identity",
    )
    _require(
        "engineering calibration only" in _text(policy["scope"], "policy/scope"),
        "policy/scope: engineering-only boundary drifted",
    )
    expected_deck = _source_binding(
        "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
        "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1",
    )
    expected_analysis = {
        "output_primitives": _source_binding(
            "tst/publication/analyze_q011_section54_outputs.py",
            "abe664ece580f70fcde022aa15583846c4050b342d2072c8cd9c181c0b9d6d22",
        ),
        "particle_vtk_reader": _source_binding(
            "tst/publication/pvtk_particles.py",
            "187254c4ed10ce20ec383e710dadfa2e03e9f45317cce53e5d644c1db2be7339",
        ),
    }
    _strict_equal(policy["active_deck_binding"], expected_deck, "policy/active_deck")
    _strict_equal(policy["analysis_bindings"], expected_analysis, "policy/analysis_bindings")
    contract = _object(
        policy["pilot_contract"],
        {
            "physical_mode",
            "cases",
            "fixed_overrides",
            "required_times_omega0_inverse",
            "required_snapshot_products",
            "mhd_w_bcc_required_variables",
            "pressure_formula",
            "particle_provenance_scalars",
            "particle_provenance_vectors",
            "terminal_particle_rule",
            "required_q017_telemetry",
            "required_shock_diagnostic_prefixes",
            "terminal_restart_policy",
            "overlay_profile_fields",
        },
        "policy/pilot_contract",
    )
    _strict_equal(contract["physical_mode"], "paper_mhd_pic_vl2_tsc", "policy/physical_mode")
    _strict_equal(
        contract["cases"],
        [
            {"case_id": name, "problem_ps_p0": value, "argv_value": argv}
            for name, value, argv in _CASES
        ],
        "policy/cases",
    )
    _strict_equal(contract["fixed_overrides"], list(_FIXED_OVERRIDES), "policy/overrides")
    _strict_equal(contract["required_times_omega0_inverse"], list(_TIMES), "policy/times")
    _strict_equal(
        contract["required_snapshot_products"],
        ["mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all"],
        "policy/products",
    )
    _strict_equal(contract["mhd_w_bcc_required_variables"], list(_MHD_FIELDS), "policy/mhd fields")
    _strict_equal(contract["particle_provenance_scalars"], _PVTK_SCALAR_LIST, "policy/pvtk scalars")
    _strict_equal(contract["particle_provenance_vectors"], ["vel"], "policy/pvtk vectors")
    _strict_equal(contract["required_q017_telemetry"], _Q017_REQUIRED_LIST, "policy/q017")
    _strict_equal(contract["required_shock_diagnostic_prefixes"], list(_SHOCK_PREFIXES), "policy/shock diagnostics")
    _strict_equal(contract["overlay_profile_fields"], ["rho", "p", "vx", "|B|"], "policy/profiles")
    _strict_equal(
        policy["execution_policy"],
        {
            "frontier_execution_authorized_by_this_record": False,
            "scheduler_commands_authorized_by_this_record": False,
            "scientific_evidence_eligible": False,
            "sun_bai_claim": False,
        },
        "policy/execution",
    )
    limitations = _list(policy["limitations"], "policy/limitations")
    _require(bool(limitations) and all(type(item) is str and item for item in limitations),
             "policy/limitations: expected nonempty text list")
    sources = {
        expected_deck["path"]: ACTIVE_DECK_PATH,
        expected_analysis["output_primitives"]["path"]: Path(output_primitives.__file__),
        expected_analysis["particle_vtk_reader"]["path"]: REPO_ROOT / expected_analysis["particle_vtk_reader"]["path"],
    }
    for relative, path in sources.items():
        binding = expected_deck if relative == expected_deck["path"] else (
            expected_analysis["output_primitives"]
            if relative == expected_analysis["output_primitives"]["path"]
            else expected_analysis["particle_vtk_reader"]
        )
        measured = _sha256_bytes(_regular_bytes(path, f"bound source {relative}"))
        _require(measured == binding["sha256"], f"bound source SHA-256 drifted: {relative}")
    return policy


def _manifest_schema(payload: bytes) -> dict[str, Any]:
    manifest = _object(
        _decode_json(payload, "pressure-pilot manifest"),
        {
            "schema_version",
            "record_type",
            "evidence_class",
            "qualification_effect",
            "active_deck_binding",
            "preregistration_binding",
            "registered_execution_preregistration_binding",
            "cases",
        },
        "pressure-pilot manifest",
    )
    _require(_int(manifest["schema_version"], "manifest/schema_version") == 1,
             "manifest/schema_version: expected 1")
    _strict_equal(manifest["record_type"], "q011_section54_pressure_pilot_bundle_manifest", "manifest/record_type")
    _strict_equal(manifest["evidence_class"], EVIDENCE_CLASS, "manifest/evidence_class")
    _strict_equal(manifest["qualification_effect"], QUALIFICATION_EFFECT, "manifest/qualification_effect")
    parsed_cases = []
    for index, raw in enumerate(_list(manifest["cases"], "manifest/cases")):
        label = f"manifest/cases[{index}]"
        case = _object(raw, {"case_id", "ps_p0", "overrides", "snapshots", "stdout", "terminal_restart"}, label)
        snapshots = []
        for snapshot_index, raw_snapshot in enumerate(_list(case["snapshots"], f"{label}/snapshots")):
            snap_label = f"{label}/snapshots[{snapshot_index}]"
            snapshot = _object(raw_snapshot, {"time", "mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all"}, snap_label)
            snapshots.append({
                "time": _float(snapshot["time"], f"{snap_label}/time"),
                **{
                    kind: _binding(snapshot[kind], f"{snap_label}/{kind}")
                    for kind in ("mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all")
                },
            })
        restart = _object(case["terminal_restart"], {"time", "manifest", "manifest_complete", "members"}, f"{label}/terminal_restart")
        members = []
        for member_index, raw_member in enumerate(_list(restart["members"], f"{label}/terminal_restart/members")):
            member_label = f"{label}/terminal_restart/members[{member_index}]"
            member = _object(raw_member, {"artifact", "complete"}, member_label)
            members.append({
                "artifact": _binding(member["artifact"], f"{member_label}/artifact"),
                "complete": _binding(member["complete"], f"{member_label}/complete"),
            })
        parsed_cases.append({
            "case_id": _text(case["case_id"], f"{label}/case_id"),
            "ps_p0": _float(case["ps_p0"], f"{label}/ps_p0"),
            "overrides": [_text(item, f"{label}/overrides") for item in _list(case["overrides"], f"{label}/overrides")],
            "snapshots": snapshots,
            "stdout": _binding(case["stdout"], f"{label}/stdout"),
            "terminal_restart": {
                "time": _float(restart["time"], f"{label}/terminal_restart/time"),
                "manifest": _binding(restart["manifest"], f"{label}/terminal_restart/manifest"),
                "manifest_complete": _binding(restart["manifest_complete"], f"{label}/terminal_restart/manifest_complete"),
                "members": members,
            },
        })
    return {
        **manifest,
        "active_deck_binding": _binding(manifest["active_deck_binding"], "manifest/active_deck_binding"),
        "preregistration_binding": _binding(manifest["preregistration_binding"], "manifest/preregistration_binding"),
        "registered_execution_preregistration_binding": _binding(
            manifest["registered_execution_preregistration_binding"],
            "manifest/registered_execution_preregistration_binding",
        ),
        "cases": parsed_cases,
    }


def _member_payload(root: Path, binding: Mapping[str, str], label: str) -> bytes:
    path = root / binding["path"]
    payload = _regular_bytes(path, label)
    _require(_sha256_bytes(payload) == binding["sha256"], f"{label}: SHA-256 drifted")
    _require(bool(payload), f"{label}: retained product is empty")
    return payload


def _expected_snapshot_paths(case_id: str, index: int) -> dict[str, str]:
    suffix = f"{index:05d}"
    prefix = f"cases/{case_id}"
    return {
        "mhd_w_bcc": f"{prefix}/bin/{case_id}.mhd_w_bcc.{suffix}.bin",
        "bmag": f"{prefix}/bin/{case_id}.bmag.{suffix}.bin",
        "prtcl_jx": f"{prefix}/bin/{case_id}.prtcl_jx.{suffix}.bin",
        "j2": f"{prefix}/bin/{case_id}.j2.{suffix}.bin",
        "prtcl_all": f"{prefix}/pvtk/{case_id}.prtcl_all.{suffix}.part.vtk",
    }


def _validate_case_identity(cases: Sequence[Mapping[str, Any]]) -> None:
    _require(len(cases) == len(_CASES), "pressure-pilot manifest must contain exactly four cases")
    for index, (case, expected) in enumerate(zip(cases, _CASES)):
        case_id, ps_p0, argv = expected
        _strict_equal(case["case_id"], case_id, f"case[{index}]/case_id")
        _strict_equal(case["ps_p0"], ps_p0, f"case[{index}]/ps_p0")
        _strict_equal(case["overrides"], [*_FIXED_OVERRIDES, f"problem/ps_p0={argv}"], f"case[{index}]/overrides")
        snapshots = case["snapshots"]
        _require(len(snapshots) == len(_TIMES), f"{case_id}: expected exactly five snapshots")
        for snap_index, (snapshot, time) in enumerate(zip(snapshots, _TIMES)):
            _strict_equal(snapshot["time"], time, f"{case_id}/snapshots[{snap_index}]/time")
            for kind, path in _expected_snapshot_paths(case_id, snap_index).items():
                _strict_equal(snapshot[kind]["path"], path, f"{case_id}/{kind} path")
        _strict_equal(case["stdout"]["path"], f"cases/{case_id}/stdout.txt", f"{case_id}/stdout path")
        restart = case["terminal_restart"]
        _strict_equal(restart["time"], 60.0, f"{case_id}/restart time")
        manifest_path = f"cases/{case_id}/rst/{case_id}.00004.rst.manifest"
        _strict_equal(restart["manifest"]["path"], manifest_path, f"{case_id}/restart manifest path")
        _strict_equal(restart["manifest_complete"]["path"], manifest_path + ".complete", f"{case_id}/restart manifest marker path")
        _require(bool(restart["members"]), f"{case_id}: restart publication has no members")


def _declared_paths(manifest: Mapping[str, Any]) -> set[str]:
    declared = {MANIFEST_NAME}
    for case in manifest["cases"]:
        declared.add(case["stdout"]["path"])
        for snapshot in case["snapshots"]:
            declared.update(snapshot[kind]["path"] for kind in ("mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all"))
        restart = case["terminal_restart"]
        declared.add(restart["manifest"]["path"])
        declared.add(restart["manifest_complete"]["path"])
        for member in restart["members"]:
            declared.add(member["artifact"]["path"])
            declared.add(member["complete"]["path"])
    expected_count = 1 + len(_CASES) * (1 + len(_TIMES) * 5)
    expected_count += sum(2 + 2 * len(case["terminal_restart"]["members"]) for case in manifest["cases"])
    _require(len(declared) == expected_count, "pressure-pilot manifest contains duplicate retained paths")
    return declared


def _actual_files(root: Path) -> set[str]:
    files = set()
    for directory, names, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        for name in names:
            _require(not (base / name).is_symlink(), "pressure-pilot bundle contains a directory symlink")
        for name in filenames:
            path = base / name
            _require(not path.is_symlink(), "pressure-pilot bundle contains a file symlink")
            _require(path.is_file(), "pressure-pilot bundle contains a special file")
            files.add(path.relative_to(root).as_posix())
    return files


def _validate_runtime_parameters(dataset: output_primitives.AthenaBinaryDataset, ps_p0: float) -> float:
    parameters = dataset.input_parameters
    for (block, key), expected in _RUNTIME_PARAMETERS.items():
        _require(block in parameters and key in parameters[block],
                 f"{dataset.source}: missing runtime parameter {block}/{key}")
        _strict_equal(parameters[block][key], expected, f"{dataset.source}: {block}/{key}")
    try:
        measured_ps_p0 = float(parameters.get("problem", {}).get("ps_p0", ""))
    except ValueError as error:
        raise PilotAnalysisError(f"{dataset.source}: problem/ps_p0 is malformed") from error
    _require(math.isfinite(measured_ps_p0) and measured_ps_p0 == ps_p0,
             f"{dataset.source}: problem/ps_p0 drifted")
    _require(dataset.root_grid_shape == (100, 20, 1), f"{dataset.source}: root-grid shape drifted")
    _require(dataset.domain_bounds == (0.0, 1200.0, 0.0, 240.0, 0.0, 1.0),
             f"{dataset.source}: domain bounds drifted")
    return float(parameters["mhd"]["gamma"])


def _binary_dataset(root: Path, binding: Mapping[str, str], time: float, ps_p0: float, fields: Sequence[str]) -> tuple[output_primitives.AthenaBinaryDataset, dict[str, output_primitives.CompositeGrid]]:
    payload = _member_payload(root, binding, f"binary product {binding['path']}")
    try:
        dataset = output_primitives.parse_athenak_binary_bytes(payload, source=binding["path"])
    except output_primitives.AnalysisError as error:
        raise PilotAnalysisError(f"malformed Athena bin {binding['path']}: {error}") from error
    _require(dataset.time == time, f"{binding['path']}: internal snapshot time drifted")
    _validate_runtime_parameters(dataset, ps_p0)
    _require(tuple(dataset.variable_names) == tuple(fields), f"{binding['path']}: variable inventory drifted")
    composites = {}
    for field in fields:
        try:
            composites[field] = output_primitives.compose_leaf_field(dataset, field)
        except output_primitives.AnalysisError as error:
            raise PilotAnalysisError(f"{binding['path']}: invalid composite field {field}: {error}") from error
    return dataset, composites


def _particle_data(payload: bytes, label: str) -> ParticleVTKData:
    with tempfile.NamedTemporaryFile(suffix=".part.vtk") as temporary:
        temporary.write(payload)
        temporary.flush()
        try:
            return read_particle_vtk(temporary.name)
        except (OSError, ValueError) as error:
            raise PilotAnalysisError(f"{label}: malformed particle VTK: {error}") from error


def _particle_report(root: Path, binding: Mapping[str, str], time: float) -> dict[str, int]:
    payload = _member_payload(root, binding, f"particle product {binding['path']}")
    lines = payload.splitlines()
    _require(len(lines) >= 2, f"{binding['path']}: particle VTK provenance header is absent")
    match = _PVTK_PROVENANCE_PATTERN.fullmatch(lines[1])
    _require(match is not None, f"{binding['path']}: particle VTK provenance header drifted")
    try:
        internal_time = float(match.group(1))
    except ValueError as error:
        raise PilotAnalysisError(f"{binding['path']}: particle VTK time is malformed") from error
    _require(math.isfinite(internal_time) and internal_time == time,
             f"{binding['path']}: particle VTK internal time drifted")
    _require(int(match.group(2)) >= 1 and int(match.group(3)) >= 0,
             f"{binding['path']}: particle VTK rank or cycle provenance is invalid")
    data = _particle_data(payload, binding["path"])
    _require(set(data.scalars) == _PVTK_SCALARS, f"{binding['path']}: particle scalar provenance drifted")
    _require(set(data.vectors) == {"vel"}, f"{binding['path']}: particle vector provenance drifted")
    count = data.points.shape[0]
    _require(data.points.shape == (count, 3) and data.vectors["vel"].shape == (count, 3),
             f"{binding['path']}: particle point/vector shape drifted")
    for name in _PVTK_SCALARS:
        _require(data.scalars[name].shape == (count,), f"{binding['path']}: particle scalar shape drifted")
    for values in (data.points, data.vectors["vel"], *(data.scalars[name] for name in ("macro_weight", "birth_time", "deltaf_f0", "deltaf_weight"))):
        _require(np.all(np.isfinite(values)), f"{binding['path']}: particle payload contains non-finite values")
    _require(np.all(data.scalars["gid"] >= 0) and np.all(data.scalars["ptag"] >= 0),
             f"{binding['path']}: particle gid or ptag is invalid")
    _require(np.unique(data.scalars["ptag"]).size == count, f"{binding['path']}: particle tags are not unique")
    _require(np.all(data.scalars["species"] == 0), f"{binding['path']}: particle species provenance drifted")
    _require(np.all(np.isin(data.scalars["cr_source"], (0, 1))),
             f"{binding['path']}: particle source provenance drifted")
    _require(np.all(data.scalars["macro_weight"] >= 0.0),
             f"{binding['path']}: particle macro weight is negative")
    retained_early = (data.scalars["cr_source"] == 1) & (data.scalars["birth_time"] < 45.0)
    if time == 60.0:
        _require(not np.any(retained_early),
                 f"{binding['path']}: retained post-removal shock-injected particle has birth_time < 45")
    return {
        "particle_count": int(count),
        "shock_injected_particle_count": int(np.count_nonzero(data.scalars["cr_source"] == 1)),
        "retained_early_shock_injected_particle_count": int(np.count_nonzero(retained_early)),
    }


def _stdout_report(root: Path, binding: Mapping[str, str], terminal_particles: int) -> dict[str, Any]:
    payload = _member_payload(root, binding, f"stdout product {binding['path']}")
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PilotAnalysisError(f"{binding['path']}: stdout is not UTF-8") from error
    telemetry = {}
    lines = text.splitlines()
    for line in lines:
        if not line.startswith("q017.telemetry."):
            continue
        match = _Q017_PATTERN.fullmatch(line)
        _require(match is not None, f"{binding['path']}: malformed Q017 telemetry line")
        name, raw = match.groups()
        _require(name not in telemetry, f"{binding['path']}: duplicate Q017 telemetry key {name}")
        value = float(raw)
        _require(math.isfinite(value) and value >= 0.0, f"{binding['path']}: invalid Q017 telemetry value")
        telemetry[name] = value
    _require(_Q017_REQUIRED <= set(telemetry), f"{binding['path']}: required Q017 telemetry is missing")
    _require(telemetry["schema_version"] == 2.0, f"{binding['path']}: Q017 schema drifted")
    _require(telemetry["mpi.ranks"] >= 1.0, f"{binding['path']}: Q017 rank count is invalid")
    _require(telemetry["cycles"] > 0.0, f"{binding['path']}: Q017 cycle count is invalid")
    _require(telemetry["mesh.active_cells"] == 2000.0, f"{binding['path']}: Q017 active-cell count drifted")
    _require(telemetry["particles.total"] == float(terminal_particles), f"{binding['path']}: Q017 particle count disagrees with t=60 pVTK")
    _require(telemetry["amr.enabled"] == 0.0, f"{binding['path']}: Q017 AMR state drifted")
    _require(telemetry["particle_memory.invalid_records"] == 0.0, f"{binding['path']}: Q017 reports invalid particle records")
    diagnostic_counts = {}
    for prefix in _SHOCK_PREFIXES:
        matching = [line for line in lines if line.startswith(prefix)]
        _require(bool(matching), f"{binding['path']}: required shock diagnostic is missing: {prefix}")
        for line in matching:
            _require(_NONFINITE_PATTERN.search(line) is None, f"{binding['path']}: shock diagnostic contains non-finite value")
            _require(bool(_NUMBER_PATTERN.findall(line)), f"{binding['path']}: shock diagnostic contains no numeric payload")
            _require(all(math.isfinite(float(value)) for value in _NUMBER_PATTERN.findall(line)),
                     f"{binding['path']}: shock diagnostic contains invalid numeric payload")
        diagnostic_counts[prefix] = len(matching)
    cfg = next(line for line in lines if line.startswith(_SHOCK_PREFIXES[0]))
    for fragment in ("mom_feedback=1", "eng_feedback=1", "pic_enable_2d3v=1", "use_delta_feedback=1"):
        _require(fragment in cfg, f"{binding['path']}: feedback diagnostic configuration drifted")
    return {"telemetry": telemetry, "shock_diagnostic_counts": diagnostic_counts}


def _verify_marker(payload: bytes, artifact: bytes, label: str) -> None:
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as error:
        raise PilotAnalysisError(f"{label}: restart completion marker is not ASCII") from error
    lines = text.splitlines()
    _require(len(lines) == 3 and lines[0] == "ATHENAK_RESTART_COMPLETE_V1",
             f"{label}: restart completion marker is malformed")
    _require(lines[1] == f"size={len(artifact)}", f"{label}: restart completion size drifted")
    expected = _fnv1a64(artifact)
    _require(lines[2] == f"fnv1a64={expected}", f"{label}: restart completion digest drifted")


def _restart_report(root: Path, case_id: str, restart: Mapping[str, Any]) -> dict[str, Any]:
    manifest_binding = restart["manifest"]
    manifest_payload = _member_payload(root, manifest_binding, f"{case_id} restart manifest")
    manifest_marker = _member_payload(root, restart["manifest_complete"], f"{case_id} restart manifest marker")
    _verify_marker(manifest_marker, manifest_payload, f"{case_id} restart manifest marker")
    decoded = _object(_decode_json(manifest_payload, f"{case_id} restart manifest"), {"schema", "members"}, f"{case_id} restart manifest")
    _strict_equal(decoded["schema"], "ATHENAK_RESTART_MANIFEST_V1", f"{case_id} restart schema")
    raw_members = _list(decoded["members"], f"{case_id} restart manifest/members")
    _require(bool(raw_members), f"{case_id}: restart manifest members are empty")
    declared = []
    expected_name = f"{case_id}.00004.rst"
    ranks = []
    shared = False
    for index, raw in enumerate(raw_members):
        label = f"{case_id} restart manifest/members[{index}]"
        member = _object(raw, {"path", "size", "fnv1a64"}, label)
        path = _relative_path(member["path"], f"{label}/path")
        size = _int(member["size"], f"{label}/size")
        digest = _text(member["fnv1a64"], f"{label}/fnv1a64")
        _require(size >= 0 and _FNV1A64_PATTERN.fullmatch(digest) is not None, f"{label}: digest metadata is invalid")
        relative = PurePosixPath(path)
        _require(relative.parts[0] == "rst" and relative.name == expected_name, f"{label}: restart member path drifted")
        if len(relative.parts) == 2:
            shared = True
        else:
            _require(len(relative.parts) == 3, f"{label}: restart member layout drifted")
            match = _RANK_DIRECTORY_PATTERN.fullmatch(relative.parts[1])
            _require(match is not None, f"{label}: restart rank directory drifted")
            ranks.append(int(match.group(1)))
        declared.append({"path": f"cases/{case_id}/{path}", "size": size, "fnv1a64": digest})
    _require(not (shared and len(declared) != 1), f"{case_id}: shared restart manifest must contain one member")
    _require(not (shared and ranks), f"{case_id}: restart manifest mixes shared and ranked members")
    if ranks:
        _require(ranks == list(range(len(ranks))), f"{case_id}: restart rank members are not contiguous")
    outer = restart["members"]
    _require(len(outer) == len(declared), f"{case_id}: restart product/member count drifted")
    for index, (member, expected) in enumerate(zip(outer, declared)):
        artifact_binding = member["artifact"]
        marker_binding = member["complete"]
        _strict_equal(artifact_binding["path"], expected["path"], f"{case_id} restart member[{index}] path")
        _strict_equal(marker_binding["path"], expected["path"] + ".complete", f"{case_id} restart member[{index}] marker path")
        artifact = _member_payload(root, artifact_binding, f"{case_id} restart member[{index}]")
        marker = _member_payload(root, marker_binding, f"{case_id} restart member[{index}] marker")
        _verify_marker(marker, artifact, f"{case_id} restart member[{index}] marker")
        _require(len(artifact) == expected["size"] and _fnv1a64(artifact) == expected["fnv1a64"],
                 f"{case_id}: restart manifest member digest drifted")
    return {"layout": "shared" if shared else "rank_sharded", "member_count": len(declared)}


def _snapshot_report(root: Path, case_id: str, ps_p0: float, snapshot: Mapping[str, Any]) -> dict[str, Any]:
    time = snapshot["time"]
    mhd_dataset, mhd = _binary_dataset(root, snapshot["mhd_w_bcc"], time, ps_p0, _MHD_FIELDS)
    _, bmag_product = _binary_dataset(root, snapshot["bmag"], time, ps_p0, _SCALAR_BIN_FIELDS["bmag"])
    _binary_dataset(root, snapshot["prtcl_jx"], time, ps_p0, _SCALAR_BIN_FIELDS["prtcl_jx"])
    _binary_dataset(root, snapshot["j2"], time, ps_p0, _SCALAR_BIN_FIELDS["j2"])
    particles = _particle_report(root, snapshot["prtcl_all"], time)
    rho = mhd["dens"].values
    pressure = (float(mhd_dataset.input_parameters["mhd"]["gamma"]) - 1.0) * mhd["eint"].values
    vx = mhd["velx"].values
    bmag = np.sqrt(mhd["bcc1"].values**2 + mhd["bcc2"].values**2 + mhd["bcc3"].values**2)
    _require(np.all(np.isfinite(pressure)) and np.all(pressure >= 0.0),
             f"{case_id} t={time}: pressure is non-finite or negative")
    _require(np.allclose(bmag, bmag_product["bmag"].values, rtol=1.0e-6, atol=1.0e-7),
             f"{case_id} t={time}: raw bmag disagrees with mhd_w_bcc")
    x1 = 0.5 * (mhd["dens"].x1_faces[:-1] + mhd["dens"].x1_faces[1:])
    profile = {
        "case_id": case_id,
        "ps_p0": ps_p0,
        "time_omega0_inverse": time,
        "x1_c_over_omega_pi": x1.tolist(),
        "rho_y_average": np.mean(rho, axis=(0, 1)).tolist(),
        "p_y_average": np.mean(pressure, axis=(0, 1)).tolist(),
        "vx_y_average": np.mean(vx, axis=(0, 1)).tolist(),
        "bmag_y_average": np.mean(bmag, axis=(0, 1)).tolist(),
    }
    return {"profile": profile, "particles": particles}


def analyze_pressure_pilot_bundle(
    root: str | Path,
    expected_manifest_sha256: str,
    *,
    authorized_publication_root: Path = AUTHORIZED_PUBLICATION_ROOT,
) -> dict[str, Any]:
    """Validate one complete four-case pilot bundle and emit overlay-ready records."""
    execution.validate_source_tranche()
    _require(_SHA256_PATTERN.fullmatch(expected_manifest_sha256) is not None,
             "expected manifest SHA-256 must contain 64 lowercase hexadecimal digits")
    bundle_root = Path(root)
    _require(bundle_root.is_absolute(), "pressure-pilot bundle root must be absolute")
    bundle_lexical = Path(os.path.abspath(bundle_root))
    try:
        bundle_root = bundle_lexical.resolve(strict=True)
        publication_root = Path(os.path.abspath(authorized_publication_root)).resolve(
            strict=True
        )
    except OSError as error:
        raise PilotAnalysisError(
            "pressure-pilot bundle or authorized publication root is unavailable"
        ) from error
    _require(
        bundle_root == bundle_lexical
        and publication_root == Path(os.path.abspath(authorized_publication_root))
        and bundle_root.parent == publication_root,
        "pressure-pilot bundle is outside the authorized publication root",
    )
    _require(bundle_root.is_dir(), "pressure-pilot bundle root must be a directory")
    policy = _load_policy()
    manifest_payload = _regular_bytes(bundle_root / MANIFEST_NAME, "pressure-pilot manifest")
    _require(_sha256_bytes(manifest_payload) == expected_manifest_sha256,
             "pressure-pilot manifest SHA-256 drifted")
    manifest = _manifest_schema(manifest_payload)
    _strict_equal(manifest["active_deck_binding"], policy["active_deck_binding"], "manifest/active deck binding")
    expected_prereg = {
        "path": PREREGISTRATION_PATH.relative_to(REPO_ROOT).as_posix(),
        "sha256": _sha256_bytes(_regular_bytes(PREREGISTRATION_PATH, "pressure-pilot preregistration")),
    }
    _strict_equal(manifest["preregistration_binding"], expected_prereg, "manifest/preregistration binding")
    expected_registered_execution = {
        "path": REGISTERED_EXECUTION_PREREGISTRATION_PATH.relative_to(
            REPO_ROOT
        ).as_posix(),
        "sha256": _sha256_bytes(
            _regular_bytes(
                REGISTERED_EXECUTION_PREREGISTRATION_PATH,
                "pressure-pilot registered-execution preregistration",
            )
        ),
    }
    _strict_equal(
        manifest["registered_execution_preregistration_binding"],
        expected_registered_execution,
        "manifest/registered-execution preregistration binding",
    )
    _validate_case_identity(manifest["cases"])
    _require(_actual_files(bundle_root) == _declared_paths(manifest),
             "pressure-pilot raw product inventory drifted")
    profiles = []
    case_summaries = []
    for case in manifest["cases"]:
        snapshot_reports = [
            _snapshot_report(bundle_root, case["case_id"], case["ps_p0"], snapshot)
            for snapshot in case["snapshots"]
        ]
        profiles.extend(report["profile"] for report in snapshot_reports)
        terminal_particles = snapshot_reports[-1]["particles"]["particle_count"]
        stdout = _stdout_report(bundle_root, case["stdout"], terminal_particles)
        restart = _restart_report(bundle_root, case["case_id"], case["terminal_restart"])
        case_summaries.append({
            "case_id": case["case_id"],
            "ps_p0": case["ps_p0"],
            "snapshot_count": len(snapshot_reports),
            "terminal_particles": snapshot_reports[-1]["particles"],
            "q017_telemetry": stdout["telemetry"],
            "shock_diagnostic_counts": stdout["shock_diagnostic_counts"],
            "terminal_restart": restart,
        })
    return {
        "schema_version": 1,
        "record_type": RESULT_RECORD_TYPE,
        "evidence_class": EVIDENCE_CLASS,
        "qualification_effect": QUALIFICATION_EFFECT,
        "scientific_evidence_eligible": False,
        "sun_bai_claim": False,
        "frontier_execution_authorized": False,
        "status": "pass_engineering_calibration_only",
        "manifest_sha256": expected_manifest_sha256,
        "overlay_profile_fields": ["rho", "p", "vx", "|B|"],
        "overlay_profiles": profiles,
        "case_summaries": case_summaries,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("bundle_root")
    parser.add_argument("expected_manifest_sha256")
    args = parser.parse_args(argv)
    try:
        result = analyze_pressure_pilot_bundle(args.bundle_root, args.expected_manifest_sha256)
    except (OSError, PilotAnalysisError, ValueError) as error:
        result = {
            "schema_version": 1,
            "record_type": RESULT_RECORD_TYPE,
            "evidence_class": EVIDENCE_CLASS,
            "qualification_effect": QUALIFICATION_EFFECT,
            "scientific_evidence_eligible": False,
            "sun_bai_claim": False,
            "frontier_execution_authorized": False,
            "status": "rejected",
            "failure": str(error),
        }
        print(json.dumps(result, indent=2, sort_keys=True))
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
