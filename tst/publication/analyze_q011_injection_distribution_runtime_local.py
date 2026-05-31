#!/usr/bin/env python3
"""Bounded local Q-011 shock-injection distribution runtime audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

if __package__:
    from .immutable_orion_tree import freeze_tree as freeze_immutable_tree
    from .immutable_orion_tree import staged_verified_frozen_tree
    from .immutable_orion_tree import verify_frozen_tree as verify_immutable_tree
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    from immutable_orion_tree import freeze_tree as freeze_immutable_tree
    from immutable_orion_tree import staged_verified_frozen_tree
    from immutable_orion_tree import verify_frozen_tree as verify_immutable_tree
    from pvtk_particles import ParticleVTKData, read_particle_vtk


REPO_ROOT = Path(__file__).resolve().parents[2]
DECK = REPO_ROOT / "inputs/tests/pic_q011_injection_distribution_runtime_local.athinput"
SOURCE = REPO_ROOT / "src/pgen/tests/pic_parallel_shock.cpp"
EXPECTED_DECK_SHA256 = "6e84e34a91e48b4933f26ee1c3354e95139ff89ddf75ec1b358a7a7d3c1d1c17"
EXPECTED_SOURCE_SHA256 = "509fb2387e3b9933f919b062962158805bacce27e928fa396cd0583e082cbbf9"
ORION_BULK_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
ARTIFACT_ROLE = "bounded_serial_host_runtime_diagnostic_only"
QUALIFICATION_EFFECT = "none"
AUDIT_METHOD = "deterministic_bounded_statistics_not_exact_rng_replay"
INVENTORY_NAME = "artifact_inventory.sha256"
FREEZE_RECEIPT_NAME = "freeze_receipt.json"
_PVTK_EXECUTION_PATTERN = re.compile(
    rb"^# vtk DataFile Version 2\.0\n"
    rb"# AthenaK particle data at time= ([^ \n]+)  "
    rb"nranks= (0|[1-9][0-9]*)  cycle=(0|[1-9][0-9]*)  variables=([^\n]+)\n"
)

_EXPECTED_DECK_VALUES = {
    ("job", "basename"): "pic_q011_injection_distribution_runtime_local",
    ("mesh", "nx1"): "8",
    ("mesh", "x1min"): "0.0",
    ("mesh", "x1max"): "8.0",
    ("mesh", "ix1_bc"): "reflect",
    ("mesh", "ox1_bc"): "inflow",
    ("mesh", "nx2"): "16",
    ("mesh", "x2min"): "0.0",
    ("mesh", "x2max"): "16.0",
    ("mesh", "ix2_bc"): "periodic",
    ("mesh", "ox2_bc"): "periodic",
    ("mesh", "nx3"): "1",
    ("meshblock", "nx1"): "8",
    ("meshblock", "nx2"): "16",
    ("meshblock", "nx3"): "1",
    ("mesh_refinement", "refinement"): "none",
    ("time", "integrator"): "rk1",
    ("time", "nlim"): "1",
    ("particles", "ppc"): "0.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "cr_distribution"): "random",
    ("particles", "deposit_qscale"): "2.0e-2",
    ("particles", "pic_physical_mode"): "paper_mhd_pic",
    ("particles", "pic_cr_initial_state"): "momentum",
    ("problem", "pgen_name"): "pic_parallel_shock",
    ("problem", "ps_rho0"): "1.0",
    ("problem", "ps_p0"): "0.10",
    ("problem", "ps_u0"): "30.0",
    ("problem", "ps_eta"): "1.0",
    ("problem", "ps_vinj_over_u0"): "3.16227766017",
    ("problem", "ps_inject_half_width_cells"): "0.5",
    ("problem", "ps_shock_speed_model"): "ideal_surface",
    ("problem", "ps_enable_injection"): "true",
    ("problem", "ps_enable_gas_subtraction"): "false",
    ("problem", "ps_enable_curvature_amr"): "false",
    ("problem", "ps_inject_seed"): "11011",
    ("problem", "ps_enable_frame_tracking"): "false",
    ("q011_injection_distribution_runtime_local", "audit_role"): ARTIFACT_ROLE,
    ("q011_injection_distribution_runtime_local", "qualification_effect"):
        QUALIFICATION_EFFECT,
    ("q011_injection_distribution_runtime_local", "frontier_authorization"):
        "not_bound",
    ("q011_injection_distribution_runtime_local", "amr_qualification"):
        "not_claimed",
    ("q011_injection_distribution_runtime_local", "gpu_qualification"):
        "not_claimed",
    ("q011_injection_distribution_runtime_local", "physical_calibration"):
        "not_claimed",
    ("q011_injection_distribution_runtime_local", "external_review"):
        "not_claimed",
    ("q011_injection_distribution_runtime_local", "audit_method"): AUDIT_METHOD,
    ("output1", "file_type"): "pvtk",
    ("output1", "variable"): "prtcl_all",
    ("output1", "dcycle"): "1",
}

_EXPLICIT_GAPS = [
    "physical gas-pressure, unit-normalization, macro-mass and downstream-ppc calibration",
    "AMR versus fine-uniform residual qualification",
    "GPU, HIP, MPI, decomposition and Frontier qualification",
    "independent raw-artifact recompute and external review",
]


class AuditError(ValueError):
    """Raised when the bounded Q-011 runtime audit fails closed."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_athinput(path: Path = DECK) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the Q-011 local deck."""
    blocks: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            if not line.endswith(">"):
                raise AuditError(f"{path}:{lineno}: malformed block header")
            name = line[1:-1].strip()
            if not name:
                raise AuditError(f"{path}:{lineno}: empty block name")
            current = blocks.setdefault(name, {})
            continue
        if current is None or "=" not in line:
            raise AuditError(f"{path}:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        if not name or not value:
            raise AuditError(f"{path}:{lineno}: empty parameter name or value")
        if name in current:
            raise AuditError(f"{path}:{lineno}: duplicate parameter {name}")
        current[name] = value
    return blocks


def validate_deck(path: Path = DECK) -> dict[str, Any]:
    """Require the reduced serial-host role and preserve every nonclaim."""
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise AuditError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )
    if _sha256(path) != EXPECTED_DECK_SHA256:
        raise AuditError(f"{path}: exact bounded Q-011 deck SHA-256 drifted")
    gamma = float(blocks["mhd"]["gamma"])
    u0 = float(blocks["problem"]["ps_u0"])
    x1min = float(blocks["mesh"]["x1min"])
    x1max = float(blocks["mesh"]["x1max"])
    shock_speed = 0.5 * (gamma - 1.0) * u0
    injection_speed = float(blocks["problem"]["ps_vinj_over_u0"]) * u0
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": _sha256(path),
        "audit_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "audit_method": AUDIT_METHOD,
        "shock_speed": shock_speed,
        "injection_speed_relative_to_surface": injection_speed,
        "shock_surface_at_injection_time": x1min,
        "clamped_surface_output_x1": x1min + max(1.0e-12 * (x1max - x1min), 1.0e-14),
        "explicit_gaps": list(_EXPLICIT_GAPS),
    }


def validate_source_contract() -> dict[str, Any]:
    """Bind the audit to the committed full-sphere surface-injection path."""
    source_sha256 = _sha256(SOURCE)
    if source_sha256 != EXPECTED_SOURCE_SHA256:
        raise AuditError("pic_parallel_shock source SHA-256 drifted")
    source = SOURCE.read_text(encoding="utf-8")
    snippets = [
        "std::mt19937_64 rng(seed);",
        "const Real mu = 2.0*Uniform01(rng) - 1.0;",
        "const Real phi = 2.0*M_PI*Uniform01(rng);",
        "dirx = mu;",
        "diry = st*std::cos(phi);",
        "dirz = st*std::sin(phi);",
        "part.x1 = xshock;",
        "xshock >= x1c - half_width && xshock < x1c + half_width",
        "part.vx = ps_shock_speed + frame_vx + vinj*dirx;",
        "part.vy = vinj*diry;",
        "part.vz = vinj*dirz;",
    ]
    missing = [snippet for snippet in snippets if snippet not in source]
    if missing:
        raise AuditError("pic_parallel_shock source contract drift:\n" + "\n".join(missing))
    return {
        "path": str(SOURCE.relative_to(REPO_ROOT)),
        "sha256": source_sha256,
        "sampler": "full_sphere_isotropic_monoenergetic_relative_to_ideal_surface",
        "placement": "single_half_open_carrier_cell_with_x1_at_clamped_surface",
    }


def _scalar(data: ParticleVTKData, name: str) -> np.ndarray:
    if name not in data.scalars:
        raise AuditError(f"Particle VTK scalar {name!r} is required")
    values = np.asarray(data.scalars[name])
    if values.shape != (data.points.shape[0],):
        raise AuditError(f"Particle VTK scalar {name!r} has the wrong shape")
    return values


def _integer_scalar(data: ParticleVTKData, name: str) -> np.ndarray:
    """Return an exact integer scalar without accepting lossy truncation."""
    values = _scalar(data, name)
    if np.issubdtype(values.dtype, np.integer):
        return values.astype(np.int64)
    _require(np.all(np.isfinite(values)), f"Particle VTK scalar {name!r} must be finite")
    rounded = np.rint(values)
    _require(
        np.array_equal(values, rounded),
        f"Particle VTK scalar {name!r} must contain exact integers",
    )
    return rounded.astype(np.int64)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def analyze_particle_payload(data: ParticleVTKData) -> dict[str, Any]:
    """Audit one untransported runtime injection payload using deterministic bounds."""
    deck = validate_deck()
    validate_source_contract()
    points = np.asarray(data.points, dtype=np.float64)
    velocity = np.asarray(data.vectors.get("vel"), dtype=np.float64)
    _require(points.ndim == 2 and points.shape[1:] == (3,), "points shape must be (n, 3)")
    _require(velocity.shape == points.shape, "velocity shape must match points")
    _require(np.all(np.isfinite(points)), "points must contain finite values")
    _require(np.all(np.isfinite(velocity)), "velocity must contain finite values")
    count = points.shape[0]
    _require(128 <= count <= 1024, "runtime injection count must remain bounded in [128, 1024]")

    tags = _integer_scalar(data, "ptag")
    species = _integer_scalar(data, "species")
    source = _integer_scalar(data, "cr_source")
    birth_time = _scalar(data, "birth_time").astype(np.float64)
    _require(np.all(np.isfinite(birth_time)), "runtime birth times must be finite")
    _require(np.array_equal(tags, np.arange(count)), "runtime tags must be serial from zero")
    _require(np.all(species == 0), "runtime injection must use species zero")
    _require(np.all(source == 1), "runtime particles must be shock_injected provenance")
    _require(np.all(birth_time == 0.0), "runtime particles must be born at source time zero")

    x1_expected = float(deck["clamped_surface_output_x1"])
    surface_residual = np.abs(points[:, 0] - x1_expected)
    _require(
        float(np.max(surface_residual)) <= 1.0e-17,
        "runtime particle x1 placement drifted from the clamped shock surface",
    )
    _require(np.all((points[:, 1] >= 0.0) & (points[:, 1] < 16.0)),
             "runtime particle x2 placement left the carrier domain")
    _require(np.all(points[:, 2] == 0.0), "thin 2D PVTK output must report x3=0")

    shock_speed = float(deck["shock_speed"])
    injection_speed = float(deck["injection_speed_relative_to_surface"])
    relative_velocity = velocity.copy()
    relative_velocity[:, 0] -= shock_speed
    relative_speed = np.linalg.norm(relative_velocity, axis=1)
    speed_residual = np.abs(relative_speed - injection_speed)
    _require(
        float(np.max(speed_residual)) <= 1.0e-5,
        "runtime relative speed drifted from the monoenergetic shell",
    )
    direction = relative_velocity / relative_speed[:, None]
    direction_mean = np.mean(direction, axis=0)
    direction_second_moment = np.mean(direction * direction, axis=0)
    _require(
        bool(np.all(np.abs(direction_mean) <= 0.12)),
        "runtime direction means exceeded the bounded isotropy diagnostic",
    )
    _require(
        bool(np.all((direction_second_moment >= 0.25) &
                    (direction_second_moment <= 0.42))),
        "runtime direction second moments exceeded the bounded isotropy diagnostic",
    )
    octant_index = (
        (direction[:, 0] >= 0.0).astype(np.int64) * 4
        + (direction[:, 1] >= 0.0).astype(np.int64) * 2
        + (direction[:, 2] >= 0.0).astype(np.int64)
    )
    octant_counts = np.bincount(octant_index, minlength=8)
    _require(bool(np.all(octant_counts > 0)), "runtime sampler did not populate all octants")

    return {
        "schema_version": 1,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "audit_method": AUDIT_METHOD,
        "particle_count": count,
        "provenance": {
            "cr_source": "shock_injected",
            "species": 0,
            "birth_time": 0.0,
            "serial_tag_min": int(tags[0]),
            "serial_tag_max": int(tags[-1]),
        },
        "shock_surface_placement": {
            "model_surface_x1": float(deck["shock_surface_at_injection_time"]),
            "clamped_output_x1": x1_expected,
            "maximum_absolute_x1_residual": float(np.max(surface_residual)),
            "single_surface_coordinate": bool(np.unique(points[:, 0]).size == 1),
        },
        "monoenergetic_full_sphere_sampler": {
            "shock_speed": shock_speed,
            "expected_relative_speed": injection_speed,
            "maximum_absolute_relative_speed_residual": float(np.max(speed_residual)),
            "direction_mean": direction_mean.tolist(),
            "direction_second_moment": direction_second_moment.tolist(),
            "octant_counts": octant_counts.astype(int).tolist(),
        },
        "explicit_gaps": list(_EXPLICIT_GAPS),
    }


def extract_runtime_artifact(
    pvtk_path: str | Path,
    *,
    artifact_root: str | Path,
    artifact_inventory_sha256: str,
    executable: str | Path,
    executable_root: str | Path,
    executable_root_inventory_sha256: str,
    expected_executable_sha256: str,
) -> dict[str, Any]:
    """Audit one retained payload only after both frozen trees verify."""
    runtime_root = Path(artifact_root)
    path = _require_frozen_tree_member(
        pvtk_path,
        root=runtime_root,
        label="Q-011 runtime PVTK",
    )
    executable_tree_root = Path(executable_root)
    executable_path = _require_frozen_tree_member(
        executable,
        root=executable_tree_root,
        label="Q-011 pinned executable",
    )
    with staged_verified_frozen_tree(
        runtime_root,
        artifact_inventory_sha256,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-011 retained runtime tree",
    ) as (runtime_tree, runtime_snapshot), staged_verified_frozen_tree(
        executable_tree_root,
        executable_root_inventory_sha256,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-011 pinned executable tree",
    ) as (executable_tree, executable_snapshot):
        staged_path = runtime_snapshot.member_path(path.relative_to(runtime_root))
        staged_executable = executable_snapshot.member_path(
            executable_path.relative_to(executable_tree_root)
        )
        executable_sha256 = _sha256(staged_executable)
        _require(
            executable_sha256 == expected_executable_sha256,
            "Q-011 pinned executable SHA-256 drifted",
        )
        execution_metadata = _read_pvtk_execution_metadata(staged_path)
        report = analyze_particle_payload(read_particle_vtk(staged_path))
        report["immutable_runtime_artifact"] = {
            "path": str(path),
            "sha256": _sha256(staged_path),
        }
        report["immutable_runtime_tree"] = runtime_tree
        report["pvtk_execution_metadata"] = execution_metadata
        report["deck"] = validate_deck()
        report["source"] = validate_source_contract()
        report["host_executable"] = {
            "path": str(executable_path),
            "sha256": executable_sha256,
        }
        report["immutable_executable_tree"] = executable_tree
        return report


def _require_frozen_tree_member(
    member: str | Path,
    *,
    root: Path,
    label: str,
) -> Path:
    """Require one canonical real file strictly below an already-verified tree."""
    candidate = Path(member)
    if not candidate.is_absolute():
        raise AuditError(f"{label} path must be absolute")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, RuntimeError, ValueError) as error:
        raise AuditError(f"{label} must remain below {root}: {error}") from error
    if resolved != candidate:
        raise AuditError(f"{label} must use a canonical path without aliases")
    if resolved == root or not resolved.is_file():
        raise AuditError(f"{label} must be a retained regular file")
    return resolved


def _read_pvtk_execution_metadata(path: Path) -> dict[str, Any]:
    """Require the bounded serial cycle-one PVTK execution header."""
    with path.open("rb") as stream:
        header = stream.read(4096)
    match = _PVTK_EXECUTION_PATTERN.match(header)
    _require(match is not None, "PVTK execution header is missing or malformed")
    time = float(match.group(1))
    nranks = int(match.group(2))
    cycle = int(match.group(3))
    variables = match.group(4).decode("ascii")
    _require(math.isfinite(time), "PVTK execution time must be finite")
    _require(nranks == 1, "bounded runtime artifact must report nranks=1")
    _require(cycle == 1, "bounded runtime artifact must report cycle=1")
    _require(variables == "prtcl_all", "bounded runtime artifact must report variables=prtcl_all")
    return {
        "time": time,
        "nranks": nranks,
        "cycle": cycle,
        "variables": variables,
    }


def write_runtime_report(
    pvtk_path: str | Path,
    output_path: str | Path,
    *,
    artifact_root: str | Path,
    artifact_inventory_sha256: str,
    executable: str | Path,
    executable_root: str | Path,
    executable_root_inventory_sha256: str,
    expected_executable_sha256: str,
) -> dict[str, Any]:
    """Write the bounded nonqualifying runtime report."""
    report = extract_runtime_artifact(
        pvtk_path,
        artifact_root=artifact_root,
        artifact_inventory_sha256=artifact_inventory_sha256,
        executable=executable,
        executable_root=executable_root,
        executable_root_inventory_sha256=executable_root_inventory_sha256,
        expected_executable_sha256=expected_executable_sha256,
    )
    Path(output_path).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def verify_frozen_tree(root: str | Path, expected_inventory_sha256: str) -> dict[str, Any]:
    """Verify exact inventory hashes, tree membership, and recursive read-only modes."""
    return verify_immutable_tree(
        root,
        expected_inventory_sha256,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-011 retained runtime tree",
    )


def freeze_tree(root: str | Path) -> dict[str, Any]:
    """Write an exact inventory and recursively remove write bits from the tree."""
    receipt = {
        "schema_version": 1,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
    }
    return freeze_immutable_tree(
        root,
        receipt,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-011 retained runtime tree",
    )


def _write_json(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser("audit")
    audit.add_argument("pvtk")
    audit.add_argument("output_json")
    audit.add_argument("--artifact-root", required=True)
    audit.add_argument("--artifact-inventory-sha256", required=True)
    audit.add_argument("--executable", required=True)
    audit.add_argument("--executable-root", required=True)
    audit.add_argument("--executable-root-inventory-sha256", required=True)
    audit.add_argument("--expected-executable-sha256", required=True)
    freeze = subparsers.add_parser("freeze-tree")
    freeze.add_argument("root")
    verify = subparsers.add_parser("verify-frozen-tree")
    verify.add_argument("root")
    verify.add_argument("--expected-inventory-sha256", required=True)
    args = parser.parse_args()
    if args.command == "audit":
        _write_json(
            write_runtime_report(
                args.pvtk,
                args.output_json,
                artifact_root=args.artifact_root,
                artifact_inventory_sha256=args.artifact_inventory_sha256,
                executable=args.executable,
                executable_root=args.executable_root,
                executable_root_inventory_sha256=args.executable_root_inventory_sha256,
                expected_executable_sha256=args.expected_executable_sha256,
            )
        )
    elif args.command == "freeze-tree":
        _write_json(freeze_tree(args.root))
    else:
        _write_json(verify_frozen_tree(args.root, args.expected_inventory_sha256))


if __name__ == "__main__":
    main()
