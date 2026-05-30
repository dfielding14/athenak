#!/usr/bin/env python3
"""Qualify deterministic serial relativistic CR typed-v2 restart behavior."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from typing import Iterator
from typing import Optional
from typing import Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputs/unit_tests/cr_relativistic_restart_runtime.athinput"
DEFAULT_CRITERIA = (
    REPO_ROOT
    / "docs/source/modules/figures/cr_tracer_relativistic_acceleration"
    / "cr_relativistic_phase6_preregistered_criteria.json"
)
EXPECTED_CRITERIA_SHA256 = "123701e93374eb61bdd5921f6d4b07eb3899a3aaabae36db632930da980fb269"
EXPECTED_CRITERION_IDS = tuple(f"P6-{number:02d}" for number in range(1, 22))
sys.path.insert(0, str(SCRIPT_DIR))
import cr_tracer_inspect as tracer_inspect  # noqa: E402


@dataclass(frozen=True)
class Checkpoint:
    """Paths for one copied or generated serial paired restart checkpoint."""

    root: Path
    number: int
    shard_relative: Path
    manifest_relative: Path
    mesh_relative: Path
    witness_relative: Path

    @property
    def shard(self) -> Path:
        return self.root / self.shard_relative

    @property
    def manifest(self) -> Path:
        return self.root / self.manifest_relative

    @property
    def mesh(self) -> Path:
        return self.root / self.mesh_relative

    @property
    def witness(self) -> Path:
        return self.root / self.witness_relative


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _validate_criteria(report: dict, criteria_path: Path) -> None:
    criteria_blob = criteria_path.read_bytes()
    criteria_sha256 = hashlib.sha256(criteria_blob).hexdigest()
    _require(
        criteria_sha256 == EXPECTED_CRITERIA_SHA256,
        f"{criteria_path}: frozen criteria SHA-256 mismatch",
    )
    criteria = json.loads(criteria_blob)
    _require(
        criteria.get("manifest_schema_version") == 1,
        f"{criteria_path}: unsupported criteria manifest schema",
    )
    _require(
        criteria.get("status") ==
        "reviewer_rebound_freeze_after_third_restart_specialist_hold_before_corrected_replay",
        f"{criteria_path}: criteria status is not the expected reviewer-rebound freeze",
    )
    _require(
        isinstance(criteria.get("criteria"), list) and criteria["criteria"],
        f"{criteria_path}: criteria list must be non-empty",
    )
    metrics = report["metrics"]
    operators = {
        "==": lambda actual, limit: actual == limit,
        ">=": lambda actual, limit: actual >= limit,
        "<=": lambda actual, limit: actual <= limit,
    }
    criterion_ids = [criterion.get("id") for criterion in criteria["criteria"]]
    _require(
        all(isinstance(criterion_id, str) and criterion_id for criterion_id in criterion_ids),
        f"{criteria_path}: every criterion requires a non-empty ID",
    )
    _require(
        len(criterion_ids) == len(set(criterion_ids)),
        f"{criteria_path}: criterion IDs must be unique",
    )
    _require(
        tuple(criterion_ids) == EXPECTED_CRITERION_IDS,
        f"{criteria_path}: criterion ID set does not match the frozen registration",
    )
    for criterion in criteria["criteria"]:
        metric = criterion["metric"]
        _require(metric in metrics, f"{criteria_path}: missing metric {metric!r}")
        operator = criterion["operator"]
        _require(operator in operators, f"{criteria_path}: unsupported operator {operator!r}")
        _require(
            operators[operator](metrics[metric], criterion["limit"]),
            f"{criterion['id']}: {metric}={metrics[metric]!r} does not satisfy "
            f"{operator} {criterion['limit']!r}",
        )
    report["criteria"] = {
        "path": str(criteria_path),
        "count": len(criteria["criteria"]),
        "status": "PASS",
        "sha256": criteria_sha256,
    }


def _manifest_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    shard_rows = 0
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if fields[0] == "shard":
            shard_rows += 1
            continue
        if len(fields) != 2:
            raise RuntimeError(f"{path}:{line_number}: malformed manifest row")
        key, value = fields
        if key in values:
            raise RuntimeError(f"{path}:{line_number}: duplicate manifest key {key!r}")
        values[key] = value
    _require(shard_rows == 1, f"{path}: Phase-6 serial qualification expected one shard")
    return values


def _checkpoint_from_shard(root: Path, shard: Path) -> Checkpoint:
    match = re.search(r"\.(\d{5})\.prst$", shard.name)
    _require(match is not None, f"{shard}: missing five-digit checkpoint number")
    manifest = root / "prst" / f"{shard.name}.manifest"
    _require(manifest.is_file(), f"{shard}: missing manifest {manifest}")
    values = _manifest_values(manifest)
    _require("mesh_restart" in values, f"{manifest}: missing mesh_restart")
    mesh_relative = Path(values["mesh_restart"])
    _require(not mesh_relative.is_absolute(), f"{manifest}: absolute mesh restart path")
    witness_relative = Path(str(mesh_relative) + ".rmeta")
    return Checkpoint(
        root=root,
        number=int(match.group(1)),
        shard_relative=shard.relative_to(root),
        manifest_relative=manifest.relative_to(root),
        mesh_relative=mesh_relative,
        witness_relative=witness_relative,
    )


def _discover_checkpoints(root: Path) -> list[Checkpoint]:
    shards = sorted((root / "prst/rank_00000000").glob("*.prst"))
    _require(shards, f"{root}: no rank-zero typed-v2 particle restart shards")
    checkpoints = [_checkpoint_from_shard(root, shard) for shard in shards]
    numbers = [checkpoint.number for checkpoint in checkpoints]
    _require(len(numbers) == len(set(numbers)), f"{root}: duplicate checkpoint numbers")
    return checkpoints


def _inspect_checkpoint(checkpoint: Checkpoint) -> dict:
    """Use the shared particle inspector to validate the complete typed-v2 bundle."""
    decoded = tracer_inspect.read_prst_file(checkpoint.shard)
    _require(
        decoded["format_version"] == "AKPRST-v2.0",
        f"{checkpoint.shard}: expected typed-v2 output, found {decoded['format_version']}",
    )
    _require(decoded["rank"] == 0, f"{checkpoint.shard}: expected serial rank zero")
    _require(checkpoint.mesh.is_file(), f"{checkpoint.shard}: missing paired mesh restart")
    _require(checkpoint.witness.is_file(), f"{checkpoint.shard}: missing paired mesh witness")
    return decoded


def _particle_signature(decoded: dict) -> tuple:
    return tuple(
        (
            particle["pgid"],
            particle["tag"],
            particle["species"],
            tuple(particle["reals"]),
        )
        for particle in decoded["particles"]
    )


def _native_mesh_state_payload(path: Path) -> bytes:
    blob = path.read_bytes()
    marker = b"<par_end>\n"
    marker_offset = blob.find(marker)
    _require(marker_offset >= 0, f"{path}: native mesh restart parameter terminator is missing")
    return blob[marker_offset + len(marker):]


def _fresh_case(root: Path, name: str) -> Path:
    case_dir = root / name
    _require(not case_dir.exists(), f"{case_dir}: scratch case already exists")
    case_dir.mkdir(parents=True)
    return case_dir


def _write_log(case_dir: Path, result: subprocess.CompletedProcess[str]) -> str:
    output = result.stdout + result.stderr
    (case_dir / "athena.log").write_text(output)
    return output


def _run_success(binary: Path, case_dir: Path, arguments: Sequence[str],
                 timeout: float, label: str) -> None:
    result = subprocess.run(
        [str(binary), *arguments],
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    output = _write_log(case_dir, result)
    if result.returncode:
        raise RuntimeError(f"{label}: Athena failed with exit code {result.returncode}:\n{output}")


def _run_expected_failure(binary: Path, case_dir: Path, arguments: Sequence[str],
                          timeout: float, label: str, expected: str) -> dict[str, str]:
    result = subprocess.run(
        [str(binary), *arguments],
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    output = _write_log(case_dir, result)
    if result.returncode == 0 or expected.lower() not in output.lower():
        raise RuntimeError(
            f"{label}: expected a nonzero exit containing {expected!r}, got "
            f"exit code {result.returncode}:\n{output}"
        )
    return {"id": label, "status": "PASS", "expected_diagnostic": expected}


def _copy_checkpoint(checkpoint: Checkpoint, destination: Path) -> None:
    for relative in (
        checkpoint.shard_relative,
        checkpoint.manifest_relative,
        checkpoint.mesh_relative,
        checkpoint.witness_relative,
    ):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(checkpoint.root / relative, target)


def _copy_unique_checkpoints(checkpoints: Sequence[Checkpoint], destination: Path) -> None:
    seen: set[int] = set()
    for checkpoint in checkpoints:
        if checkpoint.number in seen:
            continue
        _copy_checkpoint(checkpoint, destination)
        seen.add(checkpoint.number)


def _restart_arguments(particle_checkpoint: Checkpoint,
                       mesh_checkpoint: Optional[Checkpoint] = None,
                       extra: Sequence[str] = (),
                       particle_shard_relative: Optional[Path] = None) -> list[str]:
    mesh_checkpoint = mesh_checkpoint or particle_checkpoint
    particle_shard_relative = particle_shard_relative or particle_checkpoint.shard_relative
    return [
        "-r",
        str(mesh_checkpoint.mesh_relative),
        "problem/prtcl_rst_flag=1",
        f"problem/prtcl_res_file={particle_shard_relative}",
        *extra,
    ]


def _replace_manifest_value(path: Path, key: str,
                            replacement: Callable[[str], str]) -> None:
    lines = path.read_text().splitlines()
    replaced = 0
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) == 2 and fields[0] == key:
            lines[index] = f"{key} {replacement(fields[1])}"
            replaced += 1
    _require(replaced == 1, f"{path}: expected exactly one {key!r} row")
    path.write_text("\n".join(lines) + "\n")


def _mutate_payload(case_dir: Path, checkpoint: Checkpoint) -> None:
    path = case_dir / checkpoint.shard_relative
    blob = bytearray(path.read_bytes())
    _require(
        len(blob) > tracer_inspect.TYPED_V2_HEADER_BYTES,
        f"{path}: no payload byte available for corruption control",
    )
    blob[tracer_inspect.TYPED_V2_HEADER_BYTES] ^= 0x01
    path.write_bytes(blob)


def _mutate_mesh_payload(case_dir: Path, checkpoint: Checkpoint) -> None:
    path = case_dir / checkpoint.mesh_relative
    blob = bytearray(path.read_bytes())
    _require(blob, f"{path}: no mesh byte available for corruption control")
    blob[-1] ^= 0x01
    path.write_bytes(blob)


def _fabricate_checkpoint_id(case_dir: Path, checkpoint: Checkpoint) -> None:
    for relative in (checkpoint.manifest_relative, checkpoint.witness_relative):
        _replace_manifest_value(
            case_dir / relative,
            "checkpoint_id",
            lambda value: str(int(value) + 1),
        )


def _underscore_mesh_byte_count(case_dir: Path, checkpoint: Checkpoint) -> None:
    for relative in (checkpoint.manifest_relative, checkpoint.witness_relative):
        _replace_manifest_value(
            case_dir / relative,
            "mesh_byte_count",
            lambda value: f"{int(value):_}",
        )


def _rename_overwidth_shard(case_dir: Path, checkpoint: Checkpoint) -> None:
    source = case_dir / checkpoint.shard_relative
    renamed = source.with_name(
        source.name.replace(f".{checkpoint.number:05d}.prst", ".100000.prst"))
    _require(renamed != source, f"{source}: overwidth rename pattern did not match")
    source.rename(renamed)


def _mutate_future_header(case_dir: Path, checkpoint: Checkpoint) -> None:
    path = case_dir / checkpoint.shard_relative
    blob = bytearray(path.read_bytes())
    _require(len(blob) >= 10, f"{path}: header is too short for version control")
    struct.pack_into("<H", blob, 8, 3)
    path.write_bytes(blob)


def _truncate_shard(case_dir: Path, checkpoint: Checkpoint) -> None:
    path = case_dir / checkpoint.shard_relative
    blob = path.read_bytes()
    _require(blob, f"{path}: empty shard cannot be truncated")
    path.write_bytes(blob[:-1])


def _remove_witness(case_dir: Path, checkpoint: Checkpoint) -> None:
    (case_dir / checkpoint.witness_relative).unlink()


def _stale_witness(case_dir: Path, checkpoint: Checkpoint,
                   stale_checkpoint: Checkpoint) -> None:
    shutil.copy2(stale_checkpoint.witness, case_dir / checkpoint.witness_relative)


def _corrupt_manifest_count(case_dir: Path, checkpoint: Checkpoint) -> None:
    _replace_manifest_value(
        case_dir / checkpoint.manifest_relative,
        "global_count",
        lambda value: str(int(value) + 1),
    )


def _negative_control(binary: Path, input_path: Path, root: Path, timeout: float,
                      name: str, particle_checkpoint: Checkpoint, expected: str,
                      *, mesh_checkpoint: Optional[Checkpoint] = None,
                      copied_checkpoints: Sequence[Checkpoint] = (),
                      mutate: Optional[Callable[[Path, Checkpoint], None]] = None,
                      extra: Sequence[str] = (),
                      particle_only: bool = False,
                      particle_shard_relative: Optional[Path] = None) -> dict[str, str]:
    case_dir = _fresh_case(root, name)
    checkpoints = [particle_checkpoint, *copied_checkpoints]
    if mesh_checkpoint is not None:
        checkpoints.append(mesh_checkpoint)
    _copy_unique_checkpoints(checkpoints, case_dir)
    if mutate is not None:
        mutate(case_dir, particle_checkpoint)
    if particle_only:
        arguments = [
            "-i",
            str(input_path),
            "problem/prtcl_rst_flag=1",
            f"problem/prtcl_res_file={particle_checkpoint.shard_relative}",
            *extra,
        ]
    else:
        arguments = _restart_arguments(
            particle_checkpoint, mesh_checkpoint, extra, particle_shard_relative)
    return _run_expected_failure(binary, case_dir, arguments, timeout, name, expected)


def _inspector_negative_control(root: Path, name: str, checkpoint: Checkpoint,
                                expected: str,
                                mutate: Callable[[Path, Checkpoint], None]) -> dict[str, str]:
    case_dir = _fresh_case(root, name)
    _copy_checkpoint(checkpoint, case_dir)
    mutate(case_dir, checkpoint)
    copied = Checkpoint(
        root=case_dir,
        number=checkpoint.number,
        shard_relative=checkpoint.shard_relative,
        manifest_relative=checkpoint.manifest_relative,
        mesh_relative=checkpoint.mesh_relative,
        witness_relative=checkpoint.witness_relative,
    )
    try:
        _inspect_checkpoint(copied)
    except ValueError as error:
        _require(
            expected.lower() in str(error).lower(),
            f"{name}: expected inspector rejection containing {expected!r}, got {error}",
        )
        return {"id": name, "status": "PASS", "expected_diagnostic": expected}
    raise RuntimeError(f"{name}: expected inspector rejection containing {expected!r}")


@contextlib.contextmanager
def _runtime_root(work_dir: Optional[Path]) -> Iterator[Path]:
    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="cr-rel-phase6-") as temporary:
            yield Path(temporary)
        return
    root = work_dir.resolve()
    if root.exists():
        _require(not any(root.iterdir()), f"{root}: --work-dir must be empty")
    else:
        root.mkdir(parents=True)
    yield root


def _qualify(binary: Path, input_path: Path, criteria_path: Path,
             root: Path, timeout: float) -> dict:
    source_dir = _fresh_case(root, "source_uninterrupted")
    _run_success(binary, source_dir, ["-i", str(input_path)], timeout, "source_uninterrupted")
    source_checkpoints = _discover_checkpoints(source_dir)
    _require(
        len(source_checkpoints) >= 2,
        "source_uninterrupted: expected at least two paired restart checkpoints",
    )
    source_decoded = {
        checkpoint.number: _inspect_checkpoint(checkpoint)
        for checkpoint in source_checkpoints
    }
    base = source_checkpoints[0]
    reference = source_checkpoints[-1]

    continuation_dir = _fresh_case(root, "paired_continuation")
    _copy_checkpoint(base, continuation_dir)
    _run_success(
        binary,
        continuation_dir,
        _restart_arguments(base),
        timeout,
        "paired_continuation",
    )
    continuation_checkpoints = _discover_checkpoints(continuation_dir)
    continuation_by_number = {
        checkpoint.number: checkpoint for checkpoint in continuation_checkpoints
    }
    _require(
        reference.number in continuation_by_number,
        f"paired_continuation: did not regenerate checkpoint {reference.number:05d}",
    )
    continued = continuation_by_number[reference.number]
    continued_decoded = _inspect_checkpoint(continued)
    _require(
        _particle_signature(continued_decoded) == _particle_signature(source_decoded[reference.number]),
        "paired_continuation: decoded particle state differs from uninterrupted source",
    )
    _require(
        continued.shard.read_bytes() == reference.shard.read_bytes(),
        "paired_continuation: final particle shard is not byte-identical to uninterrupted source",
    )
    _require(
        _native_mesh_state_payload(continued.mesh) == _native_mesh_state_payload(reference.mesh),
        "paired_continuation: final mesh state payload is not byte-identical to uninterrupted source",
    )

    negative_controls = [
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "payload_corruption",
            base,
            "payload checksum mismatch",
            mutate=_mutate_payload,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "future_version_header_corruption",
            base,
            "schema version is unsupported",
            mutate=_mutate_future_header,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "truncated_shard",
            base,
            "shard byte count does not match manifest",
            mutate=_truncate_shard,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "missing_mesh_rmeta",
            base,
            "could not be opened",
            mutate=_remove_witness,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "stale_mesh_rmeta",
            base,
            "byte count does not match mesh checkpoint",
            mutate=lambda case_dir, checkpoint: _stale_witness(
                case_dir, checkpoint, reference
            ),
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "mismatched_mesh_particle_pair",
            base,
            "particle and mesh checkpoint numbers do not match",
            mesh_checkpoint=reference,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "manifest_count_corruption",
            base,
            "shard counts do not sum to global count",
            mutate=_corrupt_manifest_count,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "semantic_alpha_mismatch",
            base,
            "semantic configuration does not match the active input",
            extra=("particles/alpha_s=0.8",),
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "semantic_subcycle_mismatch",
            base,
            "semantic configuration does not match the active input",
            extra=("particles/subcycle_gyro_fraction=0.2",),
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "unpaired_particle_only_reload",
            base,
            "requires paired full-mesh -r and typed-v2 particle restart input",
            particle_only=True,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "mesh_payload_corruption",
            base,
            "mesh restart witness checksum does not match mesh checkpoint",
            mutate=_mutate_mesh_payload,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "fabricated_checkpoint_id",
            base,
            "mesh restart witness checkpoint ID is inconsistent",
            mutate=_fabricate_checkpoint_id,
        ),
        _negative_control(
            binary,
            input_path,
            root,
            timeout,
            "overwidth_particle_filename",
            base,
            "missing five-digit checkpoint number",
            mutate=_rename_overwidth_shard,
            particle_shard_relative=Path(
                "prst/rank_00000000/cr_relativistic_restart.100000.prst"),
        ),
    ]
    inspector_negative_controls = [
        _inspector_negative_control(
            root,
            "inspector_numeric_grammar",
            base,
            "invalid mesh_byte_count",
            _underscore_mesh_byte_count,
        ),
    ]
    rejected = {control["id"]: True for control in negative_controls}
    inspector_rejected = {control["id"]: True for control in inspector_negative_controls}

    report = {
        "schema_version": 1,
        "qualification": "cr_relativistic_restart_runtime_phase6",
        "status": "PASS",
        "binary": str(binary),
        "input": str(input_path),
        "source": {
            "checkpoint_numbers": [checkpoint.number for checkpoint in source_checkpoints],
            "checkpoint_count": len(source_checkpoints),
            "final_particle_count": source_decoded[reference.number]["count"],
            "final_format": source_decoded[reference.number]["format_version"],
        },
        "paired_continuation": {
            "reloaded_checkpoint": base.number,
            "compared_checkpoint": reference.number,
            "decoded_particle_state_equal": True,
            "particle_shard_byte_identical": True,
            "mesh_state_payload_byte_identical": True,
        },
        "negative_controls": negative_controls,
        "negative_control_count": len(negative_controls),
        "inspector_negative_controls": inspector_negative_controls,
        "inspector_negative_control_count": len(inspector_negative_controls),
        "metrics": {
            "source_checkpoint_count": len(source_checkpoints),
            "source_final_format": source_decoded[reference.number]["format_version"],
            "continuation_decoded_particle_state_equal": True,
            "continuation_particle_shard_byte_identical": True,
            "continuation_mesh_state_payload_byte_identical": True,
            "copied_artifact_negative_control_count": len(negative_controls),
            "inspector_negative_control_count": len(inspector_negative_controls),
            **{f"{name}_rejected": value for name, value in rejected.items()},
            **{f"{name}_rejected": value for name, value in inspector_rejected.items()},
        },
    }
    _validate_criteria(report, criteria_path)
    return report


def _write_json(path_text: Optional[str], report: dict) -> None:
    if path_text is None:
        return
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if path_text == "-":
        print(payload, end="")
        return
    path = Path(path_text).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(payload)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Qualify serial relativistic CR paired typed-v2 restart runtime behavior."
    )
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument(
        "--work-dir",
        type=Path,
        help="preserve isolated runtime cases in this new or empty directory",
    )
    parser.add_argument(
        "--json",
        metavar="PATH",
        help="write the machine-readable qualification report to PATH, or '-' for stdout",
    )
    parser.add_argument("--timeout", type=float, default=30.0)
    args = parser.parse_args()

    binary = args.binary.resolve()
    input_path = args.input.resolve()
    criteria_path = args.criteria.resolve()
    report = {
        "schema_version": 1,
        "qualification": "cr_relativistic_restart_runtime_phase6",
        "status": "FAIL",
        "binary": str(binary),
        "input": str(input_path),
        "criteria": str(criteria_path),
    }
    try:
        _require(binary.is_file(), f"{binary}: --binary is not a file")
        _require(input_path.is_file(), f"{input_path}: --input is not a file")
        _require(criteria_path.is_file(), f"{criteria_path}: --criteria is not a file")
        _require(args.timeout > 0.0, "--timeout must be positive")
        with _runtime_root(args.work_dir) as root:
            report = _qualify(binary, input_path, criteria_path, root, args.timeout)
    except Exception as error:  # Keep CLI failure output concise; case logs retain detail.
        report["error"] = str(error)
        _write_json(args.json, report)
        print(f"CR relativistic Phase-6 runtime qualification FAIL: {error}", file=sys.stderr)
        return 1

    _write_json(args.json, report)
    print("CR relativistic Phase-6 runtime qualification PASS")
    print(
        "  paired typed-v2 source checkpoints inspected: "
        f"{report['source']['checkpoint_count']}"
    )
    print(
        "  continuation comparison: checkpoint "
        f"{report['paired_continuation']['reloaded_checkpoint']:05d} -> "
        f"{report['paired_continuation']['compared_checkpoint']:05d}, "
        "decoded state equal and particle shard byte-identical"
    )
    print(f"  copied-artifact negative controls rejected: {report['negative_control_count']}")
    print(
        "  inspector negative controls rejected: "
        f"{report['inspector_negative_control_count']}"
    )
    print(f"  preregistered criteria passed: {report['criteria']['count']}")
    if args.json:
        print(f"  JSON report: {args.json}")
    if args.work_dir:
        print(f"  preserved runtime cases: {args.work_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
