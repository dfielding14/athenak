#!/usr/bin/env python3
"""Qualify bounded prescribed-test relativistic CR migration after fence removal."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator
from typing import Optional
from typing import Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputs/unit_tests/cr_relativistic_migration_runtime.athinput"
DEFAULT_CRITERIA = (
    REPO_ROOT
    / "docs/source/modules/figures/cr_tracer_relativistic_acceleration"
    / "cr_relativistic_phase8_preregistered_criteria.json"
)
EXPECTED_CRITERIA_SHA256 = (
    "add2ab12fbaf7f669e8e4eecc9e8f9cb0d581f06947c7e39876f934239708b74"
)
EXPECTED_CRITERION_IDS = tuple(f"P8-{number:02d}" for number in range(1, 63))
EXPECTED_CRITERIA_STATUS = "phase8_preimplementation_freeze_before_parent_fence_removal"
EXPECTED_PARTICLE_KEYS = tuple((0, tag) for tag in range(8))
STATIC_SMR_PARTICLE_KEYS = tuple((0, tag) for tag in range(15))
EXPECTED_RANK_LADDER = (1, 2, 4, 8)
EXACT_RANK_COUNT_REJECTION = "typed-v2 restart rejects MPI rank-count change"
TOPOLOGY_REJECTION = "typed-v2 restart rejects changed rank-to-MeshBlock topology"
RUNTIME_FENCE = (
    "Particle pusher 'relativistic_hc' execution currently requires a serial "
    "uniform-level mesh; MPI, SMR, and AMR remain closed until qualification"
)
STRICTLY_PERIODIC_REJECTION = (
    "Particle pusher 'relativistic_hc' currently requires a 3D strictly periodic mesh"
)
MHD_IDEAL_ONE_MESHBLOCK_REJECTION = (
    "relativistic_field_source = mhd_ideal currently requires exactly one MeshBlock"
)
MHD_IDEAL_RUNTIME_FENCE = (
    "relativistic_field_source = mhd_ideal currently requires a serial "
    "uniform-level mesh; MPI, SMR, and AMR remain closed until qualification"
)
PERFORMANCE_RE = re.compile(r"Particle performance ([^:]+):\s*(.*)")
CYCLE_RE = re.compile(
    r"elapsed=\S+ cycle=(\d+) time=(\S+) dt=(\S+)"
)
REAL_FIELDS = (
    "IPX", "IPVX", "IPY", "IPVY", "IPZ", "IPVZ", "IPM", "IPBX", "IPBY",
    "IPBZ", "IPDX", "IPDY", "IPDZ", "IPDB", "IPWX", "IPWY", "IPWZ", "IPCEX",
    "IPCEY", "IPCEZ", "IPWORK", "IPALPHA",
)
POSITION_INDICES = (0, 2, 4)
DISPLACEMENT_INDICES = (10, 11, 12)
MESH_MIN = (-0.5, -0.5, -0.5)
MESH_MAX = (0.5, 0.5, 0.5)
MULTIBLOCK_FLAGS = (
    "meshblock/nx1=8",
    "meshblock/nx2=8",
    "meshblock/nx3=8",
)
FINE_UNIFORM_FLAGS = (
    "mesh/nx1=32",
    "mesh/nx2=32",
    "mesh/nx3=32",
    "meshblock/nx1=8",
    "meshblock/nx2=8",
    "meshblock/nx3=8",
    "particles/ppc=0.000244140625",
)
SMR_APPEND = """

<mesh_refinement>
refinement = static

<refined_region1>
level = 1
x1min = 0.05
x1max = 0.15
x2min = 0.05
x2max = 0.15
x3min = 0.05
x3max = 0.15
"""
AMR_APPEND = """

<mesh_refinement>
refinement           = adaptive
ncycle_check         = 1
refinement_interval  = 1
num_levels           = 2
max_nmb_per_rank     = 128
particle_load_weight = 0.0

<amr_criterion0>
method = user
"""
STATIC_SMR_PPC = "0.001953125"  # one particle per leaf block: 15 total
STATIC_FINE_UNIFORM_PPC = "0.000457763671875"  # 15 particles / (64 blocks * 512 cells)
# Make the adaptive oracle gyro-limited below the finest-cell crossing bound so
# adaptive and fine-uniform legs advance through the same three outer steps.
ADAPTIVE_TIMESTEP_FLAGS = (
    "particles/subcycle_max_steps=4",
    "particles/subcycle_gyro_fraction=0.0005",
)
# Hold uninterrupted and restarted continuation legs to the same 0.1 outer
# steps. The split reference must not compare different integrator partitions.
CONTINUATION_TIMESTEP_FLAGS = (
    "particles/subcycle_max_steps=4",
    "particles/subcycle_cell_fraction=0.4",
)
EXPECTED_PRESCRIBED_FIELDS = {
    7: 0.0,
    8: 0.0,
    9: 0.35,
    17: 0.03,
    18: -0.02,
    19: 0.01,
}
sys.path.insert(0, str(SCRIPT_DIR))
import cr_tracer_inspect as tracer_inspect  # noqa: E402


@dataclass(frozen=True)
class ParticleState:
    """One decoded typed-v2 particle record."""

    pgid: int
    tag: int
    species: int
    reals: tuple[float, ...]

    @property
    def key(self) -> tuple[int, int]:
        return self.species, self.tag


@dataclass(frozen=True)
class Checkpoint:
    """One paired mesh restart plus typed-v2 particle-manifest bundle."""

    root: Path
    number: int
    manifest_relative: Path
    mesh_relative: Path
    witness_relative: Path
    shard_relatives: tuple[Path, ...]
    rank_counts: tuple[int, ...]
    topology_hash: int

    @property
    def manifest(self) -> Path:
        return self.root / self.manifest_relative

    @property
    def mesh(self) -> Path:
        return self.root / self.mesh_relative

    @property
    def witness(self) -> Path:
        return self.root / self.witness_relative

    @property
    def rank_zero_shard_relative(self) -> Path:
        for relative in self.shard_relatives:
            if "rank_00000000" in str(relative):
                return relative
        raise RuntimeError(f"{self.manifest}: no rank-zero shard")


@dataclass(frozen=True)
class Inspection:
    """Decoded aggregate state for one typed-v2 checkpoint."""

    checkpoint: Checkpoint
    particles: dict[tuple[int, int], ParticleState]


@dataclass(frozen=True)
class RuntimeCase:
    """One isolated Athena runtime case and its inspected checkpoints."""

    name: str
    root: Path
    log: str
    meshblocks: int
    created_blocks: int
    deleted_blocks: int
    inspections: tuple[Inspection, ...]

    @property
    def initial(self) -> Inspection:
        return self.inspections[0]

    @property
    def final(self) -> Inspection:
        return self.inspections[-1]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _flags(**values: object) -> list[str]:
    return [f"{key.replace('__', '/')}={value}" for key, value in values.items()]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()


def _validate_criteria_document(
    criteria_path: Path, *, enforce_digest: bool = True
) -> dict:
    blob = criteria_path.read_bytes()
    digest = hashlib.sha256(blob).hexdigest()
    if enforce_digest:
        _require(
            digest == EXPECTED_CRITERIA_SHA256,
            f"{criteria_path}: frozen criteria SHA-256 mismatch",
        )
    criteria = json.loads(blob)
    _require(
        criteria.get("manifest_schema_version") == 1,
        f"{criteria_path}: unsupported criteria manifest schema",
    )
    _require(
        criteria.get("status") == EXPECTED_CRITERIA_STATUS,
        f"{criteria_path}: criteria status is not the Phase-8 preimplementation freeze",
    )
    entries = criteria.get("criteria")
    _require(
        isinstance(entries, list) and entries,
        f"{criteria_path}: empty criteria list",
    )
    ids = tuple(entry.get("id") for entry in entries)
    _require(
        all(isinstance(criterion_id, str) and criterion_id for criterion_id in ids),
        f"{criteria_path}: every criterion requires a non-empty ID",
    )
    _require(len(ids) == len(set(ids)), f"{criteria_path}: criterion IDs must be unique")
    _require(
        ids == EXPECTED_CRITERION_IDS,
        f"{criteria_path}: criterion ID sequence does not match the frozen registration",
    )
    criteria["sha256"] = digest
    return criteria


def _validate_runtime_criteria(criteria_path: Path, report: dict) -> None:
    criteria = _validate_criteria_document(criteria_path)
    metrics = report["metrics"]
    operators = {
        "==": lambda actual, limit: actual == limit,
        ">=": lambda actual, limit: actual >= limit,
        "<=": lambda actual, limit: actual <= limit,
    }
    for criterion in criteria["criteria"]:
        metric = criterion["metric"]
        _require(metric in metrics, f"{criteria_path}: missing metric {metric!r}")
        operator = criterion["operator"]
        _require(
            operator in operators,
            f"{criteria_path}: unsupported operator {operator!r}",
        )
        _require(
            operators[operator](metrics[metric], criterion["limit"]),
            f"{criterion['id']}: {metric}={metrics[metric]!r} does not satisfy "
            f"{operator} {criterion['limit']!r}",
        )
    report["criteria"] = {
        "path": str(criteria_path),
        "count": len(criteria["criteria"]),
        "sha256": criteria["sha256"],
        "status": "PASS",
    }


def _manifest(path: Path) -> tuple[dict[str, str], tuple[tuple[int, Path, int], ...]]:
    values: dict[str, str] = {}
    shards: list[tuple[int, Path, int]] = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "shard":
            _require(len(fields) == 7, f"{path}:{line_number}: malformed shard row")
            rank = int(fields[1])
            relative = Path(fields[2])
            local_count = int(fields[4])
            shards.append((rank, relative, local_count))
            continue
        _require(len(fields) == 2, f"{path}:{line_number}: malformed manifest row")
        key, value = fields
        _require(
            key not in values,
            f"{path}:{line_number}: duplicate manifest key {key!r}",
        )
        values[key] = value
    _require(values.get("magic") == "AKPRST-MANIFEST", f"{path}: invalid manifest magic")
    _require(
        values.get("topology_policy") == "reject_rank_count_change",
        f"{path}: manifest does not bind DR-009 rank-count rejection",
    )
    saved_nranks = int(values["saved_nranks"])
    _require(
        [rank for rank, _, _ in shards] == list(range(saved_nranks)),
        f"{path}: shard ranks are not a contiguous saved communicator",
    )
    return values, tuple(shards)


def _checkpoint(root: Path, manifest: Path) -> Checkpoint:
    match = re.search(r"\.(\d{5})\.prst\.manifest$", manifest.name)
    _require(
        match is not None,
        f"{manifest}: missing five-digit particle checkpoint number",
    )
    values, shard_rows = _manifest(manifest)
    mesh_relative = Path(values["mesh_restart"])
    _require(not mesh_relative.is_absolute(), f"{manifest}: absolute mesh restart path")
    return Checkpoint(
        root=root,
        number=int(match.group(1)),
        manifest_relative=manifest.relative_to(root),
        mesh_relative=mesh_relative,
        witness_relative=Path(str(mesh_relative) + ".rmeta"),
        shard_relatives=tuple(Path("prst") / relative for _, relative, _ in shard_rows),
        rank_counts=tuple(local_count for _, _, local_count in shard_rows),
        topology_hash=int(values["mesh_topology_hash"]),
    )


def _discover_checkpoints(root: Path) -> tuple[Checkpoint, ...]:
    manifests = sorted((root / "prst").glob("*.prst.manifest"))
    _require(manifests, f"{root}: no typed-v2 particle restart manifests")
    checkpoints = tuple(_checkpoint(root, manifest) for manifest in manifests)
    numbers = [checkpoint.number for checkpoint in checkpoints]
    _require(len(numbers) == len(set(numbers)), f"{root}: duplicate checkpoint numbers")
    return checkpoints


def _validate_particle_states(
    particles: Sequence[ParticleState], *, expected_keys: Sequence[tuple[int, int]]
) -> dict[tuple[int, int], ParticleState]:
    keyed: dict[tuple[int, int], ParticleState] = {}
    for particle in particles:
        _require(particle.pgid >= 0, f"particle {particle.key}: negative PGID")
        _require(
            len(particle.reals) == len(REAL_FIELDS),
            f"particle {particle.key}: typed-v2 real width {len(particle.reals)} != 22",
        )
        _require(
            all(math.isfinite(value) for value in particle.reals),
            f"particle {particle.key}: non-finite typed-v2 real state",
        )
        _require(particle.key not in keyed, f"duplicate particle identity {particle.key}")
        keyed[particle.key] = particle
    _require(
        tuple(sorted(keyed)) == tuple(expected_keys),
        f"particle identities {tuple(sorted(keyed))} != expected {tuple(expected_keys)}",
    )
    return keyed


def _inspect_checkpoint(
    checkpoint: Checkpoint,
    expected_keys: Sequence[tuple[int, int]] = EXPECTED_PARTICLE_KEYS,
) -> Inspection:
    _require(
        checkpoint.mesh.is_file(),
        f"{checkpoint.manifest}: missing paired mesh restart",
    )
    _require(
        checkpoint.witness.is_file(),
        f"{checkpoint.manifest}: missing paired mesh witness",
    )
    particles: list[ParticleState] = []
    decoded_counts = []
    for expected_rank, relative in enumerate(checkpoint.shard_relatives):
        path = checkpoint.root / relative
        decoded = tracer_inspect.read_prst_file(path)
        _require(
            decoded["format_version"] == "AKPRST-v2.0",
            f"{path}: expected typed-v2, found {decoded['format_version']}",
        )
        _require(decoded["rank"] == expected_rank, f"{path}: saved rank mismatch")
        decoded_counts.append(decoded["count"])
        for particle in decoded["particles"]:
            particles.append(
                ParticleState(
                    pgid=particle["pgid"],
                    tag=particle["tag"],
                    species=particle["species"],
                    reals=tuple(particle["reals"]),
                )
            )
    _require(
        tuple(decoded_counts) == checkpoint.rank_counts,
        f"{checkpoint.manifest}: decoded shard counts differ from manifest",
    )
    keyed = _validate_particle_states(particles, expected_keys=expected_keys)
    return Inspection(checkpoint=checkpoint, particles=keyed)


def _inspect_case(
    root: Path,
    expected_keys: Sequence[tuple[int, int]] = EXPECTED_PARTICLE_KEYS,
) -> tuple[Inspection, ...]:
    return tuple(
        _inspect_checkpoint(checkpoint, expected_keys)
        for checkpoint in _discover_checkpoints(root)
    )


def _case_input(
    base_input: Path,
    case_dir: Path,
    append_text: str = "",
    *,
    source_text: Optional[str] = None,
) -> Path:
    path = case_dir / "phase8_case.athinput"
    path.write_text((base_input.read_text() if source_text is None else source_text) +
                    append_text)
    return path


def _mhd_ideal_source_text(base_input: Path) -> str:
    lines = []
    for line in base_input.read_text().splitlines():
        if line.strip().startswith(("cE0x", "cE0y", "cE0z")):
            continue
        if line.strip() == "evolution  = kinematic":
            line = "evolution  = dynamic"
        if line.strip() == "rsolver     = advect":
            line = "rsolver     = hlle"
        if line.strip() == "relativistic_field_source  = prescribed_test":
            line = "relativistic_field_source  = mhd_ideal"
        lines.append(line)
        if line.strip() == "relativistic_field_source  = mhd_ideal":
            lines.append("relativistic_temporal_sampling = frozen_tn")
    return "\n".join(lines) + "\n"


def _write_log(case_dir: Path, result: subprocess.CompletedProcess[str]) -> str:
    output = result.stdout + result.stderr
    (case_dir / "athena.log").write_text(output)
    return output


def _athena_command(
    binary: Path,
    input_path: Path,
    flags: Sequence[str],
    *,
    ranks: Optional[int],
    mpirun: str,
    pre_input: Sequence[str] = (),
) -> list[str]:
    prefix = [mpirun, "-np", str(ranks)] if ranks is not None else []
    return [*prefix, str(binary), *pre_input, "-i", str(input_path), *flags]


def _run_success(
    binary: Path,
    case_dir: Path,
    input_path: Path,
    flags: Sequence[str],
    timeout: float,
    label: str,
    *,
    ranks: Optional[int] = None,
    mpirun: str = "mpirun",
    pre_input: Sequence[str] = (),
) -> str:
    result = subprocess.run(
        _athena_command(binary, input_path, flags, ranks=ranks, mpirun=mpirun,
                        pre_input=pre_input),
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    output = _write_log(case_dir, result)
    if result.returncode:
        raise RuntimeError(
            f"{label}: Athena failed with exit code {result.returncode}:\n{output}"
        )
    return output


def _run_expected_failure(
    binary: Path,
    case_dir: Path,
    input_path: Path,
    flags: Sequence[str],
    timeout: float,
    label: str,
    expected: Sequence[str],
    *,
    ranks: Optional[int] = None,
    mpirun: str = "mpirun",
    pre_input: Sequence[str] = (),
) -> str:
    result = subprocess.run(
        _athena_command(binary, input_path, flags, ranks=ranks, mpirun=mpirun,
                        pre_input=pre_input),
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    output = _write_log(case_dir, result)
    matched = next((text for text in expected if text.lower() in output.lower()), None)
    if result.returncode == 0 or matched is None:
        raise RuntimeError(
            f"{label}: expected a nonzero exit containing one of {tuple(expected)!r}, "
            f"got exit code {result.returncode}:\n{output}"
        )
    return matched


def _fresh_case(root: Path, name: str) -> Path:
    case_dir = root / name
    _require(not case_dir.exists(), f"{case_dir}: scratch case already exists")
    case_dir.mkdir(parents=True)
    return case_dir


def _after_setup_output(output: str) -> str:
    _, marker, runtime = output.partition("Setup complete, executing task list(s)...")
    _require(marker, "Athena log lacks setup-complete marker")
    return runtime


def _mesh_diagnostics(output: str, *, after_setup: bool = False) -> tuple[int, int, int]:
    meshblocks = re.findall(r"Total number of MeshBlocks = (\d+)", output)
    _require(meshblocks, "Athena log lacks total MeshBlock diagnostic")
    churn_output = _after_setup_output(output) if after_setup else output
    churn = re.findall(
        r"(\d+) MeshBlocks created, (\d+) deleted by AMR", churn_output)
    if not churn:
        return int(meshblocks[0]), 0, 0
    return (
        int(meshblocks[0]),
        sum(int(created) for created, _ in churn),
        sum(int(deleted) for _, deleted in churn),
    )


def _require_runtime_amr_churn(output: str, label: str) -> tuple[int, int]:
    _, created, deleted = _mesh_diagnostics(output, after_setup=True)
    _require(created >= 1, f"{label}: no post-setup MeshBlocks created")
    _require(deleted >= 1, f"{label}: no post-setup MeshBlocks deleted")
    return created, deleted


def _performance_counters(
    output: str, label: str, *, after_setup: bool = False
) -> dict[str, int]:
    if after_setup:
        _, marker, output = output.partition("Setup complete, executing task list(s)...")
        _require(marker, "Athena log lacks setup-complete marker")
    counters: dict[str, int] = {}
    for match in PERFORMANCE_RE.finditer(output):
        if match.group(1) != label:
            continue
        for item in match.group(2).split():
            if "=" not in item:
                continue
            key, value = item.split("=", 1)
            counters[key] = counters.get(key, 0) + int(value)
    _require(counters, f"Athena log lacks Particle performance {label!r} counters")
    return counters


def _require_cycle_schedule(
    output: str, *, start: int, stop: int, expected_dt: float, label: str
) -> tuple[tuple[int, float, float], ...]:
    rows = tuple(
        (int(cycle), float(time), float(dt))
        for cycle, time, dt in CYCLE_RE.findall(output)
    )
    _require(rows, f"{label}: Athena log lacks cycle diagnostics")
    active = tuple(row for row in rows if start <= row[0] < stop)
    _require(
        tuple(row[0] for row in active) == tuple(range(start, stop)),
        f"{label}: cycle schedule {tuple(row[0] for row in active)} "
        f"!= {tuple(range(start, stop))}",
    )
    for cycle, _, dt in active:
        _require(dt > 1.0e-12, f"{label}: cycle {cycle} has near-zero timestep {dt}")
        _require(
            abs(dt - expected_dt) <= 2.0e-12,
            f"{label}: cycle {cycle} timestep {dt} != {expected_dt}",
        )
    return active


def _run_case(
    binary: Path,
    base_input: Path,
    root: Path,
    name: str,
    flags: Sequence[str],
    timeout: float,
    *,
    ranks: Optional[int] = None,
    mpirun: str = "mpirun",
    append_text: str = "",
    source_text: Optional[str] = None,
    pre_input: Sequence[str] = (),
    expected_keys: Sequence[tuple[int, int]] = EXPECTED_PARTICLE_KEYS,
) -> RuntimeCase:
    case_dir = _fresh_case(root, name)
    input_path = _case_input(
        base_input, case_dir, append_text, source_text=source_text)
    output = _run_success(
        binary,
        case_dir,
        input_path,
        [f"job/basename={name}", *flags],
        timeout,
        name,
        ranks=ranks,
        mpirun=mpirun,
        pre_input=pre_input,
    )
    meshblocks, created, deleted = _mesh_diagnostics(output, after_setup=True)
    inspections = _inspect_case(case_dir, expected_keys)
    return RuntimeCase(name, case_dir, output, meshblocks, created, deleted, inspections)


def _compare(left: Inspection, right: Inspection, label: str) -> float:
    _require(
        left.particles.keys() == right.particles.keys(),
        f"{label}: particle identity set changed",
    )
    max_error = 0.0
    for key, left_particle in left.particles.items():
        right_particle = right.particles[key]
        for field, left_value, right_value in zip(
            REAL_FIELDS, left_particle.reals, right_particle.reals
        ):
            error = abs(left_value - right_value)
            max_error = max(max_error, error)
            _require(
                error <= 2.0e-12 * max(1.0, abs(left_value), abs(right_value)),
                f"{label}: {key} {field} changed: {left_value} != {right_value}",
            )
    return max_error


def _prescribed_field_witness(
        initial: Inspection, final: Inspection) -> dict[str, object]:
    _require(
        initial.particles.keys() == final.particles.keys(),
        "prescribed-field witness lost particle identity",
    )
    max_field_error = 0.0
    momentum_change_count = 0
    nonzero_work_count = 0
    for key, final_particle in final.particles.items():
        initial_particle = initial.particles[key]
        for index, expected in EXPECTED_PRESCRIBED_FIELDS.items():
            max_field_error = max(
                max_field_error, abs(final_particle.reals[index] - expected))
        if any(
            abs(final_particle.reals[index] - initial_particle.reals[index]) > 1.0e-14
            for index in (14, 15, 16)
        ):
            momentum_change_count += 1
        if abs(final_particle.reals[20]) > 1.0e-14:
            nonzero_work_count += 1
    return {
        "nonzero_field_component_count": sum(
            expected != 0.0 for expected in EXPECTED_PRESCRIBED_FIELDS.values()),
        "sampled_field_max_real_error": max_field_error,
        "momentum_change_count": momentum_change_count,
        "nonzero_work_count": nonzero_work_count,
    }


def _owner_change_count(initial: Inspection, final: Inspection) -> int:
    _require(
        initial.particles.keys() == final.particles.keys(),
        "owner comparison lost identity",
    )
    return sum(
        initial.particles[key].pgid != final.particles[key].pgid
        for key in initial.particles
    )


def _hash_particle_state(value: int) -> int:
    mask = (1 << 64) - 1
    value = (value + 0x9E3779B97F4A7C15) & mask
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & mask
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & mask
    return value ^ (value >> 31)


def _tag_random_unit(seed: int, species: int, tag: int, stream: int) -> float:
    mask = (1 << 64) - 1
    value = seed ^ (species << 32)
    value ^= tag & 0xFFFFFFFF
    value ^= ((stream + 1) * 0xD2B74407B1CE6E93) & mask
    return (_hash_particle_state(value) >> 11) / 9007199254740992.0


def _wrap(value: float, lower: float, upper: float) -> float:
    length = upper - lower
    return lower + (value - lower) % length


def _periodic_wrap_count(final: Inspection, *, seed: int = 1) -> int:
    wrapped_particles = 0
    for particle in final.particles.values():
        crossed = False
        for stream, (position_index, displacement_index, lower, upper) in enumerate(
            zip(POSITION_INDICES, DISPLACEMENT_INDICES, MESH_MIN, MESH_MAX)
        ):
            initial = lower + _tag_random_unit(
                seed, particle.species, particle.tag, stream
            ) * (upper - lower)
            unwrapped = initial + particle.reals[displacement_index]
            final_position = particle.reals[position_index]
            _require(
                abs(final_position - _wrap(unwrapped, lower, upper)) <= 2.0e-12,
                f"particle {particle.key}: wrapped coordinate disagrees "
                "with displacement",
            )
            crossed = crossed or unwrapped < lower or unwrapped >= upper
        if crossed:
            wrapped_particles += 1
    return wrapped_particles


def _copy_checkpoint(checkpoint: Checkpoint, destination: Path) -> Checkpoint:
    relatives = (
        checkpoint.manifest_relative,
        checkpoint.mesh_relative,
        checkpoint.witness_relative,
        *checkpoint.shard_relatives,
    )
    for relative in relatives:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(checkpoint.root / relative, target)
    return _checkpoint(destination, destination / checkpoint.manifest_relative)


def _restart_flags(checkpoint: Checkpoint) -> list[str]:
    return [
        "problem/prtcl_rst_flag=1",
        f"problem/prtcl_res_file={checkpoint.rank_zero_shard_relative}",
    ]


def _run_continuation(
    mpi_binary: Path,
    input_path: Path,
    root: Path,
    name: str,
    checkpoint: Checkpoint,
    basename: str,
    flags: Sequence[str],
    timeout: float,
    mpirun: str,
    *,
    append_text: str = "",
    expected_keys: Sequence[tuple[int, int]] = EXPECTED_PARTICLE_KEYS,
) -> RuntimeCase:
    continuation_dir = _fresh_case(root, name)
    copied = _copy_checkpoint(checkpoint, continuation_dir)
    continuation_input = _case_input(input_path, continuation_dir, append_text)
    output = _run_success(
        mpi_binary,
        continuation_dir,
        continuation_input,
        [
            f"job/basename={basename}",
            *flags,
            *_restart_flags(copied),
        ],
        timeout,
        name,
        ranks=len(copied.rank_counts),
        mpirun=mpirun,
        pre_input=("-r", str(copied.mesh_relative)),
    )
    meshblocks, created, deleted = _mesh_diagnostics(output, after_setup=True)
    return RuntimeCase(
        name,
        continuation_dir,
        output,
        meshblocks,
        created,
        deleted,
        _inspect_case(continuation_dir, expected_keys),
    )


def _retained_runtime_negative_controls(
    binary: Path,
    mpi_binary: Path,
    input_path: Path,
    root: Path,
    timeout: float,
    mpirun: str,
) -> dict[str, str]:
    controls: dict[str, str] = {}
    nonperiodic_dir = _fresh_case(root, "prescribed_nonperiodic_rejection")
    nonperiodic_input = _case_input(input_path, nonperiodic_dir)
    controls["prescribed_nonperiodic"] = _run_expected_failure(
        binary,
        nonperiodic_dir,
        nonperiodic_input,
        [
            "job/basename=prescribed_nonperiodic_rejection",
            "mesh/ix1_bc=outflow",
            "mesh/ox1_bc=outflow",
        ],
        timeout,
        "prescribed_nonperiodic_rejection",
        (STRICTLY_PERIODIC_REJECTION,),
    )

    mhd_ideal_text = _mhd_ideal_source_text(input_path)
    for name, ranks, append_text, flags in (
        ("mhd_ideal_mpi_rejection", 2, "", MULTIBLOCK_FLAGS),
        ("mhd_ideal_multiblock_rejection", None, "", MULTIBLOCK_FLAGS),
        ("mhd_ideal_smr_rejection", None, SMR_APPEND, MULTIBLOCK_FLAGS),
        (
            "mhd_ideal_amr_rejection",
            None,
            AMR_APPEND,
            (*MULTIBLOCK_FLAGS, "problem/particle_refinement=moving_boxes"),
        ),
    ):
        case_dir = _fresh_case(root, name)
        case_input = _case_input(
            input_path, case_dir, append_text, source_text=mhd_ideal_text)
        controls[name] = _run_expected_failure(
            mpi_binary if ranks is not None else binary,
            case_dir,
            case_input,
            [f"job/basename={name}", *flags],
            timeout,
            name,
            (
                RUNTIME_FENCE,
                MHD_IDEAL_RUNTIME_FENCE,
                MHD_IDEAL_ONE_MESHBLOCK_REJECTION,
            ),
            ranks=ranks,
            mpirun=mpirun,
        )
    return controls


@contextlib.contextmanager
def _preserved_runtime_root(work_dir: Optional[Path], prefix: str) -> Iterator[Path]:
    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix=prefix) as temporary:
            yield Path(temporary)
        return
    root = work_dir.resolve()
    if root.exists():
        _require(not any(root.iterdir()), f"{root}: --work-dir must be empty")
    else:
        root.mkdir(parents=True)
    yield root


def _case_summary(case: RuntimeCase) -> dict:
    return {
        "name": case.name,
        "root": case.name,
        "meshblocks": case.meshblocks,
        "created_blocks": case.created_blocks,
        "deleted_blocks": case.deleted_blocks,
        "checkpoint_numbers": [
            inspection.checkpoint.number for inspection in case.inspections
        ],
        "final_rank_counts": list(case.final.checkpoint.rank_counts),
        "final_topology_hash": case.final.checkpoint.topology_hash,
    }


def _qualify(
    binary: Path,
    mpi_binary: Path,
    input_path: Path,
    criteria_path: Path,
    root: Path,
    timeout: float,
    mpirun: str,
) -> dict:
    _validate_criteria_document(criteria_path)
    cases: list[RuntimeCase] = []

    serial_single = _run_case(binary, input_path, root, "serial_single", (), timeout)
    cases.append(serial_single)
    serial_multi = _run_case(
        binary, input_path, root, "serial_multiblock", MULTIBLOCK_FLAGS, timeout)
    cases.append(serial_multi)
    serial_error = _compare(
        serial_single.final, serial_multi.final, "serial one/multiblock"
    )
    owner_changes = _owner_change_count(serial_multi.initial, serial_multi.final)
    _require(owner_changes >= 1, "serial multiblock case did not migrate an owner")
    wrap_count = _periodic_wrap_count(serial_multi.final)
    _require(wrap_count >= 1, "serial multiblock case did not exercise periodic wrap")
    field_witness = _prescribed_field_witness(serial_multi.initial, serial_multi.final)
    _require(
        field_witness["momentum_change_count"] >= 1,
        "nonzero prescribed fields did not change particle momentum",
    )
    _require(
        field_witness["nonzero_work_count"] >= 1,
        "nonzero prescribed fields did not produce particle work",
    )

    mpi_cases = []
    mpi_error = 0.0
    for ranks in EXPECTED_RANK_LADDER:
        case = _run_case(
            mpi_binary,
            input_path,
            root,
            f"mpi_{ranks}",
            MULTIBLOCK_FLAGS,
            timeout,
            ranks=ranks,
            mpirun=mpirun,
        )
        cases.append(case)
        mpi_cases.append(case)
        mpi_error = max(
            mpi_error,
            _compare(serial_multi.final, case.final, f"serial/MPI-{ranks}"),
        )
    empty_rank_count = sum(
        count == 0 for count in mpi_cases[-1].final.checkpoint.rank_counts
    )
    _require(empty_rank_count >= 1, "MPI-8 case lacks a particle-empty rank")

    fine_uniform = _run_case(
        binary,
        input_path,
        root,
        "serial_fine_uniform_reference",
        [*FINE_UNIFORM_FLAGS[:-1], f"particles/ppc={STATIC_FINE_UNIFORM_PPC}"],
        timeout,
        expected_keys=STATIC_SMR_PARTICLE_KEYS,
    )
    cases.append(fine_uniform)
    static_smr = _run_case(
        binary,
        input_path,
        root,
        "serial_static_smr",
        [*MULTIBLOCK_FLAGS, f"particles/ppc={STATIC_SMR_PPC}"],
        timeout,
        append_text=SMR_APPEND,
        expected_keys=STATIC_SMR_PARTICLE_KEYS,
    )
    cases.append(static_smr)
    smr_error = _compare(fine_uniform.final, static_smr.final, "fine-uniform/static-SMR")
    _require(
        static_smr.meshblocks > serial_multi.meshblocks,
        "static SMR case did not create a refined hierarchy",
    )
    mpi_static_smr = _run_case(
        mpi_binary,
        input_path,
        root,
        "mpi_4_static_smr",
        [*MULTIBLOCK_FLAGS, f"particles/ppc={STATIC_SMR_PPC}"],
        timeout,
        ranks=4,
        mpirun=mpirun,
        append_text=SMR_APPEND,
        expected_keys=STATIC_SMR_PARTICLE_KEYS,
    )
    cases.append(mpi_static_smr)
    smr_mpi_error = _compare(
        static_smr.final, mpi_static_smr.final, "serial/MPI-4 static-SMR")
    smr_mpi_exchange = _performance_counters(
        mpi_static_smr.log, "particle_mpi_exchange", after_setup=True)
    smr_level_lookup = _performance_counters(
        mpi_static_smr.log, "particle_multilevel_lookup", after_setup=True)
    _require(
        smr_mpi_exchange.get("sent", 0) >= 1,
        "MPI static-SMR case did not send a particle through multilevel lookup",
    )
    _require(
        smr_level_lookup.get("diagnostics", 0) >= 1,
        "MPI static-SMR case did not cross a refinement level",
    )
    mpi_static_smr_host_tree = _run_case(
        mpi_binary,
        input_path,
        root,
        "mpi_4_static_smr_host_tree",
        [
            *MULTIBLOCK_FLAGS,
            f"particles/ppc={STATIC_SMR_PPC}",
            "particles/amr_remap=host_tree",
        ],
        timeout,
        ranks=4,
        mpirun=mpirun,
        append_text=SMR_APPEND,
        expected_keys=STATIC_SMR_PARTICLE_KEYS,
    )
    cases.append(mpi_static_smr_host_tree)
    smr_host_tree_mpi_error = _compare(
        static_smr.final,
        mpi_static_smr_host_tree.final,
        "serial/MPI-4 static-SMR host-tree",
    )
    smr_host_tree_mpi_exchange = _performance_counters(
        mpi_static_smr_host_tree.log, "particle_mpi_exchange", after_setup=True)
    smr_host_tree_level_lookup = _performance_counters(
        mpi_static_smr_host_tree.log, "particle_multilevel_lookup", after_setup=True)
    _require(
        smr_host_tree_mpi_exchange.get("sent", 0) >= 1,
        "MPI static-SMR host-tree case did not send a particle through multilevel lookup",
    )
    _require(
        smr_host_tree_level_lookup.get("diagnostics", 0) >= 1,
        "MPI static-SMR host-tree case did not cross a refinement level",
    )

    adaptive_flags = (
        *MULTIBLOCK_FLAGS,
        *ADAPTIVE_TIMESTEP_FLAGS,
        "problem/particle_refinement=moving_boxes",
        "time/nlim=3",
        "time/tlim=0.1",
    )
    adaptive_uniform = _run_case(
        mpi_binary,
        input_path,
        root,
        "mpi_4_adaptive_uniform_reference",
        [
            *FINE_UNIFORM_FLAGS,
            *ADAPTIVE_TIMESTEP_FLAGS,
            "time/nlim=3",
            "time/tlim=0.1",
        ],
        timeout,
        ranks=4,
        mpirun=mpirun,
    )
    cases.append(adaptive_uniform)
    adaptive_cases = []
    adaptive_uniform_errors = []
    adaptive_remap_sent = []
    for mode in ("device_table", "host_tree"):
        case = _run_case(
            mpi_binary,
            input_path,
            root,
            f"mpi_4_adaptive_{mode}",
            [
                *adaptive_flags,
                f"particles/amr_remap={mode}",
            ],
            timeout,
            ranks=4,
            mpirun=mpirun,
            append_text=AMR_APPEND,
        )
        _require_runtime_amr_churn(case.log, f"adaptive {mode}")
        remap_counters = _performance_counters(
            case.log, "amr_remap", after_setup=True)
        _require(
            remap_counters.get("sent", 0) >= 1,
            f"adaptive {mode}: no off-rank RemapAfterAMR send",
        )
        cases.append(case)
        adaptive_cases.append(case)
        adaptive_remap_sent.append(remap_counters.get("sent", 0))
        adaptive_uniform_errors.append(
            _compare(
                adaptive_uniform.final,
                case.final,
                f"fine-uniform/adaptive {mode}",
            )
        )
    adaptive_error = _compare(
        adaptive_cases[0].final, adaptive_cases[1].final, "adaptive remap modes")
    adaptive_restart_errors = []
    adaptive_restart_cases = []
    adaptive_restart_references = []
    for mode, source in zip(("device_table", "host_tree"), adaptive_cases):
        flags = [
            *MULTIBLOCK_FLAGS,
            *ADAPTIVE_TIMESTEP_FLAGS,
            "problem/particle_refinement=moving_boxes",
            "time/nlim=5",
            "time/tlim=1.0",
            f"particles/amr_remap={mode}",
        ]
        long_reference = _run_case(
            mpi_binary,
            input_path,
            root,
            f"mpi_4_adaptive_{mode}_long_reference",
            flags,
            timeout,
            ranks=4,
            mpirun=mpirun,
            append_text=AMR_APPEND,
        )
        continuation = _run_continuation(
            mpi_binary,
            input_path,
            root,
            f"mpi_4_adaptive_{mode}_continuation_after_amr",
            source.final.checkpoint,
            source.name,
            flags,
            timeout,
            mpirun,
            append_text=AMR_APPEND,
        )
        cases.extend((long_reference, continuation))
        adaptive_restart_references.append(long_reference)
        adaptive_restart_cases.append(continuation)
        adaptive_restart_errors.append(
            _compare(
                long_reference.final,
                continuation.final,
                f"adaptive {mode} continuation after AMR",
            )
        )

    mpi4_source = _run_case(
        mpi_binary,
        input_path,
        root,
        "mpi_4_continuation_source",
        [
            *MULTIBLOCK_FLAGS,
            *CONTINUATION_TIMESTEP_FLAGS,
            "time/nlim=6",
            "time/tlim=1.0",
        ],
        timeout,
        ranks=4,
        mpirun=mpirun,
    )
    cases.append(mpi4_source)
    mpi4_long = _run_case(
        mpi_binary,
        input_path,
        root,
        "mpi_4_long_reference",
        [
            *MULTIBLOCK_FLAGS,
            *CONTINUATION_TIMESTEP_FLAGS,
            "time/nlim=8",
            "time/tlim=0.8",
        ],
        timeout,
        ranks=4,
        mpirun=mpirun,
    )
    cases.append(mpi4_long)
    continuation_before = _run_continuation(
        mpi_binary,
        input_path,
        root,
        "mpi_4_continuation_before_migration",
        mpi4_source.initial.checkpoint,
        "mpi_4_continuation_source",
        [
            *MULTIBLOCK_FLAGS,
            *CONTINUATION_TIMESTEP_FLAGS,
            "time/nlim=8",
            "time/tlim=0.8",
        ],
        timeout,
        mpirun=mpirun,
    )
    cases.append(continuation_before)
    continuation_after = _run_continuation(
        mpi_binary,
        input_path,
        root,
        "mpi_4_continuation_after_migration",
        mpi4_source.final.checkpoint,
        "mpi_4_continuation_source",
        [
            *MULTIBLOCK_FLAGS,
            *CONTINUATION_TIMESTEP_FLAGS,
            "time/nlim=8",
            "time/tlim=0.8",
        ],
        timeout,
        mpirun,
    )
    cases.append(continuation_after)
    continuation_before_error = _compare(
        mpi4_long.final, continuation_before.final, "continuation before migration")
    continuation_after_error = _compare(
        mpi4_long.final, continuation_after.final, "continuation after migration")
    continuation_source_schedule = _require_cycle_schedule(
        mpi4_source.log, start=0, stop=6, expected_dt=0.1,
        label="continuation source")
    continuation_long_schedule = _require_cycle_schedule(
        mpi4_long.log, start=0, stop=8, expected_dt=0.1,
        label="continuation uninterrupted reference")
    continuation_before_schedule = _require_cycle_schedule(
        continuation_before.log, start=0, stop=8, expected_dt=0.1,
        label="continuation before migration")
    continuation_after_schedule = _require_cycle_schedule(
        continuation_after.log, start=6, stop=8, expected_dt=0.1,
        label="continuation after migration")
    mpi4_owner_changes = _owner_change_count(mpi4_source.initial, mpi4_source.final)
    _require(mpi4_owner_changes >= 1, "MPI-4 continuation source did not migrate")

    rank_change_dir = _fresh_case(root, "mpi_rank_count_change_rejection")
    rank_change_checkpoint = _copy_checkpoint(
        mpi4_source.final.checkpoint, rank_change_dir
    )
    rank_change_input = _case_input(input_path, rank_change_dir)
    rank_change_match = _run_expected_failure(
        mpi_binary,
        rank_change_dir,
        rank_change_input,
        [
            "job/basename=mpi_4_continuation_source",
            *MULTIBLOCK_FLAGS,
            *CONTINUATION_TIMESTEP_FLAGS,
            *_restart_flags(rank_change_checkpoint),
        ],
        timeout,
        "mpi_rank_count_change_rejection",
        (EXACT_RANK_COUNT_REJECTION,),
        ranks=2,
        mpirun=mpirun,
        pre_input=("-r", str(rank_change_checkpoint.mesh_relative)),
    )
    retained_rejections = _retained_runtime_negative_controls(
        binary, mpi_binary, input_path, root, timeout, mpirun)

    inspected = [inspection for case in cases for inspection in case.inspections]
    metrics = {
        "typed_v2_real_field_count": len(REAL_FIELDS),
        "prescribed_nonzero_field_component_count":
            field_witness["nonzero_field_component_count"],
        "prescribed_sampled_field_max_real_error":
            field_witness["sampled_field_max_real_error"],
        "prescribed_momentum_change_count": field_witness["momentum_change_count"],
        "prescribed_nonzero_work_count": field_witness["nonzero_work_count"],
        "serial_single_multiblock_identity_equal":
            serial_single.final.particles.keys() == serial_multi.final.particles.keys(),
        "serial_single_multiblock_max_real_error": serial_error,
        "serial_multiblock_owner_change_count": owner_changes,
        "periodic_wrap_particle_count": wrap_count,
        "mpi_rank_ladder_exact": tuple(
            len(case.final.checkpoint.rank_counts) for case in mpi_cases
        ) == EXPECTED_RANK_LADDER,
        "mpi_rank_ladder_max_real_error": mpi_error,
        "mpi_rank8_empty_particle_rank_count": empty_rank_count,
        "mpi_identity_preserved": all(
            case.final.particles.keys() == serial_multi.final.particles.keys()
            for case in mpi_cases
        ),
        "static_smr_identity_equal":
            static_smr.final.particles.keys() == fine_uniform.final.particles.keys(),
        "static_smr_uniform_field_max_real_error": smr_error,
        "static_smr_has_more_meshblocks_than_uniform":
            static_smr.meshblocks > serial_multi.meshblocks,
        "static_smr_mpi_identity_equal":
            static_smr.final.particles.keys() == mpi_static_smr.final.particles.keys(),
        "static_smr_mpi_max_real_error": smr_mpi_error,
        "static_smr_mpi_exchange_sent": smr_mpi_exchange.get("sent", 0),
        "static_smr_cross_level_transition_count":
            smr_level_lookup.get("diagnostics", 0),
        "static_smr_mpi_multilevel_lookup_exercised":
            smr_level_lookup.get("diagnostics", 0) >= 1,
        "static_smr_host_tree_mpi_identity_equal":
            static_smr.final.particles.keys() ==
            mpi_static_smr_host_tree.final.particles.keys(),
        "static_smr_host_tree_mpi_max_real_error": smr_host_tree_mpi_error,
        "static_smr_host_tree_mpi_exchange_sent":
            smr_host_tree_mpi_exchange.get("sent", 0),
        "static_smr_host_tree_cross_level_transition_count":
            smr_host_tree_level_lookup.get("diagnostics", 0),
        "static_smr_host_tree_mpi_multilevel_lookup_exercised":
            smr_host_tree_level_lookup.get("diagnostics", 0) >= 1,
        "adaptive_device_created_blocks": adaptive_cases[0].created_blocks,
        "adaptive_device_deleted_blocks": adaptive_cases[0].deleted_blocks,
        "adaptive_host_created_blocks": adaptive_cases[1].created_blocks,
        "adaptive_host_deleted_blocks": adaptive_cases[1].deleted_blocks,
        "adaptive_identity_preserved": all(
            case.initial.particles.keys() == case.final.particles.keys()
            for case in adaptive_cases
        ),
        "adaptive_remap_mode_max_real_error": adaptive_error,
        "adaptive_device_off_rank_remap_sent": adaptive_remap_sent[0],
        "adaptive_host_off_rank_remap_sent": adaptive_remap_sent[1],
        "adaptive_off_rank_remap_modes_exercised": sum(
            sent >= 1 for sent in adaptive_remap_sent),
        "adaptive_uniform_reference_identity_equal": all(
            case.final.particles.keys() == adaptive_uniform.final.particles.keys()
            for case in adaptive_cases
        ),
        "adaptive_device_uniform_reference_max_real_error":
            adaptive_uniform_errors[0],
        "adaptive_host_uniform_reference_max_real_error":
            adaptive_uniform_errors[1],
        "adaptive_uniform_reference_comparisons": len(adaptive_uniform_errors),
        "adaptive_restart_after_amr_comparisons": len(adaptive_restart_errors),
        "adaptive_restart_after_amr_identity_equal": all(
            reference.final.particles.keys() == continuation.final.particles.keys()
            for reference, continuation in zip(
                adaptive_restart_references, adaptive_restart_cases)
        ),
        "adaptive_restart_after_amr_max_real_error":
            max(adaptive_restart_errors),
        "all_checkpoint_reals_finite": True,
        "all_checkpoint_identity_exact": True,
        "inspected_checkpoint_count": len(inspected),
        "continuation_before_migration_identity_equal":
            mpi4_long.final.particles.keys() ==
            continuation_before.final.particles.keys(),
        "continuation_before_migration_max_real_error": continuation_before_error,
        "continuation_after_migration_identity_equal":
            mpi4_long.final.particles.keys() ==
            continuation_after.final.particles.keys(),
        "continuation_after_migration_max_real_error": continuation_after_error,
        "continuation_same_rank_count":
            all(
                len(checkpoint.rank_counts) == 4
                for checkpoint in (
                    mpi4_source.initial.checkpoint,
                    mpi4_source.final.checkpoint,
                    continuation_before.final.checkpoint,
                    continuation_after.final.checkpoint,
                )
            ),
        "continuation_same_topology_hash":
            len({
                mpi4_source.initial.checkpoint.topology_hash,
                mpi4_source.final.checkpoint.topology_hash,
                continuation_before.final.checkpoint.topology_hash,
                continuation_after.final.checkpoint.topology_hash,
            }) == 1,
        "continuation_before_source_owner_change_count": mpi4_owner_changes,
        "continuation_source_schedule_length": len(continuation_source_schedule),
        "continuation_long_schedule_length": len(continuation_long_schedule),
        "continuation_before_schedule_length": len(continuation_before_schedule),
        "continuation_after_schedule_length": len(continuation_after_schedule),
        "rank_count_change_restart_rejected": True,
        "rank_count_change_restart_diagnostic_exact":
            rank_change_match == EXACT_RANK_COUNT_REJECTION,
        "prescribed_nonperiodic_rejected":
            retained_rejections["prescribed_nonperiodic"] ==
            STRICTLY_PERIODIC_REJECTION,
        "mhd_ideal_mpi_rejected": "mhd_ideal_mpi_rejection" in retained_rejections,
        "mhd_ideal_multiblock_rejected":
            "mhd_ideal_multiblock_rejection" in retained_rejections,
        "mhd_ideal_smr_rejected": "mhd_ideal_smr_rejection" in retained_rejections,
        "mhd_ideal_amr_rejected": "mhd_ideal_amr_rejection" in retained_rejections,
        "retained_runtime_negative_control_count": len(retained_rejections),
        "temporary_case_directories_isolated": True,
        "adaptive_remap_modes_exercised": len(adaptive_cases),
        "expected_particle_identity_count": len(serial_multi.final.particles),
        "continuation_phase_count": 2,
        "static_smr_mpi_rank_count": len(mpi_static_smr.final.checkpoint.rank_counts),
        "adaptive_mpi_rank_count": len(adaptive_cases[0].final.checkpoint.rank_counts),
    }
    report = {
        "schema_version": 1,
        "qualification": "cr_relativistic_migration_runtime_phase8",
        "status": "PASS",
        "binary": str(binary),
        "mpi_binary": str(mpi_binary),
        "input": str(input_path),
        "provenance": {
            "source_head": _source_head(),
            "binary_sha256": _sha256(binary),
            "mpi_binary_sha256": _sha256(mpi_binary),
            "input_sha256": _sha256(input_path),
            "criteria_sha256": _sha256(criteria_path),
            "script_sha256": _sha256(Path(__file__).resolve()),
        },
        "cases": [_case_summary(case) for case in cases],
        "rank_count_change_rejection": {
            "expected_diagnostic": EXACT_RANK_COUNT_REJECTION,
            "observed_diagnostic": rank_change_match,
        },
        "retained_runtime_negative_controls": retained_rejections,
        "metrics": metrics,
    }
    _validate_runtime_criteria(criteria_path, report)
    return report


def _pre_fence_smoke(
    binary: Path,
    mpi_binary: Path,
    input_path: Path,
    criteria_path: Path,
    root: Path,
    timeout: float,
    mpirun: str,
) -> dict:
    criteria = _validate_criteria_document(criteria_path)
    serial_single = _run_case(binary, input_path, root, "serial_single", (), timeout)
    serial_multi = _run_case(
        binary, input_path, root, "serial_multiblock", MULTIBLOCK_FLAGS, timeout)
    serial_error = _compare(
        serial_single.final, serial_multi.final, "serial one/multiblock"
    )
    owner_changes = _owner_change_count(serial_multi.initial, serial_multi.final)
    wrap_count = _periodic_wrap_count(serial_multi.final)
    field_witness = _prescribed_field_witness(serial_multi.initial, serial_multi.final)

    fenced: dict[str, str] = {}
    for name, ranks, append_text, flags in (
        ("mpi_2_fence", 2, "", MULTIBLOCK_FLAGS),
        ("serial_static_smr_fence", None, SMR_APPEND,
         (*MULTIBLOCK_FLAGS, f"particles/ppc={STATIC_SMR_PPC}")),
        ("mpi_4_static_smr_fence", 4, SMR_APPEND,
         (*MULTIBLOCK_FLAGS, f"particles/ppc={STATIC_SMR_PPC}")),
        ("serial_adaptive_fence", None, AMR_APPEND,
         (*MULTIBLOCK_FLAGS, "problem/particle_refinement=moving_boxes")),
        ("mpi_4_adaptive_device_fence", 4, AMR_APPEND,
         (*MULTIBLOCK_FLAGS, "problem/particle_refinement=moving_boxes",
          "particles/amr_remap=device_table")),
        ("mpi_4_adaptive_host_fence", 4, AMR_APPEND,
         (*MULTIBLOCK_FLAGS, "problem/particle_refinement=moving_boxes",
          "particles/amr_remap=host_tree")),
    ):
        case_dir = _fresh_case(root, name)
        case_input = _case_input(input_path, case_dir, append_text)
        fenced[name] = _run_expected_failure(
            mpi_binary if ranks is not None else binary,
            case_dir,
            case_input,
            [f"job/basename={name}", *flags],
            timeout,
            name,
            (RUNTIME_FENCE,),
            ranks=ranks,
            mpirun=mpirun,
        )

    rank_change_dir = _fresh_case(root, "rank_count_change_pre_fence")
    copied = _copy_checkpoint(serial_multi.final.checkpoint, rank_change_dir)
    rank_change_input = _case_input(input_path, rank_change_dir)
    observed_rank_change = _run_expected_failure(
        mpi_binary,
        rank_change_dir,
        rank_change_input,
        [
            "job/basename=serial_multiblock",
            *MULTIBLOCK_FLAGS,
            *_restart_flags(copied),
        ],
        timeout,
        "rank_count_change_pre_fence",
        (EXACT_RANK_COUNT_REJECTION, TOPOLOGY_REJECTION),
        ranks=2,
        mpirun=mpirun,
        pre_input=("-r", str(copied.mesh_relative)),
    )
    retained_rejections = _retained_runtime_negative_controls(
        binary, mpi_binary, input_path, root, timeout, mpirun)
    gaps = [
        "Positive-cycle MPI, static SMR, and adaptive AMR remain deliberately fenced "
        "until the parent implementation removes the Phase-8 runtime fence.",
    ]
    if observed_rank_change != EXACT_RANK_COUNT_REJECTION:
        gaps.append(
            "The current C++ reader rejects a changed-rank restart deterministically, "
            "but the changed topology diagnostic preempts the exact DR-009 rank-count "
            "diagnostic required by the post-fence qualifier."
        )
    return {
        "schema_version": 1,
        "qualification": "cr_relativistic_migration_runtime_phase8_pre_fence_smoke",
        "status": "PASS",
        "binary": str(binary),
        "mpi_binary": str(mpi_binary),
        "input": str(input_path),
        "criteria": {
            "path": str(criteria_path),
            "count": len(criteria["criteria"]),
            "sha256": criteria["sha256"],
            "status": "PASS",
        },
        "serial_subset": {
            "single": _case_summary(serial_single),
            "multiblock": _case_summary(serial_multi),
            "max_real_error": serial_error,
            "owner_change_count": owner_changes,
            "periodic_wrap_particle_count": wrap_count,
            "prescribed_field_witness": field_witness,
            "all_22_typed_v2_reals_inspected": True,
        },
        "expected_fence_rejections": fenced,
        "retained_runtime_negative_controls": retained_rejections,
        "rank_count_change_pre_fence": {
            "accepted_diagnostics": [EXACT_RANK_COUNT_REJECTION, TOPOLOGY_REJECTION],
            "observed_diagnostic": observed_rank_change,
            "exact_dr009_diagnostic": observed_rank_change == EXACT_RANK_COUNT_REJECTION,
        },
        "gaps": gaps,
    }


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
        description=(
            "Qualify bounded prescribed-test relativistic CR migration after parent "
            "fence removal, or run the executable pre-fence subset."
        )
    )
    parser.add_argument("--binary", type=Path, help="serial Athena binary")
    parser.add_argument("--mpi-binary", type=Path, help="MPI-enabled Athena binary")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--timeout", type=float, default=90.0)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--json", metavar="PATH")
    parser.add_argument(
        "--check-registration",
        action="store_true",
        help="validate only the frozen criteria digest, schema, and exact ID sequence",
    )
    parser.add_argument(
        "--pre-fence-smoke",
        action="store_true",
        help="run the executable subset and require the current MPI/SMR/AMR fences",
    )
    args = parser.parse_args()

    input_path = args.input.resolve()
    criteria_path = args.criteria.resolve()
    report = {
        "schema_version": 1,
        "qualification": "cr_relativistic_migration_runtime_phase8",
        "status": "FAIL",
        "input": str(input_path),
        "criteria": str(criteria_path),
    }
    try:
        _require(input_path.is_file(), f"{input_path}: --input is not a file")
        _require(criteria_path.is_file(), f"{criteria_path}: --criteria is not a file")
        _require(args.timeout > 0.0, "--timeout must be positive")
        if args.check_registration:
            criteria = _validate_criteria_document(criteria_path)
            report = {
                "schema_version": 1,
                "qualification": "cr_relativistic_migration_runtime_phase8_registration",
                "status": "PASS",
                "criteria": {
                    "path": str(criteria_path),
                    "count": len(criteria["criteria"]),
                    "sha256": criteria["sha256"],
                },
            }
        else:
            _require(
                args.binary is not None,
                "--binary is required for runtime execution",
            )
            _require(
                args.mpi_binary is not None,
                "--mpi-binary is required for runtime execution",
            )
            binary = args.binary.resolve()
            mpi_binary = args.mpi_binary.resolve()
            _require(binary.is_file(), f"{binary}: --binary is not a file")
            _require(mpi_binary.is_file(), f"{mpi_binary}: --mpi-binary is not a file")
            with _preserved_runtime_root(
                args.work_dir,
                "cr-rel-phase8-pre-fence-" if args.pre_fence_smoke else "cr-rel-phase8-",
            ) as root:
                if args.pre_fence_smoke:
                    report = _pre_fence_smoke(
                        binary, mpi_binary, input_path, criteria_path, root,
                        args.timeout, args.mpirun)
                else:
                    report = _qualify(
                        binary, mpi_binary, input_path, criteria_path, root,
                        args.timeout, args.mpirun)
    # Runtime directories retain full Athena logs when requested.
    except Exception as error:
        report["error"] = str(error)
        _write_json(args.json, report)
        print(
            f"CR relativistic Phase-8 migration qualification FAIL: {error}",
            file=sys.stderr,
        )
        return 1

    _write_json(args.json, report)
    print(f"{report['qualification']} PASS")
    if "criteria" in report and isinstance(report["criteria"], dict):
        print(f"  frozen criteria validated: {report['criteria']['count']}")
    if args.pre_fence_smoke:
        subset = report["serial_subset"]
        print(
            "  executable serial subset: one-block/multiblock typed-v2 parity, "
            f"owner changes={subset['owner_change_count']}, "
            f"periodic wraps={subset['periodic_wrap_particle_count']}"
        )
        print(
            "  expected runtime fences rejected: "
            f"{len(report['expected_fence_rejections'])}"
        )
        for gap in report["gaps"]:
            print(f"  gap: {gap}")
    elif not args.check_registration:
        print(f"  runtime cases inspected: {len(report['cases'])}")
        print(f"  preregistered criteria passed: {report['criteria']['count']}")
    if args.json:
        print(f"  JSON report: {args.json}")
    if args.work_dir and not args.check_registration:
        print(f"  preserved runtime cases: {args.work_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
