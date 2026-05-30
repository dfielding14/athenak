#!/usr/bin/env python3
"""Qualify serial Phase-7 relativistic CR diagnostic meanings and restart persistence."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
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
DEFAULT_INPUT = (
    REPO_ROOT / "inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput"
)
DEFAULT_CRITERIA = (
    REPO_ROOT
    / "docs/source/modules/figures/cr_tracer_relativistic_acceleration"
    / "cr_relativistic_phase7_preregistered_criteria.json"
)
SAMPLE_FIELDS = (
    "vx", "vy", "vz", "wx", "wy", "wz", "wmag", "gamma",
    "kinetic_energy_model", "cex", "cey", "cez", "work", "alpha_s",
    "r_larmor_over_dx_min",
)
PSAMP_METADATA = {
    "selected_species": "-1",
    "sample_stride": "1",
    "sample_offset": "0",
    "sample_count": "1",
    "field_count": str(len(SAMPLE_FIELDS)),
    "schema": "akcr_particle_output_v1",
    "mode": "relativistic_hc",
    "units": "code_model",
    "c_model": "1.00000000000000000e+00",
    "alpha_s": "1.00000000000000000e+00",
    "position_units": "code_length",
    "velocity_units": "code_velocity",
    "w_units": "code_velocity",
    "cE_units": "code_velocity_times_code_B",
    "gamma_units": "dimensionless",
    "kinetic_energy_model_units": "code_velocity_squared",
    "work_units": "code_velocity_squared",
    "alpha_s_units": "inverse_code_time_per_code_B",
    "r_larmor_over_dx_min_units": "dimensionless",
}
PSPEC_METADATA = {
    "quantity": "wmag",
    "nspecies": "1",
    "nbin": "8",
    "reduce": "1",
    "schema": "akcr_particle_output_v1",
    "mode": "relativistic_hc",
    "units": "code_model",
    "c_model": "1.00000000000000000e+00",
    "alpha_s": "1.00000000000000000e+00",
    "quantity_units": "code_velocity",
    "histogram_units": "particle_count",
}
EXPECTED_CRITERIA_SHA256 = (
    "1608e681952fa5d0133d8fc2ae1b0681fa17a09a667272340446625dccfba0d1"
)
EXPECTED_CRITERIA_IDS = tuple(f"P7-{index:02d}" for index in range(1, 44))
sys.path.insert(0, str(SCRIPT_DIR))
import cr_tracer_inspect as tracer_inspect  # noqa: E402


@dataclass(frozen=True)
class Checkpoint:
    """Paths for one serial paired mesh and particle restart checkpoint."""

    root: Path
    shard_relative: Path
    manifest_relative: Path
    mesh_relative: Path
    witness_relative: Path


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _near(actual: float, expected: float, tolerance: float, label: str) -> float:
    error = abs(actual - expected)
    scale = max(1.0, abs(expected))
    _require(
        error <= tolerance * scale,
        f"{label}: error={error:.6e}, tolerance={tolerance:.6e}",
    )
    return error


def _vector_error(actual: Sequence[float], expected: Sequence[float],
                  tolerance: float, label: str) -> float:
    error = math.sqrt(sum((a - b) ** 2 for a, b in zip(actual, expected)))
    scale = max(1.0, math.sqrt(sum(value * value for value in expected)))
    _require(
        error <= tolerance * scale,
        f"{label}: error={error:.6e}, tolerance={tolerance:.6e}",
    )
    return error


def _flags(**values: object) -> list[str]:
    return [f"{key.replace('__', '/')}={value}" for key, value in values.items()]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_success(binary: Path, case_dir: Path, arguments: Sequence[str],
                 timeout: float, label: str) -> None:
    case_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(binary), *arguments],
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    output = result.stdout + result.stderr
    (case_dir / "athena.log").write_text(output)
    if result.returncode:
        raise RuntimeError(
            f"{label}: Athena failed with exit code {result.returncode}:\n{output}")


def _manifest_values(path: Path) -> dict[str, str]:
    values = {}
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields or fields[0] == "shard":
            continue
        if len(fields) != 2 or fields[0] in values:
            raise RuntimeError(f"{path}:{line_number}: malformed manifest row")
        values[fields[0]] = fields[1]
    return values


def _checkpoint(root: Path, shard: Path) -> Checkpoint:
    manifest_relative = Path("prst") / f"{shard.name}.manifest"
    manifest = root / manifest_relative
    _require(manifest.is_file(), f"{shard}: missing manifest {manifest}")
    values = _manifest_values(manifest)
    _require("mesh_restart" in values, f"{manifest}: missing mesh_restart")
    mesh_relative = Path(values["mesh_restart"])
    _require(not mesh_relative.is_absolute(), f"{manifest}: absolute mesh restart path")
    return Checkpoint(
        root=root,
        shard_relative=shard.relative_to(root),
        manifest_relative=manifest_relative,
        mesh_relative=mesh_relative,
        witness_relative=Path(str(mesh_relative) + ".rmeta"),
    )


def _checkpoints(root: Path) -> list[Checkpoint]:
    shards = sorted((root / "prst/rank_00000000").glob("*.prst"))
    _require(shards, f"{root}: no serial particle restart checkpoints")
    return [_checkpoint(root, shard) for shard in shards]


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


def _restart_arguments(checkpoint: Checkpoint) -> list[str]:
    return [
        "-r",
        str(checkpoint.mesh_relative),
        "problem/prtcl_rst_flag=1",
        f"problem/prtcl_res_file={checkpoint.shard_relative}",
    ]


def _require_metadata(actual: dict[str, str], expected: dict[str, str],
                      label: str) -> None:
    for key, value in expected.items():
        _require(
            actual.get(key) == value,
            f"{label}: metadata {key}={actual.get(key)!r}, expected {value!r}",
        )


def _latest_output(case_dir: Path, directory: str, suffix: str) -> Path:
    paths = list((case_dir / directory).glob(f"*.{suffix}"))
    _require(paths, f"{case_dir}: missing {directory} output")
    return max(paths, key=tracer_inspect._file_number)


def _common_metadata() -> dict[str, str]:
    return {
        "schema": "akcr_particle_output_v1",
        "mode": "relativistic_hc",
        "units": "code_model",
        "c_model": "1.00000000000000000e+00",
        "alpha_s": "1.00000000000000000e+00",
    }


def _require_histogram_metadata(record: dict, expected: dict[str, str],
                                expected_vmin: float, expected_vmax: float,
                                label: str) -> None:
    metadata = record["metadata"]
    _require_metadata(metadata, {**_common_metadata(), **expected}, label)
    _near(float(metadata["vmin"]), expected_vmin, 1.0e-15, f"{label} vmin")
    _near(float(metadata["vmax"]), expected_vmax, 1.0e-15, f"{label} vmax")


def _retained_output_dictionary(case_dir: Path,
                                checkpoints: Sequence[Checkpoint]) -> None:
    ppd = tracer_inspect.summarize_ppd(case_dir)
    _require(ppd is not None, f"{case_dir}: missing ppd output")
    _require(ppd["count"] == 1, f"{case_dir}: ppd expected one particle")
    _require_metadata(
        ppd["metadata"],
        {
            **_common_metadata(),
            "columns": "species_x_y_z",
            "position_units": "code_length",
            "payload": "float32",
        },
        f"{case_dir} ppd",
    )

    histograms = tracer_inspect.validate_histograms(case_dir, [1], 8)
    df = tracer_inspect.read_df_record(histograms["df_file"], 1, 8)
    _require_histogram_metadata(
        df,
        {
            "quantity": "mu",
            "quantity_units": "dimensionless",
            "histogram_units": "particle_count",
            "nspecies": "1",
            "nbin": "8",
            "reduce": "1",
        },
        -1.0, 1.0, f"{case_dir} df",
    )
    dxh = tracer_inspect.read_dxh_record(histograms["dxh_file"], 1, 8)
    _require_histogram_metadata(
        dxh,
        {
            "quantity": "dx_dy_dz",
            "quantity_units": "code_length",
            "histogram_units": "particle_count",
            "nspecies": "1",
            "nbin": "8",
            "reduce": "1",
        },
        -1.0, 1.0, f"{case_dir} dxh",
    )
    scalar_specs = (
        (
            "drh", b"# AthenaK particle scalar displacement histogram",
            "displacement_norm", 0.0, 1.0, {},
        ),
        (
            "dparh", b"# AthenaK particle parallel displacement histogram",
            "dparallel", -1.0, 1.0,
            {"definition": "accumulated_midpoint_sum_dx_dot_Bhat"},
        ),
    )
    for name, marker, quantity, vmin, vmax, additions in scalar_specs:
        record = tracer_inspect.read_scalar_histogram_record(
            histograms[name]["file"], marker, 1, 8)
        _require_histogram_metadata(
            record,
            {
                "quantity": quantity,
                "quantity_units": "code_length",
                "histogram_units": "particle_count",
                "nspecies": "1",
                "nbin": "8",
                "reduce": "1",
                **additions,
            },
            vmin, vmax, f"{case_dir} {name}",
        )

    tracer_inspect.validate_moments(case_dir, [1])
    pmom_metadata = tracer_inspect.read_pmom_record(
        _latest_output(case_dir, "pmom", "pmom"), 1)["metadata"]
    _require_metadata(
        pmom_metadata,
        {
            **_common_metadata(),
            "quantity_basis": "physical_velocity_shadow",
            "mu_units": "dimensionless",
            "displacement_units": "code_length",
            "displacement_second_moment_units": "code_length_squared",
            "velocity_units": "code_velocity",
            "velocity_second_moment_units": "code_velocity_squared",
            "dparallel_definition": "accumulated_midpoint_sum_dx_dot_Bhat",
        },
        f"{case_dir} pmom",
    )

    joint = tracer_inspect.validate_joint_spectra(case_dir, [1], 8, 8)
    _require_metadata(
        joint["metadata"],
        {
            **_common_metadata(),
            "quantity": "mu_wmag",
            "nspecies": "1",
            "nbin1": "8",
            "nbin2": "8",
            "reduce": "1",
            "axis1_units": "dimensionless",
            "axis2_units": "code_velocity",
            "histogram_units": "particle_count",
        },
        f"{case_dir} pspec2",
    )
    for name, expected in (
            ("vmin1", -1.0), ("vmax1", 1.0), ("vmin2", 0.0), ("vmax2", 1.0)):
        _near(float(joint["metadata"][name]), expected, 1.0e-15,
              f"{case_dir} pspec2 {name}")

    restart = tracer_inspect.summarize_restart(case_dir)
    _require(restart["total"] == 1, f"{case_dir}: restart expected one particle")
    _require(
        restart["format_versions"] == {"AKPRST-v2.0": 1},
        f"{case_dir}: expected one typed-v2 restart shard",
    )
    _require(checkpoints, f"{case_dir}: expected paired restart checkpoints")
    for checkpoint in checkpoints:
        decoded = tracer_inspect.read_prst_file(
            checkpoint.root / checkpoint.shard_relative)
        _require(
            decoded["format_version"] == "AKPRST-v2.0",
            f"{case_dir}: expected typed-v2 checkpoint {checkpoint.shard_relative}",
        )
        _require(
            decoded["count"] == 1,
            f"{case_dir}: expected one particle in {checkpoint.shard_relative}",
        )
        for relative in (
                checkpoint.shard_relative, checkpoint.manifest_relative,
                checkpoint.mesh_relative, checkpoint.witness_relative):
            _require(
                (checkpoint.root / relative).is_file(),
                f"{case_dir}: missing paired restart artifact {relative}",
            )


def _sample(case_dir: Path,
            expected_metadata: dict[str, str] = PSAMP_METADATA) -> tuple[dict, dict]:
    summary = tracer_inspect.summarize_samples(case_dir)
    _require(summary is not None, f"{case_dir}: missing psamp output")
    _require(summary["total"] == 1, f"{case_dir}: expected one sampled particle")
    _require(len(summary["blocks"]) == 1, f"{case_dir}: expected one psamp shard")
    block = summary["blocks"][0]
    expected_columns = ["rank", "pgid", "tag", "species", *SAMPLE_FIELDS]
    _require(
        block["columns"] == expected_columns,
        f"{case_dir}: sample columns {block['columns']} != {expected_columns}",
    )
    _require_metadata(block["metadata"], expected_metadata, f"{case_dir} psamp")
    return block["rows"][0], block


def _spectrum(case_dir: Path) -> dict:
    spectrum = tracer_inspect.validate_spectra(case_dir, [1], 8)
    _require_metadata(spectrum["metadata"], PSPEC_METADATA, f"{case_dir} pspec")
    _near(float(spectrum["metadata"]["vmin"]), 0.0, 0.0, f"{case_dir} pspec vmin")
    _near(float(spectrum["metadata"]["vmax"]), 0.8, 1.0e-15, f"{case_dir} pspec vmax")
    return spectrum


def _occupied_bin(spectrum: dict, expected_wmag: float, label: str) -> int:
    metadata = spectrum["metadata"]
    nbin = int(metadata["nbin"])
    vmin = float(metadata["vmin"])
    vmax = float(metadata["vmax"])
    expected = int((expected_wmag - vmin) / (vmax - vmin) * nbin)
    expected = min(nbin - 1, max(0, expected))
    occupied = [
        index for index, count in enumerate(spectrum["histogram"][0]) if count]
    _require(occupied == [expected], f"{label}: bins {occupied}, expected {[expected]}")
    _require(
        spectrum["histogram"][0][expected] == 1,
        f"{label}: expected exactly one particle in bin {expected}",
    )
    return expected


def _qualify_row(row: dict, initial_w: Sequence[float], expected_ce: Sequence[float],
                 label: str, c_model: float = 1.0,
                 alpha_s: float = 1.0) -> dict[str, float]:
    fields = row["fields"]
    w = tuple(fields[name] for name in ("wx", "wy", "wz"))
    velocity = tuple(fields[name] for name in ("vx", "vy", "vz"))
    gamma = math.sqrt(1.0 + sum(value * value for value in w) / c_model**2)
    wmag = math.sqrt(sum(value * value for value in w))
    kinetic = (gamma - 1.0) * c_model**2
    initial_gamma = math.sqrt(
        1.0 + sum(value * value for value in initial_w) / c_model**2)
    delta_kinetic = kinetic - (initial_gamma - 1.0) * c_model**2
    return {
        "gamma_error": _near(fields["gamma"], gamma, 3.0e-13, f"{label} gamma"),
        "velocity_shadow_error": _vector_error(
            velocity, tuple(value / gamma for value in w), 3.0e-13,
            f"{label} v=w/gamma"),
        "wmag_error": _near(fields["wmag"], wmag, 3.0e-13, f"{label} wmag"),
        "kinetic_energy_model_error": _near(
            fields["kinetic_energy_model"], kinetic, 3.0e-13,
            f"{label} kinetic_energy_model"),
        "work_closure_error": _near(
            fields["work"], delta_kinetic, 3.0e-13,
            f"{label} delta work=delta kinetic energy"),
        "sampled_ce_error": _vector_error(
            tuple(fields[name] for name in ("cex", "cey", "cez")),
            expected_ce, 3.0e-13, f"{label} sampled cE"),
        "alpha_s_error": _near(
            fields["alpha_s"], alpha_s, 3.0e-13, f"{label} alpha_s"),
        "applicability_ratio_error": _near(
            fields["r_larmor_over_dx_min"], 0.0, 3.0e-13,
            f"{label} r_larmor_over_dx_min"),
    }


def _fields_equal(left: dict[str, float], right: dict[str, float],
                  tolerance: float, label: str) -> float:
    error = max(abs(left[name] - right[name]) for name in SAMPLE_FIELDS)
    _require(error <= tolerance, f"{label}: maximum field error={error:.6e}")
    return error


def _validate_criteria(path: Path, metrics: dict[str, object]) -> None:
    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    _require(
        digest == EXPECTED_CRITERIA_SHA256,
        f"{path}: criteria SHA-256 {digest} != frozen {EXPECTED_CRITERIA_SHA256}",
    )
    criteria = json.loads(payload)
    _require(
        criteria.get("manifest_schema_version") == 1,
        f"{path}: unsupported criteria manifest schema",
    )
    criteria_ids = tuple(item["id"] for item in criteria["criteria"])
    _require(
        criteria_ids == EXPECTED_CRITERIA_IDS,
        f"{path}: criteria IDs {criteria_ids!r} != frozen {EXPECTED_CRITERIA_IDS!r}",
    )
    registered = {item["metric"] for item in criteria["criteria"]}
    emitted = set(metrics)
    _require(
        emitted == registered,
        f"criteria completeness mismatch: unbound={sorted(emitted - registered)}, "
        f"missing={sorted(registered - emitted)}",
    )
    for item in criteria["criteria"]:
        actual = metrics[item["metric"]]
        limit = item["limit"]
        operator = item["operator"]
        accepted = (
            (operator == "<=" and actual <= limit)
            or (operator == ">=" and actual >= limit)
            or (operator == "==" and actual == limit)
        )
        _require(
            accepted,
            f"{item['id']} failed: {item['metric']}={actual!r} "
            f"{operator} {limit!r} is false",
        )


@contextlib.contextmanager
def _runtime_root(work_dir: Optional[Path]) -> Iterator[Path]:
    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="cr-rel-phase7-") as temporary:
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
    acceleration_dir = root / "acceleration_uninterrupted"
    _run_success(binary, acceleration_dir, ["-i", str(input_path)], timeout,
                 "acceleration_uninterrupted")
    acceleration_row, acceleration_block = _sample(acceleration_dir)
    acceleration_spectrum = _spectrum(acceleration_dir)
    acceleration = _qualify_row(
        acceleration_row, (0.2, 0.0, 0.0), (0.3, 0.0, 0.0), "acceleration")
    acceleration_bin = _occupied_bin(
        acceleration_spectrum, acceleration_row["fields"]["wmag"], "acceleration")
    _require(acceleration_bin == 2, f"acceleration: occupied bin {acceleration_bin} != 2")

    deceleration_dir = root / "deceleration"
    _run_success(
        binary,
        deceleration_dir,
        ["-i", str(input_path), *_flags(time__nlim=1, problem__cE0x=-0.3)],
        timeout,
        "deceleration",
    )
    deceleration_row, deceleration_block = _sample(deceleration_dir)
    deceleration_spectrum = _spectrum(deceleration_dir)
    deceleration = _qualify_row(
        deceleration_row, (0.2, 0.0, 0.0), (-0.3, 0.0, 0.0), "deceleration")
    _require(
        deceleration_row["fields"]["work"] < 0.0,
        "deceleration: expected negative accumulated work",
    )
    deceleration_bin = _occupied_bin(
        deceleration_spectrum, deceleration_row["fields"]["wmag"], "deceleration")
    _require(deceleration_bin == 1, f"deceleration: occupied bin {deceleration_bin} != 1")

    nonunit_dir = root / "nonunit_c_model"
    _run_success(
        binary,
        nonunit_dir,
        ["-i", str(input_path), *_flags(particles__c_model=2.0, time__nlim=1)],
        timeout,
        "nonunit_c_model",
    )
    nonunit_metadata = {
        **PSAMP_METADATA,
        "c_model": "2.00000000000000000e+00",
    }
    nonunit_row, _ = _sample(nonunit_dir, nonunit_metadata)
    nonunit = _qualify_row(
        nonunit_row, (0.2, 0.0, 0.0), (0.3, 0.0, 0.0),
        "nonunit c_model", c_model=2.0)
    _require(
        nonunit_row["fields"]["work"] > 0.0,
        "nonunit c_model: expected positive accumulated work",
    )

    checkpoints = _checkpoints(acceleration_dir)
    _require(
        len(checkpoints) >= 2,
        "acceleration: expected at least two restart checkpoints",
    )
    _retained_output_dictionary(acceleration_dir, checkpoints)
    continuation_dir = root / "paired_continuation"
    _copy_checkpoint(checkpoints[0], continuation_dir)
    _run_success(
        binary,
        continuation_dir,
        _restart_arguments(checkpoints[0]),
        timeout,
        "paired_continuation",
    )
    continuation_row, continuation_block = _sample(continuation_dir)
    continuation_spectrum = _spectrum(continuation_dir)
    continuation = _qualify_row(
        continuation_row, (0.2, 0.0, 0.0), (0.3, 0.0, 0.0), "paired continuation")
    restart_field_error = _fields_equal(
        acceleration_row["fields"], continuation_row["fields"], 3.0e-13,
        "paired continuation diagnostic fields",
    )
    _require(
        acceleration_spectrum["histogram"] == continuation_spectrum["histogram"],
        "paired continuation: final spectrum differs from uninterrupted source",
    )
    _require(
        acceleration_block["columns"] == continuation_block["columns"],
        "paired continuation: canonical psamp columns changed across restart",
    )

    metrics = {
        **{f"acceleration_{name}": value for name, value in acceleration.items()},
        "acceleration_sample_schema_exact": True,
        "acceleration_psamp_metadata_explicit": True,
        "acceleration_pspec_metadata_explicit": True,
        "acceleration_spectrum_bin": acceleration_bin,
        **{f"deceleration_{name}": value for name, value in deceleration.items()},
        "deceleration_work_negative": True,
        "deceleration_sample_schema_exact": True,
        "deceleration_psamp_metadata_explicit": True,
        "deceleration_pspec_metadata_explicit": True,
        "deceleration_spectrum_bin": deceleration_bin,
        **{f"nonunit_{name}": value for name, value in nonunit.items()},
        "nonunit_psamp_metadata_explicit": True,
        "nonunit_work_positive": True,
        "restart_diagnostic_field_error": restart_field_error,
        "restart_spectrum_equal": True,
        "restart_sample_schema_equal": True,
        "restart_work_closure_error": continuation["work_closure_error"],
        "retained_output_metadata_explicit": True,
        "retained_output_inspector_round_trip": True,
        "paired_restart_bundle_present": True,
        "restart_artifacts_inspected": True,
    }
    _validate_criteria(criteria_path, metrics)
    return {
        "schema_version": 1,
        "qualification": "cr_relativistic_diagnostics_runtime_phase7",
        "status": "PASS",
        "binary": str(binary),
        "binary_sha256": _sha256(binary),
        "input": str(input_path),
        "input_sha256": _sha256(input_path),
        "criteria": str(criteria_path),
        "criteria_sha256": _sha256(criteria_path),
        "metrics": metrics,
        "runtime_cases": [
            str(acceleration_dir),
            str(deceleration_dir),
            str(nonunit_dir),
            str(continuation_dir),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Qualify serial Phase-7 relativistic CR diagnostic meanings.")
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--metrics", type=Path)
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args()
    try:
        with _runtime_root(args.work_dir) as root:
            report = _qualify(
                args.binary.resolve(),
                args.input.resolve(),
                args.criteria.resolve(),
                root,
                args.timeout,
            )
    except Exception as error:  # Keep CLI failure concise; per-case logs retain detail.
        print(f"CR relativistic Phase-7 diagnostics qualification FAIL: {error}",
              file=sys.stderr)
        return 1

    if args.metrics is not None:
        path = args.metrics.resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("CR relativistic Phase-7 diagnostics qualification PASS")
    for metric, value in sorted(report["metrics"].items()):
        print(f"  {metric}: {value}")
    if args.work_dir is not None:
        print(f"  preserved runtime cases: {args.work_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
