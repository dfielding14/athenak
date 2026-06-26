#!/usr/bin/env python3
"""Run the compact Bell-current engineering-readiness matrix.

This is an operational gate for starting Bell development runs. It reuses the
Q043 decks and exact-byte raw oracle, but it is deliberately not publication
qualification and does not replace the complete 132-case matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Mapping, Sequence

from tst.publication import (
    q043_bell_current_volume_aware_deposited_current_oracle as oracle,
)


SCHEMA_VERSION = 1
RECORD_TYPE = "q043_bell_current_engineering_readiness_manifest"
SUMMARY_RECORD_TYPE = "q043_bell_current_engineering_readiness_summary"
QUALIFICATION_EFFECT = "engineering_readiness_only_not_publication_qualification"
BASELINE_C_OVER_V_CR = 1000
ENDPOINT_C_OVER_V_CR = (100, 10000)


class ReadinessError(RuntimeError):
    """Reject drifted selection, execution, or raw-output evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ReadinessError(message)


def _canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selected_cases() -> tuple[dict[str, object], ...]:
    """Return the 26 orthogonal cases used by the engineering gate."""
    selected: list[dict[str, object]] = []
    for case in oracle.expected_cases():
        resolution = str(case["resolution"])
        ppc = int(case["ppc"])
        decomposition = str(case["decomposition"])
        c_ratio = int(case["artificial_c_over_v_cr"])

        baseline_scale = (
            decomposition == "single" and c_ratio == BASELINE_C_OVER_V_CR
        )
        artificial_c_endpoint = (
            resolution == "coarse"
            and ppc == 1
            and decomposition == "single"
            and c_ratio in ENDPOINT_C_OVER_V_CR
        )
        decomposition_check = (
            resolution == "fine"
            and ppc == 4
            and decomposition != "single"
            and c_ratio == BASELINE_C_OVER_V_CR
        )
        if baseline_scale or artificial_c_endpoint or decomposition_check:
            selected.append(case)

    _validate_selection(selected)
    return tuple(selected)


def _validate_selection(cases: Sequence[Mapping[str, object]]) -> None:
    _require(len(cases) == 26, "engineering-readiness selection must contain 26 cases")
    ids = [str(case["case_id"]) for case in cases]
    _require(len(ids) == len(set(ids)), "engineering-readiness case IDs are not unique")

    baseline = [
        case
        for case in cases
        if case["decomposition"] == "single"
        and case["artificial_c_over_v_cr"] == BASELINE_C_OVER_V_CR
    ]
    _require(
        {
            (int(case["dimension"]), str(case["resolution"]), int(case["ppc"]))
            for case in baseline
        }
        == {
            (dimension, resolution, ppc)
            for dimension in oracle.DIMENSIONS
            for resolution in oracle.RESOLUTIONS
            for ppc in oracle.PPC_VALUES
        },
        "baseline dimension/resolution/PPC coverage drifted",
    )

    endpoint = [
        case
        for case in cases
        if case["resolution"] == "coarse"
        and case["ppc"] == 1
        and case["decomposition"] == "single"
        and case["artificial_c_over_v_cr"] in ENDPOINT_C_OVER_V_CR
    ]
    _require(
        {
            (int(case["dimension"]), int(case["artificial_c_over_v_cr"]))
            for case in endpoint
        }
        == {
            (dimension, c_ratio)
            for dimension in oracle.DIMENSIONS
            for c_ratio in ENDPOINT_C_OVER_V_CR
        },
        "artificial-C endpoint coverage drifted",
    )

    split = [case for case in cases if case["decomposition"] != "single"]
    _require(
        {
            (int(case["dimension"]), str(case["decomposition"])) for case in split
        }
        == {
            (dimension, decomposition)
            for dimension in oracle.DIMENSIONS
            for decomposition in oracle.DECOMPOSITIONS_BY_DIMENSION[dimension]
            if decomposition != "single"
        },
        "MPI decomposition coverage drifted",
    )


def build_manifest() -> dict[str, object]:
    """Bind the compact selection to the canonical Q043 deck bytes."""
    canonical = oracle.validate_checked_in_decks()
    by_id = {str(case["case_id"]): case for case in canonical["cases"]}
    records = []
    for case in selected_cases():
        case_id = str(case["case_id"])
        deck = oracle.CHECKED_IN_DECK_ROOT / f"{case_id}.athinput"
        canonical_record = by_id[case_id]
        _require(deck.is_file(), f"{case_id}: checked-in deck is missing")
        _require(
            _sha256(deck) == canonical_record["deck_sha256"],
            f"{case_id}: checked-in deck digest drifted",
        )
        records.append(
            {
                "case_id": case_id,
                "dimension": int(case["dimension"]),
                "resolution": str(case["resolution"]),
                "ppc": int(case["ppc"]),
                "decomposition": str(case["decomposition"]),
                "mpi_ranks": int(case["mpi_ranks"]),
                "artificial_c_over_v_cr": int(case["artificial_c_over_v_cr"]),
                "deck_path": str(deck.relative_to(oracle.REPO_ROOT)),
                "deck_sha256": _sha256(deck),
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "campaign_id": "Q043-BELL-CURRENT-ENGINEERING-READINESS-V1",
        "qualification_effect": QUALIFICATION_EFFECT,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "replaces_complete_q043_matrix": False,
        "canonical_q043_case_count": oracle.EXPECTED_CASE_COUNT,
        "engineering_case_count": len(records),
        "selection_design": {
            "baseline_scale_matrix": (
                "3 dimensions x 2 resolutions x 2 PPC at C/v_CR=1000 "
                "and single decomposition"
            ),
            "artificial_c_endpoints": (
                "C/v_CR=100 and 10000 at coarse PPC=1 single decomposition "
                "in each dimension"
            ),
            "mpi_decompositions": (
                "every non-single dimension-valid decomposition at fine PPC=4 "
                "and C/v_CR=1000"
            ),
        },
        "execution_deck_derivation": (
            "canonical_q043_deck_plus_exact_species0_velocity_parameters"
        ),
        "cases": records,
    }


def _execution_deck_bytes(case: Mapping[str, object]) -> bytes:
    """Render the exact species velocity overrides required by the pgen."""
    base = oracle.render_oracle_deck(case)
    blocks = oracle.parse_athinput_text(base)
    particles = blocks["particles"]
    marker = "\n<problem>\n"
    insertion = "".join(
        f"v{axis}0 = {particles[f'cr_v{axis}0']}\n"
        for axis in ("x", "y", "z")
    )
    _require(
        base.count(marker) == 1,
        f"{case['case_id']}: species velocity insertion point drifted",
    )
    rendered = base.replace(marker, insertion + marker, 1)
    try:
        oracle.validate_rendered_deck(case, rendered)
    except oracle.ContractError as error:
        raise ReadinessError(
            f"{case['case_id']}: derived execution deck is invalid"
        ) from error
    return rendered.encode("utf-8")


def _command(
    case: Mapping[str, object],
    *,
    executable: Path,
    execution_deck: Path,
    raw_root: Path,
    launcher: str,
) -> list[str]:
    ranks = int(case["mpi_ranks"])
    prefix = []
    if launcher:
        prefix = [
            launcher,
            "--nodes=1",
            f"--ntasks={ranks}",
            f"--ntasks-per-node={ranks}",
            "--cpus-per-task=1",
            "--gpus-per-task=1",
            "--gpu-bind=closest",
        ]
    return [
        *prefix,
        str(executable),
        "-i",
        str(execution_deck),
        "-d",
        str(raw_root),
        "time/nlim=1",
    ]


def _cycle_one_paths(
    case: Mapping[str, object], raw_root: Path
) -> dict[str, tuple[Path, ...]]:
    basename = str(case["case_id"]).replace("-", "_")
    ranks = int(case["mpi_ranks"])
    result: dict[str, tuple[Path, ...]] = {}
    for field in oracle.FIELDS:
        filename = f"{basename}.{field}.00001.bin"
        if ranks == 1:
            paths = (raw_root / "bin" / filename,)
        else:
            paths = tuple(
                raw_root / "bin" / f"rank_{rank:08d}" / filename
                for rank in range(ranks)
            )
        missing = [str(path) for path in paths if not path.is_file()]
        _require(not missing, f"{case['case_id']}: missing cycle-one {field}: {missing}")
        result[field] = paths
    return result


def _write_json(path: Path, value: object) -> None:
    path.write_bytes(_canonical_bytes(value))


def _campaign_summary(results: Sequence[Mapping[str, object]]) -> dict[str, object]:
    projected = [
        float(result["measured_guide_projected_volume_mean_j_over_c"])
        for result in results
    ]
    maximum_tolerance = max(
        float(result["representation_derived_absolute_tolerance"])
        for result in results
    )
    spread = max(projected) - min(projected)
    _require(spread <= maximum_tolerance, "engineering cross-case current spread failed")
    expected_ids = [str(case["case_id"]) for case in selected_cases()]
    observed_ids = [str(result["case_id"]) for result in results]
    _require(observed_ids == expected_ids, "engineering result order or membership drifted")
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": SUMMARY_RECORD_TYPE,
        "campaign_id": "Q043-BELL-CURRENT-ENGINEERING-READINESS-V1",
        "qualification_effect": QUALIFICATION_EFFECT,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "engineering_ready": True,
        "case_count": len(results),
        "maximum_cross_case_projected_current_spread": spread,
        "maximum_allowed_absolute_tolerance": maximum_tolerance,
        "case_ids": observed_ids,
    }


def run_campaign(
    *,
    executable: Path,
    output_root: Path,
    launcher: str = "/usr/bin/srun",
    dry_run: bool = False,
) -> dict[str, object] | None:
    """Execute all 26 cases serially inside one existing allocation."""
    _require(executable.is_file(), f"executable is unavailable: {executable}")
    _require(not output_root.exists(), f"output root already exists: {output_root}")
    output_root.mkdir(parents=True)
    _write_json(output_root / "engineering_readiness_manifest.json", build_manifest())
    results = []
    for case in selected_cases():
        case_id = str(case["case_id"])
        case_root = output_root / case_id
        raw_root = case_root / "raw"
        raw_root.mkdir(parents=True)
        execution_deck = case_root / f"{case_id}.athinput"
        execution_deck.write_bytes(_execution_deck_bytes(case))
        command = _command(
            case,
            executable=executable,
            execution_deck=execution_deck,
            raw_root=raw_root,
            launcher=launcher,
        )
        (case_root / "command.json").write_bytes(_canonical_bytes(command))
        if dry_run:
            continue
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "1",
                "MPICH_GPU_SUPPORT_ENABLED": "1",
                "MPICH_GPU_EAGER_REGISTER_HOST_MEM": "0",
                "MPICH_GPU_NO_ASYNC_COPY": "1",
                "MPICH_GPU_IPC_ENABLED": "0",
            }
        )
        with (case_root / "stdout.txt").open("wb") as stdout, (
            case_root / "stderr.txt"
        ).open("wb") as stderr:
            completed = subprocess.run(
                command,
                cwd=oracle.REPO_ROOT,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        _require(completed.returncode == 0, f"{case_id}: AthenaK exited nonzero")
        result = oracle.analyze_raw_case(case_id, _cycle_one_paths(case, raw_root))
        _require(
            result["source_local_oracle_check_pass"] is True,
            f"{case_id}: exact-byte current oracle did not pass",
        )
        _write_json(case_root / "analysis.json", result)
        results.append(result)
    if dry_run:
        return None
    summary = _campaign_summary(results)
    _write_json(output_root / "engineering_readiness_summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--print-manifest", action="store_true")
    action.add_argument("--validate-selection", action="store_true")
    action.add_argument("--run", action="store_true")
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--launcher", default="/usr/bin/srun")
    parser.add_argument("--dry-run", action="store_true")
    arguments = parser.parse_args()
    if arguments.print_manifest:
        print(json.dumps(build_manifest(), indent=2, sort_keys=True))
        return 0
    if arguments.validate_selection:
        manifest = build_manifest()
        print(
            json.dumps(
                {
                    "engineering_case_count": manifest["engineering_case_count"],
                    "status": "valid",
                },
                sort_keys=True,
            )
        )
        return 0
    _require(arguments.executable is not None, "--run requires --executable")
    _require(arguments.output_root is not None, "--run requires --output-root")
    summary = run_campaign(
        executable=arguments.executable.resolve(),
        output_root=arguments.output_root.resolve(),
        launcher=arguments.launcher,
        dry_run=arguments.dry_run,
    )
    if summary is not None:
        print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
