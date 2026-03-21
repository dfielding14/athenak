#!/usr/bin/env python3
"""Run focused 3D Clebsch audits against the diagnostics sidecar."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from collections import Counter
from pathlib import Path

SEEDS = (123, 321, 777)
SECTOR_ORDER = ("local_local", "low_high", "high_low", "high_high_cancel")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_athena_path() -> Path:
    for relpath in (
        Path("build-mpi/src/athena"),
        Path("build-codex-spectrum/src/athena"),
        Path("build/src/athena"),
    ):
        candidate = _repo_root() / relpath
        if candidate.exists():
            return candidate
    return _repo_root() / "build-codex-spectrum" / "src" / "athena"


def _run_root(athena: Path) -> Path:
    return athena.parent.parent if athena.parent.name == "src" else athena.parent


def _default_cases() -> list[dict[str, object]]:
    return [
        {
            "name": "small_unthinned",
            "run_nhigh": 6,
            "analysis_nhigh": 6,
            "k_crit": 16.0,
        },
        {
            "name": "small_thinned",
            "run_nhigh": 6,
            "analysis_nhigh": 6,
            "k_crit": 2.0,
        },
        {
            "name": "small_thinned_2x_band",
            "run_nhigh": 12,
            "analysis_nhigh": 6,
            "k_crit": 2.0,
        },
        {
            "name": "small_thinned_4x_band",
            "run_nhigh": 16,
            "analysis_nhigh": 6,
            "k_crit": 2.0,
        },
    ]


def _production_cases(repo_root: Path) -> list[dict[str, object]]:
    input_path = repo_root / "inputs" / "hydro" / "scalar_mixing_methods_paper_3d.athinput"
    return [
        {
            "name": "prod_thinned",
            "input": input_path,
            "run_nhigh": 128,
            "analysis_nhigh": 128,
            "k_crit": 2.0,
            "overrides": ["problem/turb_use_stream_function=true"],
        },
    ]


def _sector_fraction_summary(diag: dict[str, object], response_key: str,
                             analysis_nhigh: int) -> dict[str, object]:
    response = diag["sector_response"][response_key]
    totals = {}
    total_sum = 0.0
    for sector in SECTOR_ORDER:
        shell_energy = response[sector]["shell_energy"]
        kmax = min(analysis_nhigh, len(shell_energy) - 1)
        subtotal = float(sum(shell_energy[1:kmax + 1]))
        totals[sector] = subtotal
        total_sum += subtotal
    fractions = {
        sector: (totals[sector] / total_sum if total_sum > 0.0 else 0.0)
        for sector in SECTOR_ORDER
    }
    dominant_sector = max(SECTOR_ORDER, key=lambda sector: fractions[sector])
    return {
        "totals": totals,
        "fractions": fractions,
        "dominant_sector": dominant_sector,
        "analysis_nhigh": analysis_nhigh,
    }


def _run_case(athena: Path, case: dict[str, object], seed: int,
              input_path: Path, launcher: str, output_root: Path) -> dict[str, object]:
    input_override = Path(case.get("input", input_path)).resolve()
    case_run_dir = output_root / str(case["name"]) / f"seed_{seed}"
    case_run_dir.mkdir(parents=True, exist_ok=True)
    basename = f"clebsch_audit_{case['name']}_{seed}"
    diag_path = case_run_dir / f"{basename}.turb_init_diag.json"
    if diag_path.exists():
        diag_path.unlink()

    cmd = shlex.split(launcher) + [
        str(athena),
        "-i",
        str(input_override),
        f"job/basename={basename}",
        f"problem/turb_rseed={seed}",
        f"problem/turb_nhigh={case['run_nhigh']}",
        f"problem/turb_k_crit={case['k_crit']}",
    ] + list(case.get("overrides", []))
    completed = subprocess.run(
        cmd,
        cwd=case_run_dir,
        capture_output=True,
        text=True,
        check=True,
    )
    if not diag_path.exists():
        raise RuntimeError(f"Missing diagnostics JSON for {basename}")
    diag = json.loads(diag_path.read_text(encoding="utf-8"))
    fit = diag["fit"]
    quality_gate = diag["quality_gate"]
    fit_sector = _sector_fraction_summary(
        diag, "fitted_response", int(case["analysis_nhigh"]))
    raw_sector = _sector_fraction_summary(
        diag, "raw_tensor", int(case["analysis_nhigh"]))
    return {
        "case": case["name"],
        "seed": seed,
        "run_nhigh": int(case["run_nhigh"]),
        "analysis_nhigh": int(case["analysis_nhigh"]),
        "k_crit": float(case["k_crit"]),
        "phi_shell_max_used": int(diag["spectra"]["phi_shell_max_used"]),
        "seed_model_chosen": fit["seed_model_chosen"],
        "fit_starts": fit["fit_starts"],
        "deterministic_slope": float(fit["fitted_slope"]),
        "deterministic_leakage": float(fit["estimated_out_of_band_leakage"]),
        "realized_slope": float(fit["realized_slope"]),
        "realized_leakage": float(fit["realized_out_of_band_leakage"]),
        "quality_gate": quality_gate,
        "sector_raw": raw_sector,
        "sector_fitted": fit_sector,
        "init_wall_seconds": float(diag["metadata"]["init_wall_seconds"]),
        "run_dir": str(case_run_dir),
        "stdout_tail": completed.stdout.splitlines()[-12:],
    }


def _aggregate(records: list[dict[str, object]]) -> dict[str, object]:
    by_case: dict[str, list[dict[str, object]]] = {}
    for record in records:
        by_case.setdefault(str(record["case"]), []).append(record)

    summary = {}
    for case_name, items in by_case.items():
        seed_model_counts = Counter(str(item["seed_model_chosen"]) for item in items)
        dominant_raw = Counter(str(item["sector_raw"]["dominant_sector"]) for item in items)
        dominant_fitted = Counter(str(item["sector_fitted"]["dominant_sector"]) for item in items)
        summary[case_name] = {
            "run_nhigh": int(items[0]["run_nhigh"]),
            "analysis_nhigh": int(items[0]["analysis_nhigh"]),
            "k_crit": float(items[0]["k_crit"]),
            "phi_shell_max_used": sorted({int(item["phi_shell_max_used"]) for item in items}),
            "pass_fraction": float(sum(bool(item["quality_gate"]["pass"]) for item in items)) / len(items),
            "forced_continue_fraction": float(
                sum(bool(item["quality_gate"]["forced_continue"]) for item in items)
            ) / len(items),
            "mean_deterministic_slope": float(
                sum(float(item["deterministic_slope"]) for item in items) / len(items)
            ),
            "mean_realized_slope": float(
                sum(float(item["realized_slope"]) for item in items) / len(items)
            ),
            "mean_deterministic_leakage": float(
                sum(float(item["deterministic_leakage"]) for item in items) / len(items)
            ),
            "mean_realized_leakage": float(
                sum(float(item["realized_leakage"]) for item in items) / len(items)
            ),
            "mean_init_wall_seconds": float(
                sum(float(item["init_wall_seconds"]) for item in items) / len(items)
            ),
            "seed_model_counts": dict(seed_model_counts),
            "dominant_raw_sector_counts": dict(dominant_raw),
            "dominant_fitted_sector_counts": dict(dominant_fitted),
        }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--athena",
        type=Path,
        default=_default_athena_path(),
        help="Path to AthenaK built with PROBLEM=scalar_mixing.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=_repo_root() / "inputs" / "tests" / "scalar_mixing_stream_3d_audit.athinput",
        help="Audit input deck with diagnostics enabled and poor-fit continuation allowed.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_repo_root() / "tmp" / "clebsch_audit_report.json",
        help="Path to the JSON audit report.",
    )
    parser.add_argument(
        "--profile",
        choices=("small", "production"),
        default="small",
        help="Audit case profile to run.",
    )
    parser.add_argument(
        "--launcher",
        default="",
        help="Optional launcher prefix, for example 'mpiexec -n 8'.",
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        default=_repo_root() / "tmp" / "clebsch_audit_runs",
        help="Directory that will receive per-case diagnostics and outputs.",
    )
    args = parser.parse_args()

    repo_root = _repo_root()
    cases = (_default_cases() if args.profile == "small"
             else _production_cases(repo_root))
    records = []
    for case in cases:
        for seed in SEEDS:
            print(f"[audit] case={case['name']} seed={seed}")
            records.append(
                _run_case(
                    args.athena.resolve(),
                    case,
                    seed,
                    args.input.resolve(),
                    args.launcher,
                    args.run_root.resolve(),
                )
            )

    report = {
        "profile": args.profile,
        "notes": [
            "The nested cutoff cases vary run_nhigh while reporting sector totals only through analysis_nhigh.",
            "This approximates low-band cutoff sensitivity without exposing a production turb_stream_phi_nhigh knob.",
        ],
        "seeds": list(SEEDS),
        "records": records,
        "summary": _aggregate(records),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"[audit] wrote {args.output}")


if __name__ == "__main__":
    main()
