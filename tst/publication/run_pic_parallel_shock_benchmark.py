#!/usr/bin/env python3
"""Run and compare the PIC parallel-shock benchmark trio."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import re
import subprocess
import time

import matplotlib.pyplot as plt
import numpy as np

import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402

from pvtk_particles import read_particle_vtk  # noqa: E402

_CYCLE_RE = re.compile(r"\.(\d+)\.bin$")


@dataclass
class CaseSpec:
    label: str
    deck: str
    basename: str


class MissingOutputsError(RuntimeError):
    """Raised when expected benchmark artifacts are missing for analysis."""


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _files_by_cycle(bin_dir: Path, basename: str, file_id: str) -> dict[int, Path]:
    out: dict[int, Path] = {}
    for path in sorted(bin_dir.glob(f"{basename}.{file_id}.*.bin")):
        match = _CYCLE_RE.search(path.name)
        if match is None:
            continue
        out[int(match.group(1))] = path
    return out


def _latest_pvtk(pvtk_dir: Path, basename: str) -> Path | None:
    files = sorted(pvtk_dir.glob(f"{basename}.prtcl_all.*.part.vtk"))
    if not files:
        return None
    return files[-1]


def _load_pmag(pvtk_path: Path) -> np.ndarray:
    pdata = read_particle_vtk(pvtk_path)
    if "vel" in pdata.vectors:
        vel = np.asarray(pdata.vectors["vel"], dtype=float)
    elif pdata.vectors:
        vname = sorted(pdata.vectors.keys())[0]
        vel = np.asarray(pdata.vectors[vname], dtype=float)
    else:
        return np.array([], dtype=float)
    return np.linalg.norm(vel, axis=1)


def _compute_case_metrics(bin_dir: Path, pvtk_dir: Path, case: CaseSpec) -> dict:
    rho_map = _files_by_cycle(bin_dir, case.basename, "rho")
    bmag_map = _files_by_cycle(bin_dir, case.basename, "bmag")
    shared = sorted(set(rho_map) & set(bmag_map))
    if len(shared) < 2:
        raise MissingOutputsError(
                f"case={case.label}: need >=2 shared rho/bmag outputs, got {len(shared)}\n"
                f"  searched bin_dir={bin_dir}\n"
                f"  expected patterns:\n"
                f"    {case.basename}.rho.*.bin\n"
                f"    {case.basename}.bmag.*.bin\n"
                f"  found: rho={len(rho_map)} bmag={len(bmag_map)} shared={len(shared)}"
        )

    b0 = np.asarray(
            bin_convert.read_binary_as_athdf(str(bmag_map[shared[0]]))["bmag"],
            dtype=float,
    )
    b1 = np.asarray(
            bin_convert.read_binary_as_athdf(str(bmag_map[shared[-1]]))["bmag"],
            dtype=float,
    )

    amp_mean = float(np.mean(b1) / max(np.mean(b0), 1.0e-30))
    amp_p95 = float(
            np.percentile(b1, 95.0) / max(np.percentile(b0, 95.0), 1.0e-30)
    )

    rho_latest = np.asarray(
            bin_convert.read_binary_as_athdf(str(rho_map[shared[-1]]))["dens"],
            dtype=float,
    )
    if rho_latest.ndim == 3:
        rho_x = np.mean(rho_latest, axis=(0, 1))
    elif rho_latest.ndim == 2:
        rho_x = np.mean(rho_latest, axis=0)
    else:
        rho_x = rho_latest
    x = np.asarray(
            bin_convert.read_binary_as_athdf(str(rho_map[shared[-1]]))["x1v"],
            dtype=float,
    )
    grad = np.gradient(rho_x, x)
    shock_x = float(x[int(np.argmax(np.abs(grad)))])

    pvtk_path = _latest_pvtk(pvtk_dir, case.basename)
    if pvtk_path is None:
        p_mag = np.array([], dtype=float)
    else:
        p_mag = _load_pmag(pvtk_path)

    if p_mag.size > 0:
        p50 = float(np.percentile(p_mag, 50.0))
        p90 = float(np.percentile(p_mag, 90.0))
        p99 = float(np.percentile(p_mag, 99.0))
    else:
        p50 = float("nan")
        p90 = float("nan")
        p99 = float("nan")

    return {
            "label": case.label,
            "deck": case.deck,
            "basename": case.basename,
            "n_shared_cycles": len(shared),
            "first_cycle": int(shared[0]),
            "last_cycle": int(shared[-1]),
            "b_amp_mean_ratio": amp_mean,
            "b_amp_p95_ratio": amp_p95,
            "shock_x_last": shock_x,
            "n_particles_last": int(p_mag.size),
            "p50_last": p50,
            "p90_last": p90,
            "p99_last": p99,
            "pvtk_last": "" if pvtk_path is None else str(pvtk_path),
            "p_mag_last": p_mag,
    }


def _write_summary_csv(path: Path, rows: list[dict]) -> None:
    fields = [
            "label",
            "deck",
            "basename",
            "n_shared_cycles",
            "first_cycle",
            "last_cycle",
            "b_amp_mean_ratio",
            "b_amp_p95_ratio",
            "shock_x_last",
            "n_particles_last",
            "p50_last",
            "p90_last",
            "p99_last",
            "pvtk_last",
    ]
    with path.open("w", encoding="utf-8", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = {k: row.get(k, "") for k in fields}
            writer.writerow(out)


def _plot_comparison(path: Path, rows: list[dict], pinj: float) -> None:
    all_p = []
    for row in rows:
        vals = np.asarray(row["p_mag_last"], dtype=float)
        vals = vals[np.isfinite(vals)]
        vals = vals[vals > 0.0]
        if vals.size > 0:
            all_p.append(vals)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=170)
    ax0, ax1 = axes

    if all_p:
        pmin = float(min(np.min(v) for v in all_p))
        pmax = float(max(np.max(v) for v in all_p))
        if pmin == pmax:
            pmin *= 0.8
            pmax *= 1.2
        bins = np.logspace(np.log10(pmin), np.log10(pmax), 64)
        for row in rows:
            vals = np.asarray(row["p_mag_last"], dtype=float)
            vals = vals[np.isfinite(vals)]
            vals = vals[vals > 0.0]
            if vals.size == 0:
                continue
            hist, edges = np.histogram(vals, bins=bins)
            centers = np.sqrt(edges[:-1] * edges[1:])
            ax0.loglog(centers, np.maximum(hist, 1), lw=1.5, label=row["label"])
    ax0.axvline(pinj, color="k", lw=1.0, ls="--", alpha=0.7)
    ax0.set_xlabel(r"$p$")
    ax0.set_ylabel("Counts per bin")
    ax0.set_title("CR Spectrum (latest output)")
    ax0.grid(alpha=0.25, which="both")
    ax0.legend(fontsize=8, frameon=False)

    labels = [row["label"] for row in rows]
    x = np.arange(len(labels), dtype=float)
    amp = np.array([row["b_amp_p95_ratio"] for row in rows], dtype=float)
    p90 = np.array([row["p90_last"] for row in rows], dtype=float)
    ax1.bar(x - 0.17, amp, width=0.34, label=r"$B_{95}/B_{95,0}$")
    ax1.bar(x + 0.17, p90, width=0.34, label=r"$p_{90}$")
    ax1.set_xticks(x, labels, rotation=15)
    ax1.set_title("Amplification and Spectrum Proxy")
    ax1.grid(alpha=0.25, axis="y")
    ax1.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
            description="Run coarse/fine/AMR pic_parallel_shock decks and compare trends."
    )
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--nproc", type=int, default=8)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--athena-cwd", default="tst/build/src")
    parser.add_argument(
            "--output-root",
            default="tst/.codex/pic_parallel_shock_runs",
            help="Base directory for run logs and benchmark summaries.",
    )
    parser.add_argument("--run-id", default="")
    parser.add_argument("--pinj", type=float, default=3.16227766017 * 3.0)
    args = parser.parse_args()

    if not args.run and not args.run_id:
        parser.error(
                "--run-id is required when --run is not set "
                "(analyze-only mode uses existing artifacts)."
        )

    run_id = args.run_id if args.run_id else f"pic_parallel_shock_bench_{_timestamp()}"
    out_dir = (REPO_ROOT / args.output_root / run_id / "reviews" / "benchmark").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    athena_cwd = (REPO_ROOT / args.athena_cwd).resolve()
    bin_dir = athena_cwd / "bin"
    pvtk_dir = athena_cwd / "pvtk"

    cases = [
            CaseSpec(
                    label="coarse_uniform",
                    deck=str(
                            REPO_ROOT
                            / "inputs/tests/pic_parallel_shock_coarse_uniform.athinput"
                    ),
                    basename=f"{run_id}_coarse",
            ),
            CaseSpec(
                    label="fine_uniform",
                    deck=str(
                            REPO_ROOT
                            / "inputs/tests/pic_parallel_shock_fine_uniform.athinput"
                    ),
                    basename=f"{run_id}_fine",
            ),
            CaseSpec(
                    label="amr_fiducial",
                    deck=str(
                            REPO_ROOT
                            / "inputs/tests/pic_parallel_shock_amr_fiducial.athinput"
                    ),
                    basename=f"{run_id}_amr",
            ),
    ]

    run_rows: list[dict] = []
    if args.run:
        for case in cases:
            cmd = ["./athena", "-i", case.deck, f"job/basename={case.basename}"]
            if args.nproc > 1:
                full_cmd = [args.mpiexec, "-n", str(args.nproc)] + cmd
            else:
                full_cmd = cmd
            t0 = time.time()
            proc = subprocess.run(full_cmd, cwd=str(athena_cwd), check=False)
            runtime = time.time() - t0
            run_rows.append(
                    {
                            "label": case.label,
                            "deck": case.deck,
                            "basename": case.basename,
                            "command": " ".join(full_cmd),
                            "return_code": int(proc.returncode),
                            "runtime_s": float(runtime),
                    }
            )
            if proc.returncode != 0:
                raise RuntimeError(f"case failed: {case.label}")

    summary_rows: list[dict] = []
    try:
        for case in cases:
            summary_rows.append(_compute_case_metrics(bin_dir, pvtk_dir, case))
    except MissingOutputsError as exc:
        print("ERROR: missing benchmark artifacts for analysis.", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        if args.run:
            print(
                    "Hint: this --run invocation did not produce the expected outputs.",
                    file=sys.stderr,
            )
        else:
            print(
                    "Hint: use --run to generate artifacts, or pass --run-id for an "
                    "existing run in the Athena bin/pvtk outputs.",
                    file=sys.stderr,
            )
        return 2

    compare = {
            "run_id": run_id,
            "cases": [],
            "trend_checks": {},
    }
    for row in summary_rows:
        compare["cases"].append(
                {k: v for k, v in row.items() if k != "p_mag_last"}
        )

    by_label = {row["label"]: row for row in summary_rows}
    coarse = by_label["coarse_uniform"]
    fine = by_label["fine_uniform"]
    amr = by_label["amr_fiducial"]
    compare["trend_checks"] = {
            "b_amp_p95_ratio_fine_over_coarse": (
                    fine["b_amp_p95_ratio"] / max(coarse["b_amp_p95_ratio"], 1.0e-30)
            ),
            "b_amp_p95_ratio_amr_over_coarse": (
                    amr["b_amp_p95_ratio"] / max(coarse["b_amp_p95_ratio"], 1.0e-30)
            ),
            "b_amp_p95_ratio_amr_over_fine": (
                    amr["b_amp_p95_ratio"] / max(fine["b_amp_p95_ratio"], 1.0e-30)
            ),
            "p90_amr_over_fine": amr["p90_last"] / max(fine["p90_last"], 1.0e-30),
            "p90_fine_over_coarse": fine["p90_last"] / max(coarse["p90_last"], 1.0e-30),
    }

    summary_json = out_dir / "benchmark_summary.json"
    summary_csv = out_dir / "benchmark_summary.csv"
    run_json = out_dir / "benchmark_runs.json"
    fig_path = out_dir / "pic_parallel_shock_benchmark_compare.png"

    summary_json.write_text(json.dumps(compare, indent=2), encoding="utf-8")
    _write_summary_csv(summary_csv, summary_rows)
    run_json.write_text(json.dumps(run_rows, indent=2), encoding="utf-8")
    _plot_comparison(fig_path, summary_rows, args.pinj)

    print(f"run_id={run_id}")
    print(f"summary_json={summary_json}")
    print(f"summary_csv={summary_csv}")
    print(f"runs_json={run_json}")
    print(f"figure={fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
