#!/usr/bin/env python3
"""Run focused multi-seed audits of the 3D Clebsch initializer."""

from __future__ import annotations

import argparse
import json
import math
import shlex
import subprocess
from pathlib import Path

import numpy as np

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bin_convert_new  # noqa


SEEDS = (123, 321, 777)


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
    return _repo_root() / "build" / "src" / "athena"


def _run_root(athena: Path) -> Path:
    return athena.parent.parent if athena.parent.name == "src" else athena.parent


def _fit_slope(shell_energy: list[float] | tuple[float, ...], nlow: int, nhigh: int) -> float:
    xs = []
    ys = []
    for shell in range(nlow, min(nhigh, len(shell_energy) - 1) + 1):
      value = float(shell_energy[shell])
      if value <= 0.0:
          continue
      xs.append(math.log(float(shell)))
      ys.append(math.log(value))
    if len(xs) < 2:
        return float("nan")
    xmean = sum(xs) / len(xs)
    ymean = sum(ys) / len(ys)
    num = sum((x - xmean) * (y - ymean) for x, y in zip(xs, ys))
    den = sum((x - xmean) ** 2 for x in xs)
    if den <= 0.0:
        return float("nan")
    return num / den


def _leakage(shell_energy: list[float] | tuple[float, ...], nlow: int, nhigh: int) -> float:
    in_band = 0.0
    out_band = 0.0
    for shell, value in enumerate(shell_energy):
        if shell == 0:
            continue
        if nlow <= shell <= nhigh:
            in_band += float(value)
        else:
            out_band += float(value)
    total = in_band + out_band
    return out_band / total if total > 0.0 else 0.0


def _power_law_rel_err(shell_energy: list[float] | tuple[float, ...],
                       slope: float, nlow: int, nhigh: int) -> float:
    if nhigh >= len(shell_energy):
        nhigh = len(shell_energy) - 1
    ref = float(shell_energy[nlow])
    if ref <= 0.0:
        return float("inf")
    target = [0.0] * len(shell_energy)
    for shell in range(nlow, nhigh + 1):
        target[shell] = shell ** (-slope)
    norm = ref / target[nlow]
    max_rel = 0.0
    for shell in range(nlow, nhigh + 1):
        expected = norm * target[shell]
        if expected <= 0.0:
            continue
        rel = abs(float(shell_energy[shell]) - expected) / expected
        max_rel = max(max_rel, rel)
    return max_rel


def _find_hydro_output(case_run_dir: Path, basename: str) -> Path:
    matches = sorted(case_run_dir.rglob(f"{basename}.hydro_w.*.bin"))
    if not matches:
        raise RuntimeError(f"Missing hydro_w output for {basename}")
    return matches[-1]


def _fft_shell_spectrum(path: Path) -> list[float]:
    data = bin_convert_new.read_binary_as_athdf(
        str(path), quantities=["velx", "vely", "velz"], dtype=np.float64
    )
    vx = np.asarray(data["velx"], dtype=np.float64)
    vy = np.asarray(data["vely"], dtype=np.float64)
    vz = np.asarray(data["velz"], dtype=np.float64)
    ntot = float(vx.size)
    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / (ntot * ntot)
    nz, ny, nx = vx.shape
    kz = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    ky = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    kx = np.fft.fftfreq(nx) * nx
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    shell = np.ceil(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3) - 1.0e-12).astype(np.int64)
    return np.bincount(shell.ravel(), weights=energy.ravel()).astype(np.float64).tolist()


def _default_cases() -> list[dict[str, object]]:
    return [
        {
            "name": "small_unthinned",
            "run_nhigh": 6,
            "analysis_nhigh": 6,
            "k_crit": 16.0,
            "alpha": 1.0 / 3.0,
        },
        {
            "name": "small_thinned",
            "run_nhigh": 6,
            "analysis_nhigh": 6,
            "k_crit": 2.0,
            "alpha": 1.0 / 3.0,
        },
    ]


def _production_cases(repo_root: Path) -> list[dict[str, object]]:
    return [
        {
            "name": "prod_thinned",
            "input": repo_root / "inputs" / "hydro" / "scalar_mixing_methods_paper_3d.athinput",
            "run_nhigh": 128,
            "analysis_nhigh": 128,
            "k_crit": 2.0,
            "alpha": 1.0 / 3.0,
        },
    ]


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
        "problem/turb_velocity_method=clebsch",
        f"problem/turb_alpha={case['alpha']}",
        f"problem/turb_rseed={seed}",
        f"problem/turb_nhigh={case['run_nhigh']}",
        f"problem/turb_k_crit={case['k_crit']}",
        "problem/turb_dump_generator_diagnostics=true",
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
    hydro_path = _find_hydro_output(case_run_dir, basename)
    meta = diag["metadata"]
    target = meta["target"]
    spectra = diag["spectra"]
    phi1_shell = [float(x) for x in spectra["phi1_shell_energy"]]
    phi2_shell = [float(x) for x in spectra["phi2_shell_energy"]]
    velocity_shell = _fft_shell_spectrum(hydro_path)
    nlow = int(target["nlow"])
    analysis_nhigh = int(case["analysis_nhigh"])
    phi_slope = float(target["phi_slope"])
    shared_kept_modes = int(diag["catalogs"]["phi"]["kept_modes"])

    return {
        "case": str(case["name"]),
        "seed": seed,
        "run_nhigh": int(case["run_nhigh"]),
        "analysis_nhigh": analysis_nhigh,
        "k_crit": float(case["k_crit"]),
        "alpha": float(target["alpha"]),
        "phi_slope": phi_slope,
        "velocity_slope": float(target["velocity_slope"]),
        "init_wall_seconds": float(meta["init_wall_seconds"]),
        "phi1_kept_modes": shared_kept_modes,
        "phi2_kept_modes": shared_kept_modes,
        "velocity_kept_modes": shared_kept_modes,
        "phi1_shell_energy": phi1_shell,
        "phi2_shell_energy": phi2_shell,
        "velocity_shell_energy": velocity_shell,
        "phi1_exact_rel_err": _power_law_rel_err(phi1_shell, phi_slope, nlow, analysis_nhigh),
        "phi2_exact_rel_err": _power_law_rel_err(phi2_shell, phi_slope, nlow, analysis_nhigh),
        "retained_velocity_slope": _fit_slope(velocity_shell, nlow, analysis_nhigh),
        "retained_velocity_leakage": _leakage(velocity_shell, nlow, analysis_nhigh),
        "stdout_tail": completed.stdout.splitlines()[-12:],
    }


def _aggregate(records: list[dict[str, object]]) -> dict[str, object]:
    by_case: dict[str, list[dict[str, object]]] = {}
    for record in records:
        by_case.setdefault(str(record["case"]), []).append(record)

    summary: dict[str, object] = {}
    for case_name, items in by_case.items():
        summary[case_name] = {
            "run_nhigh": int(items[0]["run_nhigh"]),
            "analysis_nhigh": int(items[0]["analysis_nhigh"]),
            "k_crit": float(items[0]["k_crit"]),
            "alpha": float(items[0]["alpha"]),
            "phi_slope": float(items[0]["phi_slope"]),
            "velocity_slope": float(items[0]["velocity_slope"]),
            "mean_retained_velocity_slope": float(
                sum(float(item["retained_velocity_slope"]) for item in items) / len(items)
            ),
            "mean_retained_velocity_leakage": float(
                sum(float(item["retained_velocity_leakage"]) for item in items) / len(items)
            ),
            "max_phi_power_law_rel_err": float(max(
                max(float(item["phi1_exact_rel_err"]), float(item["phi2_exact_rel_err"]))
                for item in items
            )),
            "mean_init_wall_seconds": float(
                sum(float(item["init_wall_seconds"]) for item in items) / len(items)
            ),
            "mean_velocity_kept_modes": float(
                sum(float(item["velocity_kept_modes"]) for item in items) / len(items)
            ),
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
        default=_repo_root() / "inputs" / "tests" / "scalar_mixing_clebsch_3d_audit.athinput",
        help="Default audit input deck used for the small-profile cases.",
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
    cases = _default_cases() if args.profile == "small" else _production_cases(repo_root)
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
        "seeds": list(SEEDS),
        "records": records,
        "summary": _aggregate(records),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"[audit] wrote {args.output}")


if __name__ == "__main__":
    main()
