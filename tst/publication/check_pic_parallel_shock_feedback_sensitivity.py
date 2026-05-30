#!/usr/bin/env python3
"""Quantify feedback on/off divergence for pic_parallel_shock outputs.

This checker compares two completed runs that differ only by fluid-feedback
toggles (`couple_moments_momentum_to_mhd`, `couple_moments_energy_to_mhd`).
It reports matched-snapshot field deltas and optional threshold gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402
from artifact_lineage import write_companion_manifest  # noqa: E402

_CYCLE_RE = re.compile(r"\.(\d+)\.bin$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fobj:
        for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _input_records(bin_dir: Path, basename: str) -> list[dict[str, str]]:
    rho_map = _map_files(bin_dir, basename, "rho")
    b_map = _map_files(bin_dir, basename, "bmag")
    shared = sorted(set(rho_map) & set(b_map))
    if not shared:
        return []
    paths = [rho_map[shared[-1]], b_map[min(b_map)], b_map[shared[-1]]]
    return [{"path": str(path), "sha256": _sha256(path)} for path in paths]


def _map_files(bin_dir: Path, basename: str, file_id: str) -> dict[int, Path]:
    out: dict[int, Path] = {}
    for path in sorted(bin_dir.glob(f"{basename}.{file_id}.*.bin")):
        match = _CYCLE_RE.search(path.name)
        if match is None:
            continue
        out[int(match.group(1))] = path
    return out


def _slice_xy(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 3:
        return np.asarray(arr[arr.shape[0] // 2, :, :], dtype=float)
    if arr.ndim == 2:
        return np.asarray(arr, dtype=float)
    raise RuntimeError(f"Expected 2D/3D array, got ndim={arr.ndim}")


def _shock_x(rho_xy: np.ndarray, x1v: np.ndarray) -> float:
    rho_x = np.mean(rho_xy, axis=0)
    grad = np.gradient(rho_x, x1v)
    ishock = int(np.argmax(np.abs(grad)))
    return float(x1v[ishock])


def _stats(arr: np.ndarray, xmask: np.ndarray | None) -> dict[str, float]:
    if xmask is None:
        vals = arr.ravel()
    else:
        vals = arr[:, xmask].ravel()
    if vals.size == 0:
        return {"mean": float("nan"), "p95": float("nan"), "p99": float("nan")}
    return {
        "mean": float(np.mean(vals)),
        "p95": float(np.percentile(vals, 95.0)),
        "p99": float(np.percentile(vals, 99.0)),
    }


def _rel_delta(a: float, b: float) -> float:
    if not np.isfinite(a) or abs(a) <= 0.0:
        return float("nan")
    return float((b - a) / a)


def _amp_ratio(first: np.ndarray, last: np.ndarray, percentile: float) -> float:
    p0 = float(np.percentile(first, percentile))
    p1 = float(np.percentile(last, percentile))
    return p1 / max(p0, 1.0e-30)


def _load_last_pair(
    bin_dir: Path, basename: str
) -> tuple[int, float, np.ndarray, np.ndarray, np.ndarray]:
    rho_map = _map_files(bin_dir, basename, "rho")
    b_map = _map_files(bin_dir, basename, "bmag")
    shared = sorted(set(rho_map) & set(b_map))
    if not shared:
        raise RuntimeError(f"No shared rho/bmag cycles for basename={basename}")
    cyc = shared[-1]
    rho_data = bin_convert.read_binary_as_athdf(str(rho_map[cyc]))
    b_data = bin_convert.read_binary_as_athdf(str(b_map[cyc]))
    rho_xy = _slice_xy(np.asarray(rho_data["dens"], dtype=float))
    b_xy = _slice_xy(np.asarray(b_data["bmag"], dtype=float))
    x1v = np.asarray(rho_data["x1v"], dtype=float)
    t = float(rho_data.get("Time", 0.0))
    return cyc, t, rho_xy, b_xy, x1v


def _load_first_bmag(bin_dir: Path, basename: str) -> np.ndarray:
    b_map = _map_files(bin_dir, basename, "bmag")
    if not b_map:
        raise RuntimeError(f"No bmag outputs for basename={basename}")
    first_cycle = min(b_map)
    b_data = bin_convert.read_binary_as_athdf(str(b_map[first_cycle]))
    return _slice_xy(np.asarray(b_data["bmag"], dtype=float))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare feedback-off/on pic_parallel_shock runs with thresholds."
    )
    parser.add_argument("--off-bin-dir", required=True)
    parser.add_argument("--on-bin-dir", required=True)
    parser.add_argument("--off-basename", required=True)
    parser.add_argument("--on-basename", required=True)
    parser.add_argument(
        "--shock-buffer",
        type=float,
        default=0.05,
        help="Exclude [shock_x-buffer, shock_x+buffer] when defining up/down masks.",
    )
    parser.add_argument(
        "--out-json",
        required=True,
        help="Path for metrics JSON.",
    )
    parser.add_argument(
        "--min-b-l2-diff",
        type=float,
        default=0.0,
        help="Fail if b_l2 < threshold.",
    )
    parser.add_argument(
        "--min-rho-l2-diff",
        type=float,
        default=0.0,
        help="Fail if rho_l2 < threshold.",
    )
    parser.add_argument(
        "--min-b-downstream-mean-rel",
        type=float,
        default=0.0,
        help="Fail if |b_downstream_mean_rel| < threshold.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    off_bin = Path(args.off_bin_dir).resolve()
    on_bin = Path(args.on_bin_dir).resolve()

    off_cycle, off_time, off_rho, off_b, x1v = _load_last_pair(off_bin, args.off_basename)
    on_cycle, on_time, on_rho, on_b, _ = _load_last_pair(on_bin, args.on_basename)
    shock_x_off = _shock_x(off_rho, x1v)
    shock_x_on = _shock_x(on_rho, x1v)
    shock_x_mid = 0.5 * (shock_x_off + shock_x_on)
    down_mask = x1v < (shock_x_mid - args.shock_buffer)
    up_mask = x1v > (shock_x_mid + args.shock_buffer)

    drho = on_rho - off_rho
    db = on_b - off_b

    off_b0 = _load_first_bmag(off_bin, args.off_basename)
    on_b0 = _load_first_bmag(on_bin, args.on_basename)

    off_inputs = _input_records(off_bin, args.off_basename)
    on_inputs = _input_records(on_bin, args.on_basename)
    metrics = {
        "artifact_metadata": {
            "evidence_class": "engineering_proxy",
            "not_sun_bai_reproduction": True,
            "qualification_status": "unqualified_feedback_sensitivity_metric",
            "physical_model": "pic_parallel_shock_engineering_scaffold",
            "off_inputs": off_inputs,
            "on_inputs": on_inputs,
        },
        "off": {
            "cycle": int(off_cycle),
            "time": float(off_time),
            "shock_x": float(shock_x_off),
            "rho_all": _stats(off_rho, None),
            "rho_downstream": _stats(off_rho, down_mask),
            "rho_upstream": _stats(off_rho, up_mask),
            "b_all": _stats(off_b, None),
            "b_downstream": _stats(off_b, down_mask),
            "b_upstream": _stats(off_b, up_mask),
            "b_amp_p95_ratio_vs_first": float(_amp_ratio(off_b0, off_b, 95.0)),
            "b_amp_p99_ratio_vs_first": float(_amp_ratio(off_b0, off_b, 99.0)),
        },
        "on": {
            "cycle": int(on_cycle),
            "time": float(on_time),
            "shock_x": float(shock_x_on),
            "rho_all": _stats(on_rho, None),
            "rho_downstream": _stats(on_rho, down_mask),
            "rho_upstream": _stats(on_rho, up_mask),
            "b_all": _stats(on_b, None),
            "b_downstream": _stats(on_b, down_mask),
            "b_upstream": _stats(on_b, up_mask),
            "b_amp_p95_ratio_vs_first": float(_amp_ratio(on_b0, on_b, 95.0)),
            "b_amp_p99_ratio_vs_first": float(_amp_ratio(on_b0, on_b, 99.0)),
        },
        "delta_on_minus_off": {
            "rho_l2": float(np.sqrt(np.mean(drho * drho))),
            "rho_linf": float(np.max(np.abs(drho))),
            "b_l2": float(np.sqrt(np.mean(db * db))),
            "b_linf": float(np.max(np.abs(db))),
            "rho_downstream_mean_rel": _rel_delta(
                float(np.mean(off_rho[:, down_mask]))
                if np.any(down_mask)
                else float("nan"),
                float(np.mean(on_rho[:, down_mask]))
                if np.any(down_mask)
                else float("nan"),
            ),
            "rho_upstream_mean_rel": _rel_delta(
                float(np.mean(off_rho[:, up_mask])) if np.any(up_mask) else float("nan"),
                float(np.mean(on_rho[:, up_mask])) if np.any(up_mask) else float("nan"),
            ),
            "b_downstream_mean_rel": _rel_delta(
                float(np.mean(off_b[:, down_mask]))
                if np.any(down_mask)
                else float("nan"),
                float(np.mean(on_b[:, down_mask])) if np.any(down_mask) else float("nan"),
            ),
            "b_upstream_mean_rel": _rel_delta(
                float(np.mean(off_b[:, up_mask])) if np.any(up_mask) else float("nan"),
                float(np.mean(on_b[:, up_mask])) if np.any(up_mask) else float("nan"),
            ),
        },
    }

    gate = {
        "thresholds": {
            "min_b_l2_diff": float(args.min_b_l2_diff),
            "min_rho_l2_diff": float(args.min_rho_l2_diff),
            "min_b_downstream_mean_rel": float(args.min_b_downstream_mean_rel),
        },
        "checks": {},
        "pass": True,
    }
    b_l2 = metrics["delta_on_minus_off"]["b_l2"]
    rho_l2 = metrics["delta_on_minus_off"]["rho_l2"]
    b_down_rel = abs(metrics["delta_on_minus_off"]["b_downstream_mean_rel"])
    gate["checks"]["b_l2"] = bool(b_l2 >= args.min_b_l2_diff)
    gate["checks"]["rho_l2"] = bool(rho_l2 >= args.min_rho_l2_diff)
    gate["checks"]["b_downstream_mean_rel"] = bool(
        b_down_rel >= args.min_b_downstream_mean_rel
    )
    gate["pass"] = all(gate["checks"].values())
    metrics["gate"] = gate

    out_path = Path(args.out_json).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    write_companion_manifest(
        out_path.with_name(out_path.stem + "_artifact_manifest.json"),
        generator="check_pic_parallel_shock_feedback_sensitivity.py",
        physical_model="pic_parallel_shock_engineering_scaffold",
        inputs=[Path(record["path"]) for record in off_inputs + on_inputs],
        outputs=(out_path,),
    )

    print(f"metrics_json={out_path}")
    print(f"gate_pass={gate['pass']}")
    print(
        "delta_summary: "
        f"b_l2={b_l2:.6e} rho_l2={rho_l2:.6e} "
        f"|b_down_rel|={b_down_rel:.6e}"
    )
    return 0 if gate["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
