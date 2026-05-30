#!/usr/bin/env python3
"""Compute exploratory F01-F09 engineering metrics from one artifact run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402
from artifact_lineage import write_companion_manifest  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.particles.pic_analysis_utils import fit_exponential_growth  # noqa: E402
from scripts.particles.pic_analysis_utils import (  # noqa: E402
    fit_exponential_growth_windowed,
)
from scripts.particles.pic_analysis_utils import polarization_split  # noqa: E402


_SELECTED_SOURCES: list[dict[str, object]] = []


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fobj:
        for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_records(paths: list[Path]) -> list[dict[str, str]]:
    return [{"path": str(path), "sha256": _sha256(path)} for path in paths]


def _selected_source_paths() -> list[Path]:
    return [
        Path(str(record["path"]))
        for source in _SELECTED_SOURCES
        for record in source["files"]
    ]


def _case_dir(run_dir: Path, case_id: str) -> Path:
    return run_dir / "cases" / case_id / "files"


def _case_bin_dir(run_dir: Path, case_id: str) -> Path:
    return _case_dir(run_dir, case_id) / "tst" / "build" / "src" / "bin"


def _case_src_dir(run_dir: Path, case_id: str) -> Path:
    return _case_dir(run_dir, case_id) / "tst" / "build" / "src"


def _select_case_bin_dir(run_dir: Path, candidates: tuple[str, ...]) -> Path:
    for cid in candidates:
        p = _case_bin_dir(run_dir, cid)
        if p.is_dir():
            _record_selected_source(run_dir, "bin", candidates, cid)
            return p
    _record_selected_source(run_dir, "bin", candidates, candidates[0])
    return _case_bin_dir(run_dir, candidates[0])


def _select_case_src_dir(run_dir: Path, candidates: tuple[str, ...]) -> Path:
    for cid in candidates:
        p = _case_src_dir(run_dir, cid)
        if p.is_dir():
            _record_selected_source(run_dir, "src", candidates, cid)
            return p
    _record_selected_source(run_dir, "src", candidates, candidates[0])
    return _case_src_dir(run_dir, candidates[0])


def _record_selected_source(
    run_dir: Path,
    artifact_kind: str,
    candidates: tuple[str, ...],
    selected_case: str,
) -> None:
    case_record = run_dir / "cases" / selected_case / "case.json"
    evidence_class = "legacy_unclassified"
    physical_model = "unknown"
    if case_record.is_file():
        data = json.loads(case_record.read_text(encoding="utf-8"))
        evidence_class = str(data.get("evidence_class", evidence_class))
        physical_model = str(data.get("physical_model", physical_model))
    selected_path = (
        _case_bin_dir(run_dir, selected_case)
        if artifact_kind == "bin"
        else _case_src_dir(run_dir, selected_case)
    )
    files = (
        sorted(path for path in selected_path.rglob("*") if path.is_file())
        if selected_path.is_dir()
        else []
    )
    _SELECTED_SOURCES.append(
        {
            "artifact_kind": artifact_kind,
            "candidate_cases": list(candidates),
            "selected_case": selected_case,
            "selected_case_exists": case_record.is_file(),
            "evidence_class": evidence_class,
            "physical_model": physical_model,
            "selected_path": str(selected_path),
            "files": _file_records(files),
        }
    )


def _latest(pattern: str, root: Path) -> Path | None:
    files = sorted(root.glob(pattern))
    return files[-1] if files else None


def _mode1_complex(dataset: dict, field: str) -> complex:
    values = np.asarray(dataset[field], dtype=float)
    x_mode = np.mean(values, axis=(0, 1))
    fluc = x_mode - np.mean(x_mode)
    x1f = np.asarray(dataset["x1f"], dtype=float)
    x1v = np.asarray(dataset["x1v"], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    return np.sum(fluc * phase) / x_mode.size


def _integrate(dataset: dict, quantity: str) -> float:
    dx1 = np.diff(dataset["x1f"])
    dx2 = np.diff(dataset["x2f"])
    dx3 = np.diff(dataset["x3f"])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(np.asarray(dataset[quantity], dtype=float) * dvol))


def _freq_nonuniform(times: np.ndarray, signal: np.ndarray,
                     fmin: float, fmax: float, nsample: int = 2000) -> float:
    t = np.asarray(times, dtype=float)
    y = np.asarray(signal, dtype=float) - np.mean(signal)
    if t.size < 8:
        return 0.0
    freqs = np.linspace(fmin, fmax, nsample)
    best_f = float(freqs[0])
    best_amp = -1.0
    for freq in freqs:
        omega_t = 2.0 * np.pi * freq * t
        design = np.column_stack((np.sin(omega_t), np.cos(omega_t)))
        coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        amp = float(np.hypot(coeff[0], coeff[1]))
        if amp > best_amp:
            best_amp = amp
            best_f = float(freq)
    return best_f


def _zero_crossings(signal: np.ndarray) -> float:
    shifted = np.asarray(signal, dtype=float) - np.mean(signal)
    signs = np.sign(shifted)
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0.0
    return float(np.count_nonzero(signs[1:] != signs[:-1]))


def _growth_fit_from_series(times: np.ndarray, amp: np.ndarray) -> tuple[float, float]:
    floor = max(1.0e-30, 1.0e-4 * float(np.max(amp)))
    fit = fit_exponential_growth_windowed(
        times, amp, min_points=6, floor=floor, min_growth_factor=1.2
    )
    return float(fit["gamma"]), float(fit["r2"])


def _proxy_stability_fit(
    times: np.ndarray, amp: np.ndarray
) -> tuple[float, float, float]:
    t = np.asarray(times, dtype=float)
    a = np.asarray(amp, dtype=float)
    if t.size != a.size or t.size < 6:
        raise ValueError("proxy stability fit requires at least six samples")
    floor = max(1.0e-16, 1.0e-3 * float(np.max(a)))
    gamma, _, r2 = fit_exponential_growth(
        t, a + floor, float(t[2]), float(t[-3]), floor=1.0e-30
    )
    ratio = float((np.max(a) + floor) / max(np.min(a) + floor, 1.0e-30))
    return float(gamma), float(r2), ratio


def _extended_growth_fit(
    times: np.ndarray, amp: np.ndarray
) -> tuple[float, float, float]:
    t = np.asarray(times, dtype=float)
    a = np.asarray(amp, dtype=float)
    floor = max(1.0e-30, 1.0e-4 * float(np.max(a)))
    fit = fit_exponential_growth_windowed(
        t, a, min_points=6, floor=floor, min_growth_factor=2.0
    )
    return float(fit["gamma"]), float(fit["r2"]), float(fit["growth_factor"])


def _add_entity_metrics(run_dir: Path, out: dict[str, object]) -> None:
    for key, cid in (
        ("F01_mink", "entity_deposit_mink_publication", "entity_deposit_mink"),
        ("F01_reflect", "entity_deposit_reflect_publication", "entity_deposit_reflect"),
    ):
        pass


def _entity_case_metrics(
    bin_dir: Path,
    stem: str,
    ranks: tuple[int, ...],
) -> dict[str, float]:
    np1 = _latest("*_np1.prtcl_d.*.bin", bin_dir)
    if np1 is None:
        np1 = _latest("*.prtcl_d.*.bin", bin_dir)
    if np1 is None:
        return {}
    data1 = bin_convert.read_binary_as_athdf(str(np1))
    pd1 = np.asarray(data1["pdens"], dtype=float)
    out = {
        "np1_sum": float(np.sum(pd1)),
        "np1_max": float(np.max(pd1)),
    }
    for rank in ranks:
        rp = _latest(f"*_np{rank}.prtcl_d.*.bin", bin_dir)
        if rp is None:
            continue
        dr = bin_convert.read_binary_as_athdf(str(rp))
        pdr = np.asarray(dr["pdens"], dtype=float)
        if pdr.shape == pd1.shape:
            out[f"maxdiff_np{rank}"] = float(np.max(np.abs(pdr - pd1)))
    return out


def _collect_entity(run_dir: Path, out: dict[str, object]) -> None:
    b_mink = _select_case_bin_dir(
        run_dir, ("entity_deposit_mink_publication", "entity_deposit_mink")
    )
    b_ref = _select_case_bin_dir(
        run_dir, ("entity_deposit_reflect_publication", "entity_deposit_reflect")
    )
    out["F01_mink"] = _entity_case_metrics(b_mink, "mink", (2, 3, 4))
    out["F01_reflect"] = _entity_case_metrics(b_ref, "reflect", (2, 4))


def _collect_em(run_dir: Path, out: dict[str, object]) -> None:
    src = _select_case_src_dir(
        run_dir, ("em_vacuum_wave_publication", "em_vacuum_wave")
    )
    files = sorted(src.glob("pic_em_vacuum_*np*-errs.dat"))
    pat = re.compile(r"pic_em_vacuum_(?:pub_)?np(\d+)_n(\d+)-errs\.dat")
    vals = {}
    order = {}
    by_np = {}
    for fpath in files:
        m = pat.match(fpath.name)
        if m is None:
            continue
        nproc = int(m.group(1))
        nres = int(m.group(2))
        table = np.loadtxt(fpath)
        row = table if table.ndim == 1 else table[-1, :]
        l1 = float(row[4])
        vals[f"np{nproc}_n{nres}"] = l1
        by_np.setdefault(nproc, {})[nres] = l1

    for nproc, d in by_np.items():
        if len(d) >= 2:
            x = np.log(np.asarray(sorted(d), dtype=float))
            y = np.log(np.asarray([d[r] for r in sorted(d)], dtype=float))
            c = np.polyfit(x, y, 1)
            order[f"np{nproc}_order"] = float(-c[0])
    vals.update(order)
    out["F02_em"] = vals


def _collect_langmuir(run_dir: Path, out: dict[str, object]) -> None:
    bin_dir = _select_case_bin_dir(
        run_dir, ("langmuir_frequency_publication", "langmuir_frequency_proxy")
    )
    metrics = {}
    for rank in (1, 2):
        jx_files = sorted(bin_dir.glob(f"*np{rank}.prtcl_jx.*.bin"))
        jy_files = sorted(bin_dir.glob(f"*np{rank}.prtcl_jy.*.bin"))
        if not jx_files or not jy_files:
            continue
        cycles = {}
        for f in jx_files:
            cycles.setdefault(f.name.split(".")[-2], {})["jx"] = f
        for f in jy_files:
            cycles.setdefault(f.name.split(".")[-2], {})["jy"] = f
        t, jxv, jyv = [], [], []
        for cyc in sorted(cycles):
            if "jx" not in cycles[cyc] or "jy" not in cycles[cyc]:
                continue
            djx = bin_convert.read_binary_as_athdf(str(cycles[cyc]["jx"]))
            djy = bin_convert.read_binary_as_athdf(str(cycles[cyc]["jy"]))
            t.append(float(djx["Time"]))
            jxv.append(_integrate(djx, "prtcl_jx"))
            jyv.append(_integrate(djy, "prtcl_jy"))
        if len(t) < 8:
            continue
        t = np.asarray(t, dtype=float)
        metrics[f"np{rank}_fx"] = _freq_nonuniform(t, np.asarray(jxv), 0.05, 0.40)
        metrics[f"np{rank}_fy"] = _freq_nonuniform(t, np.asarray(jyv), 0.05, 0.40)
    out["F03_langmuir"] = metrics


def _collect_two_stream_weibel(run_dir: Path, out: dict[str, object]) -> None:
    ts_dir = _select_case_bin_dir(
        run_dir, ("two_stream_growth_publication", "two_stream_growth_proxy")
    )
    wb_dir = _select_case_bin_dir(
        run_dir, ("weibel_growth_publication", "weibel_growth_proxy")
    )

    ts = {}
    for rank in (1, 2):
        files = sorted(ts_dir.glob(f"pic_two_stream_pub_np{rank}.mhd_bcc.*.bin"))
        field = "bcc2"
        extended = bool(files)
        if not files:
            files = sorted(
                ts_dir.glob(f"pic_two_stream_growth_proxy_np{rank}.prtcl_rho.*.bin")
            )
            field = "prtcl_rho"
        if not files:
            continue
        t, a = [], []
        for f in files:
            d = bin_convert.read_binary_as_athdf(str(f))
            t.append(float(d["Time"]))
            a.append(abs(_mode1_complex(d, field)))
        if extended:
            gamma, r2, growth_factor = _extended_growth_fit(
                np.asarray(t), np.asarray(a)
            )
            ts[f"np{rank}_growth_factor"] = growth_factor
        else:
            gamma, r2, ratio = _proxy_stability_fit(np.asarray(t), np.asarray(a))
            ts[f"np{rank}_ratio"] = ratio
        ts[f"np{rank}_gamma"] = gamma
        ts[f"np{rank}_r2"] = r2
        ts[f"np{rank}_max"] = float(np.max(a))
    out["F04_twostream"] = ts

    wb = {}
    for rank in (1, 2):
        files = sorted(wb_dir.glob(f"pic_weibel_pub_np{rank}.mhd_bcc.*.bin"))
        field = "bcc2"
        extended = bool(files)
        if not files:
            files = sorted(
                wb_dir.glob(f"pic_weibel_growth_proxy_np{rank}.prtcl_jy.*.bin")
            )
            field = "prtcl_jy"
        if not files:
            continue
        t, a = [], []
        for f in files:
            d = bin_convert.read_binary_as_athdf(str(f))
            t.append(float(d["Time"]))
            a.append(abs(_mode1_complex(d, field)))
        if extended:
            gamma, r2, growth_factor = _extended_growth_fit(
                np.asarray(t), np.asarray(a)
            )
            wb[f"np{rank}_growth_factor"] = growth_factor
        else:
            gamma, r2, ratio = _proxy_stability_fit(np.asarray(t), np.asarray(a))
            wb[f"np{rank}_ratio"] = ratio
        wb[f"np{rank}_gamma"] = gamma
        wb[f"np{rank}_r2"] = r2
        wb[f"np{rank}_max"] = float(np.max(a))
    out["F04_weibel"] = wb


def _collect_bell(run_dir: Path, out: dict[str, object]) -> None:
    bdir = _select_case_bin_dir(
        run_dir, ("bell_growth_publication", "bell_growth_proxy")
    )
    tags = [
        "pic_bell_pub_serial_uncoupled",
        "pic_bell_pub_serial_coupled",
        "pic_bell_pub_mpi2_uncoupled",
        "pic_bell_pub_mpi2_coupled",
        "pic_bell_proxy_serial_uncoupled",
        "pic_bell_proxy_serial_coupled",
        "pic_bell_proxy_mpi2_uncoupled",
        "pic_bell_proxy_mpi2_coupled",
    ]
    metrics = {}
    for tag in tags:
        files = sorted(bdir.glob(f"{tag}.mhd_bcc.*.bin"))
        if not files:
            continue
        t, b = [], []
        for f in files:
            d = bin_convert.read_binary_as_athdf(str(f))
            b2_k = _mode1_complex(d, "bcc2")
            b3_k = _mode1_complex(d, "bcc3")
            bperp = np.sqrt(np.abs(b2_k) * np.abs(b2_k) +
                            np.abs(b3_k) * np.abs(b3_k))
            t.append(float(d["Time"]))
            b.append(float(bperp))
        t = np.asarray(t)
        b = np.asarray(b)
        floor = max(1.0e-30, 1.0e-6 * float(np.max(b)))
        if "coupled" in tag and "uncoupled" not in tag:
            fit = fit_exponential_growth_windowed(
                t, b, min_points=8, floor=floor, min_growth_factor=2.0
            )
            gamma, r2 = float(fit["gamma"]), float(fit["r2"])
        else:
            gamma, _, r2 = fit_exponential_growth(
                t, np.maximum(b, floor), float(t[0]), float(t[-1]), floor=1.0e-30
            )
        key = tag
        metrics[key + "_ratio"] = float(b[-1] / max(b[0], 1.0e-30))
        metrics[key + "_gamma"] = float(gamma)
        metrics[key + "_r2"] = float(r2)
    out["F05_bell"] = metrics


def _collect_mso(run_dir: Path, out: dict[str, object]) -> None:
    bdir = _select_case_bin_dir(
        run_dir,
        (
            "multispecies_backreaction_publication",
            "multispecies_backreaction_oscillation",
        ),
    )
    metrics = {}
    tags = ["uniform", "smr", "amr_publication", "amr_proxy"]
    prefixes = ["pic_mso_pub_", "pic_mso_"]
    for tag in tags:
        for pref in prefixes:
            base = f"{pref}{tag}"
            for rank in (1, 2):
                m2_files = sorted(bdir.glob(f"{base}_np{rank}.mhd_u_m2.*.bin"))
                e_files = sorted(bdir.glob(f"{base}_np{rank}.mhd_u_e.*.bin"))
                if not m2_files or not e_files:
                    continue
                t, mom, ener = [], [], []
                cycles = {}
                for f in m2_files:
                    cycles.setdefault(f.name.split(".")[-2], {})["m2"] = f
                for f in e_files:
                    cycles.setdefault(f.name.split(".")[-2], {})["e"] = f
                for cyc in sorted(cycles):
                    if "m2" not in cycles[cyc] or "e" not in cycles[cyc]:
                        continue
                    dm = bin_convert.read_binary_as_athdf(str(cycles[cyc]["m2"]))
                    de = bin_convert.read_binary_as_athdf(str(cycles[cyc]["e"]))
                    t.append(float(dm["Time"]))
                    mom.append(_integrate(dm, "mom2"))
                    ener.append(_integrate(de, "ener"))
                if len(t) < 8:
                    continue
                t = np.asarray(t)
                mom = np.asarray(mom)
                ener = np.asarray(ener)
                finite = bool(np.all(np.isfinite(mom)) and np.all(np.isfinite(ener)))
                if finite:
                    amp = float(np.max(mom) - np.min(mom))
                    turns = _zero_crossings(mom)
                    drift = float(
                        np.max(np.abs(ener - ener[0])) / max(abs(ener[0]), 1.0)
                    )
                else:
                    amp = -1.0
                    turns = -1.0
                    drift = 1.0e300
                key = f"{tag}_np{rank}".replace("amr_publication", "amr")
                metrics[key + "_amp"] = amp
                metrics[key + "_turns"] = turns
                metrics[key + "_drift"] = drift
                metrics[key + "_finite"] = 1.0 if finite else 0.0
    out["F06_mso"] = metrics


def _collect_crsi(run_dir: Path, out: dict[str, object]) -> None:
    bdir = _select_case_bin_dir(run_dir, ("crsi_deltaf_publication", "crsi_deltaf_proxy"))
    metrics = {}
    tags = [
        ("serial_on", ["pic_crsi_pub_serial_on", "pic_crsi_deltaf_serial_on"]),
        ("mpi2_on", ["pic_crsi_pub_mpi2_on", "pic_crsi_deltaf_mpi2_on"]),
        ("mpi4_on", ["pic_crsi_pub_mpi4_on", "pic_crsi_deltaf_mpi4_on"]),
    ]
    for label, bases in tags:
        files = []
        for base in bases:
            files = sorted(bdir.glob(base + ".mhd_bcc.*.bin"))
            if files:
                break
        if not files:
            continue
        t, by, bz = [], [], []
        for f in files:
            d = bin_convert.read_binary_as_athdf(str(f))
            t.append(float(d["Time"]))
            by.append(_mode1_complex(d, "bcc2"))
            bz.append(_mode1_complex(d, "bcc3"))
        t = np.asarray(t)
        right, left = polarization_split(np.asarray(by), np.asarray(bz))
        ar = np.abs(right)
        al = np.abs(left)
        gr, _ = _growth_fit_from_series(t, ar)
        gl, _ = _growth_fit_from_series(t, al)
        dom = "right" if ar[-1] >= al[-1] else "left"
        metrics[label + "_dom"] = dom
        metrics[label + "_gr"] = float(gr)
        metrics[label + "_gl"] = float(gl)
    out["F07_crsi"] = metrics


def _collect_crpai(run_dir: Path, out: dict[str, object]) -> None:
    bdir = _select_case_bin_dir(
        run_dir, ("crpai_polarization_publication", "crpai_polarization_proxy")
    )
    metrics = {}
    for label in ("prolate", "oblate"):
        for rank in (1, 2, 4):
            files = sorted(bdir.glob(f"pic_crpai_pub_{label}_np{rank}.mhd_bcc.*.bin"))
            if not files:
                files = sorted(bdir.glob(f"pic_crpai_{label}_np{rank}.mhd_bcc.*.bin"))
            if not files:
                continue
            t, by, bz = [], [], []
            for f in files:
                d = bin_convert.read_binary_as_athdf(str(f))
                t.append(float(d["Time"]))
                by.append(_mode1_complex(d, "bcc2"))
                bz.append(_mode1_complex(d, "bcc3"))
            t = np.asarray(t)
            right, left = polarization_split(np.asarray(by), np.asarray(bz))
            ar = np.abs(right)
            al = np.abs(left)
            gr, _ = _growth_fit_from_series(t, ar)
            gl, _ = _growth_fit_from_series(t, al)
            dom = "right" if ar[-1] >= al[-1] else "left"
            asym = float((ar[-1] - al[-1]) / max(ar[-1] + al[-1], 1.0e-30))
            k = f"{label}_np{rank}"
            metrics[k + "_dom"] = dom
            metrics[k + "_asym"] = asym
            metrics[k + "_gr"] = float(gr)
            metrics[k + "_gl"] = float(gl)
    out["F08_crpai"] = metrics


def _collect_aniso(run_dir: Path, out: dict[str, object]) -> None:
    bdir = _select_case_bin_dir(
        run_dir,
        ("expanding_box_anisotropy_publication", "expanding_box_anisotropy_proxy"),
    )
    metrics = {}
    for label in ("expanding", "compressing"):
        base = f"pic_box_{label}_np"
        for rank in (1, 2, 4):
            jx_files = sorted(bdir.glob(f"{base}{rank}.prtcl_jx.*.bin"))
            jy_files = sorted(bdir.glob(f"{base}{rank}.prtcl_jy.*.bin"))
            if not jx_files or not jy_files:
                continue
            cycles = {}
            for f in jx_files:
                cycles.setdefault(f.name.split(".")[-2], {})["jx"] = f
            for f in jy_files:
                cycles.setdefault(f.name.split(".")[-2], {})["jy"] = f
            t, ratio = [], []
            for cyc in sorted(cycles):
                if "jx" not in cycles[cyc] or "jy" not in cycles[cyc]:
                    continue
                djx = bin_convert.read_binary_as_athdf(str(cycles[cyc]["jx"]))
                djy = bin_convert.read_binary_as_athdf(str(cycles[cyc]["jy"]))
                jx = _integrate(djx, "prtcl_jx")
                jy = _integrate(djy, "prtcl_jy")
                if abs(jy) <= 1.0e-12:
                    continue
                t.append(float(djx["Time"]))
                ratio.append(abs(jx) / max(abs(jy), 1.0e-30))
            if len(t) < 3:
                continue
            c = np.polyfit(np.asarray(t), np.asarray(ratio), 1)
            metrics[f"{label}_np{rank}_slope"] = float(c[0])
    out["F09_aniso"] = metrics


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute exploratory engineering F01-F09 metrics summary."
    )
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    out = {}

    _collect_entity(run_dir, out)
    _collect_em(run_dir, out)
    _collect_langmuir(run_dir, out)
    _collect_two_stream_weibel(run_dir, out)
    _collect_bell(run_dir, out)
    _collect_mso(run_dir, out)
    _collect_crsi(run_dir, out)
    _collect_crpai(run_dir, out)
    _collect_aniso(run_dir, out)
    out["artifact_metadata"] = {
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_exploratory_metrics",
        "selected_sources": _SELECTED_SOURCES,
    }

    out_path = run_dir / "metrics_summary.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    metrics_artifact_manifest = run_dir / "metrics_artifact_manifest.json"
    metrics_artifact_manifest.write_text(
        json.dumps(
            {
                "evidence_class": "engineering_proxy",
                "not_sun_bai_reproduction": True,
                "metrics_summary": {
                    "path": str(out_path),
                    "sha256": _sha256(out_path),
                },
                "selected_sources": _SELECTED_SOURCES,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    write_companion_manifest(
        run_dir / "metrics_lineage_manifest.json",
        generator="compute_pic_publication_metrics.py",
        physical_model="exploratory_pic_metric_bundle",
        inputs=_selected_source_paths(),
        outputs=(out_path, metrics_artifact_manifest),
    )
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
