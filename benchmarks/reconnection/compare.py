"""Compare analyzed short runs at a common final time; no steady-rate inference.

Run analyze.py first into --analysis. Flux rate is the initial-to-final change,
compared with trapezoidal integration of the complete physical-EMF history.
Field differences restrict fine Bcc by volume averaging on the uniform mesh.
"""
import argparse
import csv
import json
from pathlib import Path

from analyze import hst, np, read_binary_as_athdf


def flux_budget(case, analysis):
    manifest = json.loads((case/"run.json").read_text())
    if manifest["status"] != "complete":
        raise ValueError(f"Incomplete case: {case}")
    with (analysis/f"{case.name}.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    first, last = rows[0], rows[-1]
    start, end = float(first["time"]), float(last["time"])
    if end <= start:
        raise ValueError(f"Need distinct snapshot times: {case}")
    history = hst(str(next(case.rglob("*.user.hst"))))
    if start < history["time"][0] or end > history["time"][-1]:
        raise ValueError(f"History does not cover the comparison interval: {case}")
    t = np.r_[start, history["time"][(history["time"] > start) & (history["time"] < end)], end]
    emf = np.interp(t, history["time"], history["x_Ez"]-history["ref_Ez"])
    physical = float(np.sum(0.5*(emf[1:]+emf[:-1])*np.diff(t))/(end-start))
    flux = (float(last["psi_ref_minus_X"])-float(first["psi_ref_minus_X"])) / (
        (end-start)*float(first["rate_normalization"]))
    return dict(case=case.name, model=manifest["physics"]["model"],
                cells_per_di=manifest["physics"]["cells_per_di"],
                cfl=manifest["arguments"]["cfl"], start=start, end=end,
                mean_flux_rate=flux, mean_physical_EMF_difference=physical,
                mean_rate_discrepancy=flux-physical,
                fixed_X_verified_every_dump=all(r["fixed_X_verified"] == "True" for r in rows))


def field_difference(coarse, fine):
    a, b = [read_binary_as_athdf(sorted(p.glob("bin/*.state.*.bin"))[-1], dtype=np.float64)
            for p in (coarse, fine)]
    if abs(a["Time"]-b["Time"]) > 1e-12:
        raise ValueError("Field comparisons require the same final time")
    result = dict(coarse=coarse.name, fine=fine.name, time=a["Time"])
    for component in (1, 2, 3):
        low, high = a[f"bcc{component}"][0], b[f"bcc{component}"][0]
        ny, nx = low.shape
        if high.shape[0] % ny or high.shape[1] % nx:
            raise ValueError("Fine mesh must be an integer refinement of coarse mesh")
        restricted = high.reshape(ny, high.shape[0]//ny, nx, high.shape[1]//nx).mean(axis=(1, 3))
        difference = np.abs(low-restricted)
        result[f"B{component}_L1"] = float(difference.mean())
        result[f"B{component}_Linf"] = float(difference.max())
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--analysis", required=True, type=Path)
    args = parser.parse_args()
    rates = [flux_budget(case, args.analysis) for case in args.cases]
    if len({(r["start"], r["end"]) for r in rates}) != 1:
        parser.error("All cases must cover the same time interval")
    fields = []
    lookup = {(r["model"], r["cells_per_di"], r["cfl"]): p for r,p in zip(rates,args.cases)}
    for model in ("current_limited", "constant", "ideal"):
        for low, high in ((2, 4), (4, 8)):
            keys = [(model, c, 0.4) for c in (low, high)]
            if all(k in lookup for k in keys):
                fields.append(field_difference(*(lookup[k] for k in keys)))
    keys = [("current_limited", 8, cfl) for cfl in (0.4, 0.2)]
    if all(k in lookup for k in keys):
        fields.append(field_difference(*(lookup[k] for k in keys)))
    result = dict(flux_budgets=rates, Bcc_differences=fields,
                  field_norms="B0=1 code units; L1 is domain mean absolute error, Linf is maximum",
                  interpretation="Startup comparison only; physical EMFs are not numerical CT EMFs")
    (args.analysis/"comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
