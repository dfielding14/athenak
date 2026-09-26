"""Analyze complete uniform 2-D campaign dumps; the right boundary is a reference, not O.

Rates use independent magnetic-flux differences. History EMFs are reconstructed
physical EMFs, not the time-integrated numerical CT flux. Energy boundary flux is
a cell-centered, extrapolated quadrature requiring resolution/cadence convergence.
Flux uses linear half-cell extrapolation to the physical positive-x boundary.
Snapshot energies, energy balance and etaJ2_cc are per unit depth; raw tot-E is
the full-volume history. Heating is a proxy. Null/X/O tests use in-plane B only.
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vis/python"))
from athena_read import hst
from bin_convert import read_binary, read_binary_as_athdf
from reconnection_flux import sheet_flux, sheet_geometry


def parameters(header):
    result, block = {}, ""
    for line in header:
        if line.startswith("<"):
            block = line.strip("<>")
        elif "=" in line:
            key, value = line.split("=", 1)
            result[f"{block}/{key.strip()}"] = value.strip()
    return result


def snapshot(path, null_tolerance):
    raw = read_binary(path)
    active = tuple(raw[f"nx{d}_mb"] for d in (1, 2, 3))
    output = tuple(raw[f"nx{d}_out_mb"] for d in (1, 2, 3))
    expected = {(i, j, 0, 0) for j in range(raw["Nx2"] // active[1])
                for i in range(raw["Nx1"] // active[0])}
    if (raw["Nx3"] != 1 or raw["Nx2"] < 2 or active != output or
            set(map(tuple, raw["mb_logical"])) != expected or raw["n_mbs"] != len(expected)):
        raise ValueError(f"{path}: need full uniform 2-D active-cell dumps")
    p = parameters(raw["header"])
    del raw
    data = read_binary_as_athdf(path, dtype=np.float64)
    x, y = data["x1v"], data["x2v"]
    if not (x[0] < 0 < x[-1] and y[0] < 0 < y[-1]):
        raise ValueError("The fixed X candidate (0,0) must be inside the domain")
    if abs(data["x1f"][0]+data["x1f"][-1])+abs(data["x2f"][0]+data["x2f"][-1]) > 1e-12:
        raise ValueError("History probes must match (0,0): require a domain centered at zero")
    b0, rho_up, di = (float(p[k]) for k in ("problem/b0", "problem/density", "mhd/d_i"))
    bx, by, bz = (data[f"bcc{d}"][0] for d in (1, 2, 3))
    rho, energy = data["dens"][0], data["ener"][0]
    vx, vy, vz = (data[f"mom{d}"][0]/rho for d in (1, 2, 3))
    j = np.searchsorted(y, 0)-1
    weight = -y[j]/(y[j+1]-y[j])
    line = lambda a: (1-weight)*a[j] + weight*a[j+1]
    at = lambda a, xp=0: float(np.interp(xp, x, line(a)))
    bx_y, bx_x = np.gradient(bx, y, x, edge_order=2)
    by_y, by_x = np.gradient(by, y, x, edge_order=2)
    bz_y, bz_x = np.gradient(bz, y, x, edge_order=2)
    axy = 0.5*(bx_x-by_y)
    det_at = lambda xp=0: -at(by_x, xp)*at(bx_y, xp)-at(axy, xp)**2
    # Simple stationary points on the midplane; no search or claim about off-plane O points.
    byline = line(by)
    roots = [float(x[n]-byline[n]*(x[n+1]-x[n])/(byline[n+1]-byline[n]))
             for n in np.where(byline[:-1]*byline[1:] < 0)[0]]
    roots += list(x[1:-1][(byline[1:-1] == 0) & (byline[:-2]*byline[2:] < 0)])
    nulls = [r for r in roots if abs(at(bx, r)) <= null_tolerance*b0]
    valid_x = max(abs(at(bx)), abs(at(by))) <= null_tolerance*b0 and det_at() < 0
    geometry = sheet_geometry(x, y, bx, by, rho, 0, 0, b0, di)
    psi = sheet_flux(x, byline, 0, data["x1f"][-1])
    b2, kinetic = bx*bx+by*by+bz*bz, 0.5*rho*(vx*vx+vy*vy+vz*vz)
    pressure = (float(p["mhd/gamma"])-1)*(energy-kinetic-0.5*b2)
    jx, jy, jz = bz_y, -bz_x, by_x-bx_y
    eta0 = float(p.get("mhd/ohmic_resistivity", 0))
    eta = np.full_like(rho, eta0)
    if p.get("mhd/resistivity_model", "constant") == "current_limited":
        if p.get("mhd/b_rec_method", "constant") != "constant":
            raise ValueError("This pilot analysis supports fixed b_rec only")
        maximum, brec = float(p["mhd/eta_max"]), float(p["mhd/b_rec"])
        q = np.sqrt(jx*jx+jy*jy+jz*jz)*di/brec/np.sqrt(np.maximum(
            rho, float(p.get("mhd/dfloor", np.finfo(np.float32).tiny))))
        qs, es = 1-np.sqrt(eta0/maximum), np.sqrt(eta0*maximum)
        low = q <= qs
        eta[low] = eta0/(1-q[low])
        eta[~low] = maximum-(qs/q[~low])*(maximum-es)
    vdotb = vx*bx+vy*by+vz*bz
    fx = (energy+pressure+0.5*b2)*vx-vdotb*bx+eta*(jy*bz-jz*by)
    fy = (energy+pressure+0.5*b2)*vy-vdotb*by+eta*(jz*bx-jx*bz)
    dx, dy = np.diff(data["x1f"]), np.diff(data["x2f"])
    integrate = lambda a: float(np.sum(a*dy[:, None]*dx[None, :]))
    outward = np.dot(dy, (3*fx[:, -1]-fx[:, -2]-3*fx[:, 0]+fx[:, 1])/2)
    outward += np.dot(dx, (3*fy[-1]-fy[-2]-3*fy[0]+fy[1])/2)
    row = dict(time=data["Time"], time_binary=data["Time"], cycle=data["NumCycles"], psi_ref_minus_X=psi,
               fixed_X_verified=bool(valid_x), bx_at_X=at(bx), by_at_X=at(by),
               hessian_det_X=det_at(), stationary_midplane=len(roots),
               nulls_midplane=len(nulls), O_candidates_midplane=sum(det_at(r)>0 for r in nulls),
               total_energy_snapshot=integrate(energy), magnetic_energy_bcc=integrate(0.5*b2),
               kinetic_energy=integrate(kinetic), internal_energy=integrate(energy-kinetic-0.5*b2),
               density_min=float(rho.min()), pressure_min=float(pressure.min()),
               finite_state=all(np.isfinite(a).all() for a in (rho,energy,pressure,vx,vy,vz,bx,by,bz)),
               X_local_di_cells=di/np.sqrt(at(rho))/max(dx.max(), dy.max()),
               outward_energy_flux_approx=float(outward), divb_max=np.nan,
               delta=geometry[0], delta_local_di=geometry[1], half_length=geometry[2],
               aspect_ratio=geometry[3], rate_normalization=b0*b0/np.sqrt(rho_up),
               depth=float(p["mesh/x3max"])-float(p["mesh/x3min"]))
    divpath = path.with_name(path.name.replace(".state.", ".divb."))
    if divpath.exists():
        div = read_binary(divpath)
        row["divb_max"] = float(max(np.max(np.abs(a)) for a in div["mb_data"]["divb"]))
    return row


def analyze(case, output, null_tolerance):
    rows = sorted((snapshot(p, null_tolerance) for p in case.glob("bin/*.state.*.bin")),
                  key=lambda r: r["time"])
    if not rows:
        raise ValueError(f"No state dumps in {case}/bin")
    t = np.array([r["time"] for r in rows])
    if np.any(np.diff(t) <= 0):
        raise ValueError(f"Duplicate/nonincreasing state times in {case}")
    for kind, names in (("user", ("q_max_edge", "frac_qstar", "frac_q1", "x_q",
                                  "etaJ2_cc", "heat_frac", "x_etaJz", "x_Ez", "ref_etaJz", "ref_Ez")),
                        ("mhd", ("tot-E",))):
        files = list(case.rglob(f"*.{kind}.hst"))
        if len(files) != 1:
            raise ValueError(f"Need one {kind} history in {case}")
        history = hst(str(files[0]))
        # Binary headers store six significant digits; history times retain full precision.
        query = t.copy()
        for n, time in enumerate(t):
            nearest = int(np.argmin(abs(history["time"]-time)))
            rounding = 0.5*10**(np.floor(np.log10(abs(time)))-5) if time else 0
            if abs(history["time"][nearest]-time) <= rounding:
                query[n] = history["time"][nearest]
        if query.min() < history["time"][0] or query.max() > history["time"][-1]:
            raise ValueError(f"{files[0]} does not cover all snapshot times")
        if kind == "user":
            t = query
            for row, time in zip(rows, t):
                row["time"] = float(time)
            history_start, history_end = history["time"][0], history["time"][-1]
            maxima = {name: float(np.max(history[name]))
                      for name in ("q_max_edge", "frac_qstar", "frac_q1")}
        for name in names:
            values = np.interp(query, history["time"], history[name], left=np.nan, right=np.nan)
            for row, value in zip(rows, values):
                row[name] = float(value)
    rate = np.gradient([r["psi_ref_minus_X"] for r in rows], t) if len(t)>1 else [np.nan]
    boundary = np.zeros(len(t))
    flux = np.array([r["outward_energy_flux_approx"] for r in rows])
    boundary[1:] = np.cumsum(0.5*(flux[1:]+flux[:-1])*np.diff(t))
    initial_energy = rows[0]["tot-E"]/rows[0]["depth"]
    for row, derivative, lost in zip(rows, rate, boundary):
        row["flux_rate_raw"] = float(derivative/row["rate_normalization"])
        row["flux_rate_verified_X"] = row["flux_rate_raw"] if row["fixed_X_verified"] else np.nan
        row["physical_EMF_difference"] = row["x_Ez"]-row["ref_Ez"]
        row["resistive_EMF_difference"] = row["x_etaJz"]-row["ref_etaJz"]
        row["rate_minus_physical_EMF"] = row["flux_rate_verified_X"]-row["physical_EMF_difference"]
        row["energy_balance_approx"] = row["tot-E"]/row["depth"]-initial_energy+lost
        row["etaJ2_cc"] /= row["depth"]
    with (output/f"{case.name}.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return dict(case=str(case.resolve()), start=t[0], end=t[-1], snapshots=len(t),
                history_start=history_start, history_end=history_end,
                max_q=maxima["q_max_edge"], max_q_occupancy=maxima["frac_q1"],
                max_qstar_occupancy=maxima["frac_qstar"],
                X_verified_every_dump=all(r["fixed_X_verified"] for r in rows),
                finite_every_dump=all(r["finite_state"] for r in rows),
                density_min=min(r["density_min"] for r in rows),
                pressure_min=min(r["pressure_min"] for r in rows),
                max_abs_divb=max(r["divb_max"] for r in rows),
                energy_balance_approx_final=rows[-1]["energy_balance_approx"],
                energy_balance_baseline_time=t[0],
                flux_reference="Fixed X candidate (0,0) to positive-x boundary at y=0; reference is not O",
                status="Pilot only; no steady-state classification; off-midplane topology not searched")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--null-tolerance", type=float, default=1e-6, help="Null tolerance in units of B0")
    args = parser.parse_args()
    if not np.isfinite(args.null_tolerance) or args.null_tolerance <= 0:
        parser.error("null-tolerance must be finite and positive")
    if len({case.name for case in args.cases}) != len(args.cases):
        parser.error("case directory names must be unique")
    args.output.mkdir(parents=True, exist_ok=True)
    summary = [analyze(case, args.output, args.null_tolerance) for case in args.cases]
    summary = [{k: None if isinstance(v, float) and not np.isfinite(v) else v for k,v in row.items()}
               for row in summary]
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
