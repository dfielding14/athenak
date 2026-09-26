"""X--O flux and energy budgets from uniform 2-D mhd_u_bcc binary dumps.

Supply verified X- and O-point x coordinates on the same sheet midplane.
This script does not track changing magnetic topology. Its cell-centered flux
quadrature is independent of the operator's edge-current history diagnostics.
"""

import argparse
import csv
import sys

import numpy as np

from bin_convert import read_binary, read_binary_as_athdf


def sheet_flux(x, by, xpoint, opoint):
    """Return Az(O)-Az(X) = -integral_X^O By dx along a fixed y line."""
    # Extend half a cell to physical domain faces (e.g. a periodic GEM O point).
    dx0, dx1 = x[1] - x[0], x[-1] - x[-2]
    x = np.r_[x[0] - dx0/2, x, x[-1] + dx1/2]
    by = np.r_[by[0] - (by[1]-by[0])/2, by, by[-1] + (by[-1]-by[-2])/2]
    if not x[0] <= min(xpoint, opoint) <= max(xpoint, opoint) <= x[-1]:
        raise ValueError("X and O points must lie within the domain")
    # Include the requested positions before quadrature, preserving linear By.
    nodes = np.unique(np.r_[x, xpoint, opoint])
    by = np.interp(nodes, x, by)
    x = nodes
    potential = np.r_[0.0, -np.cumsum(0.5 * (by[1:] + by[:-1]) * np.diff(x))]
    return np.interp(opoint, x, potential) - np.interp(xpoint, x, potential)


def sheet_geometry(x, y, bx, by, rho, xpoint, midplane, b0, di):
    """Slope thickness and connected half-maximum current length, as proxies."""
    dbx_dy = np.gradient(bx, y, axis=0)
    jz = np.gradient(by, x, axis=1) - dbx_dy

    def midline(a):
        return np.array([np.interp(midplane, y, column) for column in a.T])

    slope = abs(np.interp(xpoint, x, midline(dbx_dy)))
    delta = b0 / slope if slope > 0 else np.nan
    local_di = di / np.sqrt(np.interp(xpoint, x, midline(rho)))
    current = abs(midline(jz))
    center = int(np.argmin(abs(x-xpoint)))
    threshold = 0.5 * current[center]
    left = right = center
    while left > 0 and current[left-1] >= threshold:
        left -= 1
    while right < len(x)-1 and current[right+1] >= threshold:
        right += 1
    # A layer reaching an open boundary has no measured half-maximum length.
    length = np.nan
    if threshold > 0 and left > 0 and right < len(x)-1:
        xl = x[left-1] + (x[left]-x[left-1]) * (threshold-current[left-1]) / (
            current[left]-current[left-1])
        xr = x[right] + (x[right+1]-x[right]) * (current[right]-threshold) / (
            current[right]-current[right+1])
        length = 0.5 * (xr-xl)
    return delta, delta/local_di, length, delta/length


def measure(path, xpoint, opoint, midplane, b0, di):
    raw = read_binary(path)
    if raw["Nx3"] != 1 or raw["Nx2"] < 2 or np.any(raw["mb_logical"][:, 3] != 0):
        raise ValueError("Flux measurement requires an unrefined, unsliced 2-D dump")
    active = (raw["nx1_mb"], raw["nx2_mb"], raw["nx3_mb"])
    output = (raw["nx1_out_mb"], raw["nx2_out_mb"], raw["nx3_out_mb"])
    blocks = {(i, j, 0) for j in range(raw["Nx2"] // active[1])
              for i in range(raw["Nx1"] // active[0])}
    present = {tuple(p[:3]) for p in raw["mb_logical"]}
    if active != output or present != blocks or len(present) != raw["n_mbs"]:
        raise ValueError("Use a complete, unsharded mhd_u_bcc dump without ghost zones")
    data = read_binary_as_athdf(path, dtype=np.float64)
    x, y = data["x1v"], data["x2v"]
    if not y[0] <= midplane <= y[-1]:
        raise ValueError("The sheet midplane must lie within the output")
    by = np.array([np.interp(midplane, y, column) for column in data["bcc2"][0].T])
    flux = sheet_flux(x, by, xpoint, opoint)
    dv = np.diff(data["x2f"])[:, None] * np.diff(data["x1f"])[None, :]
    total = np.sum(dv * data["ener"][0])
    magnetic = np.sum(dv * sum(data[f"bcc{d}"][0]**2 for d in (1, 2, 3))) / 2
    kinetic = np.sum(dv * sum(data[f"mom{d}"][0]**2 for d in (1, 2, 3))
                     / data["dens"][0]) / 2
    geometry = sheet_geometry(x, y, data["bcc1"][0], data["bcc2"][0],
                              data["dens"][0], xpoint, midplane, b0, di)
    return [data["Time"], data["NumCycles"], flux, total, magnetic,
            kinetic, total - magnetic - kinetic, *geometry]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="Time sequence of full mhd_u_bcc dumps")
    parser.add_argument("--xpoint", type=float, required=True)
    parser.add_argument("--opoint", type=float, required=True)
    parser.add_argument("--midplane", type=float, default=0.0)
    parser.add_argument("--b0", type=float, required=True, help="Upstream reconnecting field")
    parser.add_argument("--rho-up", type=float, required=True, help="Upstream density")
    parser.add_argument("--di", type=float, required=True, help="Ion length at unit density")
    args = parser.parse_args()
    if not np.isfinite([args.b0, args.rho_up, args.di]).all() or min(
            args.b0, args.rho_up, args.di) <= 0:
        parser.error("b0, rho-up, and di must be finite and positive")
    rows = sorted(measure(p, args.xpoint, args.opoint, args.midplane, args.b0, args.di)
                  for p in args.files)
    data = np.asarray(rows)
    if len(data) < 2 or np.any(np.diff(data[:, 0]) <= 0):
        parser.error("At least two distinct, increasing output times are required")
    rate = np.gradient(data[:, 2], data[:, 0]) / (args.b0**2 / np.sqrt(args.rho_up))
    writer = csv.writer(sys.stdout)
    writer.writerow(["time", "cycle", "psi_O_minus_X", "total_energy", "magnetic_energy",
                     "kinetic_energy", "internal_energy", "delta_slope", "delta_local_di",
                     "current_half_length", "delta_over_half_length", "signed_normalized_rate"])
    writer.writerows([*row, value] for row, value in zip(rows, rate))


if __name__ == "__main__":
    main()
