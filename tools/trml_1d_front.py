#!/usr/bin/env python3
"""Generate simple 1D steady TRML cooling-front profiles.

The model is intentionally minimal:

  * z increases from the cold phase to the hot phase.
  * j = rho*vz is the constant normal mass flux. Condensation from hot to cold
    has j < 0.
  * q = kappa*dT/dz >= 0 is the magnitude of the conductive heat flux from hot
    to cold. The signed conductive flux is F_cond = -q.
  * H = j*c_p*T + F_cond is the thermal enthalpy plus conductive flux.
  * Steady energy balance is dH/dz = -cooling_loss.

The shooting variable is temperature rather than physical position:

    dq/dT = kappa*cooling_loss/q + j*c_p

with q(T_cold)=q(T_hot)=0.  This avoids imposing exact zero-gradient
conditions at finite z; the asymptotic tails are truncated when the profile is
written onto the requested output interval.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import re
import sys
from typing import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import root_scalar


DEFAULT_OUTDIR = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/TRML/simple/front_profiles"
)


@dataclass(frozen=True)
class TRMLParams:
    rho_hot: float
    rho_cold: float
    pres: float
    velocity: float
    t_cool_min: float
    t_cutoff_over_t_cold: float
    beta: float
    gamma: float
    x1min: float
    x1max: float
    x3min: float
    x3max: float

    @property
    def gm1(self) -> float:
        return self.gamma - 1.0

    @property
    def cp(self) -> float:
        # Ideal gas with gas constant R=1 in AthenaK's code units.
        return self.gamma / self.gm1

    @property
    def t_cold(self) -> float:
        # This is the cooling floor used by simple_TRML.cpp.
        return self.pres / self.rho_cold

    @property
    def t_hot(self) -> float:
        return self.pres / self.rho_hot

    @property
    def lx(self) -> float:
        return self.x1max - self.x1min


@dataclass
class FrontSolution:
    mode: str
    j: float
    d_eff: float
    z: np.ndarray
    rho: np.ndarray
    pres: np.ndarray
    temp: np.ndarray
    vz: np.ndarray
    scalar_cold: np.ndarray
    cooling: np.ndarray
    dtdz: np.ndarray
    f_cond: np.ndarray
    h_total: np.ndarray
    kappa: np.ndarray
    residual: np.ndarray
    z_trunc_min: float
    z_trunc_max: float
    q_max: float
    max_abs_residual: float


def parse_value(raw: str) -> str:
    return raw.strip()


def read_athinput(path: Path) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    block = None
    block_re = re.compile(r"^<([^>]+)>$")
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            match = block_re.match(line)
            if match:
                block = match.group(1).strip()
                blocks.setdefault(block, {})
                continue
            if block is None or "=" not in line:
                continue
            key, value = line.split("=", 1)
            blocks[block][key.strip()] = parse_value(value)
    return blocks


def get_float(blocks: dict[str, dict[str, str]], block: str, key: str) -> float:
    try:
        return float(blocks[block][key])
    except KeyError as exc:
        raise KeyError(f"Missing <{block}>/{key}") from exc


def params_from_input(path: Path) -> TRMLParams:
    blocks = read_athinput(path)
    return TRMLParams(
        rho_hot=get_float(blocks, "problem", "rho_hot"),
        rho_cold=get_float(blocks, "problem", "rho_cold"),
        pres=get_float(blocks, "problem", "pres"),
        velocity=get_float(blocks, "problem", "velocity"),
        t_cool_min=get_float(blocks, "problem", "t_cool_min"),
        t_cutoff_over_t_cold=get_float(blocks, "problem", "T_cutoff_over_T_cold"),
        beta=get_float(blocks, "problem", "beta"),
        gamma=get_float(blocks, "hydro", "gamma"),
        x1min=get_float(blocks, "mesh", "x1min"),
        x1max=get_float(blocks, "mesh", "x1max"),
        x3min=get_float(blocks, "mesh", "x3min"),
        x3max=get_float(blocks, "mesh", "x3max"),
    )


def state_from_t(
    temp: np.ndarray | float, j: float, params: TRMLParams, mode: str
) -> tuple[np.ndarray | float, np.ndarray | float]:
    """Return rho, P for a temperature and pressure model."""
    temp_arr = np.asarray(temp)
    if mode == "isobaric":
        pres = np.full_like(temp_arr, params.pres, dtype=float)
        rho = pres / temp_arr
    elif mode == "momentum_flux":
        # Constant normal momentum flux P + rho*vz^2 = P + j^2/rho = M.
        # The hot boundary is the reference state.
        momentum_flux = params.pres + j * j / params.rho_hot
        disc = momentum_flux * momentum_flux - 4.0 * temp_arr * j * j
        if np.any(disc < -1.0e-12):
            raise FloatingPointError(
                "Momentum-flux state became non-real; reduce |j| or check inputs."
            )
        disc = np.maximum(disc, 0.0)
        sqrt_disc = np.sqrt(disc)
        # Pick the branch connected to the requested hot reservoir density.
        # For subsonic hot normal inflow this is the dense/plus branch. If
        # |j| exceeds sqrt(P_hot*rho_hot), the branch swaps; choosing
        # continuously avoids a silent hot-boundary density jump.
        branch_sign = 1.0 if (params.pres - j*j/params.rho_hot) >= 0.0 else -1.0
        rho = (momentum_flux + branch_sign * sqrt_disc) / (2.0 * temp_arr)
        pres = rho * temp_arr
    else:
        raise ValueError(f"unknown pressure mode {mode!r}")

    if np.isscalar(temp):
        return float(np.asarray(rho)), float(np.asarray(pres))
    return rho, pres


def cooling_loss(
    temp: np.ndarray | float, j: float, params: TRMLParams, mode: str
) -> np.ndarray | float:
    """Positive radiative energy loss per unit volume per unit time.

    This matches the continuous cooling law integrated exactly in
    simple_TRML.cpp, including the active temperature interval and the
    T_cold floor.
    """
    temp_arr = np.asarray(temp)
    rho, _ = state_from_t(temp_arr, j, params, mode)
    active = (temp_arr / params.t_cold > 1.0) & (
        temp_arr / params.t_cold <= params.t_cutoff_over_t_cold
    )

    if params.beta != 0.0:
        temp_at_t_cool_min = (
            params.t_cold
            if params.beta > 0.0
            else params.t_cutoff_over_t_cold * params.t_cold
        )
        loss = (
            rho
            / params.gm1
            * temp_at_t_cool_min**params.beta
            / params.t_cool_min
            * temp_arr ** (1.0 - params.beta)
        )
    else:
        eint = rho * temp_arr / params.gm1
        loss = eint / params.t_cool_min

    loss = np.where(active, loss, 0.0)
    if np.isscalar(temp):
        return float(np.asarray(loss))
    return loss


def kappa_eff(
    temp: np.ndarray | float,
    j: float,
    params: TRMLParams,
    mode: str,
    d_eff: float,
) -> np.ndarray | float:
    rho, _ = state_from_t(temp, j, params, mode)
    return rho * params.cp * d_eff


def integrate_q(
    j: float,
    params: TRMLParams,
    mode: str,
    d_eff: float,
    t_start: float,
    t_stop: float,
    n_temp: int | None,
    rtol: float,
    atol: float,
) -> solve_ivp:
    delta_t0 = t_start - params.t_cold
    c0 = cooling_loss(t_start, j, params, mode)
    k0 = kappa_eff(t_start, j, params, mode, d_eff)
    # Local balance near T_cold: d(q^2)/dT ~= 2*kappa*C.
    q0 = math.sqrt(max(2.0 * k0 * c0 * delta_t0, 1.0e-300))

    def rhs(temp: float, y: np.ndarray) -> list[float]:
        q = max(float(y[0]), 1.0e-300)
        c = cooling_loss(temp, j, params, mode)
        kappa = kappa_eff(temp, j, params, mode, d_eff)
        return [kappa * c / q + j * params.cp]

    def q_zero(_temp: float, y: np.ndarray) -> float:
        return float(y[0])

    q_zero.terminal = True  # type: ignore[attr-defined]
    q_zero.direction = -1.0  # type: ignore[attr-defined]

    t_eval = None
    if n_temp is not None:
        # Cluster samples near both asymptotic endpoints.  The physical mapping
        # dz/dT=kappa/q is largest where q -> 0, so a linear temperature grid
        # makes the inferred front span unnecessarily resolution-sensitive.
        u = np.linspace(0.0, 1.0, n_temp)
        smooth = 0.5 * (1.0 - np.cos(np.pi * u))
        t_eval = t_start + (t_stop - t_start) * smooth

    return solve_ivp(
        rhs,
        (t_start, t_stop),
        [q0],
        method="DOP853",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        events=q_zero,
        max_step=(t_stop - t_start) / 100.0,
    )


def find_j(
    params: TRMLParams,
    mode: str,
    d_eff: float,
    t_eps: float,
    rtol: float,
    atol: float,
) -> float:
    t_start = params.t_cold * (1.0 + t_eps)
    t_stop = params.t_hot - (params.t_hot - params.t_cold) * t_eps
    if not (params.t_cold < t_start < t_stop < params.t_hot):
        raise ValueError("invalid temperature endpoints")

    if params.t_hot / params.t_cold <= params.t_cutoff_over_t_cold:
        raise ValueError(
            "The hot phase is inside the cooling-active interval; a zero-gradient "
            "hot asymptote is not available for this simplified front model."
        )

    def residual_for_mag(jmag: float) -> float:
        j = -jmag
        try:
            sol = integrate_q(
                j, params, mode, d_eff, t_start, t_stop, None, rtol, atol
            )
        except (FloatingPointError, ValueError):
            return -np.inf
        if not sol.success and sol.status != 1:
            return np.nan
        if sol.status == 1:
            # q hit zero before the hot phase: |j| is too large.
            return -abs(j * params.cp * (params.t_hot - t_stop)) - 1.0
        q_end = float(sol.y[0, -1])
        q_hot_tail = -j * params.cp * (params.t_hot - t_stop)
        return q_end - q_hot_tail

    # Condensation fronts have j<0.  Search for the magnitude that makes q
    # return to zero at T_hot.
    j_ref = max(params.rho_hot * abs(params.velocity), 1.0e-12)
    j_min = max(1.0e-14 * j_ref, 1.0e-16)
    if mode == "momentum_flux":
        # The algebraic momentum-flux closure becomes non-real once the normal
        # mass flux exceeds the hot-state isothermal momentum limit.
        j_sound = math.sqrt(params.pres * params.rho_hot)
        j_max = min(0.95 * j_sound, 100.0 * j_ref)
    else:
        # The isobaric diagnostic closure has no momentum discriminant.  Permit
        # large |j| so very rapid-cooling / high-diffusivity cases can still be
        # solved and identified as dynamically aggressive.
        j_max = 1.0e5 * j_ref
    samples = np.logspace(math.log10(j_min), math.log10(j_max), 96)

    prev_mag = None
    prev_res = None
    bracket = None
    sampled: list[tuple[float, float]] = []
    for mag in samples:
        res = residual_for_mag(float(mag))
        sampled.append((float(mag), float(res)))
        if not np.isfinite(res):
            continue
        if prev_res is not None and prev_mag is not None:
            if prev_res == 0.0:
                return -prev_mag
            if prev_res * res < 0.0:
                bracket = (prev_mag, float(mag))
                break
        prev_mag = float(mag)
        prev_res = float(res)

    if bracket is None:
        finite = [(mag, res) for mag, res in sampled if np.isfinite(res)]
        if finite:
            best_mag, best_res = min(finite, key=lambda item: abs(item[1]))
            res_min = min(res for _, res in finite)
            res_max = max(res for _, res in finite)
            summary = (
                f"searched |j|=[{finite[0][0]:.6e}, {finite[-1][0]:.6e}], "
                f"residual range=[{res_min:.6e}, {res_max:.6e}], "
                f"closest residual={best_res:.6e} at |j|={best_mag:.6e}"
            )
        else:
            summary = "all residual samples were non-finite"
        raise RuntimeError(
            f"Could not bracket j for mode={mode}; {summary}"
        )

    root = root_scalar(
        residual_for_mag,
        bracket=bracket,
        method="brentq",
        xtol=max(1.0e-13, 1.0e-10 * bracket[1]),
        rtol=1.0e-10,
    )
    if not root.converged:
        raise RuntimeError(f"j root solve did not converge for mode={mode}")
    return -float(root.root)


def solve_front(
    params: TRMLParams,
    mode: str,
    d_eff: float,
    zmin: float,
    zmax: float,
    n_z: int,
    n_temp: int,
    t_eps: float,
    rtol: float,
    atol: float,
) -> FrontSolution:
    j = find_j(params, mode, d_eff, t_eps, rtol, atol)

    t_start = params.t_cold * (1.0 + t_eps)
    t_stop = params.t_hot - (params.t_hot - params.t_cold) * t_eps
    sol = integrate_q(j, params, mode, d_eff, t_start, t_stop, n_temp, rtol, atol)
    if not sol.success:
        raise RuntimeError(f"Final q integration failed for mode={mode}: {sol.message}")

    temp_front = sol.t
    q_front = np.maximum(sol.y[0], 0.0)
    rho_front, pres_front = state_from_t(temp_front, j, params, mode)
    kappa_front = kappa_eff(temp_front, j, params, mode, d_eff)
    dzdt = kappa_front / np.maximum(q_front, 1.0e-300)
    z_front = cumulative_trapezoid(dzdt, temp_front, initial=0.0)

    t_center = math.sqrt(params.t_cold * params.t_hot)
    center_idx = int(np.argmin(np.abs(np.log(temp_front / t_center))))
    z_rel = z_front - z_front[center_idx]

    z = np.linspace(zmin, zmax, n_z)
    temp = np.interp(z, z_rel, temp_front, left=params.t_cold, right=params.t_hot)
    q = np.interp(z, z_rel, q_front, left=0.0, right=0.0)
    rho, pres = state_from_t(temp, j, params, mode)
    kappa = kappa_eff(temp, j, params, mode, d_eff)
    cooling = cooling_loss(temp, j, params, mode)
    dtdz = q / np.maximum(kappa, 1.0e-300)
    f_cond = -q
    h_total = j * params.cp * temp + f_cond
    vz = j / rho
    scalar_cold = np.clip(
        (params.t_hot - temp) / (params.t_hot - params.t_cold), 0.0, 1.0
    )

    residual = np.gradient(h_total, z, edge_order=2) + cooling
    max_abs_residual = float(np.max(np.abs(residual)))

    return FrontSolution(
        mode=mode,
        j=j,
        d_eff=d_eff,
        z=z,
        rho=rho,
        pres=pres,
        temp=temp,
        vz=vz,
        scalar_cold=scalar_cold,
        cooling=cooling,
        dtdz=dtdz,
        f_cond=f_cond,
        h_total=h_total,
        kappa=kappa,
        residual=residual,
        z_trunc_min=float(z_rel[0]),
        z_trunc_max=float(z_rel[-1]),
        q_max=float(np.max(q_front)),
        max_abs_residual=max_abs_residual,
    )


def safe_float_tag(value: float) -> str:
    text = f"{value:g}"
    return text.replace("-", "m").replace("+", "").replace(".", "p")


def default_tag(params: TRMLParams, d_eff_factor: float) -> str:
    chi = params.rho_cold / params.rho_hot
    return (
        f"front_chi{safe_float_tag(chi)}_Mrel{safe_float_tag(params.velocity)}_"
        f"tcool{safe_float_tag(params.t_cool_min)}_Deff{safe_float_tag(d_eff_factor)}"
    )


def write_profile(path: Path, params: TRMLParams, sol: FrontSolution, args: argparse.Namespace) -> None:
    data = np.column_stack(
        [
            sol.z,
            sol.rho,
            sol.pres,
            sol.temp,
            sol.vz,
            sol.scalar_cold,
            sol.cooling,
            sol.dtdz,
            sol.f_cond,
            sol.h_total,
            sol.kappa,
            sol.residual,
        ]
    )
    header = "\n".join(
        [
            "TRML 1D steady cooling-front profile",
            f"input_file = {args.input}",
            f"pressure_mode = {sol.mode}",
            "z convention: cold phase at low z, hot phase at high z",
            "j = rho*vz; j < 0 means hot-to-cold condensation",
            f"j = {sol.j:.17e}",
            f"D_eff = {sol.d_eff:.17e}",
            f"D_eff_factor = {args.d_eff_factor:.17e}",
            f"gamma = {params.gamma:.17e}",
            f"cp = {params.cp:.17e}",
            f"T_cold_floor = {params.t_cold:.17e}",
            f"T_hot = {params.t_hot:.17e}",
            f"rho_hot_input = {params.rho_hot:.17e}",
            f"rho_cold_input = {params.rho_cold:.17e}",
            f"P0_input = {params.pres:.17e}",
            f"cooling_beta = {params.beta:.17e}",
            f"t_cool_min = {params.t_cool_min:.17e}",
            f"T_cutoff_over_T_cold = {params.t_cutoff_over_t_cold:.17e}",
            f"truncated_front_z_min = {sol.z_trunc_min:.17e}",
            f"truncated_front_z_max = {sol.z_trunc_max:.17e}",
            f"q_max = {sol.q_max:.17e}",
            f"max_abs_discrete_residual = {sol.max_abs_residual:.17e}",
            "columns: z rho P T vz_lab scalar_cold cooling_loss dTdz F_cond H_total kappa_eff residual_dHdz_plus_cooling",
        ]
    )
    np.savetxt(path, data, header=header)


def plot_solution(path: Path, params: TRMLParams, sol: FrontSolution) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    ax = axes.ravel()
    ax[0].plot(sol.z, sol.temp)
    ax[0].axhline(params.t_cold, color="k", lw=0.8, ls=":")
    ax[0].axhline(params.t_hot, color="k", lw=0.8, ls=":")
    ax[0].set_ylabel("T")

    ax[1].plot(sol.z, sol.rho)
    ax[1].axhline(params.rho_cold, color="k", lw=0.8, ls=":")
    ax[1].axhline(params.rho_hot, color="k", lw=0.8, ls=":")
    ax[1].set_ylabel(r"$\rho$")

    ax[2].plot(sol.z, sol.pres)
    ax[2].axhline(params.pres, color="k", lw=0.8, ls=":")
    ax[2].set_ylabel("P")

    ax[3].plot(sol.z, sol.vz)
    ax[3].axhline(0.0, color="k", lw=0.8, ls=":")
    ax[3].set_ylabel(r"$v_z$")

    ax[4].semilogy(sol.z, np.maximum(sol.cooling, 1.0e-300))
    ax[4].set_ylabel("cooling loss")
    ax[4].set_xlabel("z")

    ax[5].plot(sol.z, sol.f_cond, label=r"$F_{\rm cond}$")
    ax[5].plot(sol.z, sol.j * params.cp * sol.temp, label=r"$j c_p T$")
    ax[5].plot(sol.z, sol.h_total, label="H")
    ax[5].plot(sol.z, sol.residual, label=r"$dH/dz+\mathcal{L}$", lw=0.8)
    ax[5].axhline(0.0, color="k", lw=0.8, ls=":")
    ax[5].set_ylabel("flux / residual")
    ax[5].set_xlabel("z")
    ax[5].legend(fontsize=8)

    fig.suptitle(
        f"{sol.mode}: j={sol.j:.4e}, D_eff={sol.d_eff:.4e}, "
        f"max|res|={sol.max_abs_residual:.3e}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_compare(path: Path, params: TRMLParams, solutions: list[FrontSolution]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    ax = axes.ravel()
    for sol in solutions:
        label = f"{sol.mode}, j={sol.j:.3e}"
        ax[0].plot(sol.z, sol.temp, label=label)
        ax[1].plot(sol.z, sol.rho, label=label)
        ax[2].plot(sol.z, sol.pres, label=label)
        ax[3].plot(sol.z, sol.vz, label=label)
        ax[4].semilogy(sol.z, np.maximum(sol.cooling, 1.0e-300), label=label)
        ax[5].plot(sol.z, sol.h_total, label=label)

    ax[0].set_ylabel("T")
    ax[1].set_ylabel(r"$\rho$")
    ax[2].set_ylabel("P")
    ax[3].set_ylabel(r"$v_z$")
    ax[4].set_ylabel("cooling loss")
    ax[5].set_ylabel("H")
    ax[4].set_xlabel("z")
    ax[5].set_xlabel("z")
    for axis in ax:
        axis.axvline(0.0, color="k", lw=0.5, ls=":")
    ax[0].legend(fontsize=8)
    fig.suptitle(
        f"TRML 1D front comparison, T_cold={params.t_cold:.3e}, T_hot={params.t_hot:.3e}"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="AthenaK athinput file")
    parser.add_argument(
        "--pressure-mode",
        choices=("isobaric", "momentum_flux", "both"),
        default="both",
    )
    parser.add_argument(
        "--D-eff-factor",
        dest="d_eff_factor",
        type=float,
        default=0.1,
        help="D_eff = factor * <problem>/velocity * (x1max-x1min)",
    )
    parser.add_argument("--zmin", type=float, default=None)
    parser.add_argument("--zmax", type=float, default=None)
    parser.add_argument("--n", type=int, default=2048, help="number of output z points")
    parser.add_argument(
        "--n-temp",
        type=int,
        default=4096,
        help="number of temperature-space integration samples",
    )
    parser.add_argument(
        "--temperature-epsilon",
        type=float,
        default=1.0e-6,
        help="fractional truncation away from asymptotic temperatures",
    )
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument("--atol", type=float, default=1.0e-11)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--tag", default=None)
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    params = params_from_input(args.input)

    modes = (
        ["isobaric", "momentum_flux"]
        if args.pressure_mode == "both"
        else [args.pressure_mode]
    )
    zmin = params.x3min if args.zmin is None else args.zmin
    zmax = params.x3max if args.zmax is None else args.zmax
    d_eff = args.d_eff_factor * abs(params.velocity) * params.lx
    tag = args.tag or default_tag(params, args.d_eff_factor)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    solutions: list[FrontSolution] = []
    failures = 0
    for mode in modes:
        try:
            sol = solve_front(
                params=params,
                mode=mode,
                d_eff=d_eff,
                zmin=zmin,
                zmax=zmax,
                n_z=args.n,
                n_temp=args.n_temp,
                t_eps=args.temperature_epsilon,
                rtol=args.rtol,
                atol=args.atol,
            )
        except Exception as exc:
            failures += 1
            print(f"{mode}: FAILED: {exc}", file=sys.stderr)
            continue
        solutions.append(sol)
        txt_path = args.output_dir / f"{tag}_{mode}.txt"
        png_path = args.output_dir / f"{tag}_{mode}.png"
        write_profile(txt_path, params, sol, args)
        if not args.no_plot:
            plot_solution(png_path, params, sol)
        print(
            f"{mode}: wrote {txt_path}  j={sol.j:.8e}  "
            f"front_span=[{sol.z_trunc_min:.4g},{sol.z_trunc_max:.4g}]  "
            f"max|res|={sol.max_abs_residual:.3e}"
        )
        if not args.no_plot:
            print(f"{mode}: wrote {png_path}")

    if len(solutions) > 1 and not args.no_plot:
        compare_path = args.output_dir / f"{tag}_compare.png"
        plot_compare(compare_path, params, solutions)
        print(f"compare: wrote {compare_path}")

    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
