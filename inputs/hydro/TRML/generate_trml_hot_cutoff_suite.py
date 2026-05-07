#!/usr/bin/env python3
"""Generate a fixed-grid TRML hot-cutoff cooling suite.

The defaults reproduce the 126-run suite:

  slopes = 0, 2, 4
  log10(T_hot_cut_off/T_cold) = 0.25, 0.50, ..., log10(chi) - 0.25
  log10(xi) = 0.0, 0.5, ..., 3.0

Relative output paths are interpreted from the repository root.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_DIR = Path("inputs/hydro/TRML/generated/hot_cutoff_scan")


@dataclass(frozen=True)
class SuiteConfig:
    out_dir: Path
    chi: float
    mach_rel: float
    slopes: tuple[float, ...]
    xi_log_min: float
    xi_log_max: float
    xi_log_step: float
    tcut_log_min: float
    tcut_log_max: float
    tcut_log_step: float
    t_peak_log_over_tcold: float
    nx: tuple[int, int, int]
    meshblock: tuple[int, int, int]
    domain_size: tuple[float, float, float]
    tlim: float
    overwrite: bool
    dry_run: bool
    particles: bool
    particle_target_count: float
    particle_uniform_by_volume: bool
    particle_pos_init_seed: int
    particle_random_seed: int
    particle_min_radius: float
    particle_inject_at_inflow: bool
    particle_inject_seed: int
    particle_max_inject_per_step: int
    particle_vtk_dt: float
    track_dt: float
    tracked_particles: int

    rho_0: float = 1.0
    pgas_0: float = 1.0
    gamma: float = 5.0 / 3.0
    mu: float = 0.6
    lambda_smooth_width: float = 0.05
    t_floor_over_tcold: float = 0.1
    cfl_number: float = 0.3
    cfl_cool: float = 0.1
    dt_cutoff: float = 1.0e-4
    dt_min_terminate: float = -1.0
    bin_dt: float = 10.0
    slice_dt: float = 1.0


@dataclass(frozen=True)
class SuiteCase:
    filename: str
    basename: str
    slope: float
    lambda_slope_lo: float
    lambda_slope_hi: float
    log10_tcut_over_tcold: float
    t_hot_cut_off: float
    cool_T_min: float
    cool_T_max: float
    log10_xi: float
    xi: float
    T_cold: float
    T_hot: float
    T_peak_over_T_cold: float
    T_peak: float
    velocity: float


def parse_csv_floats(raw: str) -> tuple[float, ...]:
    values = tuple(float(item.strip()) for item in raw.split(",") if item.strip())
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated value")
    return values


def parse_int3(raw: str) -> tuple[int, int, int]:
    parts = tuple(int(item.strip()) for item in raw.split(",") if item.strip())
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated integers")
    return parts


def parse_float3(raw: str) -> tuple[float, float, float]:
    parts = tuple(float(item.strip()) for item in raw.split(",") if item.strip())
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated floats")
    return parts


def resolve_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def fnum(value: float) -> str:
    if value == 0.0:
        value = 0.0
    return f"{value:.16g}"


def token(value: float, digits: int, fixed: bool = False) -> str:
    text = f"{value:.{digits}f}" if fixed else f"{value:.{digits}f}".rstrip("0").rstrip(".")
    if text == "-0":
        text = "0"
    return text.replace("-", "m").replace(".", "p")


def slope_token(value: float) -> str:
    if math.isclose(value, round(value), rel_tol=0.0, abs_tol=1.0e-12):
        return str(int(round(value)))
    return token(value, 3)


def inclusive_values(start: float, stop: float, step: float) -> tuple[float, ...]:
    if step <= 0.0:
        raise ValueError("scan step must be positive")
    if stop < start - 1.0e-12:
        raise ValueError("scan stop must be >= start")
    nstep = int(math.floor((stop - start) / step + 1.0e-9))
    values = [start + i * step for i in range(nstep + 1)]
    if values[-1] < stop - 1.0e-9:
        values.append(stop)
    return tuple(round(value, 12) for value in values)


def x_bounds(length: float) -> tuple[float, float]:
    half = 0.5 * length
    return -half, half


def bool_token(value: bool) -> str:
    return "true" if value else "false"


def validate_config(config: SuiteConfig) -> None:
    if config.chi <= 1.0:
        raise ValueError("chi must be > 1")
    if config.mach_rel < 0.0:
        raise ValueError("Mach_rel must be non-negative")
    if any(slope < 0.0 for slope in config.slopes):
        raise ValueError("cooling-curve suite slopes must be non-negative")
    if config.t_peak_log_over_tcold < 0.0:
        raise ValueError("T_peak must be at or above T_cold for cool_T_min=T_cold")
    if config.t_peak_log_over_tcold > config.tcut_log_min + 1.0e-12:
        raise ValueError("T_peak must be inside every active window; require "
                         "t_peak_log_over_tcold <= tcut_log_min")
    if any(n <= 0 for n in config.nx + config.meshblock):
        raise ValueError("mesh sizes must be positive")
    if any(n % m != 0 for n, m in zip(config.nx, config.meshblock)):
        raise ValueError("each mesh dimension must be divisible by its meshblock size")
    if any(length <= 0.0 for length in config.domain_size):
        raise ValueError("domain sizes must be positive")
    if config.tlim <= 0.0:
        raise ValueError("tlim must be positive")
    if config.particles:
        if config.particle_target_count <= 0.0:
            raise ValueError("particle target count must be positive")
        if config.tracked_particles <= 0:
            raise ValueError("tracked particle count must be positive")
        if config.track_dt <= 0.0:
            raise ValueError("tracked-particle output dt must be positive")
        if config.particle_vtk_dt <= 0.0:
            raise ValueError("particle VTK output dt must be positive")


def build_cases(config: SuiteConfig) -> list[SuiteCase]:
    validate_config(config)

    log10_chi = math.log10(config.chi)
    xi_logs = inclusive_values(config.xi_log_min, config.xi_log_max, config.xi_log_step)
    tcut_logs = inclusive_values(config.tcut_log_min, config.tcut_log_max,
                                 config.tcut_log_step)
    T_hot = config.pgas_0 / config.rho_0
    T_cold = T_hot / config.chi
    T_peak_over_T_cold = 10.0 ** config.t_peak_log_over_tcold
    T_peak = T_cold * T_peak_over_T_cold
    velocity = config.mach_rel * math.sqrt(config.gamma * config.pgas_0 / config.rho_0)

    cases: list[SuiteCase] = []
    chi_part = token(log10_chi, 2, fixed=True)
    mach_part = token(config.mach_rel, 2, fixed=True)

    for slope in config.slopes:
        for tcut_log in tcut_logs:
            t_hot_cut_off = T_cold * 10.0 ** tcut_log
            for xi_log in xi_logs:
                xi = 10.0 ** xi_log
                basename = (
                    f"TRML_chi{chi_part}_M{mach_part}_lam{slope_token(slope)}"
                    f"_TcutTc{token(tcut_log, 2, fixed=True)}"
                    f"_xiLog{token(xi_log, 1, fixed=True)}"
                )
                cases.append(SuiteCase(
                    filename=f"athinput.{basename}",
                    basename=basename,
                    slope=slope,
                    lambda_slope_lo=slope,
                    lambda_slope_hi=-slope,
                    log10_tcut_over_tcold=tcut_log,
                    t_hot_cut_off=t_hot_cut_off,
                    cool_T_min=T_cold,
                    cool_T_max=t_hot_cut_off,
                    log10_xi=xi_log,
                    xi=xi,
                    T_cold=T_cold,
                    T_hot=T_hot,
                    T_peak_over_T_cold=T_peak_over_T_cold,
                    T_peak=T_peak,
                    velocity=velocity,
                ))

    return cases


def render_input(config: SuiteConfig, case: SuiteCase) -> str:
    x1min, x1max = x_bounds(config.domain_size[0])
    x2min, x2max = x_bounds(config.domain_size[1])
    x3min, x3max = x_bounds(config.domain_size[2])
    ft_temp_lo = 1.01 * case.T_cold
    ft_temp_hi = 1.10 * case.T_cold
    particle_block = ""
    particle_output_block = ""
    if config.particles:
        particle_block = f"""
<particles>
particle_type = lagrangian_mc
pusher = lagrangian_mc
target_count = {fnum(config.particle_target_count)}
uniform_by_volume = {bool_token(config.particle_uniform_by_volume)}
pos_init_seed = {config.particle_pos_init_seed}
random_seed = {config.particle_random_seed}
min_radius = {fnum(config.particle_min_radius)}
inject_at_inflow = {bool_token(config.particle_inject_at_inflow)}
inject_seed = {config.particle_inject_seed}
max_inject_per_step = {config.particle_max_inject_per_step}
"""
        particle_output_block = f"""
<output8>
file_type = pvtk
id = prtcl
variable = prtcl_all
dt = {fnum(config.particle_vtk_dt)}

<output9>
file_type = trk
id = trk
variable = prtcl_all
dt = {fnum(config.track_dt)}
nparticles = {config.tracked_particles}
"""

    return f"""# AthenaK Input File: TRML hot-cutoff cooling suite
# Generated by inputs/hydro/TRML/generate_trml_hot_cutoff_suite.py
# basename = {case.basename}
# chi = {fnum(config.chi)}
# Mach_rel = {fnum(config.mach_rel)}
# slope_abs = {fnum(case.slope)}
# log10(T_hot_cut_off/T_cold) = {fnum(case.log10_tcut_over_tcold)}
# xi = {fnum(case.xi)}

<job>
basename = {case.basename}

<mesh>
nghost = 4
nx1 = {config.nx[0]}
x1min = {fnum(x1min)}
x1max = {fnum(x1max)}
ix1_bc = periodic
ox1_bc = periodic
nx2 = {config.nx[1]}
x2min = {fnum(x2min)}
x2max = {fnum(x2max)}
ix2_bc = periodic
ox2_bc = periodic
nx3 = {config.nx[2]}
x3min = {fnum(x3min)}
x3max = {fnum(x3max)}
ix3_bc = outflow
ox3_bc = user

<meshblock>
nx1 = {config.meshblock[0]}
nx2 = {config.meshblock[1]}
nx3 = {config.meshblock[2]}

<time>
evolution = dynamic
integrator = rk3
cfl_number = {fnum(config.cfl_number)}
tlim = {fnum(config.tlim)}
ndiag = 1

<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
gamma = {fnum(config.gamma)}
nscalars = 1
fofc = true
dfloor = 1.0e-4
pfloor = 1.0e-5
{particle_block}

<units>
length_cgs = 1.0
mass_cgs = 1.0
time_cgs = 1.0
mu = {fnum(config.mu)}

<problem>
rho_0 = {fnum(config.rho_0)}
pgas_0 = {fnum(config.pgas_0)}
density_contrast = {fnum(config.chi)}

velocity = {fnum(case.velocity)}
z_interface = 0.0
smoothing_thickness = 0.01

lambda_slope_lo = {fnum(case.lambda_slope_lo)}
lambda_slope_hi = {fnum(case.lambda_slope_hi)}
lambda_smooth_width = {fnum(config.lambda_smooth_width)}
T_peak_over_T_cold = {fnum(case.T_peak_over_T_cold)}
cool_T_min = {fnum(case.cool_T_min)}
cool_T_max = {fnum(case.cool_T_max)}
T_floor = {fnum(config.t_floor_over_tcold)}
heating_mode = off

xi = {fnum(case.xi)}
t_cool_start = 0.0
dt_cutoff = {fnum(config.dt_cutoff)}
dt_min_terminate = {fnum(config.dt_min_terminate)}
cfl_cool = {fnum(config.cfl_cool)}

user_srcs = true
user_dt = true
user_hist = true
ndiag = 100

user_work_in_loop = true
use_frame_tracking = true
t_frame_tracking_start = 1.0
ft_apply_every = 10
ft_mode = pd
ft_position_signal = blend
ft_position_blend = 0.7
ft_max_abs_boost = 0.1
ft_tau_int = 12.0
ft_int_max_abs = 0.0015
ft_int_leak_tau = 8.0
ft_int_unsat_only = true
ft_temp_lo = {fnum(ft_temp_lo)}
ft_temp_hi = {fnum(ft_temp_hi)}
ft_reacquire_enable = true
ft_reacquire_after = 3
ft_reacquire_expand = 1.6
ft_reacquire_max = 6.0
ft_min_confidence = 0.15
ft_conf_smooth_alpha = 0.35
ft_min_cold_mass = 1.0e-8
ft_weight_floor = 1.0e-12
ft_min_support_cells = 16
ft_min_support_volume = 1.0e-6
ft_min_component_weight = 1.0e-12
ft_use_cooling_observer = true
ft_cooling_weight = 0.35
ft_robust_peak_weight = 0.45
ft_shear_weight = 0.10
ft_scalar_weight = 0.10

analysis_on = true
analysis_file = {case.basename}_lightdiag.dat
analysis_every = 10
analysis_clobber = true
analysis_region_zmin = -0.25
analysis_region_zmax = 0.25

<initial_turb>
turb_flag = 1
tcorr = 1.0
dt_turb_update = 1.0
dedt = 1.0
spect_form = 2
expo = 0
nlow = 2
nhigh = 16
sol_fraction = 1.0
rseed = -914
constant_edot = false
z_turb_scale_height = 0.1
z_turb_center = 0.0

<output1>
file_type = log
dcycle = 1

<output2>
file_type = hst
data_format = %12.6e
dcycle = 10

<output3>
file_type = bin
variable = hydro_w
dt = {fnum(config.bin_dt)}

<output4>
file_type = rst
dcycle = 10000
dt_wall = 3600

<output5>
file_type = bin
slice_x1 = 0.0
id = slice_x
variable = hydro_w
dt = {fnum(config.slice_dt)}

<output6>
file_type = bin
slice_x2 = 0.0
id = slice_y
variable = hydro_w
dt = {fnum(config.slice_dt)}

<output7>
file_type = bin
slice_x3 = 0.0
id = slice_z
variable = hydro_w
dt = {fnum(config.slice_dt)}
{particle_output_block}
"""


def manifest_row(config: SuiteConfig, case: SuiteCase) -> dict[str, str]:
    x1min, x1max = x_bounds(config.domain_size[0])
    x2min, x2max = x_bounds(config.domain_size[1])
    x3min, x3max = x_bounds(config.domain_size[2])

    row = {
        "filename": case.filename,
        "basename": case.basename,
        "slope": fnum(case.slope),
        "lambda_slope_lo": fnum(case.lambda_slope_lo),
        "lambda_slope_hi": fnum(case.lambda_slope_hi),
        "lambda_smooth_width": fnum(config.lambda_smooth_width),
        "log10_T_hot_cut_off_over_T_cold": fnum(case.log10_tcut_over_tcold),
        "T_hot_cut_off": fnum(case.t_hot_cut_off),
        "cool_T_min": fnum(case.cool_T_min),
        "cool_T_max": fnum(case.cool_T_max),
        "log10_xi": fnum(case.log10_xi),
        "xi": fnum(case.xi),
        "chi": fnum(config.chi),
        "log10_chi": fnum(math.log10(config.chi)),
        "T_cold": fnum(case.T_cold),
        "T_hot": fnum(case.T_hot),
        "T_peak_over_T_cold": fnum(case.T_peak_over_T_cold),
        "T_peak": fnum(case.T_peak),
        "Mach_rel": fnum(config.mach_rel),
        "velocity": fnum(case.velocity),
        "gamma": fnum(config.gamma),
        "rho_0": fnum(config.rho_0),
        "pgas_0": fnum(config.pgas_0),
        "nx1": str(config.nx[0]),
        "nx2": str(config.nx[1]),
        "nx3": str(config.nx[2]),
        "mb_nx1": str(config.meshblock[0]),
        "mb_nx2": str(config.meshblock[1]),
        "mb_nx3": str(config.meshblock[2]),
        "x1min": fnum(x1min),
        "x1max": fnum(x1max),
        "x2min": fnum(x2min),
        "x2max": fnum(x2max),
        "x3min": fnum(x3min),
        "x3max": fnum(x3max),
        "tlim": fnum(config.tlim),
        "analysis_file": f"{case.basename}_lightdiag.dat",
        "particles": bool_token(config.particles),
        "particle_type": "lagrangian_mc" if config.particles else "",
        "particle_pusher": "lagrangian_mc" if config.particles else "",
        "particle_target_count": fnum(config.particle_target_count) if config.particles else "",
        "particle_uniform_by_volume": (
            bool_token(config.particle_uniform_by_volume) if config.particles else ""
        ),
        "particle_pos_init_seed": str(config.particle_pos_init_seed) if config.particles else "",
        "particle_random_seed": str(config.particle_random_seed) if config.particles else "",
        "particle_min_radius": fnum(config.particle_min_radius) if config.particles else "",
        "particle_inject_at_inflow": (
            bool_token(config.particle_inject_at_inflow) if config.particles else ""
        ),
        "particle_inject_seed": str(config.particle_inject_seed) if config.particles else "",
        "particle_max_inject_per_step": (
            str(config.particle_max_inject_per_step) if config.particles else ""
        ),
        "particle_vtk_dt": fnum(config.particle_vtk_dt) if config.particles else "",
        "track_dt": fnum(config.track_dt) if config.particles else "",
        "tracked_particles": str(config.tracked_particles) if config.particles else "",
    }
    return row


def clear_generated_files(out_dir: Path) -> None:
    for path in out_dir.glob("athinput.TRML_*"):
        if path.is_file():
            path.unlink()
    manifest = out_dir / "manifest.csv"
    if manifest.exists() and manifest.is_file():
        manifest.unlink()


def write_suite(config: SuiteConfig, cases: Iterable[SuiteCase]) -> None:
    cases = list(cases)
    out_dir = resolve_path(config.out_dir)
    if out_dir.exists() and any(out_dir.iterdir()) and not config.overwrite:
        raise SystemExit(f"{out_dir} is not empty; pass --overwrite to replace generated files")
    out_dir.mkdir(parents=True, exist_ok=True)
    if config.overwrite:
        clear_generated_files(out_dir)

    for case in cases:
        (out_dir / case.filename).write_text(render_input(config, case), encoding="utf-8")

    rows = [manifest_row(config, case) for case in cases]
    manifest_path = out_dir / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def print_summary(config: SuiteConfig, cases: list[SuiteCase]) -> None:
    out_dir = resolve_path(config.out_dir)
    first = cases[0]
    last = cases[-1]
    print("TRML hot-cutoff suite")
    print(f"  output directory: {out_dir}")
    print(f"  cases: {len(cases)}")
    print(f"  chi: {fnum(config.chi)}")
    print(f"  Mach_rel: {fnum(config.mach_rel)}")
    print(f"  slopes: {', '.join(fnum(slope) for slope in config.slopes)}")
    if config.particles:
        print(f"  particles: lagrangian_mc target={fnum(config.particle_target_count)} "
              f"tracked={config.tracked_particles} trk_dt={fnum(config.track_dt)}")
    else:
        print("  particles: disabled")
    print(f"  first: {first.filename}  T_hot_cut_off={fnum(first.t_hot_cut_off)}  "
          f"xi={fnum(first.xi)}")
    print(f"  last:  {last.filename}  T_hot_cut_off={fnum(last.t_hot_cut_off)}  "
          f"xi={fnum(last.xi)}")
    if config.dry_run:
        print("  dry run: no files written")
    else:
        print("  wrote athinput files and manifest.csv")


def parse_args() -> SuiteConfig:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--chi", type=float, default=10.0 ** 1.75)
    parser.add_argument("--mach-rel", type=float, default=0.5)
    parser.add_argument("--slopes", type=parse_csv_floats, default=(0.0, 2.0, 4.0))
    parser.add_argument("--xi-log-min", type=float, default=0.0)
    parser.add_argument("--xi-log-max", type=float, default=3.0)
    parser.add_argument("--xi-log-step", type=float, default=0.5)
    parser.add_argument("--tcut-log-min", type=float, default=0.25)
    parser.add_argument("--tcut-log-max", type=float, default=None,
                        help="Default is log10(chi) - 0.25")
    parser.add_argument("--tcut-log-step", type=float, default=0.25)
    parser.add_argument("--t-peak-log-over-tcold", type=float, default=0.25)
    parser.add_argument("--nx", type=parse_int3, default=(128, 128, 256))
    parser.add_argument("--meshblock", type=parse_int3, default=(32, 32, 64))
    parser.add_argument("--domain-size", type=parse_float3, default=(1.0, 1.0, 2.0))
    parser.add_argument("--tlim", type=float, default=100.0)
    parser.add_argument("--no-particles", dest="particles", action="store_false",
                        help="Generate hydro-only inputs without tracer particles")
    parser.set_defaults(particles=True)
    parser.add_argument("--particle-target-count", type=float, default=100000.0)
    parser.add_argument("--particle-uniform-by-volume", action="store_true",
                        help="Sample initial particles uniformly by volume instead of by mass")
    parser.add_argument("--particle-pos-init-seed", type=int, default=280496)
    parser.add_argument("--particle-random-seed", type=int, default=12345)
    parser.add_argument("--particle-min-radius", type=float, default=0.0)
    parser.add_argument("--no-particle-inflow-injection", dest="particle_inject_at_inflow",
                        action="store_false",
                        help="Disable injection of new tracers through physical inflow faces")
    parser.set_defaults(particle_inject_at_inflow=True)
    parser.add_argument("--particle-inject-seed", type=int, default=123456)
    parser.add_argument("--particle-max-inject-per-step", type=int, default=-1)
    parser.add_argument("--particle-vtk-dt", type=float, default=1.0)
    parser.add_argument("--track-dt", type=float, default=0.1)
    parser.add_argument("--tracked-particles", type=int, default=1000)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    tcut_log_max = args.tcut_log_max
    if tcut_log_max is None:
        tcut_log_max = math.log10(args.chi) - 0.25

    return SuiteConfig(
        out_dir=args.out_dir,
        chi=args.chi,
        mach_rel=args.mach_rel,
        slopes=args.slopes,
        xi_log_min=args.xi_log_min,
        xi_log_max=args.xi_log_max,
        xi_log_step=args.xi_log_step,
        tcut_log_min=args.tcut_log_min,
        tcut_log_max=tcut_log_max,
        tcut_log_step=args.tcut_log_step,
        t_peak_log_over_tcold=args.t_peak_log_over_tcold,
        nx=args.nx,
        meshblock=args.meshblock,
        domain_size=args.domain_size,
        tlim=args.tlim,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        particles=args.particles,
        particle_target_count=args.particle_target_count,
        particle_uniform_by_volume=args.particle_uniform_by_volume,
        particle_pos_init_seed=args.particle_pos_init_seed,
        particle_random_seed=args.particle_random_seed,
        particle_min_radius=args.particle_min_radius,
        particle_inject_at_inflow=args.particle_inject_at_inflow,
        particle_inject_seed=args.particle_inject_seed,
        particle_max_inject_per_step=args.particle_max_inject_per_step,
        particle_vtk_dt=args.particle_vtk_dt,
        track_dt=args.track_dt,
        tracked_particles=args.tracked_particles,
    )


def main() -> None:
    config = parse_args()
    cases = build_cases(config)
    if not cases:
        raise SystemExit("scan produced no cases")
    if not config.dry_run:
        write_suite(config, cases)
    print_summary(config, cases)


if __name__ == "__main__":
    main()
