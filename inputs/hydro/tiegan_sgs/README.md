# Tiegan SGS

This project studies fully two-dimensional, compressible isothermal turbulence and
on-the-fly subgrid-scale (SGS) filtering. The explicit-viscosity science pair is
`mach010_12288_n64_viscous.athinput` and
`mach010_16384_n64_viscous.athinput`. The older `512^2` and `1024^2` inputs are
historical ILES commissioning runs and must not be used as resolved-viscosity DNS.

See [HANDOFF.md](HANDOFF.md) for the scientific rationale, completed pilot results,
known limitations, and the staged plan for moving toward a `16384 x 16384`
production run.

## Resolved-Viscosity Production Contract

The science inputs use AthenaK's uniform isotropic kinematic shear viscosity with
`viscosity = 1.316e-7`. AthenaK's `nlow`, `npeak`, and `nhigh` are dimensionless
integer box-mode numbers. For a periodic box of length `L_box`, distinguish this
mode number `n` from physical angular wavenumber `k = 2*pi*n/L_box`:

```text
k_peak = 2*pi*n_peak/L_box
eta_est = k_peak^2 * dedt
ell_visc = (viscosity^3/eta_est)^(1/6) = 1/k_visc
n_visc = L_box*k_visc/(2*pi) = L_box/(2*pi*ell_visc)
lambda_visc = L_box/n_visc = 2*pi*ell_visc
```

For `n_peak = 64`, `dedt = 0.001`, and `L_box = 1`, this gives
`eta_est = 161.704`, `ell_visc = 1.5542e-4`, `n_visc = 1024.03`, and
`lambda_visc = 9.7654e-4`. Thus `n_visc/n_peak = 16.0`. The viscous wavelength is
sampled by `12.0` cells at `12288^2` and `16.0` cells at `16384^2`. The corresponding
`ell_visc = 1/k_visc` sampling is `1.91` and `2.55` cells because it differs from
the wavelength by `2*pi`. The `12288^2` run is the convergence reference, and the
`16384^2` result is the production candidate.

At the same `n_visc = 1024`, `n_peak = 128` would give a ratio of `8`, and
`n_peak = 96` would give `10.67`. The production pair therefore uses `n_peak = 64`,
preserving the requested factor of 16 between forcing and viscous mode while
retaining a factor of 64 between the box and forcing modes.

Use `rsolver = roe`. AthenaK's HLLC implementation is ideal-gas only and rejects an
isothermal EOS. HLLE supports isothermal hydro but is intentionally more diffusive;
it is a robustness fallback, not the baseline for a calculation in which explicit
viscosity should control the small-scale dissipation.

## Two-Dimensional Forcing And Drag

The production inputs use a global parabolic spectrum over the narrow mode annulus
`63 <= |n| <= 65`, peaked at `npeak = 64`. The sparse-annulus sampler selects
64 equal-angle complex modes from one Fourier half-plane; their conjugates supply
the other half-plane. The driver therefore evolves 64 modes without
enumerating the full Cartesian mode volume or periodically repeating a smaller
forcing realization.

The standard large-scale sink is uniform linear Rayleigh drag,
`dv/dt = -drag_rate*v`. It damps every Fourier mode at the same rate, but removes
most energy from the large scales where the inverse cascade accumulates it. The
initial calibration uses

```text
drag_rate = v_rms/L_box
dedt = drag_rate*v_rms^2 = v_rms^3/L_box
```

so drag arrests the inverse cascade near the box scale and the injection-drag
balance estimates the requested global RMS velocity. This is a calibration, not
a Mach-number thermostat.

## Two-Dimensional Contract

- Use `nx3 = 1` and `meshblock/nx3 = 1`.
- Keep the `x3` boundaries periodic. Reflecting boundaries are not needed for a
  degenerate dimension.
- Set `turb_driving/min_kz = turb_driving/max_kz = 0`.
- The turbulence driver only generates in-plane force components on an `nx3 = 1`
  mesh. With zero initial `v3`, the third momentum and kinetic energy remain zero.

## Coarsened SGS Output

Use `file_type = cbin` with `variable = hydro_sgs_2d`. A square top-hat filter of
width `coarsen_factor` is applied in each active dimension. Degenerate dimensions
remain one cell wide.

Each file contains:

| Field | Definition |
| --- | --- |
| `dens` | `bar(rho)` |
| `velx` | `tilde(v_x) = bar(rho v_x) / bar(rho)` |
| `vely` | `tilde(v_y) = bar(rho v_y) / bar(rho)` |
| `tau_xx` | `bar(rho v_x v_x) - bar(rho) tilde(v_x)^2` |
| `tau_xy` | `bar(rho v_x v_y) - bar(rho) tilde(v_x) tilde(v_y)` |
| `tau_yy` | `bar(rho v_y v_y) - bar(rho) tilde(v_y)^2` |

The coarsening factor must divide every active MeshBlock dimension. Keeping every
factor aligned with MeshBlock boundaries makes the local box filters equivalent to
a global non-overlapping square filter. The supplied inputs use powers of two.
Factor 512 therefore requires at least 512 cells in both active MeshBlock
dimensions; it cannot span several smaller MeshBlocks. Coarsened output requires
the full domain without ghost zones, slices, or `gid` selection. Compact hydro SGS
outputs do not support `compute_moments=true` because their fields are already
the final Favre quantities.

Use the default double-precision build for SGS production. The moments and Favre
subtraction use simulation precision; only the final fields are stored as float32.
The SGS kernel reads density and the two momenta once per fine cell per filter
width, accumulating all six moments together without fine-grid moment arrays or
global atomic adds. Large filters are split into 1024-cell partial reductions to
keep many GPU teams busy, then combined. Each width transfers one complete coarse
array to the host, where Favre subtraction precedes float32 serialization.
The six-field product covers the isothermal momentum closure; it does not include
the additional energy-flux moments needed for a non-isothermal closure.

For three-dimensional hydro, use `variable = hydro_sgs_3d` with `file_type = cbin`.
It writes ten final Favre fields: `dens`, `velx`, `vely`, `velz`, `tau_xx`,
`tau_xy`, `tau_xz`, `tau_yy`, `tau_yz`, and `tau_zz`. The definitions above extend
to all three velocity components. This momentum product supports Newtonian ideal
and isothermal hydro; it contains no energy-flux closure. The box filter spans all
three active dimensions and still must fit within each MeshBlock.

For Newtonian ideal-gas MHD, `variable = mhd_sgs` retains the existing 59 raw
moment fields, `mhd_sgs_1` through `mhd_sgs_59`. These are filtered moments from
which SGS stresses, EMFs, and energy-flux terms can be constructed; they are not
final Favre-subtracted labels. Field `mhd_sgs_5` is conserved total energy,
including kinetic and magnetic energy. The field order and duplicate moments are
preserved. Isothermal and relativistic MHD are rejected for this product.
For MHD, `compute_moments=true` retains the existing four raw powers per field
and uses the same optimized coarsening path.

The optimized 3D hydro and MHD coarsening paths accumulate at most eight moments
per GPU team, with 1024 fine cells per partial reduction, then combine the partial
sums and copy the whole coarse output to the host once. This bounds accumulator
register use and avoids allocating ten or 59 full-resolution moment arrays.
The standard moment groups use fixed formulas, avoiding field-selection branches
inside the fine-cell loop.
The 2D six-moment path is unchanged. These SGS reductions use no global
floating-point atomics: HIP's `-munsafe-fp-atomics` is compatible but does not
speed up the reduction itself. It is distinct from `Kokkos_ENABLE_ATOMICS_BYPASS`,
which must remain off for GPU builds.

Measure overhead with the production GPU count, MeshBlocks, output cadence, and
filesystem. This helper recognizes all three SGS products, alternates SGS-on/off
runs, and reports median timings:

```bash
python3 scripts/benchmark_sgs.py build/src/athena case.athinput \
  --launcher 'srun -n 4' --output-dir /scratch/sgs-benchmark time/nlim=1000
```

Choose enough cycles to include several scheduled snapshots beyond the initial
and final dumps. All other outputs and physical parameters stay identical. The
saved logs and `results.json` include timings, output counts, and bytes. Kokkos
profiling regions `SGS2D/load`, `SGS3D/load`, `MHD_SGS/load`, and `cbin/write`
separate calculation/transfer from file packing/writing when a profiling tool is
attached. GPU throughput must be measured on the target hardware; CPU timings do
not establish GPU overhead.

## Historical Commissioning Runs

`mach025_1024.athinput` uses:

- a `1024 x 1024 x 1` mesh with `512 x 512 x 1` MeshBlocks;
- isothermal sound speed `c_s = 1`;
- target steady-state `v_rms = 0.25` and Mach `0.25`;
- sparse global forcing over `31 <= |n| <= 33`, peaked at mode `32`, so
  `L_drive = 1/32`;
- OU correlation time `tcorr = t_eddy = L_drive/v_rms = 0.125`;
- Rayleigh `drag_rate = 0.25` and constant-Edot normalization with
  `dedt = drag_rate*v_rms^2 = 0.015625`;
- ten nominal eddy times and 20 SGS outputs per eddy time;
- full-resolution primitive-state binary output at the SGS cadence;
- square filter factors `4, 8, 16, 32, 64, 128`.

The history output should be used to measure `v_rms` and decide whether `dedt` or
`drag_rate` needs retuning before production runs.

The target Mach number is a steady-state calibration. Starting from rest, five
or ten forcing-scale eddy times are commissioning spin-up runs, not long enough
to reach the target RMS velocity: `t_drag = 1/drag_rate = 4 = 32*t_eddy`. The
earlier tiled five-eddy MPI pipeline test reached `v_rms = 0.128946` at
`t = 0.625`, with `v3 = 0`. Its repeated forcing pattern makes it unsuitable for
SGS science; the sparse-annulus setup requires a fresh spin-up and calibration.
Statistically steady production measurements should begin after several
drag/box-turnover times or from a saturated restart.

The fiducial `mach010_1024.athinput` uses the same numerical and output setup with
`v_rms = 0.1`, `t_eddy = tcorr = 0.3125`, `drag_rate = 0.1`, `dedt = 0.001`,
`tlim = 3.125`, and an SGS output cadence of `0.015625`.

`mach025_1024_mpi8_5eddy.athinput` is a shorter five-eddy target-Mach-0.25
commissioning run. It uses `256 x 512 x 1` MeshBlocks, forming a `4 x 2`
decomposition intended for eight MPI ranks.

`mach025_512_k16_mpi8_steady.athinput` is the first inverse-cascade steady-state
run. It uses a `512 x 512 x 1` mesh split into eight `128 x 256 x 1` MeshBlocks,
64 complex modes in the narrow annulus `15 <= |n| <= 17`, and `npeak = 16`.
With `drag_rate = 0.25` and `dedt = 0.015625`, the inverse-cascade friction-scale
estimate is `k_drag ~ (drag_rate^3/dedt)^(1/2) = 1`, while drag changes the
forcing-scale velocity by only `drag_rate*t_eddy = 1/16` per turnover. The run
lasts 40 forcing-scale turnover times, or five linear energy-relaxation times.
A completed eight-rank calibration reached a final `v_rms = 0.24452`; over the
last time unit its mean and standard deviation were `0.24464` and `0.00012`.
The final velocity spectrum placed `1.31%` of its energy in the box mode,
`21.2%` at `k <= 4`, and `81.8%` below the forcing band, indicating an active
inverse cascade without a box-scale condensate. The third velocity remained zero.

`mach025_512_k16_mpi8_drag0025.athinput` is a controlled comparison that changes
only the Rayleigh drag rate from `0.25` to `0.025`, leaving `dedt`, the random
seed, grid, outputs, and `tlim = 10` unchanged. Its estimated friction scale is
`k_drag ~ 0.0316`, below the box mode, and its linear drag time is `40`. This
run therefore measures transient box-scale accumulation over 40 initial
forcing-scale turnover times rather than a new drag-balanced steady state. A
completed eight-rank run reached a final `v_rms = 0.45581`; over the last time
unit its mean and standard deviation were `0.44648` and `0.00417`. The final
velocity spectrum placed `23.6%` of its energy in the box mode, `62.2%` at
`k <= 4`, and `94.1%` below the forcing band. Thus a tenfold drag reduction at
fixed `dedt` produces a strong box-scale condensate and does not preserve the
target Mach number. The third velocity remained zero.

## Explicit-Viscosity Science Pair

`mach010_12288_n64_viscous.athinput` and
`mach010_16384_n64_viscous.athinput` hold the physical parameters fixed while
increasing the linear resolution by a factor of `4/3`. Both use Mach `0.1`,
`n_peak = 64`, ordinary kinematic viscosity `1.316e-7`, Rayleigh drag `0.1`,
`dedt = 0.001`, and the same random seed. With `tcorr = t_eddy = 0.15625`, each
run lasts 40 forcing-scale turnover times (`tlim = 6.25`). This duration is not
assumed to establish large-scale stationarity from rest; the history and spectra
must demonstrate it.

SGS products are written 20 times per turnover for factors
`4, 8, 16, 32, 64, 128, 256, 512`. Full-resolution primitive states are written
once per turnover, and restarts every four turnovers. A successful run must still
demonstrate the measured enstrophy budget, a viscous dissipation rolloff separated
from the numerical cutoff, and agreement between the `12288^2` and `16384^2`
solutions over their shared resolved range before it is labeled resolved DNS.

For the isothermal hydro history file,
`v_rms = sqrt(2 * (1-KE + 2-KE + 3-KE) / mass)`. The `3-KE` column should remain
zero for these fully 2D runs.
