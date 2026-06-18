# Tiegan SGS

This project studies fully two-dimensional, compressible isothermal turbulence and
on-the-fly subgrid-scale (SGS) filtering. The explicit-viscosity science pair is
`mach010_12288_k64_viscous.athinput` and
`mach010_16384_k64_viscous.athinput`. The older `512^2` and `1024^2` inputs are
historical ILES commissioning runs and must not be used as resolved-viscosity DNS.

See [HANDOFF.md](HANDOFF.md) for the scientific rationale, completed pilot results,
known limitations, and the staged plan for moving toward a `16384 x 16384`
production run.

## Resolved-Viscosity Production Contract

The science inputs use AthenaK's uniform isotropic kinematic shear viscosity with
`viscosity = 5.196e-6`. For a narrow forcing annulus centered on Fourier mode
`k_f`, the viscous length `l_nu` uses the same inverse-angular-wavenumber definition
as the earlier 15.76-cell estimate. Its associated wavelength is `lambda_nu`:

```text
eta_est = (2*pi*k_f/L_box)^2 * dedt
l_nu = (viscosity^3/eta_est)^(1/6)
n_nu = L_box/(2*pi*l_nu)
lambda_nu = L_box/n_nu = 2*pi*l_nu
```

For `k_f = 64`, `dedt = 0.001`, and `L_box = 1`, this gives
`eta_est = 161.704`, `l_nu = 9.7660e-4`, `n_nu = 162.97`, and
`lambda_nu = 6.1361e-3`. The viscous length is sampled by `12.0` cells at
`12288^2` and `16.0` cells at `16384^2`. The associated wavelengths are `75.4`
and `100.5` cells, and `n_nu/k_f = 2.546`. The `12288^2` run is the convergence
reference, and the `16384^2` result is the production candidate.

The three desired conditions cannot all hold at `16384^2`: `k_f = 64`,
`l_nu = 16*dx`, and `n_nu/k_f = 16`. The first two imply `n_nu/k_f = 2.546`.
Using `k_f = 96` would reduce it further to `1.698`; obtaining a mode ratio of 16
with a 16-cell `l_nu` would require `k_f` near `10.2`. The production pair keeps
`k_f = 64` because it is the better of the requested `64` and `96` choices, but it
does not provide a broad forward-cascade interval under this estimate.

Use `rsolver = roe`. AthenaK's HLLC implementation is ideal-gas only and rejects an
isothermal EOS. HLLE supports isothermal hydro but is intentionally more diffusive;
it is a robustness fallback, not the baseline for a calculation in which explicit
viscosity should control the small-scale dissipation.

## Two-Dimensional Forcing And Drag

The production inputs use a global parabolic spectrum over the narrow annulus
`63 <= |k| <= 65`, peaked at `npeak = 64`. The sparse-annulus sampler selects
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

## Historical Commissioning Runs

`mach025_1024.athinput` uses:

- a `1024 x 1024 x 1` mesh with `512 x 512 x 1` MeshBlocks;
- isothermal sound speed `c_s = 1`;
- target steady-state `v_rms = 0.25` and Mach `0.25`;
- sparse global forcing over `31 <= |k| <= 33`, peaked at mode `32`, so
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
64 complex modes in the narrow annulus `15 <= |k| <= 17`, and `npeak = 16`.
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

`mach010_12288_k64_viscous.athinput` and
`mach010_16384_k64_viscous.athinput` hold the physical parameters fixed while
increasing the linear resolution by a factor of `4/3`. Both use Mach `0.1`,
`k_f = 64`, ordinary kinematic viscosity `5.196e-6`, Rayleigh drag `0.1`,
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
