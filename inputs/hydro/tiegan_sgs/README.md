# Tiegan SGS

This project studies fully two-dimensional, compressible isothermal turbulence and
on-the-fly subgrid-scale (SGS) filtering. `mach010_1024.athinput` is the project
fiducial Mach-0.1 run; `mach025_1024.athinput` is the first higher-Mach commissioning
run.

## Two-Dimensional Forcing And Drag

For `Nres = 1024`, the supplied inputs target global mode
`k_drive = sqrt(Nres) = 32`. They use a tile-local parabolic spectrum over modes
`1-3`, peaked at `npeak = 2`, and repeat the realization on `16 x 16` tiles. Thus
the global peak is `tile_nx*npeak = 32`, while the driver only evolves the small
tile-local Fourier mode set.

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

## First Run

`mach025_1024.athinput` uses:

- a `1024 x 1024 x 1` mesh with `512 x 512 x 1` MeshBlocks;
- isothermal sound speed `c_s = 1`;
- target steady-state `v_rms = 0.25` and Mach `0.25`;
- forcing peaked at global mode `32`, so `L_drive = 1/32`;
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
completed five-eddy MPI run reached `v_rms = 0.128946` at `t = 0.625`, with
`v3 = 0`. Statistically steady production measurements should begin after
several drag/box-turnover times or from a saturated restart.

The fiducial `mach010_1024.athinput` uses the same numerical and output setup with
`v_rms = 0.1`, `t_eddy = tcorr = 0.3125`, `drag_rate = 0.1`, `dedt = 0.001`,
`tlim = 3.125`, and an SGS output cadence of `0.015625`.

`mach025_1024_mpi8_5eddy.athinput` is a shorter five-eddy target-Mach-0.25
commissioning run. It uses `256 x 512 x 1` MeshBlocks, forming a `4 x 2`
decomposition intended for eight MPI ranks.

For the isothermal hydro history file,
`v_rms = sqrt(2 * (1-KE + 2-KE + 3-KE) / mass)`. The `3-KE` column should remain
zero for these fully 2D runs.

## Slice Plots

The production inputs write full-resolution primitive binaries and coarsened SGS
binaries at matching times. Plot the latest synchronized snapshot with:

```bash
python inputs/hydro/tiegan_sgs/plot_sgs_slices.py --run-dir RUN_DIRECTORY
```

This writes one three-panel full-resolution `rho`, `vx`, `vy` figure and one
six-panel state-plus-SGS figure for every available coarsening factor under
`RUN_DIRECTORY/plots`.
