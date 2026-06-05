# Tiegan SGS

This project studies fully two-dimensional, compressible isothermal turbulence and
on-the-fly subgrid-scale (SGS) filtering. `mach010_1024.athinput` is the project
fiducial Mach-0.1 run; `mach025_1024.athinput` is the first higher-Mach commissioning
run.

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
- nominal `v_rms = 0.25` and Mach `0.25`;
- forcing scale `L_drive = 0.5`, so `t_eddy = L_drive/v_rms = 2`;
- OU correlation time `tcorr = t_eddy = 2`;
- constant-Edot normalization with `dedt = v_rms^3/L_drive = 0.03125`;
- ten nominal eddy times and 20 SGS outputs per eddy time;
- square filter factors `4, 8, 16, 32, 64, 128`.

Constant Edot is a target-scale choice, not a Mach-number thermostat. In 2D, the
inverse cascade can build large-scale kinetic energy, so the history output should
be used to measure `v_rms` and decide whether the forcing or a large-scale drag
needs retuning before production runs.

The fiducial `mach010_1024.athinput` uses the same numerical and output setup with
`v_rms = 0.1`, `t_eddy = tcorr = 5`, `dedt = 0.002`, `tlim = 50`, and an SGS
output cadence of `0.25`.

For the isothermal hydro history file,
`v_rms = sqrt(2 * (1-KE + 2-KE + 3-KE) / mass)`. The `3-KE` column should remain
zero for these fully 2D runs.
