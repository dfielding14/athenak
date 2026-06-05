<!--
GitHub Pages integration checklist:
1. Copy this file to docs/source/modules/ito_tracers.md on origin/gh-pages.
2. Add modules/ito_tracers to the Modules toctree in docs/source/index.md.
3. Add a row for this page under Physics Modules in docs/source/modules/index.md.
4. Add a short link from docs/source/modules/particles.md and, optionally,
   docs/source/modules/outputs.md.
5. Rebuild with: cd docs && make html.
-->

# Module: Second-Moment Itô Mass-Flux Tracers

## Overview

`particle_type = lagrangian_ito` with `pusher = ito2` implements the Itô-2
method described by Moseley, Teyssier, and Abel in
[arXiv:2604.23041v2](https://arxiv.org/abs/2604.23041). It builds on the
Lagrangian Monte-Carlo tracer infrastructure, but replaces discrete cell-center
jumps with continuous Euler-Maruyama trajectories whose first two displacement
moments match the MC jump kernel.

This implementation intentionally does not implement Itô-3, the piecewise-skew-
uniform sampler, or any third-moment matching. Inputs requesting `pusher = ito3`
or `ito_order` other than `2` fail explicitly.

## Algorithm

For each cell and active coordinate direction, AthenaK obtains the outward MC
jump probabilities from the final Runge-Kutta-weighted mass fluxes:

```text
pR = max( F_R dt / (rho h), 0)
pL = max(-F_L dt / (rho h), 0)
C+ = pR + pL
C- = pR - pL
```

The Itô-2 drift and diffusion coefficient are:

```text
u     = h C- / dt
kappa = h^2 (C+ - C-^2) / (2 dt)
```

The cell-centered `u` and `kappa` fields are communicated across MeshBlocks and
coarse/fine interfaces, then independently interpolated to each particle with
cloud-in-cell interpolation. Each active coordinate is advanced once per fluid
timestep with Euler-Maruyama:

```text
X(n+1) = X(n) + u dt + sqrt(2 kappa dt) xi
```

where each `xi` is an independent bounded uniform draw on
`[-sqrt(3), +sqrt(3)]`, giving zero mean and unit variance. The random draw is
deterministic for a fixed tracer tag, cycle, seed, and coordinate direction.

The stored mass flux is accumulated with the weights that contribute to the
final `rk1`, `rk2`, or `rk3` state, after AMR flux correction. This differs from
a simple arithmetic average of stage fluxes and is required for the Itô moments
to match the final finite-volume update.

## Quick Start

```ini
<time>
evolution  = dynamic
integrator = rk2

<particles>
particle_type   = lagrangian_ito
pusher          = ito2
ito_order       = 2
tracer_kick_pdf = uniform
random_seed     = 12345
track_variables = density, pressure, temperature, v1

<tracer_seed1>
id              = 1
start_time      = 0.0
end_time        = 0.0
cadence         = -1.0
count_per_event = 4096
weight          = mass
region          = all
seed            = 24680
```

Run the uniform-flow example:

```bash
./build/src/athena -i inputs/particles/ito_tracers.athinput -d run_ito2
```

## Supported Configuration

| Capability | Support |
| --- | --- |
| Fluid systems | Non-relativistic Cartesian Hydro and MHD |
| Equation of state | Ideal gas and isothermal |
| Dimensions | 2D and 3D particle configurations |
| Time integrators | Dynamic `rk1`, `rk2`, and `rk3` |
| Mesh | Uniform and AMR |
| Parallelism | Serial and MPI; Kokkos device kernels |
| Boundaries | Periodic in every active dimension |
| Persistence | Restart, particle VTK, and thermodynamic history output |
| Seeding | Shared MC tracer seed schedules and field masks |

The tracer uses one particle step per fluid timestep. To remain compatible with
the existing particle migration path, a realized displacement that spans more
than one MeshBlock in any direction is rejected.

## Validation and Fail-Closed Checks

AthenaK exits with a fatal error when:

- the final mass-flux probabilities are non-finite, have negative variance, or
  have total outward probability greater than one;
- interpolated drift, diffusion, or displacement is invalid;
- a particle would move farther than the existing one-neighbor MeshBlock
  migration path can represent;
- an unsupported integrator, boundary condition, relativistic coordinate
  system, kick distribution, or Itô order is requested.

The regression suite checks:

- the drift and variance of uniform advection against the analytical MC moments;
- both RK2 and RK3 final-stage flux weighting;
- isothermal Hydro and ideal-gas MHD;
- serial restart and explicit Itô-3 rejection;
- MPI migration and restart, AMR coefficient exchange, and preservation of
  continuous subcell positions.

## Validation

Validation on June 5, 2026 passed serial and MPI release builds, the repository
style gate, the thermodynamic-history reader test, the two CPU Itô regression
tests, and the two-rank MPI/AMR regression test. Additional smoke runs passed for
serial and MPI AMR restart, ideal-gas MHD, 3D RK3 transport, explicit Itô-3
rejection, and the pre-existing MC tracer Hydro and MPI/AMR inputs.

## Source Location

| Path | Role |
| --- | --- |
| `src/particles/particles_lagrangian_ito.cpp` | Itô-2 coefficients, AMR communication, CIC interpolation, and Euler-Maruyama push. |
| `src/particles/particles_lagrangian_mc.cpp` | Shared seeding, restart, and post-AMR remapping. |
| `src/hydro/hydro_fluxes.cpp`, `src/mhd/mhd_fluxes.cpp` | Final-RK-weighted mass-flux accumulation. |
| `inputs/particles/ito_tracers*.athinput` | Uniform-flow and AMR smoke inputs. |
| `tst/test_suite/particles/test_particles_ito_*.py` | CPU and MPI regression coverage. |
