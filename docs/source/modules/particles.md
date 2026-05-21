# Module: Particles

The particle module advances Lagrangian particles on top of a hydro or MHD mesh.  This
branch adds massive dust particles with gas drag and optional backreaction on the fluid.
Drag-coupled particles currently require exactly one host fluid module, either `<hydro>`
or `<mhd>`; combined hydro+MHD two-fluid configurations are rejected until the
backreaction update is generalized for both primitive states.

## Source Location

| File | Purpose |
| --- | --- |
| `src/particles/particles.hpp` | Particle data layout, pusher and drag model enums |
| `src/particles/particles.cpp` | Particle construction, input parsing, AMR remapping |
| `src/particles/particles_pushers.cpp` | Position pushers |
| `src/particles/particles_drag.cpp` | Drag force, backreaction, pgen drag hook dispatch |
| `src/particles/particles_tasks.cpp` | Particle task registration |
| `src/bvals/bvals_part.cpp` | MeshBlock and MPI particle exchange |

## Particle Types

`<particles>/particle_type = cosmic_ray` keeps the historical passive particle layout.
`<particles>/particle_type = dust` allocates position, velocity, and mass slots for
drag-coupled massive particles:

```cpp
PGID = 0, PTAG = 1
IPX = 0, IPVX = 1
IPY = 2, IPVY = 3
IPZ = 4, IPVZ = 5
IPM = 6
```

`PGID` is the global MeshBlock ID that currently owns the particle.  `PTAG` is the
global particle tag.  `IPM` is the particle mass used for drag backreaction.

## Standard Drag

Enable the standard stopping-time law with:

```ini
<particles>
particle_type = dust
pusher = drag
ppc = 1.0

<drag_particles>
enabled = true
model = stopping_time
stopping_time = 0.1
cfl_drag = 0.05
backreaction = true
include_energy = false
interpolation = host_cell
deposition = host_cell
particle_mass = 1.0
orbital_terms = false
```

The standard force is

```{math}
\frac{d\mathbf{v}_p}{dt} = \frac{\mathbf{u}_g-\mathbf{v}_p}{t_s}.
```

`stopping_time` is in code time units.  `cfl_drag` caps the particle timestep by
`cfl_drag * stopping_time`; the default is deliberately conservative because the drag
exchange is operator split from the fluid update.  `particle_mass` is used only when a
pgen leaves `IPM <= 0`; production pgens should initialize each particle mass explicitly.

The implementation samples the host fluid primitive velocity, applies the exact
single-particle stopping-time update for the chosen step, and deposits the equal and
opposite momentum impulse into the fluid when `backreaction = true`:

```{math}
\Delta \mathbf{M}_g = -m_p\,\Delta \mathbf{v}_{p,\mathrm{drag}}.
```

The sign convention is that a positive particle velocity change gives the gas a negative
momentum impulse.  The deposited conserved momentum is divided by the cell volume because
AthenaK stores cell-averaged conserved variables.  If `include_energy = true` and the gas
EOS is ideal, the opposite of the drag change in particle kinetic energy is deposited
into the gas total energy.  For isothermal fluids, `include_energy` is ignored.

The drag exchange is operator split in the `after_timeintegrator` task list.  After the
exchange, hydro or MHD primitives are refreshed so the next step sees the updated fluid.

## Coupling Stencils

Two coupling stencils are available independently for cell-to-particle interpolation and
particle-to-cell backreaction deposition:

| Value | Meaning |
| --- | --- |
| `host_cell` or `nearest` | Use the particle's current host cell only. This is the default. |
| `cloud_in_cell` or `cic` | Use the 2, 4, or 8 nearest active cells in 1D, 2D, or 3D geometry. |

Example:

```ini
<drag_particles>
interpolation = cloud_in_cell
deposition = cloud_in_cell
```

The current `cloud_in_cell` implementation is block-local.  At a MeshBlock boundary the
stencil is clipped to active cells in the owning block and renormalized, rather than
communicating a cross-block deposition stencil.  This is still a true weighted
interpolation/deposition option inside each block, but it is not yet a full mesh-wide CIC
scheme.  Use `host_cell` when bitwise-minimal behavior is more important than smoothing
the coupling stencil.

## Orbital Terms

For local streaming-instability tests without the built-in shearing-box module, the drag
pusher can add particle epicyclic terms:

```ini
<drag_particles>
orbital_terms = true
omega0 = 1.0
qshear = 1.5
```

The particle update adds

```{math}
\dot v_x = 2\Omega v_y,\qquad
\dot v_y = -(2-q)\Omega v_x
```

before applying drag.  `omega0` has units of inverse code time and `qshear` is
dimensionless.  These orbital accelerations are external forces and are not backreacted
onto the gas.

## Pgen Drag Hook

Set `model = user` to let the problem generator define the force law:

```ini
<drag_particles>
enabled = true
model = user
stopping_time = 0.1   # optional, used only as a timestep cap by the core module
cfl_drag = 0.05
```

The pgen must enroll:

```cpp
void MyParticleDrag(MeshBlockPack *pmbp, const Real bdt);

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_particle_drag_func = MyParticleDrag;
}
```

The callback receives the local `MeshBlockPack` and the operator-split time step.  It can
read `pmbp->ppart->prtcl_rdata`, `pmbp->ppart->prtcl_idata`, and the host hydro/MHD
primitive and conserved arrays, then apply any particle and fluid updates needed by the
problem.

For common drag callbacks, use the helper API in `src/particles/particles.hpp` instead of
hand-writing raw array indexing:

```cpp
DragParticleStencil stencil;
BuildDragParticleStencil(ppart->drag_interpolation, x1, x2, x3, x1min, x2min, x3min,
                         dx1, dx2, dx3, is, ie, js, je, ks, ke, multi_d, three_d,
                         stencil);

Real vg1, vg2, vg3;
InterpolateDragVelocity(w0, m, stencil, multi_d, three_d, vg1, vg2, vg3);

DepositDragBackreaction(u0, m, stencil, inv_vol, dmx, dmy, dmz, denergy, multi_d,
                        three_d, include_energy);
```

`inputs/particles/user_drag_hook.athinput` shows a minimal pgen-enrolled callback using
these helpers through the `streaming_instability` pgen's relaxation mode.

## AMR and MPI Behavior

Particles continue to use the existing MeshBlock boundary exchange for MPI migration.
This branch also adds `Particles::RemapAfterAMR()`:

1. Copy particle positions and IDs to host.
2. Wrap periodic coordinates back into the mesh domain.
3. Find the new leaf MeshBlock containing each particle after the AMR rebuild by logical
   location using the MeshBlock tree.
4. Update `PGID`.
5. Use the existing particle MPI exchange for particles now owned by another rank.
6. Recompute per-rank and global particle counts with `Mesh::UpdateParticleCounts()`.

The remap is called after AMR rebuilds mesh metadata and physics arrays.

The remap path avoids particle host/device copies on ranks with zero local particles, but
those ranks still participate in MPI send/receive discovery so sparse-rank AMR rebuilds
remain collective-safe.

## Restarts

Restart files now append particle metadata, per-rank particle counts, `prtcl_rdata`, and
`prtcl_idata` after the MeshBlock-strided fluid state.  On restart the particle arrays are
resized, particle counts are restored, and the pgen is called with `restart = true` so
user drag callbacks and history functions can be re-enrolled without overwriting particle
state.

## Validation Inputs

| Input | Purpose |
| --- | --- |
| `inputs/tests/particle_drag_relaxation.athinput` | Dust-gas drag damping and momentum check |
| `inputs/tests/particle_drag_mhd_relaxation.athinput` | Zero-field MHD host smoke test |
| `inputs/tests/particle_drag_restart.athinput` | Restart and user drag callback re-enrollment |
| `inputs/particles/hello_drag_particles.athinput` | Minimal standard drag example |
| `inputs/particles/user_drag_hook.athinput` | Minimal pgen user drag callback example |
| `inputs/particles/streaming_instability.athinput` | Uniform-grid streaming setup |
| `inputs/particles/streaming_instability_amr.athinput` | AMR/MPI streaming setup |
| `inputs/particles/streaming_instability_amr_dynamic.athinput` | Dynamic AMR remap check |

The pytest coverage is in `tst/test_suite/particles/`.  The helper scripts
`tst/scripts/particles/streaming_eigenmode.py` and
`tst/scripts/particles/profile_drag_particles.py` regenerate the streaming eigenmode and
run short drag-deposition timing sweeps.  The profiling helper accepts `--np` and
`--extra` so the same script can time dynamic-AMR MPI cases that exercise particle
sendlist packing as well as local deposition atomics.
