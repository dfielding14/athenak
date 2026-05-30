# Star Particles

The star-particle feature is a standalone particle type for collisionless stellar
objects coupled to hydro gas through formation and accretion.  It starts from the
Gotham file-initialized idea, but does not require all stars to be known at startup.

## Enabling The Module

```ini
<particles>
particle_type = star
pusher = drift
dt = 0.05
```

Star particles use the particle task list before the fluid time integrator:

1. Form new stars from eligible gas cells.
2. Accrete gas onto existing stars.
3. Push particle positions.
4. Exchange particles that crossed MeshBlock boundaries.

Formation and accretion require `<hydro>`.  File-only star-particle runs can be used
without enabling formation or accretion.

## File Initialization

```ini
star_init = file
particle_file = inputs/particles/star_particles.tbl
particle_file_columns = mass_then_time
```

Table columns are:

```text
x y z vx vy vz mass t_create
```

For legacy Gotham-style tables with creation time before mass, use:

```ini
particle_file_columns = time_then_mass
```

## Formation Parameters

```ini
star_formation = true
star_formation_density_threshold = 5.0
star_formation_particle_mass = 0.02
star_formation_density_floor = 1.0
star_formation_interval = 1
star_formation_max_per_cycle = 1
star_remove_gas_on_formation = true
```

A cell forms a star when its conserved density is above
`star_formation_density_threshold`.  The star is placed at the cell center and receives
the local gas velocity.  If gas removal is enabled, the code removes up to
`star_formation_particle_mass` from the host cell without taking the gas below
`star_formation_density_floor`.  `star_formation_max_per_cycle` is applied per rank.

## Accretion Parameters

```ini
star_accretion = true
star_accretion_rate = 1.0
star_accretion_radius_cells = 1
star_accretion_density_floor = 1.0
star_accretion_max_fraction = 0.05
```

For each star, the code visits the host cell and neighboring cells within
`star_accretion_radius_cells`.  In each cell it removes a fraction

```text
min(star_accretion_max_fraction, star_accretion_rate*dt)
```

of the gas above `star_accretion_density_floor`, then adds that mass to the particle.
Momentum, total energy, and scalar conserved variables are removed in the same fraction
so the gas velocity and composition remain unchanged by the sink operation.

The stencil crosses MeshBlock and MPI rank boundaries.  Every rank evaluates its owned
active cells against the same replicated star snapshot, and the removed mass and
momentum are reduced back to the owning star.  Stable particle-tag order makes
overlapping stencils independent of rank-local storage order.

## Boundary Behavior

Star particles use the existing particle boundary task list.  Particles that cross an
internal MeshBlock boundary are moved to the new owning MeshBlock and, under MPI, to the
owning rank.  Particles that leave the root mesh through a non-periodic boundary are
removed from the particle arrays and from the global particle counts.  Periodic and
shear-periodic boundaries retain the wrapped-particle behavior of the base particle
module.  A particle update must cross at most one MeshBlock in each direction.  AthenaK
stops with an explicit error if `<particles>/dt` is too large for that requirement.

## Gravity

The gravity-aware particle pusher is documented separately in
`star_particle_gravity.md`.  In brief, `pusher = gravity` and a `<star_gravity>` block
allow star particles to feel exact direct star-star gravity and configured external
accelerations.  This is separate from hydro/MHD source terms: stars do not automatically
inherit a pgen's gas source term unless the pgen explicitly provides a star-particle
acceleration callback.

## Sidecar Restarts

Star particles are restartable with a particle sidecar output:

```ini
<output2>
file_type = rst
dt = 0.1

<output3>
file_type = rst_prtcl
dt = 0.1
```

The mesh state remains in `rst/<basename>.<file_number>.rst`; particle real/int arrays
and star bookkeeping totals are stored in
`rst_prtcl/<basename>.<file_number>.rst_prtcl`.  Restart with:

```bash
./athena -r rst/my_run.00010.rst \
  particles/star_init=restart \
  particles/particle_restart_file=rst_prtcl/my_run.00010.rst_prtcl
```

The sidecar loader validates time and cycle against the mesh restart, restores tags
without reassignment, recomputes particle MeshBlock ownership from the current geometry,
and refreshes the next dynamic star tag before any new formation.  Particle sidecars are
written before mesh restart files regardless of output-block order, so resumed runs
continue with the next unused sidecar number.

## History Output

Any run with particles and history output writes `<basename>.part.hst`.

| Quantity | Description |
| --- | --- |
| `npart` | Particle count |
| `part-mass` | Total particle mass |
| `formed` | Cumulative mass created through star formation |
| `accreted` | Cumulative gas mass accreted onto particles |

## Regression Coverage

The branch includes CPU regression tests for:

- File initialization preserves particle count and total mass.
- A dense hydro cell forms one star and grows it by accretion.
- A star that crosses an outflow boundary is removed cleanly.

It also includes MPI CPU regression tests for:

- File initialization with particles distributed over two ranks.
- Formation/accretion when the new particle is created on only one rank.
- Outflow-boundary removal when the particle is owned by one MPI rank and the global
  particle history is written by rank 0.

Run only these tests with:

```bash
cd tst
python run_test_suite.py --cpu --test test_suite/particles/test_star_particles_cpu.py
python run_test_suite.py --mpicpu --test test_suite/particles/test_star_particles_mpicpu.py
```

## Notes

This implementation is intentionally conservative.  The host-side formation and
accretion scan is straightforward and portable, and the task ordering keeps all gas
updates before the hydro time integrator.  More specialized star-formation criteria,
feedback, stochastic creation, or device-side large-scale formation scans can be added
on top of this branch without changing the file initialization path.  Neighboring-cell
accretion currently replicates a compact star snapshot on each rank; a distributed sink
search is a future optimization for very large star counts.
