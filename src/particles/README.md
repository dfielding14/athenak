# Star Particles

This branch adds a standalone `particle_type = star` path to the existing particle
module.  The implementation keeps the original file-initialized use case, but it also
supports creating particles from gas and growing them through local gas accretion.

## Data Model

Star particles use the existing particle arrays:

| Index | Meaning |
| --- | --- |
| `PGID` | Global MeshBlock id containing the particle |
| `PTAG` | Globally unique integer particle tag |
| `IPX`, `IPY`, `IPZ` | Position |
| `IPVX`, `IPVY`, `IPVZ` | Velocity |
| `IPMASS` | Particle mass |
| `IPT_CREATE` | Creation time |

Only one particle type is active in a run.  The cosmic-ray path is left unchanged.

## Initialization

Set `particle_type = star` in the `<particles>` block.  By default, a run starts with no
star particles:

```ini
<particles>
particle_type = star
pusher = drift
star_init = none
```

File initialization is still available:

```ini
<particles>
particle_type = star
pusher = drift
star_init = file
particle_file = inputs/particles/star_particles.tbl
particle_file_columns = mass_then_time
```

The file format is ASCII with one particle per row:

```text
# x y z vx vy vz mass t_create
0.25 0.25 0.25 0.0 0.0 0.0 1.0 0.0
```

Set `particle_file_columns = time_then_mass` for legacy tables where the last two columns
are creation time followed by mass.

## Formation

Formation is enabled with:

```ini
star_formation = true
star_formation_density_threshold = 5.0
star_formation_particle_mass = 0.02
star_formation_density_floor = 1.0
star_formation_max_per_cycle = 1
```

Every `star_formation_interval` cycles, the module scans active hydro cells.  Cells with
conserved density above `star_formation_density_threshold` can create a star at the cell
center.  If `star_remove_gas_on_formation = true`, the new particle mass is removed from
the host gas cell while preserving the local gas velocity and specific energy.
`star_formation_max_per_cycle` limits the number of new particles created on each rank
per formation pass.

## Accretion

Accretion is enabled with:

```ini
star_accretion = true
star_accretion_rate = 1.0
star_accretion_radius_cells = 1
star_accretion_density_floor = 1.0
star_accretion_max_fraction = 0.05
```

Each cycle, every star removes a fraction
`min(star_accretion_max_fraction, star_accretion_rate*dt)` of gas above the density floor
from its host cell and the surrounding `star_accretion_radius_cells` stencil.  Removed gas
mass is added to `IPMASS`.  A replicated star snapshot lets the stencil cross MeshBlock
and MPI-rank boundaries while each rank updates only the active gas cells that it owns.
Stars are processed in stable tag order so overlapping stencils are decomposition
independent.

## Boundary Handling

The particle boundary task list keeps star particles attached to the owning MeshBlock and
rank.  Particles crossing periodic or shear-periodic root boundaries are wrapped.  Particles
crossing non-periodic root boundaries, such as outflow boundaries, are removed and the
mesh-level particle counts are updated collectively.  A particle update must cross at
most one MeshBlock in each direction.  The boundary task reports an explicit error when
`<particles>/dt` is too large for that requirement.

## Gravity

Set `pusher = gravity` and add a `<star_gravity>` input block to make stars respond to
particle gravity.  The current implementation supports:

- exact direct star-star gravity with Plummer-style softening;
- replicated Barnes-Hut tree star-star gravity with monopole moments;
- `kdk` and `rk4` particle integrators;
- constant external acceleration;
- fixed softened point-mass external acceleration;
- a pgen callback through `user_star_particle_accel_func` for custom external
  accelerations;
- gravity timestep estimates through `Particles::dtnew`;
- global history diagnostics for energy, momentum, angular momentum, center of mass,
  maximum acceleration, and minimum pair separation.

Both gravity backends gather the global particle snapshot with collectives and compute
accelerations for locally owned stars.  The direct backend remains the correctness
reference.  The tree backend builds a replicated host-side octree on each rank, uses
monopole moments with opening angle `tree_theta`, and is intended as the first large-N
option before a fully distributed tree or particle-mesh method.  Set
`exact_diagnostics = false` to skip the otherwise exact `O(N^2)` tree-history pass for
potential energy and minimum pair separation.

KDK integration is split around the existing migration path: the first kick and drift run
before particle boundary exchange, and the final kick runs after exchange.  Formation and
accretion happen before the first gravity acceleration, so newly formed particles are
included in the same update.  Accretion conserves particle momentum by default using the
velocity of the removed gas.

## Outputs And Tests

Particle history output is available through ordinary history files.  With `file_type =
hst`, particle runs now write `<basename>.part.hst` containing:

| Column | Meaning |
| --- | --- |
| `npart` | Total particle count |
| `part-mass` | Total particle mass |
| `formed` | Cumulative mass formed by threshold creation |
| `accreted` | Cumulative mass accreted from gas |

Regression inputs live in `inputs/particles/`:

- `star_particle_file.athinput`
- `star_particle_formation.athinput`
- `star_particle_outflow.athinput`
- `star_gravity_binary.athinput`
- `star_gravity_tree.athinput`
- `star_gravity_external.athinput`
- `star_gravity_formation.athinput`
- `star_gravity_accretion.athinput`
- `star_gravity_restart.athinput`
- `star_gravity_restart_formation.athinput`

Pytest coverage lives in:

- `tst/test_suite/particles/test_star_particles_cpu.py`
- `tst/test_suite/particles/test_star_particles_mpicpu.py`
- `tst/test_suite/particles/test_star_gravity_cpu.py`
- `tst/test_suite/particles/test_star_gravity_mpicpu.py`

## Current Limits

- Star particles require a 3D mesh.
- Formation and accretion require a hydro module.
- Accretion currently uses a cell-stencil prescription, not a Bondi or feedback model.
- The formation/accretion scan is host-side for clarity and testability; it is suitable
  for the current feature branch, but large production runs may want a device-side
  implementation.
- Gravity, neighboring-cell accretion, and sidecar reload currently replicate global
  star snapshots or payloads on each rank.  A distributed large-count path remains
  future work.
- Star-star gravity does not include gas self-gravity or star-gas gravitational
  backreaction.

## Sidecar Restarts

Mesh restart files still store the mesh state.  Star particles are restartable through a
matching `file_type = rst_prtcl` sidecar output in `rst_prtcl/`.  Restart with
`star_init = restart` and `particle_restart_file = rst_prtcl/<basename>.00010.rst_prtcl`.

The sidecar header records a magic/version, header byte size, byte-order marker, scalar
sizes, mesh time/cycle, total particle count, particle array dimensions, maximum tag,
formed/accreted totals, and basename.  The loader validates the header and complete byte
count against the mesh restart, restores particle real/int arrays exactly, recomputes
`PGID` from the current mesh geometry so rank counts can change, and refreshes
`next_star_tag` from restored tags rather than reassigning tags.  Sidecars are written
before mesh restarts regardless of output-block order so the continued run uses the next
unused sidecar number.
