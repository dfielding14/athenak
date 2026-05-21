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
mass is added to `IPMASS`.

## Boundary Handling

The particle boundary task list keeps star particles attached to the owning MeshBlock and
rank.  Particles crossing periodic or shear-periodic root boundaries are wrapped.  Particles
crossing non-periodic root boundaries, such as outflow boundaries, are removed and the
mesh-level particle counts are updated collectively.

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

Pytest coverage lives in:

- `tst/test_suite/particles/test_star_particles_cpu.py`
- `tst/test_suite/particles/test_star_particles_mpicpu.py`

## Current Limits

- Star particles require a 3D mesh.
- Formation and accretion require a hydro module.
- Accretion currently uses a cell-stencil prescription, not a Bondi or feedback model.
- The formation/accretion scan is host-side for clarity and testability; it is suitable
  for the current feature branch, but large production runs may want a device-side
  implementation.
