# Module: Particles

## Overview
The Particles module implements Lagrangian particle tracking with pusher algorithms,
interpolation schemes, and particle-mesh coupling. The module supports two
**mutually exclusive** particle types: star particles and cosmic ray particles.

## Source Location
`src/particles/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `particles.hpp/cpp` | Core particle class | Initialization, type selection, memory management |
| `particles_pushers.cpp` | Integration algorithms | `rk4_gravity()`, `boris_lin()`, `boris_tsc()` |
| `particles_tasks.cpp` | Task registration | Particle evolution tasks |
| `bvals/bvals_part.cpp` | Boundary communication | MPI particle exchange |

## Particle Types (Mutually Exclusive)

### Star Particles
- **Purpose**: Represent stellar objects in gravitational potentials
- **Memory**: 9 real data slots per particle
- **Compatible Pushers**: `rk4_gravity`, `drift`
- **Incompatible**: Boris pushers (will cause fatal error)

### Cosmic Ray Particles
- **Purpose**: Charged particles in magnetic fields
- **Memory**: 23 real data slots per particle, including sampled fields,
  displacement/feedback diagnostics, and per-particle macro weight
- **Compatible Pushers**: `drift`, `boris_lin`, `boris_tsc`

## Configuration Parameters

### Common Parameters (`<particles>` block)

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `particle_type` | string | Yes | Either `star` or `cosmic_ray` |
| `pusher` | string | Yes | Integration algorithm (see below) |
| `pic_interp_scheme` | string | No | PIC interpolation scheme; currently only `tsc` is supported |

### Star Particle Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `star_particle_file` | string | Yes | ASCII file with initial positions |
| `grav_dx` | Real | No | Finite difference spacing for gravity (default: 1e-6) |

**Required `<potential>` block parameters:**
```ini
r_scale     # Radial scale
rho_scale   # Density scale
mass_gal    # Galaxy mass
scale_gal   # Galaxy scale length
z_gal       # Vertical scale
r_200       # Virial radius
rho_mean    # Mean density
```

### Cosmic Ray Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nspecies` | int | 1 | Number of species |
| `ppc` | Real | 1.0 | Particles per cell; must be non-negative; fractional counts are distributed by global MeshBlock ID |
| `cr_distribution` | string | center | `center` places particles on deterministic cell centers; `random` uses deterministic per-block pseudorandom positions |
| `assign_tag` | string | index_order | Tagging method |
| `deposit_moments` | bool | false | Deposit particle charge/current moments |
| `deposit_qscale` | Real | 1.0 | Root-grid macro-charge and macro-mass scale |
| `pic_background_mode` | string | coupled | `coupled`, `passive_mhd`, or `no_mhd` |
| `pic_feedback_mode` | string | mode-dependent | `coupled` or `test_particle` |
| `pic_enable_2d3v` | bool | false | Required for Boris pushers on 2D meshes; keeps `vz` and `Bz` active when `nx3=1` |
| `pic_cr_light_speed` | Real | 1.0 | Reserved; non-default values are rejected |
| `pic_max_cell_cross` | int | 2 | Particle cell-crossing timestep limit before global CFL scaling; must not exceed the smallest active MeshBlock dimension |
| `pic_theta_max` | Real | 0.3 | Boris gyro-angle timestep limit before global CFL scaling |
| `pic_random_seed` | int | 0 | Seed for deterministic `cr_distribution=random` placement |
| `pic_deltaf_mode` | string | off | Staged quiet-start control; accepted values are `off` and `on` |
| `pic_deltaf_f0` | string | empty | Required when `pic_deltaf_mode=on`; accepted staged labels are `kappa_iso` and `uniform_quiet` |

Each `<speciesN>` block defines at least `mass` and `charge`. Species masses
must be positive.

When `couple_j_deposition_mode=direct_staggered`, the deposited edge current is
constructed from each particle trajectory over the step. The implementation
requires the old-to-new shape support to shift by no more than one cell in each
active dimension and aborts if a step is too large. Reduce `time/cfl_number` or
`pic_max_cell_cross` if that guard is triggered. Direct edge-current physical
BCs support periodic, reflecting, and outflow boundaries; use
`couple_j_deposition_mode=cc_convert` for inflow-boundary coupled runs.

## Particle Data Arrays

### Index Definitions
```cpp
enum ParticlesIndex {
  // Integer indices
  PGID=0,  // Global MeshBlock ID
  PTAG=1,  // Particle tag
  PSP=2,   // Species index

  // Real data indices
  IPX=0, IPVX=1,  // X position, velocity
  IPY=2, IPVY=3,  // Y position, velocity
  IPZ=4, IPVZ=5,  // Z position, velocity

  // Cosmic ray only (indices 6-22)
  IPM=6,          // Charge/mass ratio
  IPBX=7, IPBY=8, IPBZ=9,  // B-field components
  IPDX=10, IPDY=11, IPDZ=12, IPDB=13, // Displacements
  IPEX=14, IPEY=15, IPEZ=16,          // Electric-field samples
  IPDPX=17, IPDPY=18, IPDPZ=19,       // Feedback momentum rate
  IPDE=20, IPEBDOT=21,                // Feedback energy diagnostics
  IPWT=22                             // Relative macro-particle weight
};
```

Cosmic-ray particles created from `ppc` carry `IPWT = cell_volume/root_cell_volume`
so refined-cell particles keep the same physical density represented by a
root-level `deposit_qscale`. Manually injected particles may set `IPWT=1.0`
when `deposit_qscale` is already the intended macro-particle mass/charge.

`pic_deltaf_mode=on` currently selects deterministic low-discrepancy placement
for `cr_distribution=random` tests. It does not yet implement full delta-f
particle-weight evolution or use `pic_deltaf_f0` in moment deposition.

## Pusher Algorithms

### Compatibility Matrix

| Pusher | Star Particles | Cosmic Rays | Description |
|--------|---------------|-------------|-------------|
| `drift` | yes | yes | Simple drift motion |
| `rk4_gravity` | yes | no | RK4 in gravitational potential |
| `boris_lin` | no | yes | Boris pusher with linear field interpolation |
| `boris_tsc` | no | yes | Boris pusher with TSC field interpolation |

### RK4 Gravity Pusher (Stars)
4th-order Runge-Kutta integration in gravitational potential:
```cpp
// Evaluates gravitational acceleration via finite differences
// Uses potential parameters from <potential> block
```

### Boris Pushers (Cosmic Rays)
Relativistic particle motion in electromagnetic fields:
- `boris_lin`: Linear interpolation of fields
- `boris_tsc`: Triangular-shaped cloud (TSC) interpolation

## Usage Examples

### Star Particle Simulation
```ini
<job>
problem_id = stellar_dynamics

<particles>
particle_type = star
pusher = rk4_gravity
star_particle_file = initial_stars.txt
grav_dx = 1.0e-6

<potential>
r_scale = 8.0      # kpc
rho_scale = 1.0e7  # Msun/kpc^3
mass_gal = 1.0e12  # Msun
scale_gal = 3.5    # kpc
z_gal = 0.5        # kpc
r_200 = 200.0      # kpc
rho_mean = 1.0e6   # Msun/kpc^3
```

### Cosmic Ray Simulation
```ini
<job>
problem_id = cosmic_ray_propagation

<particles>
particle_type = cosmic_ray
pusher = boris_tsc
pic_interp_scheme = tsc
nspecies = 1
ppc = 8.0
pic_cr_light_speed = 1.0
pic_max_cell_cross = 2
pic_theta_max = 0.3

<species0>
mass = 1.0
charge = 1.0
```

## Star Particle File Format
ASCII file with 8 columns (3D only):
```
# x y z vx vy vz mass property
0.1 0.2 0.3 1.0 0.5 0.0 1.0e6 0
0.2 0.3 0.4 0.8 0.6 0.1 2.0e6 1
```

## Safety Features

### Type-Specific Memory Allocation
- Star particles: Exactly 9 slots (no waste)
- Cosmic rays: 23 slots, including field samples, feedback diagnostics, and
  relative macro-particle weight

### Runtime Protection
1. **Initialization check**: Fatal error if incompatible pusher selected
2. **Pusher verification**: Boris pushers require cosmic-ray particles; RK4
   gravity requires star particles
3. **Clear error messages**: Explains exactly what went wrong

Example error:
```
### FATAL ERROR in particles.cpp at line 198
Boris pushers are incompatible with star particles; use
<particles>/particle_type=cosmic_ray for boris_lin/boris_tsc, or use
drift/rk4_gravity for star particles.
```

## Implementation Notes

### Mutual Exclusivity
- Only ONE particle type per simulation
- Types cannot be mixed in same run
- Must restart to change particle type
- This is by design for safety and efficiency

### Memory Layout
- Real data: `prtcl_rdata[nrdata][nparticles]`
- Integer data: `prtcl_idata[nidata][nparticles]`
- Kokkos Views for GPU portability
- PIC restart files store particle records, `nrdata`, and `nidata` for every
  active particle mode. Runs with deposited moments also persist the moment
  arrays, and edge-current coupling persists its staggered current arrays.
  Changing the particle record layout intentionally rejects older incompatible
  PIC restart files.

### Boundary Communication
- Full MPI support for parallel runs
- Particles automatically exchanged between MeshBlocks
- Handles AMR level transitions

## Limitations

1. Star particles require 3D simulations
2. Particle types cannot be mixed
3. No dynamic type switching
4. Boris pushers only work with cosmic rays

## Performance Considerations

- Star particles: 72 bytes per particle (9x8)
- Cosmic rays: 184 bytes per particle for real data (23x8), plus integer data
- All pushers use Kokkos parallel_for for GPU acceleration

## References

- Stone et al. (2020) for Athena++ particle framework
