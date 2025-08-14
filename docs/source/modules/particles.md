# Module: Particles

## Overview
The Particles module implements Lagrangian particle tracking with various pusher algorithms, interpolation schemes, and particle-mesh coupling. The module supports two **mutually exclusive** particle types: star particles (from brent_test branch) and cosmic ray particles (from athenak-PK branch).

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
- **Compatible Pushers**: `rk4_gravity`, `drift`, `leap_frog`, `lagrangian_tracer`, `lagrangian_mc`
- **Incompatible**: Boris pushers (will cause fatal error)
- **From**: brent_test branch

### Cosmic Ray Particles  
- **Purpose**: Charged particles in magnetic fields
- **Memory**: 14 real data slots per particle (includes B-field tracking)
- **Compatible Pushers**: All, especially `boris_lin`, `boris_tsc`
- **From**: athenak-PK branch

## Configuration Parameters

### Common Parameters (`<particles>` block)

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `type` | string | Yes | Either `star` or `cosmic_ray` |
| `pusher` | string | Yes | Integration algorithm (see below) |
| `interpolation` | string | No | For boris: `lin` or `tsc` (default: `tsc`) |

### Star Particle Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `particle_file` | string | Yes | ASCII file with initial positions |
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
| `ppc` | Real | 1.0 | Particles per cell |
| `assign_tag` | string | index_order | Tagging method |
| `mass_log_spacing` | Real | 1.0 | Mass spectrum spacing |
| `min_mass` | Real | 1.0 | Minimum particle mass |

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
  
  // Cosmic ray only (indices 6-13)
  IPM=6,          // Mass/charge ratio
  IPBX=7, IPBY=8, IPBZ=9,  // B-field components
  IPDX=10, IPDY=11, IPDZ=12, IPDB=13  // Displacements
};
```

## Pusher Algorithms

### Compatibility Matrix

| Pusher | Star Particles | Cosmic Rays | Description |
|--------|---------------|-------------|-------------|
| `drift` | ✅ | ✅ | Simple drift motion |
| `rk4_gravity` | ✅ **Best** | ⚠️ | RK4 in gravitational potential |
| `leap_frog` | ✅ | ✅ | 2nd order leapfrog |
| `lagrangian_tracer` | ✅ | ✅ | Passive tracers |
| `lagrangian_mc` | ✅ | ✅ | Monte Carlo tracers |
| `boris` (→boris_lin/tsc) | ❌ **Fatal** | ✅ **Best** | Boris for E&M fields |

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
type = star
pusher = rk4_gravity
particle_file = initial_stars.txt
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
type = cosmic_ray
pusher = boris
interpolation = tsc
nspecies = 1
ppc = 8.0
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
- Cosmic rays: 14 slots (includes B-field tracking)
- ~36% memory savings for star particles vs unified approach

### Runtime Protection
1. **Initialization check**: Fatal error if incompatible pusher selected
2. **Pusher verification**: Boris pushers verify particle type
3. **Clear error messages**: Explains exactly what went wrong

Example error:
```
### FATAL ERROR in particles.cpp at line 198
Boris pushers are incompatible with star particles!
Boris pushers require magnetic field tracking arrays not allocated for stars.
Use 'rk4_gravity' or 'drift' pusher for star particles.
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

- Star particles: 72 bytes per particle (9×8)
- Cosmic rays: 112 bytes per particle (14×8)
- Memory savings with isolation: 40 bytes per star particle
- All pushers use Kokkos parallel_for for GPU acceleration

## References

- Original brent_test implementation for star particles
- athenak-PK implementation for magnetic field tracking
- Stone et al. (2020) for Athena++ particle framework