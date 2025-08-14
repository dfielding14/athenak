# Particle System Merge Summary

## Successfully Merged Features

### From brent_test branch:
- **Star particles** (`ParticleType::star`)
- **RK4 gravity pusher** for star particles in gravitational potentials
- Gravity potential parameters (r_scale, rho_scale, m_gal, a_gal, z_gal, r_200, rho_mean)
- Star particle file loading from ASCII format

### From athenak-PK branch:
- **Magnetic field tracking** for particles (IPBX, IPBY, IPBZ indices)
- **Displacement tracking** (IPDX, IPDY, IPDZ, IPDB indices)
- **Particle mass** index (IPM)
- **Species support** (PSP index)
- **Boris pushers** (boris_lin, boris_tsc) for charged particles in magnetic fields
- Enhanced particle boundary communication functions

### From main branch:
- **Cosmic ray particles** (`ParticleType::cosmic_ray`)
- Updated AMR features
- Kokkos 4.1.00 compatibility

## Critical Issues Resolved

### 1. Array Size Conflict (FIXED)
**Problem**: Different particle types used different array sizes:
- Star particles: 9 real data slots
- Cosmic rays with B-field: 14 real data slots

**Solution**: Unified all particle types to use 14 slots to ensure compatibility and prevent buffer overruns.

### 2. Pusher-Particle Compatibility (ADDRESSED)
**Problem**: Boris pushers access magnetic field indices that star particles didn't originally have.

**Solution**: 
- Allocated full array size for all particles
- Added runtime warnings for physically questionable combinations
- Preserved all pusher functionality

### 3. Index Definitions (MERGED)
Successfully merged all particle indices:
```cpp
enum ParticlesIndex {
  PGID=0, PTAG=1, PSP=2,  // Integer indices
  IPX=0, IPVX=1, IPY=2, IPVY=3, IPZ=4, IPVZ=5,  // Position & velocity
  IPM=6,   // Mass
  IPBX=7, IPBY=8, IPBZ=9,  // Magnetic field
  IPDX=10, IPDY=11, IPDZ=12, IPDB=13  // Displacements
};
```

## Available Pushers

1. **drift** - Simple drift motion
2. **rk4_gravity** - RK4 integration in gravitational potential (for stars)
3. **leap_frog** - Leap-frog integration
4. **lagrangian_tracer** - Passive tracer particles
5. **lagrangian_mc** - Lagrangian Monte Carlo
6. **boris_lin** - Boris pusher with linear interpolation (for cosmic rays)
7. **boris_tsc** - Boris pusher with TSC interpolation (for cosmic rays)

## Usage Recommendations

### For Star Particles:
```
<particles>
type = star
pusher = rk4_gravity
particle_file = star_particles.txt
</particles>

<potential>
r_scale = 1.0
rho_scale = 1.0
mass_gal = 1.0
scale_gal = 1.0
z_gal = 1.0
r_200 = 1.0
rho_mean = 1.0
</potential>
```

### For Cosmic Ray Particles:
```
<particles>
type = cosmic_ray
pusher = boris
interpolation = tsc  # or lin
</particles>
```

## Memory Overhead

The unified array size approach adds ~56% memory overhead for star particles (14 slots vs 9 needed), but this ensures:
- Complete safety from buffer overruns
- Full compatibility between all features
- Flexibility to switch pushers at runtime
- Future extensibility

## Testing Recommendations

1. Test star particles with rk4_gravity pusher
2. Test cosmic rays with boris pushers
3. Verify boundary communication with both particle types
4. Check AMR refinement with particles
5. Validate magnetic field interpolation for cosmic rays
6. Verify gravitational potential integration for stars

## Known Limitations

1. Star particles still require 3D simulations
2. Memory overhead for star particles (acceptable trade-off for safety)
3. Some pusher-particle combinations may not be physically meaningful (warnings added)

## Compilation Status

✅ **Successfully compiles and links** with all merged features