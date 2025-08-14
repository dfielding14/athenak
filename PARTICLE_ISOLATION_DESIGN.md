# Particle System Isolation Design

## Executive Summary

The brent_test (star particles) and athenak-PK (cosmic ray particles) systems are **mutually exclusive** - they cannot be used simultaneously in the same simulation. This design decision eliminates memory overhead and ensures type safety.

## Key Design Principles

### 1. Type-Specific Memory Allocation
- **Star particles**: Allocate exactly 9 real data slots (no waste)
- **Cosmic ray particles**: Allocate 14 real data slots for magnetic field tracking
- **Memory savings**: 40 bytes per star particle (vs unified approach)

### 2. Compile-Time Configuration, Runtime Enforcement
- Particle type selected via input file at runtime
- Array sizes determined based on particle type
- Incompatible pusher-particle combinations cause **fatal errors** (not warnings)

### 3. Defense in Depth
Three layers of protection against buffer overruns:
1. **Initialization check**: Fatal error if incompatible pusher selected
2. **Runtime check**: Boris pushers verify particle type before execution
3. **Clear documentation**: Explicit warnings about incompatibility

## Array Allocation Details

### Star Particles (9 slots)
```cpp
IPX(0), IPVX(1), IPY(2), IPVY(3), IPZ(4), IPVZ(5)  // Position & velocity
[6-8]: Star-specific data (mass, age, metallicity, etc.)
```

### Cosmic Ray Particles (14 slots)
```cpp
IPX(0), IPVX(1), IPY(2), IPVY(3), IPZ(4), IPVZ(5)  // Position & velocity
IPM(6)                                               // Mass/charge ratio
IPBX(7), IPBY(8), IPBZ(9)                          // Magnetic field
IPDX(10), IPDY(11), IPDZ(12), IPDB(13)            // Displacement tracking
```

## Pusher Compatibility Matrix

| Pusher | Star Particles | Cosmic Ray Particles |
|--------|---------------|---------------------|
| drift | ✅ Yes | ✅ Yes |
| rk4_gravity | ✅ **Recommended** | ⚠️ Works but unusual |
| leap_frog | ✅ Yes | ✅ Yes |
| lagrangian_tracer | ✅ Yes | ✅ Yes |
| lagrangian_mc | ✅ Yes | ✅ Yes |
| boris_lin | ❌ **Fatal Error** | ✅ **Recommended** |
| boris_tsc | ❌ **Fatal Error** | ✅ **Recommended** |

## Error Handling

### Fatal Errors (Program terminates)
1. Attempting to use boris pushers with star particles
2. Unknown particle type specified
3. Missing pusher specification

### Warnings (Program continues)
1. Using rk4_gravity with cosmic ray particles (unusual but allowed)

## Usage Examples

### Star Particle Simulation
```ini
<particles>
type = star
pusher = rk4_gravity  # or drift, leap_frog, etc.
particle_file = stars.txt
</particles>

<potential>
r_scale = 1.0
rho_scale = 1.0  
mass_gal = 1.0
# ... other gravity parameters
</potential>
```

### Cosmic Ray Simulation
```ini
<particles>
type = cosmic_ray
pusher = boris       # becomes boris_lin or boris_tsc
interpolation = tsc  # or lin
</particles>
```

## Implementation Safety Features

### 1. Initialization Phase (particles.cpp)
```cpp
switch (particle_type) {
  case ParticleType::star:
    nrdata = 9;  // Exact allocation
    if (pusher == boris_lin || pusher == boris_tsc) {
      FATAL_ERROR("Boris pushers incompatible with stars!");
    }
    break;
  case ParticleType::cosmic_ray:
    nrdata = 14;  // Full magnetic field tracking
    break;
}
```

### 2. Runtime Phase (particles_pushers.cpp)
```cpp
case ParticlesPusher::boris_lin:
  if (particle_type != ParticleType::cosmic_ray) {
    FATAL_ERROR("Boris requires cosmic rays!");
  }
  // ... pusher code ...
```

## Benefits of This Approach

1. **Zero memory waste** for star particles
2. **Type safety** - impossible to corrupt memory with wrong pusher
3. **Clear semantics** - each simulation has one particle type
4. **Future extensibility** - easy to add new particle types
5. **Performance** - no runtime overhead from unused fields

## Limitations

1. Cannot mix particle types in single simulation
2. Must restart to change particle type
3. Some generic pushers may not be optimal for all types

## Migration Guide

### From Unified Approach
If you were using the unified 14-slot approach:
- Star simulations will use less memory (automatic)
- No code changes needed
- Boris+star combinations now properly error out

### From Original Branches
- brent_test users: Your star particles work unchanged
- athenak-PK users: Your cosmic rays work unchanged
- Mixed usage: Not supported (by design)

## Testing Checklist

- [ ] Star particles with rk4_gravity pusher
- [ ] Star particles with drift pusher
- [ ] Cosmic rays with boris_lin pusher
- [ ] Cosmic rays with boris_tsc pusher
- [ ] Error handling: boris + star (should fail)
- [ ] Warning: rk4_gravity + cosmic ray
- [ ] Memory usage comparison
- [ ] Performance benchmarks