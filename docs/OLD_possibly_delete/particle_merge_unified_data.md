# Unified Particle Data Structures

## Core Architecture

### Flexible Index System

```cpp
// athena.hpp
enum ParticlesIndex {
  // Common particle properties
  PGID=0, PTAG=1,     // Integer indices  
  IPX=0, IPVX=1,      // Real data position/velocity
  IPY=2, IPVY=3, 
  IPZ=4, IPVZ=5,
  
  // Type-specific (dynamically allocated)
  // Cosmic Ray: IPM=6, IPBX=7, IPBY=8, IPBZ=9, IPDX=10, IPDY=11, IPDZ=12, IPDB=13  
  // Star: IPT_CREATE=6, IPMASS=7, IPT_NEXT_SN=8
};

// Additional integer indices for particle type
enum ParticleIntegerIndex {
  PSP=2    // Species for cosmic rays
  NSN=2    // SN count for stars  
};
```

### Type-Aware Memory Allocation

```cpp
// particles.cpp constructor modification
void Particles::AllocateArrays() {
  switch(particle_type) {
    case ParticleType::cosmic_ray:
      // Full tracking arrays for CR physics
      nrdata = pmy_pack->pmesh->three_d ? 14 : 12;
      nidata = 3;  // PGID, PTAG, PSP
      break;
      
    case ParticleType::star:
      // Compact arrays for stellar evolution
      nrdata = 9;  
      nidata = 3;  // PGID, PTAG, NSN
      break;
  }
  
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
}
```

## Implementation Steps

### 1. Extend Particle Type Enum
```cpp
// particles.hpp
enum class ParticleType {cosmic_ray, star, mixed};
```

### 2. Add Species Support
```cpp
class Particles {
  // Cosmic ray specific
  int nspecies;         // Number of CR species
  bool track_displace;  // Enable displacement tracking
  
  // Star specific  
  bool enable_sn_feedback;
  Real sn_energy;
};
```

### 3. Constructor Branching
```cpp
Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) {
  std::string ptype = pin->GetString("particles","particle_type");
  
  if (ptype == "cosmic_ray") {
    nspecies = pin->GetOrAddInteger("particles","nspecies",1);
    track_displace = pin->GetOrAddBoolean("particles","track_displacement",false);
    InitializeCosmicRays(pin);
  } else if (ptype == "star") {
    enable_sn_feedback = pin->GetOrAddBoolean("particles","sn_feedback",true);
    InitializeStars(pin);
  }
}
```

### 4. Accessor Functions

```cpp
// Type-safe accessors prevent index confusion
inline Real& GetMass(int p) {
  return particle_type == ParticleType::cosmic_ray ? 
         prtcl_rdata(6,p) : prtcl_rdata(7,p);
}

inline Real& GetBField(int comp, int p) {
  assert(particle_type == ParticleType::cosmic_ray);
  return prtcl_rdata(7+comp,p);  // IPBX, IPBY, IPBZ
}

inline Real& GetDisplacement(int comp, int p) {
  assert(particle_type == ParticleType::cosmic_ray);
  return prtcl_rdata(10+comp,p); // IPDX, IPDY, IPDZ  
}
```

## Memory Layout

### Cosmic Ray Particles (14 elements)
```
[0-5]:  x, vx, y, vy, z, vz  (position/velocity)
[6]:    mass
[7-9]:  Bx, By, Bz            (magnetic field)  
[10-13]: dx, dy, dz, db       (displacement tracking)
```

### Star Particles (9 elements)
```
[0-5]: x, vx, y, vy, z, vz   (position/velocity)
[6]:   t_create               (formation time)
[7]:   mass                   (stellar mass)
[8]:   t_next_SN              (next supernova time)
```

## Optimization Notes

- Use template specialization for hot loops to avoid runtime branching
- Consider AoS → SoA transformation for better GPU memory coalescing
- Maintain separate kernels for each particle type to minimize divergence
- Pool memory allocation for dynamic particle creation/destruction