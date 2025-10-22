# Cosmic Ray Pusher Porting Guide

## Boris Pusher Implementation

### Core Algorithm
The Boris pusher solves charged particle motion in electromagnetic fields using a centered time-stepping scheme that preserves energy exactly.

```cpp
// particles_pushers.cpp
TaskStatus Particles::PushCosmicRays(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  Real dt = pmy_pack->pmesh->dt;
  
  // Select interpolation method
  bool use_tsc = (pusher == ParticlesPusher::boris_tsc);
  
  par_for("cr_push", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    // Get particle properties
    Real x = pr(IPX,p), y = pr(IPY,p), z = pr(IPZ,p);
    Real vx = pr(IPVX,p), vy = pr(IPVY,p), vz = pr(IPVZ,p);
    Real q_over_m = pr(IPM,p);  // charge-to-mass ratio
    
    // Interpolate fields to particle position
    Real Bx, By, Bz, Ex, Ey, Ez;
    if (use_tsc) {
      InterpolateTSC(x, y, z, Bx, By, Bz, Ex, Ey, Ez);
    } else {
      InterpolateLinear(x, y, z, Bx, By, Bz, Ex, Ey, Ez);
    }
    
    // Boris algorithm
    Real dt_half = 0.5 * dt;
    Real qdt_2m = q_over_m * dt_half;
    
    // Half electric acceleration
    vx += qdt_2m * Ex;
    vy += qdt_2m * Ey;
    vz += qdt_2m * Ez;
    
    // Magnetic rotation
    Real tx = qdt_2m * Bx;
    Real ty = qdt_2m * By;
    Real tz = qdt_2m * Bz;
    Real t2 = tx*tx + ty*ty + tz*tz;
    Real s = 2.0 / (1.0 + t2);
    
    Real vx1 = vx + vy*tz - vz*ty;
    Real vy1 = vy + vz*tx - vx*tz;
    Real vz1 = vz + vx*ty - vy*tx;
    
    vx += s*(vy1*tz - vz1*ty);
    vy += s*(vz1*tx - vx1*tz);
    vz += s*(vx1*ty - vy1*tx);
    
    // Half electric acceleration
    vx += qdt_2m * Ex;
    vy += qdt_2m * Ey;
    vz += qdt_2m * Ez;
    
    // Update position
    pr(IPX,p) += dt * vx;
    pr(IPY,p) += dt * vy;
    pr(IPZ,p) += dt * vz;
    
    // Store updated velocity and fields
    pr(IPVX,p) = vx;
    pr(IPVY,p) = vy;
    pr(IPVZ,p) = vz;
    pr(IPBX,p) = Bx;
    pr(IPBY,p) = By;
    pr(IPBZ,p) = Bz;
    
    // Update displacement tracking
    if (track_displace) {
      pr(IPDX,p) += dt * vx;
      pr(IPDY,p) += dt * vy;
      pr(IPDZ,p) += dt * vz;
      pr(IPDB,p) += dt * (vx*Bx + vy*By + vz*Bz) / sqrt(Bx*Bx + By*By + Bz*Bz);
    }
  });
}
```

## Field Interpolation Methods

### Linear Interpolation
```cpp
KOKKOS_INLINE_FUNCTION
void InterpolateLinear(Real x, Real y, Real z,
                      Real &Bx, Real &By, Real &Bz,
                      Real &Ex, Real &Ey, Real &Ez) {
  // Find cell indices
  int i = static_cast<int>((x - x1min) / dx1);
  int j = static_cast<int>((y - x2min) / dx2);
  int k = static_cast<int>((z - x3min) / dx3);
  
  // Compute weights
  Real wx = (x - x1v(i)) / dx1;
  Real wy = (y - x2v(j)) / dx2;
  Real wz = (z - x3v(k)) / dx3;
  
  // Trilinear interpolation
  Bx = (1-wx)*(1-wy)*(1-wz)*b.x1f(m,k,j,i) +
       wx*(1-wy)*(1-wz)*b.x1f(m,k,j,i+1) + ...;
  // Similar for By, Bz, Ex, Ey, Ez
}
```

### TSC (Triangular-Shaped Cloud) Interpolation
```cpp
KOKKOS_INLINE_FUNCTION
void InterpolateTSC(Real x, Real y, Real z,
                   Real &Bx, Real &By, Real &Bz,
                   Real &Ex, Real &Ey, Real &Ez) {
  // TSC uses 27-point stencil for smoother interpolation
  // Better conservation properties but more expensive
  
  // Compute TSC weights
  Real wx[3], wy[3], wz[3];
  ComputeTSCWeights(x, wx);
  ComputeTSCWeights(y, wy);
  ComputeTSCWeights(z, wz);
  
  // Apply 3x3x3 stencil
  Bx = 0.0;
  for (int ii=-1; ii<=1; ++ii) {
    for (int jj=-1; jj<=1; ++jj) {
      for (int kk=-1; kk<=1; ++kk) {
        Bx += wx[ii+1]*wy[jj+1]*wz[kk+1]*b.x1f(m,k+kk,j+jj,i+ii);
      }
    }
  }
}
```

## Pusher Selection Logic

```cpp
// particles.hpp
enum class ParticlesPusher {
  drift,         // Simple advection
  boris_lin,     // Boris with linear interpolation
  boris_tsc,     // Boris with TSC interpolation
  rk4_gravity,   // RK4 for stars
  leap_frog      // Symplectic integrator
};

// particles_pushers.cpp
TaskStatus Particles::Push(Driver *pdriver, int stage) {
  switch(pusher) {
    case ParticlesPusher::drift:
      return PushDrift(pdriver, stage);
    case ParticlesPusher::boris_lin:
    case ParticlesPusher::boris_tsc:
      return PushCosmicRays(pdriver, stage);
    case ParticlesPusher::rk4_gravity:
      return PushStars(pdriver, stage);
    default:
      return TaskStatus::fail;
  }
}
```

## Multiple Species Support

```cpp
void InitializeCosmicRays(ParameterInput *pin) {
  nspecies = pin->GetOrAddInteger("particles","nspecies",1);
  
  // Allocate per-species properties
  Kokkos::realloc(species_mass, nspecies);
  Kokkos::realloc(species_charge, nspecies);
  
  // Read species parameters
  for (int s=0; s<nspecies; ++s) {
    std::string block = "species" + std::to_string(s);
    species_mass(s) = pin->GetReal(block,"mass");
    species_charge(s) = pin->GetReal(block,"charge");
  }
  
  // Initialize particles with species tags
  par_for("init_species", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    int species = p % nspecies;  // Simple round-robin
    pi(PSP,p) = species;
    pr(IPM,p) = species_charge(species) / species_mass(species);
  });
}
```

## Input Parameters

```ini
<particles>
particle_type     = cosmic_ray
pusher           = boris_tsc
ppc              = 10          # particles per cell
nspecies         = 2           # number of species
track_displacement = true       # enable tracking

<species0>
mass             = 1.0
charge           = 1.0

<species1>  
mass             = 1836.0      # proton mass
charge           = 1.0
```

## Performance Considerations

- Boris pusher is ~2x more expensive than drift
- TSC interpolation is ~3x slower than linear but more accurate
- Displacement tracking adds 30% overhead
- Multiple species increases memory by `nspecies * sizeof(properties)`
- GPU optimization: coalesce field access, minimize divergence