# CGM Cooling Flow AMR with Metals: Complete Technical Documentation

```{warning}
This is internal technical documentation for the CGM cooling flow with metals implementation.
This page is intentionally not linked from the main documentation.
```

## Table of Contents

1. [Overview](#overview)
2. [Physical Model](#physical-model)
3. [Code Architecture](#code-architecture)
4. [Implementation Details](#implementation-details)
5. [Input Parameters](#input-parameters)
6. [Algorithms and Numerical Methods](#algorithms-and-numerical-methods)
7. [Performance Considerations](#performance-considerations)
8. [Usage Examples](#usage-examples)
9. [Initial Condition Generation](ic_generation.md)
10. [Technical Reference](technical_reference.md)

## Overview

The `cgm_cooling_flow_amr_metals.cpp` problem generator implements a sophisticated circumgalactic medium (CGM) simulation with:

- **Multi-phase gas dynamics**: Hot CGM halo with embedded cool disk
- **Metal tracking**: Passive scalar tracking of metallicity
- **Supernova feedback**: Particle-based SN injection with FIRE-3 delay times
- **Adaptive Mesh Refinement**: Density gradient-based refinement
- **Multi-component gravity**: NFW + Miyamoto-Nagai + outer halo
- **Cooling flows**: Profile-based inflow solutions
- **Rotation**: Both CGM rotation and disk Keplerian rotation

## Physical Model

### 1. Gravitational Potential

The total gravitational potential consists of three components:

$$\Phi_{\text{total}} = \Phi_{\text{NFW}} + \Phi_{\text{MN}} + \Phi_{\text{outer}}$$

#### NFW Dark Matter Halo
```cpp
Real x = r / r_s;
Real phi_NFW = -4 * M_PI * G * rho_s * r_s^2 * log(1 + x) / x;
```
- Scale radius: `r_s` (typically 10-30 kpc)
- Scale density: `rho_s` 
- Provides dominant gravity at intermediate radii

#### Miyamoto-Nagai Disk
```cpp
Real phi_MN = -G * M_gal / sqrt(R^2 + (sqrt(z^2 + z_gal^2) + a_gal)^2);
```
- Galaxy mass: `M_gal` (typically $10^{10}$-$10^{11}$ $M_\odot$)
- Scale length: `a_gal` (few kpc)
- Scale height: `z_gal` (< 1 kpc)
- Provides disk gravity

#### Outer Halo
```cpp
Real term1 = (4.0/3.0) * pow(5*R200, 1.5) * sqrt(r);
Real term2 = (1.0/6.0) * r^2;
Real phi_Outer = 4 * M_PI * G * rho_mean * (term1 + term2);
```
- Virial radius: `R200`
- Mean density: `rho_mean`
- Provides confinement at large radii

### 2. CGM Cooling Flow

The cooling flow solution is read from pre-computed 1D profiles containing:
- Density profile: $\rho(r)$
- Temperature profile: $T(r)$
- Radial velocity profile: $v_r(r)$

The profiles are loaded via `ProfileReader` class and interpolated onto the 3D grid:

```cpp
SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
```

Key features:
- Radial inflow velocity: $\vec{v} = -v_r(r) \hat{r}$
- Smooth velocity tapering near origin: $v_r \to 0$ as $r \to 0$
- Pressure from ideal gas: $P = \rho T$ (in code units)

### 3. Disk Equilibrium

The disk is initialized in hydrostatic equilibrium with rotation:

```cpp
SetEquilibriumState(u0, m, k, j, i, x1v, x2v, x3v, ...);
```

The disk follows:
- Vertical hydrostatic equilibrium: $\partial P/\partial z = -\rho \partial\Phi/\partial z$
- Radial force balance: $v_\phi^2/R = (1/\rho)\partial P/\partial R - \partial\Phi/\partial R$
- Isothermal assumption: $c_s^2 = T$ (constant)
- Density stratification: $\rho(R,z) = \rho_0(R) \exp[-(\Phi(R,z) - \Phi(R,0))/c_s^2]$

### 4. Rotation

Two rotation components are implemented:

#### CGM Rotation
```cpp
if (r < r_circ) {
    v_phi = v_circ * sin_theta;
} else {
    v_phi = v_circ * sin_theta * r_circ / r;
}
```
- Solid body rotation within $r_{circ}$
- Keplerian falloff beyond $r_{circ}$

#### Disk Rotation
Computed from radial force balance:
```cpp
v_phi = sqrt(R * max(dP_dR/rho - dPhi_dR, 0.0));
```

### 5. Metal Tracking

Metallicity is tracked as a passive scalar:
```cpp
u0(m, nhydro, k, j, i) = Z * Zsol * u0(m, IDN, k, j, i);
```
- Solar metallicity: $Z_{\odot} = 0.02$
- Initial metallicity: $Z$ (typically $1/3$ solar)
- Advected with gas flow

## Code Architecture

### Global Variables

```cpp
namespace {
  // Gravitational parameters
  Real r_scale, rho_scale, m_gal, a_gal, z_gal, r_200, rho_mean;
  
  // Rotation parameters
  Real r_circ, v_circ;
  
  // Metallicity
  Real Z;
  
  // SN feedback
  Real r_inj, e_sn, m_ej;
  
  // Profile readers
  ProfileReaderHost profile_reader_host;      // CGM profiles
  ProfileReader profile_reader;               // Device-side
  ProfileReaderHost disk_profile_reader_host; // Disk profiles
  ProfileReader disk_profile_reader;          // Device-side
  
  // AMR threshold
  Real ddens_threshold;
}
```

### Key Functions

#### Problem Generator
```cpp
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart)
```
- Reads parameters
- Loads profile files
- Initializes grid with cooling flow + disk + rotation
- Sets up user functions

#### Source Terms
```cpp
void UserSource(Mesh* pm, const Real bdt)
  ├── GravitySource(pm, bdt)  // Gravitational acceleration
  └── SNSource(pm, bdt)        // Supernova feedback
```

#### Boundary Conditions
```cpp
void UserBoundary(Mesh* pm)
```
Applies cooling flow profile at all boundaries

#### Refinement
```cpp
void RefinementCondition(MeshBlockPack* pmbp)
```
Density gradient-based AMR criterion

## Implementation Details

### 1. Profile Reading System

The profile reader uses a two-stage system for GPU compatibility:

```cpp
// Host-side loading
ProfileReaderHost profile_reader_host;
profile_reader_host.ReadProfiles(profile_file);

// Create device-accessible copy
ProfileReader profile_reader = profile_reader_host.CreateDeviceReader();
```

Profile format (ASCII):
```
# r[kpc]  T[K]  rho[g/cm^3]  v[km/s]
0.1       1e6   1e-25        -10.0
0.2       1e6   8e-26        -12.0
...
```

### 2. Supernova Feedback Implementation

#### Particle-based SN tracking
Each star particle tracks:
- Creation time: `pr(nrdata-3, p)`
- Cluster mass: `pr(nrdata-2, p)` 
- Next SN time: `pr(nrdata-1, p)`
- SN count: `pi(2, p)`

#### SN Delay Time Distribution (FIRE-3)

```cpp
KOKKOS_INLINE_FUNCTION
Real GetNthSNTime(Real M, Real t_start, Real unit_time, int n)
```

Implements piecewise power-law DTD:
- Early CC SNe (3.7-7 Myr): $\propto t^{s_1}$
- Late CC SNe (7-44 Myr): $\propto t^{s_2}$  
- Type Ia SNe (>44 Myr): $\propto t^{-1.1}$

#### Energy Injection
```cpp
if (r <= r_inj) {
    u0(m,IEN,k,j,i) += e_sn;  // Thermal energy
}
```
- Injection radius: `r_inj` (typically 50-100 pc)
- Energy per SN: $10^{51}$ ergs
- Mass ejection: 8.4 $M_\odot$ (optional)

### 3. AMR Refinement Criterion

Refinement based on normalized density gradient:

```cpp
Real d2 = SQR(u0(m,IDN,k,j,i+1) - u0(m,IDN,k,j,i-1));
if (multi_d) {d2 += SQR(u0(m,IDN,k,j+1,i) - u0(m,IDN,k,j-1,i));}
if (three_d) {d2 += SQR(u0(m,IDN,k+1,j,i) - u0(m,IDN,k-1,j,i));}
Real ddmax = sqrt(d2)/u0(m,IDN,k,j,i);

if (ddmax > ddens_threshold) {refine_flag = 1;}        // Refine
if (ddmax < 0.25*ddens_threshold) {refine_flag = -1;}  // Derefine
```

### 4. Boundary Condition Implementation

User boundaries apply the cooling flow solution at all domain edges:

```cpp
void UserBoundary(Mesh* pm) {
  // For each boundary face:
  SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
  SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
  u0(m, nhydro, k, j, i) = Z * Zsol * u0(m, IDN, k, j, i);  // Metallicity
}
```

## Input Parameters

### Required Parameters

```ini
<potential>
r_scale   = 20.0      # NFW scale radius [kpc]
rho_scale = 1e-25     # NFW scale density [g/cm^3]
mass_gal  = 1e10      # Galaxy mass [M_sun]
scale_gal = 3.0       # Disk scale length [kpc]
z_gal     = 0.3       # Disk scale height [kpc]
r_200     = 200.0     # Virial radius [kpc]
rho_mean  = 1e-28     # Mean density [g/cm^3]

<problem>
profile_file      = "cgm_profile.dat"      # CGM profile
disk_profile_file = "disk_profile.dat"     # Disk profile
r_circ = 10.0        # Rotation scale radius [kpc]
v_circ = 100.0       # Rotation velocity [km/s]
metallicity = 0.333  # Initial metallicity [solar]
ddens_max = 2.0      # AMR density gradient threshold

<SN>
r_inj = 0.1          # Injection radius [code units]
E_sn  = 1e51         # Energy per SN [ergs]
M_ej  = 8.4          # Ejecta mass [M_sun]
```

### Units Configuration

The code uses physical units via the Units class:
```cpp
Unit Length         : kpc
Unit Temperature    : K
Unit Number Density : cm^-3
Unit Velocity       : km/s
Unit Time           : Myr
```

## Algorithms and Numerical Methods

### 1. Interpolation Algorithm

Linear interpolation for profile lookups:

```cpp
KOKKOS_INLINE_FUNCTION 
Real ProfileReader::GetDensity(Real r) const {
  // Binary search for bracketing indices
  int i_low = 0, i_high = num_points - 1;
  while (i_high - i_low > 1) {
    int i_mid = (i_low + i_high) / 2;
    if (r < d_r_vals(i_mid)) {
      i_high = i_mid;
    } else {
      i_low = i_mid;
    }
  }
  
  // Linear interpolation
  Real t = (r - d_r_vals(i_low)) / (d_r_vals(i_high) - d_r_vals(i_low));
  return d_rho_vals(i_low) + t * (d_rho_vals(i_high) - d_rho_vals(i_low));
}
```

### 2. Force Calculation

Gravitational forces via finite differences:

```cpp
Real phi1l = GravPot(x1l, x2v, x3v, ...);
Real phi1r = GravPot(x1r, x2v, x3v, ...);
Real f_x1 = -(phi1r - phi1l)/(x1r - x1l);
```

Second-order accurate central differences

### 3. Hydrostatic Equilibrium Solver

For disk initialization:
1. Compute potential at cell center and faces
2. Calculate pressure gradient from density profile
3. Solve for rotation velocity from force balance
4. Apply vertical stratification

### 4. SN Injection Algorithm

1. **Timing**: Check each particle for next SN time
2. **Registration**: Atomic counter for concurrent SNe
3. **Injection**: Loop over cells within r_inj
4. **Update**: Advance particle to next SN time

## Performance Considerations

### GPU Optimization

All kernels use Kokkos lambdas for GPU execution:
```cpp
par_for("kernel_name", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
KOKKOS_LAMBDA(int m, int k, int j, int i) {
  // Kernel body
});
```

### Memory Management

- Profile data stored in device-accessible Kokkos Views
- Atomic operations for SN counting
- Team-based reductions for AMR criterion

### Load Balancing

AMR can create load imbalance:
- Disk region typically most refined
- SN injection sites temporarily refined
- Use appropriate decomposition strategy

## Usage Examples

### Basic CGM Simulation

```ini
# Basic cooling flow setup
<job>
problem_id = cgm_cooling

<mesh>
nx1 = 128
x1min = -50.0
x1max = 50.0
# ... (similar for x2, x3)

<meshblock>
nx1 = 32
nx2 = 32
nx3 = 32

<mesh_refinement>
adaptive = true
max_level = 3

<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
active_scalars = true
nscalars = 1

<cooling>
enable = true
```

### High-Resolution Disk

```ini
# Focus refinement on disk
<problem>
ddens_max = 0.5  # More aggressive refinement

<mesh_refinement>
max_level = 5    # Allow deeper refinement

# Use static refinement for disk region
<refinement1>
level = 3
x1min = -5.0
x1max = 5.0
x2min = -5.0
x2max = 5.0
x3min = -1.0
x3max = 1.0
```

### SN Feedback Study

```ini
# Enhanced feedback parameters
<SN>
r_inj = 0.2    # Larger injection region
E_sn = 2e51    # Higher energy
M_ej = 10.0    # More ejecta

<particles>
enable = true
particle_type = star
```

## Debugging and Validation

### Common Issues

1. **Profile extrapolation**: Check Rmax in profiles
2. **SN injection overflow**: Reduce e_sn or increase r_inj
3. **AMR thrashing**: Adjust ddens_threshold
4. **Boundary artifacts**: Verify profile continuity

### Validation Tests

1. **Hydrostatic equilibrium**: Turn off cooling, check drift
2. **Cooling flow**: Compare with 1D solution
3. **SN energy conservation**: Track total energy
4. **Metal conservation**: Check scalar conservation

### Diagnostic Outputs

Enable these derived variables:
```cpp
// In outputs configuration
<output1>
variable = prim
dt = 0.1
x3_slice = 0.0  # Midplane slice

<output2>
variable = user
user_vars = metallicity,sn_rate,cooling_rate
```

## Advanced Features

### 1. Multiple Metallicity Species

Extend to track individual elements:
```cpp
// Allocate multiple scalars
nscalars = 10;  // Fe, O, Si, etc.

// Track each species
u0(m, nhydro+i, k, j, i) = Z_i * Zsol_i * rho;
```

### 2. Turbulent Driving

Combine with turbulence module:
```cpp
#include "srcterms/turb_driver.hpp"

// In UserSource
TurbulenceDriver(pm, bdt);
```

### 3. Sink Particles

For star formation:
```cpp
// Check density threshold
if (rho > rho_sink) {
    CreateSinkParticle(x, v, m);
}
```

### 4. Non-Equilibrium Chemistry

Replace ideal gas with chemistry network:
```cpp
// Use chemical abundances
Real ne, nH, nHe;
ChemicalNetwork::Solve(T, rho, Z, ne, nH, nHe);
```

## Code Maintenance

### Adding New Physics

1. Extend source terms in `UserSource`
2. Add parameters to namespace
3. Update boundary conditions if needed
4. Consider AMR implications

### Performance Profiling

Key hotspots:
- Gravity calculation: ~30% runtime
- Profile interpolation: ~10% runtime  
- SN injection: Variable (0-20%)
- AMR operations: ~15% runtime

### Future Improvements

Planned enhancements:
- [ ] Magnetic fields (MHD version)
- [ ] Radiation transport
- [ ] Self-consistent cooling
- [ ] Cosmological expansion
- [ ] Subgrid turbulence model

## References

Key papers for this implementation:

1. **Cooling flows**: Stern et al. (2019, 2020)
2. **CGM structure**: Fielding et al. (2017, 2020)
3. **SN feedback**: Hopkins et al. (2018, FIRE-2), Hopkins et al. (2023, FIRE-3)
4. **AMR methods**: Berger & Colella (1989)
5. **Disk models**: Miyamoto & Nagai (1975)

## Contact

For questions about this implementation:
- Check the hidden documentation at `/docs/source/_hidden/cgm_metals/`
- Review the source code in `src/pgen/cgm_cooling_flow_amr_metals.cpp`
- Examine test cases in `inputs/cgm/`

---

*Last updated: 2024*  
*Version: AthenaK 1.0*  
*Status: Production Ready*