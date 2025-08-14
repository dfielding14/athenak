# CGM Cooling Flow with Metals and Supernovae

## Overview

The CGM (Circumgalactic Medium) cooling flow with metals implementation simulates the cooling flow in galaxy halos with metal enrichment from supernovae. This module combines gravitational potential modeling, radiative cooling, AMR (Adaptive Mesh Refinement), and supernova feedback with metal injection.

## Source Files

### Primary Implementation
- **`src/pgen/cgm_cooling_flow_amr_metals.cpp`**: Main problem generator with AMR and metal tracking
- **Lines 1-850**: Complete implementation including initialization, source terms, and refinement

### Supporting Components
- **`src/utils/profile_reader.hpp`**: Lines 1-200 - Profile data interpolation for initial conditions
- **`src/utils/sn_scheduler.hpp`**: Lines 1-56 - Supernova timing based on FIRE-3 model
- **`src/srcterms/ismcooling.hpp`**: Lines 1-50 - ISM cooling function implementation

### Initial Condition Generation
- **`CGM_ICs_external/ipynb/generate_CGM_ics.ipynb`**: Jupyter notebook for generating steady-state cooling flow profiles

## Key Components

### 1. Gravitational Potential (`GravPot` function, lines 504-532)

The gravitational potential combines three components:

#### NFW (Navarro-Frenk-White) Halo
```cpp
// Line 518: NFW component
Real phi_NFW = -4 * M_PI * G * rho_s * SQR(r_s) * log(1 + x) / x;
```
Where `x = r/r_s` is the normalized radius.

#### Miyamoto-Nagai Disk
```cpp  
// Line 521: Disk component
Real phi_MN = -G * M_gal / sqrt(R*R + SQR(sqrt(x3*x3 + z_gal*z_gal) + a_gal));
```
Models a flattened galactic disk with scale length `a_gal` and scale height `z_gal`.

#### Outer Halo
```cpp
// Lines 524-526: Outer component for large radii
Real term1 = (4.0/3.0) * pow(5*R200, 1.5) * sqrt(r);
Real term2 = (1.0/6.0) * r2;
Real phi_Outer = 4 * M_PI * G * rho_mean * (term1 + term2);
```

### 2. Initial Conditions

#### Profile Reading (lines 171-197)
The code reads two profile files:
1. **CGM profile** (`profile_file`): Contains radial profiles of density, temperature, and velocity
2. **Disk profile** (`disk_profile_file`): Contains disk density distribution

```cpp
// Lines 174-176: Load CGM profiles
profile_reader_host.ReadProfiles(profile_file);
profile_reader = profile_reader_host.CreateDeviceReader();
```

#### State Initialization (lines 270-308)
Three functions set the initial state:

1. **`SetCoolingFlowState`** (lines 270-308): Sets density, temperature, and radial velocity from profiles
   ```cpp
   // Line 278-279: Interpolate from profiles
   Real rho = profile.GetDensity(r);
   Real temp = profile.GetTemperature(r);
   Real vr = profile.GetVelocity(r);
   ```

2. **`SetRotation`** (lines 368-417): Adds rotational velocity in hydrostatic equilibrium
   ```cpp
   // Line 401: Calculate circular velocity from pressure gradient
   Real v_phi = sqrt(R*fmax(dP_dR_over_rho - dPhi_dR, 0.0));
   ```

3. **`SetEquilibriumState`** (lines 310-366): Sets disk component in vertical hydrostatic equilibrium

#### Metallicity (line 258)
```cpp
// Uniform initial metallicity (Z/Z_solar = 1/3 by default)
u0(m, nhydro, k, j, i) = Z_ * Zsol * u0(m, IDN, k, j, i);
```
Where `Zsol = 0.02` (solar metallicity).

### 3. Cooling Mechanism

The cooling is implemented via the **ISMCoolFn** in `src/srcterms/ismcooling.hpp`:

#### Temperature Regimes (lines 31-48)
1. **T < 10^4 K**: Koyama & Inutsuka (2002) formula
   ```cpp
   // Line 32: Low temperature cooling
   return (2.0e-19*exp(-1.184e5/(temp + 1.0e3)) + 2.8e-28*sqrt(temp)*exp(-92.0/temp));
   ```

2. **10^4 K < T < 10^8.15 K**: SPEX cooling tables (Schure et al. 2009)
   ```cpp
   // Lines 42-47: Tabulated cooling rate interpolation
   Real logcool = (lhd[ipps+1]*dx - lhd[ipps]*(dx - 0.04))*25.0;
   ```

3. **T > 10^8.15 K**: CGOLS fit
   ```cpp
   // Line 37: High temperature cooling
   return pow(10.0, (0.45*logt - 26.065));
   ```

### 4. Supernova Feedback (`SNSource` function, lines 534-652)

#### SN Timing (lines 564-574)
Uses FIRE-3 model from `sn_scheduler.hpp`:
```cpp
// Line 573: Calculate next SN time for cluster
pr(nrdata-1, p) = GetNthSNTime(cluster_mass, par_t_create, unit_time, sn_idx);
```

The timing follows three phases:
1. **Early core-collapse** (3.7-7 Myr): High rate
2. **Late core-collapse** (7-44 Myr): Declining rate  
3. **Type Ia** (>44 Myr): Low constant rate

#### Energy and Mass Parameters (lines 131-136)
```cpp
// Line 133: Injection volume
const Real sphere_vol = (4.0/3.0)*M_PI*std::pow(r_inj,3);
// Line 135: Energy density (default 10^51 ergs per SN)
e_sn = E_def*pmbp->punit->erg()/sphere_vol;
// Line 136: Mass density (default 8.4 solar masses per SN)
m_ej = M_def*pmbp->punit->msun()/sphere_vol;
```

#### Injection Process (lines 605-642)
**Current Implementation Status**: Incomplete - only energy injection is active.

For each cell within injection radius:
```cpp
// Lines 628-633: Calculate distance from SN center
Real dx = x1v - sn_x;
Real dy = x2v - sn_y;
Real dz = x3v - sn_z;
Real r = sqrt(dx*dx + dy*dy + dz*dz);

// Lines 635-639: Current injection (incomplete)
if (r <= dr) {
  // u0(m,IDN,k,j,i) += m_ej_;           // Mass injection - COMMENTED OUT
  u0(m,IEN,k,j,i) += e_sn_;              // Energy injection - ACTIVE
  // u0(m,nhydro,k,j,i) += m_ej_;        // Metal injection - NOT IMPLEMENTED
}
```

**Implementation Status**: Currently only energy injection is active. Mass and metal injection are planned features (TODO) for the next version. The infrastructure for metal tracking exists (metallicity scalar field) and the code reads `M_ej` from input, but SNe do not yet enrich the gas with mass or metals.

### 5. AMR Refinement (`RefinementCondition`, lines 793-839)

The refinement is based on density gradient:

```cpp
// Lines 824-827: Calculate density gradient magnitude
Real d2 = SQR(u0(m,IDN,k,j,i+1) - u0(m,IDN,k,j,i-1));
if (multi_d) {d2 += SQR(u0(m,IDN,k,j+1,i) - u0(m,IDN,k,j-1,i));}
if (three_d) {d2 += SQR(u0(m,IDN,k+1,j,i) - u0(m,IDN,k-1,j,i));}

// Lines 831-832: Refinement decision
if (team_ddmax > ddens_thresh) {refine_flag.d_view(m+mbs) = 1;}      // Refine
if (team_ddmax < 0.25*ddens_thresh) {refine_flag.d_view(m+mbs) = -1;} // Coarsen
```

The threshold `ddens_threshold` is set via input parameter `ddens_max`.

### 6. Source Terms (`UserSource`, lines 423-428)

Combines two source terms applied each timestep:
```cpp
void UserSource(Mesh* pm, const Real bdt) {
  GravitySource(pm, bdt);  // Gravitational acceleration
  SNSource(pm, bdt);       // Supernova feedback
}
```

## Input Parameters

### Required Parameters

#### `<potential>` Block (lines 119-125)
```ini
r_scale   = 11.605   # NFW scale radius [code units]
rho_scale = 2.174    # NFW scale density [code units]
mass_gal  = 6547.4   # Galaxy mass [code units]
scale_gal = 1.469    # Disk scale radius [code units]
z_gal     = 0.15     # Disk scale height [code units]
r_200     = 144.77   # R200 radius [code units]
rho_mean  = 2.58e-5  # Mean density [code units]
```

#### `<problem>` Block (lines 126-139)
```ini
r_circ      = 10.0        # Rotation radius [code units]
v_circ      = 20.0        # Rotation velocity [code units]
metallicity = 0.333       # Initial Z/Z_solar
ddens_max   = 2.0         # Density gradient threshold for AMR
profile_file      = "cgm_profiles.txt"  # CGM profile data
disk_profile_file = "disk_profiles.txt" # Disk profile data
```

#### `<SN>` Block (lines 131-136)
```ini
r_inj = 0.5    # Injection radius [code units]
E_sn  = 1e51   # Energy per SN [ergs] (optional, default 1e51)
M_ej  = 8.4    # Ejecta mass per SN [solar masses] (optional, default 8.4)
```

## Profile File Format

Both profile files use the format (lines 172-197):
```
# r T rho v
1.0  1.39  0.135  0.0
2.0  1.41  0.068  0.5
...
```
Where:
- `r`: Radius [code units]
- `T`: Temperature [code units]
- `rho`: Density [code units]
- `v`: Radial velocity [code units]

## Units System

The code uses a configurable unit system (lines 105-115):
```cpp
Unit Length         : kpc
Unit Temperature    : K  
Unit Number Density : cm^-3
Unit Velocity       : km/s
Unit Time           : Myr
```

## CGM_ICs: Initial Condition Generator

The `CGM_ICs_external/ipynb/generate_CGM_ics.ipynb` notebook generates steady-state cooling flow profiles by:

1. **Setting up the gravitational potential** (cells 4-5):
   - NFW halo with concentration from M-c relation
   - Miyamoto-Nagai disk
   - Outer halo component

2. **Solving cooling flow equations** (cell 11):
   - Shoots from sonic radius outward and inward
   - Balances cooling, gravity, and pressure forces
   - Uses Wiersma cooling tables with metallicity

3. **Generating star clusters** (cell 9):
   - Power-law mass distribution
   - Random positions within specified radius
   - Assigns circular velocities

4. **Creating disk profiles** (cells 14-19):
   - Exponential surface density
   - Vertical hydrostatic equilibrium
   - Isothermal at 10^4 K

## Key Physics

### Cooling Flow Solution
The steady-state solution balances:
- **Cooling**: Radiative losses via metal-line cooling
- **Heating**: Compression from inflowing gas
- **Gravity**: Drives the inflow
- **Pressure**: Provides support against gravity

### Metal Enrichment
- Initial uniform metallicity (default Z = 1/3 Z_solar)
- SNe inject metals with each explosion
- Metals enhance cooling rates

### AMR Strategy
- Refines regions with steep density gradients
- Captures shocks and contact discontinuities
- Essential for resolving SN bubbles and cooling interfaces

## Usage Example

```ini
<job>
problem_id = cgm_cooling_metals

<mesh>
refinement  = adaptive
max_level   = 3

<hydro>
eos         = ideal
gamma       = 1.66666667

<problem>
profile_file      = cgm_profiles.txt
disk_profile_file = disk_profiles.txt
metallicity       = 0.333
ddens_max         = 2.0
r_circ            = 10.0
v_circ            = 20.0

<potential>
# Parameters from CGM_ICs output
r_scale   = 11.605
rho_scale = 2.174
mass_gal  = 6547.4
scale_gal = 1.469
z_gal     = 0.15
r_200     = 144.77
rho_mean  = 2.58e-5

<SN>
r_inj = 0.5
E_sn  = 1e51
M_ej  = 8.4

<particles>
type = star
pusher = rk4_gravity
particle_file = star_clusters.txt
```

## Implementation Notes

1. **Kokkos Parallelization**: All kernels use `KOKKOS_LAMBDA` for GPU compatibility
2. **Profile Interpolation**: Device-accessible views for efficient GPU access
3. **Boundary Conditions**: Custom user boundaries maintain cooling flow at domain edges
4. **Memory Management**: Explicit destructor calls for profile readers (lines 847-848)
5. **Particle Integration**: Star particles track SN timing and metal production

## References

- **Cooling Tables**: Schure et al. 2009, A&A 508, 751 (SPEX cooling)
- **Low-T Cooling**: Koyama & Inutsuka 2002 (molecular cooling)
- **SN Rates**: FIRE-3 model for stellar populations
- **Profile Generation**: Based on cooling flow solutions in gravitational potentials