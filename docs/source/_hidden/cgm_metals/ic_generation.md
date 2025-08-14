# Initial Condition Generation with CGM_ICs

## Overview

The CGM_ICs repository ([https://github.com/zunyibrt/CGM_ICs](https://github.com/zunyibrt/CGM_ICs)) provides Python tools for generating 1D steady-state cooling flow solutions that serve as initial conditions for the `cgm_cooling_flow_amr_metals.cpp` problem generator in AthenaK.

## Package Structure

### Core Modules (`pysrc/`)

#### `cooling_flow.py`
Main module for integrating steady-state flow equations based on Stern et al. (2019, MNRAS 488, 2549).

Key components:
- **Base Classes**: `Cooling` and `Potential` interfaces
- **Integration Function**: `IntegrateFlowEquations()` for solving ODEs
- **Shooting Methods**: 
  - `shoot_from_R_sonic()` - Find transsonic solutions
  - `shoot_from_R_circ()` - Find subsonic solutions with rotation

Physical equations solved:
```python
# Momentum equation (steady-state Euler)
d(ln rho)/d(ln r) = (-t_cool/t_flow/gamma - v_c^2/c_s^2 + 2M^2)/(1 - M^2)

# Energy equation
d(ln T)/d(ln r) = t_cool/t_flow + (gamma-1) d(ln rho)/d(ln r)
```

Where:
- $t_{cool}/t_{flow}$ = cooling-to-flow time ratio
- $M$ = Mach number $(v/c_s)$
- $v_c$ = circular velocity
- $\gamma = 5/3$ (adiabatic index)

#### `HaloPotential.py` / `HaloPotential_new.py`
Gravitational potential implementations matching those in AthenaK:

1. **NFWPotential**: NFW dark matter halo
   ```python
   vc(r) = sqrt(4*pi*G*rho_s*r_s^2 * [ln(1+x) - x/(1+x)]/x)
   where x = r/r_s
   ```

2. **MiyamotoNagaiPotential**: Galactic disk
   ```python
   Phi = -G*M_gal/sqrt(R^2 + (sqrt(z^2 + z_gal^2) + a_gal)^2)
   ```

3. **OuterHaloPotential**: Extended halo component
   ```python
   Phi_outer = 4*pi*G*rho_mean * [(4/3)*(5*R200)^1.5 * sqrt(r) + r^2/6]
   ```

4. **CombinedPotential**: Sum of all components

#### `WiersmaCooling.py`
Wrapper for Wiersma et al. (2009) cooling tables with metallicity dependence.

Features:
- Collisional ionization equilibrium (CIE) cooling
- Metallicity scaling: $\Lambda(T, Z) = \Lambda_H(T) + (Z/Z_\odot)\Lambda_{metals}(T)$
- Temperature range: $10^4 - 10^{8.5}$ K
- Redshift evolution tables

#### `ClusterGenerator.py`
Star cluster generation for supernova feedback particles.

Parameters:
- Cluster mass function: $dN/dM \propto M^{-\alpha}$ ($\alpha = 2$)
- Mass range: $10^4 - 10^6$ $M_\odot$
- Spatial distribution: exponential disk
- Velocity profile: derived from circular velocity

### Notebooks (`ipynb/`)

#### `generate_CGM_ics.ipynb`
Main workflow for creating initial conditions:

1. **Set Cosmology**: Planck18 parameters
2. **Define Halo Properties**:
   ```python
   M_vir = 1e11  # Halo mass [M_sun]
   r_vir = 122   # Virial radius [kpc]
   c_vir = 10.5  # Concentration
   M_gal = 1e10  # Galaxy mass [M_sun]
   ```

3. **Configure Units** (matching AthenaK):
   ```python
   cu_length = kpc
   cu_temperature = K
   cu_density = g/cm^3  # for n=0.1 cm^-3 at 1 kpc^3
   cu_velocity = km/s
   cu_time = Myr
   ```

4. **Generate Cooling Flow**:
   ```python
   sol = shoot_from_R_sonic(potential, cooling, 
                           R_sonic=100*kpc, 
                           R_max=230*kpc, 
                           R_min=1*kpc)
   ```

5. **Compute Disk Profile**:
   - Surface density: $\Sigma(R) = \Sigma_0 \exp(-R/R_{gal})$
   - Vertical structure: $\rho(z) \propto \exp(-\Phi(z)/c_s^2)$
   - Temperature: $T_{disk} = 10^4$ K

6. **Export Profiles**:
   ```
   # r[code] T[code] rho[code] v[code]
   0.001  13.9  1.5e-5  -0.1
   ...
   ```

#### `steady_state_integration_example.ipynb`
Tutorial demonstrating the integration methodology and parameter space exploration.

### Cooling Tables (`cooling/`)

Pre-computed cooling functions:
- **Wiersma09_CoolingTables/**: HDF5 files for different redshifts
- **DopitaSutherland_CIE.dat**: Alternative CIE cooling
- **GnatSternberg07_CIE.dat**: High-temperature cooling
- **Kartick_CIE_cooling.table**: Custom cooling function

## Integration with AthenaK

### Profile File Format

The generated profiles use a simple ASCII format:
```
# r T rho v
1.0000e+00  1.3922e+01  1.5432e-05  -1.0234e-01
1.0100e+00  1.3918e+01  1.5398e-05  -1.0198e-01
...
```

Where all quantities are in code units:
- `r`: Radius [kpc]
- `T`: Temperature [K/7181.24]
- `rho`: Density [g/cm$^3$/1.53e-34]
- `v`: Radial velocity [km/s/9.78]

### Loading in AthenaK

The `cgm_cooling_flow_amr_metals.cpp` problem generator reads these profiles:

```cpp
// In parameter file
<problem>
profile_file = cgm_profile.dat
disk_profile_file = disk_profile.dat

// In code
ProfileReaderHost profile_reader_host;
profile_reader_host.ReadProfiles(profile_file);
```

### Workflow Example

1. **Generate ICs with Python**:
   ```bash
   cd CGM_ICs_external/ipynb
   jupyter notebook generate_CGM_ics.ipynb
   # Run all cells
   # Files saved to tables/
   ```

2. **Copy profiles to simulation directory**:
   ```bash
   cp tables/cgm_profiles.txt ../../../inputs/cgm/
   cp tables/disk_profiles.txt ../../../inputs/cgm/
   ```

3. **Run AthenaK simulation**:
   ```bash
   ./athena -i inputs/cgm/cgm_cooling_flow.athinput
   ```

## Physical Model Details

### Steady-State Assumptions

The cooling flow solutions assume:
1. **Spherical symmetry** (for CGM)
2. **Time-steady flow**: $\partial/\partial t = 0$
3. **Mass conservation**: $\dot{M} = 4\pi r^2 \rho v = \text{const}$
4. **Momentum balance**: $v(dv/dr) = -c_s^2(d\rho/dr)/\rho - GM/r^2$
5. **Energy balance**: Cooling = Compression heating

### Solution Types

#### Transsonic Solutions
- Start subsonic at large radii
- Pass through sonic point ($M = 1$)
- Become supersonic at small radii
- Physically represent cooling flows

#### Subsonic Solutions
- Remain subsonic throughout
- Can include rotation: $v_\phi = v_c(R/R_{circ})\sin(\theta)$
- Represent hot halos in quasi-equilibrium

#### Breeze Solutions
- Very slow inflow ($v \ll c_s$)
- Nearly hydrostatic
- $t_{cool} \gg t_{flow}$

### Parameter Space

Key dimensionless parameters:
1. **$t_{cool}/t_{flow}$**: Cooling efficiency
2. **$v_c^2/c_s^2$**: Gravitational vs. thermal energy
3. **$R_{circ}/r$**: Angular momentum distribution

Typical values for MW-like halos:
- $R_{sonic} \sim 50-150$ kpc
- $\dot{M} \sim 0.1-1$ $M_\odot$/yr
- $T_{vir} \sim 10^6$ K
- $\rho(r_{vir}) \sim 10^{-27}$ g/cm$^3$

## Advanced Features

### Multiple Metallicity Profiles

Generate profiles with varying metallicity:
```python
for Z in [0.1, 0.3, 1.0]:  # Solar units
    cooling = Cool.Wiersma_Cooling(Z2Zsun=Z, z=0)
    sol = shoot_from_R_sonic(potential, cooling, ...)
    save_profile(f"cgm_Z{Z}.dat", sol)
```

### Non-Equilibrium Disk

Include disk winds or fountains:
```python
# Add vertical velocity component
v_z = v_wind * exp(-R/R_wind) * (z/z_wind)
```

### Time-Dependent Solutions

For studying transient phenomena:
```python
# Vary mass accretion rate
for t in times:
    Mdot = Mdot_0 * (1 + A*sin(2*pi*t/P))
    sol = solve_with_Mdot(Mdot)
```

## Validation and Testing

### Convergence Tests

Check numerical convergence:
```python
max_steps = [0.1, 0.05, 0.01, 0.005]
for step in max_steps:
    sol = shoot_from_R_sonic(..., max_step=step)
    plot_solution(sol)
```

### Comparison with Analytics

For isothermal potential ($v_c \propto r^0$):
- Analytical: $\rho \propto r^{-2}$
- Check: `plot(r, rho*r^2)` should be flat

### Conservation Checks

Verify conserved quantities:
```python
# Bernoulli parameter
B = 0.5*v**2 + cs**2/(gamma-1) + Phi
assert np.std(B)/np.mean(B) < 0.01

# Mass flux
Mdot = 4*pi*r**2*rho*v
assert np.std(Mdot)/np.mean(Mdot) < 0.01
```

## Common Issues and Solutions

### Issue 1: No Transsonic Solution Found

**Problem**: `shoot_from_R_sonic` fails to converge

**Solutions**:
- Adjust $R_{sonic}$ (try $0.5-2.0 \times r_{vir}$)
- Modify cooling function normalization
- Check potential is smooth at $R_{sonic}$

### Issue 2: Unphysical Disk Densities

**Problem**: Disk midplane density exceeds CGM by $10^6$

**Solutions**:
- Reduce gas fraction: $M_{gas}/M_{gal}$
- Increase disk temperature
- Add turbulent pressure support

### Issue 3: Profile Extrapolation Errors

**Problem**: AthenaK crashes at boundaries

**Solutions**:
- Extend profile range ($R_{min} \to 0.1$ kpc, $R_{max} \to 5 \times r_{vir}$)
- Add buffer zones with constant extension
- Smooth profile endpoints

## Performance Optimization

### Adaptive Step Size

Use scipy's adaptive ODE solvers:
```python
from scipy.integrate import solve_ivp
sol = solve_ivp(odes, [ln_R0, ln_Rf], y0, 
                method='DOP853',  # 8th order
                rtol=1e-8, atol=1e-10)
```

### Parallel Parameter Studies

Generate multiple profiles simultaneously:
```python
from multiprocessing import Pool

def generate_profile(params):
    R_sonic, cooling = params
    return shoot_from_R_sonic(...)

with Pool() as pool:
    results = pool.map(generate_profile, parameter_list)
```

### Caching Cooling Tables

Pre-load cooling functions:
```python
cooling_cache = {}
for Z in metallicities:
    cooling_cache[Z] = Cool.Wiersma_Cooling(Z2Zsun=Z)
```

## Repository Information

- **Location**: External to AthenaK (in `.gitignore`)
- **Clone**: `git clone https://github.com/zunyibrt/CGM_ICs CGM_ICs_external`
- **Dependencies**: numpy, scipy, astropy, matplotlib
- **Python Version**: 3.7+

## References

Key papers for the IC generation methodology:

1. **Cooling flows**: Stern et al. (2019, MNRAS 488, 2549)
2. **CGM structure**: Fielding et al. (2017, MNRAS 466, 3810)
3. **Cooling functions**: Wiersma et al. (2009, MNRAS 393, 99)
4. **Disk models**: Mo, Mao & White (1998, MNRAS 295, 319)
5. **Numerical methods**: Mathews & Bregman (1978, ApJ 224, 308)

---

*This documentation is part of the hidden CGM metals implementation guide.*
*Last updated: 2024*