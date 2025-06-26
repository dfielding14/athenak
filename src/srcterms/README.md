# Turbulence Driving Methods in AthenaK

This directory contains source term implementations for AthenaK, including turbulence driving methods. The turbulence driver (`turb_driver.hpp/cpp`) implements two distinct methods for driving turbulence in simulations:

1. **Cartesian Fourier Driving** (standard method)
2. **Spherical Fourier-Bessel (SFB) Driving** (new method)

## Cartesian Fourier Driving

### Mathematical Foundation

The Cartesian method drives turbulence by applying a stochastic force field constructed from a superposition of plane waves:

```
F(x,t) = ∑_k A_k(t) exp(i k·x)
```

where:
- `k` is the wavevector in Cartesian coordinates
- `A_k(t)` is the time-dependent complex amplitude
- The sum is over all modes satisfying `k_low ≤ |k| ≤ k_high`

### Implementation Details

1. **Mode Selection**: Modes are selected within a spherical shell in k-space:
   - Integer wavevectors: `k = (k_x, k_y, k_z)` where each component is an integer
   - Constraint: `k_low² ≤ k_x² + k_y² + k_z² ≤ k_high²`

2. **Force Computation**: For each mode, the force at position `x` is:
   ```
   F_mode = amplitude × [cos(k·x) + i sin(k·x)]
   ```
   
3. **Amplitude Evolution**: Amplitudes follow an Ornstein-Uhlenbeck process:
   ```
   dA_k/dt = -A_k/t_corr + σ dW/dt
   ```
   where `t_corr` is the correlation time and `dW` is a Wiener process.

4. **Key Features**:
   - Maintains solenoidal (divergence-free) forcing when configured
   - Precomputes sine/cosine values for efficiency
   - Uses global coordinates to ensure continuity across AMR boundaries

## Spherical Fourier-Bessel (SFB) Driving

### Mathematical Foundation

The SFB method expands the force field in spherical harmonics and Bessel functions:

```
F(r,θ,φ,t) = ∑_{n,l,m} A_{nlm}(t) j_l(k_n r) Y_lm(θ,φ)
```

where:
- `j_l` are spherical Bessel functions
- `Y_lm` are spherical harmonics
- `k_n` are discrete wavenumbers determined by boundary conditions
- `n` is the radial mode index, `l` is the angular degree, `m` is the azimuthal order

### Implementation Details

1. **Basis Functions**:
   - Radial: Spherical Bessel functions `j_l(k_n r)`
   - Angular: Real spherical harmonics `Y_lm(θ,φ)`
   - Wavenumbers from boundary conditions at `r_inner` and `r_outer`

2. **Mode Selection**:
   - Radial modes: `n_min ≤ n ≤ n_max`
   - Angular modes: `l_min ≤ l ≤ l_max`, `|m| ≤ l`
   - Total number of modes scales as O(n_max × l_max²)

3. **Vector Field Construction**:
   - Poloidal component: `F_pol = ∇ × (∇ × [Φ(r,θ,φ) r̂])`
   - Toroidal component: `F_tor = ∇ × [Ψ(r,θ,φ) r̂]`
   - Where Φ and Ψ are scalar potentials expanded in SFB basis

4. **Key Features**:
   - Naturally respects spherical geometry
   - Can enforce boundary conditions at inner/outer radii
   - Suitable for spherical shells and full spheres
   - Automatically divergence-free when using vector spherical harmonics

## Configuration Parameters

### Common Parameters
- `type`: "hydro" or "mhd" - determines which fluid variables to drive
- `dedt`: Energy injection rate (energy/time/volume)
- `tcorr`: Correlation time for amplitude evolution
- `rseed`: Random seed for reproducibility

### Cartesian-Specific Parameters
- `nlow`, `nhigh`: Minimum and maximum |k| for mode selection
- `min_kx`, `max_kx`: Range for k_x component (similarly for y,z)

### SFB-Specific Parameters
- `r_inner`, `r_outer`: Inner and outer radii of driven region
- `n_min`, `n_max`: Range of radial mode indices
- `l_min`, `l_max`: Range of spherical harmonic degrees
- `driving_pol`, `driving_tor`: Weights for poloidal/toroidal components

## Algorithm Flow

1. **Initialization**:
   - Select modes based on configuration
   - Precompute basis functions (sin/cos for Cartesian, Bessel/harmonics for SFB)
   - Initialize random amplitudes

2. **Update** (every `dt_turb_update`):
   - Evolve amplitudes using OU process
   - Ensure desired power spectrum

3. **Force Application** (every timestep):
   - Compute force field at each cell center
   - Add to momentum equation source terms
   - Add corresponding energy source term

## AMR Considerations

Both methods handle Adaptive Mesh Refinement (AMR) by:
- Using global coordinates for basis function evaluation
- Ensuring phase continuity across refinement boundaries
- Computing forces independently on each MeshBlock using global position

## Performance Notes

- Cartesian method: O(N_cells × N_modes) per update
- SFB method: O(N_cells × N_modes) but with more expensive basis functions
- Both methods benefit from precomputation and vectorization
- Memory usage scales with number of modes

## References

- Cartesian driving: Stone et al. (1998), Mac Low (1999)
- SFB driving: Inspired by geophysical turbulence studies and stellar convection simulations
- Implementation follows Schmidt et al. (2009) for stochastic forcing