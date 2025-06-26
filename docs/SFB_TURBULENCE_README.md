# Spherical Fourier-Bessel (SFB) Turbulence Driving

## Overview

The Spherical Fourier-Bessel (SFB) turbulence driving method is designed for driving turbulence in systems with natural spherical symmetries, particularly useful for:
- Simulations with Adaptive Mesh Refinement (AMR)
- Cases where turbulence is needed only in the central region of the domain
- Spherical or quasi-spherical systems (e.g., star formation, galaxy clusters)

## Mathematical Foundation

### 1. Spherical Fourier-Bessel Basis

The SFB basis consists of divergence-free vector fields constructed from spherical Bessel functions and vector spherical harmonics. Each mode is characterized by three quantum numbers:
- `l` (angular momentum): l = 0, 1, 2, ..., lmax
- `m` (azimuthal): m = -l, -l+1, ..., l-1, l
- `n` (radial): n = 1, 2, ..., nmax

### 2. Vector Spherical Harmonics

We use magnetic-type vector spherical harmonics **B**ₗₘ, which are inherently divergence-free:

```
B_lm = ∇ × (r̂ ψ_lm)
```

where ψ_lm = jₗ(kₗₙr) Yₗₘ(θ,φ) is the scalar potential.

The components in spherical coordinates are:
- B_r = 0 (radial component always zero)
- B_θ = -(1/r sin θ) ∂ψ/∂φ
- B_φ = (1/r) ∂ψ/∂θ

### 3. Spherical Bessel Functions

The radial dependence uses spherical Bessel functions jₗ(x) with boundary conditions:

```
jₗ(kₗₙ r₀) = 0
```

where r₀ is the outer radius of the turbulent region. This gives discrete wavenumbers:

```
kₗₙ = xₗₙ / r₀
```

where xₗₙ is the n-th root of jₗ(x).

### 4. Velocity Field Construction

The total velocity field is a superposition of modes:

```
v(r,θ,φ) = Σ_lmn [a_lmn^(Re) B_lmn^(Re) + a_lmn^(Im) B_lmn^(Im)]
```

where a_lmn are random amplitudes drawn from a Gaussian distribution.

## Implementation Details

### Mode Selection

Modes are selected based on their wavenumber kₗₙ:
- Only modes with `kmin ≤ kₗₙ ≤ kmax` are included
- For a given l and n: kₗₙ = xₗₙ / r0_turb
- Total number of modes: ~lmax² × nmax (filtered by k-range)

### Power Spectrum

The amplitude of each mode is scaled according to the desired power spectrum:

1. **Power-law spectrum** (`spect_form = 2`):
   ```
   |a_lmn| ∝ k_ln^(-(expo+2)/2)
   ```
   This gives E(k) ∝ k^(-expo), e.g., expo = 5/3 for Kolmogorov turbulence.

2. **Parabolic spectrum** (`spect_form = 1`):
   ```
   |a_lmn| ∝ sqrt(1 - 4(k_ln - k_peak)²/(k_max - k_min)²)
   ```
   This peaks at k = k_peak.

### Regularization Near Origin

To avoid divergence at r → 0 (especially for grids ≥ 256³):

1. **l = 0 modes**: Set to zero near origin (they vanish naturally)
2. **l = 1 modes**: Use linear approximation jₗ(kr) ≈ kr/3
3. **Higher l modes**: Suppressed near origin

Regularization radius: r_min ≈ 2-3 grid cells

### Momentum Conservation

After generating the field, we ensure zero net momentum:
1. Calculate mean velocity inside r₀: <vᵢ>
2. Subtract mean from each component: vᵢ → vᵢ - <vᵢ>
3. Normalize to unit variance if requested

## Input Parameters

### Required Parameters

- `basis_type = 1`: Selects SFB basis (0 = Cartesian)
- `lmax`: Maximum angular momentum quantum number (typically 10-20)
- `nmax`: Maximum radial mode number (typically 5-10)
- `r0_turb`: Outer radius of turbulent region (in code units)

### Optional Parameters

- `kmin`: Minimum wavenumber (default: nlow × 2π/L)
- `kmax`: Maximum wavenumber (default: nhigh × 2π/L)
- `spect_form`: Power spectrum form (1 = parabolic, 2 = power-law)
- `expo`: Power-law exponent (default: 5/3)
- `kpeak`: Peak wavenumber for parabolic spectrum
- `dedt`: Energy injection rate
- `tcorr`: Correlation time for Ornstein-Uhlenbeck process
- `dt_turb_update`: Time between forcing updates

## Example Configuration

```
<turb_driving>
basis_type    = 1       # Use SFB basis
lmax          = 15      # Max l quantum number
nmax          = 10      # Max n quantum number
r0_turb       = 0.4     # Turbulence radius
nlow          = 2       # Min k ≈ 2×(2π/L)
nhigh         = 8       # Max k ≈ 8×(2π/L)
spect_form    = 2       # Power-law spectrum
expo          = 1.667   # Kolmogorov (5/3)
dedt          = 0.1     # Energy injection rate
tcorr         = 0.5     # Correlation time
dt_turb_update = 0.01   # Update interval
```

## Advantages Over Cartesian Driving

1. **Memory Efficiency**: For driving in r < r₀ « L, SFB requires far fewer modes than Cartesian
2. **Natural Spherical Symmetry**: No preferred Cartesian directions
3. **Exact Zero Divergence**: Vector spherical harmonics are analytically divergence-free
4. **Smooth Radial Cutoff**: Natural boundary at r = r₀

## Computational Considerations

### Performance

- Mode precomputation: O(lmax² × nmax × N³) one-time cost
- Per-timestep update: O(N_modes × N³_active) where N³_active ~ (r₀/L)³ × N³

### Memory Requirements

- Mode storage: ~5 × N_modes × N³_active × sizeof(Real)
- Typically 10-100× less memory than equivalent Cartesian driving

### Parallelization

- Mode computation is parallelized over grid points using Kokkos
- Each MPI rank computes modes for its local MeshBlocks
- No communication required during forcing updates

## Testing and Validation

### Divergence Check

The implementation ensures ∇·v = 0 to machine precision by construction.

### Power Spectrum Verification

Use the provided Jupyter notebook to:
1. Generate test fields
2. Compute power spectra
3. Verify scaling laws

### Momentum Conservation

The analysis script checks:
- Net momentum: <v> ≈ 0
- RMS velocity: normalized to unity if requested

## Common Issues and Solutions

1. **Divergence at origin (r → 0)**:
   - Solution: Implemented regularization for r < r_min
   - Set r_min ≈ 2-3 grid cells

2. **Poor mode coverage**:
   - Increase lmax and/or nmax
   - Adjust kmin/kmax range
   - Check mode distribution with test script

3. **Memory usage**:
   - Reduce lmax or nmax
   - Increase r0_turb to reduce active region
   - Use fewer MeshBlocks per rank

## References

1. Dubinski et al. (1995) - Turbulence in molecular clouds
2. Arfken & Weber - Mathematical Methods for Physicists (spherical harmonics)
3. Jackson - Classical Electrodynamics (vector spherical harmonics)