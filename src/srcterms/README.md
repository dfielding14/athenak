# Turbulence Driving in AthenaK

This directory contains source term implementations for AthenaK, including the turbulence driving method. The turbulence driver (`turb_driver.hpp/cpp`) implements Cartesian Fourier driving for simulating turbulence in computational domains.

## Cartesian Fourier Driving

### Mathematical Foundation

The turbulence driver applies a stochastic force field constructed from a superposition of plane waves:

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

## Configuration Parameters

### Common Parameters
- `turb_flag`: 1=decaying turbulence, 2=continuously driven
- `constant_edot`: Defaults to `true`; `dedt` is total box power (energy/time).
  With `false`, the acceleration normalization is $\int |a|^2 dV/2 = \texttt{dedt}$.
  These are this branch's existing input conventions; no volume factor is added.
- `tcorr`: Correlation time for amplitude evolution
- `rseed`: RNG seed for the OU process (non-negative values give reproducible sequences; negative values fall back to the internal default seed = 1).
- `sol_fraction`: Amplitude-space blend between solenoidal (divergence-free) and compressive components of each Fourier mode (`1.0` = purely solenoidal, `0.0` = purely compressive).

### Spectral Parameters
- `nlow`, `nhigh`: Minimum and maximum |k| for mode selection
- `spect_form`: 1=parabolic spectrum, 2=power law
- `kpeak`: Peak wavenumber for parabolic spectrum
- `expo`: Power law exponent (if spect_form=2)

### Spatial Windowing (optional)
- `x_turb_scale_height`: Gaussian scale height in x-direction (-1 = no windowing)
- `x_turb_center`: Center of Gaussian window in x-direction
- (Similarly for y and z directions)

## Algorithm Flow

1. **Initialization**:
   - Select modes based on configuration
   - Precompute basis functions (sin/cos arrays)
   - Initialize random amplitudes

2. **Update** (when the OU cadence requires it):
   - Evolve the random forcing state before the fluid/PIC integrator.

3. **Newtonian force application** (once per eligible timestep, including cycle zero):
   - Finish the fluid/PIC update, including particle feedback and orbital remapping.
   - Remove the density-weighted mean acceleration and normalize using that state.
   - Apply one full-step kick, then use the existing boundary exchange and primitive
     recovery. Each fluid uses its own density, momentum, and EOS energy guard.

For $m=\rho v$ and $\Delta v=a\Delta t$, the exact kick is

$$
m'=m+\rho\Delta v,\qquad
E'=E+m\cdot\Delta v+\tfrac12\rho|\Delta v|^2.
$$

It preserves thermal energy during the kick. With `constant_edot=true`, the
integrated work is `dedt*dt`; zero `dedt` gives no kick. This replaces the previous
RK-stage source integration, which used the linear work $m\cdot a$ correctly as
an RHS but normalized against a different, full-step map and skipped cycle zero.
The new external forcing is first-order split from fluid/PIC evolution, including
VL2. It does not change the particle pusher, Hall induction, or deposited
particle-to-gas exchange ledger. The special-relativistic path is unchanged.

## AMR Considerations

The turbulence driver handles Adaptive Mesh Refinement (AMR) by:
- **Preserving mode amplitudes** during mesh refinement (critical!)
- Using global coordinates for basis function evaluation
- Ensuring phase continuity across refinement boundaries
- Computing forces independently on each MeshBlock using global position
- Never resetting the turbulence state during AMR events

The key insight is that AMR is a numerical change, not a physical event. The Ornstein-Uhlenbeck process must continue smoothly through refinement.

## Performance Notes

- Computational cost: O(N_cells × N_modes) per update
- Memory usage scales with number of modes
- Benefits from precomputation and vectorization
- Sine/cosine values are cached for efficiency

## References

- Cartesian driving: Stone et al. (1998), Mac Low (1999)
- Ornstein-Uhlenbeck process: Uhlenbeck & Ornstein (1930), Phys. Rev. 36, 823
- AMR with turbulence: Schmidt et al. (2009), A&A 494, 127
