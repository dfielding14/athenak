# Turbulence Driving in AthenaK

This directory contains source term implementations for AthenaK, including CGM cooling and the turbulence driving method. The turbulence driver (`turb_driver.hpp/cpp`) implements Cartesian Fourier driving for simulating turbulence in computational domains.

## CGM Cooling

The CGM cooling source term uses the shared helper in `cgm_cooling.hpp` for the energy update, optional cooling timestep estimate, and `cooling_time` derived output. This keeps the radiative cooling, heating profile, and shielding mix in one place.

Cooling is operator split: `SourceTerms::CGMCoolingSplit` runs once per step in the `after_timeintegrator` task list (not inside the RK stages, where a stiff source loses half its effect to the RK averaging). In active cells, `cgm_cooling::IntegrateTemperature` integrates T exactly over dt with the net rate taken linear in T between table nodes (Townsend 2009, generalized as in Tan & Fielding 2023, App. C), clamped to [`tfloor`, `T_max`]. The task then refreshes restriction, ghost zones and primitives (`Driver::InitBoundaryValuesAndPrimitives`), so AMR and the next step's fluxes see the cooled state, and recomputes dt. The one-zone regression test is `tst/scripts/hydro/cgm_cooling_onezone.py`.

**The cooling rate depends on resolution.** `CalculateRates` blends the UV-background (PIE) and collisional (CIE) rates with the shielding optical depth of one cell, `tau = f_neutral * n_H * 1e-17 cm^2 * dx`, where `dx` is the local cell width. Neutral gas of the same density and temperature therefore cools differently at different resolutions and on either side of an AMR level boundary: one level finer halves `tau` and moves the rate toward the unshielded PIE limit. This is deliberate for now. The `cooling_time` output and the optional dt limiter use the same rate.

The explicit cooling timestep limiter is opt-in with `cgm_cooling_limit_dt = true`. When enabled, `cgm_cooling_dt_frac` sets the fraction of the local net-cooling time used for the source-term timestep. Gotham-scale production inputs keep this off by default because a single pathological cell can otherwise force a very small global timestep.

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
- `dedt`: Total power injected over the whole domain (energy/time in code units, not per volume). With `constant_edot = true` (default) each step's kick adds exactly `dedt*dt`
- `tcorr`: Correlation time for amplitude evolution
- `rseed`: Random seed for reproducibility (-1 for time-based)
- `sol_fraction`: Fraction of solenoidal modes (1.0 = incompressible)

### Spectral Parameters
- `nlow`, `nhigh`: Minimum and maximum |k| for mode selection
- `spect_form`: 1=parabolic spectrum, 2=power law
- `kpeak`: Peak wavenumber for parabolic spectrum
- `npeak`: Peak as a mode index, converted as `kpeak = npeak*2*pi`, i.e. as if the driving box (or tile) had length 1. For a box or tile of length `L` the peak lands at index `npeak*L`. The parabolic weight per mode is `|1 - 4(n - n_peak)^2/(nhigh - nlow)^2|`, an absolute value, so a misplaced peak also gives power to shells that should get none. GOTHAM's 400/384 kpc tiles put `npeak = 9` at `n = 9.375`; with `nlow = 8`, `nhigh = 10` the per-mode amplitudes of the |n| = 8, 9, 10 shells are 1.11, 0.97, 0.73, where a peak at 9 gives 0, 1, 0. Left as is so running campaigns keep their forcing; `development` converts with the tile length (a84c7ffd6)
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

2. **Update** (every `dt_turb_update`):
   - Evolve amplitudes using Ornstein-Uhlenbeck process
   - Ensure desired power spectrum
   - Remove net momentum to maintain momentum conservation

3. **Force Application** (every timestep):
   - Newtonian runs: one kick per step after the last RK stage's fluid update and source terms, normalized against the state being kicked; the energy update `dE = m.dv + rho|dv|^2/2` leaves the thermal energy unchanged
   - Special relativity: a source term in every RK stage

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
