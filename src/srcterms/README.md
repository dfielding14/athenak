# Turbulence Driving in AthenaK

This directory contains source term implementations for AthenaK, including CGM cooling and the turbulence driving method. The turbulence driver (`turb_driver.hpp/cpp`) implements Cartesian Fourier driving for simulating turbulence in computational domains.

## CGM Cooling

The CGM cooling source term uses the shared helper in `cgm_cooling.hpp` for the energy update, optional cooling timestep estimate, and `cooling_time` derived output. This keeps the radiative cooling, heating profile, and shielding mix in one place.

Cooling is operator split: `SourceTerms::CGMCoolingSplit` runs once per step in the `after_timeintegrator` task list (not inside the RK stages, where a stiff source loses half its effect to the RK averaging). In every active cell `cgm_cooling::IntegrateTemperature` integrates T exactly over dt with the net rate taken linear in T between nodes, ten per table bin (Townsend 2009, generalized as in Tan & Fielding 2023, App. C), clamped to [`tfloor`, `T_max`]; restriction, ghost zones and primitives are then refreshed. The rate is split into a temperature part (`TempPart`), tabulated once at the nodes (`TabulateNodes`), and a cell part (`CellPart`), computed once per cell; `NetRate` combines them, and `CalculateRates` (used by the timestep limiter and outputs) is the same combination at any temperature. The one-zone regression test is `tst/scripts/hydro/cgm_cooling_onezone.py`.

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
- `dedt`: Total energy injection rate over the whole domain (energy/time, not per volume) with `constant_edot = true`, the default. With `constant_edot = false` it instead fixes Σ ½|a|²·dt summed over cells without volume weights, so the forcing amplitude then depends on resolution and dt.
- `tcorr`: Correlation time for amplitude evolution
- `rseed`: Seed of the forcing's random sequence (the ran2 initialization value). The default -1 is a fixed seed, so runs are reproducible; nothing is time-based.
- `sol_fraction`: Fraction of solenoidal modes (1.0 = incompressible)

### Spectral Parameters
- `nlow`, `nhigh`: Minimum and maximum mode number |n| for mode selection, counted in waves per tile (per box without tiling)
- `spect_form`: 1=parabolic spectrum, 2=power law
- `kpeak`: Peak wavenumber for parabolic spectrum, in code units (default 4π)
- `npeak`: Peak of the parabolic spectrum as a mode number, so kpeak = npeak·2π/tile_lx (the box length without tiling); overrides `kpeak`. With `spect_form = 1` the run stops unless `nhigh > nlow` and `nlow ≤ npeak ≤ nhigh`.
- `expo`: Power law exponent (if spect_form=2)

### Spatial Windowing (optional)
- `x_turb_scale_height`: Gaussian scale height in x-direction (-1 = no windowing)
- `x_turb_center`: Center of Gaussian window in x-direction
- (Similarly for y and z directions)

### Tiled Driving (optional)
- `tile_driving`: Repeat one forcing pattern on every tile (default false)
- `tile_factor`: Tiles per dimension (default 1), the default for `tile_nx`, `tile_ny` and `tile_nz`, which can also be set individually. Each must divide the root-grid `nx1`, `nx2`, `nx3`.

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
   - Newtonian hydro, MHD and ion-neutral: one kick per step, after the last stage's fluid update and source terms and before restriction and boundary exchange (for ion-neutral, after the implicit update). The forcing is normalized against the state being kicked, then dv = a·dt is added with dE = m·dv + ρ|dv|²/2 from the pre-kick momentum m, so the thermal energy is unchanged.
   - Special relativity keeps a source term in every stage.

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
