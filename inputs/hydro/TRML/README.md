# TRML (Turbulent Radiative Mixing Layer) Simulations

This directory contains input files for Turbulent Radiative Mixing Layer (TRML) simulations.
These setups model the interaction between hot diffuse gas and cold dense gas phases
(e.g., circumgalactic medium and cool clouds), including:

- Shear-driven Kelvin-Helmholtz instability
- Radiative cooling with a power-law cooling function
- Thermal instability and turbulent mixing

## Input Files

### 1. `athinput.TRML_simple` - Basic TRML Setup
- Single-level mesh (no AMR)
- No tracer particles
- Good for initial testing and understanding the physics
- Parameters: chi=100, mach_rel_hot=0.5, xi=100

### 2. `athinput.TRML_amr` - TRML with Adaptive Mesh Refinement
- Up to 3 total AMR levels (base + 2 refinement levels), staged with `refine_level_times`
- Refinement based on supported temperature-band gas plus density-gradient or velocity criteria
- Refines near the mixing layer interface (peak temperature gas)
- Same physics parameters as simple version

### 3. `athinput.TRML_amr_particles` - TRML with AMR + Lagrangian MC Tracer Particles
- Same AMR setup as above
- Includes Monte-Carlo Lagrangian tracer particles (Genel+ 2013 method)
- Particles follow mass flux between cells stochastically
- Useful for tracking mass exchange between phases

## Key Parameters

### Dimensionless Numbers
- **chi** (`density_contrast`): Temperature/density contrast ratio = T_hot/T_cold = 100
- **mach_rel_hot**: Mach number in hot phase relative to interface = 0.5
- **xi**: Ratio of shear time to cooling time = t_shear/t_cool = 100
  - xi < 1: Cooling-dominated regime (cold gas grows)
  - xi > 1: Mixing-dominated regime
  - xi ~ 1: Balanced regime with shattering

### Cooling Function
Uses a piecewise power-law cooling function:
- `beta_lo = -2.0`: Low-temperature slope (below T_peak)
- `beta_hi = 3.0`: High-temperature slope (above T_peak)
- `T_peak_over_T_cold = 4.64`: Location of cooling peak relative to T_cold

### Frame Tracking
The simulation can apply Galilean velocity shifts to keep the mixing layer
centered in the domain. This prevents the interface from drifting into boundaries.
- Enable with `use_frame_tracking = true`

### AMR Refinement Criteria
- `density_ratio_threshold`: Refine if the fractional density gradient exceeds this
- `T_min_threshold`, `T_max_threshold`: Temperature support gate for mixing-layer gas
- `amr_min_temp_cells`, `amr_min_temp_fraction`: Minimum temperature-band support before refining
- `vel2_rms_threshold`: Refine if the three-component velocity RMS exceeds this; disabled unless > 0
- `amr_derefine_factor`: Hysteresis factor used before derefining
- `refine_level_times`: Time schedule for the allowed total levels, e.g. `0.0:1, 0.5:2, 1.0:3`

## Running the Simulations

```bash
# Simple version (no AMR, no particles)
./athena -i inputs/hydro/TRML/athinput.TRML_simple

# AMR version
./athena -i inputs/hydro/TRML/athinput.TRML_amr

# AMR + particles version
./athena -i inputs/hydro/TRML/athinput.TRML_amr_particles
```

## Output Analysis

- Use `vis/python/plot_slice.py` for 2D visualization
- Use `vis/python/read_prtcl_bin.py` for particle data (when using particles)
- Use `vis/python/analyze_particle_movement.py` for particle trajectory analysis

## References

- Fielding et al. (2020) - Multiphase gas in the circumgalactic medium
- Tan et al. (2021) - Radiative mixing layers
- Genel et al. (2013, MNRAS.435.1426G) - Monte-Carlo tracer particle method
