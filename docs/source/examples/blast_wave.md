# Example: Blast Wave

2D/3D blast wave problem demonstrating strong shocks, AMR capabilities, and both cylindrical and spherical explosions.

## Problem Generator
**Source**: `src/pgen/blast.cpp`

## Available Input Files
- **Hydro**: `inputs/hydro/blast_hydro.athinput`
- **Hydro with AMR**: `inputs/hydro/blast_hydro_amr.athinput`
- **MHD**: `inputs/mhd/blast_mhd.athinput`
- **MHD with AMR**: `inputs/mhd/blast_mhd_amr.athinput`

## Physics

A high-pressure region in pressure equilibrium expands into a low-pressure ambient medium, creating:
- Strong outward-propagating blast wave
- Contact discontinuity
- Optional: Magnetic field interaction (MHD)

Initial setup:
```
r < r_inner: P = P_high, ρ = ρ_high
r_inner < r < r_outer: smooth transition
r > r_outer: P = P_low, ρ = ρ_low
```

## Running the Simulation

```bash
# Configure & build once (blast is part of the built-in problem set)
cmake -S . -B build
cmake --build build -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)

# Run hydro blast
./build/src/athena -i inputs/hydro/blast_hydro.athinput
```

## Input File with AMR

```ini
# Blast with AMR
<problem>
inner_radius = 0.1     # Inner radius of blast
outer_radius = 0.125   # Outer radius (transition)
prat = 100.0          # Pressure ratio
drat = 1.0            # Density ratio
coordinates = cartesian

# High pressure region
pi_amb = 100.0        # Inner pressure
di_amb = 1.0          # Inner density

# Low pressure ambient
pn_amb = 1.0          # Ambient pressure
dn_amb = 1.0          # Ambient density

<mesh_refinement>
refinement = adaptive
num_levels = 3
dpres_max = 0.5       # Refine on pressure gradients
```

## Sedov-Taylor Comparison

For a strong point explosion the solution follows the Sedov–Taylor similarity law. The theoretical shock position is

```{math}
r_\text{shock}(t) = \xi\left(\frac{E_0 t^2}{\rho_0}\right)^{1/5},
```

with $\xi \approx 1.15$ for $\gamma = 1.4$. Measure the shock radius from the simulation output and compare it against this expression to verify energy conservation.
