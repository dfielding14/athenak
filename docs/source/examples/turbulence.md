# Example: Driven Turbulence

Compressible turbulence with continuous Fourier forcing. This example demonstrates AMR with turbulence and the recent fixes for div(B) preservation. The driver also supports optional “tiled” forcing that repeats the pattern across sub-domains.

## Problem Generators
- **Standard Turbulence**: `src/pgen/turb.cpp`
- **AMR Test Version**: `src/pgen/turb_amr_test.cpp`

## Available Input Files

- **Hydro (baseline)**: `inputs/hydro/turb.athinput`
- **Hydro with staged AMR demo**: `inputs/hydro/turb_amr_demo.athinput`
- **MHD staged refinement (part 1)**: `inputs/turb_timed_amr_stage1.athinput`

## Physics

Driven turbulence maintains a statistically steady state by continuously injecting energy at large scales, which cascades to small scales where it dissipates numerically.

### Fourier Driving (Cartesian)
Energy injection in Fourier space:
$$\mathbf{F}(\mathbf{k}) = A(\mathbf{k}) \cdot P(\mathbf{k}) \cdot \exp(i\phi(\mathbf{k}))$$

where:
- $A(\mathbf{k})$ = amplitude profile
- $P(\mathbf{k})$ = projection (solenoidal/compressive)
- $\phi(\mathbf{k})$ = random phases

## Running Basic Turbulence

```bash
# Configure & build once (turbulence drivers are included by default)
cmake -S . -B build
cmake --build build -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)

# Run standard turbulence
./build/src/athena -i inputs/hydro/turb.athinput
```

### Optional Tiled Driving

Set `tile_driving = true` to repeat the forcing pattern on smaller tiles of the domain. Each tile count must evenly divide the corresponding global mesh dimension (the driver checks this and aborts if the constraint is violated). Example:

```ini
<turb_driving>
tile_driving = true
tile_factor  = 2        # default replicated to tile_nx/tile_ny/tile_nz
tile_nx      = 2        # optional explicit tile counts
tile_ny      = 1
tile_nz      = 1
```

For 2-D problems (`nx3 <= 1`) the driver automatically resets `tile_nz = 1`, and likewise for `tile_ny` in 1-D setups.

## Complete Input File with AMR

```ini
# inputs/hydro/turb_amr_demo.athinput
<comment>
problem   = Driven turbulence with AMR
reference = Recent AMR+div(B) fixes

<job>
basename  = turb_amr

<mesh>
nghost    = 3          # Need 3 for PPM+AMR
nx1       = 64
x1min     = -0.5
x1max     = 0.5
ix1_bc    = periodic
ox1_bc    = periodic

nx2       = 64
x2min     = -0.5
x2max     = 0.5
ix2_bc    = periodic
ox2_bc    = periodic

nx3       = 64
x3min     = -0.5
x3max     = 0.5
ix3_bc    = periodic
ox3_bc    = periodic

<mesh_refinement>
refinement       = adaptive
num_levels       = 3        # 3 levels total
ncycle_check     = 10       # Check every 10 cycles
refinement_interval = 10

# AMR criteria
prolong_primitives = true
dens_max = 3.0             # Refine if $\rho$ > 3
ddens_max = 0.5            # Refine if $|\nabla\rho|/\rho$ > 0.5
dvel_max = 0.5             # Refine if $|\nabla v|$ > 0.5

<meshblock>
nx1       = 16         # Smaller blocks for AMR
nx2       = 16
nx3       = 16

<time>
evolution  = dynamic
integrator = rk3       # RK3 for stability
cfl_number = 0.3
nlim       = -1
tlim       = 10.0
ndiag      = 1

<hydro>
eos         = ideal
reconstruct = ppm4     # 3rd order PPM
rsolver     = hllc
gamma       = 1.001    # Nearly isothermal
fofc        = true     # First-order flux correction

<problem>
rho0        = 1.0      # Mean density
press0      = 1.0      # Mean pressure

<turb_driving>
turb_flag        = 2      # 1=finite-length forcing, 2=continuous
dedt             = 1.0    # Energy injection rate
tcorr            = 0.5    # Correlation time
nlow             = 2      # Minimum forced wavenumber
nhigh            = 4      # Maximum forced wavenumber
spect_form       = 1      # Spectrum shape
driving_type     = 0      # 0=solenoidal, 1=compressive
sol_fraction     = 1.0    # Solenoidal fraction
constant_edot    = true   # Maintain constant energy injection
dt_turb_update   = 0.01   # Update interval
rseed            = 42     # Random seed

<output1>
file_type  = hst       # History file
dt         = 0.01

<output2>
file_type  = vtk
variable   = prim
dt         = 0.1

<output3>
file_type  = vtk
variable   = turb_force  # Output forcing field
dt         = 0.5

### MHD Staged AMR Run

The input `inputs/turb_timed_amr_stage1.athinput` evolves an isothermal MHD
cubic domain with diffusion coefficients and FOFC enabled. It is intended as
the first leg of a multi-stage refinement workflow: stage 1 runs to `t=2`, writes
restart and diagnostic outputs, and subsequent stages (not included here) can be
started from the generated checkpoint.
```

## MHD Turbulence with div(B) Preservation

Recent fixes ensure div(B)=0 is maintained with AMR:

```ini
# inputs/mhd/test_divb_minimal.athinput
<mhd>
eos         = ideal
reconstruct = plm     # More stable for MHD+AMR
rsolver     = hlld    # HLLD for MHD
gamma       = 1.4

# Recent div(B) fixes for AMR
<mesh_refinement>
# Use constrained transport prolongation
prolong_primitives = false  # Prolong conserved for B-field
```

## Analysis and Diagnostics

### Power Spectrum
```python
import numpy as np
from scipy import fft

def compute_spectrum(data):
    """Compute kinetic energy power spectrum"""
    # FFT of velocity field
    vk = fft.fftn(data.velocity)

    # Power spectrum
    Ek = 0.5 * np.abs(vk)**2

    # Shell average
    k_bins = np.arange(0, kmax)
    spectrum = shell_average(Ek, k_bins)

    return k_bins, spectrum

# Expected: E(k) ∝ k^(-5/3) for Kolmogorov
```

### Monitoring div(B) with AMR
```python
def check_divB(bx, by, bz, dx):
    """Check divergence of B"""
    divB = np.gradient(bx, dx, axis=0)
    divB += np.gradient(by, dx, axis=1)
    divB += np.gradient(bz, dx, axis=2)

    max_divB = np.max(np.abs(divB))
    rms_divB = np.sqrt(np.mean(divB**2))

    print(f"Max |div(B)|: {max_divB:.2e}")
    print(f"RMS div(B): {rms_divB:.2e}")

    return divB
```

## Recent Fixes and Improvements

### AMR + Turbulence Fixes (2024)
1. **Forcing consistency**: Energy injection maintained across refinement levels
2. **Phase coherence**: Random phases preserved during regridding
3. **Load balancing**: Improved for turbulent flows

### div(B) Preservation Fixes
1. **Constrained transport**: Properly implemented at refinement boundaries
2. **Flux correction**: Ensures conservation at coarse-fine interfaces
3. **Prolongation**: Magnetic field prolongation preserves divergence

## Performance Optimization

### GPU Acceleration
```bash
# Build for GPU
cmake -S . -B build-cuda \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ARCH_AMPERE80=ON
cmake --build build-cuda -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)

# Optimal MeshBlock size for GPU
<meshblock>
nx1 = 32  # Larger blocks for GPU efficiency
nx2 = 32
nx3 = 32
```

### Scaling Tests
```bash
# Weak scaling
mpirun -np 1 ./athena -i turb_64.athinput
mpirun -np 8 ./athena -i turb_128.athinput
mpirun -np 64 ./athena -i turb_256.athinput

# Strong scaling
for np in 1 2 4 8 16 32 64; do
    mpirun -np $np ./athena -i turb_256.athinput
done
```

## Common Issues and Solutions

### Energy Decay
- Check `dedt` parameter is set correctly
- Verify `tcorr` is appropriate for flow timescale
- Ensure forcing wavenumbers (`nlow`, `nhigh`) are in inertial range

### AMR Oscillations
- Use `fofc = true` for first-order flux correction
- Reduce `refinement_interval`
- Adjust refinement criteria thresholds

### div(B) Growth
- Ensure using latest code with div(B) fixes
- Use `prolong_primitives = false` for MHD
- Check boundary conditions are compatible

## Validation Tests

### Kolmogorov Spectrum
Expected: `E(k) ∝ k^(-5/3)` in inertial range

### Sonic Mach Number
```
M_s = v_rms / c_s ≈ 1-10 for supersonic turbulence
```

### Magnetic Energy (MHD)
```
E_mag / E_kin ≈ 0.5-2.0 for equipartition
```

## References

1. **Compressible Turbulence**: Federrath, C., "The Physics of Supersonic Turbulence", MNRAS (2013)

2. **AMR+Turbulence**: Schmidt, W. et al., "Turbulence Production in AMR Simulations", A&A (2009)

## See Also
- [Source Terms Module](../modules/srcterms.md)
- [AMR Documentation](../modules/mesh.md)
- [MHD Module](../modules/mhd.md)
