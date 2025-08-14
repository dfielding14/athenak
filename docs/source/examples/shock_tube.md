# Example: Shock Tube (Sod Problem)

The classic Sod shock tube test - a fundamental benchmark for hydrodynamics codes demonstrating shock capturing, contact discontinuities, and rarefaction waves.

## Problem Generator
**Source**: `src/pgen/tests/shock_tube.cpp`

## Available Input Files
- **Standard Sod**: `inputs/hydro/sod.athinput`
- **MHD Brio-Wu**: `inputs/mhd/bw.athinput`
- **Relativistic**: `inputs/srhydro/mb1.athinput`
- **GRMHD**: `inputs/grmhd/mub1-gr.athinput`

## Physics

The shock tube problem consists of two uniform states separated by a discontinuity:

```
Initial Condition (t=0):
|  Left State  |  Right State |
|  ρ=1.0       |  ρ=0.125    |
|  P=1.0       |  P=0.1      |
|  v=0.0       |  v=0.0      |
```

This evolves into three waves:
1. **Shock wave** propagating right
2. **Contact discontinuity** moving right
3. **Rarefaction wave** propagating left

## Running the Simulation

### Basic Hydrodynamic Version
```bash
# Build with shock_tube problem
cmake -B build -DPROBLEM=shock_tube
make -C build -j8

# Run standard Sod test
./build/src/athena -i inputs/hydro/sod.athinput
```

### Complete Input File
```ini
# inputs/hydro/sod.athinput
<comment>
problem   = Sod shock tube
reference = Sod, G.A., JCP 27, 1-31 (1978)

<job>
basename  = sod        # problem ID

<mesh>
nghost    = 2
nx1       = 256        # resolution
x1min     = -0.5
x1max     = 0.5
ix1_bc    = outflow
ox1_bc    = outflow

nx2       = 1          # 1D problem
x2min     = -0.5
x2max     = 0.5
ix2_bc    = periodic
ox2_bc    = periodic

nx3       = 1
x3min     = -0.5
x3max     = 0.5
ix3_bc    = periodic
ox3_bc    = periodic

<meshblock>
nx1       = 256

<time>
evolution  = dynamic
integrator = rk2       # 2nd-order Runge-Kutta
cfl_number = 0.4
nlim       = -1
tlim       = 0.25      # end time
ndiag      = 1

<hydro>
eos         = ideal
reconstruct = plm      # piecewise linear
rsolver     = hllc     # HLLC Riemann solver
gamma       = 1.4      # ideal gas

<problem>
shock_dir  = 1         # 1,2,3 for x,y,z
xshock     = 0.0       # position of discontinuity

# Left state
dl         = 1.0       # density
pl         = 1.0       # pressure
ul         = 0.0       # velocity

# Right state
dr         = 0.125
pr         = 0.1
ur         = 0.0

<output1>
file_type  = vtk
variable   = prim      # output primitive variables
dt         = 0.025     # output interval
```

## MHD Version (Brio-Wu)

For magnetohydrodynamics, use the Brio-Wu shock tube:

```bash
# Run MHD version
./build/src/athena -i inputs/mhd/bw.athinput
```

Key differences:
- Includes magnetic field: `B_x = 0.75`, `B_y = ±1.0`
- Seven waves instead of three
- Tests constrained transport for div(B)=0

## Relativistic Versions

### Special Relativistic Hydrodynamics
```bash
./build/src/athena -i inputs/srhydro/mb1.athinput
```

### General Relativistic MHD
```bash
./build/src/athena -i inputs/grmhd/mub1-gr.athinput
```

## Analyzing Results

### Expected Solution Structure
```python
import numpy as np
import matplotlib.pyplot as plt

# Load VTK output
data = load_vtk("sod.block0.out1.00010.vtk")

# Plot density profile
plt.figure(figsize=(10, 6))
plt.plot(data.x, data.density, 'b-', label='Density')
plt.axvline(x=-0.3, color='r', linestyle='--', label='Rarefaction')
plt.axvline(x=0.0, color='g', linestyle='--', label='Contact')
plt.axvline(x=0.2, color='m', linestyle='--', label='Shock')
plt.xlabel('Position')
plt.ylabel('Density')
plt.legend()
plt.title('Sod Shock Tube at t=0.25')
```

## Convergence Testing

Test numerical convergence by varying resolution:

```bash
# Low resolution
sed -i 's/nx1 = 256/nx1 = 128/' inputs/hydro/sod.athinput
./athena -i inputs/hydro/sod.athinput
mv sod.*.vtk results_128/

# Medium resolution  
sed -i 's/nx1 = 128/nx1 = 256/' inputs/hydro/sod.athinput
./athena -i inputs/hydro/sod.athinput
mv sod.*.vtk results_256/

# High resolution
sed -i 's/nx1 = 256/nx1 = 512/' inputs/hydro/sod.athinput
./athena -i inputs/hydro/sod.athinput
mv sod.*.vtk results_512/
```

## Variations and Extensions

### Different Riemann Solvers
```ini
<hydro>
rsolver = llf     # Most diffusive, most robust
rsolver = hlle    # Good balance
rsolver = hllc    # Most accurate (default)
rsolver = roe     # Can be unstable
```

### Different Reconstruction
```ini
<hydro>
reconstruct = dc      # 1st order (very diffusive)
reconstruct = plm     # 2nd order (default)
reconstruct = ppm4    # 3rd order PPM (4th order accurate)
reconstruct = ppmx    # PPM with extrema preserving limiter
reconstruct = wenoz   # 5th order WENO-Z
```

### Multi-dimensional
Change to 2D by modifying:
```ini
<mesh>
nx2 = 256
<meshblock>
nx2 = 32
```

## Common Issues

### Oscillations Near Discontinuities
- Reduce CFL number: `cfl_number = 0.2`
- Use more diffusive solver: `rsolver = llf`
- Use lower-order reconstruction: `reconstruct = plm`

### Wrong Wave Speeds
- Check equation of state: `gamma = 1.4`
- Verify initial conditions match literature

## References

1. **Original Sod Problem**: Sod, G.A., "A Survey of Several Finite Difference Methods for Systems of Nonlinear Hyperbolic Conservation Laws", JCP 27, 1-31 (1978)

2. **MHD Version**: Brio, M. & Wu, C.C., "An Upwind Differencing Scheme for the Equations of Ideal Magnetohydrodynamics", JCP 75, 400-422 (1988)

3. **Relativistic Tests**: Marti & Müller, "Numerical Hydrodynamics in Special Relativity", Living Rev. Relativity 6, 7 (2003)

## See Also
- [Hydro Module Documentation](../modules/hydro.md)
- [Riemann Solvers](../modules/riemann_solvers.md)
- [Problem Generator API](../modules/pgen.md)