# Configuration

AthenaK uses `.athinput` files for runtime configuration.

## Input File Structure

```ini
<job>
basename = simulation_name

<mesh>
nx1 = 256        # Grid resolution
x1min = 0.0      # Domain boundaries
x1max = 1.0
ix1_bc = periodic # Boundary conditions
ox1_bc = periodic

<time>
cfl_number = 0.4  # CFL stability factor
tlim = 10.0       # End time

<hydro>
gamma = 1.4       # Adiabatic index
reconstruct = plm # Reconstruction method
rsolver = hllc    # Riemann solver
```

## Common Parameters

### Mesh Block
- `nx1, nx2, nx3` - Resolution in each direction
- `x1min, x1max` - Domain extent
- `ix1_bc, ox1_bc` - Inner/outer boundary conditions

### Time Integration
- `cfl_number` - CFL factor (0.1-0.4)
- `tlim` - Simulation end time
- `nlim` - Maximum number of cycles
- `integrator` - rk1, rk2, rk3, rk4, imex2, imex3, imex+

### Physics Options

#### Hydrodynamics
- `gamma` - Adiabatic index (ideal gas EOS)
- `eos` - ideal, isothermal
- `reconstruct` - dc, plm, ppm4, ppmx, wenoz
- `rsolver` - llf, hlle, hllc, roe, advect

#### MHD
- `eos` - ideal, isothermal
- `rsolver` - llf, hlle, hlld, advect
- `viscosity` - Kinematic viscosity coefficient
- `conductivity` - Thermal conductivity coefficient
- `ohmic_resistivity` - Resistivity coefficient

### Boundary Conditions
- `periodic` - Periodic
- `outflow` - Zero gradient
- `reflect` - Reflection
- `inflow` - Fixed inflow
- `diode` - One-way flow
- `vacuum` - Vacuum BC
- `shear_periodic` - Shearing box
- `user` - User-defined

## Output Configuration

```ini
<output1>
file_type = vtk
dt = 0.1
variable = cons  # or prim

<output2>
file_type = bin  # binary format
dt = 1.0
single_file_per_rank = true  # for parallel I/O
```

## Example Configurations

### Shock Tube
```ini
<problem>
shock_dir = 1    # Direction (1=x, 2=y, 3=z)
xshock = 0.5     # Discontinuity position
dl = 1.0         # Left density
pl = 1.0         # Left pressure
dr = 0.125       # Right density
pr = 0.1         # Right pressure
```

### Blast Wave
```ini
<problem>
pradius = 0.1     # Blast radius
pamb = 1.0        # Ambient pressure
rhoamb = 1.0      # Ambient density
prat = 100.0      # Pressure ratio
```